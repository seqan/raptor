rule bifrost_batch_update:
    input:
        initial_files = "input_data/initial_files.txt",
        batches_dir = "input_data/batches"
    output:
        bench = "bifrost/{args}/bench.tsv"
    benchmark:
        repeat("bifrost/{args}/batch_update.bench.tsv", config["repeat"])
    log:
        "bifrost/{args}/batch_update.log"
    #conda:
        #("../envs/bifrost.yaml")
    threads:
        config["threads"]
    params:
        output_directory = "bifrost/{args}",
        bifrost_path = config["tools"]["bifrost"]["path"],
        ubt_cmd = lambda wildcards, output: get_usr_bin_time_cmd(output.bench, "build"),
        ubt_cmd_qry = lambda wildcards, output: get_usr_bin_time_cmd(output.bench, "query"),
        size_cmd = config["size_cmd"],
        query_file = config["query_file"],
        query_args = lambda wildcards: config["tools"]["bifrost"]["query_args"],
        args = lambda wildcards: str2args(wildcards.args),
        fixed_args = config["tools"]["bifrost"]["fixed_args"],
        initial_batch = config["initial_batch"] if "initial_batch" in config else 0,
        incremental_batch = config["incremental_batch"] if "incremental_batch" in config else 0,
        first_batch = config["first_batch"] if "first_batch" in config else 0,
        last_batch = config["last_batch"] if "last_batch" in config else "",
        scripts_path = os.path.join(workflow.basedir, "scripts"),
    shell: 
        """
        if [[ {params.first_batch} -eq 0 ]]; then
            BATCH=0

            {params.ubt_cmd} {params.bifrost_path} build \
                                            --input-ref-file "{input.initial_files}" \
                                            --output-file {params.output_directory}/${{BATCH}}.index \
                                            --threads {threads} \
                                            --verbose --colors {params.fixed_args} {params.args} > {log} 2>&1

            {params.size_cmd} {params.output_directory}/${{BATCH}}.index.* | sed "s/^/${{BATCH}}\tsize\t/" >> {output.bench}
            echo "${{BATCH}}\tobs\tInitial Build" >> {output.bench}
            first_batch=1
        else
            first_batch={params.first_batch}
        fi

        if [[ -z "{params.last_batch}" ]]; then
            last_batch=$(ls {input.batches_dir}/insert_batch_* | wc -l)
        else
            last_batch={params.last_batch}
        fi

        for BATCH in $(seq ${{first_batch}} ${{last_batch}}); do
            echo "${{BATCH}}" >> {log}
            
            input="{params.output_directory}/$((BATCH - 1)).index"
            output="{params.output_directory}/${{BATCH}}.index"
            
            {params.ubt_cmd} {params.bifrost_path} update --input-graph-file "${{input}}.gfa" \
                                                --input-color-file "${{input}}.color.bfg" \
                                                --input-index-file "${{input}}.bfi" \
                                                --input-ref-file "{input.batches_dir}/insert_batch_${{BATCH}}.txt" \
                                                --output-file "${{output}}" \
                                                --threads {threads} \
                                                --no-compress-out \
                                                --verbose {params.fixed_args} {params.args} >> {log} 2>&1

            #SIZE
            {params.size_cmd} ${{output}}.* | sed "s/^/${{BATCH}}\tsize\t/" >> {output.bench}

            # QUERY
            if [[ -n "{params.query_file}" ]]; then
                {params.ubt_cmd_qry} {params.bifrost_path} query --input-graph-file "${{output}}.gfa" \
                                                                 --input-color-file "${{output}}.color.bfg" \
                                                                 --input-index-file "${{output}}.bfi" \
                                                                 --input-query-file "{params.query_file}"  \
                                                                 --output-file "{params.output_directory}/${{BATCH}}.query" \
                                                                 --threads {threads} \
                                                                 {params.query_args} --verbose >> {log} 2>&1 

                awk 'NR>1{{m=0; for(i=1;i<=NF;i++){{if($i==1){{print $1"\t"i-2; m=1}}}};if(m==0){{print $1"\t"}}}}' "{params.output_directory}/${{BATCH}}.query.tsv" > "{params.output_directory}/${{BATCH}}.eval"
                
                # AVG. QUERY HITS
                total_queries=$(wc -l "{params.query_file}"| cut -d' ' -f 1)
                total_matches=$(cut -f 2 "{params.output_directory}/${{BATCH}}.eval" | sed -r '/^\s*$/d' | wc -l)
                awk -v batch=${{BATCH}} -v q=${{total_queries}} -v m=${{total_matches}} 'BEGIN{{print batch"\tavg_query\t"(m/(q/4))}}' >> {output.bench}

                # EVAL
                if [[ "{config[create_eval_file]}" == "True" ]]; then
                    {params.scripts_path}/query_eval.py {params.output_directory}/${{BATCH}}.eval {params.initial_batch} {params.incremental_batch} | sed "s/^/${{BATCH}}\t/" >> {output.bench}
                fi

                #rm -rf "{params.output_directory}/${{BATCH}}.query" "{params.output_directory}/${{BATCH}}.eval"
            fi
            
            # Delete index
            #rm -f "${{input}}".*
        done
        # Delete final index
        #rm -f "${{output}}".*
        """