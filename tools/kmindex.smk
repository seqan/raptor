rule kmindex_batch_update:
    input:
        initial_files = "input_data/initial_files.txt",
        batches_dir = "input_data/batches"
    output:
        bench = "kmindex/{args}/bench.tsv"
    benchmark:
        repeat("kmindex/{args}/batch_update.bench.tsv", config["repeat"])
    log:
        "kmindex/{args}/batch_update.log"
    conda:
        ("../envs/kmindex.yaml")
    threads:
        config["threads"]
    params:
        output_directory = "kmindex/{args}",
        kmindex_path = config["tools"]["kmindex"]["path"],
        ubt_cmd = lambda wildcards, output: get_usr_bin_time_cmd(output.bench, "build"),
        ubt_cmd_qry = lambda wildcards, output: get_usr_bin_time_cmd(output.bench, "query"),
        size_cmd = config["size_cmd"],
        query_file = config["query_file"],
        query_args = lambda wildcards: config["tools"]["kmindex"]["query_args"],
        args = lambda wildcards: str2args(wildcards.args),
        fixed_args = config["tools"]["kmindex"]["fixed_args"],
        initial_batch = config["initial_batch"] if "initial_batch" in config else 0,
        incremental_batch = config["incremental_batch"] if "incremental_batch" in config else 0,
        first_batch = config["first_batch"] if "first_batch" in config else 0,
        last_batch = config["last_batch"] if "last_batch" in config else "",
        scripts_path = os.path.join(workflow.basedir, "scripts"),
    shell: 
        """
        if [[ {params.first_batch} -eq 0 ]]; then
            BATCH=0
            awk -F"/" '{{split($NF,f,"."); print f[1]" : "$0}}' "{input.initial_files}" > "{params.output_directory}/${{BATCH}}.list"
            {params.ubt_cmd} {params.kmindex_path} build --fof "{params.output_directory}/${{BATCH}}.list" \
                                        --run-dir "@inplace" \
                                        --index "{params.output_directory}/index" \
                                        --register-as ${{BATCH}} \
                                        --threads {threads} \
                                        {params.fixed_args} {params.args} > {log} 2>&1
            
            {params.size_cmd} {params.output_directory}/index | sed "s/^/${{BATCH}}\tsize\t/" >> {output.bench}
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
            
            output="{params.output_directory}/${{BATCH}}"
            
            # BUILD INCREMENT
            awk -F"/" '{{split($NF,f,"."); print f[1]" : "$0}}' "{input.batches_dir}/insert_batch_${{BATCH}}.txt" > "{params.output_directory}/${{BATCH}}.list"
            {params.ubt_cmd} {params.kmindex_path} build --fof "{params.output_directory}/${{BATCH}}.list" \
                                    --run-dir "@inplace" \
                                    --index "{params.output_directory}/index" \
                                    --register-as "${{BATCH}}" \
                                    --threads {threads} \
                                    {params.fixed_args} {params.args} >> {log} 2>&1

            # MERGE
            #{params.ubt_cmd} {params.kmindex_path} merge --index "{params.output_directory}/index" \
            #                                             --new-name "{params.output_directory}/${{BATCH}}"
            #                                             --new-name merged_${{BATCH}} \
 
            #SIZE
            {params.size_cmd} "{params.output_directory}/index" | sed "s/^/${{BATCH}}\tsize\t/" >> {output.bench}

            # QUERY
            if [[ -n "{params.query_file}" ]]; then

                # QUERY
                {params.ubt_cmd_qry} {params.kmindex_path} query {params.query_args} \
                                                                --index "{params.output_directory}/index" \
                                                                --fastx "{params.query_file}" \
                                                                --output "{params.output_directory}/${{BATCH}}.query_tmp" \
                                                                --format jsonl \
                                                                --threads {threads} \
                                                                --aggregate >> {log} 2>&1
                                                                #--merge >> {log} 2>&1 #for query2


                # jq converts to .jsonl to tsv: query <tab> sample1,sample2,...,sampleN
                # awk concatenate results among jsons (repeated queries are printed in all files)
                jq -r '.query as $query | .samples | keys | join(",") as $keys | [$query, $keys] | join("\t")' {params.output_directory}/${{BATCH}}.query_tmp/*.jsonl | 
                    awk -F "\t" '{{if($2!=""){{t=$2","}}else{{t=""}};con[$1]=t con[$1]}}END{{for(c in con){{print c"\t"substr(con[c],0,length(con[c])-1)}}}}' > {params.output_directory}/${{BATCH}}.query_tmp/all.tsv

                # AVG. QUERY HITS
                cut -f 2 {params.output_directory}/${{BATCH}}.query_tmp/all.tsv | awk -v batch=${{BATCH}} -F"," '{{s+=NF}}END{{print batch"\tavg_query\t"s/NR}}' >> {output.bench}

                # EVAL
                if [[ "{config[create_eval_file]}" == "True" ]]; then
                    cat {params.output_directory}/${{BATCH}}.query_tmp/all.tsv | tr ',' '\t' > {params.output_directory}/${{BATCH}}.eval
                    {params.scripts_path}/query_eval.py {params.output_directory}/${{BATCH}}.eval {params.initial_batch} {params.incremental_batch} | sed "s/^/${{BATCH}}\t/" >> {output.bench}
                fi

                #rm -rf "{params.output_directory}/${{BATCH}}.query_tmp"
            fi

        done
        #rm -rf {params.output_directory}/index
        """
