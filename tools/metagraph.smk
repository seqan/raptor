rule metagraph_batch_update:
    input:
        initial_files = "input_data/initial_files.txt",
        batches_dir = "input_data/batches"
    output:
        bench = "metagraph/{args}/bench.tsv"
    benchmark:
        repeat("metagraph/{args}/batch_update.bench.tsv", config["repeat"])
    log:
        "metagraph/{args}/batch_update.log"
    threads:
        config["threads"]
    conda:
        ("../envs/metagraph.yaml")
    params:
        output_directory = "metagraph/{args}",
        metagraph_path = config["tools"]["metagraph"]["path"],
        ubt_cmd = lambda wildcards, output: get_usr_bin_time_cmd(output.bench, "build"),
        ubt_cmd_qry = lambda wildcards, output: get_usr_bin_time_cmd(output.bench, "query"),
        size_cmd = config["size_cmd"],
        fixed_args = config["tools"]["metagraph"]["fixed_args"],
        query_args = config["tools"]["metagraph"]["query_args"],
        query_file = config["query_file"],
        initial_batch = config["initial_batch"] if "initial_batch" in config else 0,
        incremental_batch = config["incremental_batch"] if "incremental_batch" in config else 0,
        first_batch = config["first_batch"] if "first_batch" in config else 1,
        last_batch = config["last_batch"] if "last_batch" in config else "",
        scripts_path = os.path.join(workflow.basedir, "scripts"),
    shell: 
        """
        if [[ -z "{params.last_batch}" ]]; then
            last_batch=$(ls {input.batches_dir}/insert_batch_* | wc -l)
        else
            last_batch={params.last_batch}
        fi

        for BATCH in $(seq {params.first_batch} ${{last_batch}}); do
            echo "${{BATCH}}" >> {log}

            output="{params.output_directory}/${{BATCH}}.index.dbg"
            cat {input.initial_files} > {params.output_directory}/${{BATCH}}.list
            seq 1 ${{BATCH}} | xargs -I {{}} cat {input.batches_dir}/insert_batch_{{}}.txt >> {params.output_directory}/${{BATCH}}.list
           
            # BUILD
            {params.ubt_cmd} {params.metagraph_path} build \
                                                    -p {threads} {params.fixed_args} \
                                                    -o ${{output}} \
                                                    < {params.output_directory}/${{BATCH}}.list >> {log} 2>&1
            # ANNOTATE
            {params.ubt_cmd} {params.metagraph_path} annotate -p {threads} -i "${{output}}" -o {params.output_directory}/${{BATCH}}.column.annodbg --anno-filename < {params.output_directory}/${{BATCH}}.list >> {log} 2>&1

            # SIZE
            {params.size_cmd} "${{output}}" "{params.output_directory}/${{BATCH}}.column.annodbg" | sed "s/^/${{BATCH}}\tsize\t/" >> {output.bench}

            # QUERY
            if [[ -n "{params.query_file}" ]]; then
                {params.ubt_cmd_qry} {params.metagraph_path} query {params.query_args} \
                                                                -p {threads} \
                                                                -i "${{output}}" \
                                                                -a {params.output_directory}/${{BATCH}}.column.annodbg \
                                                                {params.query_file} > {params.output_directory}/${{BATCH}}.query 2>> {log}
                # AVG. QUERY HITS
                cut -f 3 "{params.output_directory}/${{BATCH}}.query" | awk -v batch=${{BATCH}} -F":" '{{s+=NF}}END{{print batch"\tavg_query\t"s/NR}}' >> {output.bench}
                
                # EVAL
                if [[ "{config[create_eval_file]}" == "True" ]]; then
                    cut -f 2,- "{params.output_directory}/${{BATCH}}.query" | tr ':' '\t' | awk -F"\t" '{{print $1"\t"$2; for(i=3;i<=NF;i++) print $1"\t"$i;}}' | sed 's|/.*/||' > {params.output_directory}/${{BATCH}}.eval
                    {params.scripts_path}/query_eval.py {params.output_directory}/${{BATCH}}.eval {params.initial_batch} {params.incremental_batch} | sed "s/^/${{BATCH}}\t/" >> {output.bench}
                fi

                rm "{params.output_directory}/${{BATCH}}.query"
            fi
            
            # Delete tmp files
            rm -f "${{output}}" {params.output_directory}/${{BATCH}}.column.annodbg {params.output_directory}/${{BATCH}}.list
        done
        
        # Delete final index
        rm -f "${{output}}"
        """
