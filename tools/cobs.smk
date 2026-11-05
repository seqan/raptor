rule cobs_batch_update:
    input:
        initial_files = "input_data/initial_files.txt",
        batches_dir = "input_data/batches"
    output:
        bench = "cobs/{args}/bench.tsv"
    benchmark:
        repeat("cobs/{args}/batch_update.bench.tsv", config["repeat"])
    log:
        "cobs/{args}/batch_update.log"
    threads:
        config["threads"]
    conda:
        ("../envs/cobs.yaml")
    params:
        output_directory = "cobs/{args}",
        cobs_path = lambda wildcards: config["tools"]["cobs"]["path"],
        ubt_cmd = lambda wildcards, output: get_usr_bin_time_cmd(output.bench, "build"),
        ubt_cmd_qry = lambda wildcards, output: get_usr_bin_time_cmd(output.bench, "query"),
        size_cmd = config["size_cmd"],
        query_file = config["query_file"],
        args = lambda wildcards: str2args(wildcards.args),
        fixed_args = config["tools"]["cobs"]["fixed_args"],
        query_args = config["tools"]["cobs"]["query_args"],
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
            
            output="{params.output_directory}/${{BATCH}}.index.cobs_compact"
            cat {input.initial_files} > {params.output_directory}/${{BATCH}}.list
            seq 1 ${{BATCH}} | xargs -I {{}} cat {input.batches_dir}/insert_batch_{{}}.txt >> {params.output_directory}/${{BATCH}}.list

            # BUILD FROM SCRATCH
            {params.ubt_cmd} {params.cobs_path} compact-construct --file-type list \
                        --threads {threads} \
                        {params.fixed_args} {params.args} \
                        {params.output_directory}/${{BATCH}}.list "${{output}}"  >> {log} 2>&1
            xargs -I {{}} rm -f {{}}.cobs_cache < {params.output_directory}/${{BATCH}}.list

            #SIZE
            {params.size_cmd} "${{output}}" | sed "s/^/${{BATCH}}\tsize\t/" >> {output.bench}

            # QUERY
            if [[ -n "{params.query_file}" ]]; then
                {params.ubt_cmd_qry} {params.cobs_path} query {params.query_args} \
                                                            --index "${{output}}" \
                                                            --file "{params.query_file}" \
                                                            --threads {threads} > "{params.output_directory}/${{BATCH}}.query" 2>> {log}
                
                # AVG. QUERY HITS
                grep "^\*" "{params.output_directory}/${{BATCH}}.query" | cut -f 2 | awk -v batch=${{BATCH}} '{{s+=$1}}END{{print batch"\tavg_query\t"s/NR}}' >> {output.bench}

                # EVAL
                if [[ "{config[create_eval_file]}" == "True" ]]; then
                    awk '{{if($1 ~ /^\*/){{printf "\\n"$1"\\t"}}else{{printf $1"\\t"}}}}END{{print ""}}' {params.output_directory}/${{BATCH}}.query | tail -n+2 | tr -d '*' > {params.output_directory}/${{BATCH}}.eval
                    {params.scripts_path}/query_eval.py {params.output_directory}/${{BATCH}}.eval {params.initial_batch} {params.incremental_batch} | sed "s/^/${{BATCH}}\t/" >> {output.bench}
                fi

                rm "{params.output_directory}/${{BATCH}}.query"
            fi
            
            # Delete index
            rm -f "${{output}}" {params.output_directory}/${{BATCH}}.list
        done
        """