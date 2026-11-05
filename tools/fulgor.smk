rule fulgor_batch_update:
    input:
        initial_files = "input_data/initial_files.txt",
        batches_dir = "input_data/batches"
    output:
        bench = "fulgor/{args}/bench.tsv"
    benchmark:
        repeat("fulgor/{args}/batch_update.bench.tsv", config["repeat"])
    log:
        "fulgor/{args}/batch_update.log"
    threads:
        config["threads"]
    conda:
        ("../envs/fulgor.yaml")
    params:
        output_directory = "fulgor/{args}",
        fulgor_path = lambda wildcards: config["tools"]["fulgor"]["path"],
        ubt_cmd = lambda wildcards, output: get_usr_bin_time_cmd(output.bench, "build"),
        ubt_cmd_qry = lambda wildcards, output: get_usr_bin_time_cmd(output.bench, "query"),
        size_cmd = config["size_cmd"],
        query_file = config["query_file"],
        args = lambda wildcards: str2args(wildcards.args),
        fixed_args = config["tools"]["fulgor"]["fixed_args"],
        query_args = config["tools"]["fulgor"]["query_args"],
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
            
            output="{params.output_directory}/${{BATCH}}.index"
            cat {input.initial_files} > {params.output_directory}/${{BATCH}}.list
            seq 1 ${{BATCH}} | xargs -I {{}} cat {input.batches_dir}/insert_batch_{{}}.txt >> {params.output_directory}/${{BATCH}}.list
           
            # BUILD FROM SCRATCH
            {params.ubt_cmd} {params.fulgor_path} build -d "{params.output_directory}" \
                                                        -l "{params.output_directory}/${{BATCH}}.list" \
                                                        -o "${{output}}" \
                                                        -t {threads} \
                                                        {params.fixed_args} {params.args} >> {log} 2>&1

            #SIZE
            {params.size_cmd} "${{output}}.fur" | sed "s/^/${{BATCH}}\tsize\t/" >> {output.bench}

            # QUERY
            if [[ -n "{params.query_file}" ]]; then
                {params.ubt_cmd_qry} {params.fulgor_path} pseudoalign {params.query_args} \
                                                                    -i "${{output}}.fur" \
                                                                    -q "{params.query_file}" \
                                                                    -o {params.output_directory}/${{BATCH}}.query \
                                                                    -t {threads} >> {log} 2>&1

                # AVG. QUERY HITS
                cut -f 2 "{params.output_directory}/${{BATCH}}.query" | awk -v batch=${{BATCH}} '{{s+=$1}}END{{print batch"\tavg_query\t"s/NR}}' >> {output.bench}

                # EVAL
                if [[ "{config[create_eval_file]}" == "True" ]]; then
                    awk -F"\t" '{{print $1"\t"$3; for(i=4;i<=NF;i++) print $1"\t"$i;}}' {params.output_directory}/${{BATCH}}.query > {params.output_directory}/${{BATCH}}.eval
                    {params.scripts_path}/query_eval.py {params.output_directory}/${{BATCH}}.eval {params.initial_batch} {params.incremental_batch} | sed "s/^/${{BATCH}}\t/" >> {output.bench}
                fi

                rm "{params.output_directory}/${{BATCH}}.query"
            fi
            
            # Delete index
            rm -f "${{output}}.fur" "{params.output_directory}/${{BATCH}}.list"
        done
        """
