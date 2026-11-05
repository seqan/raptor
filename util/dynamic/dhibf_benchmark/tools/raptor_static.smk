rule raptor_static_batch_update:
    input:
        initial_files = "input_data/initial_files.txt",
        batches_dir = "input_data/batches"
    output:
        bench = "raptor_static/{args}/bench.tsv"
    benchmark:
        repeat("raptor_static/{args}/batch_update.bench.tsv", config["repeat"])
    log:
        "raptor_static/{args}/batch_update.log"
    threads:
        config["threads"]
    params:
        output_directory = "raptor_static/{args}",
        raptor_static_path = lambda wildcards: config["tools"]["raptor_static"]["path"],
        ubt_cmd = lambda wildcards, output: get_usr_bin_time_cmd(output.bench, "build"),
        ubt_cmd_qry = lambda wildcards, output: get_usr_bin_time_cmd(output.bench, "query"),
        size_cmd = config["size_cmd"],
        query_file = config["query_file"],
        query_args = lambda wildcards: config["tools"]["raptor_static"]["query_args"],
        args = lambda wildcards: str2args(wildcards.args),
        fixed_args_prepare = config["tools"]["raptor_static"]["fixed_args"]["prepare"],
        fixed_args_layout = config["tools"]["raptor_static"]["fixed_args"]["layout"],
        fixed_args_build = config["tools"]["raptor_static"]["fixed_args"]["build"],
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

            cat {input.initial_files} > {params.output_directory}/${{BATCH}}.list
            seq 1 ${{BATCH}} | xargs -I {{}} cat {input.batches_dir}/insert_batch_{{}}.txt >> {params.output_directory}/${{BATCH}}.list
     
            output="{params.output_directory}/${{BATCH}}.index"
            #{params.ubt_cmd} {params.raptor_static_path} prepare --threads {threads} --input {params.output_directory}/${{BATCH}}.list --output {params.output_directory}/${{BATCH}}_prepare {params.fixed_args_prepare} > {log} 2>&1
            #{params.ubt_cmd} {params.raptor_static_path} layout --threads {threads} --input {params.output_directory}/${{BATCH}}_prepare/minimiser.list --output {params.output_directory}/${{BATCH}}.layout {params.fixed_args_layout} {params.args} >> {log} 2>&1
            {params.ubt_cmd} {params.raptor_static_path} layout --threads {threads} --input {params.output_directory}/${{BATCH}}.list --output {params.output_directory}/${{BATCH}}.layout {params.fixed_args_layout} {params.args} > {log} 2>&1
            {params.ubt_cmd} {params.raptor_static_path} build --threads {threads} --input {params.output_directory}/${{BATCH}}.layout --output ${{output}} {params.fixed_args_build} >> {log} 2>&1
            
            # SIZE
            {params.size_cmd} "${{output}}" | sed "s/^/${{BATCH}}\tsize\t/" >> {output.bench}

            # QUERY
            if [[ -n "{params.query_file}" ]]; then
                {params.ubt_cmd_qry} {params.raptor_static_path} search {params.query_args} \
                                                                --threads {threads} \
                                                                --index "${{output}}" \
                                                                --query "{params.query_file}" \
                                                                --output {params.output_directory}/${{BATCH}}.query >> {log} 2>&1

                # AVG. QUERY HITS
                grep -v "^#" "{params.output_directory}/${{BATCH}}.query" | cut -f 2 | awk -v batch=${{BATCH}} -F"," '{{s+=NF}}END{{print batch"\tavg_query\t"s/NR}}' >> {output.bench}
                
                # EVAL
                if [[ "{config[create_eval_file]}" == "True" ]]; then
                    grep -v "^#" {params.output_directory}/${{BATCH}}.query | tr ',' '\t' | awk -F"\t" '{{print $1"\t"$2; for(i=3;i<=NF;i++) print $1"\t"$i;}}' > {params.output_directory}/${{BATCH}}.eval
                    {params.scripts_path}/query_eval.py {params.output_directory}/${{BATCH}}.eval {params.initial_batch} {params.incremental_batch} | sed "s/^/${{BATCH}}\t/" >> {output.bench}
                fi

                rm "{params.output_directory}/${{BATCH}}.query"
            fi

            # Delete files
            # rm -rf "{params.output_directory}/${{BATCH}}_prepare"
            rm -f "${{output}}" "{params.output_directory}/${{BATCH}}.list" "{params.output_directory}/${{BATCH}}.layout"
        done
        """
