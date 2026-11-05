rule pac_initial_build:
    input:
        initial_files = "input_data/initial_files.txt"
    output:
        initial_index = directory("pac/{args}/initial.index")
    benchmark:
        repeat("pac/{args}/initial_build.bench.tsv", config["repeat"])
    log:
        "pac/{args}/initial_build.log"
    threads:
        config["threads"]
    params:
        pac_path = config["tools"]["pac"]["path"],
        args = lambda wildcards: str2args(wildcards.args),
        fixed_args = config["tools"]["pac"]["fixed_args"]
    shell: 
        """      
        {params.pac_path} -f {input.initial_files} -d {output.initial_index} -c {threads} {params.fixed_args}  >> {log} 2>&1
        """

rule pac_batch_update:
    input:
        initial_index = "pac/{args}/initial.index",
        batches_dir = "input_data/batches"
    output:
        batch_update_bench = "pac/{args}/bench.tsv"
    benchmark:
        repeat("pac/{args}/batch_update.bench.tsv", config["repeat"])
    log:
        "pac/{args}/batch_update.log"
    threads:
        config["threads"]
    params:
        output_directory = "pac/{args}",
        pac_path = config["tools"]["pac"]["path"],
        scripts_path = os.path.join(workflow.basedir, "scripts"),
        ubt_cmd = lambda wildcards, output: get_usr_bin_time_cmd(output.batch_update_bench, "build"),
        ubt_cmd_qry = lambda wildcards, output: get_usr_bin_time_cmd(output.batch_update_bench, "query"),
        size_cmd = config["size_cmd"],
        fixed_args = config["tools"]["pac"]["fixed_args"],
        query_file = config["query_file"]
    shell: 
        """
        BATCHES=$(ls {input.batches_dir} | wc -l)
        for BATCH in $(seq 1 ${{BATCHES}}); do
            echo "${{BATCH}}" >> {log}

            if [[ BATCH -eq 1 ]]; then
                input={input.initial_index}
            else
                input="{params.output_directory}/$((BATCH - 1)).index"
            fi
            output="{params.output_directory}/${{BATCH}}.index"

            # BUILD
            {params.ubt_cmd} {params.pac_path} -f {input.batches_dir}/insert_batch_${{BATCH}}.txt \
                                               -l ${{input}} \
                                               -d ${{output}} \
                                               -c {threads} {params.fixed_args} >> {log} 2>&1

            # Workaround to make further updates (is that correct?)
            cp ${{input}}/MainIndex ${{output}}/MainIndex 

            # SIZE
            {params.size_cmd} "${{output}}" | sed "s/^/${{BATCH}}\tsize\t/" >> {output.batch_update_bench}

            # QUERY
            {params.ubt_cmd_qry} {params.pac_path} -l "${{output}}" -q "{params.query_file}" -c {threads} -o {params.output_directory}/${{BATCH}}.query.gz >> {log} 2>&1

        done
        """
