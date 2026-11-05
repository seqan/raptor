rule mantis_initial_build:
    input:
        initial_files = "input_data/initial_files.txt"
    output:
        initial_index = directory("mantis/{args}/initial.index")
    benchmark:
        repeat("mantis/{args}/initial_build.bench.tsv", config["repeat"])
    log:
        "mantis/{args}/initial_build.log"
    threads:
        config["threads"]
    params:
        output_directory = "mantis/{args}",
        mantis_path = config["tools"]["mantis"]["path"],
        scripts_path = os.path.join(workflow.basedir, "scripts"),
        args = lambda wildcards: str2args(wildcards.args),
        fixed_args_squeakr = config["tools"]["mantis"]["fixed_args"]["squeakr"],
        fixed_args_build = config["tools"]["mantis"]["fixed_args"]["build"]
    shell: 
        """      
        rm -rf {params.output_directory}/squeakr
        mkdir -p {params.output_directory}/squeakr

        # Run fasta_to_fastq and squeakr in parallel
        xargs -P {threads} -I {{}} {params.scripts_path}/squeakr.sh {{}} {params.output_directory}/squeakr '{params.fixed_args_squeakr}' < {input.initial_files} > {params.output_directory}/squeakr/initial_list.tmp 2>> {log}

        {params.mantis_path} build \
                            -t {threads} \
                            -i {params.output_directory}/squeakr/initial_list.tmp \
                            -o {output.initial_index} {params.fixed_args_build} >> {log} 2>&1

        {params.mantis_path} mst \
                            -p {output.initial_index} \
                            -t {threads} \
                            -d  >> {log} 2>&1
        """

rule mantis_batch_update:
    input:
        initial_index = directory("mantis/{args}/initial.index"),
        batches_dir = "input_data/batches"
    output:
        batch_update_bench = "mantis/{args}/bench.tsv"
    benchmark:
        repeat("mantis/{args}/batch_update.bench.tsv", config["repeat"])
    log:
        "mantis/{args}/batch_update.log"
    threads:
        config["threads"]
    params:
        output_directory = "mantis/{args}",
        mantis_path = config["tools"]["mantis"]["path"],
        scripts_path = os.path.join(workflow.basedir, "scripts"),
        ubt_cmd = lambda wildcards, output: get_usr_bin_time_cmd(output.batch_update_bench, "build"),
        ubt_cmd_qry = lambda wildcards, output: get_usr_bin_time_cmd(output.batch_update_bench, "query"),
        size_cmd = config["size_cmd"],
        fixed_args_squeakr = config["tools"]["mantis"]["fixed_args"]["squeakr"],
        fixed_args_build = config["tools"]["mantis"]["fixed_args"]["build"],
        fixed_args_query = config["tools"]["mantis"]["fixed_args"]["query"],
        query_file = config["query_file"]
    shell: 
        """
        BATCHES=$(ls {input.batches_dir} | wc -l)
        for BATCH in $(seq 1 ${{BATCHES}}); do
            echo "${{BATCH}}" >> {log}

            # Run fasta_to_fastq and squeakr in parallel
            {params.ubt_cmd} xargs -P {threads} -I {{}} {params.scripts_path}/squeakr.sh {{}} {params.output_directory}/squeakr '{params.fixed_args_squeakr}' < {input.batches_dir}/insert_batch_${{BATCH}}.txt > {params.output_directory}/squeakr/${{BATCH}}_list.tmp 2>> {log}

            input="{params.output_directory}/$((BATCH - 1)).index"
            output="{params.output_directory}/${{BATCH}}.index"
            if [[ BATCH -eq 1 ]]; then
                input={input.initial_index}
            fi
            
            tmp="{params.output_directory}/${{BATCH}}.tmp.index"
            {params.ubt_cmd} {params.mantis_path} build \
                            -t {threads} \
                            -i {params.output_directory}/squeakr/${{BATCH}}_list.tmp \
                            -o ${{tmp}} {params.fixed_args_build} >> {log} 2>&1

            {params.ubt_cmd} {params.mantis_path} mst -d \
                            -p ${{tmp}} \
                            -t {threads} >> {log} 2>&1

            {params.ubt_cmd} {params.mantis_path} merge -rm \
                                -i1 ${{input}} \
                                -i2 ${{tmp}} \
                                -t {threads} \
                                -o ${{output}} >> {log} 2>&1

            # SIZE
            {params.size_cmd} "${{output}}" | sed "s/^/${{BATCH}}\tsize\t/" >> {output.batch_update_bench}

            # QUERY
            {params.ubt_cmd_qry} {params.mantis_path} query \
                                        -p "${{output}}" \
                                        -o {params.output_directory}/${{BATCH}}.query \
                                        {params.fixed_args_query} \
                                        "{params.query_file}" >> {log} 2>&1
            rm -rf "${{tmp}}"
        done
        """
