rule kmindex_initial_build:
    input:
        initial_files = "input_data/initial_files.txt"
    output:
        index = directory("kmindex/{args}/index")
    benchmark:
        repeat("kmindex/{args}/initial_build.bench.tsv", config["repeat"])
    log:
        "kmindex/{args}/initial_build.log"
    threads:
        config["threads"]
    params:
        output_directory = "kmindex/{args}",
        kmindex_path = config["tools"]["kmindex"]["path"],
        args = lambda wildcards: str2args(wildcards.args),
        fixed_args = config["tools"]["kmindex"]["fixed_args"]
    shell: 
        """
        awk -F"/" '{{split($NF,f,"."); print f[1]" : "$0}}' "{input.initial_files}" > "{params.output_directory}/initial.list"
        {params.kmindex_path} build --fof "{params.output_directory}/initial.list" \
                                    --run-dir "{params.output_directory}/initial.index" \
                                    --index "{params.output_directory}/index" \
                                    --register-as initial \
                                    --threads {threads} \
                                    {params.fixed_args} {params.args} > {log} 2>&1
        """

rule kmindex_batch_update:
    input:
        index = "kmindex/{args}/index",
        batches_dir = "input_data/batches"
    output:
        batch_update_bench = "kmindex/{args}/bench.tsv"
    benchmark:
        repeat("kmindex/{args}/batch_update.bench.tsv", config["repeat"])
    log:
        "kmindex/{args}/batch_update.log"
    threads:
        config["threads"]
    params:
        output_directory = "kmindex/{args}",
        kmindex_path = lambda wildcards: config["tools"]["kmindex"]["path"],
        ubt_cmd = lambda wildcards, output: get_usr_bin_time_cmd(output.batch_update_bench, "build"),
        ubt_cmd_qry = lambda wildcards, output: get_usr_bin_time_cmd(output.batch_update_bench, "query"),
        size_cmd = config["size_cmd"],
        query_file = config["query_file"],
        args = lambda wildcards: str2args(wildcards.args),
        fixed_args = config["tools"]["kmindex"]["fixed_args"],
        query_args = config["tools"]["kmindex"]["query_args"]
    shell: 
        """
        BATCHES=$(ls {input.batches_dir} | wc -l)
        for BATCH in $(seq 1 ${{BATCHES}}); do
            echo "${{BATCH}}" >> {log}
            
            # BUILD INCREMENT
            awk -F"/" '{{split($NF,f,"."); print f[1]" : "$0}}' "{input.batches_dir}/insert_batch_${{BATCH}}.txt" > "{params.output_directory}/${{BATCH}}.list"
            {params.ubt_cmd} {params.kmindex_path} build --fof "{params.output_directory}/${{BATCH}}.list" \
                                    --run-dir "{params.output_directory}/${{BATCH}}.index" \
                                    --index "{input.index}" \
                                    --register-as "${{BATCH}}" \
                                    --threads {threads} \
                                    {params.fixed_args} {params.args} >> {log} 2>&1

            #SIZE
            {params.size_cmd} "{input.index}" | sed "s/^/${{BATCH}}\tsize\t/" >> {output.batch_update_bench}

            # QUERY
            {params.ubt_cmd_qry} {params.kmindex_path} query {params.query_args} \
                                                             --index "{input.index}" \
                                                             --fastx "{params.query_file}" \
                                                             --output "{params.output_directory}/${{BATCH}}.query" \
                                                             --format matrix \
                                                             --threads {threads} \
                                                             --aggregate >> {log} 2>&1

            # remove query subfolders since results are agg. in the tsv files (to save space)
            find {params.output_directory}/${{BATCH}}.query/ -maxdepth 1 -mindepth 1 -type d -exec rm -rf '{{}}' \;

        done
        rm -rf {params.output_directory}/*.index {params.output_directory}/index
        """