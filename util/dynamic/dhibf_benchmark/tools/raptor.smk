rule raptor_batch_update:
    input:
        initial_files = "input_data/initial_files.txt",
        batches_dir = "input_data/batches"
    output:
        bench = "raptor/{args}/bench.tsv"
    benchmark:
        repeat("raptor/{args}/batch_update.bench.tsv", config["repeat"])
    log:
        "raptor/{args}/batch_update.log"
    threads:
        config["threads"]
    params:
        output_directory = "raptor/{args}",
        raptor_path = lambda wildcards: config["tools"]["raptor"]["path"],
        ubt_cmd = lambda wildcards, output: get_usr_bin_time_cmd(output.bench, "build"),
        ubt_cmd_qry = lambda wildcards, output: get_usr_bin_time_cmd(output.bench, "query"),
        size_cmd = config["size_cmd"],
        query_file = config["query_file"],
        query_args = lambda wildcards: config["tools"]["raptor"]["query_args"],
        args = lambda wildcards: str2args(wildcards.args),
        fixed_args_prepare = config["tools"]["raptor"]["fixed_args"]["prepare"],
        fixed_args_layout = config["tools"]["raptor"]["fixed_args"]["layout"],
        fixed_args_build = config["tools"]["raptor"]["fixed_args"]["build"],
        initial_batch = config["initial_batch"] if "initial_batch" in config else 0,
        incremental_batch = config["incremental_batch"] if "incremental_batch" in config else 0,
        first_batch = config["first_batch"] if "first_batch" in config else 0,
        last_batch = config["last_batch"] if "last_batch" in config else "",
        scripts_path = os.path.join(workflow.basedir, "scripts"),
    shell: 
        """
        if [[ {params.first_batch} -eq 0 ]]; then
            BATCH=0
            # prepare
            #{params.ubt_cmd} {params.raptor_path} prepare --threads {threads} --input {input.initial_files} --output {params.output_directory}/${{BATCH}}_prepare {params.fixed_args_prepare} > {log} 2>&1
            #{params.ubt_cmd} {params.raptor_path} layout --threads {threads} --input {params.output_directory}/${{BATCH}}_prepare/minimiser.list --output {params.output_directory}/${{BATCH}}.layout {params.fixed_args_layout} {params.args} >> {log} 2>&1
            # prepare

            # no prepare
            {params.ubt_cmd} {params.raptor_path} layout --threads {threads} --input {input.initial_files} --output {params.output_directory}/${{BATCH}}.layout {params.fixed_args_layout} {params.args} >> {log} 2>&1
            # no prepare
            
            {params.ubt_cmd} {params.raptor_path} build --threads {threads} --input {params.output_directory}/${{BATCH}}.layout --output {params.output_directory}/${{BATCH}}.index {params.fixed_args_build} >> {log} 2>&1
            {params.size_cmd} {params.output_directory}/${{BATCH}}.index | sed "s/^/${{BATCH}}\tsize\t/" >> {output.bench}
            echo "${{BATCH}}\tobs\tInitial Build" >> {output.bench}
            rm {params.output_directory}/${{BATCH}}.layout
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
            
            # UPDATE - delete
            if [[ -f "{input.batches_dir}/delete_batch_${{BATCH}}.txt" ]]; then
                {params.ubt_cmd_qry} {params.raptor_path} inspect --index "${{input}}" \
                                                                  --format json --user-bins |
                                                          jq -r '.user_bins.paths | to_entries[] | "\(.value[0])\t\(.key)"' > {params.output_directory}/$((BATCH - 1)).user_bins.tsv 2>> {log}

                # Make list of ub to delete "--delete 0 --delete 1"                
                to_delete=$(join <(sort "{input.batches_dir}/delete_batch_${{BATCH}}.txt") <(sort {params.output_directory}/$((BATCH - 1)).user_bins.tsv) -o "2.2" | sed 's/^/--delete /' | tr '\n' ' ')
                {params.ubt_cmd} {params.raptor_path} update delete --threads {threads} \
                                                                    --index "${{input}}" \
                                                                    --output "${{output}}_deleted" \
                                                                    ${{to_delete}} >> {log} 2>&1
                rm -f "${{input}}" {params.output_directory}/$((BATCH - 1)).user_bins.tsv
                input="${{output}}_deleted"
            fi

            # UPDATE - insert
            # prepare
            #{params.ubt_cmd} {params.raptor_path} prepare --threads {threads} --input "{input.batches_dir}/insert_batch_${{BATCH}}.txt" --output {params.output_directory}/${{BATCH}}_prepare {params.fixed_args_prepare} >> {log} 2>&1
            #message=$({params.ubt_cmd} {params.raptor_path} update insert --threads {threads} --index "${{input}}" --output "${{output}}" --insert {params.output_directory}/${{BATCH}}_prepare/minimiser.list --stats 2>&1)
            # prepare

            # no prepare
            message=$({params.ubt_cmd} {params.raptor_path} update insert --threads {threads} --index "${{input}}" --output "${{output}}" --insert "{input.batches_dir}/insert_batch_${{BATCH}}.txt" --stats 2>&1)
            # no prepare

            echo "${{message}}" >> {log}

            has_partial_rebuild=false
            has_full_rebuild=false
            if grep "Partial Rebuilds: [1-9]" <<< "${{message}}"; then
                has_partial_rebuild=true
            fi
            if grep "Full Rebuilds: [1-9]" <<< "${{message}}"; then
                has_full_rebuild=true
            fi
            if ${{has_partial_rebuild}} && ${{has_full_rebuild}}; then
                message="Mixed Rebuild"
            elif ${{has_partial_rebuild}}; then
                message="Partial Rebuild"
            elif ${{has_full_rebuild}}; then
                message="Full Rebuild"
            else
                message="Inplace"
            fi

            # OBS
            echo "${{BATCH}}\tobs\t${{message}}" >> {output.bench}

            # SIZE
            {params.size_cmd} "${{output}}" | sed "s/^/${{BATCH}}\tsize\t/" >> {output.bench}

            # QUERY
            if [[ -n "{params.query_file}" ]]; then
                {params.ubt_cmd_qry} {params.raptor_path} search {params.query_args} \
                                                                --threads {threads} \
                                                                --index "${{output}}" \
                                                                --query "{params.query_file}" \
                                                                --output {params.output_directory}/${{BATCH}}.query >> {log} 2>&1

                # AVG. QUERY HITS
                grep -v "^#" "{params.output_directory}/${{BATCH}}.query" | cut -f 2 | awk -v batch=${{BATCH}} -F"," '{{s+=NF}}END{{print batch"\tavg_query\t"s/NR}}' >> {output.bench}
                
                # EVAL
                if [[ "{config[create_eval_file]}" == "True" ]]; then
                    grep -v "^#" {params.output_directory}/${{BATCH}}.query | tr ',' '\t' | awk -F"\t" '{{print $1"\t"$2; for(i=3;i<=NF;i++) print $1"\t"$i;}}' > {params.output_directory}/${{BATCH}}.eval_tmp
                    
                    # Re-match user bin id with linked "file" id
                    {params.raptor_path} inspect --index "${{output}}" \
                                                 --format json --user-bins |
                                                 jq -r '.user_bins.paths | to_entries[] | "\(.value[0])\t\(.key)"' > {params.output_directory}/${{BATCH}}.user_bins.tsv 2>> {log}
                    join -1 2 -2 2 <(sort -t$'\t' -k 2,2 {params.output_directory}/${{BATCH}}.eval_tmp) <(sort -t$'\t' -k 2,2 {params.output_directory}/${{BATCH}}.user_bins.tsv) -t$'\t' -o "1.1,2.1" -a 1 | sed -e 's|/.*/bin_||' -e 's|.fa$||' > {params.output_directory}/${{BATCH}}.eval
                    {params.scripts_path}/query_eval.py {params.output_directory}/${{BATCH}}.eval {params.initial_batch} {params.incremental_batch} | sed "s/^/${{BATCH}}\t/" >> {output.bench}
                fi

                rm -f "{params.output_directory}/${{BATCH}}.query" "{params.output_directory}/${{BATCH}}.eval" "{params.output_directory}/${{BATCH}}.eval_tmp"
            fi

            # Delete index
            rm -f "${{input}}"
        done
        
        # Delete final index
        rm -f "${{output}}"
        
        # Delete minimisers (with prepare only)
        # prepare
        rm -rf {params.output_directory}/*_prepare
        # prepare
        """
