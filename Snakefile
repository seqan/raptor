workdir: config["workdir"]
include: "utils.py"

if "raptor" in config["tools"]:
    include: "tools/raptor.smk"
if "raptor_static" in config["tools"]:
    include: "tools/raptor_static.smk"    
if "metagraph" in config["tools"]:    
    include: "tools/metagraph.smk"
if "fulgor" in config["tools"]:    
    include: "tools/fulgor.smk"
if "cobs" in config["tools"]:    
    include: "tools/cobs.smk"
if "kmindex" in config["tools"]:    
    include: "tools/kmindex.smk"
if "bifrost" in config["tools"]:    
    include: "tools/bifrost.smk"

def input_build():
    out = []
    for tool in config["tools"]:
        # define path to output
        path = tool + "/"
        # build product of all arguments (single or range)
        for args in args_product(config["tools"][tool]["args"] if "args" in config["tools"][tool] else None):
            out.append(path + join_args(sorted(args)))
    import pprint
    pprint.pprint(out)
    return out

rule all:
    input:
        summaries = expand("{i}/{ext}", i=input_build(), ext=["summary.tsv"]),
    output:
        report = "report.tsv",
        plot = "report.html"
    params:
        report_unit = config["report_unit"],
        scripts_path = os.path.join(workflow.basedir, "scripts"),
        wd = config["workdir"]
    conda:
        ("envs/plot.yaml")
    shell:
        """
        head -n 1 {input.summaries[0]} > {output.report}
        awk 'FNR>1' {input.summaries} >> {output.report}
        {params.scripts_path}/plot.py {output.report} {params.report_unit} {workflow.config_settings.configfiles}
        """

rule summarize:
    input:
        bench = "{tool}/{args}/bench.tsv",
    output:
        summary = "{tool}/{args}/summary.tsv"
    params:
        report_unit = config["report_unit"],
        scripts_path = os.path.join(workflow.basedir, "scripts"),
    shell: 
        """
        python3 {params.scripts_path}/eval.py \
                --input {input.bench} \
                --output {output.summary} \
                --unit {params.report_unit} \
                --label "{wildcards.tool} {wildcards.args}"
        """

rule pre_process_input_files:
    output:
        batches_dir = directory("input_data/batches"),
        initial_files = "input_data/initial_files.txt"
    log:
        "input_data/pre_process_input_files.log"
    params:
        initial_files = config["initial_files"] if "initial_files" in config else "",
        update_files_dir = config["update_files_dir"] if "update_files_dir" in config else "",
        input_files_list = config["input_files_list"] if "input_files_list" in config else "",
        initial_batch = config["initial_batch"] if "initial_batch" in config else 0,
        incremental_batch = config["incremental_batch"] if "incremental_batch" in config else 0,
        delete_batch = config["delete_batch"] if "delete_batch" in config else 0,
        max_batches = config["max_batches"] if "max_batches" in config else 0
    shell:
        """
        mkdir -p "input_data/batches"
        if [ ! -z {params.initial_files} ]; then
            cp {params.initial_files} {output.initial_files}
            cp {params.update_files_dir}*.txt {output.batches_dir}/
            for filename in {output.batches_dir}/*.txt; do
                mv $filename "{output.batches_dir}/insert_batch_${{filename##*/}}"
            done
        else
            LINES=$(wc -l < "{params.input_files_list}")
            if [[ "${{LINES}}" -lt "{params.initial_batch}" ]]; then
                echo "Not enough lines in {params.input_files_list} to create initial index. Found ${{LINES}}, required {params.initial_batch}." >> {log}
                exit 1
            fi
            head -n "{params.initial_batch}" "{params.input_files_list}" > "input_data/initial_files.txt"

            if [[ {params.delete_batch} -gt 0 ]]; then
                tail -n "{params.delete_batch}" "input_data/initial_files.txt" > "input_data/batches/delete_batch_1.txt"
            fi

            tail -n +$(({params.initial_batch} + 1)) "{params.input_files_list}" > "input_data/insert_files.txt"
            awk -vc=1 -vl={params.incremental_batch} -vm={params.max_batches} 'NR>1 && (NR-1)%l==0{{++c}}{{if(m == 0 || c <= m) print $0 > "input_data/batches/insert_batch_"c".txt"}}' input_data/insert_files.txt
            rm -rf input_data/insert_files.txt

            if [[ {params.delete_batch} -gt 0 ]]; then
                for i in $(seq 2 $(ls "input_data/batches/" | wc -l)); do
                    tail -n "{params.delete_batch}" "input_data/batches/insert_batch_$((i-1)).txt" > "input_data/batches/delete_batch_$i.txt"
                done
            fi
        fi
        """