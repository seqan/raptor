#!/usr/bin/env python3

import csv
import re
import argparse

def format_memory(build_memory_kb, unit):
    memory_mb = build_memory_kb / 1024
    if unit == 'GiB':
        return f"{memory_mb / 1024:.2f}"
    else:
        return f"{memory_mb:.2f}"

def format_size(size_bytes, unit):
    size_mb = size_bytes / (1024 * 1024)
    if unit == 'GiB':
        return f"{size_mb / 1024:.2f}"
    else:
        return f"{size_mb:.2f}"

def process_file(input_file, output_file, unit, label):
    summ = {}
    with open(input_file, 'r') as infile:
        for line in infile:
            fields = line.rstrip().split("\t")
            batch = fields[0]
            type = fields[1]

            if batch not in summ:
                summ[batch] = {"label": label, "obs": "", "index_size_bytes": 0, "build_time_seconds": 0, "build_memory_kb": 0, "query_time_seconds": 0, "query_memory_kb": 0, "build+query_time_seconds": 0, "build+query_memory_kb": 0, "TP": 0, "FP": 0, "FN": 0, "TN": 0, "avg_query": 0}
            
            if type == "obs":
                # concatenate observations with ;
                summ[batch]["obs"] += (";" + fields[2]) if summ[batch]["obs"] else fields[2]
            elif type == "size":
                # Sum sizes
                summ[batch]["index_size_bytes"]+=int(fields[2])
            elif type == "build":
                # Sum time
                summ[batch]["build_time_seconds"]+=float(fields[3])
                summ[batch]["build+query_time_seconds"]+=float(fields[3])
                # Get peak memory build
                if int(fields[4]) > summ[batch]["build_memory_kb"]:
                    summ[batch]["build_memory_kb"] = int(fields[4])
                # Get peak memory build+query
                if int(fields[4]) > summ[batch]["build+query_memory_kb"]:
                    summ[batch]["build+query_memory_kb"] = int(fields[4])
            elif type == "query":
                summ[batch]["query_time_seconds"]+=float(fields[3])
                summ[batch]["build+query_time_seconds"]+=float(fields[3])
                summ[batch]["query_memory_kb"] = int(fields[4])
                # Get peak memory build+query
                if int(fields[4]) > summ[batch]["build+query_memory_kb"]:
                    summ[batch]["build+query_memory_kb"] = int(fields[4])
            elif type in ["TP", "FP", "FN", "TN"]:
                summ[batch][type]+=int(fields[2])
            elif type == "avg_query":
                summ[batch][type]+=float(fields[2])

    with open(output_file, 'w') as outfile:
        outfile.write("\t".join(["Label", 
                        "Run Number",
                        "Build time (sec)",
                        f"Build peak memory ({unit})",
                        f"Index size ({unit})",
                        "Query time (sec)",
                        f"Query peak memory ({unit})",
                        "Build+Query time (sec)",
                        f"Build+Query peak memory ({unit})",
                        "Obs",
                        "TP", "FP", "FN", "TN", "avg_query"]))
        outfile.write("\n")
        for batch, v in summ.items():
            outfile.write("\t".join([v["label"],
                            batch,
                            str(v["build_time_seconds"]),
                            format_memory(v["build_memory_kb"], unit),
                            format_size(v["index_size_bytes"], unit),
                            str(v["query_time_seconds"]),
                            format_memory(v["query_memory_kb"], unit),
                            str(v["build+query_time_seconds"]),
                            format_memory(v["build+query_memory_kb"], unit),
                            v["obs"].rstrip(),
                            str(v["TP"]), str(v["FP"]), str(v["FN"]), str(v["TN"]), str(v["avg_query"])
                            ]))
            outfile.write("\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Summarize update benchmark metrics.')
    parser.add_argument('--input', help='Path to the input file')
    parser.add_argument('--output', help='Path to the output file')
    parser.add_argument('--label', help='Label of bench file (default: output)', default="output")
    parser.add_argument('--unit', choices=['MiB', 'GiB'], default='GiB', help='Unit to format memory and size (default: GiB)')
    args = parser.parse_args()

    process_file(args.input, args.output, args.unit, args.label)
