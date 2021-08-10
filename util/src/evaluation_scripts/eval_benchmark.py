import numpy as np
import pandas as pd

import os
import re
import time

#Adjust parameters in generate_table_unpartitioned (bin and read count)
OUT_PATH = '<output directory of benchmark with bin_number directory, e.g. /dev/shm/username/1024'
#stores results in `OUTPATH/table.csv`

def process_output(path, bin_count, read_count, path_to_fp, path_to_fn):
    with open(os.path.join(OUT_PATH, path_to_fn), "w") as fn_f:
        with open(os.path.join(OUT_PATH, path_to_fp), "w") as fp_f:
            with open(path) as f:
                tp = 0
                fp = 0
                fn = 0
                for line in f:
                    try:
                        [x, y, thr, count] = line.strip().split('\t')
                    except ValueError:
                        try:
                            [x, y] = line.strip().split('\t')
                            thr = -1
                            count = -1
                        except ValueError:
                            fn += 1
                            continue
                    [read_id, bins] = [int(x), [int(e) for e in y[:-1].split(',') if e != '']]
                    true_id = (read_id % read_count) // (read_count // bin_count)
                    if true_id in bins:
                        tp += 1
                        if len(bins) != 1:
                            fp += len(bins) - 1
                            fp_f.write("False positive: read:{} bins:{} threshold:{} count:{}\n".format(read_id, bins, thr, count))
                    else:
                        fn += 1
                        fn_f.write("False negative: read:{} bins:{} threshold:{} count:{}\n".format(read_id, bins, thr, count))
                        fp += len(bins)
                        if len(bins) != 0:
                            fp_f.write("False positive: read:{} bins:{} threshold:{} count:{}\n".format(read_id, bins, thr, count))
                return [tp,fp,fn]

def process_time_log(path):
    time = []
    mem = []
    with open(path) as f:
        for i, line in enumerate(f):
            if i == 4:
                time = line.strip().split(' ')[-1].split(':')
                if len(time) == 3:
                    time = '{}:{:02d}'.format(int(time[0]) * 60 + int(time[1]), int(time[2]))
                elif len(time) == 2:
                    if int(time[0]) == 0:
                        time = '{}:{}'.format(time[0], time[1])
                    else:
                        time = '{}:{:02d}'.format(time[0], round(float(time[1])))
            if i == 9:
                mem = round(int(line.strip().split(' ')[-1]) / 1024)
    return [time, mem]

def process_internal_time(path):
    with open(path) as f:
        line = f.readline()
        line = f.readline()
        [ibf_time, read_time, compute_time] = line.strip().split('\t')
        return [ibf_time, read_time, compute_time]

def get_param_list(path = OUT_PATH):
    result = []
    for entry in os.listdir(path):
        match = re.match('(\d+)_(\d+)_(\d+\w).out$', entry)
        if match is not None:
            result.append( (match.group(1), match.group(2), match.group(3)) )
    return result

# (how many bins: int), (how many reads: int)
def generate_table_unpartitioned(bin_count=1024, read_count=1048576):
    data = []
    params = get_param_list()
    format_string = 'Processing file {{:>{}}} of {}...'.format(len(str(len(params))), len(params))
    for i, (window_size, kmer_size, ibf_size) in enumerate(params):
        print(format_string.format(i + 1), end='', flush=True)
        path_to_build_log = os.path.join(OUT_PATH, '{}_{}_{}_build.log'.format(window_size, kmer_size, ibf_size))
        path_to_query_log = os.path.join(OUT_PATH, '{}_{}_{}_query.log'.format(window_size, kmer_size, ibf_size))
        path_to_output = os.path.join(OUT_PATH, '{}_{}_{}.out'.format(window_size, kmer_size, ibf_size))
        path_to_internal_time = os.path.join(OUT_PATH, '{}_{}_{}.out.time'.format(window_size, kmer_size, ibf_size))
        path_to_fp = os.path.join(OUT_PATH, 'fp_{}_{}_{}.txt'.format(window_size, kmer_size, ibf_size))
        path_to_fn = os.path.join(OUT_PATH, 'fn_{}_{}_{}.txt'.format(window_size, kmer_size, ibf_size))

        [build_time, build_ram] = process_time_log(path_to_build_log)
        [query_time, query_ram] = process_time_log(path_to_query_log)
        [query_ibf_time, query_reads_time, query_compute_time] = process_internal_time(path_to_internal_time)
        [tp, fp, fn] = process_output(path_to_output, bin_count, read_count, path_to_fp, path_to_fn)
        data.append([build_time, build_ram, query_time, query_ibf_time, query_reads_time, query_compute_time, query_ram, fp, fn])
        print('Done', flush=True)
    df = pd.DataFrame(data, columns = [('Construct', 'Time'),  ('Construct', 'RAM'), ('Search','Overall'),  ('Search','IBF I/O'), ('Search','Reads I/O'),  ('Search','Compute'), ('Search','RAM'),  ('Search','FP'), ('Search','FN')])
    df.index = [str(x).replace(',', ';').replace("'", '').replace('m', '') for x in params]
    df.columns = pd.MultiIndex.from_tuples(df.columns, names=['','w k size'])
    return df

path_to_csv = os.path.join(OUT_PATH, 'table.csv')
df = generate_table_unpartitioned()
df.to_csv(path_to_csv)
