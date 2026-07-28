<!--
SPDX-FileCopyrightText: 2006-2026, Knut Reinert & Freie Universität Berlin
SPDX-FileCopyrightText: 2016-2026, Knut Reinert & MPI für molekulare Genetik
SPDX-License-Identifier: BSD-3-Clause
-->

# Simulated (even/RefSeq): k=20, w=20
# HumanRNA 1%: k=20, w=24

## cobs
### build
cobs compact-construct --file-type list --threads 24 --term-size ${k} --num-hashes 2 --false-positive-rate 0.05 ${input_fof} ${index}
### query
cobs query --load-complete --threshold 0.8 --index ${index} --file ${query} --threads 24 > ${results}

## fulgor
### build
fulgor build -d ${prefix} -l ${input_fof} -o ${index} -t 24 -k ${k} -m ${w} -g 512
### query
fulgor pseudoalign -i ${index}.fur -q ${query} -o ${results} -t 24

## metagraph
### build
metagraph build  -p 24 --kmer-length 20 -o ${index} < ${input_fof}
metagraph annotate -p 24 -i ${index} -o ${output_annotate} --anno-filename < ${input_fof}
## query
metagraph query --min-kmers-fraction-label 0.8 -p 24 -i ${index} -a ${index} ${query} > ${results}

## raptor_dynamic
### initial build
raptor layout --threads 24 --input ${input_fof} --output ${layout} --window ${w} --kmer ${k} --fpr 0.05 --hash 2 --disable-estimate-union --empty-bin-fraction 0.1
raptor build --threads 24 --input ${layout} --output ${index}
### batch updates
raptor update insert --threads 24 --index ${prev_output} --output ${index} --insert ${input_fof_batch}
### query
raptor search --error 2 --query_length 250 --threads 24 --index ${index} --query ${query} --output ${results}

## raptor_static
raptor layout --threads 24 --input ${input_fof} --output ${layout} --window ${w} --kmer ${k} --fpr 0.05 --hash 2 --disable-estimate-union --empty-bin-fraction 0
raptor build --threads 24 --input ${layout} --output ${index}
raptor search --error 2 --query_length 250 --threads 24 --index ${index} --query ${query} --output ${results}

# HumanRNA Full samples (k=20, w=24)

## raptor_dynamic
### initial build
raptor prepare --threads 24 --input ${input_fof} --output ${prepare} --window 24 --kmer 20 --use-filesize-dependent-cutoff
raptor layout --threads 24 --input ${prepare}/minimiser.list --output ${layout} --window 24 --kmer 20 --fpr 0.05 --hash 2 --disable-estimate-union --empty-bin-fraction 0.1
raptor build --threads 24 --input ${layout} --output ${index}
### batch updates
raptor prepare --threads 24 --input ${input_fof_batch} --output ${prepare_batch} --window 24 --kmer 20 --use-filesize-dependent-cutoff
raptor update insert --threads 24 --index ${prev_output} --output ${index} --insert ${prepare_batch}/minimiser.list
### query
raptor search --error 2 --query_length 250 --threads 24 --index ${index} --query ${query} --output ${results}
