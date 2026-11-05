#!/usr/bin/env bash

bname="$(basename ${1})"
outdir="${2}"
params_squeakr="${3}"

fq="${outdir}/${bname%.*}.fq"
sq=${fq}.squeakr

>&2 /home/mi/pirov/code/dhibf_benchmark/scripts/fasta_to_fastq --input ${1} --output ${fq}
>&2 /srv/public/seiler/dynamic_2025-12/squeakr/squeakr count ${params_squeakr} -o ${sq} ${fq}

echo ${sq}