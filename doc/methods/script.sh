#!/usr/bin/env bash
set -Exeuo pipefail

if [[ ! -f "example_data.tar.gz" ]]; then
## [01_introduction_snippet_1]
wget https://ftp.imp.fu-berlin.de/pub/seiler/raptor/example_data.tar.gz
tar xfz example_data.tar.gz
## [01_introduction_snippet_1]
fi

## [01_introduction_snippet_2]
seq -f "example_data/64/bins/bin_%02g.fasta" 0 1 63 > all_bin_paths.txt
raptor build --input all_bin_paths.txt --kmer 19 --window 23 --fpr 0.05 --threads 2 --output raptor.index
## [01_introduction_snippet_2]

## [01_introduction_snippet_3]
raptor search --error 2 --index raptor.index --query example_data/64/reads/mini.fastq --threads 2 --output search.output
## [01_introduction_snippet_3]

## [01_introduction_snippet_4]
seq -f "example_data/64/bins/bin_%02g.fasta" 0 1 63 > all_bin_paths.txt
raptor prepare --input all_bin_paths.txt --kmer 19 --window 23 --output precomputed_minimisers1 --threads 2
## [01_introduction_snippet_4]

## [01_introduction_snippet_5]
raptor build --input precomputed_minimisers1/minimiser.list --fpr 0.05 --threads 2 --output minimiser_raptor.index
## [01_introduction_snippet_5]

## [02_layout_snippet_1]
raptor layout --input-file all_bin_paths.txt
## [02_layout_snippet_1]

## [02_layout_snippet_2]
seq -f "example_data/1024/bins/bin_%04g.fasta" 0 1 1023 > all_bin_paths.txt
## [02_layout_snippet_2]

## [02_layout_snippet_3]
raptor layout --input-file all_bin_paths.txt
## [02_layout_snippet_3]

## [02_layout_snippet_4]
raptor layout --input-file all_bin_paths.txt \
              --kmer-size 17 \
              --num-hash-functions 4 \
              --false-positive-rate 0.25 \
              --output-filename binning.layout
## [02_layout_snippet_4]

## [02_layout_snippet_5]
raptor layout --input-file all_bin_paths.txt \
              --kmer-size 16 \
              --num-hash-functions 3 \
              --false-positive-rate 0.1 \
              --threads 2 \
              --output-filename binning2.layout
## [02_layout_snippet_5]

## [03_index_snippet_1]
raptor build --output raptor.index --input all_bin_paths.txt
## [03_index_snippet_1]

## [03_index_snippet_2]
echo -e ">chr1\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGCGTTCATTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" > mini_1.fasta
echo -e ">chr2\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGCGTCATTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" > mini_2.fasta
echo -e ">chr3\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" > mini_3.fasta
echo -e "mini_1.fasta\nmini_2.fasta\nmini_3.fasta" > all_paths.txt
## [03_index_snippet_2]

## [03_index_snippet_3]
raptor build --output raptor.index --input all_paths.txt
## [03_index_snippet_3]

## [03_index_snippet_4]
raptor build --kmer 19 --fpr 0.01 --output raptor.index --input all_bin_paths.txt
## [03_index_snippet_4]

## [03_index_snippet_5]
raptor build --kmer 4 --hash 3 --fpr 0.01 --output raptor2.index --input all_paths.txt
## [03_index_snippet_5]

## [03_index_snippet_6]
raptor prepare --kmer 19 --window 23 --output precomputed_minimisers2 --input all_paths.txt
raptor build --output minimiser_raptor.index --input precomputed_minimisers2/minimiser.list
## [03_index_snippet_6]

## [03_index_snippet_7]
raptor prepare --kmer 4 --window 6 --output precomputed_minimisers3 --input all_paths.txt
raptor build --output minimiser.index --input precomputed_minimisers3/minimiser.list
## [03_index_snippet_7]

## [03_index_snippet_8]
raptor build --input binning.layout \
             --fpr 0.25  \
             --threads 2 \
             --output hibf.index
## [03_index_snippet_8]

## [03_index_snippet_9]
raptor build --input binning.out --output hibf.index
raptor build --input binning2.layout --kmer 16 --hash 3 --fpr 0.1 --output hibf2.index
## [03_index_snippet_9]

## [04_search_snippet_1]
seq -f "example_data/64/bins/bin_%02g.fasta" 0 1 63 > all_bin_paths.txt
raptor build --input all_bin_paths.txt --output raptor.index
raptor search --index raptor.index --error 2 --query example_data/64/reads/mini.fastq --output search.output
## [04_search_snippet_1]

## [04_search_snippet_2]
echo -e ">query1\nCGCGTTCATT\n>query2\nCGCGTCATTA" > query.fasta
raptor search --index raptor2.index --query query.fasta --output search2.output
## [04_search_snippet_2]

## [04_search_snippet_3]
raptor search --index raptor2.index --query query.fasta --error 1 --output search3.output
raptor search --index raptor2.index --query query.fasta --threshold 0.9 --output search4.output
## [04_search_snippet_3]

## [04_search_snippet_4]
raptor search --index minimiser.index --query query.fasta --tau 0.9 --p_max 0.9 --output search5.output
## [04_search_snippet_4]

## [04_search_snippet_5]
raptor search --index hibf.index --query example_data/1024/reads/mini.fastq --output search6.output
## [04_search_snippet_5]
