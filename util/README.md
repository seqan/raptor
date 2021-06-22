# Raptor Utilities

This subfolder contains scripts and documentation related to the experiments run in the paper.

## Prerequisites

The base requirements are the same as for the main Raptor application:

* CMake >= 3.8
* GCC 7, 8, 9 or 10 (most recent minor version)
* git

Refer to the [SeqAn3 Setup Tutorial](https://docs.seqan.de/seqan/3-master-user/setup.html) for more in-depth
information.

To cluster a data set by taxonomy, you will also need the following tools:

* [genome_updater](https://github.com/pirovc/genome_updater)
* [TaxSBP](https://github.com/pirovc/taxsbp)

In case, the bash expression `some/path/bin_{00..63}.fasta` or similar does not work for your shell, an alternative is
to use `$(BIN_NUMBER=64; seq -f "some/path/bin_%0${#BIN_NUMBER}g.fasta" 0 1 $((BIN_NUMBER-1)))` instead. Replace
`BIN_NUMBER` as needed.

## Compiling Applications

From within the `util` directory:
```bash
mkdir build
cd build
cmake ..
make -j2 install
```

`-j` determines the number of threads to use for building, `-j2` will use two threads.

### Options

To build Mason for simulating a data set, you can pass `-DRAPTOR_UTILITY_BUILD_MASON=1` to `cmake`.<br>
To build DREAM-Yara for mapping reads, you can pass `-DRAPTOR_UTILITY_BUILD_DREAM_YARA=1` to `cmake`.

In above example, the `cmake` call would change to:
```bash
cmake -DRAPTOR_UTILITY_BUILD_MASON=1 -DRAPTOR_UTILITY_BUILD_DREAM_YARA=1 ..
```
Please note that, in contrast to most other `cmake` options, these two options must be passed to `cmake` before the
directory (`..` in this case).
Enabling the DREAM-Yara build may increase the build time significantly.

## Content

### Included Applications

After building all applications, the `bin` directory within the build directory will contain a number of executables:
```console
$ tree bin
bin/
├── apply_taxsbp             # Splits a data set according to the clustering obtained from TaxSBP
├── count_minimiser          # Counts minimiser in a list of files (individual + combined)
├── dream_yara_build_filter  # [DREAM-Yara] Builds the IBF
├── dream_yara_indexer       # [DREAM-Yara] Builds the FM-Indices
├── dream_yara_mapper        # [DREAM-Yara] Maps reads
├── generate_reads           # Simulates reads from an artifical data set
├── generate_reads_refseq    # Simulates reads from any data set
├── mason_genome             # [Mason] Generates a random genome
├── mason_variator           # [Mason] Generates haplotypes of a genome
└── split_sequence           # Splits one big genome into smaller genomes
```

### Included bash scripts

`src/bash_scripts` contains a variety of bash scripts:
```console
$ tree src/bash_scripts/
src/bash_scripts/
├── benchmark.sh           # Runs benchmarks for Raptor
├── count_minimisers.sh    # Counts a variety of different minimizers for a data set
├── dream_yara.sh          # Runs benchmarks for DREAM-Yara
├── original_dream_yara.sh # Runs benchmarks for the original DREAM-Yara
├── run_squeakr.sh         # Runs Squeakr on input files
└── simulate.sh            # Simulates an artificial data set
```

[DREAM-Yara](https://github.com/seqan/dream_yara/tree/raptor_ibf) refers to the version that uses the IBF and minimizers
from Raptor and is the version that is built in the included applications.<br>
[Original DREAM-Yara](https://github.com/seqan/dream_yara/tree/master) refers to the previous version and is not built
as part of the included applications.

### Included evaluation scripts

`src/evaluation_scripts` hosts various python scripts used to generate tables in the CSV format:
```console
$ tree src/evaluation_scripts/
src/evaluation_scripts/
├── eval_benchmark.py  # Gathers information from Raptor benchmarks
├── eval_yara.py       # Gathers information from DREAM-Yara benchmarks
├── get_counts.py      # Consolidates multiple count_minimiser runs
├── get_frequencies.py # Computes minimiser frequencies from count_minimiser runs
└── get_thresholds.py  # Gathers minimiser thresholds from Raptor benchmarks
```

## How-to

### Generate artificial data set

`src/bash_scripts/simulate.sh` can be used to simulate an artificial data set. Variables in upper case are user-provided
and may be changed.

#### Notes

* `BINARY_DIR` is the absolute path to the `build/bin` directory
* `OUT_DIR` is the absolute path to the output directory
* `LENGTH % BIN_NUMBER` must be 0
* `(LENGTH / BIN_NUMBER) % HAPLOTYPE_COUNT` must be 0
* `READ_COUNT % BIN_NUMBER` must be 0

The easiest way to achieve the last three requirements is to set `LENGTH`, `BIN_NUMBER`, and `READ_COUNT` to a power
of two.

#### Example

```bash
OUT_DIR=/some/path
BIN_NUMBER=1024
ERRORS=2
READ_LENGTHS="100 150 250"
```
will result in
```
/some/path/           # OUT_DIR
└── 1024              # BIN_NUMBER
    ├── bins          # Contains 1024 FASTA files
    ├── info          # Documents mutations introduced for the haplotypes
    ├── reads_e2_100  # Reads with 2 errors and length 100
    ├── reads_e2_150  # Reads with 2 errors and length 150
    └── reads_e2_250  # Reads with 2 errors and length 250
```

The `reads_*` directories will contain one FASTQ file per bin as well as an `all.fastq` representing the
concatenation of all reads. The other additional file, `all10.fastq`, contains the content of `all.fastq` ten times.

### Download and cluster the current NCBI RefSeq

#### Download the NCBI RefSeq

We recommend using the [genome_updater](https://github.com/pirovc/genome_updater).
To download all archaea and bacteria, use:
```bash
./genome_updater.sh -g "archaea,bacteria" \
                    -d "refseq" \
                    -l "Complete Genome" \
                    -f "genomic.fna.gz,assembly_report.txt" \
                    -o "RefSeqCG_arc_bac" -b "v1" \
                    -a -m -u -r -p -t 32
```
You may have to adapt `./genome_updater.sh` to point to the correct location. `-t 32` enables the use of 32 threads for
downloading, decrease the number as needed.
In case individual downloads fail, the command can be rerun with an additional `-i` flag to only gather missing files.

The output will look similar to
```bash
$ tree -L 2
.
├── assembly_summary.txt -> v1/assembly_summary.txt
└── v1
    ├── <timestamp>.log
    ├── <timestamp>_url_downloaded.txt
    ├── <timestamp>_url_failed.txt
    ├── assembly_summary.txt
    ├── files
    ├── taxdump.tar.gz
    ├── updated_assembly_accession.txt
    └── updated_sequence_accession.txt
```

Please extract the `taxdump.tar.gz` by running `tar -xf taxdump.tar.gz` in the `v1` directory.

#### Running TaxSBP

[TaxSBP](https://github.com/pirovc/taxsbp) will generate a taxonomic clustering of the NBCI RefSeq.

For the following bash commands, we will assume that we are in a sibling directory of `v1`, like so:
```bash
$ tree -L 1
.
├── assembly_summary.txt -> v1/assembly_summary.txt
├── clustering # This will be our working directory for the following steps
└── v1
```

First, we will need to process the `updated_sequence_accession.txt` generated by `genome_updater` in the previous step:
```bash
awk 'BEGIN {FS="\t";OFS="\t"}{if($4!="na"){ print $4,$5,$6 }}' ../v1/updated_sequence_accession.txt > seqinfo.txt
```
This will extract RefSeq accession, assembly accession and taxonomic ID.

Since TaxSBP clusters approximately into the given number of bins, we will run it multiple times until we get the
desired number of bins.

```bash
goal=64; i=$goal; printf "Trying $i bins"\\r; while [[ $(taxsbp -i seqinfo.txt -n ../v1/taxdump/nodes.dmp -t -b $i | tail -n1 | cut -f 6) -ge $goal ]]; do ((i-=1)); printf "Trying $i bins"\\r; done; echo "Use -b $i to get at most $goal bins"
```
For 64 bins, the result may be that we need to run TaxSBP with 63 bins.

```bash
goal=1024; i=$goal; printf "Trying $i bins"\\r; while [[ $(taxsbp -i seqinfo.txt -n ../v1/taxdump/nodes.dmp -t -b $i | tail -n1 | cut -f 6) -ge $goal ]]; do ((i-=1)); printf "Trying $i bins"\\r; done; echo "Use -b $i to get at most $goal bins"
```
For 1024 bins, the result may be that we need to run TaxSBP with 1009 bins.

With this knowledge, we can run TaxSBP to get the desired clustering:

```bash
taxsbp -i seqinfo.txt -n ../v1/taxdump/nodes.dmp -b 63 -o 64.taxsbp
taxsbp -i seqinfo.txt -n ../v1/taxdump/nodes.dmp -b 1009 -o 1024.taxsbp
```

#### Splitting the RefSeq according to TaxSBP

We will use `apply_taxsbp` which can be found in the `bin` directory within the build directory:

```bash
mkdir clustered_64
<build_directory>/bin/apply_taxsbp --input ../v1/files/ \
                                   --output clustered_64/ \
                                   --taxsbp 64.taxsbp \
                                   --genome_update ../v1/updated_sequence_accession.txt \
                                   --assembly_summary ../v1/assembly_summary.txt \
                                   --threads 32
```

```bash
mkdir clustered_1024
<build_directory>/bin/apply_taxsbp --input ../v1/files/ \
                                   --output clustered_1024/ \
                                   --taxsbp 1024.taxsbp \
                                   --genome_update ../v1/updated_sequence_accession.txt \
                                   --assembly_summary ../v1/assembly_summary.txt \
                                   --threads 32
```

This will create the directories `clustered_64` and `clustered_1024` containing the clustered NCBI RefSeq.

#### Generate reads

We will use `generate_reads_refseq` which can be found in the `bin` directory within the build directory:

```bash
mkdir clustered_64/reads_e2
seq -f "clustered_64/bin_%02g.fasta.gz" 0 1 63 > all_bin_paths.txt
for read_length in 100 150 250
do
    <build_directory>/bin/generate_reads_refseq --output clustered_64/reads_e2/$read_length \
                                                --errors 2 \
                                                --read_length $read_length \
                                                --number_of_reads 1048576 \
                                                --threads 32 \
                                                all_bin_paths.txt
done
for read_length in 100 150 250
do
    cat clustered_64/reads_e2/$read_length/*.fastq > clustered_64/reads_e2/$read_length/all.fastq
done
```

```bash
mkdir clustered_1024/reads_e2
seq -f "clustered_1024/bin_%04g.fasta.gz" 0 1 1023 > all_bin_paths.txt
for read_length in 100 150 250
do
    <build_directory>/bin/generate_reads_refseq --output clustered_1024/reads_e2/$read_length \
                                                --errors 2 \
                                                --read_length $read_length \
                                                --number_of_reads 1048576 \
                                                --threads 32 \
                                                all_bin_paths.txt
done
for read_length in 100 150 250
do
    cat clustered_1024/reads_e2/$read_length/*.fastq > clustered_1024/reads_e2/$read_length/all.fastq
done
```

### Download the NCBI RefSeq used in the paper

The data can be downloaded here:
 * [NCBI RefSeq](https://doi.org/10.5281/zenodo.4647988)
 * [NCBI RefSeq clustered into 64 bins](https://doi.org/10.5281/zenodo.4650188)
 * [NCBI RefSeq clustered into 1024 bins](https://doi.org/10.5281/zenodo.4651078)
 * [Queries](https://doi.org/10.5281/zenodo.4651379)

In case there are multiple parts, i.e. files matching `<common_prefix>.part.00`, combine them into a single file first:
```bash
cat <common_prefix>.part* > <common_prefix>
```
`<common_prefix>` refers to the shared part of the file name, e.g., `Clustered_64_RefSeqCG_arc_bac.tar.zst`.

`tar.zst` files can be decompressed with `tar -I zstd -xf <filename>`.

### Build a DREAM-Yara index

#### Build the IBF

In this example, we will build an index for (31,19)-minimizers and a size of 2 GiB for 1024 bins:

```bash
<build_directory>/dream_yara_build_filter --output-file ibf.filter \
                                          --kmer-size 19 \
                                          --window-size 31 \
                                          --bloom-size 2G \
                                          --threads 32 \
                                          --num-hash 2 \
                                          --version-check 0 \
                                          clustered_1024/bin_{0000..1023}.fasta.gz
```

You may need to adjust the path to `clustered_1024`.

#### Build the FM-Indices

The FM-Indices only need to be built once for a data set, independent of which IBF parameters are used:

```bash
mkdir fm_indices
<build_directory>/dream_yara_indexer --output-prefix fm_indices/ \
                                     --threads 32 \
                                     --version-check 0 \
                                     clustered_1024/bin_{0000..1023}.fasta.gz
```

### Download the DREAM-Yara index used in the paper

The data can be downloaded here:
<pre>
<a href="https://ftp.imp.fu-berlin.de/pub/raptor/">https://ftp.imp.fu-berlin.de/pub/raptor/</a>
├── <a href="https://ftp.imp.fu-berlin.de/pub/raptor/1024_bins">1024_bins</a>
│   ├── <a href="https://ftp.imp.fu-berlin.de/pub/raptor/1024_bins/19_19_16G.filter">19_19_16G.filter</a>
│   ├── <a href="https://ftp.imp.fu-berlin.de/pub/raptor/1024_bins/23_19_4G.filter">23_19_4G.filter</a>
│   ├── <a href="https://ftp.imp.fu-berlin.de/pub/raptor/1024_bins/31_19_2G.filter">31_19_2G.filter</a>
│   └── <a href="https://ftp.imp.fu-berlin.de/pub/raptor/1024_bins/fm_indices.tar.zst">fm_indices.tar.zst</a>
└── <a href="https://ftp.imp.fu-berlin.de/pub/raptor/64_bins">64_bins</a>
    ├── <a href="https://ftp.imp.fu-berlin.de/pub/raptor/64_bins/19_19_16G.filter">19_19_16G.filter</a>
    ├── <a href="https://ftp.imp.fu-berlin.de/pub/raptor/64_bins/23_19_4G.filter">23_19_4G.filter</a>
    ├── <a href="https://ftp.imp.fu-berlin.de/pub/raptor/64_bins/31_19_2G.filter">31_19_2G.filter</a>
    └── <a href="https://ftp.imp.fu-berlin.de/pub/raptor/64_bins/fm_indices.tar.zst">fm_indices.tar.zst</a>
</pre>

`tar.zst` files can be decompressed with `tar -I zstd -xf <filename>`.

### Map reads to NCBI RefSeq via DREAM-Yara

Reads can now be mapped, in this example, we will map reads of length 250 to 1024 bins:
```bash
<build_directory>/dream_yara_mapper --bloom-filter ibf.filter \
                                    --output-file results.sam \
                                    --threads 32 \
                                    --error-rate 0$(bc -l <<< "2/250") \
                                    --verbose \
                                    --version-check 0 \
                                    fm_indices/ \
                                    clustered_1024/reads_e2/250/all.fastq
```

The `error-rate` is given as floating point number between 0 and 1, e.g. `-error-rate 0.05` for 5 per cent.
`--error-rate 0$(bc -l <<< "ERRORS/READ_LENGTH")` can be used to compute the error rate for `ERRORS` many errors in
reads of length `READ_LENGTH`. In above example, we search for `2` errors in reads of length `250`.

### Automation

The bash script `src/bash_scripts/dream_yara.sh` can be used to automate the building of indices and
mapping reads in DREAM-Yara. Adjustable parameters are found in the beginning and end of the script.

## Mantis and COBS
The following is a short guide on how we used Mantis and COBS in our paper.

### COBS
COBS requires an input directory `<input_directory>` containing FASTA or FASTQ files and an output filename
`<cobs_output_name>` to store the index.

Build COBS on the artificial data set:
```bash
cobs compact-construct --num-hashes 2 -k 19 -t 32 <input_directory> <cobs_output_name>.cobs_compact
```

Build COBS on the real data set:
```bash
cobs compact-construct -k 20 <input_directory> <cobs_output_name>.cobs_compact
```

Query COBS:
```bash
cobs query -i <cobs_output_name>.cobs_compact -f <query_file>
```

### Mantis
Mantis is based on squeakr files. Therefore, they need to be constructed beforehand. The squeakr files contain all
k-mers which have count values greater than or equal to a given cutoff. For the artificial data set, the default cutoff
of 1 was used. For the real data set, the cutoffs according to the Mantis paper were used. These cutoffs depend on the
file sizes, so the cutoffs differed between files. The bash script at `src/bash_scripts/run_squeakr.sh` adapts the
cutoff for each file individually.

Run squeakr on the artificial data set:
```bash
squeakr count -e -k 19 -n -o <out_name> <fastq-file>
```

Run squear on the real data set
```bash
squeakr count -e -k 20 -n -c <cutoff_based_on_file_size> -o <out_name> <fastq-file>
```

Mantis is then built using the squeakr files. The input is a file containing the paths to each squeakr file.

Build Mantis:
```
mantis build -s 34 -i <file_with_list_of_squeakr_files> -o <output_directory>
mantis mst -p <output_directory>
```

Query Mantis:
```
mantis query -k <k-mer_size> -p <output_directory> -o <output_file> <query_file>
```
