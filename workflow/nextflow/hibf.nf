#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.reads = "data/queries.fq"
params.data = "https://ftp.imp.fu-berlin.de/pub/seiler/raptor/small_genomes.tar.gz"
params.tmax = 192 /* Is best chosen as the square root of the number of bins/files/sequences to index. */
params.kmer_size = 20
params.number_of_hash_functions = 2 /* [2,3,4] has worked well in practice */
params.fpr = 0.05
params.errors = 2
params.outdir = "hibf_result"

println "\nResult files will be written to $params.outdir\n"


workflow {
    def query_ch = Channel.fromPath(params.reads)

    DOWNLOAD(params.data)
    UNTAR(DOWNLOAD.out)

    DATA_PREP(UNTAR.out)
    LAYOUT(DATA_PREP.out, params.kmer_size, params.fpr, params.tmax)
    INDEX(LAYOUT.out, params.kmer_size, params.fpr, params.number_of_hash_functions)
    SEARCH(INDEX.out, query_ch, params.fpr, params.errors)
}

process DOWNLOAD {
    tag "Only needed for the small_genomes.tar.gz example on Github"

    input:
    val data

    output:
    path "small_genomes.tar.gz"

    script:
    """
    set -x
    wget $data
    """
}

process UNTAR {
    tag "Only needed for the small_genomes.tar.gz example on Github"

    input:
    path data

    output:
    path "${data.getSimpleName()}"

    script:
    """
    set -x
    tar zxf $data
    """
}

process DATA_PREP {
    tag "Create an info file about the data to index"

    input:
    val data

    output:
    path 'input_data.tsv'

    script:
    """
    set -x # commands are printed to stderr for logging
    find $data -type f | grep -E '*\\.f[a,q,na,asta,astq].*' > input_data.tsv
    wc -l input_data.tsv
    """
}

process LAYOUT {
    tag "Creating a layout on the input_data Database"

    input:
      path refseq
      val kmer_size
      val fpr
      val tmax

    output:
      path "hibf.layout"

    script:
    """
    set -x # commands are printed to stdout for logging
    raptor layout --input-file $refseq \
                  --kmer-size "$kmer_size" \
                  --threads $task.cpus \
                  --tmax $tmax \
                  --num-hash-functions 3 \
                  --false-positive-rate "$fpr" \
                  --output-filename hibf.layout \
                  --rearrange-user-bins
    """
}

process INDEX {
    tag "Building index"
    publishDir params.outdir

    input:
      path layout
      val kmer_size
      val fpr
      val num_hash_fn

    output:
      path 'hibf.index'

    script:
    """
    raptor build --kmer "$kmer_size"  \
                 --window "$kmer_size"  \
                 --hash $num_hash_fn  \
                 --fpr "$fpr"  \
                 --threads $task.cpus  \
                 --output hibf.index  \
                 --hibf $layout
    """
}

process SEARCH {
    tag "search"
    publishDir params.outdir

    input:
      path index
      path reads
      val fpr
      val errors

    output:
      path 'hibf_search.out'

    script:
    """
    raptor search --index $index \
                  --query $reads \
                  --output hibf_search.out \
                  --error $errors \
                  --threads $task.cpus \
                  --p_max 0.4 \
                  --tau 0.99 \
    """
}
