# HIBF in Nextflow

If you want to include the HIBF in Nextflow, you can adapt this sample workflow `hibf.nf`.

It is required to have the application [raptor](https://github.com/seqan/raptor) and
[chopper](https://github.com/seqan/chopper) available in your `$PATH`.
We will work on docker images meanwhile.

### The input data

Given the above tool requirements, the workflow is running on two inputs

1. The database `data/small_genomes/`
2. The query file `data/queries.fq`

You can swap these input data with anything that fits your needs, but be aware of the following:

1. The database is always a directory where each file that ends in [fa,fna,fasta,fq,fastq], possibly
   suffixed with `.gz`, is treated as an individual bin/sample to index. The search result will contain
   a membership answer for each query sequence in each of these files. Each file may contain multiple
   sequences
2. The query file can be either in FASTA or FASTQ format. The file may contain multiple sequences.

As an example, the directory `data/small_genomes` contains the 300 smallest RefSeq genomes (as of Jan. 2022)
and the file `data/queries.fq` is a simulated read file, generated from the former mentioned genomes.

### The parameters

Some details on the input parameters you can adjust

* `params.reads` The input queries to search for in the index
* `params.data` The database directory.
* `params.tmax` An index layout parameter. It is best chosen as the square root of the number of
                bins/files/sequences to index.
* `params.kmer_size` The kmer size the algorithms uses to abstract the data for efficiency. A kmer
                     size of 20-35 works well with DNA.
* `params.number_of_hash_functions` The number of hash functions of the bloom filter in the data
                                    structures. [2,3,4] has worked well in practice.
* `params.fpr` The false positive rate you can cope with in your analysis. The higher the FPR the
               more hits in the index will be reported that are not true hits.
* `params.errors` How many errors are allowed when searching a query in the index.
* `params.outdir` The directory name of the output.
