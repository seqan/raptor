# Index with Raptor {#tutorial_index}

You will learn how to construct a Raptor index of large collections of nucleotide sequences.

\tutorial_head{Easy, 30 min, \ref tutorial_first_steps, }

[TOC]

# Index of large collections of nucleotide sequences

ToDo ...

## EXAMPLES

`raptor build --kmer 19 --window 23 --size 8m --output raptor.index all_bin_paths.txt`

`raptor build --kmer 19 --window 23 --compute-minimiser --output precomputed_minimisers all_bin_paths.txt`

`raptor build --size 8m --output minimiser_raptor.index all_minimiser_paths.txt`
