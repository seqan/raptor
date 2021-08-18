# 2.0.0

## Features
* The index additionally stores the original file paths, window size, k-mer size, number of parts, and the data layout
  ([\#34](https://github.com/seqan/raptor/pull/34), [\#47](https://github.com/seqan/raptor/pull/47),
   [\#57](https://github.com/seqan/raptor/pull/57), [\#58](https://github.com/seqan/raptor/pull/58)):
  * `--window`, `--kmer`, `--parts`, and `--compressed` are no longer supported when running `raptor search`.
  * Indices built with an older version of raptor cannot be used with version 2.0.0.
  * Old indices can be converted by running `raptor upgrade`.

* New output format ([\#47](https://github.com/seqan/raptor/pull/47)):
  * Output includes a header (lines starting with '\#').
  * The header assigns an integer to each input file, which are then used in the results:
    ```
    #0	/dev/shm/test/data/bin1.fa
    #1	/dev/shm/test/data/bin2.fa
    #2	/dev/shm/test/data/bin3.fa
    #3	/dev/shm/test/data/bin4.fa
    #QUERY_NAME	USER_BINS
    query1	0
    query2	0,1
    query3	2,3
    ```
  * The list of bins does no longer contain a trailing comma
    (the last line would have looked like this before 2.0.0: `query3	2,3,`).

* Added support for some parts of the [SOCKS Interface](https://gitlab.ub.uni-bielefeld.de/gi/socks)
  ([\#35](https://github.com/seqan/raptor/pull/35)).

## Misc
* Passing a list of bins is no longer supported. Use a file containing paths to the files instead
  ([\#35](https://github.com/seqan/raptor/pull/35)).

## Bug fixes
Except the last entry in this list, these bugs are expected to not be present in version 1.1.0.
* Validation error when using `--compute-minimiser` ([\#25](https://github.com/seqan/raptor/pull/25)).
* Typos/Errors in `util` scripts ([\#32](https://github.com/seqan/raptor/pull/32),
  [\#39](https://github.com/seqan/raptor/pull/39), [\#40](https://github.com/seqan/raptor/pull/40),
  [\#48](https://github.com/seqan/raptor/pull/48), [\#51](https://github.com/seqan/raptor/pull/51)).
* Partitioned IBF not working ([\#50](https://github.com/seqan/raptor/pull/50),
  [\#53](https://github.com/seqan/raptor/pull/53)).
* Rare segmentation fault in threshold computation ([\#54](https://github.com/seqan/raptor/pull/54)).

# 1.1.0

## Features
* Raptor accepts a text file containing the path to a bin on each line
  ([\#16](https://github.com/seqan/raptor/pull/16)).

## Bug fixes
* Threshold option not working ([\#14](https://github.com/seqan/raptor/pull/14)).

# 1.0.1

## Bug fixes

* Reduced the number of open file handles ([\#10](https://github.com/seqan/raptor/pull/10)).
