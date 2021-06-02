# 2.0.0

## Features
* The index now also stores the window and k-mer size ([\#34](https://github.com/seqan/raptor/pull/34)):
  * This means that old indices cannot be used with a newer raptor version.
  * Old indices can be converted by pre-release todo.

## Misc
* Passing a list of bins is no longer supported. Use a file containg paths to the files instead
  ([\#35](https://github.com/seqan/raptor/pull/35)).

# 1.1.0

## Features
* Raptor accepts a text file containing the path to a bin on each line
  ([\#16](https://github.com/seqan/raptor/pull/16)).

## Bug fixes
* Threshold option not working ([\#14](https://github.com/seqan/raptor/pull/14)).

# 1.0.1

## Bug fixes

* Reduced the number of open file handles ([\#10](https://github.com/seqan/raptor/pull/10)).
