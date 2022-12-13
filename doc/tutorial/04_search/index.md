# Search with Raptor {#tutorial_search}

You will learn how to search a Raptor index using one or more sequences queries.

\tutorial_head{Easy, 30 min, \ref tutorial_index, }

[TOC]

# Search an index using one or more sequences queries

## A first search

Now we want to have a look at Raptor's second method: `raptor search`.

We have already learned how to build an index in the previous tutorial on the `raptor index` and now we want to learn
how to search in it.
But before that, we delve deeper into the various possibilities and parameters. We can start a first search with the
indexes we have already created.

For this we just need the index and a query to run
`raptor search --index raptor.index --query query.fasta --output search.output`.

A result from such a call can look as follows:
```text
#0	example_data/64/bins/bin_00.fasta
#1	example_data/64/bins/bin_01.fasta
...
#62	example_data/64/bins/bin_62.fasta
#63	example_data/64/bins/bin_63.fasta
#QUERY_NAME	USER_BINS
0	0
1	0
2	0
3	0
4	0
16384	1
...
1015812	62
1032192	63
1032193	63
1032194	63
1032195	63
1032196	63
```
The output starts with a header section (lines starting with `#`). The header maps a number to each input file.
After the header section, each line of the output consists of the read ID (in this example these are numbers) and the
corresponding bins in which they were found.

\assignment{Assignment 1: Create a first search on our created example files and their raptor index}
From the first assignment in the last tutorial your directory should look like this:
```bash
tmp$ ls
all_paths.txt   mini_1.fasta    mini_2.fasta    mini_3.fasta    raptor.index    raptor2.index    ...
```
Lets search for the queries `CGCGTTCATT` amd `CGCGTCATT` with the first two indexes, with creating a `search.output` and
a `search2.output`.
\endassignment

\solution
Your `query.fasta` should look like:
```fasta
>query1
CGCGTTCATT
>query2
CGCGTCATT
```
And you should have run
```bash
raptor search --index raptor.index --query query.fasta --output search.output
raptor search --index raptor2.index --query query.fasta --output search2.output
```
Your `search.output` should look like:
```text
#0      mini_1.fasta
#1      mini_2.fasta
#2      mini_3.fasta
#QUERY_NAME     USER_BINS
query1  0,1,2
query2  0,1,2
```
and the other `search2.output` like this:
```text
#0      mini_1.fasta
#1      mini_2.fasta
#2      mini_3.fasta
#QUERY_NAME     USER_BINS
query1  0
query2
```
For the `search.output` Raptor will not output usefull results, as the index was build by the default kmer size `20`,
which is smaller than the length of the queries. This should not happen in normal cases.
For `search2.output` lets have a look again on the indexed sequences and the queries:
```fasta
>chr1
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGCGTTCATTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>query2                              ||||||||||
                                     CGCGTTCATT
>chr2
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGCGTCATTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>query2                              |||||||||
                                     CGCGTCATT
>chr3
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
```
So, with no errors allowed, we would expect a match from query1 to the first sequence and from query2 to the second
sequence. Actually, this should also appear in the output, we assume that due to the smallness of the example, the
second match does not appear.
\endsolution

\note
The Raptor search also has a help page, which can be accessed as usual by typing `raptor search -h` or
`raptor search --help`. And also an advanced help page: `raptor search -hh` or `raptor search --advanced-help`.

## General Idea & main parameters

\image html Raptor_no_title.svg width=90%

With `raptor index` we stored a representative transformation of the k-mer content of the database that is divided up
into a number of bins, typically a few hundred to a few thousand.
In `raptor search` it will combine the binning bitvectors of all representative k-mers in a query into a counting vector
and apply a thresholding step to determine the membership of a query in a bin.

\note
The term representative indicates that the k-mer content could be transformed by a function which reduces its size and
distribution, e.g. using minimizers.

There are also a few parameters for the search, which we will now explain.

Instead of an exact search, we can also allow errors in the query, this we achieve with a positive value of `--error`.
An error of 1 means either a deletion, an insertion or a substitution of a base pair (edit distance).

Raptor needs a median query size, which it calculates by default. With the parameter `--query_length` we can specify this
instead and thus avoid this calculation.

### Thresholding

If a query is searched for in a bin, this search is performed according to a given thershold. Depending on the setting,
we can influence this through various parameters (`--threshold`, `--tau`, `--p_max`, `--cache-thresholds`).

Thus the threshold can be either user-defined or computed for a given amount of allowed **errors** in the query. For the
latter, we rely on the well-known k-mer lemma. However, when winnowing minimizers are used, we apply a probabilistic
model that accounts for how errors might destroy relevant k-mers. The model returns a threshold value given how many
relevant k-mers are in the query. We have also further refined the model to incorporate the effect of the false positive
rate on the threshold. The thresholding step is a convenience for the user missing in other tools. If no minimizers are
used, our thresholding ensures that the (H)IBF gives no false negative answers.

With the parameter `--threshold` we can set a threshold ourselves, otherwise the **k-mer counting lemma** or the
probabilistic model is applied as described above. This limit value is to be seen as a percentage and thus to be
indicated as a number between `0` and `1`.

\note
**K-mer counting lemma**: For a given `k` and number of errors `e`, there are \f$k_p = |p| - k + 1\f$ many k-mers in the
pattern `p` and an approximate occurrence of `p` in text `T` has to share at least \f$t = (k_p - k \cdot e)\f$ k-mers. \n
**Example:** For the query ACGTT und `k=2` we have \f$(5 - 2 + 1) = 4\f$ k-mers (AC CG GT TT). With an error `= 1` we
have \f$4 - (2 \cdot 1) = 2\f$ k-mers which have to exist for matching the query. \n
Alternative: Given a threshold of `50%` \f$0.5 \cdot 4 = 2\f$ k-mers have to exist in a bin for matching the query.

\assignment{Assignment 2: Search with an error or a threshold.}
From the first assignment in the last tutorial we want to use the more usefull `raptor2.index` and do two searches on it.

Lets search again for the queries `CGCGTTCATT` and `CGCGTCATT` one time with an error of `1` and one time with the
threshold of `0.9`, with creating a `search3.output` and a `search4.output`.
\endassignment

\solution
You should have run
```bash
raptor search --index raptor2.index --query query.fasta --error 1 --output search3.output
raptor search --index raptor2.index --query query.fasta --threshold 0.9 --output search4.output
```
Your `search3.output` should look like:
```text
#0      mini_1.fasta
#1      mini_2.fasta
#2      mini_3.fasta
#QUERY_NAME     USER_BINS
query1  0,1
query2  0,1
```
and the other `search4.output` like this:
```text
#0      mini_1.fasta
#1      mini_2.fasta
#2      mini_3.fasta
#QUERY_NAME     USER_BINS
query1  0
query2  1
```
For the `search3.output` lets have a look again on the first two sequences and the queries:
```fasta
>chr1
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGCGTTCATTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>query2                              ||||||||||
                                     CGCGTTCATT
>chr2                                |||| |||||
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGCG-TCATTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>query2                              |||| |||||
                                     CGCG-TCATT
>chr1                                |||| |||||
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGCGTTCATTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
A
```
As you can see, the queries now also match with an error the other sequence, which means that our search result makes
sense.
For the `search4.output` we lose these matches again, because the threshold throws out these results.
\endsolution

#### Minimizers and thresholding

As written before, when winnowing minimizers are used, we apply a probabilistic model that accounts for how errors might
destroy relevant k-mers. The model returns a threshold value given how many relevant k-mers are in the query. We have
also further refined the model to incorporate the effect of the false positive rate on the threshold.

With the parameter `--threshold` we can switch off this model and use a simpler method. Alternatively, we can use the parameters `--tau` and `--p_max` to influence the probabilistic model.

In order to understand these parameters in more detail, we will go into the effects of minimisers a little more in the
following.

We use the counts of the `(w, k)`-minimizers and the probabilistic threshold derived to decide whether a query is in a
bin or not. However, having a relatively high false positive rate will affect this thresholding.

Using minimizers can largely reduce the index size and query time, but had a negative impact on the accuracy. The larger
the difference of window size (`w`) to k-mer size (`k`), the higher is the compression, which lowers the threshold for
the query membership. This results in an increasing number of false positives, thus decreasing the accuracy. False
negatives are rare and only slightly affected by the choice of `(w, k)`.

So for `(23, 19)`-minimizers, we would expect a compression by a factor of `3`. For `(32, 18)`-minimizers, we would
expect a compression factor of `8`. So using a larger `w` is beneficial for saving space. On the other hand, the
threshold for determining a match becomes smaller (roughly by the compression factor). We want to avoid having a
threshold of only `1` or `2`, since then, a few false positives k-mer matches would result in a false positive answer.

For a false positive rate of `0.05`, we encounter reads that have, e.g., `20` `(38, 20)`-minimizers. Those reads have a
chance of almost `40 %` to have one false positive count and a chance of almost `19 %` to have two. This has a small
effect on completely random hits. However, hits that match the query with more than the allowed errors could reach the
threshold due to one or two false positive minimizers. To counter this, we introduce a correction term for the
thresholds. We add a constant `c` to the threshold `t(x)`, which is determined by increasing `c` until the probability
`bp(x, c)` drops below a threshold `pmax`.

For example, for `p = 0.05` and `pmax = 0.15`, we have a correction of `+1` for `14 ≤ x ≤ 16` and a correction of `+2`
for `17 ≤ x ≤ 33`. For `p = 0.02`, we have a correction of `+1` for all `x`. The value for `pmax` was experimentally
determined, such that the benchmark yielded no false positives and no false negatives. It is set to `pmax = 0.15`. Those
corrections are precomputed and incorporated into the thresholding.

| \f$x\f$ | \f$\#x\f$ | \f$\%x\f$ | \f$t(x)\f$ | \f$c_{0.05}\f$ | \f$b_{0.05}(x,1)\f$ | \f$b_{0.05}(x,2)\f$ | \f$b_{0.05}(x,3)\f$ |
|:---:|--------:|------:|:------:|:----------:|:---------------:|:---------------:|:---------------:|
|  14 |       6 |   0.1 |    5   |      1     |       35.9      |       12.3      |       2.6       |
|  15 |     214 |   0.1 |    6   |      1     |       36.6      |       13.5      |       3.1       |
|  16 |   2,059 |   0.2 |    6   |      1     |       37.1      |       14.6      |       3.6       |
|  17 |  11,081 |   1.1 |    8   |      2     |       37.4      |       15.8      |      *4.1*      |
|  18 |  36,651 |   3.5 |    9   |      2     |       37.6      |       16.8      |      *4.7*      |
|  19 |  83,748 |   8.0 |    9   |      2     |       37.7      |       17.9      |      *5.3*      |
|  20 | 139,864 |  13.3 |   10   |      2     |       37.7      |       18.9      |      *6.0*      |
|  21 | 179,962 |  17.2 |   11   |      2     |       37.6      |       19.8      |      *6.6*      |
|  22 | 185,842 |  17.7 |   12   |      2     |       37.5      |       20.7      |      *7.3*      |
|  23 | 158,032 |  15.1 |   12   |      2     |       37.2      |       21.5      |      *7.9*      |
|  24 | 113,696 |  10.8 |   13   |      2     |       36.9      |       22.3      |      *8.6*      |
|  25 |  70,089 |   6.7 |   14   |      2     |       36.5      |       23.1      |      *9.3*      |
|  26 |  37,540 |   3.6 |   15   |      2     |       36.1      |       23.7      |     *10.0*      |
|  27 |  18,040 |   1.7 |   15   |      2     |       35.6      |       24.3      |     *10.7*      |
|  28 |   7,535 |   0.7 |   16   |      2     |       35.0      |       24.9      |     *11.4*      |
|  29 |   2,790 |   0.3 |   17   |      2     |       34.5      |       25.4      |     *12.0*      |
|  30 |     961 |  <0.1 |   18   |      2     |       33.9      |       25.9      |     *12.7*      |
|  31 |     343 |  <0.1 |   18   |      2     |       33.3      |       26.3      |     *13.4*      |
|  32 |      88 |  <0.1 |   19   |      2     |       32.6      |       26.6      |     *14.0*      |
|  33 |      28 |  <0.1 |   20   |      2     |       32.0      |       26.9      |     *14.6*      |
|  34 |       5 |  <0.1 |   22   |      3     |       31.3      |       27.2      |      15.3       |
|  35 |       2 |  <0.1 |   23   |      3     |       30.6      |       27.4      |      15.8       |

Table: **Exemplary threshold distribution.** The values are for `(38, 20)`-minimizers, `2` errors and read length `250`.
Shown are the distribution `#x` (`\%x`) of the number (percentage) of minimizers reaching from `x = 14` to `x = 35` and
the threshold `t(x)` using the probabilistic model. The threshold `t(x)` incorporates the correction term \f$c_p\f$. On
the right, you see the probability \f$b_p(x, a)\f$ of having a false positive answers from the IBF for
`p = 0.05`.

- \f$x\f$: number of minimizers in a query
- \f$a\f$: number of false positive answers
- \f$b_p(x, a)\f$: probability of returning a false positive answers when querying `x` minimizers in a Bloom filter of a
                   false positive rate
- \f$t(x)\f$: threshold for given `x`
- \f$c\f$: (constant) correction that is added to the threshold
- \f$pmax\f$: threshold for correction term; i.e., increase c until \f$b_p(x, a) < pmax\f$

Now, with all these explanations, let's look at the parameters `--tau` and `--p_max`. Like the threshold parameter, both
can be specified as a percentage (i.e. as a number between 0 and 1).

\note
The threshold is then calculated automatically, so the parameter `--threshold` is not needed here and would disable the
probabilistic model.

With the commulative probability `--tau` we influence how many percent of the cases that are possible we want to cover,
i.e. where in the binomial distribution of `#x` (*bell curve* in table) that is created above is cut off.

With the parameter `--pmax` you can limit the percentage of false positives as described, e.g. a pmax value of `0.15`
means that you want to have a probability of no more than `15%` false positives. That means the other way round, with a
higher pmax value the results are more sensitive.

Using `--cache-thresholds` stores the computed thresholds with an unique name next to the index. In the next search call
using this option, the stored thresholds are re-used. This saves time because the model does not have to be calculated
again.
Two files are stored:
- `threshold_*.bin`: Depends on pattern, window, kmer/shape, errors, and tau.
- `correction_*.bin`: Depends on pattern, window, kmer/shape, p_max, and false positive rate.

\assignment{Assignment 3: Search with minimizers.}
We want to use the `minimiser.index` from the index tutorial assignment again and use the same queries `CGCGTTCATT` and
`CGCGTCATT` to search in it.

Lets search now with a tau of `0.9` and a pmax of `0.9`, with creating a `search5.output`.
\endassignment

\solution
You should have run
```bash
raptor search --index minimiser.index --query query.fasta --tau 0.9 --p_max 0.9 --output search5.output
```
Your `search5.output` should look like:
```text
#0      mini_1.fasta
#1      mini_2.fasta
#2      mini_3.fasta
#QUERY_NAME     USER_BINS
query1  0
query2
```
\endsolution

## Others

### Parallelization

Raptor also supports parallelisation in search. By specifying `--threads`, for example, the querys are processed
simultaneously.

### Time

If you want to know the duration of raptor search including *time for reading the index*, *time for reading reads* and
*the computation time* for benchmarking, you can use `--time` to write a timing file.

### HIBF

If a HIBF was previously calculated as an index, this must also be specified for the search with `--hibf` again.

For the HIBF we look up the k-mers of a query and apply a refined threshold on each level.

Notably, having a threshold allows us to skip traversing lower level IBFs whose upper level merged bins already do not
exceed the threshold. In practice, this means that we only have to access a fraction of the HIBF.

If no minimizers are used, our thresholding ensures that the (H)IBF gives no false negative answers. Using minimizers,
there are few false negatives, because the minimizer compression is not lossless.

\assignment{Assignment 4: Search with hibf.}
We want to use the `hibf.index` from the index tutorial assignment again. <!-- and use X as X-parameter. -->

Lets search now for the `1024/reads/mini.fastq` querries with creating a `search6.output`.
\endassignment

\solution
You should have run
```bash
raptor search --hibf --index hibf.index --query 1024/reads/mini.fastq --output search6.output
```
Your `search6.output` should look like:
```text
#0      1024/bins/bin_0712.fasta
#1      1024/bins/bin_0406.fasta
...
#1021   1024/bins/bin_0533.fasta
#1022   1024/bins/bin_0624.fasta
#QUERY_NAME     USER_BINS
0
1
...
1047555
1047556
```
\note
You can also calculate the second hibf index B and you will see that they will give slightly different results
(`diff search6.output search7.output`).
\endsolution
