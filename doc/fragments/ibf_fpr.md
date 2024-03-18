<!--
SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
SPDX-License-Identifier: CC-BY-4.0
-->

Sets an upper bound for Bloom Filter false positives.<br>
**Recommendation**: default value (`0.05`)<br>
&nbsp;&nbsp;• A lower `fpr` limits the number of false-positive results, but increases index size.<br>
&nbsp;&nbsp;• A higher `fpr` can help to reduce memory consumption in cases where false-positive k-mers have little effect.<br>
See also: [Bloom Filter Calculator](https://hur.st/bloomfilter/).
