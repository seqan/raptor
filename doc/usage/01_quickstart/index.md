# Quickstart {#usage_quickstart}

[TOC]

* reference usage sections
* copy/paste examples with many/all options
* present subcommands
* present workflow ibf/hibf
* present (w,k) minimiser, (k,k) canonical k-mers

## What problems can be solved (short 5 sentences)

* Approximate Membership Query (AMQ)
* Input sequences = bins
* Given a query sequence: Show all bins that possible contain that sequence.
  * Either given a number of errors
  * Or given a percentage (e.g., 0.7 -> 70% of minimiser/k-mers of query must be present in bin)
* Does not model a colored de Bruijn graph:
  * ✔️ Given a k-mer: Tell to which color (bin) it belongs to
  * ❌ Given a color (bin): Show all k-mers that belong to this color (bin)

## HIBF vs IBF

* comparison table

\htmlonly
<!--
    Including an SVG like this allows using CSS variables inside the SVG, e.g.,
    fill="var(--page-foreground-color)"
    fig_4_many_user_bins.svg#svg refers to id="svg" in the svg file
    viewBox of svg element is copied from svg file
    svg file needs to be added to HTML_EXTRA_FILES
-->
<svg style="height: auto; width: 100%;" viewBox="0 0 342.22876 198.804504">
  <use href="/fig_4_many_user_bins.svg#svg" width="100%" height="100%"></use>
</svg>
\endhtmlonly

## Workflow

<div class="tabbed">

- <b class="tab-title">HIBF</b>
[`raptor prepare`] -> `raptor layout` -> `raptor build` -> `raptor search`

- <b class="tab-title">IBF</b>
[`raptor prepare`] -> `raptor build` -> `raptor search`

</div>

## (w,k) minimiser vs (k,k) canonical k-mers

* comparison table
  * k-mers: exact threshold when giving number of errors
  * w,k: compression. Reduces number of values to store (and hence size of index) by roughly \f$\frac{w - k + 2}{2}\f$

\htmlonly
<svg style="height: auto; width: 100%;" viewBox="0 0 307.29294 196.59375">
  <use href="/fig_6_acc.svg#svg" width="100%" height="100%"></use>
</svg>
\endhtmlonly
