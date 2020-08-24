Slant.R
=======

Overview
--------

``slanter`` contains a set of functions for reordering data, and generating hierarchical clustering
for ordered data, for improved visualization.

See the [published R package](https://CRAN.R-project.org/package=slanter) or the [latest github
version documentation](https://tanaylab.github.io/slanter/index.html) for details. Specifically, the
[meristems vignette](https://tanaylab.github.io/slanter/articles/meristems.html) explains why and
how to use this package.

Installation
------------

To install it, use:

``` r
install.packages('slanter')
```

Usage
-----

In general, if your data is a similarity matrix (each entry is a non-negative value that indicates
how similar a pair of elements is to each other, higher is better), then use `slanter::sheatmap` as
a drop-in replacement for `pheatmap::pheatmap`, and enjoy.

The lower level function `slanter::slanted_orders` will compute the visualization order and
`slanter::oclust` will compute hierarchical clustering that is consistent with this order. See the
[published R package](https://CRAN.R-project.org/package=slanter) or the latest github version
documentation [reference section](https://tanaylab.github.io/slanter/reference/index.html) for the
list of available functions.
