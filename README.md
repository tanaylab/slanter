Slant.R
=======

Overview
--------

``slanter`` contains a set of functions for reordering data, and generating hierarchical clustering
for ordered data, for improved visualization.

See the [R package](https://CRAN.R-project.org/package=slanter) for details. Specifically, the
meristems vignette explains why and how to use this package.

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
[R package](https://CRAN.R-project.org/package=slanter) reference manual for a description of all
the available functions.
