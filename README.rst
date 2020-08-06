Slant.R
=======

Installation and Usage
----------------------

Download ``slant.r``. Inside R, ``source('slant.r')``.

Use ``sheatmap`` as a drop-in replacement to ``pheatmap``. In general, any argument to ``sheatmap``
(except ``clustering_callback``) is forwarded to ``pheatmap``, so all the diagram aspects can be
controlled in the usual way.

Lower-level functions are also available; see the `code documentation <slant.r>`_ for details.

The Problem
===========

Given a matrix of some similarity data (that is, where data[i,j] is higher if the entity represented
by row i is "more similar" to the entity represented by column j), how best to visualize it?

The default R solution is to use ``pheatmap``. This has the advantage it also reorders the data so
that most-similar rows (or columns) are near each other. This is very useful for visualizing the
clustering of the entries into some groups.

However, while the resulting clustering is "optimal" (for some clustering algorithm), the resulting
visualization is not necessarily the best one. Intuitively, all ``pheatmap`` cares about (if
anything) is the clustering tree, but not its order. In one dimensions, consider clustering the
numbers 1 through 4. The ``pheatmap`` is as happy with the tree ``(4,3), (1,2)`` as it is with the
more legible tree ``(1,2), (3,4)``. In two dimensions, this results in illegible, and at times
downright misleading, heatmaps.

It is best to demonstrate this with an example. All the example code below assumes we first:

.. code:: r

    read('data.Rda')
    source('slant.r')

Then, writing:

.. code:: r

    pheatmap::pheatmap(data, show_rownames=F, show_colnames=F)

Gives us:

.. image:: large_clustered_pheatmap.png

At a first glance, this seems to indicate there are two very sharply distinguished sub-groups, with
a third mixed group. A closer look (below) will reveal a different picture.

Slanted Matrices
----------------

A "slanted matrix" is a matrix which is reordered such that its highest values are as close to the
diagonal as possible. In this package, this is implemented by the ``slanted_orders`` function, which
repeatedly sorts the rows by the index of their "diagonal" column, and then the columns by the index
of their "diagonal" rows, until the system stabilizes (this is pretty fast).

.. note::

    * This works for any size similarity matrix, not necessarily the square symmetric matrix one
      gets from computing correlations.

    * The code actually works on the square of the data. This works OK-ish if the similarity data
      is a correlation and one takes -1 (perfect negative correlation) to be as strong a similarity
      indicator as +1 (perfect correlation). If your data isn't like that, just make sure it is all
      non-negative, where larger values are "more similar".

The ``sheatmap`` function wraps this (with additional functionality described below). We'll begin
with ignoring the clustering, writing:

.. code:: r

    sheatmap(data, show_rownames=F, show_colnames=F, cluster_rows=F, cluster_cols=F)

Which gives us:

.. image:: large_unclustered_sheatmap.png

**This is exactly the same data as before.**

Note that if we do not perform any clustering, ``pheatmap`` will preserve the original data order.
Sometimes this is exactly what you want, but when the data has no a-priori order, writing:

.. code:: r

    pheatmap::pheatmap(data, show_rownames=F, show_colnames=F, cluster_rows=F, cluster_cols=F)

Gives us:

.. image:: large_unclustered_pheatmap.png

Which most definitely is *not* what we want to see.

Slanted Clustering
------------------

What if we still want to see the cluster structure of our data? The ``sheatmap``
function provides three options to do this, with different trade-offs.

Reordering
..........

By default (or if you specify ``clusters='reorder'``), if clustering is done,
then ``sheatmap`` will preserve it. However, it will use the fact that clustering
does not specify (total) order, to reorder the clustering tree to best fit the
slanted order. For example, writing:

.. code:: r

    sheatmap(data, show_rownames=F, show_colnames=F)

Gives us:

.. image:: large_reordered_sheatmap.png

Here the clustering tree can make the claim it is "optimal" in some sense - it is identical to the
original ``pheatmap`` clustering tree. The order is the best we can have under this constraint, so
we do not get the full benefits of slanting, but we get at least a hint that there is a smooth
gradient between the groups rather than a sharp divide.

Modifying
.........

A stronger option is to allow ``sheatmap`` to modify the clustering tree. By specifying
``clusters='modify'``, it maintains the "ideal" slanted order, and generate a compatible clustering
tree which is the "closest possible" to the original unconstrained clustering. For example, writing:

.. code:: r

    sheatmap(data, show_rownames=F, show_colnames=F, clusters='modify')

Will give us:

.. image:: large_modified_sheatmap.png

We can see that while in theory the approach seems promising, in practice it tends to generate
low-quality clustering trees. However, YMMV.

Replacing
.........

Finally, ``sheatmap`` allows to simply build a brand-new clustering tree, tailored to the
"ideal" slanted order. By specifying ``clusters='replace'``, the code will discard the
original clustering (if any), and will invoke ``oclust`` (see below) to get a new one.
For example, writing:

.. code:: r

    sheatmap(data, show_rownames=F, show_colnames=F, clusters='replace')

Will give us:

.. image:: large_replaced_sheatmap.png

Here we see not only the clear gradient, but also that it splits naturally to four parts
(that smoothly transition from one to the next). We can highlight this by using ``cutree``:

.. code:: r

    sheatmap(data, show_rownames=F, show_colnames=F, clusters='replace', cutree_rows=4, cutree_cols=4)

Which gives us:

.. image:: large_cut_replaced_sheatmap.png

Which is pretty good, even though it doesn't have the claim of "full optimality" that the original
unrestricted clustering offers; instead it only promises "constrained optimality" of the clustering
(subject to the slanting order).

Ordered Clustering
------------------

In general hierarchical clustering of entities (rows or columns in our case) tries to create a
binary tree such that the more similar two entities are, the closer they are in the tree. The twist
in ordered clustering (as implemented by ``oclust``) is that the entities are ordered E = { e_1, ...
e_N } and each group of entities clustered under any tree node must be a contiguous range of
entities G = { e_start, ..., e_stop }.

Visually this means that if you draw the tree on top of the ordered entities, there would be no edge
crossings. This makes ordered clustering a natural addition to the slanted matrix visualization.
Normally, an unconstrained hierarchical clustering is done first, and the matrix visualization order
is chosen to be compatible with the clustering tree. Here we do the opposite - we first order the
entities for the slanted matrix visualization, and then find a compatible clustering tree to go with
it.

In principle it is possible to adapt any clustering method to include an ordering constraint. We
chose to adapt `Ward's method <https://en.wikipedia.org/wiki/Ward%27s_method>`_. In this method, the
algorithm starts with each element on its own, and merges the two elements such that the total
variance within the merge node is minimized. It recursively merges groups of elements (minimizing
the variance each time) until obtaining the full clustering tree.

The general Ward's method can pick any two element groups to combine at each step. The ordered
``oclust`` variant is restricted to only choosing adjacent element groups, trusting the slanting
order to make similar elements adjacent to each other. This has the secondary effect that it is very
fast, which allowing for a practical pure R implementation. The full Ward's method is slower, so the
R implementation resorts to using FORTRAN code.

TODO
====

Convert this to a proper R package.

LICENSE (MIT)
=============

This code is available under the `MIT license <LICENSE.rst>`_.
