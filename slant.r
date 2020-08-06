#' Compute rows and columns orders which move high values close to the diagonal.
#'
#' For a matrix expressing the cross-similarity between two (possibly different)
#' sets of entities, this produces better results than clustering (e.g. as done
#' by \code{pheatmap}). This is because clustering does not care about the order
#' of each two sub-partitions. That is, clustering is as happy with \code{((2, 1),
#' (4, 3))} as it is with the more sensible \code{((1, 2), (3, 4))}. As a result,
#' visualizations of similarities using naive clustering can be misleading.
#'
#' @param data A rectangular matrix.
#' @param order_rows Whether to reorder the rows.
#' @param order_cols Whether to reorder the columns.
#' @param same_order Whether to apply the same order to both rows and columns.
#' @param max_spin_count How many times to retry improving the solution before giving up.
#' @return A list with two keys, \code{rows} and \code{cols}, which contain the order.
#'
#' @export
slanted_orders <- function(data, ..., order_rows=T, order_cols=T,
                           same_order=F, max_spin_count=10) {

    wrapr::stop_if_dot_args(substitute(list(...)), 'slanted_orders')
    rows_count <- dim(data)[1]
    cols_count <- dim(data)[2]

    rows_indices <- as.vector(1:rows_count)
    cols_indices <- as.vector(1:cols_count)

    rows_permutation <- rows_indices
    cols_permutation <- cols_indices

    if (same_order) {
        stopifnot(rows_count == cols_count)
        permutation <- rows_indices
    }

    if (order_rows || order_cols) {
        squared_data <- data * data
        epsilon <- min(squared_data[squared_data > 0]) / 10

        reorder_phase <- function() {
            spinning_rows_count <- 0
            spinning_cols_count <- 0
            spinning_same_count <- 0
            was_changed <- T
            error_rows <- NULL
            error_cols <- NULL
            error_same <- NULL
            while (was_changed) {
                was_changed <- F
                ideal_index <- NULL
                if (order_rows) {
                    sum_indexed_rows <- rowSums(sweep(squared_data, 2, cols_indices, `*`))
                    sum_squared_rows <- rowSums(squared_data)
                    ideal_row_index <- (sum_indexed_rows + epsilon) / (sum_squared_rows + epsilon)

                    if (same_order) {
                        ideal_index <- ideal_row_index
                    } else {
                        error <- ideal_row_index - rows_indices
                        new_error_rows <- sum(error * error)
                        new_rows_permutation <- order(ideal_row_index)
                        new_changed <- any(new_rows_permutation != rows_indices)
                        if (is.null(error_rows) || new_error_rows < error_rows) {
                            error_rows <- new_error_rows
                            spinning_rows_count <- 0
                        } else {
                            spinning_rows_count <- spinning_rows_count + 1
                        }
                        if (new_changed && spinning_rows_count < max_spin_count) {
                            was_changed <- T
                            squared_data <<- squared_data[new_rows_permutation,]
                            rows_permutation <<- rows_permutation[new_rows_permutation]
                        }
                    }
                }

                if (order_cols) {
                    sum_indexed_cols <- colSums(sweep(squared_data, 1, rows_indices, `*`))
                    sum_squared_cols <- colSums(squared_data)
                    ideal_col_index <- (sum_indexed_cols + epsilon) / (sum_squared_cols + epsilon)

                    if (same_order) {
                        if (!is.null(ideal_index)) {
                            ideal_index <- (ideal_index + ideal_col_index) / 2
                        } else {
                            ideal_index <- ideal_col_index
                        }
                    } else {
                        error <- ideal_col_index - cols_indices
                        new_error_cols <- sum(error * error)
                        new_cols_permutation <- order(ideal_col_index)
                        new_changed <- any(new_cols_permutation != cols_indices)
                        if (is.null(error_cols) || new_error_cols < error_cols) {
                            error_cols <- new_error_cols
                            spinning_cols_count <- 0
                        } else {
                            spinning_cols_count <- spinning_cols_count + 1
                        }
                        if (new_changed && spinning_cols_count < max_spin_count) {
                            was_changed <- T
                            squared_data <<- squared_data[,new_cols_permutation]
                            cols_permutation <<- cols_permutation[new_cols_permutation]
                        }
                    }
                }

                if (!is.null(ideal_index)) {
                    error <- ideal_index - rows_indices
                    new_error_same <- sum(error * error)
                    new_permutation <- order(ideal_index)
                    new_changed <- any(new_permutation != rows_indices)
                    if (is.null(error_same) || new_error_same < error_same) {
                        error_same <- new_error_same
                        spinning_same_count <- 0
                    } else {
                        spinning_same_count <- spinning_same_count + 1
                    }
                    if (new_changed && spinning_same_count < max_spin_count) {
                        was_changed <- T
                        squared_data <<- squared_data[new_permutation,new_permutation]
                        permutation <<- permutation[new_permutation]
                        rows_permutation <<- permutation
                        cols_permutation <<- permutation
                    }
                }
            }
        }

        discount_outliers <- function() {
            row_indices_matrix <- matrix(rep(rows_indices, each=cols_count),
                                         nrow=rows_count, ncol=cols_count, byrow=T)
            col_indices_matrix <- matrix(rep(cols_indices, each=rows_count),
                                         nrow=rows_count, ncol=cols_count, byrow=F)

            rows_per_col <- rows_count / cols_count
            cols_per_row <- cols_count / rows_count

            ideal_row_indices_matrix <- col_indices_matrix * rows_per_col
            ideal_col_indices_matrix <- row_indices_matrix * cols_per_row

            row_distance_matrix <- row_indices_matrix - ideal_row_indices_matrix
            col_distance_matrix <- col_indices_matrix - ideal_col_indices_matrix

            weight_matrix <- (1 + abs(row_distance_matrix)) * (1 + abs(col_distance_matrix))
            squared_data <<- squared_data / weight_matrix
        }

        reorder_phase()
        discount_outliers()
        reorder_phase()
    }

    return (list(rows=rows_permutation, cols=cols_permutation))
}

#' Reorder data rows and columns to move high values close to the diagonal.
#'
#' Given a matrix expressing the cross-similarity between two (possibly different)
#' sets of entities, this uses \code{slanted_orders} to compute the "best" order
#' for visualizing the matrix, then returns the reordered data. Commonly used in:
#' \code{pheatmap(slanted_reorder(data),cluster_rows=F,cluster_cols=F)}, and of course
#' by \code{sheatmap}.
#'
#' @param data A rectangular matrix.
#' @param order_rows Whether to reorder the rows.
#' @param order_cols Whether to reorder the columns.
#' @param same_order Whether to apply the same order to both rows and columns.
#' @return A matrix of the same shape whose rows and columns are a permutation of the input.
#'
#' @export
slanted_reorder <- function(data, ..., order_rows=T, order_cols=T, same_order=F) {
    wrapr::stop_if_dot_args(substitute(list(...)), 'slanted_reorder')
    orders <- slanted_orders(data,
                             order=order, order_rows=order_rows, order_cols=order_cols,
                             same_order=same_order)
    return (data[orders$rows, orders$cols])
}

#' Accelerated version of \code{cor} for matrices.
accelerated_cor <- function(data, ..., method='pearson') {
    wrapr::stop_if_dot_args(substitute(list(...)), 'accelerated_cor')
    if (exists('tgs_cor') && method %in% c('pearson', 'spearman')) {
        return (tgs_cor(data, spearman=(method == 'spearman')))
    }
    return (cor(data, method=method))
}

#' Plot a heatmap with values as close to the diagonal as possible.
#'
#' Given a matrix expressing the cross-similarity between two (possibly
#' different) sets of entities, this will reorder it to move the high values
#' close to the diagonal, for better visualizations.
#'
#' If you want to cluster your data, you have two options:
#'
#' If you do not specify a clustering, and just request one (e.g.,
#' ``cluster_rows=T``), then \code{sheatmap} will invoke \code{oclust} to
#' generate the "best" clustering that is also compatible with the slanted order.
#' This is the default.
#'
#' Otherwise, \code{sheatmap} will obey any clustering it was given (e.g.,
#' ``cluster_rows=hclust(...)``. Note that any given #' clustering tree allows
#' for multiple orders, since in each node, the two sub-trees can be flipped,
#' without modifying the tree. Thus, \code{sheatmap} will use the
#' "best" order that is compatible with the given clustering.
#'
#' In addition, you can give this function any of the \code{pheatmap} flags,
#' and it will just pass them on. This allows full control over the diagram's
#' features.
#'
#' @param data A rectangular matrix
#' @param annotation_row Optional data frame describing each row.
#' @param annotation_col Optional data frame describing each column.
#' @param order_rows Whether to reorder the rows.
#' @param order_cols Whether to reorder the columns.
#' @param same_order Whether to apply the same order to both rows and columns.
#' @param cluster_rows Whether to cluster the rows, or the clustering to use.
#' @param cluster_cols Whether to cluster the columns, or the clustering to use.
#' @param clustering_distance_cols The default method for computing column distances (by default, \code{euclidian}).
#' @param clustering_distance_rows The default method for computing row distances (by default, \code{euclidian}).
#' @param clustering_method The default method to use for clustering (by default, \code{ward.D2}).
#' @param clustering_callback Is not supported.
#' @param ... Additional flags to pass to \code{pheatmap}.
#' @return Whatever \code{pheatmap} returns.
#'
#' @export
sheatmap <- function(data, ...,
                     annotation_col=NULL,
                     annotation_row=NULL,
                     order_rows=T,
                     order_cols=T,
                     same_order=F,
                     cluster_rows=T,
                     cluster_cols=T,
                     clustering_distance_rows='euclidian',
                     clustering_distance_cols='euclidian',
                     clustering_method='ward.D2',
                     clustering_callback=NA) {
    stopifnot(is.na(clustering_callback))  # Not implemented.

    ideal_orders <-
        slanted_orders(data, order_rows=order_rows, order_cols=order_cols, same_order=same_order)

    if (class(cluster_rows) == 'hclust') {
        cluster_rows <- best_reordered_hclust(cluster_rows, ideal_orders$rows)
        rows_order <- cluster_rows$order
        cluster_rows <- pre_ordered_hclust(cluster_rows)
    } else if (cluster_rows) {
        rows_distances <- dist(data, method=clustering_distance_rows)
        rows_order = ideal_orders$row
        cluster_rows <- oclust(rows_distances, order=rows_order, method=clustering_method)
    } else {
        rows_order <- 1:dim(data)[1]
    }

    if (class(cluster_cols) == 'hclust') {
        cluster_cols <- best_reordered_hclust(cluster_cols, ideal_orders$cols)
        cols_order <- cluster_cols$order
        cluster_cols <- pre_ordered_hclust(cluster_cols)
    } else if (cluster_cols) {
        cols_distances <- dist(data, method=clustering_distance_cols)
        cols_order = ideal_orders$col
        cluster_cols <- oclust(cols_distances, order=cols_order, method=clustering_method)
    } else {
        cols_order <- 1:dim(data)[2]
    }

    data <- data[rows_order, cols_order]

    if (!is.null(annotation_row)) {
        annotation_row <- reorder_frame(annotation_row, orders$rows)
    }
    if (!is.null(annotation_col)) {
        annotation_col <- reorder_frame(annotation_col, orders$cols)
    }

    return (pheatmap::pheatmap(data, annotation_row=annotation_row, annotation_col=annotation_col,
                               cluster_rows=cluster_rows, cluster_cols=cluster_cols, ...))
}

#' Reorder the rows of a frame.
#'
#' If you expect \code{data[order,]} to "just work", you haven't been using R
#' for very long. It is this sort of thing that makes me *hate* coding in R.
#'
#' @param frame A data frame to reorder the rows of.
#' @param order An array containing indices permutation to apply to the rows.
#' @return The data frame with the new row orders.
#'
#' @export
reorder_frame <- function(data, order) {
    row_names <- rownames(data)
    if (ncol(data) == 1) {
        vec <- t(data[1])
        data[1] <- vec[order]
    } else {
        data <- data[order,]
    }
    rownames(data) <- row_names[order]
    return (data)
}

#' Given a clustering of some data, and some ideal order we'd like to use to
#' visualize it, reorder (but do not modify) the clustering to be as consistent
#' as possible with this ideal order.
#'
#' @param cluster The existing clustering of the data.
#' @param ideal_order The order we'd like to see the data in.
#' @return A reordered clustering which obeys, wherever possible, the ideal order.
#'
#' @export
best_reordered_hclust <- function(clusters, ideal_order) {
    old_of_mid <- clusters$order
    mid_of_old <- Matrix::invPerm(old_of_mid)

    merge_of_old <- clusters$merge
    merge_of_mid <- permute_merge(merge_of_old, mid_of_old)

    tree_of_mid <- tree_of_merge(merge_of_mid)

    ideal_old_of_new <- ideal_order
    ideal_new_of_old <- Matrix::invPerm(ideal_old_of_new)
    ideal_new_of_mid <- ideal_new_of_old[old_of_mid]

    mid_of_new <- best_tree_compatible_order(tree_of_mid, ideal_new_of_mid)
    old_of_new <- old_of_mid[mid_of_new]

    clusters$order <- old_of_new

    return (clusters)
}

# Given a clustering which specified some data order, given we reorder the
# data ourselves, return a clustering that applies to the reordered data.
pre_ordered_hclust <- function(clusters) {
    old_of_new <- clusters$order
    new_of_old <- Matrix::invPerm(old_of_new)
    clusters$merge <- permute_merge(clusters$merge, new_of_old)
    clusters$order <- 1:length(clusters$order)
    return (clusters)
}

# Given an hclust merge array, return a new one that applies the same clustering
# if the data has been reordered.
permute_merge <- function(merge, new_of_old) {
    merges_count <- dim(merge)[1]
    for (merge_index in 1:merges_count) {
        for (entry_index in 1:2) {
            if (merge[merge_index, entry_index] < 0) {
                merge[merge_index, entry_index] <- -new_of_old[-merge[merge_index, entry_index]]
            }
        }
    }

    return (merge)
}

# Given an hclust merge array, return a different tree representation of the clustering.
tree_of_merge <- function(merge) {
    low_indices <- c()
    low_nodes <- c()
    mid_indices <- c()
    high_nodes <- c()
    high_indices <- c()
    flips <- c()

    merges_count <- dim(merge)[1]
    for (merge_index in 1:merges_count) {
        low <- merge[merge_index, 1]
        high <- merge[merge_index, 2]

        if (low < 0) {
            low_node <- NA
            low_index <- -low
            low_mid_index <- -low + 1
        } else {
            low_node <- low
            low_index <- low_indices[low_node]
            low_mid_index <- high_indices[low_node] + 1
        }

        if (high < 0) {
            high_node <- NA
            high_mid_index <- -high
            high_index <- -high
        } else {
            high_node <- high
            high_mid_index <- low_indices[high_node]
            high_index <- high_indices[high_node]
        }

        stopifnot(low_mid_index == high_mid_index)

        low_indices <- append(low_indices, low_index)
        low_nodes <- append(low_nodes, low_node)
        mid_indices <- append(mid_indices, low_mid_index)
        high_nodes <- append(high_nodes, high_node)
        high_indices <- append(high_indices, high_index)
        flips <- append(flips, F)
    }

    return (list(low_indices=low_indices,
                 low_nodes=low_nodes,
                 mid_indices=mid_indices,
                 high_nodes=high_nodes,
                 high_indices=high_indices,
                 flips=flips))
}

# Given a tree (computed from clusters merge array), and a target order for the nodes, compute the
# order that is closest to the ideal, but still compatible with the tree.
best_tree_compatible_order <- function(tree, ideal_new_of_old) {
    count <- length(ideal_new_of_old)

    indices <- as.vector(1:count)

    new_of_old <- indices

    did_change <- T

    reorder_node <- function(root_node) {
        low_node <- tree$low_nodes[root_node]
        if (!is.na(low_node)) {
            reorder_node(low_node)
        }

        high_node <- tree$high_nodes[root_node]
        if (!is.na(high_node)) {
            reorder_node(high_node)
        }

        old_low_indices <- tree$low_indices[root_node] : (tree$mid_indices[root_node] - 1)
        old_high_indices <- tree$mid_indices[root_node] : tree$high_indices[root_node]
        old_full_indices <- tree$low_indices[root_node] : tree$high_indices[root_node]

        new_low_indices <- new_of_old[old_low_indices]
        new_high_indices <- new_of_old[old_high_indices]

        ideal_new_low_indices <- ideal_new_of_old[old_low_indices]
        ideal_new_high_indices <- ideal_new_of_old[old_high_indices]

        curr_low_diff <- new_low_indices - ideal_new_low_indices
        curr_high_diff <- new_high_indices - ideal_new_high_indices

        if (tree$flips[root_node]) {
            low_delta <- -length(old_high_indices)
            high_delta <- length(old_low_indices)
        } else {
            low_delta <- length(old_high_indices)
            high_delta <- -length(old_low_indices)
        }

        flip_low_diff <- curr_low_diff + low_delta
        flip_high_diff <- curr_high_diff + high_delta

        curr_error <- sum(curr_low_diff * curr_low_diff) + sum(curr_high_diff * curr_high_diff)
        flip_error <- sum(flip_low_diff * flip_low_diff) + sum(flip_high_diff * flip_high_diff)

        if (flip_error < curr_error) {
            did_change <<- T
            tree$flips[root_node] <<- !tree$flips[root_node]

            new_of_old[old_low_indices] <<- new_of_old[old_low_indices] + low_delta
            new_of_old[old_high_indices] <<- new_of_old[old_high_indices] + high_delta
        }
    }

    while (did_change) {
        did_change <- F
        reorder_node(length(tree$mid_indices))
    }

    old_of_new <- Matrix::invPerm(new_of_old)

    return (old_of_new)
}

#' Hierarchically cluster ordered data.
#'
#' Given a distance matrix for sorted objects, compute a hierarchical clustering preserving this
#' order. That is, this is similar to \code{hclust} with the constraint that the result's order is
#' always \code{1:N}.
#'
#' If an \code{order} is specified, assumes that the data will be re-ordered by this order.
#' That is, the indices in the returned \code{hclust} object will refer to the post-reorder
#' data locations, **not** to the current data locations.
#'
#" Currently, the only methods supported are \code{'ward.D'} and \code{'ward.D2'}.
#'
#' This can be applied to the results of \code{slanted_reorder\code}, to give
#' a "plausible" clustering for the data.
#'
#' @param dist A distances object (as created by \code{dist}).
#' @param method The clustering method to use (only \code{'ward.D'} and \code{'ward.D2'} are supported}).
#' @param order If specified, assume the data will be re-ordered by this order.
#' @param members Optionally, the number of members in each element of the distances (by default, one each).
#' @return A clustering object (as created by \code{hclust}).
#'
#' @export
oclust <- function(dist, method='ward.D2', order=NULL, members=NULL) {
    distances <- as.matrix(dist)
    stopifnot(dim(distances)[1] == dim(distances)[2])
    entities_count <- dim(distances)[1]

    if (method == 'ward.D2') {
        distances <- distances * distances
        sqrt_height <- T
    } else {
        stopifnot(method %in% c('ward.D', 'ward.D2'))
        sqrt_height <- F
    }

    if (!is.null(order)) {
        distances <- distances[order, order]
    }

    diag(distances) <- Inf

    merge <- matrix(0, nrow=entities_count - 1, ncol=2)
    height <- rep(0, entities_count - 1)
    merged_height <- rep(0, entities_count)
    groups <- -(1:entities_count)
    if (is.null(members)) {
        members <- rep(1, entities_count)
    }

    for (merge_index in 1:(entities_count - 1)) {
        adjacent_distances <- pracma::Diag(distances, 1)

        low_index <- which.min(adjacent_distances)
        high_index <- low_index + 1

        grouped_indices <- groups[c(low_index, high_index)]

        merged_indices <- which(groups %in% grouped_indices)

        groups[merged_indices] <- merge_index

        merge[merge_index,] <- grouped_indices

        delta_height <- adjacent_distances[low_index]
        if (sqrt_height) {
            delta_height <- sqrt(delta_height)
        }
        height[merge_index] <- max(merged_height[merged_indices]) + delta_height

        merged_height[merged_indices] <- height[merge_index]

        a_index <- merged_indices[1]
        b_index <- merged_indices[length(merged_indices)]

        a_members <- members[a_index]
        b_members <- members[b_index]

        members[merged_indices] <- a_members + b_members

        a_b_distance_value <- distances[a_index, b_index]  # d(a, b)
        a_b_distance_scaled <- members * a_b_distance_value  # |C| * d(a, b)

        a_c_distance_slice <- distances[a_index, ] # d(a, c)
        a_c_scale <- rep(a_members, entities_count) + members # |A| + |C|
        a_c_distance_scaled <- a_c_distance_slice * a_c_scale # (|A| + |C|) * d(a, c)

        b_c_distance_slice <- distances[b_index, ] # d(b, c)
        b_c_scale <- rep(b_members, entities_count) + members # |B| + |C|
        b_c_distance_scaled <- b_c_distance_slice * b_c_scale # (|B| + |C|) * d(b, c)

        a_b_c_scale <- members + a_members + b_members  # |A| + |B| + |C|

        # Ward: ( (|A| + |C|) * d(a,c) + (|B| + |C|) * d(b, c) - |C| * d(a, b) ) / ( |A| + |B| + |C| )
        merged_distance <-
            (a_c_distance_scaled + b_c_distance_scaled - a_b_distance_scaled) / a_b_c_scale

        distances[,merged_indices] <- merged_distance
        distances[merged_indices,] <- rep(merged_distance, each=length(merged_indices))
    }

    hclust <- list(merge=merge, height=height)

    hclust$method <- 'oclust'
    hclust$order <- 1:entities_count
    class(hclust) <- 'hclust'

    return (hclust)
}
