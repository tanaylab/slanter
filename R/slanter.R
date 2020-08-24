#' Compute rows and columns orders which move high values close to the diagonal.
#'
#' For a matrix expressing the cross-similarity between two (possibly different) sets of entities,
#' this produces better results than clustering (e.g. as done by \code{pheatmap}). This is because
#' clustering does not care about the order of each two sub-partitions. That is, clustering is as
#' happy with \code{((2, 1), (4, 3))} as it is with the more sensible \code{((1, 2), (3, 4))}. As a
#' result, visualizations of similarities using naive clustering can be misleading.
#'
#' @param data A rectangular matrix containing non-negative values.
#' @param order_rows Whether to reorder the rows.
#' @param order_cols Whether to reorder the columns.
#' @param squared_order Whether to reorder to minimize the l2 norm (otherwise minimizes the l1 norm).
#' @param same_order Whether to apply the same order to both rows and columns.
#' @param max_spin_count How many times to retry improving the solution before giving up.
#' @return A list with two keys, \code{rows} and \code{cols}, which contain the order.
#'
#' @export
#'
#' @examples
#' slanter::slanted_orders(cor(t(mtcars)))
slanted_orders <- function(data, order_rows=TRUE, order_cols=TRUE,
                           squared_order=TRUE, same_order=FALSE, max_spin_count=10) {

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
        stopifnot(min(data) >= 0)
        if (squared_order) {
            data <- data * data
        }
        epsilon <- min(data[data > 0]) / 10

        reorder_phase <- function() {
            spinning_rows_count <- 0
            spinning_cols_count <- 0
            spinning_same_count <- 0
            was_changed <- TRUE
            error_rows <- NULL
            error_cols <- NULL
            error_same <- NULL
            while (was_changed) {
                was_changed <- FALSE
                ideal_index <- NULL
                if (order_rows) {
                    sum_indexed_rows <- rowSums(sweep(data, 2, cols_indices, `*`))
                    sum_squared_rows <- rowSums(data)
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
                            was_changed <- TRUE
                            data <<- data[new_rows_permutation,]
                            rows_permutation <<- rows_permutation[new_rows_permutation]
                        }
                    }
                }

                if (order_cols) {
                    sum_indexed_cols <- colSums(sweep(data, 1, rows_indices, `*`))
                    sum_squared_cols <- colSums(data)
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
                            was_changed <- TRUE
                            data <<- data[,new_cols_permutation]
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
                        was_changed <- TRUE
                        data <<- data[new_permutation,new_permutation]
                        permutation <<- permutation[new_permutation]
                        rows_permutation <<- permutation
                        cols_permutation <<- permutation
                    }
                }
            }
        }

        discount_outliers <- function() {
            row_indices_matrix <- matrix(rep(rows_indices, each=cols_count),
                                         nrow=rows_count, ncol=cols_count, byrow=TRUE)
            col_indices_matrix <- matrix(rep(cols_indices, each=rows_count),
                                         nrow=rows_count, ncol=cols_count, byrow=FALSE)

            rows_per_col <- rows_count / cols_count
            cols_per_row <- cols_count / rows_count

            ideal_row_indices_matrix <- col_indices_matrix * rows_per_col
            ideal_col_indices_matrix <- row_indices_matrix * cols_per_row

            row_distance_matrix <- row_indices_matrix - ideal_row_indices_matrix
            col_distance_matrix <- col_indices_matrix - ideal_col_indices_matrix

            weight_matrix <- (1 + abs(row_distance_matrix)) * (1 + abs(col_distance_matrix))
            data <<- data / weight_matrix
        }

        reorder_phase()
        discount_outliers()
        reorder_phase()
    }

    return (list(rows=rows_permutation, cols=cols_permutation))
}

#' Reorder data rows and columns to move high values close to the diagonal.
#'
#' Given a matrix expressing the cross-similarity between two (possibly different) sets of entities,
#' this uses \code{slanted_orders} to compute the "best" order for visualizing the matrix, then
#' returns the reordered data. Commonly used in \code{pheatmap(slanted_reorder(data), ...)}, and of
#' course \code{sheatmap} does this internally for you.
#'
#' @param data A rectangular matrix to reorder, of non-negative values (unless \code{order_data} is specified).
#' @param order_data An optional matrix of non-negative values of the same size to use for computing the orders.
#' @param order_rows Whether to reorder the rows.
#' @param order_cols Whether to reorder the columns.
#' @param squared_order Whether to reorder to minimize the l2 norm (otherwise minimizes the l1 norm).
#' @param same_order Whether to apply the same order to both rows and columns.
#' @return A matrix of the same shape whose rows and columns are a permutation of the input.
#'
#' @export
#'
#' @examples
#' slanter::slanted_reorder(cor(t(mtcars)))
slanted_reorder <- function(data, order_data=NULL, order_rows=TRUE, order_cols=TRUE,
                            squared_order=TRUE, same_order=FALSE) {
    if (is.null(order_data)) {
        order_data <- data
    }
    stopifnot(all(dim(order_data) == dim(data)))

    orders <- slanted_orders(order_data,
                             squared_order=squared_order, same_order=same_order,
                             order_rows=order_rows, order_cols=order_cols)

    return (data[orders$rows, orders$cols])
}

#' Plot a heatmap with values as close to the diagonal as possible.
#'
#' Given a matrix expressing the cross-similarity between two (possibly different) sets of entities,
#' this will reorder it to move the high values close to the diagonal, for a better visualization.
#'
#' If you have an a-priori order for the rows and/or columns, you can prevent reordering either or
#' both by specifying \code{order_rows=FALSE} and/or \code{order_cols=FALSE}. Otherwise,
#' \code{slanted_orders} is used to compute the "ideal" slanted order for the data.
#'
#' By default, the rows and columns are ordered independently from each other. If the matrix is
#' asymmetric but square (e.g., a matrix of weights of a directed graph such as a
#' K-nearest-neighbors graph), then you can can specify \code{same_order=TRUE} to force both rows and
#' columns to the same order.
#'
#' There are four options for controlling clustering:
#'
#' * By default, \code{sheatmap} will generate a clustering tree using \code{oclust}, to generate
#'   the "best" clustering that is also compatible with the slanted order.
#'
#' * Request that \code{sheatmap} will use the same \code{hclust} as
#'   \code{pheatmap} (e.g., \code{oclust_rows=FALSE}). In this case, the tree is reordered to
#'   be the "most compatible" with the target slanted order. That is, \code{sheatmap} will invoke
#'   \code{reorder_hclust} so that, for each node of the tree, the order of the two sub-trees will
#'   be chosen to best match the target slanted order. The end result need not be identical to the
#'   slanted order, but is as close as possible given the \code{hclust} clustering tree.
#'
#' * Specify an explicit clustering (e.g., \code{cluster_rows=hclust(...)}). In this case,
#'   \code{sheatmap} will again merely reorder the tree but will not modify it.
#'
#" * Disable clustering altogether (e.g., \code{cluster_rows=FALSE}).
#'
#' In addition, you can give this function any of the \code{pheatmap} flags, and it will just pass
#' them on. This allows full control over the diagram's features.
#'
#' Note that \code{clustering_callback} is not supported. In addition, the default
#' \code{clustering_method} here is \code{ward.D2} instead of \code{complete}, since the only
#' methods supported by \code{oclust} are \code{ward.D} and \code{ward.D2}.
#'
#' @param data A rectangular matrix to plot, of non-negative values (unless \code{order_data} is specified).
#' @param order_data An optional matrix of non-negative values of the same size to use for computing the orders.
#' @param annotation_row Optional data frame describing each row.
#' @param annotation_col Optional data frame describing each column.
#' @param order_rows Whether to reorder the rows. Otherwise, use the current order.
#' @param order_cols Whether to reorder the columns. Otherwise, use the current order.
#' @param squared_order Whether to reorder to minimize the l2 norm (otherwise minimizes the l1 norm).
#' @param same_order Whether to apply the same order to both rows and columns (if reordering both).
#' @param cluster_rows Whether to cluster the rows, or the clustering to use.
#' @param cluster_cols Whether to cluster the columns, or the clustering to use.
#' @param oclust_cols Whether to use \code{oclust} instead of \code{hclust} for the columns (if
#'        clustering them).
#' @param oclust_rows Whether to use \code{oclust} instead of \code{hclust} for the rows (if
#'        clustering them).
#' @param clustering_distance_cols The default method for computing column distances (by default,
#'        \code{euclidian}).
#' @param clustering_distance_rows The default method for computing row distances (by default,
#'        \code{euclidian}).
#' @param clustering_method The default method to use for hierarchical clustering (by default,
#'        \code{ward.D2} and *not* \code{complete}).
#' @param clustering_callback Is not supported.
#' @param ... Additional flags to pass to \code{pheatmap}.
#' @return Whatever \code{pheatmap} returns.
#'
#' @export
#'
#' @examples
#' slanter::sheatmap(cor(t(mtcars)))
#' slanter::sheatmap(cor(t(mtcars)), oclust_rows=FALSE, oclust_cols=FALSE)
#' pheatmap::pheatmap(cor(t(mtcars)))
sheatmap <- function(data, ...,
                     order_data=NULL,
                     annotation_col=NULL,
                     annotation_row=NULL,
                     order_rows=TRUE,
                     order_cols=TRUE,
                     squared_order=TRUE,
                     same_order=FALSE,
                     cluster_rows=TRUE,
                     cluster_cols=TRUE,
                     oclust_rows=TRUE,
                     oclust_cols=TRUE,
                     clustering_distance_rows='euclidian',
                     clustering_distance_cols='euclidian',
                     clustering_method='ward.D2',
                     clustering_callback=NA) {
    stopifnot(is.na(clustering_callback))  # Not implemented.
    stopifnot(clustering_method %in% c('ward.D', 'ward.D2'))

    if (is.null(order_data)) {
        order_data <- data
    }
    stopifnot(all(dim(order_data) == dim(data)))

    ideal_orders <-
        slanted_orders(order_data, order_rows=order_rows, order_cols=order_cols,
                       same_order=same_order)

    rows_order <- NULL

    if (class(cluster_rows) == 'logical' && cluster_rows) {
        rows_distances <- stats::dist(data, method=clustering_distance_rows)
        if (oclust_rows) {
            rows_order <- ideal_orders$row
            cluster_rows <- oclust(rows_distances, order=rows_order, method=clustering_method)
        } else {
            cluster_rows <- stats::hclust(rows_distances, method=clustering_method)
        }
    }

    if (is.null(rows_order)) {
        if (class(cluster_rows) == 'hclust') {
            cluster_rows <- reorder_hclust(cluster_rows, ideal_orders$rows)
            rows_order <- cluster_rows$order
            cluster_rows <- pre_ordered_hclust(cluster_rows)
        } else {
            rows_order <- ideal_orders$row
        }
    }

    cols_order <- NULL

    if (class(cluster_cols) == 'logical' && cluster_cols) {
        cols_distances <- stats::dist(data, method=clustering_distance_cols)
        if (oclust_cols) {
            cols_order <- ideal_orders$col
            cluster_cols <- oclust(cols_distances, order=cols_order, method=clustering_method)
        } else {
            cluster_cols <- stats::hclust(cols_distances, method=clustering_method)
        }
    }

    if (is.null(cols_order)) {
        if (class(cluster_cols) == 'hclust') {
            cluster_cols <- reorder_hclust(cluster_cols, ideal_orders$cols)
            cols_order <- cluster_cols$order
            cluster_cols <- pre_ordered_hclust(cluster_cols)
        } else {
            cols_order <- ideal_orders$col
        }
    }

    data <- data[rows_order, cols_order]

    if (!is.null(annotation_row)) {
        annotation_row <- reorder_frame(annotation_row, rows_order)
    }
    if (!is.null(annotation_col)) {
        annotation_col <- reorder_frame(annotation_col, cols_order)
    }

    return (pheatmap::pheatmap(data, annotation_row=annotation_row, annotation_col=annotation_col,
                               cluster_rows=cluster_rows, cluster_cols=cluster_cols, ...))
}

#' Reorder the rows of a frame.
#'
#' You'd expect \code{data[order,]} to "just work". It doesn't for data frames with a single column,
#' which happens for annotation data, hence the need for this function. Sigh.
#'
#' @param frame A data frame to reorder the rows of.
#' @param order An array containing indices permutation to apply to the rows.
#' @return The data frame with the new row orders.
#'
#' @export
#'
#' @examples
#' df <- data.frame(foo=c(1, 2, 3))
#' df[c(1,3,2),]
#' slanter::reorder_frame(df, c(1,3,2))
reorder_frame <- function(frame, order) {
    row_names <- rownames(frame)
    if (ncol(frame) == 1) {
        vec <- t(frame[1])
        frame[1] <- vec[order]
    } else {
        frame <- frame[order,]
    }
    rownames(frame) <- row_names[order]
    return (frame)
}

#' Given a clustering of some data, and some ideal order we'd like to use to visualize it, reorder
#' (but do not modify) the clustering to be as consistent as possible with this ideal order.
#'
#' @param clusters The existing clustering of the data.
#' @param order The ideal order we'd like to see the data in.
#' @return A reordered clustering which is consistent, wherever possible, the ideal order.
#'
#' @export
#'
#' @examples
#' clusters <- hclust(dist(mtcars))
#' clusters$order
#' clusters <- slanter::reorder_hclust(clusters, 1:length(clusters$order))
#' clusters$order
reorder_hclust <- function(clusters, order) {
    old_of_new <- order
    new_of_old <- Matrix::invPerm(old_of_new)

    merge <- clusters$merge
    merges_count <- dim(merge)[1]
    merge_data <- array(list(), merges_count)

    for (merge_index in 1:merges_count) {
        a_index <- merge[merge_index, 1]
        b_index <- merge[merge_index, 2]

        if (a_index < 0) {
            a_indices <- c(-a_index)
            a_center <- new_of_old[-a_index]
        } else {
            a_data <- merge_data[[a_index]]
            a_indices <- a_data$indices
            a_center <- a_data$center
        }

        if (b_index < 0) {
            b_indices <- c(-b_index)
            b_center <- new_of_old[-b_index]
        } else {
            b_data <- merge_data[[b_index]]
            b_indices <- b_data$indices
            b_center <- b_data$center
        }

        a_members <- length(a_indices)
        b_members <- length(b_indices)

        merged_center <-
            (a_members * a_center + b_members * b_center) / (a_members + b_members)

        if (a_center < b_center) {
            merged_indices <- c(a_indices, b_indices)
        } else {
            merged_indices <- c(b_indices, a_indices)
        }

        merge_data[[merge_index]] <- list(indices=merged_indices, center=merged_center)
    }

    clusters$order <- merge_data[[merges_count]]$indices

    return (clusters)
}

# Given a clustering which specified some data order, given we reorder the data ourselves, return a
# clustering that applies to the reordered data.
pre_ordered_hclust <- function(clusters) {
    old_of_new <- clusters$order
    new_of_old <- Matrix::invPerm(old_of_new)
    clusters$merge <- permute_merge(clusters$merge, new_of_old)
    clusters$order <- 1:length(clusters$order)
    return (clusters)
}

# Given an hclust merge array, return a new one that applies the same clustering if the data has
# been reordered.
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

#' Hierarchically cluster ordered data.
#'
#' Given a distance matrix for sorted objects, compute a hierarchical clustering preserving this
#' order. That is, this is similar to \code{hclust} with the constraint that the result's order is
#' always \code{1:N}.
#'
#' If an \code{order} is specified, assumes that the data will be re-ordered by this order. That is,
#' the indices in the returned \code{hclust} object will refer to the post-reorder data locations,
#' **not** to the current data locations.
#'
#" Currently, the only methods supported are \code{ward.D} and \code{ward.D2}.
#'
#' This can be applied to the results of \code{slanted_reorder}, to give a "plausible"
#' clustering for the data.
#'
#' @param distances A distances object (as created by \code{stats::dist}).
#' @param method The clustering method to use (only \code{ward.D} and \code{ward.D2} are supported).
#' @param order If specified, assume the data will be re-ordered by this order.
#' @param members Optionally, the number of members for each row/column of the distances (by default, one each).
#' @return A clustering object (as created by \code{hclust}).
#'
#' @export
#'
#' @examples
#' clusters <- slanter::oclust(dist(mtcars), order=1:dim(mtcars)[1])
#' clusters$order
oclust <- function(distances, method='ward.D2', order=NULL, members=NULL) {

    distances <- as.matrix(distances)
    stopifnot(dim(distances)[1] == dim(distances)[2])
    entities_count <- dim(distances)[1]

    if (method == 'ward.D2') {
        distances <- distances * distances
        sqrt_height <- TRUE
    } else {
        stopifnot(method %in% c('ward.D', 'ward.D2'))
        sqrt_height <- FALSE
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
