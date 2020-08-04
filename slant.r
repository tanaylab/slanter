#' Compute rows and columns orders which move high values close to the diagonal.
#'
#' For a matrix expressing the cross-similarity between two (possibly different)
#' sets of entities, this produces better results than clustering (e.g. as done
#' by \code{pheatmap}). This is because clustering does not care about the order
#' of each two sub-partitions. That is, clustering is as happy with \code{((2, 1),
#' (4, 3))} as it is with the more sensible \code{((1, 2), (3, 4))}. As a result,
#' visualizations of similarities using naive clustering can be misleading.
#'
#' Note the code actually works on the square of the data. This works OK-ish if
#' the similarity data is a correlation and one takes -1 (perfect negative
#' correlation) to be as strong a similarity indicator as +1 (perfect
#' correlation). If your data isn't like that, just make sure it is all
#' non-negative, where larger values are "more similar".
#'
#' @param data A rectangular matrix.
#' @param order_rows Whether to reorder the rows.
#' @param order_cols Whether to reorder the columns.
#' @return A list with two keys, \code{rows} and \code{cols}, which contain the order.
#'
#' @export
slanted_orders <- function(data, ..., order_rows=T, order_cols=T) {
    wrapr::stop_if_dot_args(substitute(list(...)), 'slanted_orders')
    rows_count <- dim(data)[1]
    cols_count <- dim(data)[2]

    rows_indices <- as.vector(1:rows_count)
    cols_indices <- as.vector(1:cols_count)

    rows_permutation <- rows_indices
    cols_permutation <- cols_indices

    if (order_rows || order_cols) {
        squared_data <- data * data

        reorder_phase <- function() {
            same_changes <- 0
            total_changes <- 1
            while (total_changes > 0 && same_changes < 10) {
                previous_changes <- total_changes
                total_changes <- 0
                if (order_rows) {
                    sum_indexed_rows <- rowSums(sweep(squared_data, 2, cols_indices, `*`))
                    sum_squared_rows <- rowSums(squared_data)
                    ideal_row_index <- sum_indexed_rows / sum_squared_rows

                    new_rows_permutation <- order(ideal_row_index)
                    row_changes <- sum(new_rows_permutation != rows_indices)
                    total_changes <- total_changes + row_changes

                    squared_data <<- squared_data[new_rows_permutation,]
                    rows_permutation <<- rows_permutation[new_rows_permutation]
                }

                if (order_cols) {
                    sum_indexed_cols <- colSums(sweep(squared_data, 1, rows_indices, `*`))
                    sum_squared_cols <- colSums(squared_data)
                    ideal_col_index <- sum_indexed_cols / sum_squared_cols

                    new_cols_permutation <- order(ideal_col_index)
                    col_changes <- sum(new_cols_permutation != cols_indices)
                    total_changes <- total_changes + row_changes

                    squared_data <<- squared_data[,new_cols_permutation]
                    cols_permutation <<- cols_permutation[new_cols_permutation]
                }

                if (total_changes == previous_changes) {
                    same_changes <- same_changes + 1
                } else {
                    same_changes <- 0
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
#' \code{pheatmap(slanted_reorder(data),cluster_rows=F,cluster_cols=F)}.
#'
#' @param data A rectangular matrix.
#' @return A matrix of the same shape whose rows and columns are a permutation of the input.
#'
#' @export
slanted_reorder <- function(data, ..., order_rows=T, order_cols=T) {
    wrapr::stop_if_dot_args(substitute(list(...)), 'slanted_reorder')
    orders <- slanted_orders(data, order_rows=order_rows, order_cols=order_cols)
    return (data[orders$rows, orders$cols])
}

#' Cluster ordered data.
#'
#' Given a distance matrix for sorted objects, compute a hierarchical clustering preserving this
#' order. That is, this is similar to `hclust` with the constraint that the result's order is
#' always `1:N`. This can be applied to the results of `slanted_reorder`, to give a "plausible"
#' clustering for the data.
#'
#' @param dist A distances object (as created by `dist`).
#' @param method The method of computing the clusters. Valid values are `agglomerative`
#'        (bottom-up, the default) or `divisive` (top-down).
#' @param aggregation How to aggregate distances between clusters; valid values are `mean`
#'        (the default), `min` and `max`.
#' @return A clustering object (as created by `hclust`).
#'
#' @export
oclust <- function(dist, ...,
                   method=c('agglomerative', 'divisive'), aggregation=c('mean', 'min', 'max')) {
    wrapr::stop_if_dot_args(substitute(list(...)), 'oclust')
    aggregation <- switch(match.arg(aggregation), mean='mean', min='min', max='max')
    method <- switch(match.arg(method), agglomerative='agglomerative', divisive='divisive')

    distances <- as.matrix(dist)
    stopifnot(dim(distances)[1] == dim(distances)[2])
    rows_count <- dim(distances)[1]

    if (method == 'agglomerative') {
        hclust <- bottom_up(distances, aggregation)
    } else {
        hclust <- top_down(distances, aggregation)
    }

    hclust$labels <- rownames(data)
    hclust$method <- sprintf('oclust.%s.%s', method, aggregation)
    if (!is.null(attr(dist, 'method'))) {
        hclust$dist.method <- attr(dist, 'method')
    }
    if (!is.null(attr(dist, 'Labels'))) {
        hclust$labels <- attr(dist, 'Labels')
    }
    hclust$order <- 1:rows_count

    class(hclust) <- 'hclust'
    return (hclust)
}

#' TODO: This can be made to run much faster.
bottom_up <- function(distances, aggregation) {
    aggregate <- switch(aggregation, mean=mean, min=min, max=max)

    rows_count <- dim(distances)[1]
    diag(distances) <- Inf

    merge <- matrix(0, nrow=rows_count - 1, ncol=2)
    height <- rep(0, rows_count - 1)
    merged_height <- rep(0, rows_count)
    groups <- -(1:rows_count)

    for (merge_index in 1:(rows_count - 1)) {
        adjacent_distances <- pracma::Diag(distances, 1)

        low_index <- which.min(adjacent_distances)
        high_index <- low_index + 1

        grouped_indices <- sort(groups[c(low_index, high_index)])

        merged_indices <- which(groups %in% grouped_indices)
        groups[merged_indices] <- merge_index
        merge[merge_index,] <- grouped_indices

        height[merge_index] <- max(merged_height[merged_indices]) + adjacent_distances[low_index]
        merged_height[merged_indices] <- height[merge_index]

        merged_distances <- apply(distances[,merged_indices], 1, aggregate)
        distances[,merged_indices] <- merged_distances
        distances[merged_indices,] <- rep(merged_distances, each=length(merged_indices))

        distances[merged_indices, merged_indices] <- Inf
    }

    return (list(merge=merge, height=height))
}

top_down <- function(distances, aggregation) {
    aggregate <- switch(aggregation, mean=cumsum, min=cummin, max=cummax)

    rows_count <- dim(distances)[1]

    merge <- matrix(0, nrow=rows_count - 1, ncol=2)
    height <- rep(0, rows_count - 1)
    accumulator <- list(merge=merge, height=height, merge_index=rows_count-1)

    accumulator <- top_down_divide(accumulator, distances, aggregate, aggregation, 1:rows_count)

    return (list(merge=accumulator$merge, height=accumulator$height))
}

top_down_divide <- function(accumulator, distances, aggregate, aggregation, indices_range) {
    rows_count <- dim(distances)[1]
    split_count <- length(indices_range)
    stopifnot(split_count > 1)

    effective_distances <- distances[indices_range, rev(indices_range)]
    effective_distances <- apply(effective_distances, 2, aggregate)
    effective_distances <- t(apply(t(effective_distances), 2, aggregate))  # TODO
    effective_distances <- effective_distances[1:split_count, split_count:1]

    candidate_distances <- pracma::Diag(effective_distances, 1)
    if (aggregation == 'mean') {
        candidate_distances <- candidate_distances / 1:length(candidate_distances)
        candidate_distances <- candidate_distances / length(candidate_distances):1
    }

    split_position <- which.max(candidate_distances)
    split_index <- split_position + min(indices_range) - 1
    low_range <- min(indices_range):split_index
    high_range <- (split_index + 1):max(indices_range)
    stopifnot(length(low_range) < split_count)
    stopifnot(length(high_range) < split_count)

    merge_index <- accumulator$merge_index

    if (length(low_range) == 1) {
        low_index <- -min(low_range)
        low_height <- 0
    } else {
        low_index <- accumulator$merge_index - 1
        accumulator$merge_index <- low_index
        accumulator <- top_down_divide(accumulator, distances, aggregate, aggregation, low_range)
        low_height <- accumulator$height[low_index]
    }

    if (length(high_range) == 1) {
        high_index <- -min(high_range)
        high_height <- 0
    } else {
        high_index <- accumulator$merge_index - 1
        accumulator$merge_index <- high_index
        accumulator <- top_down_divide(accumulator, distances, aggregate, aggregation, high_range)
        high_height <- accumulator$height[high_index]
    }

    accumulator$height[merge_index] <- candidate_distances[split_position] + max(low_height,
                                                                                 high_height)
    accumulator$merge[merge_index,] <- sort(c(low_index, high_index))

    return (accumulator)
}

#' Enhanced version of `dist`.
#'
#' This also allows using the correlations as the basis for the distances.
#' If the method is `pearson`, `kendall` or `spearman`, then the distances
#' will be `2 - cor(t(data), method)`. Otherwise, `dist` will be used.
enhanced_dist <- function(data, method) {
    if (method == 'pearson' || method == 'kendall' || method == 'spearman') {
        if (method == 'pearson') {
            distances <- 2 - tgs_cor(t(data))
        } else {
            distances <- 2 - cor(t(data), method=method)
        }
        distances_attributes <- attributes(distances)
        distances_attributes$method <- method
        attributes(distances) <- distances_attributes
        return (distances)
    }
    return (dist(data, method))
}

#' Plot a heatmap with values as close to the diagonal as possible.
#'
#' Given a matrix expressing the cross-similarity between two (possibly different) sets of
#' entities, this uses \code{slanted_reorder} to move the high values close to the diagonal, then
#' computes an order-preserving clustering for visualizing the matrix with a dendrogram tree, and
#' passes all this to `pheatmap`.
#'
#' @param data A rectangular matrix
#' @param order The default for whether to order rows and columns.
#' @param order_rows Whether to reorder the rows.
#' @param order_cols Whether to reorder the columns.
#' @param cluster The default for whether to cluster rows and columns.
#' @param cluster_rows Whether to cluster the rows (specify `order_rows=F` if giving an `hclust`).
#' @param cluster_cols Whether to cluster the columns (specify `order_cols=F` if giving an `hclust`).
#' @param distance_function The function for computing distance matrices (by default, `enhanced_dist`).
#' @param clustering_distance The default method for computing distances (by default, `pearson`).
#' @param clustering_distance_rows The method for computing distances between rows.
#' @param clustering_distance_cols The method for computing distances between columns.
#' @param clustering_method The default method to use for clustering the ordered data (by default, `agglomerative`).
#' @param clustering_method_rows The method to use for clustering the ordered rows.
#' @param clustering_method_cols The method to use for clustering the ordered columns.
#' @param clustering_aggregation The default aggregation method of cluster distances (by default, `mean`).
#' @param clustering_method_rows How to aggregate distances of clusters of rows.
#' @param clustering_method_cols How to aggregate distances of clusters of columns.
#' @param ... Additional flags to pass to `pheatmap`.
#' @return Whatever `pheatmap` returns.
sheatmap <- function(data, ...,
                     order=T,
                     order_rows=NULL,
                     order_cols=NULL,
                     cluster=F,
                     cluster_rows=NULL,
                     cluster_cols=NULL,
                     distance_function=enhanced_dist,
                     clustering_distance='pearson',
                     clustering_distance_rows=NULL,
                     clustering_distance_cols=NULL,
                     clustering_method='agglomerative',
                     clustering_method_rows=NULL,
                     clustering_method_cols=NULL,
                     clustering_aggregation='mean',
                     clustering_aggregation_rows=NULL,
                     clustering_aggregation_cols=NULL) {
    if (is.null(order_rows)) { order_rows = order }
    if (is.null(order_cols)) { order_cols = order }
    if (is.null(cluster_rows)) { cluster_rows = cluster }
    if (is.null(cluster_cols)) { cluster_cols = cluster }
    if (is.null(clustering_distance_rows)) { clustering_distance_rows=clustering_distance }
    if (is.null(clustering_method_rows)) { clustering_method_rows = clustering_method }
    if (is.null(clustering_aggregation_rows)) { clustering_aggregation_rows = clustering_aggregation }
    if (is.null(clustering_distance_cols)) { clustering_distance_cols=clustering_distance }
    if (is.null(clustering_method_cols)) { clustering_method_cols = clustering_method }
    if (is.null(clustering_aggregation_cols)) { clustering_aggregation_cols = clustering_aggregation }

    data <- slanted_reorder(data, order_rows=order_rows, order_cols=order_cols)

    if (cluster_rows == T) {
        rows_distances <- distance_function(data, method=clustering_distance_rows)
        cluster_rows <- oclust(rows_distances,
                               method=clustering_method_rows,
                               aggregation=clustering_aggregation_rows)
    }

    if (cluster_cols == T) {
        cols_distances <- distance_function(t(data), method=clustering_distance_cols)
        cluster_cols <- oclust(cols_distances,
                               method=clustering_method_cols,
                               aggregation=clustering_aggregation_cols)
    }

    return (pheatmap::pheatmap(data, cluster_rows=cluster_rows, cluster_cols=cluster_cols, ...))
}
