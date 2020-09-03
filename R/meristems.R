#' Sample RNA data of similarity between batches of 1000 cells of tomato meristem cells.
#'
#' This is a simple matrix where each entry is the similarity (correlation) between a pair of
#' batches. Negative correlations were changed to zero to simplify the analysis.
#'
#' @docType data
#'
#' @usage data(meristems)
#'
#' @format A simple square matrix.
#'
#' @keywords datasets
#'
#' @examples
#' data(meristems)
#' similarity <- meristems
#' similarity[similarity < 0] = 0
#' slanter::sheatmap(meristems, order_data=similarity, show_rownames=FALSE, show_colnames=FALSE)
"meristems"
