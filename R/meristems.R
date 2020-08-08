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
#' \donttest{slanter::sheatmap(showrownames=F, show_colnames=F)}
"meristems"
