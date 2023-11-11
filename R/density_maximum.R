#' inv_normalise1
#' @description A density plot takes a numeric variable to represent a smooth distribution curve over time. The peak of the density plot shows the maximum concentration of numeric data.
#'
#' @param x the variable
#' @return The peak of distribution plot
#' @export
#'
density_maximum <- function(data) {
  
  max_idx <- which.max(density(data)$y)
  max_x <- density(data)$x[max_idx]
  
}
