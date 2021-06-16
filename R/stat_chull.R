#' @title stat_chull
#'
#' @description The source code from the documents of
#' @description official extension mechanism provided in ggplot2
#' @description \url{https://cran.r-project.org/web/packages/ggplot2/vignettes/extending-ggplot2.html}
#'
#' @param mapping Set of aesthetic mappings.
#' @param data The data to be displayed in this layer.
#' @param stat The statistical transformation to use on the data for this layer.
#'
#' @param inherit.aes whether inherit aesthetic mapping from `ggplot`
#' @param ... additional parameters
#'
#' @return give the convex hull of a set of points in ggplot2.
#'
#' @importFrom dplyr group_by mutate
#'
#' @importFrom ggplot2 ggproto aes Stat
#' @importFrom ggplot2 Stat

#' @format NULL
#'
StatChull <- ggproto("StatChull", Stat,
                     required_aes = c("x", "y"),
                     compute_group = function(data, scales) {
                       data[chull(data$x, data$y), , drop = FALSE]
                     })
#'
#'
#' @export

stat_chull <- function(mapping = NULL, data = NULL, geom = "polygon",
                       position = "identity", na.rm = FALSE, show.legend = NA,
                       inherit.aes = TRUE, ...) {
  layer(
    stat = StatChull, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

