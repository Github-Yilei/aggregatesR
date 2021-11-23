#' @title binFreqTable
#'
#' @description The bin frequency table is the table that represent the
#' frequency for a range of values in a particular variable.
#' If we want to create a bin frequency table then table function with cut
#' and breaks function can be used.
#' On the other hand, it may not flexible when creating a grouped bin frequency table.
#'
#' @param x a numeric vector, matrix or data frame.
#' @param step Step size among the range of value.
#' @return The bin frequency table
#' @export
#'
#' @examples path <- system.file("extdata", "chr.tsv", package = "aggregatesR" )
#' @examples mydf <- read.delim(path, header = TRUE, sep = "\t")
#' @examples freqtable <- tapply(mydf$value, mydf$group, function(x) binFreqTable(x, step = 1000000))
#'
#'
binFreqTable <- function(x, step) {
  bins = seq(0, ceiling(max(x)/step)*step, by = step)
  freq = hist(x, breaks=bins, include.lowest = TRUE, plot = FALSE)

  ranges = paste(head(freq$breaks, -1), freq$breaks[-1], sep = " - ")

  return(data.frame(range = ranges, frequency = freq$counts))

}

