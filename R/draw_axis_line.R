#' @title draw_axis_line
#'
#' @description The source code from the url:
#' \url{https://stackoverflow.com/questions/17753101/center-x-and-y-axis-with-ggplot2}
#'
#' @param length_x the length of x axis, a numeric
#' @param length_y the length of y axis, a numeric
#' @param tick_step the length of tick's step, a numeric
#' @param lab_step the length of label's step, a numeric
#'
#' @return an object of ggplot that we can add ggplot layer by "+"
#'
#' @importFrom ggplot2 aes ggplot geom_segment geom_text theme_void
#'
#' @export
#'
#' @examples library(ggplot2)
#' @examples p <- draw_axis_line(20, 4)
#' @examples p + geom_point(aes(x =1, y =1))
draw_axis_line <- function(length_x, length_y,
                           tick_step = NULL, lab_step = NULL){
  axis_x_begin <- -1*length_x
  axis_x_end <- length_x

  axis_y_begin  <- -1*length_y
  axis_y_end    <- length_y

  if (missing(tick_step))
    tick_step <- 1

  if (is.null(lab_step))
    lab_step <- 2

  # axis ticks data


  tick_x_frame <- data.frame(ticks = seq(axis_x_begin, axis_x_end,
                                         by = tick_step))

  tick_y_frame <-  data.frame(ticks = seq(axis_y_begin, axis_y_end,
                                          by = tick_step))

  # axis labels data
  lab_x_frame <- subset(data.frame(lab = seq(axis_x_begin, axis_x_end,
                                             by = lab_step), zero = 0),
                        lab != 0)

  lab_y_frame <- subset(data.frame(lab = seq(axis_y_begin, axis_y_end,
                                             by = lab_step),zero = 0),
                        lab != 0)

  # set tick length
  tick_x_length = 15/length(tick_x_frame$ticks)/2
  tick_y_length = 8/length(tick_y_frame$ticks)/2

  # set zero point

  data <- data.frame(x = 0, y = 0)
  p <- ggplot(data = data) +

    # draw axis line
    geom_segment(y = 0, yend = 0,
                 x = axis_x_begin,
                 xend = axis_x_end,
                 size = 0.5) +
    geom_segment(x = 0, xend = 0,
                 y = axis_y_begin,
                 yend = axis_y_end,
                 size = 0.5) +
    # x ticks
    geom_segment(data = tick_x_frame,
                 aes(x = ticks, xend = ticks,
                     y = 0, yend = 0 + tick_x_length)) +
    # y ticks
    geom_segment(data = tick_y_frame,
                 aes(x = 0, xend = 0 + tick_y_length,
                     y = ticks, yend = ticks)) +

    # labels
    geom_text(data=lab_x_frame, aes(x=lab, y=zero, label=lab), vjust = 1.5) +
    geom_text(data=lab_y_frame, aes(x=zero, y=lab, label=lab), hjust= 1.5) +
    theme_void()
  return(p)
}
