% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Lasso2plot.R
\name{Lasso2plot}
\alias{Lasso2plot}
\title{Extract and Visualize Lasso Information from glmnet package.}
\usage{
Lasso2plot(mod, format = "cv")
}
\arguments{
\item{format}{A list from Lasso or cv Lasso fit.}

\item{data}{A model from  from glmnet package.}
}
\value{
Return a ggplot2 list
}
\description{
\code{Lasso2plot} Extract and Return a ggplot2 list of Lasso Information from glmnet package.
}
\examples{
Not run:
x = matrix(rnorm(100 * 20), 100, 20)
y = rnorm(100)
mod <- glmnet(x, y)
mod <- cv.glmnet(x, y)
p <- Lasso2plot(mod, format = "cv")
}
