#' inv_normalise1
#' @description Often in large GWAS the quantitative phenotype is forced to follow a Normal distribution by quantile normalisation (also inverse-Normal transfromation).
#' @description This way the effect of phenotypic outliers is diminished while keeping the ranking of the trait values constant. This also harmonizes the trait distributions across multiple cohorts by forcing them to look similar.
#' @description There are two basic ways of INT. One is to apply INT directly to the phenotype, and the INT-transformed phenotype is then regressed on genotypes and covariates. This is so called direct method.
#' @description In the in-direct method, the phenotype is first regressed on covariates to obtain residuals. Then the residuals of the regression in ascending order.
#' @description Further details are available on the supplements of FTO genotype is associated with phenotypic variability of body mass index \url{https://www.nature.com/articles/nature11401}
#'
#' @param x the phenotypes
#' @return The normalized phenotype
#' @export
#'
inv_normalise1 <- function(x) {
  return( qnorm( (rank(x, na.last = "keep") - 0.5) / sum(!is.na(x))))
  }

#' inv_normalise2
#' @description inv_normalise2 is more aggressive that is suitable for the high skewed dataset.
#' @param x the phenotypes
#' @return The normalizedphenotype
#'
#' @export

inv_normalise2 <- function(x) {
  set.seed(13)
  numPhenos = length(which(!is.na(x)))
  quantilePheno = (rank(x, na.last="keep", ties.method="random")-0.5)/numPhenos
  phenoIRNT = qnorm(quantilePheno)
  return(phenoIRNT);
}
