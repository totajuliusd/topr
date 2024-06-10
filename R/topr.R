#'  topr
#'
#' @description A package for viewing and annotating genetic association data
#' @name topr
#' @section topr functions:
#' The main plotting functions are:
#' * \code{\link{manhattan}} to create Manhattan plot of association results
#' * \code{\link{regionplot}} to create regional plots of association results for smaller genetic regions
#' @examples
#' library(topr)
#' # Create a manhattan plot using
#' manhattan(CD_UKBB)
#'
#' # Create a regional plot
#' regionplot(CD_UKBB, gene="IL23R")
#' @docType package
"_PACKAGE"
