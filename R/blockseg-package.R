##' blockseg package
##'
##' This package is designed to segment a matrix in blocks with constant values.
##'
##' @section Features: Package for the segmentation of the rows and columns inducing a grid.
##'
##' @section Algorithm: \code{\linkS4class{blockSeg}}, \code{\linkS4class{stab.blockSeg}}
##'
##' @section Technical remarks: Display of the result with \code{\link{plot,blockSeg-method}} and \code{\link{plot,stab.blockSeg-method}} and the evolution
##' with \code{\link{predict,blockSeg-method}} and \code{\link{evolution,stab.blockSeg-method}}.
##'
##' @name blockseg-package
##' @docType package
##' @author Julien Chiquet \email{julien.chiquet@@gmail.com}
##' @author Vincent Brault \email{vincent.brault@@univ-grenoble-alpes.fr}
##'
##' @references BRAULT V, CHIQUET J. and LEVY-LEDUC C.Efficient block boundaries estimation in
##' block-wise constant matrices: An application to HiC data, Electronic Journal of Statistics,
##' Volume 11, Number 1 (2017), 1570-1599 <doi:10.1214/17-EJS1270>.
##'
##' @import methods
##' @import Matrix
##' @import ggplot2
##' @importFrom grDevices gray rgb
##' @importFrom graphics abline axis par title
##' @importFrom utils setTxtProgressBar txtProgressBar
##' @importFrom parallel mclapply detectCores
##' @importFrom reshape2 melt
##' @importFrom stats predict residuals deviance approx rnorm setNames
##' @useDynLib blockseg
##' @importFrom Rcpp sourceCpp
NULL
