suppressPackageStartupMessages({
  library(reshape2)
})


#' Clamp a vector between a minimum and maximum value.
#'
#' @param x
#' @param min
#' @param max
#' @return
clamp <- function(x, min, max)
  ifelse(x > min, ifelse(x < max, x, max), min)

#' N-dimensional N-linear interpolation.
#'
#' Interpolation is calculated as a weighted sum of the 2^n
#' points nearest each data point.
#'
#' @param data an MxN array. (i.e. long format) of hte places to interpolate at.
#' @param grid An N-dimensional array. (i.e. wide format) of the data to
#' interpolate from.
#' @param ref A list of N vectors that establish grid coordinates.
#' Expects that \code{length(ref[[i]]) == dim(grid)[[i]]. Also expects
#' ascending values.
#' @param rule The rule (see \code{\link{approx}}). 1 means values outside the box
#' are NA, 2 means extending a constant value.
#' @return M interpolated values.
#' @author Peter Meilstrup
interp.nd <- function (data, grid, ref = lapply(dim(grid), seq_len),
                       rule=1)
{
  n <- vapply(ref, length, 0)
  ndim <- length(n)
  all(dim(grid) == n) || stop("Reference must match dims of grid")
  ncol(data) == ndim || stop("Dimensions of interpolated data must match")
  nout <- nrow(data)
  output <- vector("numeric", nout)
  loc <- sapply(1:ndim, function(dim) {
    x <- ref[[dim]]
    y <- 1:length(x)
    clamp(approx(x, y, data[,dim], rule=rule)$y, 1, length(x))
  })
  lower <- floor(loc)
  lower <- sapply(seq_along(ref),
                  function(dim) clamp(lower[, dim],
                                      1, length(ref[[dim]]) - 1))
  remainder <- loc - lower
  #weights shall be calculated for each of 2^N dimensions.
  neighbors <- as.matrix(do.call(expand.grid, rep(list(c(0,1)), ndim)))
  apply(neighbors, 1, function(neighbor) {
    neighbor_rep <- rep(1, nout) %o% neighbor
    lookup <- lower + neighbor_rep
    weights <-  (neighbor_rep * (remainder) + (1-neighbor_rep) * (1-remainder))
    weights <- Reduce(`*`, unlist(recursive=FALSE, apply(weights, 2, list)))
    output <<- output + weights*grid[lookup]
  })
  output
}
