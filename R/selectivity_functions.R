# Selectivity function library
# Shared across all notebook components

#' Standard logistic selectivity
#' @param age numeric vector of ages
#' @param a50 age at 50% selectivity
#' @param delta slope width (50% to 95%)
logistic_sel <- function(age, a50, delta) {
  1 / (1 + exp(-log(19) * (age - a50) / delta))
}

#' Double logistic selectivity (3-parameter)
#' @param age numeric vector of ages
#' @param p1 ascending slope (>0)
#' @param p2 horizontal shift
#' @param p3 descending slope (>0)
dbl_logistic_sel <- function(age, p1, p2, p3) {
  gamma1 <- p1 + p2
  gamma2 <- 2 * p1 + p2 + p3
  asc  <- 1 / (1 + exp(-log(19) * (age - gamma1) / p1))
  desc <- 1 - 1 / (1 + exp(-log(19) * (age - gamma2) / p3))
  pmin(asc * desc * 0.95^(-2), 1.0)
}

#' Lognormal prior specification from median and CV
#' @param median desired median on natural scale
#' @param cv coefficient of variation on natural scale
lognorm_from_median_cv <- function(median, cv) {
  sigma <- sqrt(log(cv^2 + 1))
  mu    <- log(median)
  list(mu = mu, sigma = sigma)
}

#' Spline-based selectivity (normalized to max = 1)
#' @param age numeric vector of ages
#' @param beta spline coefficient vector
#' @param B B-spline basis matrix
spline_sel <- function(age, beta, B) {
  log_sel <- B %*% beta
  sel <- exp(log_sel - max(log_sel))
  as.numeric(sel)
}
