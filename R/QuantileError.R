#' Quantile functions for generating AFT error distributions
#'
#' A small set of quantile (inverse CDF) functions used by \code{\link{generate.aft}}
#' to generate model errors \eqn{\varepsilon} from uniforms \eqn{U \sim \mathrm{Unif}(0,1)}.
#'
#' Each function has the common interface \code{f(u, sd.y)} where:
#' \itemize{
#'   \item \code{u} is a numeric vector with values in \eqn{(0,1)}
#'   \item \code{sd.y} is a positive scale parameter
#' }
#'
#' @param u Numeric vector of probabilities in \eqn{(0,1)}.
#' @param sd.y Positive numeric scale parameter.
#'
#' @return Numeric vector of quantiles of the same length as \code{u}.
#'
#' @details
#' \itemize{
#'   \item \code{Q.norm}: Normal errors, \eqn{\varepsilon = \Phi^{-1}(u; 0, sd.y)}.
#'   \item \code{Q.weib}: Extreme-value (Gumbel-type) errors on the log scale derived
#'     from Weibull survival times. Specifically, if \eqn{T} is Weibull, then
#'     \eqn{\log(T)} follows an extreme-value distribution; here \code{sd.y = 1/k}
#'     where \eqn{k} is the Weibull shape parameter. A centering constant
#'     \eqn{\psi(1)} (digamma at 1) is subtracted so the location is centered.
#'   \item \code{Q.Cauchy}: Cauchy errors, \eqn{\varepsilon = sd.y \tan(\pi(u-1/2))}.
#' }
#'
#' @examples
#' u <- seq(0.01, 0.99, length.out = 5)
#' Q.norm(u, sd.y = 1)
#' Q.weib(u, sd.y = 1)
#' Q.Cauchy(u, sd.y = 1)
#'
#' # Use a quantile function in the AFT generator
#' \dontrun{
#' dat <- generate.aft(n.data = 200, beta = c(0.5, -0.25), sd.y = 1, F.inv = Q.Cauchy)
#' }
#'
#' @name aft_quantile_functions
#' @aliases Q.norm Q.weib Q.Cauchy
NULL

#' @rdname aft_quantile_functions
#' @export
Q.norm <- function(u, sd.y) qnorm(u, 0, sd.y)   # normal errors

#' @rdname aft_quantile_functions
#' @export
Q.weib <- function(u, sd.y) {
  # t ~ Weibull, log(t) ~ extreme value; sd.y = 1/k where k is Weibull shape parameter
  sd.y * (log(-log(1 - u)) - digamma(1))
}

#' @rdname aft_quantile_functions
#' @export
Q.Cauchy <- function(u, sd.y) sd.y * tan(pi * (u - 0.5))
