#' Weight functions for R-estimators (score functions and derivatives)
#'
#' A collection of **rank-based score/weight functions** \eqn{a(u)} and their
#' corresponding derivatives \eqn{A(u)} used in R-estimation and related estimating
#' equations. These functions are typically evaluated at scaled ranks
#' \eqn{u \in [0,1]} (including generalized ranks under censoring).
#'
#' The functions are grouped into four families:
#' \itemize{
#'   \item **Normal-score (van der Waerden)**: \code{\link{a.norm}}, \code{\link{A.norm}}
#'   \item **Weibull/log-type**: \code{\link{a.weib}}, \code{\link{A.weib}}
#'   \item **F-quantile-based**: \code{\link{a.F}}, \code{\link{A.F}}
#'   \item **Cauchy-type (sine-based)**: \code{\link{a.Cauchy}}, \code{\link{A.Cauchy}}
#' }
#'
#' @param x Numeric vector with values in \eqn{[0,1]}. Interpreted as scaled ranks.
#' @param m Numeric. Meaning depends on the weight family:
#' \itemize{
#'   \item For \code{a.F}/\code{A.F}: positive shape parameter used via \code{df1 = 2*m}.
#'   \item For \code{a.Cauchy}/\code{A.Cauchy}: positive Cauchy scale parameter (often \eqn{\gamma}).
#'   \item For \code{a.norm}/\code{A.norm} and \code{a.weib}/\code{A.weib}: unused (kept for a common signature).
#' }
#' @param n Numeric. Meaning depends on the weight family:
#' \itemize{
#'   \item For \code{a.norm}/\code{A.norm} and \code{a.weib}/\code{A.weib}: positive finite-sample scaling parameter
#'         (often a sample size or effective sample size used in rank scaling).
#'   \item For \code{a.F}/\code{A.F}: positive shape parameter used via \code{df2 = 2*n}
#'         (also appears in the closed-form expressions).
#'   \item For \code{a.Cauchy}/\code{A.Cauchy}: unused (kept for a common signature).
#' }
#'
#' @details
#' **Normal-score family**
#'
#' Uses the finite-sample mapping
#' \deqn{z(u) = \Phi^{-1}\left(\frac{n u + 0.5}{n+1}\right)}
#' and returns:
#' \itemize{
#'   \item \code{a.norm(u)} = \eqn{z(u)}
#'   \item \code{A.norm(u)} = \eqn{-\frac{n+1}{n}\,\phi(z(u))}
#' }
#' where \eqn{\Phi^{-1}} is the standard normal quantile and \eqn{\phi} its density.
#'
#' **Weibull/log-type family**
#'
#' Uses \eqn{t(u) = 1 - n u/(n+1)} and returns:
#' \itemize{
#'   \item \code{a.weib(u)} = \eqn{-1 - \log(t(u))}
#'   \item \code{A.weib(u)} = \eqn{\frac{n+1}{n}\,t(u)\,\log(t(u))}
#' }
#'
#' **F-quantile family**
#'
#' Uses the F-quantile \eqn{F^{-1}(u)} computed by \code{qf(u, df1 = 2*m, df2 = 2*n)}.
#' Since \eqn{F^{-1}(1)=\infty}, this implementation **guards against \code{x == 1}**
#' by substituting \code{F.inv = 0} when \code{x == 1}, returning a finite value.
#' In typical rank workflows, \code{x} should avoid being exactly 1.
#'
#' **Cauchy family**
#'
#' Sine-based scores parameterized by a positive scale \code{m}:
#' \itemize{
#'   \item \code{a.Cauchy(u)} = \eqn{-\sin(2\pi u)/m}
#'   \item \code{A.Cauchy(u)} = \eqn{-\sin^2(\pi u)/m}
#' }
#'
#' @return A numeric vector of the same length as \code{x}.
#'
#' @examples
#' u <- seq(0.01, 0.99, length.out = 5)
#'
#' # Normal scores
#' a.norm(u, m = NA, n = 100)
#' A.norm(u, m = NA, n = 100)
#'
#' # Weibull/log-type
#' a.weib(u, m = NA, n = 100)
#' A.weib(u, m = NA, n = 100)
#'
#' # F-based
#' a.F(u, m = 2, n = 5)
#' A.F(u, m = 2, n = 5)
#'
#' # Cauchy-type
#' a.Cauchy(u, m = 1, n = NA)
#' A.Cauchy(u, m = 1, n = NA)
#'
#' @name weight_functions
#' @aliases a.norm A.norm a.weib A.weib a.F A.F a.Cauchy A.Cauchy
NULL

# Normal weights

#' @rdname weight_functions
#' @export
A.norm <- function(x, m, n) {
  -(n + 1) / n * dnorm(qnorm((n * x + 0.5) / (n + 1), 0, 1), 0, 1)
}

#' @rdname weight_functions
#' @export
a.norm <- function(x, m, n) {
  qnorm((n * x + 0.5) / (n + 1), 0, 1)
}


# Weibull weights

#' @rdname weight_functions
#' @export
A.weib <- function(x, m, n) {
  (n + 1) / n * (1 - n * x / (n + 1)) * log(1 - n * x / (n + 1))
}

#' @rdname weight_functions
#' @export
a.weib <- function(x, m, n) {
  -1 - log(1 - n * x / (n + 1))
}

# F-quantile-based weights

#' @rdname weight_functions
#' @export
A.F <- function(x, m, n) {
  F.inv <- ifelse(x == 1, 0, qf(x, 2 * m, 2 * n))
  A <- ifelse(
    x == 1,
    0,
    -m^m * n^n * (F.inv)^m / (beta(m, n) * (n + m * F.inv)^(m + n))
  )
  return(A)
}

#' @rdname weight_functions
#' @export
a.F <- function(x, m, n) {
  F.inv <- ifelse(x == 1, 0, qf(x, 2 * m, 2 * n))
  a <- ifelse(x == 1, n, m * n * (F.inv - 1) / (n + m * F.inv))
  return(a)
}

# Cauchy-type weights (sine-based)

#' @rdname weight_functions
#' @export
A.Cauchy <- function(x, m, n) {
  -(sin(pi * x))^2 / m
}

#' @rdname weight_functions
#' @export
a.Cauchy <- function(x, m, n) {
  -(sin(2 * pi * x)) / m
}
