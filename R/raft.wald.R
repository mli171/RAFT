#' Wald test and coefficient covariance for RAFT (rank-based AFT) estimators
#'
#' Computes an estimated variance–covariance matrix for the RAFT coefficient
#' estimator \eqn{\hat\beta} and forms a Wald chi-square test for either:
#' \itemize{
#'   \item the global null \eqn{H_0: \beta = 0} (default), or
#'   \item linear restrictions \eqn{H_0: \Gamma \beta = b}.
#' }
#'
#' The variance–covariance estimate \eqn{\widehat\Omega} is obtained using a
#' perturbation approach commonly attributed to Huang-type methods in the AFT
#' estimating-equation literature: the score covariance is square-rooted and used
#' to define small offsets to the estimating equation; the resulting changes in
#' \eqn{\hat\beta} are used to approximate sensitivity and hence \eqn{\Omega}.
#'
#' @param est.raft.res A list returned by \code{raft()}.
#' @param half.width Numeric. Half-width parameter, typically controlling step
#' size. Default is 0.5.
#' @param n.iter.max Integer. Maximum iterations limit. Default is 1000.
#' @param tol Numeric tolerance. Default is 10e-12.
#' @param var.type Character scalar specifying which score-variance estimator to use
#'   \itemize{
#'     \item \code{"martingale"} (default): ???
#'     \item \code{"conservative"}: ???
#'     \item \code{"difference"}: ???
#'   }
#' @param Gamma Optional matrix (or vector coerced to a 1-row matrix) specifying
#'   linear restrictions for the Wald test. If \code{NULL}, tests \eqn{H_0:\beta=0}.
#' @param b Optional numeric vector specifying the right-hand side in \eqn{H_0:\Gamma\beta=b}.
#'   Only used when \code{Gamma} is provided. If \code{NULL}, \eqn{b} is treated as 0.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{omega}: estimated variance–covariance matrix \eqn{\widehat\Omega} of \eqn{\hat\beta}.
#'   \item \code{beta.var}: perturbation-based sensitivity matrix used to form \code{omega}.
#'   \item \code{test}: Wald chi-square test statistic.
#'   \item \code{p.value}: p-value from \code{pchisq(test, df, lower.tail = FALSE)}.
#'   \item \code{df}: degrees of freedom (\code{ncol(x)} if \code{Gamma} is \code{NULL}, else \code{nrow(Gamma)}).
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1234)
#'
#' data = generate.aft(
#'   n.data = 100,
#'   beta   = c(0.5, 1),
#'   sd.y   = 1,
#'   mu.c   = 1.5,
#'   sd.c   = 2,
#'   alpha  = 0,
#'   gamma  = 0,
#'   F.inv  = Q.norm
#' )
#'
#' raft.res.NW = raft(
#'   x=data$x,
#'   y=data$y,
#'   delta=data$delta,
#'   A=A.norm,
#'   a1=a.norm(1, 0, n.data),
#'   m=NULL,
#'   n=n.data
#' )
#'
#' waldtest.res = raft.wald(est.rank.res = raft.res.NW)
#' }
#'
#' @seealso \code{\link{raft}}, \code{\link{raft.score.test}}.
#'
#' @export
raft.wald = function(est.raft.res, half.width=0.5, n.iter.max=1000, tol=10^-12, var.type='martingale', Gamma=NULL, b=NULL) {

  #
  #   implements the huang method for calculating Omega, the variance-covariance matrix of beta_hat
  #

  beta=est.raft.res$beta
  y=est.raft.res$y
  x=est.raft.res$x
  delta=est.raft.res$delta
  n.data=length(y)
  n.lambda.est=est.raft.res$n.lambda
  n.gamma.est=est.raft.res$n.gamma
  Gamma.est=est.raft.res$Gamma
  Lambda.est=est.raft.res$Lambda
  Gamma.ginv.est=est.raft.res$Gamma.ginv
  Lambda.ginv.est=est.raft.res$Lambda.ginv
  b.est=est.raft.res$b
  A=est.raft.res$A
  a1=est.raft.res$a1
  m=est.raft.res$m
  n=est.raft.res$n

  n.param=ncol(x)
  x.beta=colSums( beta*t(x) )
  resid=y-x.beta

  var.res=calculate.score.var(resid=resid, delta=delta, x=x, A=A, a1=a1, m=m, n=n)
  if (var.type=='martingale') v=var.res$Sigma
  else if (var.type=='conservative') v=var.res$Sigma.1
  else if (var.type=='difference') v=var.res$Sigma.1 - var.res$Sigma.2

  v=v/n.data
  C=expm::sqrtm(v)


  beta.var.pos=beta.var.neg=beta.var.min=beta.var=matrix(0,n.param,n.param)

  for (k in 1:n.param) {


    pos.res=try(estimate.e.rank.aft(y=y, x=x, delta=delta, Gamma=Gamma.est, Lambda=Lambda.est, Gamma.ginv=Gamma.ginv.est,
                                    Lambda.ginv=Lambda.ginv.est, b=b.est, beta=beta, n.gamma=n.gamma.est, n.lambda=n.lambda.est,
                                    tol=tol, half.width=half.width, n.iter.max=n.iter.max,
                                    A=A, a1=a1, m=m, n=n, offset=C[,k]), TRUE)
    neg.res=try(estimate.e.rank.aft(y=y, x=x, delta=delta, Gamma=Gamma.est, Lambda=Lambda.est, Gamma.ginv=Gamma.ginv.est,
                                    Lambda.ginv=Lambda.ginv.est, b=b.est, beta=beta, n.gamma=n.gamma.est, n.lambda=n.lambda.est,
                                    tol=tol, half.width=half.width, n.iter.max=n.iter.max,
                                    A=A, a1=a1, m=m, n=n, offset=-C[,k]), TRUE)
    use.pos=!( class(pos.res)=='try-error' )
    use.neg=!( class(neg.res)=='try-error' )
    if (use.neg&use.pos) beta.var[,k]=0.5*(pos.res$beta-neg.res$beta)
    else if (use.neg&!use.pos) beta.var[,k]=beta-neg.res$beta
    else if (!use.neg&use.pos) beta.var[,k]=pos.res$beta-beta
    else stop('both negative and positive solutions failed')
  }
  omega=beta.var %*% t(beta.var)

  #
  #   calculate wald test
  #

  if (is.null(Gamma)) {
    test=beta %*% solve(omega, beta)
    df=n.param
    omega.z=omega
  }
  else {
    if( !('matrix' %in% class(Gamma)) ) Gamma=matrix(Gamma,nrow=1)
    z.gamma=Gamma %*% beta
    if (!is.null(b)) z.gamma=z.gamma-b
    omega.z=tcrossprod(Gamma %*% omega, Gamma)
    test=z.gamma %*% solve(omega.z, z.gamma)
    df=nrow(Gamma)
  }
  p.value=pchisq(test, df=df, lower.tail=FALSE)
  res=list(omega=omega, beta.var=beta.var, test=test, p.value=p.value, df=df)
  return(res)
}
