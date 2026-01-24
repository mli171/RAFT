#' Score test for rank-based AFT estimating equations under right-censoring
#'
#' Computes a chi-square score test for hypotheses
#' defined through linear constraints on the coefficient vector from a rank-based
#' accelerated failure time (AFT) model with right-censored outcomes.
#'
#' This function is intended to be called on the result object returned by
#' \code{raft()} and uses the score vector and a sandwich-type covariance
#' estimate for the score to form a quadratic form test statistic.
#'
#' @param est.rank.res A list returned by \code{raft()}.
#' @param var.type Character scalar specifying the variance estimator used for the
#'   score. One of:
#'   \itemize{
#'     \item \code{"martingale"} (default): ???
#'     \item \code{"conservative"}: ???
#'     \item \code{"difference"}: ???
#'   }
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{test}: the chi-square test statistic.
#'   \item \code{df}: degrees of freedom.
#'   \item \code{p.value}: p-value from the chi-square reference distribution.
#'   \item \code{v}: variance-covariance matrix used in the test.
#'   \item \code{W.n}: score vector after applying constraints (if applicable).
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
#' scoretest.res = raft.score.test(est.rank.res = raft.res.NW)
#' }
#'
#' @seealso \code{\link{raft}}, \code{\link{raft.wald}}.
#'
#' @export
raft.score.test = function(est.rank.res, var.type='martingale') {
  #
  #   calculates score test given output est.rank.res from call to estimate.e.rank.aft
  #
  y=est.rank.res$y
  x=est.rank.res$x
  delta=est.rank.res$delta
  Gamma=est.rank.res$Gamma
  Lambda=est.rank.res$Lambda
  Gamma.ginv=est.rank.res$Gamma.ginv
  Lambda.ginv=est.rank.res$Lambda.ginv
  n.gamma=est.rank.res$n.gamma
  n.lambda=est.rank.res$n.lambda
  beta=est.rank.res$beta
  if (is.null(beta)) beta=est.rank.res$beta.r
  b=est.rank.res$b
  n.data=length(y)
  wt=est.rank.res$wt
  A=est.rank.res$A
  a1=est.rank.res$a1
  m=est.rank.res$m
  n=est.rank.res$n

  #
  #   calculate score and variance-covariance of score function
  #

  x.beta=colSums( beta*t(x) )
  resid=y-x.beta
  W.n.beta=e.score(resid=resid, delta=delta, x=x, A=A, a1=a1, m=m, n=n)$score

  var.res=calculate.score.var(resid=resid, delta=delta, x=x, A=A, a1=a1, m=m, n=n)
  if (var.type=='martingale') v=var.res$Sigma
  else if (var.type=='conservative') v=var.res$Sigma.1
  else if (var.type=='difference') v=var.res$Sigma.1 - var.res$Sigma.2

  c.matrix=rbind(Gamma,Lambda)
  v=tcrossprod(c.matrix %*% v, c.matrix)

  #   calculate test statistics
  #

  if (is.null(Gamma)) {
    sigma.inv=solve(v,W.n.beta)            #   W.n.beta=W.n.lambda
    test=n.data*W.n.beta %*% sigma.inv
    df=n.lambda
    W.n=W.n.beta
  }
  else {
    W.n.gamma=Gamma %*% W.n.beta
    sigma.inv=solve(v)[1:n.gamma,1:n.gamma]
    test=n.data*t(W.n.gamma) %*% sigma.inv %*% W.n.gamma
    df=n.gamma
    W.n=W.n.gamma
  }


  p.value=pchisq(test, df=df, lower.tail=FALSE)
  res=list(test=test, df=df, p.value=p.value, v=v, W.n=W.n)
  return(res)

}
