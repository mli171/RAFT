#' Fit a rank-based AFT model
#'
#' Estimates regression parameters for the accelerated failure time (AFT) model (censored
#' linear regression model) for right-censored data using R-estimation with expected (imputed) ranks.
#' Can estimate either unconstrained parameters (if \code{Gamma} is NULL) or constrained
#' parameters that satisfy \eqn{\Gamma \beta = b} for user-supplied values of \code{Gamma}
#' and \code{b}.  The function raft will check inputs, remove observations that have missing values of
#' covariates, centers \code{y} and \code{x} and calculates the matrix \code{Lambda} which
#' describes the part of parameter space that is orthgonal to the constraints.
#'
#' The returned object is designed to be used directly by inference functions
#' \code{raft.score.test()} and \code{raft.wald()}.
#'
#' @param beta Optional numeric vector of initial values for the coefficient vector.
#'   If \code{NULL}, initialized to a zero vector of length \code{ncol(x)}.
#' @param y Numeric vector of observed (censored) outcomes (log-failure-times for the AFT).
#' @param x Covariate matrix or vector. A data frame is coerced to a matrix. A vector
#'   converted to a single-column matrix.
#' @param delta Numeric/Integer vector of event indicators: \code{1} for uncensored
#'   observations and \code{0} for right-censored observations.
#' @param Gamma Optional matrix defining linear constraints/hypotheses of the form
#'   \eqn{\Gamma \beta = b}. If provided, \code{Gamma} must have full row rank and
#'   \code{b} must be supplied.
#' @param b Optional numeric vector defining the right-hand side in \eqn{\Gamma \beta = b}.
#'   Required if \code{Gamma} is provided.
#' @param tol Numeric tolerance used for uniroot solver.
#' @param half.width Numeric half-width of initial interval used by uniroot to find solutions of score equations. Default value is 0.5.
#' @param n.iter.max Integer maximum iterations for iterative root-finder for multivariate beta.
#' @param a Optional weight function \eqn{a(u,m,n)} (see package weight functions).
#' @param A Optional anti-derivative of weight function \eqn{a(u, m, n)} (see package weight functions).
#' @param m Optional numeric parameter used in a(u,m,n) and A(u,m,n).  Default value is NULL (appropriate if m not used in a(u,m,n) or A(u,m,n))
#' @param n Optional numeric parameter used in a(u,m,n) and A(u,m,n).  Default value is NULL (appropriate if m not used in a(u,m,n) or A(u,m,n))
#'
#' @return
#' A list containing the fitted parameter estimate(s), constraint matrices, and
#' the original data required for downstream inference (e.g., \code{raft.score.test()},
#' \code{raft.wald()}). The list has the following components:
#'
#' \describe{
#'   \item{beta}{Numeric vector of length \code{ncol(x)} giving the regression parameter
#'     estimates when there are no constraints (i.e., \code{is.null(Gamma)= TRUE }).
#'     Otherwise \code{NULL}.}
#'   \item{beta.r}{Numeric vector of length \code{ncol(x)} giving the regression parameter
#'     estimates under constraints (i.e., when \code{is.null(Gamma)=FALSE}). Otherwise \code{NULL}}
#'   \item{Gamma}{Matrix defining null hypothesis \eqn{\Gamma\beta=b}. \code{NULL} for simple Null \eqn{\beta=0}.}
#'   \item{Lambda}{Complementary matrix with rows spanning the null space of \eqn{\Gamma^T} (via \code{t(MASS:Null(t(Gamma)))}).}
#'   \item{n.gamma}{Integer. Number of rows of \code{Gamma} (0 if \code{Gamma} is \code{NULL}).}
#'   \item{n.lambda}{Integer. Number of rows of \code{Lambda} (0 if \code{Lambda} is \code{NULL}).}
#'   \item{Gamma.ginv}{Generalized inverse of \code{Gamma} (via \code{MASS::ginv()}),
#'     or \code{NULL} if \code{Gamma} is \code{NULL}.}
#'   \item{Lambda.ginv}{Generalized inverse of \code{Lambda} (via \code{MASS::ginv()}),
#'     or \code{NULL} if \code{Lambda} is \code{NULL}.}
#'   \item{b}{Right-hand side vector in \eqn{\Gamma\beta=b}. \code{NULL} if \code{Gamma} is \code{NULL}.}
#'   \item{y}{Numeric vector of observed outcomes used in fitting (after any deletions/centering performed by \code{raft()}).}
#'   \item{x}{Design matrix used in fitting (after any deletions/centering performed by \code{raft()}).}
#'   \item{delta}{Censoring indicator vector (1 = uncensored, 0 = right-censored).}
#'   \item{a}{Weight function \eqn{A(u,m,n)}.}
#'   \item{A}{Anti-derivative of weight function \eqn{a(u,m,n)}.}
#'   \item{m}{Parameter used in A(u,m,n) and a(u,m,n).  Set to \code{NULL} if not needed}
#'   \item{n}{Optional parameter used in A(u,m,n) and a(u,m,n).  Set to \code{NULL} if not needed.}
#' }
#'
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
#' }
#'
#' @seealso \code{\link{raft.score.test}}, \code{\link{raft.wald}}.
#' @export
raft = function(beta=NULL, y, x, delta, Gamma=NULL, b=NULL, tol=10^-12, half.width=0.5, n.iter.max=100,
                A=NULL, a=NULL, m=NULL, n=NULL){
  #
  #  checks input for problems, sets up Gamma, Lambda etc. and then calls estimate.e.rank.aft() to estimate beta or beta_r.
  #  output is used as input for hypothesis testing modules
  #

  if ('data.frame' %in% class(x) ) x=as.matrix(x)
  if (!is.matrix(x)) {x = matrix(x, ncol=1)}
  if (is.null(beta)) beta = rep(0,ncol(x))
  if (!is.null(Gamma)) {
    if( !('matrix' %in% class(Gamma)) ) Gamma=matrix(Gamma,nrow=1)
    Gamma.rank=qr(Gamma)$rank
    if( Gamma.rank!=nrow(Gamma)) {stop("Gamma matrix provided does not have full row rank")}
    if (is.null(b)) {stop("b must be specified if Gamma is specified")}
  }

  # remove missing values
  pNA.x = which(apply(x, 1, function(x) any(is.na(x))))
  pNA.y = which( is.na(y) )
  pNA = unique(c(pNA.x, pNA.y))
  if(length(pNA) > 0){
    warnings("\n\n Missing values deleted!")
    y = y[-pNA]
    x = x[-pNA,]
  }

  # center x and y
  x = scale(x, center=TRUE, scale=FALSE)
  y=y-mean(y)

  n.data = nrow(x)
  n.param = ncol(x)
  if(is.null(Gamma))  n.test = 0
  else n.test=nrow(Gamma)

  # set up Gamma and Lambda matrices
  if (is.null(Gamma)) Lambda=diag(n.param)
  else if (Gamma.rank==n.param) Lambda=NULL
  else Lambda=t( MASS::Null( t(Gamma) ) )

  if (is.null(Gamma)) Gamma.ginv=NULL
  else Gamma.ginv=MASS::ginv(Gamma)

  if (is.null(Lambda)) Lambda.ginv=NULL
  else Lambda.ginv=MASS::ginv(Lambda)

  n.gamma=ifelse(is.null(Gamma),0,nrow(Gamma))
  n.lambda=ifelse(is.null(Lambda),0,nrow(Lambda))

  if (is.null(a)) a1=1
  else a1=a(1,m,n)

  est.res=estimate.e.rank.aft(y=y, x=x, delta=delta, Gamma=Gamma, Lambda=Lambda, Gamma.ginv=Gamma.ginv,
                              Lambda.ginv=Lambda.ginv, b=b, beta=beta, n.gamma=n.gamma, n.lambda=n.lambda,
                              tol=tol, half.width=half.width, n.iter.max=n.iter.max,
                              A=A, a1=a1, m=m, n=n, offset=NULL)
  return(est.res)
}
