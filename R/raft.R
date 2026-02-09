#' Fit a rank-based AFT model
#'
#' Convenience wrapper that validates inputs, removes missing values, centers data,
#' constructs the constraint matrices \code{Gamma} and \code{Lambda}, and then calls
#' internal function \code{estimate.e.rank.aft()} to compute a rank-based
#' accelerated failure time (AFT) estimator under right-censoring.
#'
#' The returned object is designed to be used directly by inference utilities such
#' as \code{raft.score.test()} and \code{raft.wald()}.
#'
#' @param beta Optional numeric vector of initial values for the coefficient vector.
#'   If \code{NULL}, initialized to a zero vector of length \code{ncol(x)}.
#' @param y Numeric vector of observed outcomes (typically log-times in AFT settings).
#' @param x Covariate matrix or vector. A data frame is coerced to a matrix. A vector
#'   is treated as a single-column design matrix.
#' @param delta Numeric/Integer vector of event indicators: \code{1} for uncensored
#'   observations and \code{0} for right-censored observations.
#' @param Gamma Optional matrix defining linear constraints/hypotheses of the form
#'   \eqn{\Gamma \beta = b}. If provided, \code{Gamma} must have full row rank and
#'   \code{b} must be supplied.
#' @param b Optional numeric vector defining the right-hand side in \eqn{\Gamma \beta = b}.
#'   Required if \code{Gamma} is provided.
#' @param tol Numeric tolerance.
#' @param half.width Numeric half-width.
#' @param n.iter.max Integer maximum iterations.
#' @param A Optional derivative weight function \eqn{A(u)} (see package weight functions).
#' @param a Numeric constant used by some score families (passed through).
#' @param m Optional numeric parameter for certain score families (passed through).
#' @param n Optional numeric parameter for certain score families (passed through).
#'
#' @return
#' A list containing the fitted coefficient estimate(s), constraint matrices, and
#' the original data required for downstream inference (e.g., \code{raft.score.test()},
#' \code{raft.wald()}). The list has the following components:
#'
#' \describe{
#'   \item{beta}{Numeric vector of length \code{ncol(x)} giving the coefficient
#'     estimate when there are no \code{Gamma} constraints (i.e., \code{n.gamma == 0}).
#'     Otherwise \code{NULL}.}
#'   \item{beta.r}{Numeric vector of length \code{ncol(x)} giving the coefficient
#'     estimate under constraints (when \code{n.gamma > 0}), computed as
#'     \eqn{\Gamma^{+} b + \Lambda^{+}\lambda}. If \code{n.lambda == 0} (fully
#'     determined by constraints), this equals \eqn{\Gamma^{+} b}. Otherwise \code{NULL}
#'     when \code{n.gamma == 0}.}
#'   \item{n.gamma}{Integer. Number of rows of \code{Gamma} (0 if \code{Gamma} is \code{NULL}).}
#'   \item{n.lambda}{Integer. Number of rows of \code{Lambda} (0 if \code{Lambda} is \code{NULL}).}
#'   \item{Gamma}{Constraint matrix defining \eqn{\Gamma\beta=b}. \code{NULL} if no constraints are specified.}
#'   \item{Lambda}{Complementary matrix spanning the null space of \eqn{\Gamma^\top}.
#'     When \code{Gamma} is \code{NULL}, \code{Lambda} is the identity matrix of
#'     dimension \code{ncol(x)}. May be \code{NULL} if \code{Gamma} has full row rank
#'     equal to \code{ncol(x)}.}
#'   \item{Gamma.ginv}{Generalized inverse of \code{Gamma} (via \code{MASS::ginv()}),
#'     or \code{NULL} if \code{Gamma} is \code{NULL}.}
#'   \item{Lambda.ginv}{Generalized inverse of \code{Lambda} (via \code{MASS::ginv()}),
#'     or \code{NULL} if \code{Lambda} is \code{NULL}.}
#'   \item{b}{Right-hand side vector in \eqn{\Gamma\beta=b}. \code{NULL} if \code{Gamma} is \code{NULL}.}
#'   \item{y}{Numeric vector of observed outcomes used in fitting (after any deletions/centering performed by \code{raft()}).}
#'   \item{x}{Design matrix used in fitting (after any deletions/centering performed by \code{raft()}).}
#'   \item{delta}{Event indicator vector used in fitting (1 = uncensored, 0 = right-censored).}
#'   \item{A}{Derivative weight function \eqn{A(u)} passed through from \code{raft()}.}
#'   \item{a}{Numeric constant passed through from \code{raft()}.}
#'   \item{m}{Optional score/weight parameter passed through from \code{raft()}.}
#'   \item{n}{Optional score/weight parameter passed through from \code{raft()}.}
#' }
#'
#'
#' @examples
#' \dontrun{
#' set.seed(1234)
#'
#' library('KMsurv')
#' data(larynx, package='KMsurv')
#' time=larynx$time
#' stage=larynx$stage
#' age=larynx$age
#' delta=larynx$delta
#' x=model.matrix(~stage+age+0)[,2:5]
#'
#' raft.res=raft(
#'   beta = NULL, y=time, x=x, delta=delta,
#'   tol = 10^-12, half.width = 0.5, n.iter.max = 100)
#' }
#'
#' @seealso \code{\link{raft.score.test}}, \code{\link{raft.wald}}.
#' @export
raft = function(beta=NULL, y, x, delta, Gamma=NULL, b=NULL, tol=10^-12, half.width=0.5, n.iter.max=100,
                A=NULL, a=0, m=NULL, n=NULL){
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

  est.res=estimate.e.rank.aft(y=y, x=x, delta=delta, Gamma=Gamma, Lambda=Lambda, Gamma.ginv=Gamma.ginv,
                              Lambda.ginv=Lambda.ginv, b=b, beta=beta, n.gamma=n.gamma, n.lambda=n.lambda,
                              tol=tol, half.width=half.width, n.iter.max=n.iter.max,
                              A=A, a1=a, m=m, n=n, offset=NULL)
  return(est.res)
}
