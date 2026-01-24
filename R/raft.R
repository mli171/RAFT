raft = function(beta=NULL, y, x, delta, Gamma=NULL, b=NULL, tol=10^-12, half.width=0.5, n.iter.max=100,
                A=NULL, a1=0, m=NULL, n=NULL){
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
                              A=A, a1=a1, m=m, n=n, offset=NULL)
  return(est.res)
}
