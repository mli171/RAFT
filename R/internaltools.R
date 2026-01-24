#' @keywords internal
#' @noRd
estimate.e.rank.aft = function(y, x, delta, Gamma=NULL, Lambda=NULL, Gamma.ginv=NULL, Lambda.ginv=NULL, b=NULL, beta=NULL,
                               n.gamma, n.lambda, test=FALSE, tol, half.width, n.iter.max,
                               A=NULL, a1=NULL, m=NULL, n=NULL, offset=NULL) {
  #
  #  estimates parameters beta (or beta_r)
  #
  n.data=nrow=(x)
  n.var=ncol(x)

  if (n.lambda==0) {
    beta.r=as.vector( Gamma.ginv %*% b )
    res=list(beta=NULL, beta.r=beta.r, n.gamma=n.gamma, n.lambda=n.lambda, Gamma=Gamma, Lambda=Lambda, Gamma.ginv=NULL,
             Lambda.ginv=NULL, b=b, y=y, x=x, delta=delta, A=A, a1=a1, m=m, n=n, offset=offset)
    return(res)
  }

  n.param=n.lambda
  if( is.null(beta) ) lambda=rep(0,n.lambda)
  else lambda=as.vector( Lambda %*% beta )
  par=lambda
  #
  #   parameter estimation
  #
  if( n.param==1 ) {
    res=uniroot(f=fcn.e.score, interval=c(-5,5), extendInt='yes', tol=tol, y=y, delta=delta, x.cov=x, Gamma=Gamma,
                Lambda=Lambda, Gamma.ginv=Gamma.ginv, Lambda.ginv=Lambda.ginv, b=b, n.gamma=n.gamma, n.lambda=n.lambda,
                A=A, a1=a1, m=m, n=n, offset=offset )
    par=res$root
  }
  else {
    res=sus(fcn=fcn.e.score, beta=par, half.width=half.width, tol.diff=tol, n.iter.max=n.iter.max, y=y, delta=delta, x.cov=x, Gamma=Gamma,
            Lambda=Lambda, Gamma.ginv=Gamma.ginv, Lambda.ginv=Lambda.ginv, b=b, n.gamma=n.gamma, n.lambda=n.lambda,
            A=A, a1=a1, m=m, n=n, offset=offset )
    par=res$beta
    #        print( res$n.iter)
    #        print( 'using nleqslv' )
    #	    res=nleqslv(x=param0, fn=fcn.e.score, jac=jac, method='Broyden', y=y, delta=delta, x.cov=x, Gamma=Gamma, Lambda=Lambda, b=b, n.gamma=n.gamma, n.lambda=n.lambda,
    #	  	                control=list(xtol=10^-12, maxit=100))
    #		par=res$x
  }

  beta=NULL
  beta.r=NULL


  lambda=par[1:n.lambda]
  if (n.gamma==0) {
    beta=as.vector(Lambda.ginv %*% lambda)
  }
  else {
    beta.r=as.vector( Gamma.ginv %*% b + Lambda.ginv %*% lambda )
  }



  res=list(beta=beta, beta.r=beta.r, n.gamma=n.gamma, n.lambda=n.lambda, Gamma=Gamma, Lambda=Lambda,
           Gamma.ginv=Gamma.ginv, Lambda.ginv=Lambda.ginv,
           y=y, x=x, delta=delta, b=b, A=A, a1=a1, m=m, n=n)
  return(res)
}

#' @keywords internal
#' @noRd
fcn.e.score = function(x, y, x.cov, delta, Gamma, Lambda, Gamma.ginv, Lambda.ginv, b, n.gamma, n.lambda,
                       A, a1, m, n, offset) {
  #
  #   calculates the part of the score function that is being solved (as determined by Gamma, Lambda, b)
  #
  par=x
  x=x.cov

  n.p=length(par)
  W.n=rep(0,n.p)

  if (is.null(Gamma)) {
    beta=par
  }
  else if (is.null(Lambda)) {
    beta=as.vector( Gamma.ginv %*% b )
  }
  else {
    lambda=par
    beta=as.vector( Gamma.ginv %*% b + Lambda.ginv %*% lambda )
  }


  x.beta=colSums( beta*t(x) )

  resid=y - x.beta


  e.score.res=e.score(resid=resid, x=x, delta=delta, A=A, a1=a1, m=m, n=n)
  W.n=e.score.res$score
  S=Lambda %*% W.n
  if (!is.null(offset)) S=S-offset

  return(S)
}



#' @keywords internal
#' @noRd
e.score = function(resid, x, delta, A, a1, m, n) {


  resid.max=max(resid)
  j.max=which( resid==resid.max )
  delta[j.max]=1


  km.res=km(resid,delta)
  s.all=km.res$km.all[km.res$t.rank]
  s.star=NULL
  s.minus=NULL

  if (is.null(A)) {
    s.star=km.res$km.star[km.res$t.rank]
    e.rank=ifelse(delta==1, 1-s.star, (1-0.5*s.all) )
  }
  else {
    s.minus=km.res$km.minus.all[km.res$t.rank]
    delta.A=A(1-s.all,m,n)-A(1-s.minus,m,n)
    delta.F=s.minus-s.all
    exp.A.cens=ifelse(s.all>0,(A(1,m,n)-A(1-s.all,m,n))/s.all,a1)
    e.rank=ifelse(delta==1, delta.A/delta.F, exp.A.cens)
  }

  e.rank.uncentered=e.rank
  e.rank=e.rank - mean(e.rank)
  score = -colMeans(e.rank*x)


  res=list(score=score, e.rank=e.rank, s.all=s.all, s.star=s.star, s.minus=s.minus, e.rank.uncentered=e.rank.uncentered)
  return(res)
}



#' @keywords internal
#' @noRd
calculate.score.var = function(resid, delta, x, A, a1, m, n)	{
  #
  #    x.beta=colSums( beta*t(x) )
  #	resid=y-x.beta
  n.obs=length(resid)
  n.param=ncol(x)
  ord=order(resid)

  X.bar=colMeans(x)

  x=x[ord,,drop=FALSE]
  resid=resid[ord]
  delta=delta[ord]

  j.max=which( resid==resid[n.obs] )
  delta[j.max]=1
  #
  km.res=km(resid, delta)
  km=km.res$km
  km.all=km.res$km.all
  km.minus=km.res$km.minus
  fail.counts=km.res$fail.counts
  tot.counts=km.res$tot.counts
  n.unique.times=length(km.minus)

  h0=0
  h1=rep(0,n.param)
  h2=Sigma=Sigma.2=matrix(0,ncol=n.param, nrow=n.param)

  if (is.null(A)) w=km.minus/2        # factor of 1/2 replaces 0.25 for wilcoxon score
  else {
    delta.A=A(1-km,m,n)-A(1-km.minus,m,n)
    delta.F=km.minus-km
    exp.A.cens=ifelse(km>0,(A(1,m,n)-A(1-km,m,n))/km,a1)
    w=ifelse(fail.counts==0, 0, delta.A/delta.F - exp.A.cens)
  }

  i.stop=n.obs
  for (k in n.unique.times:1) {
    i.start=i.stop-tot.counts[k]+1
    for (i in i.start:i.stop) {
      h0=h0+1
      h1=h1+x[i,]
      h2=h2+outer(x[i,],x[i,])
    }
    i.stop=i.start-1
    Sigma=Sigma + w[k]^2*fail.counts[k]*( h2/h0 - outer(h1,h1)/h0^2 )
    h.temp=h1/h0-X.bar
    Sigma.2=Sigma.2 + w[k]^2*fail.counts[k]*outer( h.temp, h.temp )
  }

  Sigma=Sigma/n.obs
  Sigma.2=Sigma.2/n.obs

  wt=rep(1,n.obs)
  wt[delta==0]=1-km.all[delta==0]^2
  Sigma.1= t(x) %*% (wt*x) /(12*n.obs)


  res=list(Sigma=Sigma, Sigma.1=Sigma.1, Sigma.2=Sigma.2)

  return(res)

}




#' @keywords internal
#' @noRd
f.k = function(x, fcn, beta0, R.k, ...) {
  #
  #   utility function for sus, selects the kth equation and kth component to solve
  #
  args=list(...)
  beta = beta0 + x*R.k
  args=c(args,list(x=beta))
  res.f=do.call(fcn,args)
  value=sum( R.k*res.f )
  return(value)
}

#' @keywords internal
#' @noRd
sus = function(fcn, beta, half.width=0.5, tol.diff=10^-16, n.iter.max=100, ... ) {
  #
  #   sequential univariate solver (sus)
  #
  np=length(beta)
  diff=half.width
  tol=tol.diff
  n.iter=1
  del.beta=beta
  step.size=0.5
  R=diag(np)
  while ((n.iter<=n.iter.max)&(diff>tol.diff)) {
    diff=0
    for (k in 1:np) {
      interval=c(-half.width, half.width)
      res=uniroot( f=f.k, interval=interval, extendInt='yes', tol=tol, R.k=R[,k], fcn=fcn, beta0=beta, ...)
      del.beta[k]=res$root
      diff=diff+abs(res$root)
    }
    del.beta=del.beta/sqrt( sum(del.beta^2) )
    interval=c(-half.width, half.width)
    res=uniroot( f=f.k, interval=interval, extendInt='yes', tol=tol, R.k=del.beta, fcn=fcn, beta0=beta, ...)
    diff=sum( abs(res$root*del.beta) )
    beta=beta + res$root*del.beta
    n.iter=n.iter+1
    half.width=max( half.width/2, 10^-4 )
  }
  check=fcn(beta,...)
  res=list(diff=diff,beta=beta,check=check,n.iter=n.iter)
  return(res)
}


#' @keywords internal
#' @noRd
km = function(t, delta) {
  #
  #    calculates kaplan-meier using table function
  #
  n.obs=length(t)
  t.order=order(t)
  t.rank=rep(NA,n.obs)
  t.rank[t.order]=1:n.obs

  one=rep(1,n.obs)
  delta.2=c(delta,-one)
  tab=table( x=delta.2, y=c(t,t) )
  n.unique.times=dim(tab)[2]
  if (sum(delta)<n.obs) {
    tab=tab[3:1,]          #  first row=failures, 2nd row = censored obs, 3rd row = tot counts
    fail.counts=tab[1,]
    cens.counts=tab[2,]
    tot.counts=tab[3,]
  }
  if (sum(delta)==n.obs) {
    tab=tab[2:1,]
    fail.counts=tab[1,]
    cens.counts=rep(0,n.unique.times)
    tot.counts=tab[2,]
  }
  km.times=as.numeric( colnames(tab) )
  #  fail.counts=tab[1,]
  #  cens.counts=tab[2,]
  #  tot.counts=tab[3,]
  risk.set=n.obs - c(0, cumsum(tot.counts[-n.unique.times]) )
  names(risk.set)=names(tot.counts)
  #
  km=cumprod(1 - fail.counts/risk.set)
  km.minus=c(1,km[ -n.unique.times ] )
  names(km.minus)=names(km)

  km.all=rep(km, times=tot.counts)
  km.minus.all=rep(km.minus, times=tot.counts)
  km.star=0.5*(km.all + km.minus.all)
  km.all.times=rep(km.times, times=tot.counts)

  res=list(km.all=km.all, km.minus.all=km.minus.all, km=km, km.minus=km.minus, km.star=km.star,
           km.all.times=km.all.times, t.rank=t.rank, t.order=t.order, fail.counts=fail.counts, tot.counts=tot.counts)
  return(res)
}
