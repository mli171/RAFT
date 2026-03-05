## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE
)

## ----include=TRUE-------------------------------------------------------------
library(RAFT)

## ----eval=FALSE---------------------------------------------------------------
# devtools::install_local("RAFT_1.0.tar.gz", dependencies = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# remotes::install_github("mli171/RAFT", build_vignettes = TRUE, dependencies = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# browseVignettes("RAFT")

## ----eval=FALSE---------------------------------------------------------------
# vignette("vignette", package = "RAFT")

## ----eval=TRUE----------------------------------------------------------------
   library('KMsurv')
   data(larynx, package='KMsurv')
   time=larynx$time
   stage=factor(larynx$stage)
   age=larynx$age
   delta=larynx$delta
   x=model.matrix(~stage+age+0)[,2:5]

## ----eval=TRUE----------------------------------------------------------------
raft.res=raft(beta = NULL, y=time, x=x, delta=delta, Gamma = NULL, b = NULL, tol = 10^-12,
     half.width = 0.5, n.iter.max = 100, A = NULL, a = NULL, m = NULL, n = NULL)

## ----eval=TRUE----------------------------------------------------------------
raft.res$beta

## ----eval=TRUE----------------------------------------------------------------
wald.res=raft.wald(est.raft.res=raft.res, half.width = 0.5, n.iter.max = 1000, tol = 10^-12,
          var.type = "martingale", Gamma = NULL, b = NULL, alpha=0.05)
wald.res$test
wald.res$p.value

## ----eval=TRUE----------------------------------------------------------------
wald.res$marginal.ci

## ----eval=TRUE----------------------------------------------------------------
wald.res$omega

## ----eval=TRUE----------------------------------------------------------------
G=c(0,1,-1,0)
b=0
wald.res.2=raft.wald(est.raft.res=raft.res, half.width = 0.5, n.iter.max = 1000, tol = 10^-12,
          var.type = "martingale", Gamma = G, b = b, omega=wald.res$omega)
wald.res.2$test
wald.res.2$p.value
wald.res.2$df

## ----eval=TRUE----------------------------------------------------------------
G=matrix( c(1,0,-1,0,0,1,-1,0), byrow=TRUE, nrow=2)
b=c(0,0)
raft.res.r=raft(beta = NULL,y=time, x=x, delta=delta, Gamma = G, b = b, tol = 10^-12,
     half.width = 0.5, n.iter.max = 100, A = NULL, a = NULL, m = NULL, n = NULL)
raft.res.r$beta.r

## ----eval=TRUE----------------------------------------------------------------
raft.score.res=raft.score.test(raft.res.r)
raft.score.res$test
raft.score.res$p.value

## ----eval=TRUE----------------------------------------------------------------
A.norm=function(u,m,n) -(n+1)/n*dnorm( qnorm( (n*u+0.5)/(n+1), 0, 1), 0, 1)
a.norm=function(u,m,n) qnorm( (n*u+0.5)/(n+1),0,1)

## ----eval=TRUE----------------------------------------------------------------
n.data = length(time)
raft.res.Norm = raft(y=time, x=x, delta=delta, A=A.norm, a=a.norm, m=NULL, n=n.data)
raft.res.Norm$beta
raft.res$beta

## ----eval=TRUE----------------------------------------------------------------
wald.raft.Norm=raft.wald(est.raft.res=raft.res.Norm, half.width = 0.5, n.iter.max = 1000, tol = 10^-12, var.type = "martingale", Gamma = G, b = b)
wald.raft.Norm$test
wald.raft.Norm$p.value
wald.raft.Norm$omega
wald.raft.Norm$marginal.ci
wald.res$test
wald.res$p.value
wald.res$omega
wald.res$marginal.ci

## ----eval=TRUE----------------------------------------------------------------
A.weib=function(u,m,n) (n+1)/n*( 1 - n*u/(n+1) )*log(1 - n*u/(n+1) )
a.weib=function(u,m,n) -1 - log(1 - n*u/(n+1) )

## ----eval=TRUE----------------------------------------------------------------
raft.res.Weib = raft(y=time, x=x, delta=delta, A=A.weib, a=a.weib, m=NULL, n=n.data)
raft.res.Weib$beta

## ----eval=TRUE----------------------------------------------------------------
wald.raft.Weib=raft.wald(est.raft.res=raft.res.Weib, half.width = 0.5, n.iter.max = 1000, tol = 10^-12, var.type = "martingale", Gamma = NULL, b = NULL, omega=NULL, alpha=0.05)
wald.raft.Weib$test
wald.raft.Weib$p.value
wald.raft.Weib$omega
wald.raft.Weib$marginal.ci
wald.res$test
wald.res$p.value
wald.res$omega
wald.res$marginal.ci

## ----eval=TRUE----------------------------------------------------------------
A.F=function(u,m,n) {
    F.inv=ifelse(u==1,0,qf(u,2*m,2*n))
    A=ifelse(u==1,0, -m^m*n^n*(F.inv)^m/(beta(m,n)*(n+m*F.inv)^(m+n)))
	return(A)
    }
a.F=function(u,m,n) {
    F.inv=ifelse(u==1,0,qf(u,2*m,2*n))
    a=ifelse(u==1,n,m*n*(F.inv-1)/(n+m*F.inv))
	return(a)
    }

## ----eval=TRUE----------------------------------------------------------------
A.Cauchy=function(u,m,n) -(sin(pi*u))^2
a.Cauchy=function(u,m,n) -(sin(2*pi*u))

## ----eval=TRUE----------------------------------------------------------------
raft.res.CW=raft(y=time, x=x, delta=delta, A=A.Cauchy, a=a.Cauchy, m=NULL, n=NULL)
raft.res.CW$beta
raft.res$beta
wald.res.CW=raft.wald(raft.res.CW)
wald.res.CW$test
wald.res.CW$p.value
wald.res.CW$omega
wald.res.CW$marginal.ci
wald.res$test
wald.res$p.value
wald.res$omega
wald.res$marginal.ci


