# RAFT

## R-Estimation with Right-Censored Data

Authors: Glen Satten (GSatten@emory.edu), Mo Li (mo.li@louisiana.edu), Ni Zhao (nzhao10@jhu.edu), Robert Strawderman (Robert_Strawderman@URMC.Rochester.edu)

## Overview
This R package provides tools for R-estimation in linear models with right-censored outcomes, extending classical rank-based R-estimators to seetings where the response is subject to censoring. The goal is to enable robust, rank-based estimation of regression effects without relying on fully parametric assumptions for the error distribution, while properly accounting for censoring.

### Key algorithmic components

- a
- b
- c
- d

### Intended use cases

- Regression modeling when outcomes are right-censored (e.g., time-to-event endpoints)
- Robust estimation under departures from normal errors
- Analyses where rank-based methods are preferred for stability or interpretability
- Method development and reproducible research aligned with estimating-equation theory

## Package download and installation
You can install the version of RAFT from Github:

```{r}
# install.packages("remotes")
remotes::install_github("mli171/RAFT", build_vignettes = TRUE, dependencies = TRUE)
```

## Open the Vignette in R for more details

```{r}
browseVignettes("RAFT")
```

## RAFT

The main function in RAFT package is:

```{r}
raft()
```

## An example of using the caft function

Next, we apply **RAFT** to the simulated dataset.

### Generate right censored data with specified quantile function.

Here, we use normal quantile function as an example. User can change to others, i.e, Weibull, Cauchy, etc. Quantile function examples can be found in QuantileError.R script.

```{r}
n.data  = 100
beta.T  = c(0.5, 1)
sd.y.T  = 1
mu.c.T  = 1.5
sd.c.T  = 2
alpha.T = 0
gamma.T = 0

# normal errors
Q.norm=function(u, sd.y) qnorm(u,0,sd.y)

set.seed(1234)

data = generate.aft(
  n.data = n.data,
  beta   = beta.T,
  sd.y   = sd.y.T,
  mu.c   = mu.c.T,
  sd.c   = sd.c.T,
  alpha  = alpha.T,
  gamma  = gamma.T,
  F.inv  = Q.norm
)
```

### Apply our RAFT approach for parameter estimation

```{r}
raft.res.NW = raft(
  x=data$x, 
  y=data$y, 
  delta=data$delta, 
  A=A.norm, 
  a1=a.norm(1, 0, n.data), 
  m=NULL, 
  n=n.data
)
raft.res.NW
```

### Score test

```{r}
scoretest.res = raft.score.test(est.rank.res = raft.res.NW)
```

### Wald test

```{r}
waldtest.res = raft.wald(est.raft.res = raft.res.NW)  
```

## References

Satten GA, Li M, Zhao N, Strawderman RL (2026). “R-Estimation with Right-Censored Data.” *arXiv preprint*, arXiv:2601.06685 [stat.ME, stat.AP], submitted January 10, 2026.

