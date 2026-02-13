# RAFT

## R-Estimation with Right-Censored Data

Authors: Glen Satten (GSatten@emory.edu), Mo Li (mo.li@louisiana.edu), Ni Zhao (nzhao10@jhu.edu), Robert Strawderman (Robert_Strawderman@URMC.Rochester.edu)

## Overview
This R package provides tools for R-estimation in linear models with right-censored outcomes, extending classical rank-based R-estimators to seetings where the response is subject to censoring. The goal is to enable robust, rank-based estimation of regression effects without relying on fully parametric assumptions for the error distribution, while properly accounting for censoring.

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


## Parameter estimation

```{r}
raft.res = raft(beta = NULL, y=time, x=x, delta=delta, Gamma = NULL, b = NULL, 
                tol = 10^-12, half.width = 0.5, n.iter.max = 100, A = NULL, 
                a = NULL, m = NULL, n = NULL)
```

## Score test

```{r}
raft.score.res = raft.score.test(est.rank.res = raft.res)
```

## Wald test

```{r}
raft.res.Weib = raft.wald(est.raft.res = raft.res)  
```

## References

Satten GA, Li M, Zhao N, Strawderman RL (2026). “R-Estimation with Right-Censored Data.” *arXiv preprint*, arXiv:2601.06685 [stat.ME, stat.AP], submitted January 10, 2026.

