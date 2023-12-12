# sgtbrms
Extensions for brms to fit SGT model

This is a proof-of-concept extension and wrapper code to fit Skew generalised T familiy regression model. It shall work will brms and require Rstan with Stan >= 2.30. 

## Installation

```r
remotes::install_github('trinhdhk/sgtbrms')
```

## Use

An example model of how it works.

```r
x <- rnorm(1000, 0, 1)
y <- 3*x + 3*rt(1000,3,0)
model <- brm_sgt(y ~ x, family=sgt(), prior=c(sgt_default_prior()), data=data.frame(x=x, y=y), chains=1, iter=1000)
# or equivalently
# model <- brm(y ~ x, family=sgt(), prior=c(sgt_default_prior()), stanvar=expose_sgt_stanvar(), data=data.frame(x=x, y=y), chains=1, iter=1000) 
```


I provide wrappers for `make_stancode`, `make_standata` and `brm` with suffix `_sgt`. These are the same calls but with `sgt_default_prior()` and `expose_sgt_stanvar()` added. You can make your own by calling:

```r 
X_sgt <- sgtbrms:::.__brm_wrap__('X')
```
where X is a function in `brms`.
The special case Skew T family where `p = 2` is also provided with `skew_t()`. 

All code is based on the implementation in the R package `sgt` by Carter Davis: [https://cran.r-project.org/package=sgt](https://cran.r-project.org/package=sgt).

Trinh Dong, 2023.

