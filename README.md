# sgtbrms
Extensions for brms to fit SGT model

This is a proof-of-concept extension and wrapper code to fit Skew generalised T familiy regression model. It shall work will brms and rqeuire Rstan with Stan >= 2.30. 

An example model of how it works.

```r
x <- rnorm(1000, 0, 1)
y <- 3*x + 3*rt(1000,3,0)
model <- brm_sgt(y ~ x, family=sgt(), prior=c(sgt_default_prior()), data=data.frame(x=x, y=y), chains=1, iter=1000)
```

I provide two wrappers for `make_stancode` and `brm` with suffix `_sgt`. Will add `make_standata_sgt` very soon. At the moment you can make your own by calling `r make_standata_sgt <- sgtbrms:::.__brm_wrap__('make_standata')`.

All code is based on the implementation is R package `sgt` by Carter Davis.

Trinh Dong, 2023.

