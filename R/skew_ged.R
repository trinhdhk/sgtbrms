#' @rdname sgt_default_prior
#' @export
skew_ged_default_prior <- function(params = 'all', exclude=FALSE){
  which <- c('sigma', 'lambdap1half', 'p', 'q')
  if ('all' %in% params) {
    if (exclude) stop("Are you sure?") else params <- which
  }
  include <- if (exclude) (!which%in%params) else (which%in%params)
  sgt_default_prior(which[include][which!='q'], exclude=FALSE)
}

#' @rdname sgt
#' @export
skew_ged <- function(link='identity', link_sigma='log', link_p='log', link_lambdap1half='logit'){
  brms::custom_family(
    name='skew_ged',
    dpars=c('mu', 'sigma', 'lambdap1half', 'p'),
    lb = c(NA, 0, 0, 0),
    ub = c(NA, NA, 1, NA),
    links=c(mu=link, sigma=link_sigma, lambdap1half=link_lambdap1half, p=link_p),
    posterior_predict = posterior_predict_skew_ged,
    posterior_epred = posterior_epred_skew_ged,
    log_lik = log_lik_skew_ged
  )
}

#' @rdname brms-methods
log_lik_skew_ged <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i=i)
  lambdap1half <- brms::get_dpar(prep, "lambdap1half", i=i)
  q <- Inf
  p <- brms::get_dpar(prep, "p", i=i)
  y <- prep$data$Y[i]
  mapply(sgt::dsgt,
         x=y,
         mu=mu,
         sigma=sigma,
         lambda=lambdap1half*2-1,
         p=p,
         q=q,
         mean.cent=FALSE, var.adj=FALSE, log=TRUE)
}

#' @rdname brms-methods
posterior_predict_skew_ged <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i=i)
  lambdap1half <- brms::get_dpar(prep, "lambdap1half", i=i)
  q <- Inf
  p <- brms::get_dpar(prep, "p", i=i)
  mapply(sgt::rsgt,
         n=1,
         mu=mu,
         sigma=sigma,
         lambda=lambdap1half*2-1,
         p=p,
         q=q,
         mean.cent=FALSE, var.adj=FALSE)
}

#' @rdname brms-methods
posterior_epred_skew_ged <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  sigma <- brms::get_dpar(prep, "sigma")
  lambdap1half <- brms::get_dpar(prep, "lambdap1half")
  q <- 1e9
  p <- brms::get_dpar(prep, "p", i=i)
  lambda <- lambdap1half * 2 - 1
  mu - (2 * sigma * lambda * q^(1/p) * beta(2/p, q - 1/p))/beta(1/p, q)
}
