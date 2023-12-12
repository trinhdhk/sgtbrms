#' @rdname sgt_default_prior
#' @export
skew_t_default_prior <- function(params = 'all', exclude=FALSE){
  which <- c('sigma', 'lambdap1half', 'p', 'q')
  if ('all' %in% params) {
    if (exclude) stop("Are you sure?") else params <- which
  }
  include <- if (exclude) which(!which%in%params) else which(which%in%params)
  sgt_default_prior(params[params!='p'], exclude=FALSE)
}

#' @rdname sgt
#' @export
skew_t <- function(link_mu='identity', link_sigma='identity', link_q='log', link_lambdap1half='logit'){
  brms::custom_family(
    name='skew_t',
    dpars=c('mu', 'sigma', 'lambdap1half', 'q'),
    lb = c(NA, 0, 0, 0),
    ub = c(NA, NA, 1, NA),
    links=c(mu='identity', sigma='log', lambdap1half='logit', q='log'),
    posterior_predict = posterior_predict_skew_t,
    posterior_epred = posterior_epred_skew_t,
    log_lik = log_lik_skew_t
  )
}

#' @rdname brms-methods
log_lik_skew_t <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i=i)
  lambdap1half <- brms::get_dpar(prep, "lambdap1half", i=i)
  p <- 2
  q <- brms::get_dpar(prep, "q", i=i)
  y <- prep$data$Y[i]
  mapply(sgt::dsgt,
         x=y,
         mu=mu,
         sigma=sigma,
         lambda=lambdap1half*2-1,
         p=p,
         q=q,
         mean.cent=FALSE, var.adj=FALSE)
}

#' @rdname brms-methods
posterior_predict_skew_t <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i=i)
  lambdap1half <- brms::get_dpar(prep, "lambdap1half", i=i)
  p <- 2
  q <- brms::get_dpar(prep, "q", i=i)
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
posterior_epred_skew_t <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  sigma <- brms::get_dpar(prep, "sigma")
  lambdap1half <- brms::get_dpar(prep, "lambdap1half")
  q <- brms::get_dpar(prep, "q")
  p <- 2
  lambda <- lambdap1half * 2 - 1
  mu - (2 * sigma * lambda * q^(1/p) * beta(2/p, q - 1/p))/beta(1/p, q)
}
