#' @rdname sgt
#' @param link_pq (string, default: logm1) link function for pq = p*q, must be > 1.
#' @export
constrained_sgt <- function(link='identity', link_sigma='log', link_q='log', link_pq='logm1', link_lambdap1half='logit'){
  brms::custom_family(
    name='constrained_sgt',
    dpars=c('mu', 'sigma', 'lambdap1half', 'q', 'pq'),
    lb = c(NA, 0, 0, 0, 1),
    ub = c(NA, NA, 1, NA, NA),
    links=c(mu=link, sigma=link_sigma, lambdap1half=link_lambdap1half, q=link_q, pq=link_pq),
    posterior_predict = posterior_predict_constrained_sgt,
    posterior_epred = posterior_epred_constrained_sgt,
    log_lik = log_lik_constrained_sgt
  )
}

#' @rdname sgt_default_prior
#' @export
constrained_sgt_default_prior <- function(params = 'all', exclude=FALSE){
  default_priors <-
    list(
      brms::prior(student_t(3, 0, 2.5), class='sigma', lb=0),
      brms::prior(beta(2,2), class='lambdap1half', lb=0, ub=1),
      brms::prior(gamma(3, .1), class='q', lb=0),
      brms::prior(gamma(3, .1), class='pq', lb=1)
    )

  which <- c('sigma', 'lambdap1half', 'q', 'pq')
  if ('all' %in% params) {
    if (exclude) stop("Are you sure?") else params <- which
  }
  include <- if (exclude) which(!which%in%params) else which(which%in%params)
  out <- default_priors[include]
  if (!length(out)) return(brms::empty_prior())
  do.call(c, args=out)
}

#' @rdname brms-methods
log_lik_constrained_sgt <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i=i)
  lambdap1half <- brms::get_dpar(prep, "lambdap1half", i=i)
  q <- brms::get_dpar(prep, "q", i=i)
  pq <- brms::get_dpar(prep, "pq", i=i)
  p <- pq/q
  y <- prep$data$Y[i]
  mapply(sgt::dsgt,
         x=y,
         mu=mu,
         sigma=sigma,
         lambda=lambdap1half*2-1,
         p=p,
         q=q,
         mean.cent=FALSE, var.adj=FALSE,log=TRUE)
}

#' @rdname brms-methods
posterior_predict_constrained_sgt <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i=i)
  lambdap1half <- brms::get_dpar(prep, "lambdap1half", i=i)
  q <- brms::get_dpar(prep, "q", i=i)
  pq <- brms::get_dpar(prep, "pq", i=i)
  p <- pq/q
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
posterior_epred_constrained_sgt <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  sigma <- brms::get_dpar(prep, "sigma")
  lambdap1half <- brms::get_dpar(prep, "lambdap1half")
  q <- brms::get_dpar(prep, "q")
  pq <- brms::get_dpar(prep, "pq")
  p <- pq/q
  lambda <- lambdap1half * 2 - 1
  mu - (2 * sigma * lambda * q^(1/p) * beta(2/p, q - 1/p))/beta(1/p, q)
}
