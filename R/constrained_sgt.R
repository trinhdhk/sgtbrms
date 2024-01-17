#' @rdname sgt
#' @param link_logpq (character, default: identity) link function for logpq = log(p*q)-1.
#' @param link_logq (character, default: identity) link function for log(q).
#' @export
constrained_sgt <- function(link='identity', link_sigma='log', link_logq='identity', link_logpq='identity', link_lambdap1half='logit'){
  brms::custom_family(
    name='constrained_sgt',
    dpars=c('mu', 'sigma', 'lambdap1half', 'logpq', 'logq'),
    lb = c(NA, 0, 0, NA, NA),
    ub = c(NA, NA, 1, NA, NA),
    links=c(mu=link, sigma=link_sigma, lambdap1half=link_lambdap1half, logpq=link_logpq, logq=link_logq),
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
      brms::prior(student_t(3, 0, 2.5), class='logpq'),
      brms::prior(student_t(3, 0, 2.5), class='logq')
    )

  which <- c('sigma', 'lambdap1half', 'logpq','logq')
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
  logq <- brms::get_dpar(prep, "logq", i=i)
  logpq <- brms::get_dpar(prep, "logpq", i=i)
  q <- exp(logq)
  pq <- exp(logpq) + 1
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
  logq <- brms::get_dpar(prep, "logq", i=i)
  logpq <- brms::get_dpar(prep, "logpq", i=i)
  pq <- exp(logpq) + 1
  q <- exp(logq)
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
  logq <- brms::get_dpar(prep, "logq")
  logpq <- brms::get_dpar(prep, "logpq", i=i)
  pq <- exp(logpq) + 1
  q <- exp(logq)
  p <- pq/q
  lambda <- lambdap1half * 2 - 1
  mu - (2 * sigma * lambda * q^(1/p) * beta(2/p, q - 1/p))/beta(1/p, q)
}
