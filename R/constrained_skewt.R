#' @rdname sgt
#' @param link_qmhalf (string, default: log) link function for qmhalf = q-1/2. Must be > 0.
#' @export
constrained_skew_t <- function(link='identity', link_sigma='log', link_qmhalf='log', link_lambdap1half='logit'){
  brms::custom_family(
    name='constrained_skew_t',
    dpars=c('mu', 'sigma', 'lambdap1half', 'qmhalf'),
    lb = c(NA, 0, 0, 0),
    ub = c(NA, NA, 1, NA),
    links=c(mu=link, sigma=link_sigma, lambdap1half=link_lambdap1half, qmhalf=link_qmhalf, pq=link_q),
    posterior_predict = posterior_predict_constrained_skew_t,
    posterior_epred = posterior_epred_constrained_skew_t,
    log_lik = log_lik_constrained_skew_t
  )
}

#' @rdname sgt_default_prior
#' @export
constrained_skew_t_default_prior <- function(params = 'all', exclude=FALSE){
  default_priors <-
    list(
      brms::prior(student_t(3, 0, 2.5), class='sigma', lb=0),
      brms::prior(beta(2,2), class='lambdap1half', lb=0, ub=1),
      brms::prior(gamma(3, .1), class='qmhalf', lb=0)
    )

  which <- c('sigma', 'lambdap1half', 'qmhalf')
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
  qmhalf <- brms::get_dpar(prep, "qmhalf", i=i)
  q <- qmhalf + 1/2
  p <- 2
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
  qmhalf <- brms::get_dpar(prep, "qmhalf", i=i)
  q <- qmhalf + 1/2
  p <- 2
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
  qmhalf <- brms::get_dpar(prep, "qmhalf")
  q <-  qmhalf + 1/2
  p <- 2
  lambda <- lambdap1half * 2 - 1
  mu - (2 * sigma * lambda * q^(1/p) * beta(2/p, q - 1/p))/beta(1/p, q)
}
