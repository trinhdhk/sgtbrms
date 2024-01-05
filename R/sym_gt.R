#' @rdname sgt
#' @export
sym_gt <- function(link='identity', link_sigma='log', link_p='log', link_q='log'){
  brms::custom_family(
    name='sym_gt',
    dpars=c('mu', 'sigma', 'p', 'q'),
    lb = c(NA, 0, 0, 0),
    ub = c(NA, NA, NA, NA),
    links=c(mu=link, sigma=link_sigma, p=link_p, q=link_q),
    posterior_predict = posterior_predict_sym_gt,
    posterior_epred = posterior_epred_sym_gt,
    log_lik = log_lik_sym_gt
  )
}

#' @rdname sgt_default_prior
#' @export
sym_gt_default_prior <- function(params = 'all', exclude=FALSE){
  default_priors <-
    list(
      brms::prior(student_t(3, 0, 2.5), class='sigma', lb=0),
      brms::prior(gamma(3, .1), class='p', lb=0),
      brms::prior(gamma(3, .1), class='q', lb=0)
    )

  which <- c('sigma', 'p', 'q')
  if ('all' %in% params) {
    if (exclude) stop("Are you sure?") else params <- which
  }
  include <- if (exclude) which(!which%in%params) else which(which%in%params)
  out <- default_priors[include]
  if (!length(out)) return(brms::empty_prior())
  do.call(c, args=out)
}

#' @rdname brms-methods
log_lik_sym_gt <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i=i)
  p <- brms::get_dpar(prep, "p", i=i)
  q <- brms::get_dpar(prep, "q", i=i)
  y <- prep$data$Y[i]
  mapply(sgt::dsgt,
         x=y,
         mu=mu,
         sigma=sigma,
         lambda=0,
         p=p,
         q=q,
         mean.cent=FALSE, var.adj=FALSE,log=TRUE)
}

#' @rdname brms-methods
posterior_predict_sym_gt <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i=i)
  p <- brms::get_dpar(prep, "p", i=i)
  q <- brms::get_dpar(prep, "q", i=i)
  mapply(sgt::rsgt,
         n=1,
         mu=mu,
         sigma=sigma,
         lambda=0,
         p=p,
         q=q,
         mean.cent=FALSE, var.adj=FALSE)
}

#' @rdname brms-methods
posterior_epred_sym_gt <- function(prep) {
  brms::get_dpar(prep, "mu")
}
