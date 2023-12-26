#' Skew generalised T family for brms model
#' @description Create a custom family object for brms model
#' @param link (character, default: identity) link function for the mode of SGT distribution
#' @param link_sigma (character, default: log)  link function for the scale of SGT distribution
#' @param link_p,link_q (character, default: log)  link function for the kurtosis parameters of SGT distribution
#' @param link_lambdap1half (character, default: logit) link function for the skew parameter of SGT distribution.
#' Note that the actual skew parameter of SGT is lambda where -1 < lambda < 1.
#' A proper link function shall be \emph{tanh}.
#' However, due to brms's limit options for link function, we lambda to lambdap1half = \eqn{\frac{lambda+1}{2}} and use the logit link.
#' Hence the name "lambda p(lus) 1 then half".
#' @details
#' The constrained_ versions constraint p*q > 1 to ensure identifiability of the mean,
#' but not as "generalised" due to the exclusion of the Cauchy-related branch.
#' The mu of the param is still the mode though.
#'
#' @return a `brms::custom_family` object
#' @importFrom brms custom_family
#' @seealso \link[brms]{dsgt}
#' @export
#' @examples
#' x <- rnorm(1000, 0, 1)
#' y <- 3*x + 3*rt(1000,3,0)
#' model <- brm_sgt(y ~ x, family=sgt(), prior=c(sgt_default_prior()), data=data.frame(x=x, y=y), chains=1, iter=1000)
sgt <- function(link='identity', link_sigma='log', link_p='log', link_q='log', link_lambdap1half='logit'){
  brms::custom_family(
    name='sgt',
    dpars=c('mu', 'sigma', 'lambdap1half', 'p', 'q'),
    lb = c(NA, 0, 0, 0, 0),
    ub = c(NA, NA, 1, NA, NA),
    links=c(mu=link, sigma=link_sigma, lambdap1half=link_lambdap1half, p=link_p, q=link_q),
    posterior_predict = posterior_predict_sgt,
    posterior_epred = posterior_epred_sgt,
    log_lik = log_lik_sgt
  )
}

#' Default priors for sgt family
#' @description
#' This function returns default set of prior for sgt fit. By default, returns prior for all parameters.
#' If you set some parameter sub-model, please exclude that from the param list, otherwise brms will complain.
#' @param params string vector (default="all") of params in sgt family, namely `sigma`, `lambdap1half`, `p`, `q`, a special value "all" will return priors for all.
#' @param exclude (boolean, default=FALSE) whether to include or exclude parameters in `params`.
#' @details
#' Importance notice that if you model lambdap1half as a sub-model, be aware that on logit scale, the non-informative prior
#' for intercept is logistic(0, 1) and NOT flat.
#'
#' @return object of class `brmsprior`
#' @export
#' @seealso \link[brms]{prior}
sgt_default_prior <- function(params = 'all', exclude=FALSE){
  default_priors <-
    list(
      brms::prior(student_t(3, 0, 2.5), class='sigma', lb=0),
      brms::prior(beta(2,2), class='lambdap1half', lb=0, ub=1),
      brms::prior(gamma(3, .1), class='p', lb=0),
      brms::prior(gamma(3, .1), class='q', lb=0)
    )

  which <- c('sigma', 'lambdap1half', 'p', 'q')
  if ('all' %in% params) {
    if (exclude) stop("Are you sure?") else params <- which
  }
  include <- if (exclude) which(!which%in%params) else which(which%in%params)
  out <- default_priors[include]
  if (!length(out)) return(brms::empty_prior())
  do.call(c, args=out)
}

#' Downstream methods for sgt family
#' @rdname brms-methods
#' @description Some prediction methods for brms sgt family. For internal use of brms methods only
#' @seealso \link[brms]{log_lik} \link[brms]{posterior_predict} \link[brms]{posterior_epred}
#' @param i iteration
#' @param prep a brms prepared object
log_lik_sgt <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i=i)
  lambdap1half <- brms::get_dpar(prep, "lambdap1half", i=i)
  p <- brms::get_dpar(prep, "p", i=i)
  q <- brms::get_dpar(prep, "q", i=i)
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
posterior_predict_sgt <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i=i)
  lambdap1half <- brms::get_dpar(prep, "lambdap1half", i=i)
  p <- brms::get_dpar(prep, "p", i=i)
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
posterior_epred_sgt <- function(prep) {
  mu <- brms::get_dpar(prep, "mu")
  sigma <- brms::get_dpar(prep, "sigma")
  lambdap1half <- brms::get_dpar(prep, "lambdap1half")
  q <- brms::get_dpar(prep, "q")
  p <- brms::get_dpar(prep, "p")
  lambda <- lambdap1half * 2 - 1
  mu - (2 * sigma * lambda * q^(1/p) * beta(2/p, q - 1/p))/beta(1/p, q)
}
