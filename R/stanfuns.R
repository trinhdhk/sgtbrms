#' Expose stanvars of SGT and Skew T for brm call
#' @description
#' This function exposes brms stanvar object for brm call. Adding stanvar to the conventional call in brms is analogous to an _sgt call.
#' @param raw (boolean, default=FALSE) whether to expose raw Stan string, or brms stanvar object. The former can be used for Rstan object.
#' @return stanvar object, or if raw=TRUE, a string
#' @export
expose_sgt_stanvar <- function(raw=FALSE){
  stan_code <- paste(readLines(system.file('stan', 'sgt.stan', package = 'sgtbrms')), collapse = '\n')
  if (raw) return(stan_code)
  brms::stanvar(scode=stan_code, block='functions')
}
