.onLoad <- function(libname, pkgname){
  stan_ver <- rstan::stan_version() |> strsplit(split = '\\.')
  comp <- as.integer(stan_ver[[1]][1] >= 2) && as.integer(stan_ver[[1]][2] >= 30)
  if (!comp) {
    cli::cli_alert_danger('Stan version < 2.30. Compilation will likely fail.')
    cli::cli_alert_info('Please update rstan or install cmdstanr. Ensure that `rstan::stan_version()` or `cmdstanr::cmdstan_version()`>= 2.30')
    cli::cli_alert_warning('Note that post-hoc methods in brms still require `rstan` compiled models.')
  }
}
