.onLoad <- function(libname, pkgname){
  stan_ver <- rstan::stan_version() |> strsplit(split = '\\.')
  comp <- as.integer(stan_ver[[1]][1] >= 2) && as.integer(stan_ver[[1]][2] >= 30)
  if (!comp) stop('Stan version < 2.30. Compilation would fail. Please update rstan or install cmdstan. Ensure that Stan version > 2.30')
}
