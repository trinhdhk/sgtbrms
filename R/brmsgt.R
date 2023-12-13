.__brm_wrap__ <- function(brm_fn = '', .report_call = FALSE){
  fn <- get(brm_fn, envir=asNamespace('brms'))
  old_arg <- rlang::fn_fmls(fn)
  old_arg$family <- quote(sgt())
  old_arg$prior <- quote(sgt_default_prior())
  rlang::new_function(args=old_arg,
                      body = eval(enquote(substitute({
                        args <- rlang::fn_fmls_syms()
                        args_sub <- rlang::fn_fmls()
                        if (family$name == "skew_t" && (!length(match.call()$prior)))
                          args$prior <- quote(skew_t_default_prior())
                        sgt_vars <- brms::stanvar(scode=paste(readLines(system.file('stan', 'sgt.stan', package = 'sgtbrms')), collapse = '\n'), block = 'functions')
                        args$stanvars <- quote(c(sgt_vars, stanvars))
                        if (.report_call) cat('Call to: ', deparse1(rlang::call2(brm_fn, !!!args, .ns='brms')), '\n')
                        rlang::eval_tidy(rlang::call2(brm_fn, !!!args, .ns='brms'))
                      }))))
}



#' Wrapper functions fitting brm model for sgt and skew_t families
#' @description These functions callback to corresponding functions in brms package but with added stan codes for SGT and skew_t families
#' @inheritParams brms::brm
#' @seealso \link[brms]{brm}
#' @export
#' @return same as original function
brm_sgt <- .__brm_wrap__('brm')

#' @rdname brm_sgt
#' @inheritParams brms::make_stancode
#' @export
make_stancode_sgt <- .__brm_wrap__('make_stancode')

#' @rdname brm_sgt
#' @inheritParams brms::make_standata
#' @export
make_standata_sgt <- .__brm_wrap__('make_standata')

