int sgn(real x){
  return (x >= 0 ? 1 : -1);
}

real sgt_lpdf(real x, real mu, real sigma, real lambdap1half, real p, real q){
  real lambda = lambdap1half*2 - 1;
  // real p = 2;
  return (log(p)-log(2)-log(sigma)-log(q)/p-lbeta(1/p,q)-(1/p+q)*log1p(pow(abs(x-mu),p)/(q*pow(sigma,p)*pow(1+lambda*sgn(x-mu),p))));
}

real sgtv_lpdf(vector x, real mu, real sigma, real lambdap1half, real p, real q){
  real lambda = lambdap1half*2 - 1;
  real lpdf = 0;
  for (i in 1:num_elements(x)){
    lpdf += sgt_lpdf(x[i], mu, sigma, lambdap1half, p, q);
  } 
  return lpdf;
}

real sgt_invcdf(real prob, real mu, real sigma, real lambdap1half, real p, real q){
  real lambda = lambdap1half*2 - 1;
  // real p = 2;
  int flip = prob > (1-lambda)/2;
  real lam = lambda;
  real pr = prob;
  real out;

  if (flip){
    pr = 1 - pr;
    lam = - lam;
  }

  out = sigma*(lam-1)*(1/(q*inv_inc_beta(1/p,q, 1-2*pr/(1-lam)))-1/q)^(-1/p);
  out = flip ? -out : out;
  return (out + mu);
}

real sgt_rng(real mu, real sigma, real lambdap1half, real p, real q){
  // real p = 2;
  real prob = uniform_rng(0,1);
  return (sgt_invcdf(prob, mu, sigma, lambdap1half, p, q));
}

real sgt_lcdf(real x, real mu, real sigma, real lambdap1half, real p, real q){
  real lambda = lambdap1half*2 - 1;
  real quant = x - mu;
  int flip = quant > 0;
  real lam = lambda;
  real out;
  if (flip){
    lam = -lam;
    quant = -quant;
  }
  out = log_sum_exp(
    (1-lam)/2,
    log(lam-1)-log(2)+beta_lcdf(1/(1+q*(sigma*(1-lam)/(-quant))^p) | 1/p, q)
  );
  if (flip) out = 1 - out;
  return(out);
}


real sgt_lccdf(real x, real mu, real sigma, real lambdap1half, real p, real q){
  return(log(1-exp(sgt_lcdf(x | mu, sigma, lambdap1half, p, q))));
}

// Skew T distribution

real skew_t_lpdf(real x, real mu, real sigma, real lambdap1half, real q){
  return(sgt_lpdf(x | mu, sigma, lambdap1half, 2, q));
}

real skew_t_lcdf(real x, real mu, real sigma, real lambdap1half, real q){
  return(sgt_lcdf(x |  mu, sigma, lambdap1half, 2, q));
}

real skew_t_lccdf(real x, real mu, real sigma, real lambdap1half, real q){
  return(sgt_lccdf(x | mu, sigma, lambdap1half, 2, q));
}

real skew_t_rng(real mu, real sigma, real lambdap1half, real q){
  return(sgt_rng(mu, sigma, lambdap1half, 2, q));
}

// Constrained SGT distribution

real constrained_sgt_lpdf(real x, real mu, real sigma, real lambdap1half, real logpq, real logq){
  real q = exp(logq);
  real p = exp(log_sum_exp(logpq,0) - logq);
  return(sgt_lpdf(x | mu, sigma, lambdap1half, p, q));
}

real constrained_sgt_lcdf(real x, real mu, real sigma, real lambdap1half, real logpq, real logq){
  real q = exp(logq);
  real p = exp(log_sum_exp(logpq,0) - logq);
  return(sgt_lcdf(x |  mu, sigma, lambdap1half, p, q));
}

real constrained_sgt_lccdf(real x, real mu, real sigma, real lambdap1half, real logpq, real logq){
  real q = exp(logq);
  real p = exp(log_sum_exp(logpq,0) - logq);
  return(sgt_lccdf(x | mu, sigma, lambdap1half, p, q));
}

real constrained_sgt_rng(real mu, real sigma, real lambdap1half, real logpq, real logq){
  real q = exp(logq);
  real p = exp(log_sum_exp(logpq,0) - logq);
  return(sgt_rng(mu, sigma, lambdap1half, p, q));
}

// Constrained skew T distribution

real constrained_skew_t_lpdf(real x, real mu, real sigma, real lambdap1half, real qmhalf){
  real q = qmhalf + 1.0/2;
  return(sgt_lpdf(x | mu, sigma, lambdap1half, 2, q));
}

real constrained_skew_t_lcdf(real x, real mu, real sigma, real lambdap1half, real qmhalf){
  real q = qmhalf + 1.0/2;
  return(sgt_lcdf(x |  mu, sigma, lambdap1half, 2, q));
}

real constrained_skew_t_lccdf(real x, real mu, real sigma, real lambdap1half, real qmhalf){
  real q = qmhalf + 1.0/2;
  return(sgt_lccdf(x | mu, sigma, lambdap1half, 2, q));
}

real constrained_skew_t_rng(real mu, real sigma, real lambdap1half, real qmhalf){
  real q = qmhalf + 1.0/2;
  return(sgt_rng(mu, sigma, lambdap1half, 2, q));
}


// Symmetric GT distribution

real sym_gt_lpdf(real x, real mu, real sigma, real p, real q){
  return(sgt_lpdf(x | mu, sigma, 0.5, p, q));
}

real sym_gt_lcdf(real x, real mu, real sigma, real p, real q){
  return(sgt_lcdf(x |  mu, sigma, 0.5, p, q));
}

real sym_gt_lccdf(real x, real mu, real sigma, real p, real q){
  return(sgt_lccdf(x | mu, sigma, 0.5, p, q));
}

real sym_gt_rng(real mu, real sigma, real p, real q){
  return(sgt_rng(mu, sigma, 0.5, p, q));
}
