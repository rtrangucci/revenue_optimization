functions {
  int neg_binomial_2_log_safe_rng(real eta, real phi) {
    real phi_div_exp_eta;
    real gamma_rate;
    phi_div_exp_eta = phi/exp(eta);
    gamma_rate = gamma_rng(phi, phi_div_exp_eta);
    if (gamma_rate >= exp(20.79))
      return -9;
    return poisson_rng(gamma_rate);
  }
}
data {
  int<lower=1> N;
  vector[N] traps;
  vector[N] live_in_super;
  vector[N] log_sq_foot;
  int complaints[N];
}
parameters {
  real alpha;
  real<upper=0> beta;
  real beta_sq_foot;
  real beta_super;
 real<lower=0> inv_prec;
}
transformed parameters {
  real prec = inv(inv_prec);
}
model {
  beta ~ normal(0, 0.5);
  alpha ~ normal(0, 1);
  beta_super ~ normal(0, 1);
  beta_sq_foot ~ normal(0, 1);
  inv_prec ~ normal(0, 1);
  
  complaints ~ neg_binomial_2_log(alpha + beta * traps + beta_super * live_in_super
                               + beta_sq_foot * log_sq_foot, prec);
} 
generated quantities {
  int y_rep[N];
  
  for (n in 1:N) 
    y_rep[n] = neg_binomial_2_log_rng(alpha + beta * traps[n] + beta_super * live_in_super[n]
                                       + beta_sq_foot * log_sq_foot[n], prec);
  
}

