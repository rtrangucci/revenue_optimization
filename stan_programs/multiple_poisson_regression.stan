functions {
  int poisson_log_safe_rng(real eta) {
    real pois_rate = exp(eta);
    if (pois_rate >= exp(20.79))
      return -9;
    return poisson_rng(pois_rate);
  }
}
data {
  int<lower=1> N;
  vector[N] traps;
  vector[N] live_in_super;
  vector[N] sq_foot;
  int complaints[N];
}
parameters {
  real alpha;
  real<upper=0> beta;
  real beta_super;
  real beta_sq_foot;
}
model {
  beta ~ normal(0, 0.5);
  alpha ~ normal(0, 10);
  beta_super ~ normal(0, 0.5);
  beta_sq_foot ~ normal(0, 0.5);
  
  complaints ~ poisson_log(alpha + beta * traps + beta_super * live_in_super + beta_sq_foot * sq_foot);
} 
generated quantities {
  vector[N] pp_y;
  
  for (n in 1:N) 
    pp_y[n] = poisson_log_safe_rng(alpha + beta * traps[n] + beta_super * live_in_super[n]
                                   + beta_sq_foot * sq_foot[n]);
}
