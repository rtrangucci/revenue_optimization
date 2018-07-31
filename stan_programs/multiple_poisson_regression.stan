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
  vector[N] log_sq_foot;
  int complaints[N];
}
parameters {
  real alpha;
  real beta;
  real beta_super;
}
model {
  beta ~ normal(-0.25, 1);
  alpha ~ normal(0, 10);
  beta_super ~ normal(0, 0.5);
  
  complaints ~ poisson_log(alpha + beta * traps + beta_super * live_in_super + log_sq_foot);
} 
generated quantities {
  int y_rep[N];
  
  for (n in 1:N) 
    y_rep[n] = poisson_log_safe_rng(alpha + beta * traps[n] + beta_super * live_in_super[n]
                                   + log_sq_foot[n]);
}
