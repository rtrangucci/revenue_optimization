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
  int complaints[N];
  vector[N] traps;
}
parameters {
  real alpha;
  real<upper=0> beta;
}
model {
  beta ~ normal(0, 0.5);
  alpha ~ normal(4, 2);
  
  complaints ~ poisson_log(alpha + beta * traps);
} 
generated quantities {
  vector[N] pp_y;
  
  for (n in 1:N) 
    pp_y[n] = poisson_log_safe_rng(alpha + beta * traps[n]);
}
