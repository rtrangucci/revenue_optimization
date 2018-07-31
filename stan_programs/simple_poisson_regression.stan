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
  int<lower=0> complaints[N];
  int<lower=0> traps[N];
}
transformed data {
  vector[N] traps_num = to_vector(traps);
}
parameters {
  real alpha;
  real beta;
}
model {
  beta ~ normal(-0.25, 1);
  alpha ~ normal(4, 2);
  complaints ~ poisson_log(alpha + beta * traps_num);
  /*
  for (n in 1:N) {
    complaints[n] ~ poisson(exp(alpha + beta * traps_num[n]));
  }
  */
} 
generated quantities {
  int y_rep[N];
  
  for (n in 1:N) 
    y_rep[n] = poisson_log_safe_rng(alpha + beta * traps_num[n]);
}
