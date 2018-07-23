data {
  int<lower=1> N;
  real<lower=0> mean_traps;
}
model {
} 
generated quantities {
  int traps[N];
  int complaints[N];
  real alpha = normal_rng(4, 2);
  real beta = -fabs(normal_rng(0, 0.5));
  
  for (n in 1:N)  {
    traps[n] = poisson_rng(mean_traps);
    complaints[n] = poisson_log_rng(alpha + beta * traps[n]);
  }
}
