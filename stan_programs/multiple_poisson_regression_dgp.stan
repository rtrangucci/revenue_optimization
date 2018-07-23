data {
  int<lower=1> N;
}
model {
} 
generated quantities {
  vector[N] sq_foot;
  int live_in_super[N];
  int traps[N];
  int complaints[N];
  real alpha = normal_rng(0.5, 0.25);
  real beta = -fabs(normal_rng(0, 0.5));
  real beta_sq_foot = normal_rng(0, 0.5);
  real beta_super = normal_rng(-0.5, 0.5);
  
  for (n in 1:N) {
    sq_foot[n] = normal_rng(1, 1);
    live_in_super[n] = bernoulli_rng(0.5);
    traps[n] = poisson_rng(8);
    complaints[n] = poisson_log_rng(alpha + beta_sq_foot * sq_foot[n] 
                               + beta * traps[n] 
                               + beta_super * live_in_super[n]);
  }
}
