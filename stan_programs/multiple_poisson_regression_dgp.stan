data {
  int<lower=1> N;
  vector[N] traps;
  vector[N] live_in_super;
  vector[N] sq_foot;
}
model {
} 
generated quantities {
  vector[N] y_gen;
  real alpha = normal_rng(0.5, 0.25);
  real beta = -fabs(normal_rng(0, 0.5));
  real beta_sq_foot = normal_rng(0, 0.5);
  real beta_super = normal_rng(-0.5, 0.5);
  
  for (n in 1:N) 
    y_gen[n] = poisson_log_rng(alpha + beta_sq_foot * sq_foot[n] 
                               + beta * traps[n] 
                               + beta_super * live_in_super[n]);
}
