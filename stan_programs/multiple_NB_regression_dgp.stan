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
  real alpha = normal_rng(4, 2);
  real beta = -fabs(normal_rng(0, 0.5));
  real beta_super = normal_rng(1, 0.5);
  real beta_sq_foot = normal_rng(1, 0.5);
  real inv_prec = fabs(normal_rng(0,1));
  
  for (n in 1:N) 
    y_gen[n] = neg_binomial_2_log_rng(alpha + beta * traps[n] + beta_super * live_in_super[n]
                                       + beta_sq_foot * sq_foot[n], inv(inv_prec));
}
