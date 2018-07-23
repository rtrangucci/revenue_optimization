data {
  int<lower=1> N;
  vector[N] traps;
}
model {
} 
generated quantities {
  vector[N] y_gen;
  real alpha = normal_rng(4, 2);
  real beta = -fabs(normal_rng(0, 0.5));
  
  for (n in 1:N) 
    y_gen[n] = poisson_log_rng(alpha + beta * traps[n]);
}
