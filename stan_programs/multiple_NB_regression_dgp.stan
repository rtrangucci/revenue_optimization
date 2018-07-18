data {
  int<lower=1> N;
  vector[N] price;
  vector[N] rel_year;
  vector[N] rating;
}
model {
} 
generated quantities {
  vector[N] y_gen;
  real alpha = normal_rng(4, 2);
  real beta = -fabs(normal_rng(0, 0.5));
  real beta_year = normal_rng(1, 0.5);
  real beta_rating = normal_rng(1, 0.5);
  real inv_prec = fabs(normal_rng(0,1));
  
  for (n in 1:N) 
    y_gen[n] = neg_binomial_2_log_rng(alpha + beta * price[n] + beta_year * rel_year[n]
                                       + beta_rating * rating[n], inv(inv_prec));
}
