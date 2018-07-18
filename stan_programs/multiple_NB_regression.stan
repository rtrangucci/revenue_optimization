functions {
  int neg_binomial_2_log_safe_rng(real eta, real phi) {
    real phi_div_exp_eta;
    real gamma_rate;
    phi_div_exp_eta = phi/exp(eta);
    gamma_rate = gamma_rng(phi, phi_div_exp_eta);
    if (gamma_rate >= exp(20.79))
      return -9;
    return poisson_rng(gamma_rate);
  }
}
data {
  int<lower=1> N;
  int qty_sld[N];
  vector[N] price;
  vector[N] rel_year;
  vector[N] rating;
}
parameters {
  real alpha;
  real<upper=0> beta;
  real beta_year;
  real beta_rating;
  real<lower=0> inv_phi;
}
transformed parameters {
  real phi = inv(inv_phi);
}
model {
  beta ~ normal(0, 0.5);
  alpha ~ normal(0, 10);
  beta_year ~ normal(1, 0.5);
  beta_rating ~ normal(1, 0.5);
  inv_phi ~ normal(0, 1);
  
  qty_sld ~ neg_binomial_2_log(alpha + beta * price + beta_year * rel_year
                               + beta_rating * rating,phi);
} 
generated quantities {
  vector[N] pp_y;
  
  for (n in 1:N) 
    pp_y[n] = neg_binomial_2_log_safe_rng(alpha + beta * price[n] + beta_year * rel_year[n]
                                          + beta_rating * rating[n], phi);
}

