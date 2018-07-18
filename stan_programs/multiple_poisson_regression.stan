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
}
model {
  beta ~ normal(0, 0.5);
  alpha ~ normal(0, 10);
  beta_year ~ normal(1, 0.5);
  beta_rating ~ normal(1, 0.5);
  
  qty_sld ~ poisson_log(alpha + beta * price + beta_year * rel_year + beta_rating * rating);
} 
generated quantities {
  vector[N] pp_y;
  
  for (n in 1:N) 
    pp_y[n] = poisson_log_safe_rng(alpha + beta * price[n] + beta_year * rel_year[n]
                                   + beta_rating * rating[n]);
}
