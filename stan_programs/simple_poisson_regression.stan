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
}
parameters {
  real alpha;
  real<upper=0> beta;
}
model {
  beta ~ normal(0, 0.5);
  alpha ~ normal(4, 2);
  
  qty_sld ~ poisson_log(alpha + beta * price);
} 
generated quantities {
  vector[N] pp_y;
  
  for (n in 1:N) 
    pp_y[n] = poisson_log_safe_rng(alpha + beta * price[n]);
}
