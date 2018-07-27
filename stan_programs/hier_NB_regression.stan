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
  int<lower=1> K;
  int complaints[N];
  vector[N] traps;
  int<lower=1> J;
  int<lower=1, upper=J> building_idx[N];
  matrix[J,K] building_data;
  vector[N] log_sq_foot;
}
parameters {
  real alpha;
  real<lower=0> sigma_alpha;
  vector[J] alphas;
  real beta;
  real<lower=0> inv_prec;
  vector[K] zeta;
}
transformed parameters {
  real prec = inv(inv_prec);
}
model {
  beta ~ normal(0, 1);
  zeta ~ normal(0, 1);
  alphas ~ normal(alpha + building_data * zeta, sigma_alpha);
  sigma_alpha ~ normal(0, 1);
  alpha ~ normal(0, 1);
  inv_prec ~ normal(0, 1);
  
  complaints ~ neg_binomial_2_log(alphas[building_idx] + beta * traps 
                               + log_sq_foot,
                               prec);
} 
generated quantities {
  int y_rep[N];
  
  for (n in 1:N) 
    y_rep[n] = neg_binomial_2_log_safe_rng(alphas[building_idx[n]] + beta * traps[n]
                                          + log_sq_foot[n],
                                          prec);
}
