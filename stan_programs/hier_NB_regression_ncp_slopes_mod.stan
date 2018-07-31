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
  real<lower=0> sigma_mu;
  real<lower=0> sigma_kappa;
  vector[J] mu_raw;
  vector[J] kappa_raw;
  real beta;
  real<lower=0> inv_phi;
  vector[K] zeta;
  vector[K] gamma;
}
transformed parameters {
  vector[J] mu = alpha + building_data * zeta + sigma_mu * mu_raw;
  vector[J] kappa = beta + building_data * gamma + sigma_kappa * kappa_raw;
  real phi = inv(inv_phi);
}
model {
  beta ~ normal(-0.25, 1);
  mu_raw ~ normal(0,1) ;
  kappa_raw ~ normal(0,1) ;
  sigma_mu ~ normal(0, 1);
  sigma_kappa ~ normal(0, 1);
  alpha ~ normal(0, 1);
  zeta ~ normal(0, 1);
  gamma ~ normal(0, 1);
  inv_phi ~ normal(0, 1);
  
  complaints ~ neg_binomial_2_log(mu[building_idx] + kappa[building_idx] .* traps 
                               + log_sq_foot,
                               phi);
} 
generated quantities {
  int y_rep[N];
  
  for (n in 1:N) 
    y_rep[n] = neg_binomial_2_log_safe_rng(mu[building_idx[n]] + kappa[building_idx[n]] * traps[n]
                                          + log_sq_foot[n],
                                          phi);
}
