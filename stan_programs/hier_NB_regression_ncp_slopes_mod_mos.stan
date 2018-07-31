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
  int<lower=1> M;
  int<lower=1> K;
  int complaints[N];
  vector[N] traps;
  int<lower=1> J;
  int<lower=1, upper=J> building_idx[N];
  matrix[J,K] building_data;
  vector[N] log_sq_foot;
  int<lower=1> mo_idx[N];
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
  real<lower=0> sigma_mo;
  vector[M] mo_raw;
  real<lower=0,upper=1> rho_raw;
}
transformed parameters {
  vector[J] mu = alpha + building_data * zeta + sigma_mu * mu_raw;
  vector[J] kappa = beta + building_data * gamma + sigma_kappa * kappa_raw;
  vector[M] mo = sigma_mo * mo_raw;
  real rho = 2.0 * rho_raw - 1.0;
  real phi = inv(inv_phi);
  mo[1] /= sqrt(1 - rho^2);
  for (m in 2:M) 
    mo[m] += rho * mo[m-1];
}
model {
  beta ~ normal(-0.25, 1);
  mu_raw ~ normal(0,1);
  kappa_raw ~ normal(0,1);
  mo_raw ~ normal(0,1);
  sigma_mu ~ normal(0, 1);
  sigma_kappa ~ normal(0, 1);
  sigma_mo ~ normal(0, 1);
  alpha ~ normal(0, 1);
  zeta ~ normal(0, 1);
  gamma ~ normal(0, 1);
  inv_phi ~ normal(0, 1);
  rho_raw ~ beta(10, 5);
  
  complaints ~ neg_binomial_2_log(mu[building_idx] + kappa[building_idx] .* traps 
                                 + mo[mo_idx] + log_sq_foot, phi);
} 
generated quantities {
  int y_rep[N];
  
  for (n in 1:N) 
    y_rep[n] = neg_binomial_2_log_safe_rng(mu[building_idx[n]] + kappa[building_idx[n]] * traps[n]
                                          + mo[mo_idx[n]] + log_sq_foot[n],
                                          phi);
}
