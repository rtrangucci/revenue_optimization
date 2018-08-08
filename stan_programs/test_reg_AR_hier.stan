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
  matrix[N,K] pred_mat;
  vector[N] log_sq_foot;
  int<lower=1> mo_idx[N];
  int<lower=1> J;
  int<lower=1> building_idx[N];
  matrix[J,K] building_attr;
}
parameters {
  real alpha;
  real beta;
  real<lower=0> inv_phi;
  vector[K] zeta;
  real<lower=0> sigma_mo;
  vector[M] mo_raw;
  real<lower=0,upper=1> rho_raw;
  real<lower=0> sigma_mu;
  vector[J] mu_raw;
}
transformed parameters {
  real phi = inv(inv_phi);
  vector[J] mu = building_attr * zeta + mu_raw * sigma_mu;
  vector[M] mo = sigma_mo * mo_raw;
  real rho = 2.0 * rho_raw - 1.0;
  mo[1] /= sqrt(1 - rho^2);
  for (m in 2:M) 
    mo[m] += rho * mo[m-1];
}
model {
  beta ~ normal(-0.25, 1);
  zeta ~ normal(0, 1);
  alpha ~ normal(0, 10);
  inv_phi ~ normal(0, 1);
  sigma_mo ~ normal(0, 1);
  mo_raw ~ normal(0,1);
  rho_raw ~ beta(10, 5);
  mu_raw ~ normal(0, 1);
  sigma_mu ~ normal(0, 1);
  
  complaints ~ neg_binomial_2_log(alpha + beta * traps + mo[mo_idx] 
                                 + mu[building_idx]
                               + log_sq_foot,
                               phi);
} 
