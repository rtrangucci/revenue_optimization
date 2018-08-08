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
  matrix[N,K] pred_mat;
  vector[N] log_sq_foot;
}
parameters {
  real alpha;
  real beta;
  real<lower=0> inv_phi;
  vector[K] zeta;
}
transformed parameters {
  real phi = inv(inv_phi);
}
model {
  beta ~ normal(-0.25, 1);
  zeta ~ normal(0, 1);
  alpha ~ normal(0, 1);
  inv_phi ~ normal(0, 1);
  
  complaints ~ neg_binomial_2_log(alpha + pred_mat * zeta + beta * traps 
                               + log_sq_foot,
                               phi);
} 
