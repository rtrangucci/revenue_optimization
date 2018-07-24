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
  matrix[J,K] meta_data;
  vector[N] sq_foot;
  vector[J] sq_foot_pred;
  int<lower=1> mo_idx[N];
}
transformed data {
  int N_hypo_traps = 20;
  int hypo_traps[N_hypo_traps];
  for (i in 1:N_hypo_traps)
    hypo_traps[i] = i;
}
parameters {
  real alpha;
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
  vector[J] std_alphas;
  vector[J] std_betas;
  real beta;
  real<lower=0> inv_prec;
  vector[K] zeta;
  vector[K] gamma;
  real<lower=0> sigma_mos;
  vector[M] std_mos;
}
transformed parameters {
  vector[J] alphas = alpha + meta_data * zeta + sigma_alpha * std_alphas;
  vector[J] betas = beta + meta_data * gamma + sigma_beta * std_betas;
  vector[M] mo = sigma_mos * std_mos;
  real prec = inv(inv_prec);
}
model {
  beta ~ normal(0, 1);
  std_alphas ~ normal(0,1) ;
  std_betas ~ normal(0,1) ;
  std_mos ~ normal(0,1) ;
  sigma_alpha ~ normal(0, 1);
  sigma_beta ~ normal(0, 1);
  sigma_mos ~ normal(0, 1);
  alpha ~ normal(0, 1);
  zeta ~ normal(0, 1);
  gamma ~ normal(0, 1);
  inv_prec ~ normal(0, 1);
  
  complaints ~ neg_binomial_2_log(alphas[building_idx] + betas[building_idx] .* traps 
                                 + mo[mo_idx] + sq_foot,
                               prec);
} 
generated quantities {
  matrix[J,N_hypo_traps] pp_y;
  matrix[J,N_hypo_traps] pp_rev;
  
  for (j in 1:J)
    for (i in 1:N_hypo_traps) {
      int ys[M];
      for (m in 1:M)
        ys[m] = neg_binomial_2_log_safe_rng(alphas[j] + betas[j] * hypo_traps[i]
                                                + mo[m] + sq_foot_pred[j],
                                                prec);
      pp_y[j,i] = sum(ys);
      pp_rev[j,i] = pp_y[j,i]/10.0 * (-100);
    }
}
