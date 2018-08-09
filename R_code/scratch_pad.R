library(lubridate)
library(dplyr)
library(rstan)
library(zoo)

tpois <- function(mu, lower, upper) {
  samp <- rpois(1, mu)
  while (samp < lower | samp > upper)
    samp <- rpois(1, mu)
  return(samp)
}

rtpois <- function(N, mu, lower, upper) {
  samps <- rep(NA_real_)
  for (n in 1:N)
    samps[n] <- tpois(mu, lower, upper)
  return(samps)
}

gen_data <- function(N_buildings, N_wks, obs_noise = 20, seed=123, n_changes = NULL, ar = NULL) {
  set.seed(seed)
  building_inds <- as.vector(sapply(1:N_buildings, rep, N_wks))
  building_data <- data.frame(building_id = 1:N_buildings,
                              floors = rtpois(N_buildings, 10, 3, 20),
                              sq_footage_p_floor = rnorm(N_buildings) * 7 * sqrt(3500) + 5000,
                              #sq_footage_p_floor = sample(c(465,672,1028),N_buildings,replace=T),
                              # live_in_super = sample(c(0,1),N_buildings, replace=T),
                              live_in_super = c(rep(0,ceiling(0.7*N_buildings)),
                                                rep(1,N_buildings - ceiling(0.7 * N_buildings)))[sample(1:N_buildings,N_buildings,replace=F)],
                              monthly_average_rent = rnorm(N_buildings) * 7 * sqrt(3500) + 3500,
                              average_tenant_age = rnorm(N_buildings) * 7 + 50,
                              age_of_building = rtpois(N_buildings, 50, 20, 70)) %>% 
    mutate(
      total_sq_foot = sq_footage_p_floor * floors
    )
  
  wk_inds <- rep(1:N_wks, N_buildings)
  dates <- as_date('2015-01-15') + (0:(N_wks-1))*30
  #n_trap_changes <- rtpois(N_buildings,5,4,10)
  n_trap_changes <- rtpois(N_buildings,20,15,30)
  if (!is.null(n_changes)) 
    n_trap_changes <- rep(n_changes, N_buildings)
  starting_traps <- rtpois(N_buildings,9,4,8)
  traps <- sapply(n_trap_changes,function(x) sort(sample(2:N_wks,x, replace=F)))
  building_traps <- c()
  phi_buildings <- 0.5
  sig_buildings <- 0.05
  for (n in 1:N_buildings) {
    changes <- traps[[n]]
    starting_trap <- starting_traps[n]
    level_changes <- 2*rbinom(length(changes), 1, 0.5) - 1
    levels <- cumsum(c(starting_trap, level_changes))
    df_levels <- data.frame(trap = levels, ind = c(1,changes))
    p_series <- rep(NA_real_,N_wks)
    p_series[df_levels$ind] <- df_levels$trap
    p_series <- zoo::na.locf(p_series)
    building_traps <- c(building_traps, p_series)
  }
  building_traps[building_traps < 0] = 0
  building_attr <- scale(building_data[1:N_buildings,c('live_in_super','age_of_building',
                                                       'monthly_average_rent',
                                                       'average_tenant_age')],scale=F)
  building_attr <- sweep(building_attr,MARGIN = 2, STATS = c(1,10,1e3,10),FUN='/')
  # Add AR 1
  ar_coef <- 0.90
  mo_1 <- rnorm(1)/sqrt(1 - ar_coef^2) * 0.55
  wk_noise <- rep(NA_real_,N_wks)
  wk_noise[1] <- mo_1
  for (n in 2:N_wks) {
    wk_noise[n] = wk_noise[n-1] * ar_coef + rnorm(1) * 0.55
  }
  if (!is.null(ar)) {
    wk_noise = ar
  }
  overall_mu <- log(2.5)
  sd_alpha <- 0.25
  overall_beta <- -0.2
  sd_beta <- 0.05
  gamma <- c(-0.15,0.15,0,0)
  alphas <- (overall_mu + sd_alpha * rnorm(N_buildings)) + log(building_data$total_sq_foot/1e4) + 
    building_attr %*% c(log(0.75),log(1.2),log(0.9),log(0.8))
  betas <- overall_beta + building_attr %*% gamma + sd_beta * rnorm(N_buildings)
  mus <- alphas[building_inds] + betas[building_inds] * building_traps + wk_noise[wk_inds]
  df_pred <- data.frame(mus = mus, building_id = building_inds, wk_ind = wk_inds, date = dates[wk_inds], trap = building_traps)
  df_pred <- df_pred %>% left_join(building_data, by = 'building_id') %>%
    mutate(month = lubridate::month(date))
  for (n in 1:nrow(df_pred))
    df_pred$complaints[n] <- rnbinom(n = 1, mu = exp(df_pred$mus[n]), size = obs_noise)
  
  df_ret <- df_pred %>% select(complaints, building_id, trap, date, 
                               live_in_super, age_of_building, total_sq_foot,
                               average_tenant_age, monthly_average_rent,
                               floors)
  return(list(df = df_ret,
              df_pred = df_pred,
              betas = betas,
              alphas = alphas,
              ar = wk_noise,
              n_trap_changes = n_trap_changes,
              building_data = building_attr,
              log_sq_foot_pred = building_data$total_sq_foot)) 
}

#set.seed(12)
# df_test <- gen_data(10, 12, obs_noise = 2, seed=123)
# ar_vec <- c(df_test$ar[1:10],df_test$ar[4:3])
# df_test_2 <- gen_data(10, 12, obs_noise = 2, seed=123, n_changes = NULL, ar = ar_vec*2)
df_test <- gen_data(10, 36, obs_noise = 10, seed=123)
#save(df_test$df)
saveRDS(df_test$df, file = 'data/pest_data_longer.RDS')

stan_dat <- with(df_test$df,
                 list(
                   complaints = complaints,
                   traps = trap,
                   pred_mat = cbind(live_in_super, age_of_building/10,
                                    average_tenant_age/10, monthly_average_rent/1000),
                   N = length(complaints),
                   K = 4,
                   log_sq_foot = log(total_sq_foot/1e4),
                   mo_idx = as.integer(as.factor(date)),
                   M = 36,
                   J = length(unique(building_id)),
                   building_idx = building_id,
                   log_sq_foot_pred = log(df_test$log_sq_foot_pred/1e4),
                   M_forward = 12
                 ))
stan_dat$building_data <- df_test$building_data
test_mod_comp <- stan_model('stan_programs/test_reg.stan')
test_mod_fit <- sampling(test_mod_comp, cores = 4, chains = 4, data = stan_dat)
mcmc_recover_hist(as.matrix(test_mod_fit, pars = c('zeta', 'alpha', 'beta', 'phi')),
                  true = c(log(0.75),log(1.2),log(0.8),log(0.9), log(2.5), -0.21, 2))

test_mod_comp <- stan_model('stan_programs/test_reg_AR.stan')
test_mod_fit <- sampling(test_mod_comp, cores = 4, chains = 4, data = stan_dat)
mcmc_recover_hist(as.matrix(test_mod_fit, pars = c('zeta', 'alpha', 'beta', 'phi','rho','sigma_mo')),
                  true = c(log(0.75),log(1.2),log(0.8),log(0.9), log(2.5), -0.21, 2,0.95,0.55))

test_mod_comp <- stan_model('stan_programs/test_reg_AR_hier.stan')
test_mod_fit <- sampling(test_mod_comp, data = stan_dat, chains = 3, cores = 3, iter = 1000)
mcmc_recover_hist(as.matrix(test_mod_fit, pars = c('zeta', 'alpha', 'beta', 'phi','rho','sigma_mo','sigma_mu')),
                  true = c(log(0.75),log(1.2),log(0.8),log(0.9), log(2.5), -0.21, 2,0.6,0.55,0.25))
# with building attr
mcmc_recover_hist(as.matrix(test_mod_fit, pars = c('zeta', 'alpha', 'beta', 'phi','rho','sigma_mo','sigma_mu','gamma','sigma_kappa')),
                  true = c(log(0.75),log(1.2),log(0.9),log(0.8), log(2.5), -0.2, 10,0.9,0.55,0.25,-0.15,0.15,0.05))
saveRDS(stan_dat, 'data/pest_data_longer_stan_dat.RDS')
plot(colMeans(rstan::extract(test_mod_fit)$mo),df_test$ar)
abline(0,1)

plot(colMeans(rstan::extract(test_mod_fit)$kappa),df_test$betas)
abline(0,1)

plot(colMeans(rstan::extract(test_mod_fit)$mu),df_test$mu)
abline(0,1)


gen_beta <- function(x) {
  alpha <- 10
  bet <- 5
  a <- -1
  c <- 1
  lp <- (alpha - 1) * log(x - a) + (bet - 1) * log(c - x) - (alpha + bet - 1) * log(c - a) - lbeta(alpha, bet)
  return(exp(lp))
}
mcmc_hist(
  as.matrix(test_mod_fit, pars = "rho"), freq=F,
  binwidth = 0.01
) + overlay_function(fun = gen_beta)
print(test_mod_fit)
plot(df_test_2$df$trap,log(df_test_2$df$complaints/df_test_2$df$total_sq_foot*1e4))
hist(df_test$alphas)
hist(df_test$betas)
plot(df_test$ar,type='l')
plot(df_test_2$ar,type='l')
group_by(df_test_2$df,date) %>% summarise(m_c = mean(complaints)) %>% plot(type='l')

tt <- glm(complaints ~ offset(log(total_sq_foot/1e4)) + trap, data = df_test_2$df, family = quasipoisson)
summary(tt)
df_test_2$df$resid <- resid(lm(log1p(complaints) ~ offset(log(total_sq_foot/1e4)) + trap, data = df_test_2$df))
group_by(df_test_2$df,date) %>% summarise(m_r = mean(resid)) %>% plot(type='l')
summary(tt)

tt <- glm(complaints ~ trap, data = df_test_2$df, family = quasipoisson)
summary(tt)
hist(df_test_2$n_trap_changes)

test <- readRDS('data/building_data_20180727_time_trend.RDS')
df_test_2$df_pred$building_id <- unique(test$building_id)[df_test_2$df_pred$building_id]

#hist(df_test$alphas)
df_test_2 <- gen_data(100, 50, obs_noise = 12, seed=12)
df <- df_test_2$df
df %>% mutate(idx = as.integer(as.factor(title)), opt_price = -1/df_test_2$betas[idx]) -> df_mod
unique(df_mod[,c('title','opt_price')])
summary(df$units_sold)
with(dplyr::filter(df_test$df, building_id == 1), plot(date, complaints, type = 'l'))
#saveRDS(rename(df_test$df,traps=trap), 'data/building_data_20180727.RDS')
saveRDS(rename(df_test_2$df_pred,traps=trap), 'data/building_data_20180803_alt.RDS')


tt <- glm(complaints ~ offset(log(total_sq_foot/1e4)) + trap, data = test, family = quasipoisson)
summary(tt)

group_by(test,date) %>% summarise(m_c = mean(complaints)) %>% plot(type='l')

  phi_buildings <- 0.5
  sig_buildings <- 0.05
  for (n in 1:N_buildings) {
    changes <- traps[[n]]
    starting_trap <- starting_traps[n]
    level_changes <- 2*rbinom(length(changes), 1, 0.5) - 1
    levels <- cumsum(c(starting_trap, level_changes))
    df_levels <- data.frame(trap = levels, ind = c(1,changes))
    p_series <- rep(NA_real_,N_wks)
    p_series[df_levels$ind] <- df_levels$trap
    p_series <- zoo::na.locf(p_series)
    building_traps <- c(building_traps, p_series)
  }
  building_attr <- scale(building_data[1:N_buildings,c('live_in_super','age_of_building',
                                                             'monthly_average_rent',
                                                             'average_tenant_age')],scale=F)
  building_attr <- sweep(building_attr,MARGIN = 2, STATS = c(1,10,1e3,10),FUN='/')
  # Add AR 1
  mo_1 <- rnorm(1) * 0.15 - 0.2
  ar_coef <- 0.95
  wk_noise <- rep(NA_real_,N_wks)
  wk_noise[1] <- mo_1
  for (n in 2:N_wks) {
    wk_noise[n] = wk_noise[n-1] * ar_coef + rnorm(1) * 0.15
  }
  if (!is.null(ar)) {
    wk_noise = ar
  }
  overall_mu <- log(2.5)
  sd_alpha <- 0.25
  overall_beta <- 0.75
  sd_beta <- 0.05
  gamma <- c(0.15,-0.15,0,0)
  alphas <- (overall_mu + sd_alpha * rnorm(N_buildings)) + log(building_data$total_sq_foot/1e4) + 
    building_attr %*% c(log(0.75),log(1.2),log(0.9),log(0.8))
  betas <- -exp(overall_beta + building_attr %*% gamma + sd_beta * rnorm(N_buildings))/10
  mus <- alphas[building_inds] + betas[building_inds] * building_traps + wk_noise[wk_inds]
  df_pred <- data.frame(mus = mus, building_id = building_inds, wk_ind = wk_inds, date = dates[wk_inds], trap = building_traps)
  df_pred <- df_pred %>% left_join(building_data, by = 'building_id') %>%
    mutate(month = lubridate::month(date))
  for (n in 1:nrow(df_pred))
    df_pred$complaints[n] <- rnbinom(n = 1, mu = exp(df_pred$mus[n]), size = obs_noise)

  df_ret <- df_pred %>% select(complaints, building_id, trap, date, 
                               live_in_super, age_of_building, total_sq_foot,
                               average_tenant_age, monthly_average_rent,
                               floors)
  return(list(df = df_ret,
              df_pred = df_pred,
              betas = betas,
              alphas = alphas,
              ar = wk_noise,
              n_trap_changes = n_trap_changes)) 
}

#set.seed(12)
df_test <- gen_data(10, 24, obs_noise = 2, seed=123)
ar_vec <- c(df_test$ar[1:10], df_test$ar[4:3])
df_test_2 <- gen_data(10, 24, obs_noise = 2, seed=123, n_changes = NULL, ar = ar_vec*2)
plot(df_test_2$df$trap,log(df_test_2$df$complaints/df_test_2$df$total_sq_foot*1e4))
hist(df_test$alphas)
hist(df_test$betas)
plot(df_test$ar,type='l')
plot(df_test_2$ar,type='l')
group_by(df_test_2$df,date) %>% summarise(m_c = mean(complaints)) %>% plot(type='l')

tt <- glm(complaints ~ offset(log(total_sq_foot/1e4)) + trap, data = df_test_2$df, family = quasipoisson)
summary(tt)
df_test_2$df$resid <- resid(lm(log1p(complaints) ~ offset(log(total_sq_foot/1e4)) + trap, data = df_test_2$df))
group_by(df_test_2$df,date) %>% summarise(m_r = mean(resid)) %>% plot(type='l')
summary(tt)

tt <- glm(complaints ~ trap, data = df_test_2$df, family = quasipoisson)
summary(tt)
hist(df_test_2$n_trap_changes)

test <- readRDS('data/building_data_20180727_time_trend.RDS')
df_test_2$df_pred$building_id <- unique(test$building_id)[df_test_2$df_pred$building_id]

#hist(df_test$alphas)
df_test_2 <- gen_data(100, 50, obs_noise = 12, seed=12)
df <- df_test_2$df
df %>% mutate(idx = as.integer(as.factor(title)), opt_price = -1/df_test_2$betas[idx]) -> df_mod
unique(df_mod[,c('title','opt_price')])
summary(df$units_sold)
with(dplyr::filter(df_test$df, building_id == 1), plot(date, complaints, type = 'l'))
#saveRDS(rename(df_test$df,traps=trap), 'data/building_data_20180727.RDS')

df_test$df_pred$trap[df_test$df_pred$trap < 0] <- 0
df_test_2$df_pred$trap[df_test_2$df_pred$trap < 0] <- 0
saveRDS(rename(df_test_2$df_pred,traps=trap), 'data/building_data_20180807_alt.RDS')


tt <- glm(complaints ~ offset(log(total_sq_foot/1e4)) + trap, data = test, family = quasipoisson)
summary(tt)

group_by(test,date) %>% summarise(m_c = mean(complaints)) %>% plot(type='l')
