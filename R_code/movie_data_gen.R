library(lubridate)
library(dplyr)
library(rstan)
library(zoo)

mod <- stan_model('neg_bi_2_rng.stan')
nms <- expose_stan_functions(mod)

cov_ar <- function(N, alpha, phi) {
  x <- 1:N
  dim_x <- N
  sq <- alpha^2 / (1 - phi^2)
  cov_ret <- matrix(NA_real_, dim_x, dim_x)
  for (i in 1:(dim_x - 1)) {
    cov_ret[i,i] <- sq 
    for (j in (i + 1):dim_x) {
      cov_ret[i, j] = sq * phi^(j - i)
      cov_ret[j,i] = cov_ret[i, j]
    }
  }
  cov_ret[dim_x, dim_x] = sq
  return(cov_ret)
}

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

gen_data <- function(N_movies, N_wks, obs_noise = 20, seed=123) {
  set.seed(seed)
  movie_inds <- as.vector(sapply(1:N_movies, rep, N_wks))
  ttl_data <- read.csv(file = '~/Downloads/imdb-250-2010-2009-2008.csv', sep=',',stringsAsFactors = F)[,c(1,4,5)] %>%
    rename(title = Title, release_year = year, rating = Rating) %>% mutate(movie = 1:250)
  
  wk_inds <- rep(1:N_wks, N_movies)
  dates <- as_date('2015-01-04') + (1:N_wks)*7
  n_price_changes <- rtpois(N_movies,3,1,5)
  starting_prices <- rtpois(N_movies,9,6,15)
  prices <- sapply(n_price_changes,function(x) sort(sample(2:N_wks,x, replace=F)))
  movie_prices <- c()
  phi_movies <- 0.5
  sig_movies <- 0.05
  movie_eta <- matrix(rt(N_wks * N_movies, df = 4),N_wks,N_movies)
  K <- cov_ar(N_wks, sig_movies, phi_movies)
  L_K <- t(chol(K))
  movie_series <- L_K %*% movie_eta
  movie_series <- as.vector(movie_series)
  for (n in 1:N_movies) {
    changes <- prices[[n]]
    starting_price <- starting_prices[n]
    level_changes <- 2*rbinom(length(changes), 1, 0.5) - 1
    levels <- cumsum(c(starting_price, level_changes))
    df_levels <- data.frame(price = levels, ind = c(1,changes))
    p_series <- rep(NA_real_,N_wks)
    p_series[df_levels$ind] <- df_levels$price
    p_series <- zoo::na.locf(p_series)
    movie_prices <- c(movie_prices, p_series)
  }
  overall_mu <- log(40)
  sd_alpha <- 0.05
  overall_beta <- -0.2
  sd_beta <- 0.05
  alphas <- overall_mu + sd_alpha * rnorm(N_movies)
  betas <- overall_beta + sd_beta * rnorm(N_movies)
  wk_noise <- rnorm(N_wks) * 0.10
  mus <- alphas[movie_inds] + betas[movie_inds] * (movie_prices - mean(movie_prices)) + wk_noise[wk_inds] + movie_series
  df_pred <- data.frame(mus = mus, movie = movie_inds, wk_ind = wk_inds, date = dates[wk_inds], price = movie_prices)
  df_pred <- df_pred %>% left_join(ttl_data, by = 'movie') %>%
    mutate(mus = mus + (release_year - mean(release_year)) * 0.05) %>%
    mutate(month = lubridate::month(date),
      mus = mus + ifelse(month != 12, sin(month) / 6, 0.25), 
      mus = mus + ifelse(date == '2015-12-27',0.25,0)) 
  for (n in 1:nrow(df_pred))
    df_pred$units_sold[n] <- neg_bi_rng(df_pred$mus[n], obs_noise)

  df_ret <- df_pred %>% select(units_sold, title, price, date, release_year, rating)
  return(df_ret)
}

df_test <- gen_data(100, 50, obs_noise = 1/2, seed=12)
summary(df_test$units_sold)
with(filter(df_test, title == unique(df_test$title)[100]), plot(date, units_sold, type = 'l'))
saveRDS(df_test, 'movie_data_tt.RDS')
