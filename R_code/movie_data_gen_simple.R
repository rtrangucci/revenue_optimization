library(lubridate)
library(dplyr)
library(rstan)
library(zoo)

mod <- stan_model('neg_bi_2_rng.stan')
nms <- expose_stan_functions(mod)

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
    rename(title = Title, release_year = year, rating = Rating) %>% arrange(title) %>% mutate(movie = 1:250)
  
  wk_inds <- rep(1:N_wks, N_movies)
  dates <- as_date('2015-01-04') + (1:N_wks)*7
  n_price_changes <- rtpois(N_movies,3,1,5)
  starting_prices <- rtpois(N_movies,9,6,15)
  prices <- sapply(n_price_changes,function(x) sort(sample(2:N_wks,x, replace=F)))
  movie_prices <- c()
  phi_movies <- 0.5
  sig_movies <- 0.05
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
  release_dates <- scale(ttl_data[1:N_movies,c('release_year','rating')])
  overall_mu <- log(40)
  sd_alpha <- 0.05
  overall_beta <- 0.75
  sd_beta <- 0.05
  gamma <- c(-0.15,-0.15)
  alphas <- overall_mu + sd_alpha * rnorm(N_movies)
  betas <- -exp(overall_beta + release_dates %*% gamma + sd_beta * rnorm(N_movies))/10
  wk_noise <- rnorm(N_wks) * 0.10
  mus <- alphas[movie_inds] + betas[movie_inds] * (movie_prices - mean(movie_prices)) + wk_noise[wk_inds]# + movie_series
  df_pred <- data.frame(mus = mus, movie = movie_inds, wk_ind = wk_inds, date = dates[wk_inds], price = movie_prices)
  df_pred <- df_pred %>% left_join(ttl_data, by = 'movie') %>%
    mutate(mus = mus + (release_year - mean(release_year)) * 0.05) %>%
    mutate(month = lubridate::month(date))#,
 #     mus = mus + ifelse(month != 12, sin(month) / 6, 0.25), 
 #     mus = mus + ifelse(date == '2015-12-27',0.25,0)) 
  for (n in 1:nrow(df_pred))
    df_pred$units_sold[n] <- neg_bi_rng(df_pred$mus[n], obs_noise)

  df_ret <- df_pred %>% select(units_sold, title, price, date, release_year, rating)
  return(list(df = df_ret,
              betas = betas)) 
}

df_test <- gen_data(100, 50, obs_noise = 2, seed=12)
df_test_2 <- gen_data(100, 50, obs_noise = 12, seed=12)
df <- df_test_2$df
df %>% mutate(idx = as.integer(as.factor(title)), opt_price = -1/df_test_2$betas[idx]) -> df_mod
unique(df_mod[,c('title','opt_price')])
summary(df$units_sold)
with(filter(df, title == unique(df$title)[1]), plot(date, units_sold, type = 'l'))
saveRDS(df_test_2, 'movie_data_tt_new_set_20180717.RDS')
