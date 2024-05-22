library(tidyverse)

source('R/mcmc-sampling-fxns.R')


malaria_parms <- function(rnot = 1.1,
                          num_intros = 1,
                          distribution = "pois",
                          overdispersion=1,
                          date = NA,
                          rnot_dist=NA,
                          reporting_rate = 1,
                          secondary_trans = 0,
                          county_month = NA,
                          inform_prior=TRUE){
  return(as.list(environment()))
}


num_counties <- 2000
sec_trans <- rep(0, num_counties)
# sec_trans[3] <- 1
num_imports <- rep(3, num_counties)
cty_mth <- paste0('sarasota', seq(1,num_counties))
# rnot_dist_vec <- list(estimate=c(shape= 1.4, rate= 1))
rnot_dist_df <- tibble(shape = rep(1.4, num_counties), 
                        rate = rep(1, num_counties))
rgamma(10000, 1.4, 1) |> hist()

subs_parms(list(rnot = NA,
                rnot_dist = rnot_dist_df,
                num_intros = num_imports,
                distribution="nbinom",
                date='2023-08-15',
                secondary_trans = sec_trans,
                county_month = cty_mth),
           malaria_parms()) -> curr_parms

est_alphas <- mcmc_zika_rnot(zika_parms = curr_parms,
                             alpha_tuning = .1,
                             rnot_tuning = .1,
                             rr_tuning = .1,
                             burnin = 100000,
                             N = 400000,
                             thin=10)

est_alphas$samples |> head()

# est_alphas$samples[,1] |> plot(type='l')
est_alphas$samples[,2] |> plot(type='l') ## Alpha
# est_alphas$samples[,3] |> plot(type='l')
est_alphas$samples[,4] |> plot(type='l') ## Reproduction number
# est_alphas$samples[,5] |> plot(type='l') ## Reproduction number



est_alphas_df <- est_alphas$samples %>% as_tibble()
colnames(est_alphas_df) <- c('loglik', 'alpha', 'rr', paste0('R0_', 1:num_counties))

est_alphas_df |> 
  select(-loglik,-alpha,-rr) |> 
  gather(key, value) |> 
  group_by(key) |> 
  summarize(avg_r = mean(value),
            lb_r = quantile(value, probs = 0.025),
            ub_r = quantile(value, probs = 0.975)) |> 
  ggplot(aes(key, avg_r)) +
    geom_point() +
    geom_errorbar(aes(ymin = lb_r, ymax = ub_r))
