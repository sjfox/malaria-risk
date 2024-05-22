## Script to run malaria estimation
library(MASS)
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

get_gamma_parms <- function(rnots){
  rnots <- unlist(rnots)
  fit <- try(MASS::fitdistr(as.numeric(unlist(rnots)), "gamma"), silent = T)
  if(class(fit) == "try-error"){
    # browser()
    fit <- list(estimate=c(shape=mean(as.numeric(unlist(rnots))), rate=1))
  } 
  tibble(shape = fit$estimate["shape"], 
         rate = fit$estimate["rate"])
}

## Process the county prior r0 distributions into gamma functions
load('processed-data/2023-08-09_county-monthly-rnot-estimates.rda')
county_parm_samples |> 
  select(fips,month,rnot) |> 
  nest(rnots = rnot) |> 
  mutate(gamma_parms = map(rnots, get_gamma_parms)) |> 
  unnest(gamma_parms) |> 
  select(-rnots) |> 
  mutate(county_month = paste0(fips, '-',month)) -> county_month_rnot_priors

## Estimate importations for each county and month
read_and_clean_annual_imports <- function(file_path){
  read_csv(file_path, 
           col_select = c(-1,-3,-4),
           col_types = c('ccnnnnnnnnnnnnnnn')) |> 
    mutate(fips = str_pad(FIPS,width=5, pad = '0')) |> 
    select(-FIPS) |> 
    gather(month, imports, -fips) |> 
    mutate(month = month.abb[as.numeric(str_remove(month, 'Est Imports Month '))])
}
set.seed(34542)
list.files('processed-data/annual-county-imports', full.names = T) |> 
  map(read_and_clean_annual_imports) |> 
  bind_rows() |> 
  group_by(fips, month) |> 
  summarize(imports = sum(imports)) |> 
  mutate(imports = rpois(n(), lambda = imports)) -> county_month_imports

## Gather secondary transmission for each county and month
tibble::tribble(
                                                                                                    ~county,     ~state,   ~fips,        ~date, ~single_case, ~all_cases,
                                                                                                  "Cameron",    "Texas", "48061", "2023-06-23",           1L,         1L,
                                                                                                 "Sarasota",  "Florida", "12115", "2023-06-26",           1L,         7L,
   "Prince George's (It just says nations capital region so I chose the highest import county in maryland)", "Maryland", "24033", "2023-08-18",           1L,         1L,
                                                                                                   "Saline", "Arkansas", "05125", "2023-10-04",           1L,         1L
  ) |> 
  mutate(month = month.abb[month(ymd(date))]) |> 
  select(fips, month, single_case) -> sec_trans_df

## Join into master tibble for running estimation
county_month_rnot_priors |> 
  left_join(county_month_imports, by = c('fips', 'month')) |> 
  left_join(sec_trans_df, by = c('fips', 'month')) |> 
  mutate(imports = ifelse(is.na(imports), 0 , imports),
         single_case = ifelse(is.na(single_case), 0 , single_case)) -> county_month_full_df

## Only need counties with imports (no value added to no imports)
## Right now cheating with imports numbers
county_month_full_df |> 
  mutate(imports = ifelse(imports ==0 & single_case == 1, 1, imports)) |> 
  filter(imports !=0, shape != 0) -> county_month_reduced_df


subs_parms(list(rnot = NA,
                rnot_dist = county_month_reduced_df |> select(shape,rate),
                num_intros = county_month_reduced_df$imports,
                distribution = "nbinom",
                date = '2024-01-01',
                secondary_trans = county_month_reduced_df$single_case,
                county_month = county_month_reduced_df$county_month),
           malaria_parms()) -> curr_parms

est_alphas <- mcmc_zika_rnot(zika_parms = curr_parms,
                             alpha_tuning = .05,
                             rnot_tuning = .05,
                             rr_tuning = .05,
                             burnin = 100000,
                             N = 500000,
                             thin=40)

est_alphas$samples[,2] |> plot(type='l') 
est_alphas$samples[,2] |> summary() 

# est_alphas$samples |> as_tibble() |> 
#   select(-1, -2, -3) |> 
#   gather(key, value) |> 
#   group_by(key) |> 
#   summarize(avg_val = mean(value),
#             hi_rnot = quantile(value, probs = 0.975)) |> 
#   arrange(desc(hi_rnot)) -> temp
# temp
# dim(est_alphas$samples)
# length(curr_parms$county_month)


county_rnot_samps <- est_alphas$samples |> as_tibble() |> 
  select(-1, -2, -3)
colnames(county_rnot_samps) <- curr_parms$county_month
county_rnot_samps |> 
  gather(key, rnot) |> 
  separate(key, into = c('zip', 'month'), sep = '-') -> county_rnot_samps


est_alphas$samples |> as_tibble() |> 
  select(2) |> 
  rename(alpha = V2) |> 
  mutate(samp = seq_along(alpha)) -> alpha_samps

est_alphas$samples |> as_tibble() |> 
  select(1) |> 
  rename(loglik = V1) |> 
  mutate(samp = seq_along(loglik)) -> ll_samps

est_alphas$samples |> as_tibble() |> 
  select(3) |> 
  rename(rr = V3) |> 
  mutate(samp = seq_along(rr)) -> rr_samps



save(rr_samps, ll_samps, alpha_samps, county_rnot_samps, file = 'processed-data/malaria-risk-mcmc-samples.rda')



