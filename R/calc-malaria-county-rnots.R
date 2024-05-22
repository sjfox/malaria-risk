library(sf)
library(tidyverse)
library(cowplot)
library(weathermetrics)
theme_set(theme_cowplot())

set.seed(450123)
# Read and clean the county temperature data ------------------------------
## This is using my temperature code
# county_month_temperature <- read_csv('processed-data/county-historic-temps.csv')

## Using Mathew's temperature data base

temperature_df <- read_csv('processed-data/county-average-temps.csv')
temperature_df |>
  filter(Year %in% c(2013:2022)) |>
  select(-Year) |>
  group_by(FIPS) |>
  summarize_all(.funs=mean) |>
  gather(month, avg_temp, -FIPS) |>
  rename(fips= FIPS) |>
  mutate(avg_temp = fahrenheit.to.celsius(avg_temp)) -> county_month_temperature
county_month_temperature

# Read in the economic - mosquito relationship ----------------------------
## This file from Downgrading disease transmission risk paper
## Gathers scam distribution from perkins et al paper
load('processed-data/parms_fxns_r0.RData')

## County per capita GDP info from
## US Census 2021
load('processed-data/county_gdp.rda')

# Read in mosquito abundance ----------------------------------------------
load('processed-data/county_mosq_prob_clean.rda')
## Abundance conversion from perkins et al
##  abundance = –ln(1 – occurrence probability)
get_mosq_abundance <- function(x){
  -log(1-x)
}
county_mosq_prob_clean |> 
  # select(-Anopheles_barberi,-Anopheles_bradleyi, -Anopheles_hermsi, -Anopheles_punctipennis, -Anopheles_crucians) |> 
  # mutate(Anopheles_crucians = ifelse(startsWith(fips, '12'), .95, 0)) |> 
  mutate(mosq_prob = 1-(apply((1-across(!starts_with("fips"))), 1, prod)),
         mosq_species_present = rowSums(ifelse(across(!starts_with("fips")) > 0.4, 1, 0))) |> 
  mutate(mosq_prob = ifelse(mosq_prob >.99, .99, mosq_prob)) |> 
  mutate(mosq_abundance_lambda = get_mosq_abundance(mosq_prob)) |>
  select(fips, mosq_prob, mosq_abundance_lambda, mosq_species_present) -> county_mosq_abundance
  # mutate_at(.vars = vars(-fips), .funs = get_mosq_abundance) |> 
  # mutate(mosq_abundance_lambda = rowSums(across(starts_with("2010"))),
         # mosq_species_present = rowSums(ifelse(across(starts_with("2010")) > 0.001, 1, 0))) |> 
  # select(fips, mosq_abundance_lambda, mosq_species_present) -> county_mosq_abundance

  
  
# Sampling from parameter distributions -----------------------------------
draw_param_vals <- function(param_list,
                           draw_sample = F,
                           num_samps = 1000){
  if(draw_sample){
    ## This is if you want to draw from parameter space
    map(1:num_samps, ~as_tibble(matrix(data = rnorm(length(param_list$parms), param_list$parms, param_list$parms_sd), nrow =1))) |> 
      bind_rows(.id = 'samp')
  } else{
    map(1:num_samps, ~as_tibble(matrix(data = param_list$parms, nrow =1))) |> 
      bind_rows(.id = 'samp')
  }
}

calc_briere_value <- function(temp, c, tm, t0){
  ## For Briere, parms are ordered as c, tm, t0
  result <- c * temp * (temp - t0) * sqrt(tm - temp )
  ifelse(is.nan(result) | result < 0, 0, result)
}

calc_quadratic_value <- function(temp, q, r, s){
  ## For quadratic, parms are ordered as q, r, s
  result <- q * temp^2 + r * temp + s
  ifelse(is.nan(result) | result < 0, 0, result)
}

# All Mosquito parameters --------------------------------------------
## Drawn from: https://onlinelibrary.wiley.com/doi/10.1111/ele.12015
## Mosquito biting rate (a)
br <- list(fxn = 'briere',
            ## For Briere, parms are ordered as c, tm, t0
            parms = c(0.000203, 42.3, 11.7),
            parms_sd = c(0.0000576, 3.53, 2.47))


## Daily adult survival probability (e^-mu) or e^-mortality_rate
dasp <- list(fxn = 'quadratic',
           ## For quadratic, parms are ordered as q, r, s
           parms = c(-0.000828, 0.0367, 0.522),
           parms_sd = c(0.0000519, 0.00239, 0.0235))

## parasite development rate (1/n) or (1/extrinsic incubation period)
pdr <- list(fxn = 'briere',
                    ## For Briere, parms are ordered as c, tm, t0
                    parms = c(0.000111, 34.4, 14.7),
                    parms_sd = c(0.0000161, 0.000176, 1.48))




# Human recovery rate -----------------------------------------------------
## Human recovery duration: 34.4 days (95% CI 22, 46 days)
## https://malariajournal.biomedcentral.com/articles/10.1186/s12936-015-0580-z
human_rr_mean <- 34.4
human_rr_sd <- (34.4-22)/1.96


# Mosquito to human transmission probability ------------------------------
## mosquito to human transmission probs 
## 9.2% (4.5%-16.0%) = https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5230737/

# Function to calculate the squared difference between observed quantiles and beta distribution quantiles
squared_difference <- function(params, mean_val, q_lower, q_upper) {
  alpha <- params[1]
  beta <- params[2]
  quantiles <- qbeta(c(0.025, 0.975), alpha, beta)
  diff <- c(mean_val - alpha/(alpha+beta), q_lower - quantiles[1], q_upper - quantiles[2])
  return(sum(diff^2))
}

# Estimate beta distribution parameters using optim
estimate_beta_parameters <- function(mean, q_lower, q_upper) {
  # Set initial parameter values
  initial_params <- c(1, 1)
  
  # Define the objective function to minimize
  objective <- function(params) squared_difference(params, mean_val, q_lower, q_upper)
  
  # Use optim to minimize the objective function
  result <- optim(initial_params, objective, method = "L-BFGS-B", lower = c(0, 0))
  
  # Extract the estimated parameters
  estimated_params <- result$par
  
  # Return the estimated parameters
  return(estimated_params)
}

# Input values
mean_val <- 0.092
q_lower <- 0.045
q_upper <- 0.160
##Estimate parameters
## Gives beta distribution parameters
mosq_human_prob_parms <- estimate_beta_parameters(mean_val, q_lower, q_upper)
mosq_human_prob_parms

## Confirm they look correct
# rbeta(10000, shape1 = mosq_human_prob_parms[1], shape2 = mosq_human_prob_parms[2]) ->test
# mean(test)
# quantile(test, probs = c(0.025, 0.975))
# hist(test)

# Human to mosquito transmission probability ------------------------------
## Assumes 20% human prevalence from this paper
## https://www.nature.com/articles/ncomms7054
## Averages distributions for two mosquitoes analyzed within (including parameter uncertainty)
## Drawn from supplementary table 1
## Assume that gambiae intercept has typo for upper bound (should be .29 instead of .19)

## Ordered, mean, lower, and upper of 95% credibility interval
get_mean_sd <- function(parms){
  mean_val <- parms[1]
  ## Standard deviation is average of either side
  sd_val <- mean(c((parms[1]-parms[2])/1.96,(parms[3]-parms[1])/1.96))
  return(c(mean_val, sd_val))
}
gambiae_intercept <- get_mean_sd(c(0.22, 0.11, 0.29))
gambiae_slope <- get_mean_sd(c(-0.27, -0.36, -0.18))
funestus_intercept <- get_mean_sd(c(0.15, 0.11, 0.19))
funestus_slope <- get_mean_sd(c(-0.13, -0.21, -0.06))

## Get normal distribution parameters for average distribution at 20% prevelance
((0.2*rnorm(10000, gambiae_slope[1], gambiae_slope[2]) + 
  rnorm(10000, gambiae_intercept[1], gambiae_intercept[2])) +
  (0.2*rnorm(10000, funestus_slope[1], funestus_slope[2]) + 
     rnorm(10000, funestus_intercept[1], funestus_intercept[2]))) / 2 -> human_mosq_probs

## Parameters of normal distribution for human to mosquito transmission probability
human_mosq_prob_parms <- c(mean(human_mosq_probs), sd(human_mosq_probs))
human_mosq_prob_parms

# Sample all parameters and calculate R0 -------------------------------
## Maximum of 1,000 samps based on length of scam list right now
num_samps <- 1000

draw_param_vals(br, draw_sample = T, num_samps=num_samps) |> 
  rename(br_c=V1, br_tm=V2, br_t0=V3) |> 
  left_join(draw_param_vals(dasp, draw_sample = F, num_samps=num_samps) |> ## Don't sample from vector competence, because uncertainty estimates too high
              rename(dasp_q=V1, dasp_r=V2, dasp_s=V3),
            by='samp') |> 
  left_join(draw_param_vals(pdr, draw_sample = T, num_samps=num_samps) |> 
              rename(pdr_c=V1, pdr_tm=V2, pdr_t0=V3),
            by='samp') |> 
  mutate(samp = as.integer(samp),
         human_rr = rnorm(n(), 
                          mean = human_rr_mean, 
                          sd = human_rr_sd),
         mosq_human_trans_prob = rbeta(n(), 
                                       shape1 = mosq_human_prob_parms[1], 
                                       shape2 = mosq_human_prob_parms[2]),
         human_mosq_trans_prob = rnorm(n(),
                                       mean = human_mosq_prob_parms[1],
                                       sd = human_mosq_prob_parms[2])) -> parm_samps

# Gathers samples for all necessary parameters and calculates R0
county_month_temperature |> 
  left_join(county_mosq_abundance, by = 'fips') |> 
  filter(!is.na(mosq_abundance_lambda)) |>
  left_join(county_gdp, by = 'fips') |> 
  filter(!is.na(gdp_percap)) |> 
  slice(rep(1:n(), each = num_samps)) |> 
  mutate(mosq_abundance_samp = rpois(n(), mosq_abundance_lambda)) |> 
  group_by(fips,month) |> 
  mutate(samp = seq_along(avg_temp)) |>
  group_by(samp) |> 
  mutate(gdp_scaling = exp(scam::predict.scam(scam.est.list[[unique(samp)]],
                                          newdata=data.frame(econ=log(gdp_percap))))) |> 
  ungroup() |> 
  left_join(parm_samps, by = 'samp') |> 
  mutate(biting_rate = calc_briere_value(avg_temp, br_c, br_tm, br_t0),
         adult_survival = calc_quadratic_value(avg_temp, dasp_q, dasp_r, dasp_s),
         parasite_development = calc_briere_value(avg_temp, pdr_c, pdr_tm, pdr_t0)) |>  
  mutate(rnot = (gdp_scaling * 
                   mosq_abundance_samp * 
                   mosq_human_trans_prob * 
                   human_mosq_trans_prob * 
                   biting_rate^2 *
                   exp( log(adult_survival) * 1/parasite_development)) / 
           (-log(adult_survival) * (1/human_rr))) -> county_parm_samples

## Some NAs produced frmo really low temperatures
county_parm_samples |>
  mutate(rnot = ifelse(is.na(rnot) & avg_temp < 1, 0, rnot)) -> county_parm_samples

county_parm_samples |> 
  group_by(fips, month) |> 
  summarize(rnot_mean = mean(rnot),
            rnot_lo = quantile(rnot,probs = 0.025),
            rnot_hi = quantile(rnot, probs = 0.975)) -> county_month_rnot_summary


save(county_parm_samples, county_month_rnot_summary, 
     file = 'processed-data/2023-08-09_county-monthly-rnot-estimates.rda')


# Convert to singular R0 value based on avg of 6 max months for each samp --------
# From perkins paper
# Because temperature values were available for each location on a monthly basis, 
## we computed monthly values of R0 for each location and then used the mean of the 
## highest six monthly values of R0 as a singular estimate of R0 for each location. 
## This approach was broadly consistent with the way in which a temperature suitability 
## index was used to inform mosquito occurrence probabilities by Kraemer and co-authors7.

# get_top_k_avg <- function(df, k){
#   sort(df$rnot, decreasing = T) |> 
#     head(k) |> 
#     mean()
# }

county_parm_samples |> 
  group_by(fips, month) |> 
  # mutate(id = seq_along(avg_temp)) |> 
  mutate(id=samp) |> 
  ungroup() |> 
  select(fips, month, id, rnot) |>
  # nest(data = c(month, rnot)) |> 
  # mutate(rnot = map(data, get_top_k_avg, k = 6) |> unlist()) -> county_single_rnot_samples
  group_by(id, fips) |> 
  summarize(rnot = mean(head(sort(rnot, decreasing=T), 6))) |> 
  ungroup() -> county_single_rnot_samples

county_single_rnot_samples |> 
  group_by(fips) |> 
  summarize(rnot_mean = mean(rnot),
            rnot_lo = quantile(rnot, probs = 0.025),
            rnot_hi = quantile(rnot, probs = 0.975)) -> county_single_rnot_summary


save(county_single_rnot_summary, county_single_rnot_samples, 
     file = 'processed-data/2023-08-09_county-single-rnot-estimates.rda')


# Add in code to create and save non-crucians rnots -----------------------






# Not used anymore --------------------------------------------------------

## Not needed with perkins et al formulation
# # Egg to adult survival probability
# etasp <- list(fxn = 'quadratic',
#                ## For quadratic, parms are ordered as q, r, s
#                parms = c(-0.00924, 0.453, -4.77),
#                parms_sd = c(0.00123, 0.0618, 0.746))
# 
# ## Mosquito development rate (1/larval development time)
# mdr <- list(fxn = 'briere',
#                     ## For Briere, parms are ordered as c, tm, t0
#                     parms = c(0.000111, 34, 14.7),
#                     parms_sd = c(0.00000954, 0.000106, 0.831))
# 
# # Eggs laid per adult female per day
# elpafpd <- list(fxn = 'quadratic',
#               ## For quadratic, parms are ordered as q, r, s
#               parms = c(-0.153, 8.61, -97.7),
#               parms_sd = c(0.0307, 1.69, 22.6))
# ## vector competence (b*c)
# vc <- list(fxn = 'quadratic',
#            ## For quadratic, parms are ordered as q, r, s
#            parms = c(-0.54, 25.2, -206),
#            parms_sd = c(0.18, 9.04, 108))
