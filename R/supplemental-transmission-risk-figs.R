## Supplemental figures for the transmission rate
library(tidyverse)
library(sf)
library(cowplot)
theme_set(theme_cowplot())
load('processed-data/county_sf.rda')


# Create figure showing parameters for transmission risks ----------------
load('processed-data/parms_fxns_r0.RData')
load('processed-data/county_gdp.rda')


get_gdp <- function(scam_est, gdp_vec){
  tibble(gdp = gdp_vec,
         gdp_multiplier = exp(scam::predict.scam(scam_est,
                                newdata=data.frame(econ=log(gdp_vec)))))
}

gdp_vec <- seq(1000, max(county_gdp$gdp_percap), length.out = 1000)
scam.est.list |> 
  map(get_gdp, gdp_vec = gdp_vec) |> 
  bind_rows(.id = 'id') |> 
  ggplot(aes(log(gdp), gdp_multiplier, group = id)) +
    geom_line(alpha = .1) + 
    background_grid(major = 'xy') +
    labs(x = 'log(GDP)', y = expression(italic(R[0])~Scale~Factor)) -> scale_factor_plot
scale_factor_plot

## Sampling from parameter distributions
draw_param_val <- function(temp, ## Temperature value of interest
                           param_list,
                           draw_sample = F){
  if(draw_sample){
    ## This is if you want to draw from parameter space
    vals <- rnorm(length(param_list$parms),
                  param_list$parms,
                  param_list$parms_sd)  
  } else{
    vals <- param_list$parms
  }
  if(param_list$fxn == 'briere'){
    ## For Briere, parms are ordered as c, tm, t0
    result <- vals[1] * temp * (temp - vals[3]) * sqrt(vals[2] - temp )
    if(is.nan(result)) {
      return(0)
    }else if(result<0){
      return(0)
    } else{
      return(result)
    }
  } else if(param_list$fxn == 'quadratic'){
    ## For quadratic, parms are ordered as q, r, s
    result <- vals[1] * temp^2 + vals[2] * temp + vals[3]
    if(is.nan(result)) {
      return(0)
    }else if(result<0){
      return(0)
    } else{
      return(result)
    }
  } else{
    stop("Don't recognize the function type.")
  }
}

## All Mosquito parameters 
## Drawn from: https://onlinelibrary.wiley.com/doi/10.1111/ele.12015
## Mosquito biting rate (a)
br <- list(fxn = 'briere',
           ## For Briere, parms are ordered as c, tm, t0
           parms = c(0.000203, 42.3, 11.7),
           parms_sd = c(0.0000576, 3.53, 2.47))

tibble(temperature = rep(0:50, 100),
       biting_rate = map(temperature, 
                         draw_param_val, 
                         param_list = br,
                         draw_sample = T) |> unlist()) |> 
  ggplot(aes(temperature, biting_rate)) +
    geom_point(alpha = .3) +
    stat_smooth(se=F) +
    background_grid(major = 'xy') +
    labs(x = expression('Temperature ('*degree*C*')'), y = 'Adult survival') ->br_plot
br_plot



## Daily adult survival probability (e^-mu) or e^-mortality_rate
dasp <- list(fxn = 'quadratic',
             ## For quadratic, parms are ordered as q, r, s
             parms = c(-0.000828, 0.0367, 0.522),
             parms_sd = c(0.0000519, 0.00239, 0.0235))

tibble(temperature = rep(0:40, 1),
       adult_survival = map(temperature, 
                            draw_param_val, 
                            param_list = dasp,
                            draw_sample = F) |> unlist()) |> 
  ggplot(aes(temperature, -log(adult_survival))) +
  geom_line() +
  background_grid(major = 'xy') +
  labs(x = expression('Temperature ('*degree*C*')'), y = 'Daily mortality rate') ->mortality_rate_plot
mortality_rate_plot


## parasite development rate (1/n) or (1/extrinsic incubation period)
pdr <- list(fxn = 'briere',
            ## For Briere, parms are ordered as c, tm, t0
            parms = c(0.000111, 34.4, 14.7),
            parms_sd = c(0.0000161, 0.000176, 1.48))

tibble(temperature = rep(10:35, 100),
       parasite_development = map(temperature, 
                                  draw_param_val, 
                                  param_list = pdr,
                                  draw_sample = T) |> unlist()) |> 
  ggplot(aes(temperature, parasite_development)) +
  geom_point(alpha = .3) +
  background_grid(major = 'xy') +
  stat_smooth(se=F) +
  labs(x = expression('Temperature ('*degree*C*')'), y = 'Parasite development rate') -> pdr_plot
pdr_plot


## Combine them into a single plot
plot_grid(scale_factor_plot, 
          br_plot, 
          mortality_rate_plot, 
          pdr_plot, 
          nrow = 2, 
          align = 'hv', 
          axis = 'tblr',
          labels = 'AUTO') -> county_varying_parm_plot
county_varying_parm_plot
save_plot('figs/county_varying_parm_plot.png', 
          county_varying_parm_plot,
          base_height = 6.5,
          base_asp = 1.3,
          bg = 'white')



# Create figure showing the monthly mean and upperbound R0 estimates --------------------
load('processed-data/2023-08-09_county-monthly-rnot-estimates.rda')
county_sf |>
  inner_join(county_month_rnot_summary , 
             by = 'fips', multiple = 'all') |> 
  filter(abbr != 'AK') |> 
  mutate(month = factor(month, levels= month.abb)) |> 
  ggplot() +
  geom_sf(aes(fill = rnot_mean)) +
  facet_wrap(~month, nrow = 4) +
  scale_fill_gradient(low = 'white', high = 'darkblue', na.value = 'lightgrey') +
  labs(fill = expression(Mean~italic(R[0]))) +
  theme_map() -> rnot_mean_month_map
rnot_mean_month_map
save_plot('figs/rnot_mean_month_map.png', 
          rnot_mean_month_map, 
          base_height = 10, 
          base_asp = 1.3, 
          bg = 'white')

county_sf |>
  inner_join(county_month_rnot_summary , 
             by = 'fips', multiple = 'all') |> 
  filter(abbr != 'AK') |> 
  mutate(month = factor(month, levels= month.abb)) |> 
  ggplot() +
  geom_sf(aes(fill = rnot_hi)) +
  facet_wrap(~month,nrow=4) +
  scale_fill_gradientn(colours = c('white', 'darkblue', 'yellow', 'red'), 
                       values = scales::rescale(c(0, 1, 1.01, max(county_month_rnot_summary$rnot_hi))),
                       na.value = 'lightgrey') +
  labs(fill = expression('97.5%'~italic(R[0]))) +
  theme_map() -> rnot_hi_month_map
rnot_hi_month_map
save_plot('figs/rnot_hi_month_map.png', 
          rnot_hi_month_map, 
          base_height = 10, 
          base_asp = 1.3, 
          bg = 'white')

# Create single rnot upperbound plot --------------------------------------
load('processed-data/2023-08-09_county-single-rnot-estimates.rda')

county_sf |>
  inner_join(county_single_rnot_summary ,
             by = 'fips') |>
  # inner_join(county_month_rnot_summary |>
  #               group_by(fips) |>
  #               mutate(rnot_mean = mean(head(sort(rnot_mean, decreasing = T), 6))), by = 'fips') |>
  filter(abbr != 'AK') |> 
  ggplot() +
  geom_sf(aes(fill = rnot_hi)) +
  # scale_fill_gradient(low = 'white', high = 'darkblue', na.value = 'lightgrey') +
  scale_fill_gradientn(colours = c('white', 'darkblue', 'yellow', 'red'),
                       values = scales::rescale(c(0, 1, 1.01, max(county_single_rnot_summary$rnot_hi))),
                       na.value = 'lightgrey') +
  labs(fill = expression('97.5%'~italic(R[0]))) +
  theme_map() -> single_rnot_hi_map
single_rnot_hi_map
save_plot('figs/single_rnot_hi_map.png', 
          single_rnot_hi_map, 
          base_height = 5, 
          base_asp = 1.8, 
          bg = 'white')

load('processed-data/county_gdp.rda')

county_sf |>
  inner_join(county_month_temperature, 
             by = 'fips') |> 
  filter(abbr != 'AK') |> 
  mutate(month = factor(month, levels= month.abb)) |> 
  ggplot() +
  geom_sf(aes(fill = avg_temp)) +
  scale_fill_gradient(low = 'white', high = 'darkblue', na.value = 'lightgrey') +
  theme_map() + 
  facet_wrap(~month, nrow = 4) -> temp_month_map
temp_month_map
save_plot('figs/temp_month_map.png', 
          temp_month_map, 
          base_height = 10, 
          base_asp = 1.3, 
          bg = 'white')

county_sf |>
  # distinct(abbr,full,county,fips) |> 
  inner_join(county_mosq_abundance, 
             by = 'fips') |> 
  # group_by(abbr) |> filter(mosq_prob == max(mosq_prob), mosq_prob !=0) -> temp
  filter(abbr != 'AK', abbr != 'HI') |> 
  # mutate(mosq_pres = mosq_prob>0.4) |> 
  ggplot() +
  geom_sf(aes(fill = -log(1-mosq_prob))) +
  scale_fill_gradient(low = 'white', high = 'darkblue', na.value = 'lightgrey') +
  # # scale_fill_gradientn(colours = c('white', 'darkblue', 'yellow', 'red'), 
  # #                      values = scales::rescale(c(0, 1, 1.01, max(county_single_rnot_summary$rnot_hi))),
  # #                      na.value = 'lightgrey') +
  # labs(fill = expression('97.5%'~italic(R[0]))) +
  theme_map() 
