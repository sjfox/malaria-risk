##########################################################################################################
# Code to estimate the overdispersion parameter of negative binomial from given R0
#   and the proportion of P infections that cause Q secondary infection 
# e.g. 20% of the most infectious malaria infections are responsible for 80% of the secondary infections
# Based on this paper: https://onlinelibrary.wiley.com/doi/10.1111/tbed.14655#pane-pcw-figures
##########################################################################################################

library(tidyverse)

fit_dispersion <- function(dispersion=0.16, rnot = 1.1){
  P <- 0.2 # For malaria from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3128496/
  Q <- 0.8 # and https://www.nature.com/articles/438293a
  xs <- 0:100 # X = 100
  cum_sum_vals <- cumsum(xs * dnbinom(xs, mu = rnot, size = dispersion) / rnot) # 1 - Q
  x_val <- xs[ which(cum_sum_vals<(1-Q))[length( which(cum_sum_vals<(1-Q)) )] ] # get upper limit
  
  return( abs(1-P - pnbinom(q = x_val, mu = rnot, size = dispersion)) )
}

get_optimal_dispersion <- function(rnot){
  max_search = 10
  if(rnot==0){ # When 0 neg binom could return anything, so set as the max of interval below
    min_disp = max_search
  }else{
    fit_disp_param <- optimise(fit_dispersion, interval = c(0.0,max_search), rnot = rnot)
    min_disp = fit_disp_param$minimum
  } # end if
  return(min_disp)
} # end function

tibble(potential_rnots = seq(0.0, 3.65, 0.01) ) |> 
  mutate(dispersion = map(potential_rnots, get_optimal_dispersion) |> unlist()) -> combos

combos |> 
  # filter(dispersion < ) |> 
  filter(potential_rnots < 1) |> 
  ggplot(aes(potential_rnots, dispersion )) +
  geom_line()+
  theme_bw()


combos |> filter(potential_rnots==0.5) |> 
  summarize(mean(dispersion))
  mutate(is_trans = trans>0) |> 
  summarize(frac_trans = sum(is_trans)/n())
