## Code to estimate an overdispersion parameter for R0
## Based on 20/80 rule from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3128496/
## and https://www.nature.com/articles/438293a

## Converts that to overdispersion parameter for a specific R0 using the method from:
## https://onlinelibrary.wiley.com/doi/10.1111/tbed.14655

library(tidyverse)

fit_dispersion <- function(dispersion, rnot = 1.1){
  if(rnot == 0){
    return(100)
  }
  P <- 0.2
  Q <- 0.8
  xs <- 0:100
  cum_sum_vals <- cumsum(xs * dnbinom(xs, mu = rnot, size = dispersion) / rnot)
  x_val <- xs[which(cum_sum_vals < 1- Q)[length(which(cum_sum_vals < 1- Q))]]
  abs(1-P - pnbinom(q = x_val, mu = rnot, size = dispersion))
}

get_optimal_dispersion <- function(rnot){
  fit_disp_param <- optimise(fit_dispersion, interval = c(0.0,10), rnot = rnot)
  fit_disp_param$minimum
}

tibble(potential_rnots = seq(0.00, 3, length = 1000)) |> 
  mutate(dispersion = map(potential_rnots, get_optimal_dispersion) |> unlist()) -> combos
combos

combos |> 
  # filter(dispersion < ) |> 
  ggplot(aes(potential_rnots, dispersion)) +
    geom_line()


