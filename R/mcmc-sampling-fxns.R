############################################
## File holding all code necessary for running the mcmc
############################################

require(Rcpp)

sourceCpp("cpp/cpp_fitting_fxns.cpp")

# draw_zika_rnots <- function(gamma_parms){
#   num_rnots <- nrow(rnot_dist)
#   rnot_sample_inds <- sample(x = seq(1, ncol(rnot_dist)), size = num_rnots, replace = T)
#   as.numeric(rnot_dist[1:num_rnots, rnot_sample_inds])
# }


lprior <- function(alpha, parms){
  if(parms$inform_prior){
    sum(dgamma(parms$rnot, shape = parms$rnot_dist$shape, rate = parms$rnot_dist$rate, log = T))
  } else{
    0
  }
}

llprior <- function(alpha, parms2){
  # print('here')
  loglik <- try(scaling_loglike_cpp(alpha = alpha, params = parms2))
  if(class(loglik) == "try-error"){browser()}
  prior <- lprior(alpha, parms2)
  if(is.na(loglik + prior) | is.infinite(loglik + prior)) browser()
  # browser()
  # print(paste0("loglik: ", loglik, " prior: ", prior))
  return(loglik + prior)
}


draw_new_alpha <- function(alpha, tuning){
  runif(1, min = 0, max = 2)
  # plogis(qlogis(alpha) + rnorm(1, sd = tuning))
}
draw_new_reporting_rate <- function(rr, tuning){
  # plogis(qlogis(rr) + rnorm(1, sd = tuning))
  ## Asymptomatic rate:
  ## 70.9% range from 58.5% to 92.6%
  ## https://www.thelancet.com/journals/lanam/article/PIIS2667-193X(21)00165-4/fulltext#seccesectitle0016
  ## Set standard deviation so the max and min observed are 3 standard deviations from the mean
  mal_sd <- mean(c(abs(.926 - .709)/3, abs(.585 - .709)/3))
  
  ## reporting rate 
  reporting <- 1 - rnorm(1, mean = .709, sd = mal_sd)
  while(reporting<0){
    reporting <- 1 - rnorm(1, mean = .709, sd = mal_sd)
  }
  reporting
}

draw_new_rnots <- function(rnots, tuning){
  new_rnots <- exp(log(rnots) + rnorm(n = length(rnots), sd = tuning))
  ifelse(new_rnots<1e-16, 1e-16, new_rnots)
}


## MCMC pseudocode
mcmc_zika_rnot <- function (zika_parms,
                            alpha_tuning,
                            rnot_tuning,
                            rr_tuning,
                            burnin = 1000,
                            N= 10000,
                            thin = 1){
  # browser()
  accept <- 0
  
  ## rnot_cols stores # of columns necessary for storing posterior rnots
  ## Only estimating posterior proabbilities for necessary counties
  ## All other county posteriors can be retrospectively estimated
  ## straight from their prior distributions
  num_rnots <- nrow(zika_parms$rnot_dist)
  
  ###### Create matrix for saving
  ## adding 2 extra columns into estimate, alpha posterior and loglike
  saved_samps <- matrix(data = 0, nrow = (N-burnin)/thin, ncol = num_rnots+3)
  
  ###### Draw Rnots and alpha
  curr_rnots <- runif(num_rnots, min = 0, max = 2)
  ## Assumption that alpha is between 0 and 1, could change in future iterations
  curr_alpha <- runif(1)
  
  ## test sampling reporting rate
  curr_rr <- runif(1)
  
  ## Make sure duplicated county/month Rnots are the same
  ## Make the last instances equal the first ones
  if(anyDuplicated(zika_parms$county_month)){
    curr_rnots[duplicated(zika_parms$county_month)] <- curr_rnots[duplicated(zika_parms$county_month, fromLast = TRUE)]
  }
  ###### Calc log like + log prior
  curr_parms <- subs_parms(list(rnot = curr_rnots, reporting_rate = curr_rr), zika_parms)
  # browser()
  curr_llprior <- llprior(curr_alpha, curr_parms)
  
  ##### No longer save first results
  # saved_samps[1, ] <- c(curr_llprior, curr_alpha, curr_rnots)
  
  for( ii in 2:N){
    ###### Draw Proposed Rnots and alpha
    proposed_rnots <- draw_new_rnots(curr_rnots, rnot_tuning)
    ## link the duplicated months here as well
    if(anyDuplicated(zika_parms$county_month)){
      proposed_rnots[duplicated(zika_parms$county_month)] <- proposed_rnots[duplicated(zika_parms$county_month, fromLast = TRUE)]
    }
    
    
    ## Assumption that alpha is between 0 and 1, could change in future iterations
    proposed_alpha <- draw_new_alpha(curr_alpha, alpha_tuning)
    
    proposed_rr <- draw_new_reporting_rate(curr_rr, rr_tuning)
    
    ###### Calc log like
    proposed_parms <- subs_parms(list(rnot = proposed_rnots, reporting_rate = proposed_rr), zika_parms)
    proposed_llprior <- llprior(proposed_alpha, proposed_parms)
    # print(ii)
    mh_prob <- proposed_llprior - curr_llprior
    if(is.na(mh_prob) | is.infinite(mh_prob)) browser()
    if(mh_prob >= 0 | (runif(1) <= exp(mh_prob))){
      curr_alpha <- proposed_alpha
      curr_rnots <- proposed_rnots
      curr_llprior <- proposed_llprior
      accept <- accept + 1
    }
    curr_rr <- proposed_rr
    
    if(N > burnin & ii %% thin == 0){
      saved_samps[(ii - burnin)/thin, ] <- c(curr_llprior, curr_alpha, curr_rr, curr_alpha*curr_rnots)
    }
    
  }
  return(list(samples = saved_samps,
              aratio = accept/N))
}
