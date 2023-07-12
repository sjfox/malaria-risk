## This file compiles the gdp data for each county in US
library(tidyverse)
library(tidycensus)

# Gather the per capita GDP data for counties -----------------------------
## only needed if need to search for correct variable
# temp <- load_variables(2020, "acs5", cache = F)
# temp$label[grepl('income', temp$label)]
# temp |> 
#   filter(concept == 'PER CAPITA INCOME IN THE PAST 12 MONTHS (IN 2020 INFLATION-ADJUSTED DOLLARS)') |> 
#   pull(name) -> gdp_percap_variable
gdp_percap_variable <- 'B19301_001'
## Get county-level GDP
county_gdp <- get_acs('county', variables = gdp_percap_variable)
county_gdp |> 
  select(fips = GEOID,
         gdp_percap = estimate) -> county_gdp

save(county_gdp, file = 'processed-data/county_gdp.rda')


# No longer needed in this file
# econ = seq(6.5,10.5,0.1)
# plot(
#   econ,exp(scam::predict.scam(scam.est.list[[1]],newdata=data.frame(econ=econ))),
#   type='l',col=rgb(0,0,0,.05),ylim = c(-3,2), xlab='',ylab='')
# for(ii in 2:length(scam.est.list)){
#   lines(econ,exp(scam::predict.scam(scam.est.list[[ii]],newdata=data.frame(econ=econ))),col=rgb(0,0,0,.1))
# }
# 
# 
# 
# 
# 
# econ = seq(6.5,10.5,0.1)
# plot(
#   econ,predict(scam.est.list[[1]],newdata=data.frame(econ=econ)),
#   type='l',col=rgb(0,0,0,.05),ylim=c(-1,6),xlab='',ylab='')
# for(ii in 2:length(scam.est.list)){
#   lines(econ,predict(scam.est.list[[ii]],newdata=data.frame(econ=econ)),col=rgb(0,0,0,.1))
# }
# mtext('ln(econ. index)',1,cex=.7,line=2)
# mtext('ln(m multiplier)',2,cex=.7,line=2)
# mtext('d',3,at=6.5,cex=.7)
