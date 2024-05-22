# Malaria transmission risk estimation

## Description of code
- To run the code first you have to run `convert-mosq-raster-to-county.R` which takes all of the mosquito occurrence probabilities and averages them at the county level. It then saves the county mosquito probabilities into a nice usable format.
- Run `get-county-gdp.R` to gather all of the county per capita GDP estimates from ACD
- Run `calc-malaria-county-rnots.R` to calculate the R0s based on all of the data and parameter estimates gathered.
- Run `posterior-malaria-rnots.R` to calculate the posterior distributions for R0 based on the prior distribution and the importations estimated.
- `calc-county-temperatures.R` is not needed if using Mathew's provided temperature data


## Data needed for running
- You need the `anopholes_maps/` folder in raw-data, which has the compressed zip of mosquito occurrence provided by the malaria atlas project. You can unzip them yourself or the code can do it for you.
- Right now you don't need the `uscounties_malaria_case.csv` file in raw-data, but eventually that will probably be used.

