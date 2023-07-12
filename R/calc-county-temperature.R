#################################################
## All data can be downloaded from the repository specified in the README.md
#################################################


####################################################
## Processes the historic temperature data
library(tidyverse)
library(terra)
library(exactextractr)
library(sf)

compile_temp_df <- function(raster_location, us_counties){
  ## Takes in the raster location and compiles all of the average county temperatures from that raster
  temperature_raster <- rast(raster_location) # "data/tmean_2-5m_bil/tmean1.bil")
  
  us_counties |> 
    st_transform(st_crs(temperature_raster)) |> 
    mutate(avg_temp = exact_extract(temperature_raster, geometry, 'mean')/10) |> 
    mutate(month = month.abb[as.numeric(strsplit(strsplit(raster_location, split = ".", fixed = T)[[1]], split="n", fixed=T)[[1]][3])]) |> 
    as_tibble() |> 
    dplyr::select(-geometry)
}

us_counties <- tigris::counties() |> 
  dplyr::select(fips = GEOID) 

raster_locations <- list.files(path="raw-data/tmean_2-5m_bil", pattern = "*.bil", full.names = T)
county_temperatures <- purrr::map(raster_locations, compile_temp_df, us_counties = us_counties)



county_temperatures %>% 
  bind_rows() %>% 
  mutate(month = factor(month, levels = month.abb)) %>%
  arrange(month, fips) |> 
  write_csv(file= "processed-data/county-historic-temps.csv")
