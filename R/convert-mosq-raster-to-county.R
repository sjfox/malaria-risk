library(terra)
library(tidyverse)
library(tigris)
library(sf)
library(exactextractr)
library(sfheaders)



## Get the county shapefile for clipping mosquito rasters
tigris_counties_loc <- 'processed-data/us-counties-tigris.rda'
if(file.exists(tigris_counties_loc)){
  load(tigris_counties_loc)
} else{
  us_counties <- tigris::counties() |> 
    select(fips = GEOID) 
  save(us_counties, file = tigris_counties_loc)
}



clip_anoph_data <- function(an_map_loc, us_counties){
  ## Unzips the folder downloaded from https://data.malariaatlas.org/
  ## Only downloaded 2010 ones, becase 2016 and 2017 don't impact united states
  ## Clips occurrence probability to US counties
  ## Returns county fips and probability of mosquito occurrence data frame
  ## Maps from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3349467/
  
  
  # browser()
  new_folder_loc <- str_replace(string = an_map_loc, pattern = '.zip', replace = '')
  anoph_species <- str_replace(new_folder_loc, pattern = 'raw-data/anopholes_maps/', replace = '')
  if(!dir.exists(new_folder_loc)){
    unzip(an_map_loc, exdir = new_folder_loc)  
  }
  tif_file <- list.files(new_folder_loc, pattern = '.tif', full.names = T)
  if(length(tif_file) > 1){
    browser()
  } else{
    an_tiff <- rast(tif_file)
  }
  
  us_counties |> 
    st_transform(crs(an_tiff)) |> 
    mutate(anoph_prob = exact_extract(an_tiff, geometry, 'mean', default_value = 0)) |> 
    as_tibble() |> 
    mutate(anoph_prob = ifelse(is.nan(anoph_prob), 0 , anoph_prob)) |> 
    select(-geometry) -> county_prob
  names(county_prob)[2] <- anoph_species
  county_prob
}

## All zipped anopholes maps
an_maps <- list.files('raw-data/anopholes_maps', pattern = '.zip', full.names = T)
mosq_probs <- map(an_maps, clip_anoph_data, us_counties = us_counties)

## Combine into single mosquito probability tibble
mosq_probs |> 
  reduce(.f = left_join, by = 'fips') -> county_mosq_prob_clean


## Get county shaprefiles for future plotting
county_sf <- sf_polygon(obj = usmap::us_map('counties'),
           x = 'x', y = 'y', polygon_id = "group", keep = T) |> 
  select(fips, abbr,full,county)
  

save(county_mosq_prob_clean, file = 'processed-data/county_mosq_prob_clean.rda')
save(county_sf, file = 'processed-data/county_sf.rda')


# Plotting code not needed - not sure if it still works, but gives idea of how to plot ----
# county_sf |> 
#   group_by(fips) |> 
#   summarize(geometry = st_union(geometry)) |> 
#   left_join(county_mosq_prob |> 
#               as_tibble() |> 
#               select(-geometry), by = 'fips') |> 
#   mutate_at(.vars = vars(freeborni, pseudopunctipennis,quadrimaculatus),
#             .funs = ~ifelse(is.nan(.), 0, .)) -> county_mosq_prob_clean
# 
# county_mosq_prob_clean |> 
#   gather(species,probability, freeborni, pseudopunctipennis, quadrimaculatus) |> 
#   ggplot() +
#   geom_sf(aes(fill = probability), color="black", size=0.25) +
#   facet_wrap(~species) +
#   scale_fill_gradient(low = 'white', high = 'darkblue')
# 
# county_mosq_prob_clean |> 
#   gather(species,probability, freeborni, pseudopunctipennis, quadrimaculatus) |> 
#   mutate(abundance = -log(1-probability)) |> 
#   group_by(fips) |> 
#   summarize(abundance = sum(abundance,na.rm=T)) |> 
#   ggplot() +
#   geom_sf(aes(fill = abundance), color="black", size=0.25) +
#   scale_fill_gradient(low = 'white', high = 'darkblue')

