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
    mutate(anoph_prob = exact_extract(an_tiff, geometry, fun = 'mean')) |> 
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



## Map drawn from the mosquitoes of the southeast book
## Focused only on regions of US where crucians complex is underrepresented from other anopheles species (mainly in the range of the bradleyi species)
## Selected from: https://www.mapchart.net/usa-counties.html removed all instances of spaces in names (__) and ___
crucians_ctys <- tibble(county = c("Terrebonne__LA", "Lafourche__LA", "Plaquemines__LA", "Kleberg__TX", "Willacy__TX", 
                                   "Cameron__TX", "Hidalgo__TX", "Jim_Hogg__TX", "Brooks__TX", "Starr__TX", "Zapata__TX", 
                                   "Duval__TX", "Jim_Wells__TX", "Nueces__TX", "Kenedy__TX", "San_Patricio__TX", "Bee__TX", 
                                   "Live_Oak__TX", "Goliad__TX", "Refugio__TX", "Aransas__TX", "Calhoun__TX", "Matagorda__TX", 
                                   "Jackson__TX", "Wharton__TX", "Victoria__TX", "Webb__TX", "La_Salle__TX", "Dimmit__TX", 
                                   "McMullen__TX", "Karnes__TX", "DeWitt__TX", "Lavaca__TX", "Colorado__TX", "Austin__TX", 
                                   "Waller__TX", "Harris__TX", "Fort_Bend__TX", "Brazoria__TX", "Galveston__TX", "Chambers__TX", 
                                   "Jefferson__TX", "Liberty__TX", "San_Jacinto__TX", "Polk__TX", "Tyler__TX", "Hardin__TX", 
                                   "Jasper__TX", "Newton__TX", "Orange__TX", "Montgomery__TX", "Cameron__LA", "Calcasieu__LA", 
                                   "Acadia__LA", "Jefferson_Davis__LA", "Vermilion__LA", "St_Martin__LA", "Lafayette__LA", 
                                   "Iberia__LA", "St_Mary__LA", "Assumption__LA", "St_James__LA", "St_John_the_Baptist__LA", 
                                   "St_Charles__LA", "Jefferson__LA", "Orleans__LA", "St_Bernard__LA", "Livingston__LA", 
                                   "Ascension__LA", "Tangipahoa__LA", "St_Tammany__LA", "Washington__LA", "Hancock__MS", 
                                   "Pearl_River__MS", "Harrison__MS", "Stone__MS", "George__MS", "Jackson__MS", "Mobile__AL", 
                                   "Baldwin__AL", "Santa_Rosa__FL", "Escambia__FL", "Okaloosa__FL", "Walton__FL", "Iberville__LA", 
                                   "Beauregard__LA", "Allen__LA", "Evangeline__LA", "St_Landry__LA", "West_Baton_Rouge__LA", 
                                   "East_Baton_Rouge__LA", "Holmes__FL", "Washington__FL", "Bay__FL", "Calhoun__FL", "Gulf__FL", 
                                   "Jackson__FL", "Liberty__FL", "Franklin__FL", "Gadsden__FL", "Leon__FL", "Wakulla__FL", 
                                   "Madison__FL", "Jefferson__FL", "Decatur__GA", "Thomas__GA", "Grady__GA", "Brooks__GA", 
                                   "Lowndes__GA", "Echols__GA", "Hamilton__FL", "Taylor__FL", "Lafayette__FL", "Dixie__FL", 
                                   "Monroe__FL", "Miami_Dade__FL", "Collier__FL", "Broward__FL", "Palm_Beach__FL", "Hendry__FL", 
                                   "Lee__FL", "Martin__FL", "St_Lucie__FL", "Indian_River__FL", "Brevard__FL", "Volusia__FL", 
                                   "Flagler__FL", "St_Johns__FL", "Duval__FL", "Nassau__FL", "Levy__FL", "Citrus__FL", "Hernando__FL", 
                                   "Pasco__FL", "Pinellas__FL", "Hillsborough__FL", "Manatee__FL", "Sarasota__FL", "Charlotte__FL", 
                                   "DeSoto__FL", "Hardee__FL", "Highlands__FL", "Glades__FL", "Okeechobee__FL", "Osceola__FL", 
                                   "Polk__FL", "Lake__FL", "Orange__FL", "Seminole__FL", "Sumter__FL", "Marion__FL", "Putnam__FL", 
                                   "Alachua__FL", "Gilchrist__FL", "Suwannee__FL", "Columbia__FL", "Union__FL", "Bradford__FL", 
                                   "Clay__FL", "Baker__FL", "Camden__GA", "Glynn__GA", "McIntosh__GA", "Liberty__GA", "Bryan__GA", 
                                   "Long__GA", "Wayne__GA", "Brantley__GA", "Charlton__GA", "Chatham__GA", "Effingham__GA", 
                                   "Jasper__SC", "Beaufort__SC", "Hampton__SC", "Colleton__SC", "Charleston__SC", "Dorchester__SC", 
                                   "Berkeley__SC", "Williamsburg__SC", "Georgetown__SC", "Horry__SC", "Marion__SC", "Florence__SC", 
                                   "Columbus__NC", "Brunswick__NC", "New_Hanover__NC", "Pender__NC", "Onslow__NC", "Duplin__NC", 
                                   "Jones__NC", "Craven__NC", "Carteret__NC", "Pamlico__NC", "Beaufort__NC", "Hyde__NC", "Dare__NC", 
                                   "Tyrrell__NC", "Washington__NC", "Chowan__NC", "Perquimans__NC", "Pasquotank__NC", "Camden__NC", 
                                   "Currituck__NC", "Gates__NC", "Bertie__NC", "Martin__NC", "Pitt__NC", "Lenoir__NC", "Bladen__NC", 
                                   "Sampson__NC", "Virginia_Beach__VA", "Norfolk__VA", "Portsmouth__VA", "Chesapeake__VA", 
                                   "Suffolk__VA", "Isle_of_Wight__VA", "Surry__VA", "James_City__VA", "Newport_News__VA", 
                                   "Hampton__VA", "Poquoson__VA", "York__VA", "Gloucester__VA", "Mathews__VA", "Middlesex__VA", 
                                   "Lancaster__VA", "Northumberland__VA", "Richmond_Co___VA", "Westmoreland__VA", "Essex__VA", 
                                   "King_and_Queen__VA", "New_Kent__VA", "King_William__VA", "King_George__VA", "Stafford__VA", 
                                   "Prince_William__VA", "Fairfax_Co___VA", "Charles__MD", "Prince_George_s__MD", "St_Mary_s__MD", 
                                   "Calvert__MD", "Northampton__VA", "Accomack__VA", "Somerset__MD", "Dorchester__MD", "Wicomico__MD", 
                                   "Worcester__MD", "Talbot__MD", "Queen_Anne_s__MD", "Kent__MD", "Cecil__MD", "Harford__MD", 
                                   "Baltimore_County__MD", "Anne_Arundel__MD", "Baltimore_City__MD", "Sussex__DE", "Kent__DE", 
                                   "Caroline__MD", "New_Castle__DE", "Cumberland__NJ", "Cape_May__NJ", "Salem__NJ", "Gloucester__NJ", 
                                   "Atlantic__NJ", "Ocean__NJ", "Monmouth__NJ", "Middlesex__NJ", "Burlington__NJ", "Camden__NJ", 
                                   "Maverick__TX", "Zavala__TX", "Frio__TX", "Atascosa__TX", "Ware__GA", "Clinch__GA", "Pierce__GA", 
                                   "Appling__GA", "Tattnall__GA", "Evans__GA", "Bulloch__GA", "Screven__GA", "Allendale__SC", 
                                   "Bamberg__SC", "Orangeburg__SC", "Clarendon__SC", "Sumter__SC", "Lee__SC", "Darlington__SC", 
                                   "Dillon__SC", "Robeson__NC"))

crucians_ctys |> 
  separate(col = 'county', into =c('county', 'state'), sep = '__') |> 
  mutate(county=str_replace(county, '_', ' ')) |> 
  mutate(state=str_replace(state, '_', '')) |> 
  mutate(county=str_replace(county, '_', ' ')) |> 
  mutate(county=str_replace(county, '_', ' ')) |> 
  mutate(county=str_replace(county, ' Co', '')) |> 
  mutate(county=str_replace(county, ' s', "'s")) |> 
  mutate(county = ifelse(county == 'Miami Dade', 'Miami-Dade', county)) |>
  mutate(county = ifelse(county == 'Baltimoreunty', 'Baltimore', county))  |> 
  inner_join(county_sf |> 
              distinct(fips, county,abbr) |> 
              mutate(county = str_replace(county, ' County', '')) |>
              filter(county!= 'Richmond city', county != 'Fairfax city') |> 
              mutate(county = ifelse(county == 'Baltimore city', 'Baltimore City', county)) |> 
              mutate(county = str_replace(county, ' city', '')) |>
              mutate(county = str_replace(county, ' Parish', '')) |> 
              mutate(county = str_replace(county, '\\.', '')), 
             by = c('state' = 'abbr', 'county')) -> crucians_ctys_clean

## Should be same number of rows as crucians
nrow(crucians_ctys) == nrow(crucians_ctys_clean)

crucians_addition <- crucians_ctys_clean |> 
  select(fips) |> 
  mutate(Anopheles_crucians = 0.95)



# Combine into single dataframe -------------------------------------------
county_mosq_prob_clean |> 
  left_join(crucians_addition, by = 'fips') |> 
  mutate_at(.vars = vars(-fips), .funs = ~replace_na(.,0)) -> county_mosq_prob_clean


# county_mosq_prob_clean |>
#   filter(fips == '12115') |>
#   gather(key, value, -fips) |>
#   arrange(desc(value))
# county_mosq_prob_clean |>
#   filter(fips == '12049') |>
#   gather(key, value, -fips) |>
#   arrange(desc(value))
# 1 12115 2010_Anopheles_quadrimacul  0.434

## Get county shaprefiles for future plotting
county_sf <- sf_polygon(obj = usmap::us_map('counties'),
           x = 'x', y = 'y', polygon_id = "group", keep = T) |> 
  select(fips, abbr,full,county)
  

save(county_mosq_prob_clean,inat_pts, file = 'processed-data/county_mosq_prob_clean.rda')
save(county_sf, file = 'processed-data/county_sf.rda')



# Add code to save the non-crucians one -----------------------------------





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


# ## Estimate the missing anopheles malaria vectors from inaturalist
# inat <- read_csv('raw-data/observations-350305.csv')
# inat_pts <- st_as_sf(inat, coords = c('longitude','latitude'), crs = crs(us_counties))
# inat_pts |>
#   mutate(fips = us_counties$fips[as.integer(st_intersects(inat_pts, us_counties))]) |>
#   select(fips, scientific_name) -> inat_pts
# 
# anoph_additions <- c('Anopheles punctipennis',
#                      'Anopheles crucians',
#                      # 'Anopheles barberi',
#                      'Anopheles bradleyi')
# 
# 
# inat_pts |>
#   filter(scientific_name %in% anoph_additions) |>
#   distinct(fips, scientific_name) |>
#   mutate(scientific_name = str_replace(scientific_name, ' ', replacement = '_'),
#          prob = 0.95) |>
#   spread(scientific_name, prob) |>
#   mutate_at(.vars = vars(-fips), .funs = ~replace_na(.,0)) -> county_anoph_additions
# 
# county_sf |> 
#   left_join(county_anoph_additions, by = 'fips') |> 
#   filter(abbr != 'AK', abbr != 'HI') |> 
#   ggplot() +
#   geom_sf(aes(fill=Anopheles_punctipennis)) +
#   scale_fill_gradient(low = 'white', high = 'darkred', na.value = 'white')

