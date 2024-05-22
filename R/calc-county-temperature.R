#################################################
## All data can be downloaded from the repository specified in the README.md
#################################################


####################################################
## Processes the historic temperature data
## No longer used

library(tidyverse)
library(weathermetrics)


# If we were to use, the fips codes are not correct -----------------------
## Readme here gives information: https://www.ncei.noaa.gov/pub/data/cirs/climdiv/

# read_tsv('https://www.ncei.noaa.gov/pub/data/cirs/climdiv/climdiv-tmpccy-v1.0.0-20230707',
#          col_names = F) -> df
#   
# df |> 
#   mutate(fips = substr(X1, 1, 5),
#          element = substr(X1, 6, 7),
#          year = substr(X1, 8, 11),
#          Jan = substr(X1, 12, 18),
#          Feb = substr(X1, 19, 25),
#          Mar = substr(X1, 26, 32),
#          Apr = substr(X1, 33, 39),
#          May = substr(X1, 40, 46),
#          Jun = substr(X1, 47, 53),
#          Jul = substr(X1, 54, 60),
#          Aug = substr(X1, 61, 67),
#          Sep = substr(X1, 68, 74),
#          Oct = substr(X1, 75, 81),
#          Nov = substr(X1, 82, 88),
#          Dec = substr(X1, 89, 95)
#          ) |> 
#   select(-X1) |> 
#   mutate_at(vars(-fips,-element, -year), as.numeric) |> 
#   mutate_at(vars(-fips,-element, -year), fahrenheit.to.celsius) |> 
#   filter(element == '02') |> 
#   select(-element) -> temperature_df
# 
# temperature_df <- read_csv('processed-data/fixed_tmp_data.csv')
# 
# temperature_df |> 
#   rename(fips=FIPS, year = Year) |> 
#   filter(year %in% c(2013:2022)) |>
#   select(-year) |>
#   group_by(fips) |>
#   summarize_all(.funs=mean) |> 
#   gather(month, avg_temp, -fips) |> 
#   filter(month == 'May') |> 
#   left_join(temperature_df |> 
#               rename(fips=FIPS, year = Year) |> 
#                filter(year %in% c(2023)) |>
#                select(-year) |>
#                group_by(fips) |>
#                summarize_all(.funs=mean) |> 
#                gather(month, avg_temp_2023, -fips) |> 
#                filter(month == 'May'), by = 'fips') |> 
#   select(fips, avg_temp, avg_temp_2023) |> 
#   mutate(relative_diff = avg_temp_2023/avg_temp) -> year_comparison
# 
# 
# load('processed-data/county_sf.rda')
# 
# county_sf |> 
#   left_join(year_comparison, by = 'fips') |> 
#   ggplot() +
#   geom_sf(aes(fill = relative_diff)) +
#   scale_fill_gradientn(colours = c('white', 'darkblue', 'yellow', 'red'), 
#                        values = scales::rescale(c(min(year_comparison$relative_diff), 1, 1.01, max(year_comparison$relative_diff))),
#                        na.value = 'lightgrey') +
#   theme_map() 


# library(tidyverse)
# library(terra)
# library(exactextractr)
# library(sf)
# 
# compile_temp_df <- function(raster_location, us_counties){
#   ## Takes in the raster location and compiles all of the average county temperatures from that raster
#   temperature_raster <- rast(raster_location) # "data/tmean_2-5m_bil/tmean1.bil")
#   
#   us_counties |> 
#     st_transform(st_crs(temperature_raster)) |> 
#     mutate(avg_temp = exact_extract(temperature_raster, geometry, 'mean')/10) |> 
#     mutate(month = month.abb[as.numeric(strsplit(strsplit(raster_location, split = ".", fixed = T)[[1]], split="n", fixed=T)[[1]][3])]) |> 
#     as_tibble() |> 
#     dplyr::select(-geometry)
# }
# 
# us_counties <- tigris::counties() |> 
#   dplyr::select(fips = GEOID) 
# 
# raster_locations <- list.files(path="raw-data/tmean_2-5m_bil", pattern = "*.bil", full.names = T)
# county_temperatures <- purrr::map(raster_locations, compile_temp_df, us_counties = us_counties)
# 
# 
# 
# county_temperatures %>% 
#   bind_rows() %>% 
#   mutate(month = factor(month, levels = month.abb)) %>%
#   arrange(month, fips) |> 
#   write_csv(file= "processed-data/county-historic-temps.csv")
