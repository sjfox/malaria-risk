library(tidyverse)
library(cowplot)
library(sf)
theme_set(theme_cowplot())


load('processed-data/county_sf.rda')

rf_preds <- read_csv('processed-data/malaria_probability_next_import_comes_from_county.csv', 
                     col_types = 'ccn')
rf_preds <- rf_preds |> 
  select(-`...1`) |> 
  mutate(FIPS = str_pad(FIPS, 5, pad = "0")) |> 
  rename(fips = FIPS)



county_sf |> 
  left_join(rf_preds, by = 'fips') |> 
  filter(abbr != 'AK', abbr != 'HI') |> 
  ggplot() +
    geom_sf(aes(fill = log(`Next Import Probability`))) +
    theme_map() +
    scale_fill_gradient(low = 'white', high = 'darkblue') +
    labs(fill = 'Predicted\nImportat\nProbability') -> pred_imports
pred_imports
save_plot('figs/pred_imports.png', pred_imports,base_height = 6, base_asp = 1.6, bg = 'white')


imports_2017 <- read_csv('raw-data/Number_of_Reported_Malaria_Cases_by_County__United_States__2017.csv')
rf_preds |> 
  rename(import_prob = `Next Import Probability`) |> 
  filter(substr(fips,1,2) == '48') |> 
  arrange(desc(import_prob)) |> 
  left_join(imports_2017 |> 
              select(fips = FIPS,
                     malaria_cases = MAL_FREQ_2017), by = 'fips') |> 
  mutate(reported_cases = !is.na(malaria_cases)) |> 
  group_by(reported_cases) |> 
  summarize(total_cases = sum(malaria_cases),
            total_prob = sum(import_prob)) |> 
  mutate(total_prob = total_prob/sum(total_prob))


