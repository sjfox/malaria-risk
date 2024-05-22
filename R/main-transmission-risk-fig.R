## Analyze R0 estimates
library(tidyverse)
library(cowplot)
library(sf)
theme_set(theme_cowplot())


load('processed-data/county_sf.rda')

# Make the mosquito occurrence maps ----------------------------------------
load('processed-data/county_mosq_prob_clean.rda')

county_mosq_prob_clean |> 
  summarize_at(.vars = vars(-fips), .funs = sum) |> 
  gather(mosq_species, value) |> 
  filter(value!=0) |> 
  pull(mosq_species) -> species_present
  
county_mosq_prob_clean |> 
  select(fips, all_of(species_present)) |> 
  rename(`2010_Anopheles_crucians` = Anopheles_crucians) |> 
  gather(mosq_species, probability, -fips) |> 
  mutate(mosq_species = str_replace(mosq_species, '2010_Anopheles_', '')) |> 
  mutate(mosq_species = paste0(toupper(substring(mosq_species, 1,1)), substring(mosq_species, 2), sep="")) |> 
  group_by(fips) |> 
  filter(probability == max(probability), probability>0) -> dominant_mosq_county

get_mosq_abundance <- function(x){
  -log(1-x)
}
county_mosq_prob_clean |> 
  mutate(mosq_prob = 1-(apply((1-across(!starts_with("fips"))), 1, prod)),
         mosq_species_present = rowSums(ifelse(across(!starts_with("fips")) > 0.4, 1, 0))) |> 
  mutate(mosq_prob = ifelse(mosq_prob >.99, .99, mosq_prob)) |> 
  mutate(mosq_abundance_lambda = get_mosq_abundance(mosq_prob)) |>
  select(fips, mosq_prob, mosq_abundance_lambda, mosq_species_present) -> county_mosq_abundance

krig_legend <- function(minQ, maxQ, kvalues){
  # browser()
  vals <- seq(minQ, maxQ, by = 1.00)
  kvals <- kvalues[!is.na(kvalues)]
  dat <- as.data.frame(tidyr::expand_grid(vals, kvals))
  dat$kvals <- as.character(dat$kvals)
  dat$vals <- as.character(dat$vals)
  
  ggplot(dat, aes(x = kvals, y = vals, fill = kvals, alpha = vals, group = kvals)) +
    geom_raster() +
    scale_y_discrete(expand = c(0,0), name = "Abundance", breaks = vals) +
    scale_x_discrete(expand = c(0,0)) +
    coord_fixed(ratio = 1) +
    labs(x = NULL, y = '') +
    scale_alpha_discrete(range = c(0, 1)) +
    scale_fill_brewer(type = 'qual', palette = 2, na.translate = F) +
    theme_cowplot() %+replace% theme(legend.position = "none",
                                     axis.ticks.x = element_blank(),
                                     axis.line = element_blank(),
                                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
}

mosq_plot_df <- county_sf |> 
                  left_join(county_mosq_abundance, by = 'fips') |> 
                  left_join(dominant_mosq_county, by = 'fips') |> 
                  filter(abbr != 'AK', abbr != 'HI')
mosq_plot_df |> 
  ggplot() +
    geom_sf(aes(fill = mosq_species, alpha = mosq_abundance_lambda)) +
    scale_fill_brewer(type = 'qual', palette = 2,na.translate = F) +
    scale_alpha_continuous(range = c(0, 1)) +
    theme_map() %+replace% theme(legend.position = 'none') -> mosq_abundance_map
krig_legend(minQ = 0,
            maxQ = max(mosq_plot_df$mosq_abundance_lambda),
            kvalues = unique(mosq_plot_df$mosq_species)) -> mosq_abundance_map_legend

mosq_abundance_map_combined <- plot_grid(mosq_abundance_map, 
                                         mosq_abundance_map_legend, rel_widths = c(0.8,0.2))
mosq_abundance_map_combined

save_plot('figs/mosq_abundance_map_combined2.png', 
          mosq_abundance_map_combined,
          base_height = 4.5,
          base_asp = 1.8,
          bg='white')


# Make R0 figure ----------------------------------------------------------
load('processed-data/2023-08-09_county-monthly-rnot-estimates.rda')
load('processed-data/2023-08-09_county-single-rnot-estimates.rda')
load('processed-data/malaria-risk-mcmc-samples.rda')

set.seed(808)
alpha_samps |> 
  sample_n(1000, replace = F) |> 
  mutate(samp = seq_along(alpha))-> alpha_reduced_samp

county_parm_samples |> 
  left_join(alpha_reduced_samp, by = 'samp') |> 
  group_by(fips, month) |> 
  summarize(rnot_mean = mean(rnot*alpha),
            rnot_lo = quantile(rnot*alpha,probs = 0.025),
            rnot_hi = quantile(rnot*alpha, probs = 0.975)) -> post_county_rnot_summary

## Pull counties that have at least one month where 97.5% quantile is above 1
post_county_rnot_summary |> 
  filter(rnot_hi>1) |> 
  distinct(fips) |> 
  mutate(rnot_above_one = T) -> county_risk_above1


county_sf |> 
  left_join(post_county_rnot_summary |> 
              filter(month == 'Aug'), by = 'fips') |> 
  filter(abbr != 'AK', abbr != 'HI') |> 
  ggplot() +
    geom_sf(aes(fill=rnot_mean)) +
    geom_sf(data = county_sf |>
                    inner_join(county_risk_above1, by = 'fips'),
            color = 'black', fill = NA, linewidth = .3) +
    scale_fill_gradient(low = 'white', high = 'darkred', na.value = 'white') +
    labs(fill = expression(Mean~italic(R[0]))) +
    guides(color="none", linewidth = 'none') +
    theme_map() -> rnot_mean_map
rnot_mean_map

save_plot('figs/rnot_mean_map2.png', 
          rnot_mean_map,
          base_height = 4.5,
          base_asp = 1.8,
          bg='white')


plot_grid(mosq_abundance_map, 
          # mosq_abundance_map_legend,
          rnot_mean_map + theme(legend.position = 'none'), 
          labels = 'AUTO',
          # get_legend(rnot_mean_map),
          # align = 'hv', axis = 'tblr',
          nrow = 2) -> transmission_risk_figure

plot_grid(transmission_risk_figure, plot_grid(mosq_abundance_map_legend,
          get_legend(rnot_mean_map + 
                       theme(legend.box.margin = margin(0, 0, 0, 50))), nrow = 2), 
          nrow = 1, rel_widths = c(1, .25)) -> temp

  
save_plot(filename = 'figs/transmission-risk-figure.png', 
            plot = temp,
            base_height = 7, 
            base_asp = 1.1,
          bg = 'white')

alpha_samps |> 
  ggplot(aes(alpha)) + 
  geom_density() +
  labs(x = 'Alpha', y = 'Density')

