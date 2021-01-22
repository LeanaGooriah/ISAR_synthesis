# visual inspection of model fits
# ISAR data

library(brms)
library(tidyverse)

# load the model fits
load('./results/univariate_models.Rdata')

# study-levels
study_levels <- Sn_lnorm$data %>%
  as_tibble() %>%
  distinct(Study.x) %>%
  mutate(level = Study.x) %>%
  nest(data = c(level))

# meta data
env_file0 <- read_csv("./data/env_file_august.xlsx - new.csv")

DEM_filter <- tibble(study_ID = c('Andrade_2002', 'Bell_2017', 'Dasilva_2019', 'Jonsson_2011', 'Macdonald_2018', 'Pereira_2017', 'Perillo_2017', 'Werden_2012',
                                  'Werner_Zalewski_2006'))

elev <- read.delim("./data/env_file_a.csv", sep = ',') %>%
  as_tibble() %>%
  select(study_ID, Study, island_code, elevation_m_asl, DEM_elevation_max) %>%
  rename(Study.x = Study) %>%
  # filter to studies in the model
  filter(Study.x %in% study_levels$Study.x)

elev_als <- elev %>%
  filter(study_ID %in% DEM_filter$study_ID) %>%
  mutate(elevation = as.numeric(elevation_m_asl))

elev_dem <- elev %>%
  filter(!study_ID %in% DEM_filter$study_ID) %>%
  mutate(elevation = as.numeric(DEM_elevation_max))

elev_range <- bind_rows(elev_als,
                        elev_dem) %>%
  group_by(Study.x) %>%
  # filter(!is.na(elevation_m_asl)) %>%
  summarise(max_elev = max(elevation, na.rm = TRUE),
            min_elev = min(elevation, na.rm = TRUE),
            elev_range = max_elev - min_elev) %>%
  ungroup()

meta <- env_file0 %>%
  rename(Study.x = Study) %>%
  distinct(Study.x, Taxa, Type_of_island, n_island_cat, number_of_islands) %>%
  # filter to studies in the model
  filter(Study.x %in% study_levels$Study.x) %>%
  left_join(elev_range) %>%
  ungroup()

# residuals
Sn_resid <- resid(Sn_lnorm,  type = 'pearson', method = 'posterior_predict') %>%
  as_tibble() %>%
  bind_cols(Sn_lnorm$data) %>%
  left_join(meta)

Spie_resid <- resid(Spie_lnorm,  type = 'pearson', method = 'posterior_predict') %>%
  as_tibble() %>%
  bind_cols(Spie_lnorm$data) %>%
  left_join(meta)

# Sn
png('./figures/Sn_resid.png',width = 1024, height = 768)
par(mfrow=c(2,3))
with(Sn_resid, boxplot(Estimate ~ Study.x));abline(h = 0, lty = 2)
with(Sn_resid, boxplot(Estimate ~ Taxa));abline(h = 0, lty = 2)
with(Sn_resid, boxplot(Estimate ~ Type_of_island));abline(h = 0, lty = 2)
with(Sn_resid, plot(Estimate ~ larea));abline(h = 0, lty = 2)
with(Sn_resid, plot(Estimate ~ number_of_islands));abline(h = 0, lty = 2)
with(Sn_resid, plot(Estimate ~ log10(elev_range)));abline(h = 0, lty = 2)
dev.off()

# Spie
png('./figures/Spie_resid.png',width = 1024, height = 768)
par(mfrow=c(2,3))
with(Spie_resid, boxplot(Estimate ~ Study.x));abline(h = 0, lty = 2)
with(Spie_resid, boxplot(Estimate ~ Taxa));abline(h = 0, lty = 2)
with(Spie_resid, boxplot(Estimate ~ Type_of_island));abline(h = 0, lty = 2)
with(Spie_resid, plot(Estimate ~ larea));abline(h = 0, lty = 2)
with(Spie_resid, plot(Estimate ~ number_of_islands));abline(h = 0, lty = 2)
with(Spie_resid, plot(Estimate ~ log10(elev_range)));abline(h = 0, lty = 2)
dev.off()

cowplot::plot_grid(pp_check(Sn_lnorm) +
  scale_x_continuous(name = expression(paste(S[n])),
                     trans = 'log2') +
    theme_minimal() +
    theme(legend.position = c(1,1),
          legend.justification = c(1,1),
          panel.border = element_rect(fill = NA, colour = 'black')),
pp_check(Spie_lnorm) +
  scale_x_continuous(name = expression(paste(S[PIE])),
                     trans = 'log2')+
  theme_minimal() +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        panel.border = element_rect(fill = NA, colour = 'black')),
  labels = 'auto', align = 'hv')

ggsave('./figures/posterior_density_check.png',
       width = 250, height = 130, units = 'mm')
