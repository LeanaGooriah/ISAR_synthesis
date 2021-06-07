# Code for model analysis of ISARs using brms
# Default, weakly regularising priors were used for all parameters

# load necessary packages

library(tidyverse)
library(brms)
library(Rcpp)
library(StanHeaders)

############################
# Set path and directories #
############################
work_dir <- getwd() # first set working directory in Menu/Session/Set working directory/to Project

data_path <- paste0(work_dir, "./data")
result_path <- paste0(work_dir, "./results")

# 1) load the data (full dataset)
#load("C:/Users/leana/Desktop/ISAR_datasets/diversity_global_synthesis_all_data.Rdata")

all_ISAR <- readRDS('./data/diversity_global_synthesis_all_data.Rds')
all_ISAR <- all_ISAR %>% rename(Study = Study.x)

# log and mean centre island area to use as a predictor
temp = all_ISAR %>%
  filter(!is.na(island_area_ha))

temp = temp %>%
  filter(island_area_ha != 0)

temp$larea <- log(temp$island_area_ha) - mean(log(temp$island_area_ha))

## univariate responses (same as Chase et al 2020 Nature analyses)
Sn = temp %>%
  filter(index=='S_n') %>%
  select(Study, value, larea, island_area_ha, island_code, Sn = value)

Spie = temp %>%
  filter(index=='ENS') %>%
  select(Study, value, larea,island_area_ha, island_code, Spie = value)


#  non-varying intercept and slope (with study-level variation)
Sn_lnorm <- brm(Sn ~ larea + (larea | Study),
                family = lognormal(),
                data = Sn,
                cores = 4, chains = 4, iter = 3000)

pp_check(Sn_lnorm) +
  scale_x_continuous(trans = 'log2')


Spie_lnorm <- brm(Spie ~ larea + (larea | Study),
                  family = lognormal(),
                  data = Spie,
                  cores = 4, chains = 4, iter = 3000)

pp_check(Spie_lnorm) +
  scale_x_continuous(trans = 'log2')


fixef(Sn_lnorm)
bayes_R2(Sn_lnorm)
fixef(Spie_lnorm)
bayes_R2(Spie_lnorm)

save(Sn_lnorm, Spie_lnorm, file = (paste(result_path, "/univariate_models_all_data.Rdata",sep="")))


# 2) load the data (prop sampling only)
all_ISAR_prop <- readRDS('./data/diversity_global_synthesis_prop_sampling.Rds')

# log and mean centre island area to use as a predictor
temp2 = all_ISAR_prop %>%
  filter(!is.na(island_area_ha))

temp2 = temp2 %>%
  filter(island_area_ha != 0)

temp2$larea <- log(temp2$island_area_ha) - mean(log(temp2$island_area_ha))

## univariate responses (same as Chase et al 2020 Nature analyses)
Sn_prop = temp2 %>%
  filter(index=='S_n') %>%
  select(Study, value, larea, island_area_ha, island_code, Sn = value)

Spie_prop = temp2 %>%
  filter(index=='ENS') %>%
  select(Study, value, larea,island_area_ha, island_code, Spie = value)


#  non-varying intercept and slope (with study-level variation)
Sn_lnorm_prop <- brm(Sn ~ larea + (larea | Study),
                     family = lognormal(),
                     data = Sn_prop,
                     cores = 4, chains = 4, iter = 3000)

pp_check(Sn_lnorm_prop) +
  scale_x_continuous(trans = 'log2')


Spie_lnorm_prop <- brm(Spie ~ larea + (larea | Study),
                       family = lognormal(),
                       data = Spie_prop,
                       cores = 4, chains = 4, iter = 3000)

pp_check(Spie_lnorm_prop) +
  scale_x_continuous(trans = 'log2')


bayes_R2(Sn_lnorm_prop)
bayes_R2(Spie_lnorm_prop)
fixef(Sn_lnorm_prop)
fixef(Spie_lnorm_prop)

save(Sn_lnorm_prop, Spie_lnorm_prop, file = (paste(result_path, "/univariate_models_prop.Rdata",sep="")))


# 3) load the data (islands with N > 10)
all_ISAR_N10 <- readRDS('./data/diversity_global_synthesis_N10.Rds')
all_ISAR_N10 <- all_ISAR_N10 %>% rename(Study = Study.x)

# log and mean centre island area to use as a predictor
temp3 = all_ISAR_N10 %>%
  filter(!is.na(island_area_ha))

temp3 = temp3 %>%
  filter(island_area_ha != 0)

temp3$larea <- log(temp3$island_area_ha) - mean(log(temp3$island_area_ha))

## univariate responses (same as Chase et al 2020 Nature analyses)
Sn_N10 = temp3 %>%
  filter(index=='S_n') %>%
  select(Study, value, larea, island_area_ha, island_code, Sn = value)

Spie_N10 = temp3 %>%
  filter(index=='ENS') %>%
  select(Study, value, larea,island_area_ha, island_code, Spie = value)


#  non-varying intercept and slope (with study-level variation)
Sn_lnorm_N10 <- brm(Sn ~ larea + (larea | Study),
                    family = lognormal(),
                    data = Sn_N10,
                    cores = 4, chains = 4, iter = 3000)

pp_check(Sn_lnorm_N10) +
  scale_x_continuous(trans = 'log2')


Spie_lnorm_N10 <- brm(Spie ~ larea + (larea | Study),
                      family = lognormal(),
                      data = Spie_N10,
                      cores = 4, chains = 4, iter = 3000)

pp_check(Spie_lnorm_N10) +
  scale_x_continuous(trans = 'log2')


bayes_R2(Sn_lnorm_N10)
bayes_R2(Spie_lnorm_N10)
fixef(Sn_lnorm_N10)
fixef(Spie_lnorm_N10)
save(Sn_lnorm_N10, Spie_lnorm_N10, file = (paste(result_path, "/univariate_models_N10.Rdata",sep="")))


# 4) load the data (islands with N > 20)

all_ISAR_N20 <- readRDS('./data/diversity_global_synthesis_N20.Rds')
all_ISAR_N20 <- all_ISAR_N20 %>% rename(Study = Study.x)

# log and mean centre island area to use as a predictor
temp4 = all_ISAR_N20 %>%
  filter(!is.na(island_area_ha))

temp4 = temp4 %>%
  filter(island_area_ha != 0)

temp4$larea <- log(temp4$island_area_ha) - mean(log(temp4$island_area_ha))

## univariate responses (same as Chase et al 2020 Nature analyses)
Sn_N20 = temp4 %>%
  filter(index=='S_n') %>%
  select(Study, value, larea, island_area_ha, island_code, Sn = value)

Spie_N20 = temp4 %>%
  filter(index=='ENS') %>%
  select(Study, value, larea,island_area_ha, island_code, Spie = value)


#  non-varying intercept and slope (with study-level variation)
Sn_lnorm_N20 <- brm(Sn ~ larea + (larea | Study),
                    family = lognormal(),
                    data = Sn_N20,
                    cores = 4, chains = 4, iter = 3000)

pp_check(Sn_lnorm_N20) +
  scale_x_continuous(trans = 'log2')


Spie_lnorm_N20 <- brm(Spie ~ larea + (larea | Study),
                      family = lognormal(),
                      data = Spie_N20,
                      cores = 4, chains = 4, iter = 3000)

pp_check(Spie_lnorm_N20) +
  scale_x_continuous(trans = 'log2')

fixef(Sn_lnorm_N20)
bayes_R2(Sn_lnorm_N20)
bayes_R2(Spie_lnorm_N20)
fixef(Spie_lnorm_N20)

save(Sn_lnorm_N20, Spie_lnorm_N20, file = (paste(result_path, "/univariate_models_N20.Rdata",sep="")))

