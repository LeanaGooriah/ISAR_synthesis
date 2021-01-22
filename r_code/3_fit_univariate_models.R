# Code for model analysis of ISARs using brms
# Default, weakly regularising priors were used for all parameters

# load necessary packages

# install.packages("tidyverse")
# install.packages("brms")
library(tidyverse)
library(brms)

# load the data
work_dir <- getwd()
data_path <- paste(work_dir,"/data",sep="")
load(paste(result_path, "/diversity_global_synthesis_with_area_NEW.Rdata",sep=""))

# log and mean centre island area to use as a predictor
temp = all_ISAR %>%
  filter(!is.na(island_area_ha))

temp = temp %>%
  filter(island_area_ha != 0)

temp$larea <- log(temp$island_area_ha) - mean(log(temp$island_area_ha))

## univariate responses (same as Chase et al 2020 Nature analyses)
Sn = temp %>%
  filter(index=='S_n') %>%
  select(Study.x, value, larea, island_area_ha, island_code, Sn = value)

Spie = temp %>%
  filter(index=='ENS') %>%
  select(Study.x, value, larea,island_area_ha, island_code, Spie = value)


#  non-varying intercept and slope (with study-level variation)
Sn_lnorm <- brm(Sn ~ larea + (larea | Study.x),
                family = lognormal(),
                data = Sn,
                cores = 4, chains = 4)

pp_check(Sn_lnorm) +
  scale_x_continuous(trans = 'log2')

Spie_lnorm <- brm(Spie ~ larea + (larea | Study.x),
                family = lognormal(),
                data = Spie,
                cores = 4, chains = 4)

pp_check(Spie_lnorm) +
  scale_x_continuous(trans = 'log2')

save(Sn_lnorm, Spie_lnorm, file = (paste(getwd(), "/results/univariate_models.Rdata",sep="")))



