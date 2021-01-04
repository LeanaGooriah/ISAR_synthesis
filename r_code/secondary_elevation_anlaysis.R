# code to model study-level variation as a function of elevation range
library(brms)
library(tidyverse)

# load the model fit
load('~/Dropbox/ISAR Meta-analysis/new_Sept2020/results/univariate_models.Rdata')

# load the meta data: this is not very clean, and the join with the DEM elevation data is not working
env_file0 <- read_csv("~/Dropbox/ISAR Meta-analysis/new_Sept2020/data/env_file_august.xlsx - new.csv")

# investigate elevation data: use DEM for all studies except: Andrade_2002, Bell_2017, Dasilva_2019, 
#Jonsson_2011, Macdonald_2018, Pereira_2017, Perillo_2017, Werden_2012
# Werner_Zalewski_2006

DEM_filter <- tibble(study_ID = c('Andrade_2002', 'Bell_2017', 'Dasilva_2019', 
                                  'Jonsson_2011', 'Macdonald_2018', 'Pereira_2017', 'Perillo_2017', 'Werden_2012',
                                  'Werner_Zalewski_2006'))

elev <- read.delim("~/Dropbox/ISAR Meta-analysis/new_Sept2020/data/env_file_alban.csv", sep = ';') %>% 
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
  distinct(Study.x, Taxa, Type_of_island, n_island_cat) %>% 
  # filter to studies in the model
  filter(Study.x %in% study_levels$Study.x) %>% 
  left_join(elev_range) %>% 
  ungroup() %>% 
  mutate(simple_taxa = ifelse(Taxa=='Beetles', 'Invertebrates', Taxa))


Sn_coefs <- coef(Sn_lnorm)$Study.x[,,'larea'] %>% 
  as_tibble() %>% 
  mutate(Study.x = rownames(coef(Sn_lnorm)$Study.x[,,'larea'])) %>% 
  left_join(meta) %>% 
  filter(elev_range!=-Inf) %>% 
  mutate(cElevRange = elev_range - mean(elev_range),
         cElevRange_log = log(elev_range+1) - mean(log(elev_range+1)))

Spie_coefs <- coef(Spie_lnorm)$Study.x[,,'larea'] %>% 
  as_tibble() %>% 
  mutate(Study.x = rownames(coef(Sn_lnorm)$Study.x[,,'larea'])) %>% 
  left_join(meta) %>% 
  filter(elev_range!=-Inf) %>% 
  mutate(cElevRange = elev_range - mean(elev_range),
         cElevRange_log = log(elev_range+1) - mean(log(elev_range+1)))

# fit second stage model (include uncertainty in study-level slope estimate)
Sn_elev_log <- brm(bf(Estimate | se(Est.Error) ~ cElevRange_log + (1 | Study.x)),
                   data = Sn_coefs, 
                   cores = 4)

# Spie ~  elev_range on log-scale
Spie_elev_log <- brm(bf(Estimate | se(Est.Error) ~ cElevRange_log + (1 | Study.x)),
                     data = Spie_coefs, 
                     cores = 4)

plot(Sn_elev_log)
plot(Spie_elev_log)

pp_check(Sn_elev_log)
pp_check(Spie_elev_log)

# get coefficients for plotting model results
Sn_elev_fit <- cbind(Sn_coefs %>% 
                       mutate(Slope = Estimate,
                              Slope_error = Est.Error,
                              S_Q2.5 = Q2.5,
                              S_Q97.5 = Q97.5) %>% 
                       select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                     fitted(Sn_elev_log, re_formula = NA)) %>% 
  as_tibble()

Sn_elev_global <- fixef(Sn_elev_log)

Spie_elev_fit <- cbind(Spie_coefs %>% 
                         mutate(Slope = Estimate,
                                Slope_error = Est.Error,
                                S_Q2.5 = Q2.5,
                                S_Q97.5 = Q97.5) %>% 
                         select(-Estimate, -Est.Error, -Q2.5, -Q97.5),
                       fitted(Spie_elev_log, re_formula = NA)) %>% 
  as_tibble()

Spie_elev_global <- fixef(Spie_elev_log)

# elevation range
Sn_elev_plot <- ggplot() +
  facet_wrap(~factor(metric,
                     levels = c('Sn', 'Spie'),
                     labels = c(expression(paste(S[n])),
                                expression(paste(S[PIE])))),
             ncol = 2,
             labeller = label_parsed) +
  geom_point(data = Sn_elev_fit %>% 
               mutate(metric = 'Sn'),
             aes(x = elev_range + 1, y = Slope)) +
  geom_linerange(data = Sn_elev_fit %>% 
                   mutate(metric = 'Sn'),
                 aes(x = elev_range + 1, ymin = S_Q2.5, ymax = S_Q97.5)) +
  geom_ribbon(data = Sn_elev_fit %>% 
                mutate(metric = 'Sn'),
              aes(x = elev_range + 1, ymin = Q2.5, ymax = Q97.5),
              alpha = 0.5) +
  geom_line(data = Sn_elev_fit %>% 
              mutate(metric = 'Sn'),
            aes(x = elev_range + 1, y = Estimate)) +
  # add regression coefficient and uncertainty interval
  annotate('text', x = 1, y = Inf, hjust = 0.1, vjust = 1.4,
           label = paste("beta == ", #[Frag.~size]
                         round(Sn_elev_global['cElevRange_log','Estimate'],3),
                         "  (",
                         round(Sn_elev_global['cElevRange_log','Q2.5'],3),
                         " - ",
                         round(Sn_elev_global['cElevRange_log','Q97.5'],3),
                         ")"),  
           parse = T#, size = 2
  ) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(y = 'Study-level slope estimate',
       x = '') +
  theme_minimal() +
  theme(#panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, colour = 'black'),
    legend.key = element_blank(),
    legend.position = 'none', 
    strip.text = element_text(hjust = 0, size = 12),
    # axis.title.y = element_blank(),
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    aspect.ratio = 1)

# Spie
Spie_elev_plot <- ggplot() +
  facet_wrap(~factor(metric,
                     levels = c('Sn', 'Spie'),
                     labels = c(expression(paste(S[n])),
                                expression(paste(S[PIE])))),
             ncol = 2,
             labeller = label_parsed) +
  geom_point(data = Spie_elev_fit %>% 
               mutate(metric = 'Spie'),
             aes(x = elev_range + 1, y = Slope)) +
  geom_linerange(data = Spie_elev_fit %>% 
                   mutate(metric = 'Spie'),
                 aes(x = elev_range + 1, ymin = S_Q2.5, ymax = S_Q97.5)) +
  geom_ribbon(data = Spie_elev_fit %>% 
                mutate(metric = 'Spie'),
              aes(x = elev_range + 1, ymin = Q2.5, ymax = Q97.5),
              alpha = 0.5) +
  geom_line(data = Spie_elev_fit %>% 
              mutate(metric = 'Spie'),
            aes(x = elev_range + 1, y = Estimate)) +
  # add regression coefficient and uncertainty interval
  annotate('text', x = 1, y = Inf, hjust = 0.1, vjust = 1.4,
           label = paste("beta == ", #[Frag.~size]
                         round(Spie_elev_global['cElevRange_log','Estimate'],3),
                         "  (",
                         round(Spie_elev_global['cElevRange_log','Q2.5'],3),
                         " - ",
                         round(Spie_elev_global['cElevRange_log','Q97.5'],3),
                         ")"),  
           parse = T#, size = 2
  ) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(y = '',
       x = '') +
  theme_minimal() +
  theme(#panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, colour = 'black'),
    legend.key = element_blank(),
    legend.position = 'none', 
    strip.text = element_text(hjust = 0, size = 12),
    # axis.title.y = element_blank(),
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    aspect.ratio = 1)

cowplot::plot_grid(Sn_elev_plot,
                   Spie_elev_plot,
                   labels = 'auto',
                   align = 'hv') +
  cowplot::draw_label(y = 0.04, label = 'Elevation range (m)')


ggsave('~/Dropbox/ISAR Meta-analysis/new_Sept2020/figures/fig4_log_scale.png', 
       width = 250, height = 130, units = 'mm')


# check rhats of main models
Sn_rhat <- rhat(Sn_lnorm)
Spie_rhat <- rhat(Spie_lnorm)
bayesplot::mcmc_rhat(Sn_rhat)
bayesplot::mcmc_rhat(Spie_rhat)
