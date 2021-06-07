# code to get posterior samples for each study-level slope from multivariate model fit to
# ISAR data

library(brms)
library(tidyverse)
library(ggridges)

# load the model fit
load("./results/univariate_models_all_data.Rdata")

# study-levels
study_levels <- Sn_lnorm$data %>%
  as_tibble() %>%
  distinct(Study) %>%
  mutate(level = Study) %>%
  nest(data = c(level))

study_sample_posterior <- study_levels %>%
  mutate(Spie_global = purrr::map(data, ~posterior_samples(Spie_lnorm,
                                                           pars = 'b_larea',
                                                           fixed = TRUE) %>% unlist() %>% as.numeric()),
         Spie_study = purrr::map(data, ~posterior_samples(Spie_lnorm,
                                                          pars = paste('r_Study[', as.character(.x$level), ',larea]', sep=''),
                                                          fixed = TRUE) %>% unlist() %>% as.numeric()),
         Sn_global = purrr::map(data, ~posterior_samples(Sn_lnorm,
                                                         pars = 'b_larea',
                                                         fixed = TRUE) %>% unlist() %>% as.numeric()),
         Sn_study = purrr::map(data, ~posterior_samples(Sn_lnorm,
                                                        pars = paste('r_Study[', as.character(.x$level), ',larea]', sep=''),
                                                        fixed = TRUE) %>% unlist() %>% as.numeric()))

# load the meta data
env_file0 <- read_csv("./data/env_file_UTF-8.csv")

meta <- env_file0 %>%
  distinct(Study, Taxa, island_type2, n_island_cat) %>%
  # filter to studies in the model
  filter(Study %in% study_levels$Study) %>%
  mutate(simple_taxa = ifelse(Taxa=='Beetles', 'Invertebrates', Taxa))

# check that you have no NAs
study_sample_posterior <- left_join(study_sample_posterior,
                                    meta,
                                    by = 'Study') %>%
  unnest(cols = c(Spie_global, Spie_study,
                  Sn_global, Sn_study)) %>%
  select(-data)

# simplify island type: true (ocean / archipelego) or other

study_sample_posterior
study_sample_posterior <- study_sample_posterior %>%
  mutate(true_other = ifelse(island_type2=='volcanic', 'True','Other'))


# group posterior distributions for each taxa
post_taxa <-
  ggplot() +
  facet_wrap(~factor(metric,
                     levels = c('Rarefied richness',
                                'Evenness'),
                     labels = c(expression(paste(S[n])),
                                expression(paste(S[PIE])))),
             ncol = 2,
             scales = 'free_x', labeller = label_parsed) +
  # density of posteriors of Spie study-level slopes grouped by taxa
  geom_density_ridges_gradient(data = study_sample_posterior %>%
                                 mutate(metric='Evenness'),
                               aes(x = Spie_global + Spie_study,
                                   y = simple_taxa,
                                   # col="#fdbe85",
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9,
                               linetype = 0) +
  # repeat for Sn
  geom_density_ridges_gradient(data = study_sample_posterior %>%
                                 mutate(metric='Rarefied richness'),
                               aes(x = Sn_global + Sn_study,
                                   y = simple_taxa,
                                   # col = "#fdbe85",
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9,
                               linetype = 0) +
  geom_vline(xintercept = 0, lty = 2) +
  # global estimates
  geom_rect(data = study_sample_posterior %>%
              mutate(metric = 'Evenness',
                     lower = quantile(Spie_global, probs = c(0.025)),
                     upper = quantile(Spie_global, probs = c(0.975))) %>%
              distinct(metric, lower, upper),
            aes(xmin = lower, xmax = upper, ymin= -Inf, ymax = Inf),
            alpha = 0.5) +
  geom_rect(data = study_sample_posterior %>%
              mutate(metric = 'Rarefied richness',
                     lower = quantile(Sn_global, probs = c(0.025)),
                     upper = quantile(Sn_global, probs = c(0.975))) %>%
              distinct(metric, lower, upper),
            aes(xmin = lower, xmax = upper, ymin= -Inf, ymax = Inf),
            alpha = 0.5) +
  geom_vline(data = study_sample_posterior %>%
               mutate(metric = 'Evenness'),
             aes(xintercept = median(Spie_global))) +
  geom_vline(data = study_sample_posterior %>%
               mutate(metric = 'Rarefied richness'),
             aes(xintercept = median(Sn_global))) +
  # add point for median of posterior distribution
  geom_point(data = study_sample_posterior %>%
               mutate(metric='Rarefied richness'),
             aes(x = Sn_global + Sn_study,
                 y = simple_taxa),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 2.5, shape = 17) +
  # add point for median of posterior distribution
  geom_point(data = study_sample_posterior %>%
               mutate(metric='Evenness'),
             aes(x = Spie_global + Spie_study,
                 y = simple_taxa),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 2.5, shape = 17) +
  geom_text(data = study_sample_posterior %>%
              group_by(simple_taxa) %>%
              summarise(n_study = n_distinct(Study)) %>%
              ungroup() %>%
              distinct(simple_taxa, n_study, .keep_all = T) %>%
              mutate(metric = 'Rarefied richness'),
            aes(x=0.4, y=simple_taxa,
                label=paste('n[study] == ', n_study)),
            size=4.1,
            nudge_y = 0.1, parse = T) +
  labs(y = 'Taxon group',
       x = 'Study-level slopes') +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#bdc9e1', '#74a9cf', '#2b8cbe',
                               '#74a9cf', '#bdc9e1')) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = 'black'),
        legend.key = element_blank(),
        legend.position = 'none',
        strip.text = element_text(hjust = 0, size = 12),
        # axis.title.y = element_blank(),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        aspect.ratio = 1)

# group posterior distributions for island types
post_type <-
  ggplot() +
  facet_wrap(~factor(metric,
                     levels = c('Sn', 'Spie'),
                     labels = c(expression(paste(S[n])),
                                expression(paste(S[PIE])))),
             ncol = 2,
             scales = 'free_x', labeller = label_parsed) +
  # density of posteriors of Spie study-level slopes grouped by taxa
  geom_density_ridges_gradient(data = study_sample_posterior %>%
                                 mutate(metric = 'Spie'),
                               aes(x = Spie_global + Spie_study,
                                   y = true_other,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, linetype = 0) +
  # repeat for Sn
  geom_density_ridges_gradient(data = study_sample_posterior %>%
                                 mutate(metric = 'Sn'),
                               aes(x = Sn_global + Sn_study,
                                   y = true_other,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.025, 0.25, 0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, linetype = 0) +
  # global estimates
  geom_rect(data = study_sample_posterior %>%
              mutate(metric = 'Spie',
                     lower = quantile(Spie_global, probs = c(0.025)),
                     upper = quantile(Spie_global, probs = c(0.975))) %>%
              distinct(metric, lower, upper),
            aes(xmin = lower, xmax = upper, ymin= -Inf, ymax = Inf),
            alpha = 0.5) +
  geom_rect(data = study_sample_posterior %>%
              mutate(metric = 'Sn',
                     lower = quantile(Sn_global, probs = c(0.025)),
                     upper = quantile(Sn_global, probs = c(0.975))) %>%
              distinct(metric, lower, upper),
            aes(xmin = lower, xmax = upper, ymin= -Inf, ymax = Inf),
            alpha = 0.5) +
  geom_vline(data = study_sample_posterior %>%
               mutate(metric = 'Spie'),
             aes(xintercept = median(Spie_global))) +
  geom_vline(data = study_sample_posterior %>%
               mutate(metric = 'Sn'),
             aes(xintercept = median(Sn_global))) +
  # add point for median of posterior distribution
  geom_point(data = study_sample_posterior %>%
               mutate(metric = 'Sn'),
             aes(x = Sn_global + Sn_study,
                 y = true_other),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 2.5, shape = 17) +
  # add point for median of posterior distribution
  geom_point(data = study_sample_posterior %>%
               mutate(metric = 'Spie'),
             aes(x = Spie_global + Spie_study,
                 y = true_other),
             stat = ggstance:::StatSummaryh,
             fun.x = median,
             size = 2.5, shape = 17) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_text(data = study_sample_posterior %>%
              group_by(true_other) %>%
              summarise(n_study = n_distinct(Study)) %>%
              ungroup() %>%
              distinct(true_other, n_study, .keep_all = T) %>%
              mutate(metric = 'Sn'),
            aes(x=0.4, y=true_other,
                label=paste('n[study] == ', n_study)),
            size=4.2,
            nudge_y = 0.1, parse = T) +
  theme_bw() +
  labs(y = 'Island type',
       x = expression(paste('Study-level slopes'))
       # subtitle = expression(paste('Posterior samples of study-level ', S[std], ' fragment area slopes'))#,
       # tag = 'c'
  ) +
  scale_y_discrete(labels = scales::wrap_format(12), expand = c(0.05,0,0.1,0)) +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#bdc9e1', '#74a9cf', '#2b8cbe',
                               '#74a9cf', '#bdc9e1')) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(fill = NA, colour = 'black'),
        legend.key = element_blank(),
        legend.position = 'none',
        strip.text = element_text(hjust = 0, size = 12),
        # axis.title.y = element_blank(),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        aspect.ratio = 1)


## dummy plot for creating separate legend
three_grey_legend <- ggplot() +
  # facet_grid(continent ~ ., scale = 'free') +
  geom_density_ridges_gradient(data = study_sample_posterior %>%
                                 mutate(metric = 'Sn'),
                               aes(x = Sn_global + Sn_study,
                                   y = true_other,
                                   fill = stat(quantile)
                               ),
                               quantiles = c(0.75, 0.975),
                               calc_ecdf = T,
                               scale = 0.9, linetype = 0) +
  scale_fill_manual(name = 'Posterior probability',
                    values = c('#bdc9e1', '#74a9cf', '#2b8cbe'),
                    labels = c('< 5%', '< 45%',  '50%')) +
  theme(panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        # legend.justification = c(1, 1),
        legend.background = element_blank(),
        #legend.text = element_text(size = 5, face = 'plain'),
        # legend.title = element_text(size = 6, face = 'plain'),
        legend.margin = margin(),
        legend.box.spacing = unit(c(0,0,0,0), units = 'mm'),
        # legend.key.size = unit(2, units = 'mm'),
        plot.margin = unit(c(0,0,0,0), units = 'mm')) #+

source('./functions/gg_legend.R')
legend <- gg_legend(three_grey_legend)

cowplot::plot_grid(post_taxa,
                   post_type,
                   legend,
                   nrow = 3,
                   rel_heights = c(1, 1, 0.05),
                   labels = c('a', 'b', ''))

ggsave("./figures/Fig4.png",
       width = 250, height = 250, units = 'mm')
