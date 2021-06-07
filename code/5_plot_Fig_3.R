# load packages

require(dplyr)
require(tidyr)
library(devtools)
library(tidyr)
library(ggplot2)
library(brms)


#load necessary R data

all_ISAR <- readRDS("./data/diversity_global_synthesis_all_data.Rds")
all_ISAR <- all_ISAR %>% rename(Study = Study.x)
load("./results/univariate_models_all_data.Rdata")

# data that models were fit to:
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

# Sn and Spie: study (archipelago) level estimates
Sn_study_coefs <- coef(Sn_lnorm)
Spie_study_coefs <- coef(Spie_lnorm)
study_coefs <- tibble(Study = rownames (Sn_study_coefs$Study),
                      Spie_intercept = Spie_study_coefs$Study[, , "Intercept"][,'Estimate'],
                      Spie_slope = Spie_study_coefs$Study[, , "larea"][,'Estimate'],
                      Spiemin = Spie_study_coefs$Study[, , "larea"][,'Q2.5'],
                      Spiemax = Spie_study_coefs$Study[, , "larea"][,'Q97.5'],
                      Sn_intercept = Sn_study_coefs$Study[, , "Intercept"][,'Estimate'],
                      Sn_slope = Sn_study_coefs$Study[, , "larea"][,'Estimate'],
                      Snmin = Sn_study_coefs$Study[, , "larea"][,'Q2.5'],
                      Snmax = Sn_study_coefs$Study[, , "larea"][,'Q97.5'])

study_coefs <- study_coefs %>% mutate(hypothesis = ifelse(Snmin > 0 & Spiemin <= 0, "Rare species",
                                                          ifelse(Snmin > 0 & Spiemin > 0, "Rare and common", "Sampling effects")))

study_coefs %>%
  group_by(hypothesis) %>%
  summarise(n())

# global parameters
Sn_global <- fixef(Sn_lnorm)

Sn_fitted <- cbind(Sn_lnorm$data,
                   fitted(Sn_lnorm, re_formula = NA)) %>%
  as_tibble() %>%
  inner_join(Sn %>%
               distinct(Study, larea, island_area_ha))

Spie_global <- fixef(Spie_lnorm)

Spie_fitted <- cbind(Spie_lnorm$data,
                     fitted(Spie_lnorm, re_formula = NA)) %>%
  as_tibble() %>%
  inner_join(Spie %>%
               distinct(Study, larea, island_area_ha))

# load the meta data

env_file <- read_csv("./data/env_file_UTF-8.csv")


env <- env_file %>%
  select(Study, island_code,island_area_ha, island_type2, Taxa) %>%
  filter(island_area_ha!=0) %>%
  mutate(log_a = log(island_area_ha) - mean(log(island_area_ha)),
         # put beetles back in with the other invertebrates
         simple_taxa = ifelse(Taxa=='Beetles', 'Invertebrates', Taxa))

env_tot <- env %>%
  rename(Study = Study) %>%
  group_by(Study) %>%
  summarise(xmin = min(island_area_ha),
            xmax = max(island_area_ha),
            cxmin = min(log_a),
            cxmax = max(log_a),
            Type_of_island = unique(island_type2),
            Taxa = unique(Taxa),
            simple_taxa = unique(simple_taxa))

hyp_colour = c('Rare and common' = '#7570b3',
               'Rare species' = '#d95f02',
               'Sampling effects' = '#1b9e77')

Spie_regPlot <-
  ggplot() +
  # data
  geom_point(data = left_join(Spie %>%
                                mutate(Study = Study),
                              env_tot %>%
                                distinct(Study, Taxa, simple_taxa)),
             aes(x = island_area_ha, y = Spie, colour = simple_taxa),
             size = 1, alpha = 0.25) +
  geom_segment(data = left_join(study_coefs %>%
                                  rename(Study = Study),
                                env_tot),
               aes(group = Study,
                   colour = simple_taxa,
                   x = xmin,
                   xend = xmax,
                   y = exp(Spie_intercept + Spie_slope * cxmin),
                   yend = exp(Spie_intercept + Spie_slope * cxmax)),
               size = 0.5) +
  # fixed effect
  geom_line(data = Spie_fitted,
            aes(x = island_area_ha,
                y = Estimate),
            size = 1) +
  # fixed effect uncertainty
  geom_ribbon(data = Spie_fitted,
              aes(x = island_area_ha,
                  ymin = Q2.5,
                  ymax = Q97.5),
              alpha = 0.3) +
  # add regression coefficient and uncertainty interval
  annotate('text', x = 0.02, y = Inf, hjust = 0.1, vjust = 1.4,
           label = paste("beta == ", #[Frag.~size]
                         round(Spie_global['larea','Estimate'],2),
                         "  (",
                         round(Spie_global['larea','Q2.5'],2),
                         " - ",
                         round(Spie_global['larea','Q97.5'],2),
                         ")"),
           parse = T#, size = 2
  ) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(1, 2, 4, 8, 16,32,64)) +
  # scale_color_manual(name = '', values = hyp_colour, guide = F) +
  scale_colour_viridis_d(guide = F) +
  labs(x = '',
       y = expression(paste(S[PIE], ' (evenness)'))#,
       # tag = 'd'
  ) +
  theme_bw() +
  theme(legend.position = c(0.5, 0.9),
        legend.direction = 'horizontal',
        legend.background = element_blank()#,
        # text = element_text(size = 7),
        # plot.tag = element_text(size = 8, face = 'bold')
  )

Sn_regPlot <-
  ggplot() +
  # data
  geom_point(data = left_join(Sn %>%
                                mutate(Study = Study),
                              env_tot %>%
                                distinct(Study, Taxa, simple_taxa)),
             aes(x = island_area_ha, y = Sn, colour = simple_taxa),
             size = 1, alpha = 0.25) +
  geom_segment(data = left_join(study_coefs %>%
                                  rename(Study = Study),
                                env_tot),
               aes(group = Study,
                   colour = simple_taxa,
                   x = xmin,
                   xend = xmax,
                   y = exp(Sn_intercept + Sn_slope * cxmin),
                   yend = exp(Sn_intercept + Sn_slope * cxmax)),
               size = 0.5) +
  # fixed effect
  geom_line(data = Sn_fitted,
            aes(x = island_area_ha,
                y = Estimate),
            size = 1) +
  # fixed effect uncertainty
  geom_ribbon(data = Sn_fitted,
              aes(x = island_area_ha,
                  ymin = Q2.5,
                  ymax = Q97.5),
              alpha = 0.3) +
  # add regression coefficient and uncertainty interval
  annotate('text', x = 0.02, y = Inf, hjust = 0.1, vjust = 1.4,
           label = paste("beta == ", #[Frag.~size]
                         round(Sn_global['larea','Estimate'],2),
                         "  (",
                         round(Sn_global['larea','Q2.5'],2),
                         " - ",
                         round(Sn_global['larea','Q97.5'],2),
                         ")"),
           parse = T#, size = 2
  ) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log', breaks = c(1, 2, 4, 8, 16,32,64, 128, 256)) +
  # scale_color_manual(name = '', values = hyp_colour, guide = F) +
  scale_colour_viridis_d(guide = F) +
  labs(x = '',
       y = expression(paste(S[n], ' (rarefied richness)'))#,
       # tag = 'd'
  ) +
  theme_bw() +
  theme(legend.position = c(0.5, 0.9),
        legend.direction = 'horizontal',
        legend.background = element_blank()#,
        # text = element_text(size = 7),
        # plot.tag = element_text(size = 8, face = 'bold')
  )

# separate colour legend
col_legend <- ggplot() +
  geom_point(data = left_join(Sn %>%
                                mutate(Study = Study),
                              env_tot %>%
                                distinct(Study, Taxa, simple_taxa)),
             aes(x = island_area_ha, y = Sn, colour = simple_taxa),
             size = 1, alpha = 0.25)+
  geom_segment(data = left_join(study_coefs %>%
                                  rename(Study = Study),
                                env_tot),
               aes(group = Study,
                   colour = simple_taxa,
                   x = xmin,
                   xend = xmax,
                   y = exp(Sn_intercept + Sn_slope * cxmin),
                   yend = exp(Sn_intercept + Sn_slope * cxmax)),
               size = 0.5) +
  # scale_color_manual(name = 'Hypothesis', values = hyp_colour) +
  scale_color_viridis_d(name = 'Taxa') +
  theme_bw() +
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.background = element_blank()) +
  guides(colour = guide_legend(nrow = 1))

source('./functions/gg_legend.R')
gg_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}
col_legend2 <- gg_legend(col_legend)

cowplot::plot_grid(cowplot::plot_grid(Sn_regPlot,
                                      Spie_regPlot, labels = 'auto'),
                   col_legend2,
                   rel_heights = c(1, 0.05),
                   ncol = 1) +
  cowplot::draw_label(y = 0.1,
                      label = 'Island area (ha)')
ggsave("./figures/Fig3.png", width = 200, height = 90, units = 'mm')
