# load necessary libraries
library(tidyverse)
library(RCurl)
library(png)
library(mapproj)

############################
# Set path and directories #
############################

work_dir <- getwd() # first set working directory in Menu/Session/Set working directory/to Project
data_path <- paste(work_dir, "/data", sep = "")
figures_path <- paste(work_dir, "/figures", sep = "")


# load the dataset coordinates data

env <-  read_csv(paste0(data_path, "/env_file_UTF-8.csv"))
coord_data <- data.frame(
  unique(env[, c("study_ID","number_of_islands","island_type","Taxa")]),
  aggregate(env[, c("island_lat","island_lon")], by = list(study_ID = env$study_ID), mean, na.rm = TRUE)[, -1L]
)
colnames(coord_data) <- c("Dataset_ID","N_islands","island_type","Taxa","lat","long")

# load the map
world <- map_data('world') %>%
  as_tibble()

world <- world %>%
  rename(
    Longitude = long,
    Latitude = lat
  )

map_taxa <- ggplot() +
  geom_polygon(data=world,
               aes(Longitude, Latitude, group = group), colour=NA, fill='#CCCCCC', size=0) +

  geom_point(data = coord_data,
             aes(x = long, y = lat, shape = island_type, colour = Taxa, size = N_islands)) +
  coord_map('mollweide', ylim = c(-60, 90), xlim = c(-180, 180)) +
  scale_x_continuous(breaks = seq(-180, 180, by = 30)) +
  scale_y_continuous(breaks = c(0, -23.5, 23.5, -60, 60)) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = 'black', size = 0.1),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = 'bottom',
        legend.direction = 'horizontal',
        legend.text=element_text(size=13),
        legend.title = element_text(color="#636363", size = 15, face = "bold"),

        plot.margin = unit(c(0,0,0,0), units = 'mm')) +
  guides(colour = guide_legend(title.position = "top", title.hjust = 0.5,nrow = 2, override.aes = list(size=7)), shape =  guide_legend(title.position = "top", title.hjust = 0.5, nrow = 2,override.aes = list(size=7)), size =  guide_legend(title.position = "top",title.hjust = 0.5)) +
  scale_colour_manual(name="Taxa",values = c("Beetles"= "#d95f02", "Herpetofauna" = "#31a354",
                                             "Mammals" = "#386cb0", "Birds"= "#984ea3",
                                             "Invertebrates"= "#fecc5c", "Plants" = "#66c2a4")) +
  scale_shape_manual(name="Island Type",values = c("Atoll"= 8, "Barrier Island" = 15,
                                                   "True Island" = 18, "Forest Island"= 17 ,
                                                   "Lake Island"= 16 ))+
  labs(shape="island_type", colour="Taxa", size = "N_islands")

ggsave((paste(figures_path, "/dataset_map.png",sep="")),
       width = 340, height = 200, units = 'mm')
