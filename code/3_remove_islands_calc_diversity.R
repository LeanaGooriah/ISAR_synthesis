#############################################################
# CODE TO REMOVE ISLANDS ACCORDING TO NUMBER OF INDIVIDUALS #
#############################################################
rm(list = ls())
require(dplyr)
require(tidyr)
library(devtools)
library(vegan)
library(tidyr)
library(readr)
# install_github('MoBiodiv/mobr') # installing version saved in renv.lock is preferred
require(mobr) # version 1.0

work_dir <- getwd() # first set working directory in Menu/Session/Set working directory/to Project


# main data path
data_path <- paste0(work_dir, "/data/isar_datasets")


# data paths to save updated datasets according to number of individuals
data_path1 <- paste(data_path,"/ISAR_datasets_N10",sep="")
data_path2 <- paste(data_path,"/ISAR_datasets_N20",sep="")

# Read all data file names
filenames <- list.files(data_path, pattern="*.csv*", full.names = F)

# Make list of study ids
filename_roots <- gsub(".csv","",filenames) # remove .csv from filenames

study_ids <- unique(filename_roots)
study_ids


###################
# 1.Read in data  #
###################

n_files <- length(unique(study_ids))
data_out <- list()

for (i in 1:n_files){
  data_file <- paste(data_path,"/", study_ids[i], ".csv", sep ="")
  data_out[[i]] <- read.csv(data_file, header = TRUE, stringsAsFactors = F)
}

#################################
# 2. Loop through all
# csv files and remove islands
# where total N is < 10 and < 20
#################################

div_out_N10 <- list()
div_out_N20 <- list()

for (i in 1:n_files){
  gamma_tab <- data_out[[i]]
  gamma_tab[is.na(gamma_tab)] <- 0
  gamma_tab$island_name <- as.character(gamma_tab$island_name)
  class(gamma_tab) <- ("data.frame")
  # calculate total number of individuals per island
  gamma_div <- calc_biodiv(gamma_tab[,-1],
                           groups = gamma_tab$island_name,
                           index = "N")
  gamma_div$island_name <- gamma_div$group
  gamma_div <- subset(gamma_div, select=-c(effort,group,index))

  gamma_tab_N <- merge(gamma_tab,gamma_div, by = "island_name")

  # remove islands with less than 10 individuals

  gamma_tab_N10 <- subset(gamma_tab_N, value>10)

  gamma_tab_N10$value <- NULL
  gamma_tab_N10 <- as.data.frame(gamma_tab_N10)

  write.csv(gamma_tab_N10, paste(data_path1,"/", study_ids[i], ".csv",sep = ""), row.names = F)

  # removes islands with less than 20 individuals

  gamma_tab_N20 <- subset(gamma_tab_N, value>20)

  gamma_tab_N20$value <- NULL
  gamma_tab_N20 <- as.data.frame(gamma_tab_N20)

  write.csv(gamma_tab_N20, paste(data_path2,"/", study_ids[i], ".csv",sep = ""), row.names = F)

  # estimate reference n for rarefaction and extrapolations for islands with N >10
  gamma_tab_N10 <- as.data.frame(gamma_tab_N10)
  n_island_names1 <- rowSums(gamma_tab_N10[,-1])
  r <- 2
  n_ref1<- round(min(r * n_island_names1[n_island_names1 > 0]))


  gamma_div_N10 <- calc_biodiv(gamma_tab_N10[,-1],
                               groups = gamma_tab_N10$island_name,
                               index = c("S_n","S_asymp","S_PIE","S","N"),
                               effort = n_ref1,
                               extrapolate = T,
                               return_NA = F)

  gamma_div_N10 <- subset(gamma_div_N10, select=-c(effort))


  # estimate reference n for rarefaction and extrapolations for islands with N >20
  n_island_names2 <- rowSums(gamma_tab_N20[,-1])
  r <- 2
  n_ref2<- round(min(r * n_island_names2[n_island_names2 > 0]))


  gamma_div_N20 <- calc_biodiv(gamma_tab_N20[,-1],
                               groups = gamma_tab_N20$island_name,
                               index = c("S_n","S_asymp","S_PIE","S","N"),
                               effort = n_ref2,
                               extrapolate = T,
                               return_NA = F)

  gamma_div_N20 <- subset(gamma_div_N20, select=-c(effort))

  # convert data to wide format
  gamma_div_wide_N10 <- spread(gamma_div_N10, index, value)

  gamma_div_wide_N20 <- spread(gamma_div_N20, index, value)

  # calculate Spie using vegan

  group <- gamma_tab_N10$island_name
  group <- as.data.frame(group)

  ENS <- diversity(gamma_tab_N10[,-1], index = "invsimpson")
  ENS <- as.data.frame(ENS)

  ens_data <- merge(group, ENS, by = "row.names")
  ens_data <- ens_data[,-1]

  gamma_div_N10 <- merge(gamma_div_wide_N10, ens_data, by = "group")

  #########################################################################
  group <- gamma_tab_N20$island_name
  group <- as.data.frame(group)

  ENS <- diversity(gamma_tab_N20[,-1], index = "invsimpson")
  ENS <- as.data.frame(ENS)
  ENS$group <- gamma_tab_N20$island_name

  #ens_data1 <- merge(group, ENS1, by = "row.names")
  #ens_data1 <- ens_data1[,-1]

  gamma_div_N20 <- merge(gamma_div_wide_N20, ENS, by = "group")


  # remove islands with no diversity
  gamma_div_N10 <- subset(gamma_div_N10, S > 0)

  gamma_div_N20 <- subset(gamma_div_N20, S > 0)

  # Covert data to long format
  gamma_div_N10<- gather(gamma_div_N10, index, value, S_asymp, S_PIE, S_n, N, ENS)
  gamma_div_N20<- gather(gamma_div_N20, index, value, S_asymp, S_PIE, S_n, N, ENS)

  # Add study iDs

  gamma_div_N10$Study <- study_ids[i]
  gamma_div_N20$Study <- study_ids[i]

  div_out_N10[[i]] <- bind_rows(gamma_div_N10)
  div_out_N20[[i]] <- bind_rows(gamma_div_N20)

}


div_out_N10 <- bind_rows(div_out_N10)
div_out_N10 = div_out_N10 %>% as_tibble()

div_out_N20 <- bind_rows(div_out_N20)
div_out_N20 = div_out_N20 %>% as_tibble()


# load env file
env_file <- read_csv("./data/env_file_UTF-8.csv")
env_file$unique_id <- paste(env_file$Study,env_file$island_code, sep ="_")

# create unique ids to merge with env_file

div_out_N10$unique_id <- paste(div_out_N10$Study,div_out_N10$group,sep ="_")
div_out_N20$unique_id <- paste(div_out_N20$Study,div_out_N20$group,sep ="_")

# Merge diversity data of islands with > 10 individuals with the env file

all_ISAR_N10 <-  left_join(div_out_N10, env_file,
                           by = "unique_id")


# check to see whether some studies have less than 3 islands

all_ISAR_N10 %>%
  filter(index=='S_n') %>%
  group_by(Study.x) %>%
  tally() %>%
  filter(n<3)


# Merge diversity data of islands with > 20 individuals with the env file

all_ISAR_N20 <-  left_join(div_out_N20, env_file,
                           by = "unique_id")

# check to see whether some studies have less than 3 islands

all_ISAR_N20 %>%
  filter(index=='S_n') %>%
  group_by(Study.x) %>%
  tally() %>%
  filter(n<3)

# save data with all diversity metrics for :

# 1) data with islands where N>10
# 2) data with islands where N>20

saveRDS(all_ISAR_N10, file = "./data/diversity_global_synthesis_N10.Rds")
saveRDS(all_ISAR_N20, file = "./data/diversity_global_synthesis_N20.Rds")


