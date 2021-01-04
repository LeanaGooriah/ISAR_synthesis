#####################################################################################
# Author : Leana Gooriah
# Email : leana.gooriah@gmail.com
#####################################################################################
# Description : This script can be used to calculate the following diversity indices:
# Stotal, S_n, S_PIE (following 'mobr' terminology)

#####################################################################################
# 1. load data
# 2. calculate Stotal & S_PIE 
# 3. calculate reference n 
# 4. calculate S_n 

# LOAD R PACKAGES
rm(list = ls())
require(dplyr)
require(tidyr)
library(devtools)
library(vegan)
library(tidyr)
install_github('MoBiodiv/mobr') # to install most recent version of mobr
require(mobr) # version 1.0

############################
# Set path and directories #
############################
work_dir <- getwd() # first set working directory in Menu/Session/Set working directory/to Project 

data_path <- paste(work_dir,"/ISAR_datasets",sep="")

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


div_out <- list()

for (i in 1:n_files){
  
  gamma_tab <- data_out[[i]]
  gamma_tab[is.na(gamma_tab)] <- 0
  gamma_tab$island_name <- as.character(gamma_tab$island_name)
  class(gamma_tab) <- ("data.frame")
  # estimate reference n for rarefaction and extrapolations
  n_island_names <- rowSums(gamma_tab[,-1])
  r <- 2
  n_ref <- round(min(r * n_island_names[n_island_names > 0]))
  
  
  gamma_div <- calc_biodiv(gamma_tab[,-1], 
                           groups = gamma_tab$island_name,
                           index = c("S_n","S_asymp","S_PIE","S"),
                           effort = n_ref,
                           extrapolate = T,
                           return_NA = F)
  
  gamma_div <- subset(gamma_div, select=-c(effort))
  gamma_div_wide <- spread(gamma_div, index, value)

  
  # calculate Spie using vegan 
  
  group <- gamma_tab$island_name
  group <- as.data.frame(group)
  
  ENS <- diversity(gamma_tab[,-1], index = "invsimpson")
  ENS <- as.data.frame(ENS)
  
  ens_data <- merge(group, ENS, by = "row.names")
  ens_data <- ens_data[,-1]
  
  gamma_div <- merge(gamma_div_wide, ens_data, by = "group")
  # Add study iDs  
  
  gamma_div <- subset(gamma_div, S > 0)
  gamma_div<- gather(gamma_div, index, value, S_asymp, S_PIE, S_n, ENS)
  
  gamma_div$Study <- study_ids[i]
  
  div_out[[i]] <- bind_rows(gamma_div)
  
}

div_out <- bind_rows(div_out)
div_out = div_out %>% as_tibble()


library(readr)
env_file <- read_csv("ISAR_meta_analysis/env_file.csv")
env_file$unique_id <- paste(env_file$Study,env_file$island_code, sep ="_")


div_out$unique_id <- paste(div_out$Study,div_out$group,sep ="_")

  
# Join env file with diversity file  
  
  all_ISAR <-  left_join(div_out, env_file,
             by = "unique_id")
  
  save(all_ISAR, file = (paste(data_path, "/diversity_global_synthesis_with_area.Rdata",sep="")))
