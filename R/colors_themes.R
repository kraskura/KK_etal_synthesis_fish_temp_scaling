
# colors 
# if neccessary : 
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("ggtree")

# install.packages("lme4", type = "source")
# install.packages("evolvability", type = "source")
library(ggh4x)
library(tidyverse)
library(dplyr)

library(colorBlindness)
library(patchwork)
library(weathermetrics)

library(lme4)
library(emmeans)
library(car)

library(evolvability) # Almer function
library(ape)
library(rotl)
library(ggtree)
library(Matrix)

# library(kableExtra)

library(pryr)
library(viridisLite)
library(TDbook)
library(ggimage)

library(ggformat2) # from github curstom built to format figures. 
library(weathermetrics)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(forcats)
library(patchwork)
library(reshape2)
library(emmeans)
library(systemfonts)
library(flextable)
library(here)


# displayAvailablePalette(color="white")
# colorBlindness::SteppedSequential5Steps

cols.amr<-c(colorBlindness::SteppedSequential5Steps[c(21, 22, 23, 24,  25)]) 
cols.rmr<-c(colorBlindness::SteppedSequential5Steps[c(11, 12, 13,14,  15)]) # "#6B990F" "#A3CC51" "#E5FFB2"
cols.fas<-c(colorBlindness::SteppedSequential5Steps[c(1,2, 3,4, 5)])  # "#990F0F" "#CC5151" "#FFB2B2"
cols.as<-c(colorBlindness::SteppedSequential5Steps[c(16,17,  18, 19, 20)]) 


cols<-c(cols.amr[1],cols.rmr[1],cols.amr[3],cols.rmr[3], cols.amr[5],cols.rmr[5], cols.as[2], cols.fas[2])
        # AMR -rmr- AMR dark - rmr dark - light - as - fas

# DATA SETS -------------------
message("import datasets")

data.list<-get_data_temp(data.amr = here("Data", "Fish_AMR_temp_dataset_mar2022.csv"),
                                 data.rmr = here("Data","Fish_RMR_temp_dataset_mar2022.csv"),
                                 ecology.data = here("Data", "Kraskura_species_ecologies_mar2022.csv"),
                                 onlyTop.above = TRUE,
                                 save.FishBase.species.data = F,
                                 calc_mass_specific = FALSE)

# data.list.NO.SALMON<-get_data_temp(data.amr = here("Data", "Fish_AMR_temp_dataset_mar2022.csv"),
#                                  data.rmr = here("Data","Fish_RMR_temp_dataset_mar2022.csv"),
#                                  ecology.data = here("Data", "Kraskura_species_ecologies_mar2022.csv"),
#                                  onlyTop.above = TRUE, save.FishBase.species.data = F,
#                                  calc_mass_specific = FALSE)

data.amrAC<-data.frame(data.list[1])
data.rmrAC<-data.frame(data.list[2])
data.amrAM<-data.frame(data.list[3])
data.rmrAM<-data.frame(data.list[4])
data.amrER<-data.frame(data.list[5])
data.rmrER<-data.frame(data.list[6])

data.asAC<-data.frame(data.list[7])
data.fasAC<-data.frame(data.list[8])
data.asAM<-data.frame(data.list[9])
data.fasAM<-data.frame(data.list[10])
data.asER<-data.frame(data.list[11])
data.fasER<-data.frame(data.list[12])

data.amr<-data.frame(data.list[13])
data.rmr<-data.frame(data.list[14])
dataMR<-data.frame(data.list[15])
data.as<-data.frame(data.list[16])
data.fas<-data.frame(data.list[17])

# the warm ones
data.amr.test<-rbind(data.amrAC, data.amrAM)
data.rmr.test<-rbind(data.rmrAC, data.rmrAM)
data.fas.test<-rbind(data.fasAC, data.fasAM)
data.as.test<-rbind(data.asAC, data.asAM)
data.fas.test<-data.fas.test[c(!is.na(data.fas.test$FAS) & is.finite(data.fas.test$FAS)) , ]
data.as.test<-data.as.test[c(!is.na(data.as.test$lnAS) & is.finite(data.as.test$lnAS)) , ]

# Get model data set specific phylogentics model matrixes -------
data.rmrER<-droplevels(data.rmrER)
data.rmr.test<-droplevels(data.rmr.test)
data.amrER<-droplevels(data.amrER)
data.amr.test<-droplevels(data.amr.test)
data.asER<-droplevels(data.asER)
data.as.test<-droplevels(data.as.test)
data.fasER<-droplevels(data.fasER)
data.fas.test<-droplevels(data.fas.test)

# test categories
data.fas$test_category3 <- "warm"
data.fas[data.fas$test_category == "ecol_relev", "test_category3"] <- "optimal"
data.as$test_category3 <- "warm"
data.as[data.as$test_category == "ecol_relev", "test_category3"] <- "optimal"
data.amr$test_category3 <- "warm"
data.amr[data.amr$test_category == "ecol_relev", "test_category3"] <- "optimal"
data.rmr$test_category3 <- "warm"
data.rmr[data.rmr$test_category == "ecol_relev", "test_category3"] <- "optimal"

