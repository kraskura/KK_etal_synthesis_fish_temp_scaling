

# Import libraries ----------
library(tidyverse)
library(dplyr)
library(reshape2)


# source the data --------
source("/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Codes/Metadata-temperature/JEB/get_scaling_data_temp.R")
source("/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Codes/MR-fish-metadata/deltaIC_BICdelta_ggformat.R")

data.list<-get_scaling_data_temp(data.amr = "/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Data/MR-fish-metadata-data/Fish_AMR_temp_dataset_mar2022.csv",
                                 data.rmr = "/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Data/MR-fish-metadata-data/Fish_RMR_temp_dataset_mar2022.csv",
                                 ecology.data = "/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Data/MR-fish-metadata-data/Kraskura_species_ecologies_mar2022.csv", 
                                 onlyTop.above = TRUE, 
                                 exp_rmr = 0.81, exp_amr_warm = 0.81, 
                                 exp_rmr_warm = 0.84, exp_amr = 0.82, 
                                 exp_as = 0.82, exp_as_warm = 0.78)

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


# Data where FAS is greater than 12. why?
data.fas12<-data.fas[(data.fas$FAS>12),]
# write.csv(file= "/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Data/MR-fish-metadata-data/dataFAS_largerThan12_mar22-2022_ALLtempCATEG.csv", data.fas12, row.names=FALSE)

# studies that has FAS > 12
unique(data.fas12$study_ID)
mean(data.fas$FAS[data.fas$study_ID==4])
mean(data.fas$FAS[data.fas$study_ID==6])
mean(data.fas$FAS[data.fas$study_ID==15])
mean(data.fas$FAS[data.fas$study_ID==20])
mean(data.fas$FAS[data.fas$study_ID==61])
mean(data.fas$FAS[data.fas$study_ID==410])
mean(data.fas$FAS[data.fas$study_ID==405])

length(data.fas$FAS[data.fas$study_ID==4])
length(data.fas$FAS[data.fas$study_ID==6])
length(data.fas$FAS[data.fas$study_ID==15])
length(data.fas$FAS[data.fas$study_ID==20])
length(data.fas$FAS[data.fas$study_ID==61])
length(data.fas$FAS[data.fas$study_ID==410])
length(data.fas$FAS[data.fas$study_ID==405])

k<-(8.62*10^(-5)) # Boltzmann's constant - units eV/K
# Boltz<-1.38 * 10^23; units J/K
E<-0.63 # activation energy MTE
# cols <- c("#d55e00", "#0072b2", "#009e73")

cols.as<-c("#265F73", "#00929A", "#00C5A3")
cols.rmr<-c("#C70039", "#FF6D7C", "#FFA3AC")
cols.amr<-c("#00749F","#00A8D6", "#9CE9FF")
cols.fas<-c("#C94F00", "#F37121", "#FFD784")

# the warm ones 
data.amr.test<-rbind(data.amrAC, data.amrAM)
data.rmr.test<-rbind(data.rmrAC, data.rmrAM)
data.fas.test<-data.amr.test[c(!is.na(data.amr.test$FAS) & is.finite(data.amr.test$FAS)) , ]
data.as.test<-data.amr.test[c(!is.na(data.amr.test$lnAS) & is.finite(data.amr.test$lnAS)) , ]

nrow(data.amr.test)+nrow(data.amrER)
nrow(data.rmr.test)+nrow(data.rmrER)

data.fas$test_category3<-"ecol_relev"
data.fas$test_category3[data.fas$test_category=="acute"] <- "warm"
data.fas$test_category3[data.fas$test_category=="acclim"] <- "warm"

data.as$test_category3<-"ecol_relev"
data.as$test_category3[data.as$test_category=="acute"] <- "warm"
data.as$test_category3[data.as$test_category=="acclim"] <- "warm"

data.amr$test_category3<-"ecol_relev"
data.amr$test_category3[data.amr$test_category=="acute"] <- "warm"
data.amr$test_category3[data.amr$test_category=="acclim"] <- "warm"

data.rmr$test_category3<-"ecol_relev"
data.rmr$test_category3[data.rmr$test_category=="acute"] <- "warm"
data.rmr$test_category3[data.rmr$test_category=="acclim"] <- "warm"


# check for any NAs
table(is.na(dataMR))
levels(dataMR$test_category)
levels(dataMR$species)
unique(dataMR$species)
# unique(data.rmr$species)
# unique(data.amr$species)
levels(factor(dataMR$MR_type2))
length(levels(dataMR$species)) # 94 species aug 2022
nrow(dataMR) # 10287 (when AMR, RMR, both together) mar 2022




# 1. test categories -------
summaryFAS <-  data.fas %>%
  dplyr:::group_by(test_category3) %>%
  summarise(size_mean = mean(BW_g), size_min = min(BW_g), size_max = max(BW_g), 
            temp_mean = mean(tempTest), temp_min = min(tempTest), temp_max = max(tempTest),
            n_species = length(unique(species)), 
            n = length(BW_g),
            n_trial = length(unique(trial)), 
            n_studies = length(unique(study_ID)))%>% 
  as.data.frame()
summaryFAS$temp_sum<-paste(round(summaryFAS$temp_mean, 2), " (", round(summaryFAS$temp_min, 2), ",", round(summaryFAS$temp_max,2), ")", sep="")
summaryFAS$size_sum<-paste(round(summaryFAS$size_mean, 2), " (", round(summaryFAS$size_min, 2), ",", round(summaryFAS$size_max,2), ")", sep="")
summaryFAS$MR<-"FAS"

summaryAS <-  data.as %>%
  dplyr:::group_by(test_category3) %>%
  summarise(size_mean = mean(BW_g), size_min = min(BW_g), size_max = max(BW_g), 
            temp_mean = mean(tempTest), temp_min = min(tempTest), temp_max = max(tempTest),
            n_species = length(unique(species)), 
            n = length(BW_g),
            n_trial = length(unique(trial)), 
            n_studies = length(unique(study_ID)))%>% 
  as.data.frame()
summaryAS$temp_sum<-paste(round(summaryAS$temp_mean, 2), " (", round(summaryAS$temp_min, 2), ",", round(summaryAS$temp_max,2), ")", sep="")
summaryAS$size_sum<-paste(round(summaryAS$size_mean, 2), " (", round(summaryAS$size_min, 2), ",", round(summaryAS$size_max,2), ")", sep="")
summaryAS$MR<-"AS"

summaryAMR <- data.amr %>%
  dplyr:::group_by(test_category3) %>%
  summarise(size_mean = mean(BW_g), size_min = min(BW_g), size_max = max(BW_g), 
            temp_mean = mean(tempTest), temp_min = min(tempTest), temp_max = max(tempTest),
            n_species = length(unique(species)), 
            n = length(BW_g),
            n_trial = length(unique(trial)), 
            n_studies = length(unique(study_ID)))%>% 
  as.data.frame()
summaryAMR$temp_sum<-paste(round(summaryAMR$temp_mean, 2), " (", round(summaryAMR$temp_min, 2), ",", round(summaryAMR$temp_max,2), ")", sep="")
summaryAMR$size_sum<-paste(round(summaryAMR$size_mean, 2), " (", round(summaryAMR$size_min, 2), ",", round(summaryAMR$size_max,2), ")", sep="")
summaryAMR$MR<-"AMR"

summaryRMR <- data.rmr %>%
  dplyr:::group_by(test_category3) %>%
  summarise(size_mean = mean(BW_g), size_min = min(BW_g), size_max = max(BW_g), 
            temp_mean = mean(tempTest), temp_min = min(tempTest), temp_max = max(tempTest),
            n_species = length(unique(species)), 
            n = length(BW_g),
            n_trial = length(unique(trial)), 
            n_studies = length(unique(study_ID)))%>% 
  as.data.frame()
summaryRMR$temp_sum<-paste(round(summaryRMR$temp_mean, 2), " (", round(summaryRMR$temp_min, 2), ",", round(summaryRMR$temp_max,2), ")", sep="")
summaryRMR$size_sum<-paste(round(summaryRMR$size_mean, 2), " (", round(summaryRMR$size_min, 2), ",", round(summaryRMR$size_max,2), ")", sep="")
summaryRMR$MR<-"RMR"

sum_table1<-rbind(summaryFAS, summaryAS, summaryAMR, summaryRMR)

write.csv(file = "/Users/kristakraskura/Desktop/BOX/UCSB/Research/Metabolic_scaling/ms-AMR-RMR-Temperature/thesis/data files/table1_summary.csv", sum_table1, row.names=FALSE)

# 2. overall, all data -----
length(levels(dataMR$trial))
length(levels(dataMR$study_ID))
length(levels(dataMR$species))
nrow(dataMR)

length(levels(data.as$trial))
length(levels(data.as$study_ID))
length(levels(data.as$species))
nrow(data.as)

# how many fish total? 
length(unique(dataMR$fish_ID))
nrow(dataMR[!duplicated(dataMR$fish_ID),])
nrow(dataMR[c(!duplicated(dataMR$fish_ID)) & dataMR$BW_g<1000,])


data.fas %>% 
  group_by(test_category3) %>% 
  summarise(min_fas = min(FAS), max_fas = max(FAS), mean_fas = mean(FAS), n = length(FAS), se = sd(FAS)/(sqrt(n)))


# 3. ecologies -------

data.fas %>% 
  group_by(DemersPelag, test_category3) %>% 
  summarise(min_fas = min(FAS), max_fas = max(FAS), mean_fas = mean(FAS), n = length(FAS))
data.fas[(data.fas$FAS>12),]

data.fas %>% 
  group_by(Climate, test_category3) %>% 
  summarise(min_fas = min(FAS), max_fas = max(FAS), mean_fas = mean(FAS), n = length(FAS))

data.fas %>% 
  group_by(BodyShapeI, test_category3) %>% 
  summarise(min_fas = min(FAS), max_fas = max(FAS), mean_fas = mean(FAS), n = length(FAS))

data.fas %>% 
  group_by(salintyComb, test_category3) %>% 
  summarise(min_fas = min(FAS), max_fas = max(FAS), mean_fas = mean(FAS), n = length(FAS))

# tropical
nrow(data.fas[(data.fas$Climate == "Tropical" & data.fas$FAS<3 & data.fas$test_category3 == "ecol_relev"),]) /
  nrow(data.fas[(data.fas$Climate == "Tropical" & data.fas$test_category3 == "ecol_relev"),])

# subtropical
nrow(data.fas[(data.fas$Climate == "Subtropical" & data.fas$FAS<3 & data.fas$test_category3 == "ecol_relev"),]) /
  nrow(data.fas[(data.fas$Climate == "Subtropical" & data.fas$test_category3 == "ecol_relev"),])

# temperates
nrow(data.fas[(data.fas$Climate == "Temperate" & data.fas$FAS<3 & data.fas$test_category3 == "ecol_relev"),]) /
  nrow(data.fas[(data.fas$Climate == "Temperate" & data.fas$test_category3 == "ecol_relev"),])

# polar
nrow(data.fas[(data.fas$Climate == "Polar" & data.fas$FAS<3 & data.fas$test_category3 == "ecol_relev"),]) /
nrow(data.fas[(data.fas$Climate == "Polar" & data.fas$test_category3 == "ecol_relev"),])


# how many fish < 1 kg?
nrow(data.fas[(data.fas$Climate == "Polar" & data.fas$test_category3 == "ecol_relev"),])
data.fas.test


# overall sums 
# summaryFAS <-  data.fas %>%
#   dplyr:::group_by(test_category3) %>%
#   summarise(size_mean = mean(BW_g), size_min = min(BW_g), size_max = max(BW_g), 
#             temp_mean = mean(tempTest), temp_min = min(tempTest), temp_max = max(tempTest),
#             n_species = length(unique(species)), 
#             n = length(BW_g),
#             n_trial = length(unique(trial)), 
#             n_studies = length(unique(study_ID)))%>% 
#   as.data.frame()
# summaryFAS$temp_sum<-paste(round(summaryFAS$temp_mean, 2), " (", round(summaryFAS$temp_min, 2), ",", round(summaryFAS$temp_max,2), ")", sep="")
# summaryFAS$size_sum<-paste(round(summaryFAS$size_mean, 2), " (", round(summaryFAS$size_min, 2), ",", round(summaryFAS$size_max,2), ")", sep="")
# summaryFAS$MR<-"FAS"
# 
# summaryAS <-  data.as %>%
#   dplyr:::group_by(test_category3) %>%
#   summarise(size_mean = mean(BW_g), size_min = min(BW_g), size_max = max(BW_g), 
#             temp_mean = mean(tempTest), temp_min = min(tempTest), temp_max = max(tempTest),
#             n_species = length(unique(species)), 
#             n = length(BW_g),
#             n_trial = length(unique(trial)), 
#             n_studies = length(unique(study_ID)))%>% 
#   as.data.frame()
# summaryAS$temp_sum<-paste(round(summaryAS$temp_mean, 2), " (", round(summaryAS$temp_min, 2), ",", round(summaryAS$temp_max,2), ")", sep="")
# summaryAS$size_sum<-paste(round(summaryAS$size_mean, 2), " (", round(summaryAS$size_min, 2), ",", round(summaryAS$size_max,2), ")", sep="")
# summaryAS$MR<-"AS"
# 
# summaryAMR <- data.amr %>%
#   dplyr:::group_by(test_category3) %>%
#   summarise(size_mean = mean(BW_g), size_min = min(BW_g), size_max = max(BW_g), 
#             temp_mean = mean(tempTest), temp_min = min(tempTest), temp_max = max(tempTest),
#             n_species = length(unique(species)), 
#             n = length(BW_g),
#             n_trial = length(unique(trial)), 
#             n_studies = length(unique(study_ID)))%>% 
#   as.data.frame()
# summaryAMR$temp_sum<-paste(round(summaryAMR$temp_mean, 2), " (", round(summaryAMR$temp_min, 2), ",", round(summaryAMR$temp_max,2), ")", sep="")
# summaryAMR$size_sum<-paste(round(summaryAMR$size_mean, 2), " (", round(summaryAMR$size_min, 2), ",", round(summaryAMR$size_max,2), ")", sep="")
# summaryAMR$MR<-"AMR"
# 
summaryRMRecol.sal <- data.rmr %>%
  dplyr:::group_by(test_category3, salintyComb) %>%
  summarise(size_mean = mean(BW_g), size_min = min(BW_g), size_max = max(BW_g),
            temp_mean = mean(tempTest), temp_min = min(tempTest), temp_max = max(tempTest),
            n_species = length(unique(species)),
            n = length(BW_g),
            n_trial = length(unique(trial)),
            n_studies = length(unique(study_ID)))%>%
  as.data.frame()

summaryRMRecol.morph <- data.rmr %>%
  dplyr:::group_by(test_category3, BodyShapeI) %>%
  summarise(size_mean = mean(BW_g), size_min = min(BW_g), size_max = max(BW_g),
            temp_mean = mean(tempTest), temp_min = min(tempTest), temp_max = max(tempTest),
            n_species = length(unique(species)),
            n = length(BW_g),
            n_trial = length(unique(trial)),
            n_studies = length(unique(study_ID)))%>%
  as.data.frame()
# summaryRMR$temp_sum<-paste(round(summaryRMR$temp_mean, 2), " (", round(summaryRMR$temp_min, 2), ",", round(summaryRMR$temp_max,2), ")", sep="")
# summaryRMR$size_sum<-paste(round(summaryRMR$size_mean, 2), " (", round(summaryRMR$size_min, 2), ",", round(summaryRMR$size_max,2), ")", sep="")
# summaryRMR$MR<-"RMR"
# 
# rbind(summaryFAS, summaryAS, summaryAMR, summaryRMR)


# summaries
summary(data.rmrER$tempTest)
min(data.rmrER$BW_g)  
max(data.rmrER$BW_g)  
mean(data.rmrER$BW_g)

summary(data.amrER$tempTest)
min(data.amrER$BW_g)  
max(data.amrER$BW_g)  
mean(data.amrER$BW_g)

summary(data.fasER$tempTest)
min(data.fasER$BW_g)  
max(data.fasER$BW_g)  
mean(data.fasER$BW_g)

# warm
summary(data.rmr.test$tempTest)
min(data.rmr.test$BW_g)  
max(data.rmr.test$BW_g)  
mean(data.rmr.test$BW_g)

summary(data.amr.test$tempTest)
min(data.amr.test$BW_g)  
max(data.amr.test$BW_g)  
mean(data.amr.test$BW_g)

summary(data.fas.test$tempTest)
min(data.fas.test$BW_g)  
max(data.fas.test$BW_g)  
mean(data.fas.test$BW_g)

# 4. species ecology and morphology categorisation --------
data.amr$categ_species<-paste(data.amr$test_category3, data.amr$species, sep="-")
unique(data.amr$categ_species)
data.rmr$categ_species<-paste(data.rmr$test_category3, data.rmr$species, sep="-")
unique(data.rmr$categ_species)

data.amr.sum.ecol<-data.amr[!duplicated(data.amr$categ_species),c("species","test_category3", "salintyComb", "Climate", "BodyShapeI", "DemersPelag")]
data.amr.sum.ecol$MR<-"MMR"
unique(data.amr.sum.ecol$species)
# names(data.amr.sum.ecol)
data.rmr.sum.ecol<-data.rmr[!duplicated(data.rmr$categ_species),c("species","test_category3", "salintyComb", "Climate", "BodyShapeI", "DemersPelag")]
data.rmr.sum.ecol$MR<-"RMR"
unique(data.rmr.sum.ecol$species)

data.sum.ecol.gr<-rbind(data.amr.sum.ecol, data.rmr.sum.ecol)
data.sum.ecol.gr$categ_mr<-paste(data.sum.ecol.gr$test_category3, data.sum.ecol.gr$MR, sep="-")
data.sum.ecol.gr2<-as.data.frame(dcast(data.sum.ecol.gr, species ~ categ_mr))
colnames(data.sum.ecol.gr2)<-c("species","ecol_relev_MMR", "ecol_relev_RMR", "warm_MMR", "warm_RMR")
nrow(data.sum.ecol.gr2[!is.na(data.sum.ecol.gr2$ecol_relev_MMR),])
nrow(data.sum.ecol.gr2[!is.na(data.sum.ecol.gr2$ecol_relev_RMR),])
data.sum.ecol.gr2$ecol_relev_MMR[!is.na(data.sum.ecol.gr2$ecol_relev_MMR)]<-"Y"
data.sum.ecol.gr2$ecol_relev_RMR[!is.na(data.sum.ecol.gr2$ecol_relev_RMR)]<-"Y"
data.sum.ecol.gr2$warm_MMR[!is.na(data.sum.ecol.gr2$warm_MMR)]<-"Y"
data.sum.ecol.gr2$warm_RMR[!is.na(data.sum.ecol.gr2$warm_RMR)]<-"Y"
data.sum.ecol.gr2<-left_join(data.sum.ecol.gr2, dataMR[!duplicated(dataMR$species), c("Climate", "salintyComb", "BodyShapeI", "DemersPelag", "species")], by = "species")

write.csv(file = "/Users/kristakraskura/Desktop/BOX/UCSB/Research/Metabolic_scaling/ms-AMR-RMR-Temperature/thesis/data files/table2Sup_Ecolgroups_summary.csv", data.sum.ecol.gr2, row.names=FALSE)









