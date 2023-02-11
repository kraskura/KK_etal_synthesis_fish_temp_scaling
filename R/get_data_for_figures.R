
# Sys.Date()
date.imports<-"2023-01-23"
getModelParameters<-"NO" 
# ^ this is computationally heavy task, to get saved data: "NO" or 2
# to get compute now: "YES" or 1 

setwd("/Users/kristakraskura/Github_repositories/KK_etal_synthesis_fish_temp_scaling/")
source("./R/phylo_mixed model.R")

# 2. DATA update with mass specific values  --------
data.list<-get_scaling_data_temp(data.amr = "/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Data/MR-fish-metadata-data/Fish_AMR_temp_dataset_mar2022.csv",
                                 data.rmr = "/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Data/MR-fish-metadata-data/Fish_RMR_temp_dataset_mar2022.csv",
                                 ecology.data = "/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Data/MR-fish-metadata-data/Kraskura_species_ecologies_mar2022.csv", 
                                 onlyTop.above = TRUE,
                                 calc_mass_specific = TRUE,
                                 exp_rmr = round(fixef(rmr_mod_ER)[2],3),
                                 exp_amr = round(AMR.slopes$lnBWg.trend[3],3), # at 20ºC
                                 exp_as = round(fixef(as_mod_ER)[2],3),
                                 exp_rmr_warm = round(fixef(rmr_mod_W)[2],3),
                                 exp_amr_warm = round(fixef(amr_mod_W)[2],3),
                                 exp_as_warm = round(fixef(as_mod_W)[2],3))

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

k<-(8.62*10^(-5)) # Boltzmann's constant
E<-0.63 # activation energy MTE

data.amr$test_category2<-"acclim"
data.amr$test_category2[data.amr$test_category=="acute"] <- "acute"
data.amr$test_category3<-"ecol_relev"
data.amr$test_category3[data.amr$test_category=="acute"] <- "warm"
data.amr$test_category3[data.amr$test_category=="acclim"] <- "warm"
data.rmr$test_category2<-"acclim"
data.rmr$test_category2[data.rmr$test_category=="acute"] <- "acute"
data.rmr$test_category3<-"ecol_relev"
data.rmr$test_category3[data.rmr$test_category=="acute"] <- "warm"
data.rmr$test_category3[data.rmr$test_category=="acclim"] <- "warm"
data.fas$test_category2<-"acclim"
data.fas$test_category2[data.fas$test_category=="acute"] <- "acute"
data.fas$test_category3<-"ecol_relev"
data.fas$test_category3[data.fas$test_category=="acute"] <- "warm"
data.fas$test_category3[data.fas$test_category=="acclim"] <- "warm"
data.as$test_category2<-"acclim"
data.as$test_category2[data.as$test_category=="acute"] <- "acute"
data.as$test_category3<-"ecol_relev"
data.as$test_category3[data.as$test_category=="acute"] <- "warm"
data.as$test_category3[data.as$test_category=="acclim"] <- "warm"

# these datasets are > Topt temperatures
data.amr.test<-rbind(data.amrAC, data.amrAM)
data.rmr.test<-rbind(data.rmrAC, data.rmrAM)
data.fas.test<-data.amr.test[c(!is.na(data.amr.test$FAS) & is.finite(data.amr.test$FAS)) , ]
data.as.test<-data.amr.test[c(!is.na(data.amr.test$lnAS) & is.finite(data.amr.test$lnAS)) , ]

data.amrER$tempTestK1000_inC<-((1000/data.amrER$tempTestK1000))-273.15



# Ecology groups -----
ecology_data.FAS <- read.csv(file = "./Data_exports/Ecologies/ecologies_FAS.csv")
ecology_data.RMR <- read.csv(file = "./Data_exports/Ecologies/ecologies_RMR.csv")
ecology_data.AMR <- read.csv(file = "./Data_exports/Ecologies/ecologies_MMR.csv")

# only those that have warm 
ecology_data.AMRd<-ecology_data.AMR[ ave(1:nrow(ecology_data.AMR), ecology_data.AMR$ecology_subgroup, FUN=length) > 1 , ]
ecology_data.FASd<-ecology_data.FAS[ ave(1:nrow(ecology_data.FAS), ecology_data.FAS$ecology_subgroup, FUN=length) > 1 , ]
ecology_data.RMRd<-ecology_data.RMR[ ave(1:nrow(ecology_data.RMR), ecology_data.RMR$ecology_subgroup, FUN=length) > 1 , ]
ecology_data.AMRd$ecol_temp_cat<-factor(ecology_data.AMRd$ecol_temp_cat)
ecology_data.RMRd$ecol_temp_cat<-factor(ecology_data.RMRd$ecol_temp_cat)
ecology_data.FASd$ecol_temp_cat<-factor(ecology_data.FASd$ecol_temp_cat)

order.ecology_data.FASd<- ecology_data.FASd 
order.ecology_data.FASd$ecol_temp_cat<- factor(order.ecology_data.FASd$ecol_temp_cat,
                                               levels = order.ecology_data.FASd$ecol_temp_cat[order(order.ecology_data.FASd$slope)])

# warm only:
ecology_data.RMRd_w<-ecology_data.RMRd[grepl("warm", as.character(ecology_data.RMRd$ecol_temp_cat)), ]
ecology_data.AMRd_w<-ecology_data.AMRd[grepl("warm", as.character(ecology_data.AMRd$ecol_temp_cat)), ]
ecology_data.FASd_w<-ecology_data.FASd[grepl("warm", as.character(ecology_data.FASd$ecol_temp_cat)), ]

# ecol relev only:
ecology_data.RMRd_er<-ecology_data.RMRd[grepl("ecol_relev", as.character(ecology_data.RMRd$ecol_temp_cat)), ]
ecology_data.AMRd_er<-ecology_data.AMRd[grepl("ecol_relev", as.character(ecology_data.AMRd$ecol_temp_cat)), ]
ecology_data.FASd_er<-ecology_data.FASd[grepl("ecol_relev", as.character(ecology_data.FASd$ecol_temp_cat)), ]

order.ecology_data.FASd.W <- order.ecology_data.FASd[grepl("warm", as.character(order.ecology_data.FASd$ecol_temp_cat)), ]
order.ecology_data.FASd.ER <- order.ecology_data.FASd[grepl("ecol_relev", as.character(order.ecology_data.FASd$ecol_temp_cat)), ]

ecology_data.RMRd_er$DIFF_mmr_rmr<-ecology_data.RMRd_er$slope - ecology_data.AMRd_er$slope
ecology_data.RMRd_er[ecology_data.RMRd_er$DIFF_mmr_rmr > 0,] # bMMR < bRMR
ecology_data.RMRd_er[ecology_data.RMRd_er$DIFF_mmr_rmr < 0,] # bMMR > bRMR

ecology_data.RMRd_w$DIFF_mmr_rmr<-ecology_data.RMRd_w$slope - ecology_data.AMRd_w$slope
ecology_data.RMRd_w[ecology_data.RMRd_w$DIFF_mmr_rmr > 0,] # bMMR < bRMR
ecology_data.RMRd_w[ecology_data.RMRd_w$DIFF_mmr_rmr < 0,] # bMMR > bRMR

ecology_data.AMRd_w[ecology_data.AMRd_w$Ecol == "salintyComb",]
ecology_data.RMRd_w[ecology_data.RMRd_w$Ecol == "salintyComb",]
ecology_data.RMRd_er[ecology_data.RMRd_er$Ecol == "salintyComb",]
ecology_data.AMRd_er[ecology_data.AMRd_er$Ecol == "salintyComb",]

ecology_data.RMRd.g<-ecology_data.RMRd[c(ecology_data.RMRd$scalingGood=="fusiform" | 
                                           ecology_data.RMRd$ecology_subgroup=="Subtropical"|
                                           ecology_data.RMRd$ecology_subgroup=="Temperate"|
                                           ecology_data.RMRd$ecology_subgroup=="Tropical"|
                                           ecology_data.RMRd$ecology_subgroup=="benthopelagic"|
                                           ecology_data.RMRd$ecology_subgroup=="demersal"|
                                           ecology_data.RMRd$ecology_subgroup=="reef-associated"|
                                           ecology_data.RMRd$ecology_subgroup=="Marine"|
                                           ecology_data.RMRd$ecology_subgroup=="Marine; freshwater; brackish"),]

ecology_data.AMRd.g<-ecology_data.AMRd[c(ecology_data.AMRd$ecology_subgroup=="fusiform" | 
                                           ecology_data.AMRd$ecology_subgroup=="Subtropical"|
                                           ecology_data.AMRd$ecology_subgroup=="Temperate"|
                                           ecology_data.AMRd$ecology_subgroup=="Tropical"|
                                           ecology_data.AMRd$ecology_subgroup=="benthopelagic"|
                                           ecology_data.AMRd$ecology_subgroup=="demersal"|
                                           ecology_data.AMRd$ecology_subgroup=="reef-associated"|
                                           ecology_data.AMRd$ecology_subgroup=="Marine"|
                                           ecology_data.AMRd$ecology_subgroup=="Marine; freshwater; brackish"),]

ecology_data.FASd.g<-ecology_data.FASd[c(ecology_data.FASd$ecology_subgroup=="fusiform" | 
                                           ecology_data.FASd$ecology_subgroup=="Subtropical"|
                                           ecology_data.FASd$ecology_subgroup=="Temperate"|
                                           ecology_data.FASd$ecology_subgroup=="Tropical"|
                                           ecology_data.FASd$ecology_subgroup=="benthopelagic"|
                                           ecology_data.FASd$ecology_subgroup=="demersal"|
                                           ecology_data.FASd$ecology_subgroup=="reef-associated"|
                                           ecology_data.FASd$ecology_subgroup=="Marine"|
                                           ecology_data.FASd$ecology_subgroup=="Marine; freshwater; brackish"),]

# SPECIES -----
species.fas<-read.csv("./Data_exports/Species/species_FAS.csv")
species.rmr<-read.csv("./Data_exports/Species/species_RMR.csv")
species.amr<-read.csv("./Data_exports/Species/species_MMR.csv")

species.fas.w<-species.fas[species.fas$test_category3=="warm", ]
species.fas.er<-species.fas[species.fas$test_category3=="ecol_relev", ]
species.rmr.w<-species.rmr[species.rmr$test_category3=="warm", ]
species.rmr.er<-species.rmr[species.rmr$test_category3=="ecol_relev", ]
species.amr.w<-species.amr[species.amr$test_category3=="warm", ]
species.amr.er<-species.amr[species.amr$test_category3=="ecol_relev", ]

