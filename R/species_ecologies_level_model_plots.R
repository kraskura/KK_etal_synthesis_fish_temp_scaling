

options(dplyr.summarise.inform = FALSE)
# *****************************************************************************
#  not needed for final version
# k<-(8.62*10^(-5)) # Boltzmann's constant
# E<-0.63 # activation energy MTE

# *******************************************************************************

# for emmeans library
# Ecologically Relevant conditions: 
emm_options(pbkrtest.limit = 7000)
emm_options(lmerTest.limit = 7000)

# *******************************************************************************
# *******************************************************************************
# GLOBAL DATA (with mass specific values) ---------

# 1) for phylo mixed models
source(here("R", "phylo_mixed_model.R")) # get final model outputs
source(here("R/colors_themes.R"))
# 2) for non-phylo mixed models 
# source("./R/nonPhylo_mixed_models.R")

# *********************************************************************************
# **********************************************************************************
# ECOLOGY - SPECIFIC (estimate scaling) ----
# 
summarise_ecologiesLM<-function(scalingGood=FALSE,
                                NAMEecol,
                                lmdata.AMRW,
                                summarylm.AMR,
                                lmdata.RMRW,
                                summarylm.RMR,
                                lmdata.FASW,
                                summarylm.FAS){
  
  lmdata.AMR$ecol_temp_cat<-paste(lmdata.AMR[,1], lmdata.AMR$test_category3, sep="_")
  summarylm.AMR$ecol_temp_cat<-paste(summarylm.AMR[,1], summarylm.AMR$test_category3, sep="_")
  
  lmdata.RMR$ecol_temp_cat<-paste(lmdata.RMR[,1], lmdata.RMR$test_category3, sep="_")
  summarylm.RMR$ecol_temp_cat<-paste(summarylm.RMR[,1], summarylm.RMR$test_category3, sep="_")
  
  lmdata.FAS$ecol_temp_cat<-paste(lmdata.FAS[,1], lmdata.FAS$test_category3, sep="_")
  summarylm.FAS$ecol_temp_cat<-paste(summarylm.FAS[,1], summarylm.FAS$test_category3, sep="_")
  
  
  lmdata.AMRW<-dcast(melt(lmdata.AMR[, c("ecol_temp_cat", "estimate", "term")], id.vars=c("ecol_temp_cat", "estimate")),
                     ecol_temp_cat ~ value + variable, value.var=c("estimate"))
  lmdata.AMRWciL<-dcast(melt(lmdata.AMR[, c("ecol_temp_cat",  "term", "conf.low")], id.vars=c("ecol_temp_cat",  "conf.low")),
                       ecol_temp_cat ~ value + variable, value.var=c("conf.low"))
  lmdata.AMRWciH<-dcast(melt(lmdata.AMR[, c("ecol_temp_cat",  "term", "conf.high")], id.vars=c("ecol_temp_cat",  "conf.high")),
                        ecol_temp_cat ~ value + variable, value.var=c("conf.high"))
  lmdata.AMRW<-cbind(lmdata.AMRW, lmdata.AMRWciL[,c(2:4)],lmdata.AMRWciH[,c(2:4)])
  colnames(lmdata.AMRW)<-c("ecol_temp_cat", "intercept", "slope", "temp", "intercept.ciL", "slope.ciL", "temp.ciL", "intercept.ciH", "slope.ciH", "temp.ciH")
  
  
  lmdata.RMRW<-dcast(melt(lmdata.RMR[, c("ecol_temp_cat", "estimate", "term")], id.vars=c("ecol_temp_cat", "estimate")),
                     ecol_temp_cat ~ value + variable, value.var=c("estimate"))
  lmdata.RMRWciL<-dcast(melt(lmdata.RMR[, c("ecol_temp_cat",  "term", "conf.low")], id.vars=c("ecol_temp_cat",  "conf.low")),
                        ecol_temp_cat ~ value + variable, value.var=c("conf.low"))
  lmdata.RMRWciH<-dcast(melt(lmdata.RMR[, c("ecol_temp_cat",  "term", "conf.high")], id.vars=c("ecol_temp_cat",  "conf.high")),
                        ecol_temp_cat ~ value + variable, value.var=c("conf.high"))
  lmdata.RMRW<-cbind(lmdata.RMRW, lmdata.RMRWciL[,c(2:4)],lmdata.RMRWciH[,c(2:4)])
  colnames(lmdata.RMRW)<-c("ecol_temp_cat", "intercept", "slope", "temp", "intercept.ciL", "slope.ciL", "temp.ciL", "intercept.ciH", "slope.ciH", "temp.ciH")
  
  lmdata.FASW<-dcast(melt(lmdata.FAS[, c("ecol_temp_cat", "estimate", "term")], id.vars=c("ecol_temp_cat", "estimate")),
                     ecol_temp_cat ~ value + variable, value.var=c("estimate"))
  lmdata.FASWciL<-dcast(melt(lmdata.FAS[, c("ecol_temp_cat",  "term", "conf.low")], id.vars=c("ecol_temp_cat",  "conf.low")),
                        ecol_temp_cat ~ value + variable, value.var=c("conf.low"))
  lmdata.FASWciH<-dcast(melt(lmdata.FAS[, c("ecol_temp_cat",  "term", "conf.high")], id.vars=c("ecol_temp_cat",  "conf.high")),
                        ecol_temp_cat ~ value + variable, value.var=c("conf.high"))
  lmdata.FASW<-cbind(lmdata.FASW, lmdata.FASWciL[,c(2:4)],lmdata.FASWciH[,c(2:4)])
  colnames(lmdata.FASW)<-c("ecol_temp_cat", "intercept", "slope", "temp", "intercept.ciL", "slope.ciL", "temp.ciL", "intercept.ciH", "slope.ciH", "temp.ciH")
  
  
  AMR_sum<-merge(lmdata.AMRW, summarylm.AMR, by = "ecol_temp_cat", all.x = TRUE)
  RMR_sum<-merge(lmdata.RMRW, summarylm.RMR, by = "ecol_temp_cat", all.x = TRUE)
  FAS_sum<-merge(lmdata.FASW, summarylm.FAS, by = "ecol_temp_cat", all.x = TRUE)
  
  AMR_sum$MR<-"AMR"
  RMR_sum$MR<-"RMR"
  FAS_sum$MR<-"FAS"
  sum<-rbind(AMR_sum, RMR_sum, FAS_sum)
  sum$Ecol<-NAMEecol
  sum$n_data_n_species<-paste(sum$n, " (",  sum$species_n, ")", sep="")
  sum<-sum %>% mutate_if(is.numeric, ~round(., 2))
  
  sum$scalingGood<-NA
  for(i in 1:nrow(sum)){
    if(sum$n[i]<100 | c(sum$size_max[i] <= sum$size_min[i]*10^1)){
      sum$scalingGood[i]<-0
    }else{
      sum$scalingGood[i]<-1
    }
  }
  
  return(sum)
}


# *************************
# demersal and pelagic 
lmdata.AMR <- data.amr %>%
  dplyr:::group_by(DemersPelag, test_category3) %>%
  group_modify(~ broom::tidy(lm(lnAMR ~ lnBWg + tempTest, data = .x), conf.int = TRUE)) %>% 
  as.data.frame()

summarylm.AMR<-  data.amr %>%
  dplyr:::group_by(DemersPelag, test_category3) %>%
  summarise(size_min = min(BW_g), size_max = max(BW_g), n = length(BW_g),
            temp_min = min(tempTest), temp_max = max(tempTest), species_n=length(unique(species)))%>% 
  as.data.frame()

lmdata.RMR <- data.rmr %>%
  dplyr:::group_by(DemersPelag, test_category3) %>%
  group_modify(~ broom::tidy(lm(lnRMR ~ lnBWg + tempTest, data = .x), conf.int = TRUE))%>% 
  as.data.frame()

summarylm.RMR<-  data.rmr %>%
  dplyr:::group_by(DemersPelag, test_category3) %>%
  summarise(size_min = min(BW_g), size_max = max(BW_g), n = length(BW_g),
            temp_min = min(tempTest), temp_max = max(tempTest), species_n=length(unique(species)))%>% 
  as.data.frame()

lmdata.FAS <- data.fas %>%
  dplyr:::group_by(DemersPelag, test_category3) %>%
  group_modify(~ broom::tidy(lm(lnFAS ~ lnBWg + tempTest, data = .x), conf.int = TRUE))%>% 
  as.data.frame()

summarylm.FAS<-  data.fas %>%
  dplyr:::group_by(DemersPelag, test_category3) %>%
  summarise(size_min = min(BW_g), size_max = max(BW_g), n = length(BW_g),
            temp_min = min(tempTest), temp_max = max(tempTest), species_n=length(unique(species)))%>% 
  as.data.frame()

DemersPelagLM <- summarise_ecologiesLM(scalingGood = FALSE, "DemersPelag",  lmdata.AMRW, summarylm.AMR, lmdata.RMRW, summarylm.RMR, lmdata.FASW, summarylm.FAS)


# *************************
# Climate
lmdata.AMR <- data.amr %>%
  dplyr:::group_by(Climate, test_category3) %>%
  group_modify(~ broom::tidy(lm(lnAMR ~ lnBWg + tempTest, data = .x), conf.int = TRUE)) %>% 
  as.data.frame()

summarylm.AMR<-  data.amr %>%
  dplyr:::group_by(Climate, test_category3) %>%
  summarise(size_min = min(BW_g), size_max = max(BW_g), n = length(BW_g),
            temp_min = min(tempTest), temp_max = max(tempTest), species_n=length(unique(species)))%>% 
  as.data.frame()

lmdata.RMR <- data.rmr %>%
  dplyr:::group_by(Climate, test_category3) %>%
  group_modify(~ broom::tidy(lm(lnRMR ~ lnBWg + tempTest, data = .x), conf.int = TRUE))%>% 
  as.data.frame()

summarylm.RMR<-  data.rmr %>%
  dplyr:::group_by(Climate, test_category3) %>%
  summarise(size_min = min(BW_g), size_max = max(BW_g), n = length(BW_g),
            temp_min = min(tempTest), temp_max = max(tempTest), species_n=length(unique(species)))%>% 
  as.data.frame()

lmdata.FAS <- data.fas %>%
  dplyr:::group_by(Climate, test_category3) %>%
  group_modify(~ broom::tidy(lm(lnFAS ~ lnBWg + tempTest, data = .x), conf.int = TRUE))%>% 
  as.data.frame()

summarylm.FAS<-  data.fas %>%
  dplyr:::group_by(Climate, test_category3) %>%
  summarise(size_min = min(BW_g), size_max = max(BW_g), n = length(BW_g),
            temp_min = min(tempTest), temp_max = max(tempTest), species_n=length(unique(species)))%>% 
  as.data.frame()

ClimateLM <- summarise_ecologiesLM(scalingGood = FALSE, "Climate", lmdata.AMRW, summarylm.AMR, lmdata.RMRW, summarylm.RMR, lmdata.FASW, summarylm.FAS)


# *************************
# salintyComb

lmdata.AMR <- data.amr %>%
  dplyr:::group_by(salintyComb, test_category3) %>%
  group_modify(~ broom::tidy(lm(lnAMR ~ lnBWg + tempTest, data = .x), conf.int = TRUE)) %>% 
  as.data.frame()

summarylm.AMR<-  data.amr %>%
  dplyr:::group_by(salintyComb, test_category3) %>%
  summarise(size_min = min(BW_g), size_max = max(BW_g), n = length(BW_g),
            temp_min = min(tempTest), temp_max = max(tempTest), species_n=length(unique(species)))%>% 
  as.data.frame()

lmdata.RMR <- data.rmr %>%
  dplyr:::group_by(salintyComb, test_category3) %>%
  group_modify(~ broom::tidy(lm(lnRMR ~ lnBWg + tempTest, data = .x), conf.int = TRUE))%>% 
  as.data.frame()

summarylm.RMR<-  data.rmr %>%
  dplyr:::group_by(salintyComb, test_category3) %>%
  summarise(size_min = min(BW_g), size_max = max(BW_g), n = length(BW_g),
            temp_min = min(tempTest), temp_max = max(tempTest), species_n=length(unique(species)))%>% 
  as.data.frame()

lmdata.FAS <- data.fas %>%
  dplyr:::group_by(salintyComb, test_category3) %>%
  group_modify(~ broom::tidy(lm(lnFAS ~ lnBWg + tempTest, data = .x), conf.int = TRUE))%>% 
  as.data.frame()

summarylm.FAS<-  data.fas %>%
  dplyr:::group_by(salintyComb, test_category3) %>%
  summarise(size_min = min(BW_g), size_max = max(BW_g), n = length(BW_g),
            temp_min = min(tempTest), temp_max = max(tempTest), species_n=length(unique(species)))%>% 
  as.data.frame()

salintyCombLM <- summarise_ecologiesLM(scalingGood = FALSE, "salintyComb", lmdata.AMRW, summarylm.AMR, lmdata.RMRW, summarylm.RMR, lmdata.FASW, summarylm.FAS)

# *************************
# BodyShapeI
lmdata.AMR <- data.amr %>%
  dplyr:::group_by(BodyShapeI, test_category3) %>%
  group_modify(~ broom::tidy(lm(lnAMR ~ lnBWg + tempTest, data = .x), conf.int = TRUE), .groups = ) %>% 
  as.data.frame()

summarylm.AMR<-  data.amr %>%
  dplyr:::group_by(BodyShapeI, test_category3) %>%
  summarise(size_min = min(BW_g), size_max = max(BW_g), n = length(BW_g),
            temp_min = min(tempTest), temp_max = max(tempTest), species_n=length(unique(species)))%>% 
  as.data.frame()

lmdata.RMR <- data.rmr %>%
  dplyr:::group_by(BodyShapeI, test_category3) %>%
  group_modify(~ broom::tidy(lm(lnRMR ~ lnBWg + tempTest, data = .x), conf.int = TRUE))%>% 
  as.data.frame()

summarylm.RMR<-  data.rmr %>%
  dplyr:::group_by(BodyShapeI, test_category3) %>%
  summarise(size_min = min(BW_g), size_max = max(BW_g), n = length(BW_g),
            temp_min = min(tempTest), temp_max = max(tempTest), species_n=length(unique(species)))%>% 
  as.data.frame()

lmdata.FAS <- data.fas %>%
  dplyr:::group_by(BodyShapeI, test_category3) %>%
  group_modify(~ broom::tidy(lm(lnFAS ~ lnBWg + tempTest, data = .x), conf.int = TRUE))%>% 
  as.data.frame()

summarylm.FAS<-  data.fas %>%
  dplyr:::group_by(BodyShapeI, test_category3) %>%
  summarise(size_min = min(BW_g), size_max = max(BW_g), n = length(BW_g),
            temp_min = min(tempTest), temp_max = max(tempTest), species_n=length(unique(species)))%>% 
  as.data.frame()

BodyShapeILM <- summarise_ecologiesLM(scalingGood = FALSE, "BodyShapeI", lmdata.AMRW, summarylm.AMR, lmdata.RMRW, summarylm.RMR, lmdata.FASW, summarylm.FAS)

colnames(BodyShapeILM)[colnames(BodyShapeILM) == "BodyShapeI"] <- "ecology_subgroup"
colnames(ClimateLM)[colnames(ClimateLM) == "Climate"] <- "ecology_subgroup"
colnames(salintyCombLM)[colnames(salintyCombLM) == "salintyComb"] <- "ecology_subgroup"
colnames(DemersPelagLM)[colnames(DemersPelagLM) == "DemersPelag"] <- "ecology_subgroup"

## Data curating and creation for plotting

ecology_data<-rbind(BodyShapeILM, ClimateLM, DemersPelagLM, salintyCombLM)
ecology_data.AMR<-ecology_data[ecology_data$MR=="AMR", ]
ecology_data.AMR$ecol_temp_cat<-factor(ecology_data.AMR$ecol_temp_cat, level = c(as.character(ecology_data.AMR$ecol_temp_cat)))

ecology_data.RMR<-ecology_data[ecology_data$MR=="RMR", ]
ecology_data.RMR$ecol_temp_cat<-factor(ecology_data.RMR$ecol_temp_cat, level = c(as.character(ecology_data.RMR$ecol_temp_cat)))

ecology_data.FAS<-ecology_data[ecology_data$MR=="FAS", ]
ecology_data.FAS$ecol_temp_cat<-factor(ecology_data.FAS$ecol_temp_cat, level = c(as.character(ecology_data.FAS$ecol_temp_cat)))

# only those that have warm temperatures
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
ecology_data.RMRd_er<-ecology_data.RMRd[grepl("optimal", as.character(ecology_data.RMRd$ecol_temp_cat)), ]
ecology_data.AMRd_er<-ecology_data.AMRd[grepl("optimal", as.character(ecology_data.AMRd$ecol_temp_cat)), ]
ecology_data.FASd_er<-ecology_data.FASd[grepl("optimal", as.character(ecology_data.FASd$ecol_temp_cat)), ]

order.ecology_data.FASd.W <- order.ecology_data.FASd[grepl("warm", as.character(order.ecology_data.FASd$ecol_temp_cat)), ]
order.ecology_data.FASd.ER <- order.ecology_data.FASd[grepl("optimal", as.character(order.ecology_data.FASd$ecol_temp_cat)), ]

ecology_data.RMRd_er$DIFF_mmr_rmr<-ecology_data.RMRd_er$slope - ecology_data.AMRd_er$slope
# ecology_data.RMRd_er[ecology_data.RMRd_er$DIFF_mmr_rmr > 0,] # bMMR < bRMR
# ecology_data.RMRd_er[ecology_data.RMRd_er$DIFF_mmr_rmr < 0,] # bMMR > bRMR
# 
ecology_data.RMRd_w$DIFF_mmr_rmr<-ecology_data.RMRd_w$slope - ecology_data.AMRd_w$slope
# ecology_data.RMRd_w[ecology_data.RMRd_w$DIFF_mmr_rmr > 0,] # bMMR < bRMR
# ecology_data.RMRd_w[ecology_data.RMRd_w$DIFF_mmr_rmr < 0,] # bMMR > bRMR

# ecology_data.AMRd_w[ecology_data.AMRd_w$Ecol == "salintyComb",]
# ecology_data.RMRd_w[ecology_data.RMRd_w$Ecol == "salintyComb",]
# ecology_data.RMRd_er[ecology_data.RMRd_er$Ecol == "salintyComb",]
# ecology_data.AMRd_er[ecology_data.AMRd_er$Ecol == "salintyComb",]

# g for 'good' for scaling purposes 
ecology_data.RMRd.g<-ecology_data.RMRd[c(ecology_data.RMRd$ecology_subgroup=="fusiform" | 
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

# saving 
write.csv(file = "./Data_exports/Ecologies/ecologies_FAS.csv", ecology_data.FAS, row.names=FALSE)
write.csv(file = "./Data_exports/Ecologies/ecologies_RMR.csv", ecology_data.RMR, row.names=FALSE)
write.csv(file = "./Data_exports/Ecologies/ecologies_MMR.csv", ecology_data.AMR, row.names=FALSE)


# ************************************************************************************
# ************************************************************************************
# SPECIES - SPECIFIC (estimate scaling) -----------------------
summarise_speciesLM.AMR<-function(lmdata.AMR, summarylm.AMR){
  
  lmdata.AMR$species_temp_cat<-paste(lmdata.AMR[,1], lmdata.AMR$test_category3, sep="_")
  summarylm.AMR$species_temp_cat<-paste(summarylm.AMR[,1], summarylm.AMR$test_category3, sep="_")
  
  
  lmdata.AMRW<-dcast(melt(lmdata.AMR[, c("species_temp_cat", "estimate", "term")], id.vars=c("species_temp_cat", "estimate")),
                     species_temp_cat ~ value + variable, value.var=c("estimate"))
  lmdata.AMRWciL<-dcast(melt(lmdata.AMR[, c("species_temp_cat",  "term", "conf.low")], id.vars=c("species_temp_cat",  "conf.low")),
                        species_temp_cat ~ value + variable, value.var=c("conf.low"))
  lmdata.AMRWciH<-dcast(melt(lmdata.AMR[, c("species_temp_cat",  "term", "conf.high")], id.vars=c("species_temp_cat",  "conf.high")),
                        species_temp_cat ~ value + variable, value.var=c("conf.high"))
  lmdata.AMRW<-cbind(lmdata.AMRW, lmdata.AMRWciL[,c(2:3)],lmdata.AMRWciH[,c(2:3)])
  colnames(lmdata.AMRW)<-c("species_temp_cat", "intercept", "slope", "intercept.ciL", "slope.ciL",  "intercept.ciH", "slope.ciH" )
  

  AMR_sum<-merge(lmdata.AMRW, summarylm.AMR, by = "species_temp_cat", all.x = TRUE)

  AMR_sum$MR<-"AMR"
  AMR_sum<-AMR_sum %>% mutate_if(is.numeric, ~round(., 2))
  
  return(AMR_sum)
}
  
summarise_speciesLM.RMR<-function(lmdata.RMR, summarylm.RMR){
  
  lmdata.RMR$species_temp_cat<-paste(lmdata.RMR[,1], lmdata.RMR$test_category3, sep="_")
  summarylm.RMR$species_temp_cat<-paste(summarylm.RMR[,1], summarylm.RMR$test_category3, sep="_")

  
  lmdata.RMRW<-dcast(melt(lmdata.RMR[, c("species_temp_cat", "estimate", "term")], id.vars=c("species_temp_cat", "estimate")),
                     species_temp_cat ~ value + variable, value.var=c("estimate"))
  lmdata.RMRWciL<-dcast(melt(lmdata.RMR[, c("species_temp_cat",  "term", "conf.low")], id.vars=c("species_temp_cat",  "conf.low")),
                        species_temp_cat ~ value + variable, value.var=c("conf.low"))
  lmdata.RMRWciH<-dcast(melt(lmdata.RMR[, c("species_temp_cat",  "term", "conf.high")], id.vars=c("species_temp_cat",  "conf.high")),
                        species_temp_cat ~ value + variable, value.var=c("conf.high"))
  lmdata.RMRW<-cbind(lmdata.RMRW, lmdata.RMRWciL[,c(2:3)],lmdata.RMRWciH[,c(2:3)])
  colnames(lmdata.RMRW)<-c("species_temp_cat", "intercept", "slope",  "intercept.ciL", "slope.ciL", "intercept.ciH", "slope.ciH")
  
 
  RMR_sum<-merge(lmdata.RMRW, summarylm.RMR, by = "species_temp_cat", all.x = TRUE)
  RMR_sum$MR<-"RMR"
 
  RMR_sum<-RMR_sum %>% mutate_if(is.numeric, ~round(., 2))
  return(RMR_sum)
}

summarise_speciesLM.FAS<-function(lmdata.FAS, summarylm.FAS){
  
  lmdata.FAS$species_temp_cat<-paste(lmdata.FAS[,1], lmdata.FAS$test_category3, sep="_")
  summarylm.FAS$species_temp_cat<-paste(summarylm.FAS[,1], summarylm.FAS$test_category3, sep="_")
  
 
  lmdata.FASW<-dcast(melt(lmdata.FAS[, c("species_temp_cat", "estimate", "term")], id.vars=c("species_temp_cat", "estimate")),
                     species_temp_cat ~ value + variable, value.var=c("estimate"))
  lmdata.FASWciL<-dcast(melt(lmdata.FAS[, c("species_temp_cat",  "term", "conf.low")], id.vars=c("species_temp_cat",  "conf.low")),
                        species_temp_cat ~ value + variable, value.var=c("conf.low"))
  lmdata.FASWciH<-dcast(melt(lmdata.FAS[, c("species_temp_cat",  "term", "conf.high")], id.vars=c("species_temp_cat",  "conf.high")),
                        species_temp_cat ~ value + variable, value.var=c("conf.high"))
  lmdata.FASW<-cbind(lmdata.FASW, lmdata.FASWciL[,c(2:3)],lmdata.FASWciH[,c(2:3)])
  colnames(lmdata.FASW)<-c("species_temp_cat", "intercept", "slope", "intercept.ciL", "slope.ciL",  "intercept.ciH", "slope.ciH")
  
  FAS_sum<-merge(lmdata.FASW, summarylm.FAS, by = "species_temp_cat", all.x = TRUE)
  
 
  FAS_sum$MR<-"FAS"
  FAS_sum<-FAS_sum %>% mutate_if(is.numeric, ~round(., 2))

  return(FAS_sum)
}


lmdata.FAS <- data.fas %>%
  dplyr:::group_by(species, test_category3) %>%
  group_modify(~ broom::tidy(lm(lnFAS ~ lnBWg , data = .x), conf.int = TRUE))%>% 
  as.data.frame()

summarylm.FAS<-  data.fas %>%
  dplyr:::group_by(species, test_category3) %>%
  summarise(size_min = min(BW_g), size_max = max(BW_g), n = length(BW_g),
            temp_min = min(tempTest), temp_max = max(tempTest))%>% 
  as.data.frame()

lmdata.RMR <- data.rmr %>%
  dplyr:::group_by(species, test_category3) %>%
  group_modify(~ broom::tidy(lm(lnRMR ~ lnBWg , data = .x), conf.int = TRUE))%>% 
  as.data.frame()

summarylm.RMR <- data.rmr %>%
  dplyr:::group_by(species, test_category3) %>%
  summarise(size_min = min(BW_g), size_max = max(BW_g), n = length(BW_g),
            temp_min = min(tempTest), temp_max = max(tempTest))%>% 
  as.data.frame()

lmdata.AMR <- data.amr %>%
  dplyr:::group_by(species, test_category3) %>%
  group_modify(~ broom::tidy(lm(lnAMR ~ lnBWg , data = .x), conf.int = TRUE))%>% 
  as.data.frame()

summarylm.AMR <- data.amr %>%
  dplyr:::group_by(species, test_category3) %>%
  summarise(size_min = min(BW_g), size_max = max(BW_g), n = length(BW_g),
            temp_min = min(tempTest), temp_max = max(tempTest))%>% 
  as.data.frame()

species.fas <- summarise_speciesLM.FAS(lmdata.FAS, summarylm.FAS)
species.rmr <- summarise_speciesLM.RMR(lmdata.RMR, summarylm.RMR)
species.amr <- summarise_speciesLM.AMR(lmdata.AMR, summarylm.AMR)

# Data curating for plotting
species.fas.w<-species.fas[species.fas$test_category3=="warm", ]
species.fas.er<-species.fas[species.fas$test_category3=="optimal", ]
species.rmr.w<-species.rmr[species.rmr$test_category3=="warm", ]
species.rmr.er<-species.rmr[species.rmr$test_category3=="optimal", ]
species.amr.w<-species.amr[species.amr$test_category3=="warm", ]
species.amr.er<-species.amr[species.amr$test_category3=="optimal", ]

# saving 
write.csv(file = "./Data_exports/Species/species_FAS.csv", species.fas, row.names=FALSE)
write.csv(file = "./Data_exports/Species/species_RMR.csv", species.rmr, row.names=FALSE)
write.csv(file = "./Data_exports/Species/species_MMR.csv", species.amr, row.names=FALSE)

# *****************************************************************************************
# *****************************************************************************************
# FIGURES ------------

## SUPPLEMENT Ecologies: MMR and RMR trends (warm, optimal)  -------
ecol_amr_sum1<-ggplot(ecology_data.AMRd.g, aes(color = test_category3)) +
  geom_segment(aes(x = log(size_min),
                   xend = log(size_max),
                   y = intercept + slope*log(size_min),
                   yend = intercept + slope*log(size_max)) ,
               show.legend = FALSE )+
  ylim(x = -7, 12 )+
  xlim(x = -7, 12)+
  scale_color_manual(values = c("black", cols.amr[2]))
ggformat(ecol_amr_sum1,
         x_title=expression(italic(ln)*Body~mass~(g)),
         y_title=expression(italic(ln)*MMR~(mg~O[2]~h^-1)),
         print = F)

ecol_amr_sum1<-ecol_amr_sum1 +
  ggtitle("ECOLOGIES") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.margin = margin(unit(c(0,0,-2.5,0), "cm")),
    # axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    # axis.text.y = element_blank(),
    axis.text.x = element_blank()
    )


ecol_rmr_sum1<-ggplot(ecology_data.RMRd.g, aes(color = test_category3)) +
  geom_segment(aes(x = log(size_min),
                   xend = log(size_max),
                   y = intercept + slope*log(size_min), 
                   yend = intercept + slope*log(size_max)),
               show.legend = FALSE )+
  ylim(x = -7,12 )+
  xlim(x = -7, 12)+
  scale_color_manual(values = c("black", cols.rmr[2]))
ggformat(ecol_rmr_sum1, x_title=expression(italic(ln)*Body~mass~(g)),
         y_title=expression(italic(ln)*RMR~(mg~O[2]~h^-1)), print = F)
ecol_rmr_sum1<-ecol_rmr_sum1 +
  theme(plot.margin = margin(unit(c(-2.5,0,0,0), "cm")))


# Species specific groupings: 
species_amr_sum1<-ggplot(species.amr, aes(color = test_category3)) +
  geom_segment(aes(x = log(size_min),
                   xend = log(size_max),
                   y = intercept + slope*log(size_min),
                   yend = intercept + slope*log(size_max)) ,
               show.legend = FALSE )+
  ylim(x = -7, 12 )+
  xlim(x = -7, 12)+
  scale_color_manual(values = c("black", cols.amr[2]))
ggformat(species_amr_sum1, x_title=expression(italic(ln)*Body~mass~(g)),
         y_title=expression(italic(ln)*MMR~(mg~O[2]~h^-1)), print = F)

species_amr_sum1<-species_amr_sum1 +
  ggtitle("SPECIES") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.margin = margin(unit(c(0,0.5,-2.5,0), "cm")),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank()
    )    

species_rmr_sum1<-ggplot(species.rmr, aes(color = test_category3)) +
  geom_segment(aes(x = log(size_min),
                   xend = log(size_max),
                   y = intercept + slope*log(size_min),
                   yend = intercept + slope*log(size_max)) ,
               show.legend = FALSE )+
  ylim(x = -7, 12 )+
  xlim(x = -7, 12)+
  scale_color_manual(values = c("black", cols.rmr[2]))
ggformat(species_rmr_sum1, x_title=expression(italic(ln)*Body~mass~(g)),
         y_title=expression(italic(ln)*RMR~(mg~O[2]~h^-1)), print = F )

species_rmr_sum1<-species_rmr_sum1 +
  theme(plot.margin = margin(unit(c(-2.5,0.5,0,0), "cm")), 
    axis.title.y = element_blank(),
    axis.text.y = element_blank())


grid2<-cowplot:::plot_grid(ecol_amr_sum1, species_amr_sum1,
                          ecol_rmr_sum1, species_rmr_sum1,
                          nrow = 2,
                          ncol = 2,
                          labels = c("B", "C", "D", "E"),
                          label_size = 13, 
                          label_x = c(0.24, 0.05),
                          label_y = c(0.87, 0.88),
                          rel_widths = c(1.3, 1), 
                          rel_heights = c(1, 1.3))
ggsave(filename = paste("./Figures/Suppl_Figurex.png", sep=""), grid2,
       width = 6, height = 6, units = "in")

  

 ## Ecologies: FAS trends (warm, optimal) -------
rmr.WARM.ER.good<-
  ggplot(data=ecology_data.RMRd.g,
                         aes(x=slope,
                             y=factor(ecol_temp_cat, level = c(as.character(ecol_temp_cat))),
                             color=interaction(test_category3, MR),
                             fill=interaction(test_category3, MR), 
                             group= ecology_subgroup,
                             size = MR))+
  scale_y_discrete(limits=rev,
    labels = c("","Any salinity", 
               "", "Marine",
               "","Reef-associated",
              "",  "Demersal",
              "", "Benthopelagic", 
               "","Tropical",
               "","Temperate",
               "","Subtropical",
               "","Fusiform"))+
  geom_text(aes(label = slope,  x = 0.25, color=interaction(test_category3, MR)),
            family = "Helvetica",  size=3.5, hjust = 1)+
  geom_text(data = ecology_data.AMRd.g, aes(label = slope,  x = 0.3),
            family = "Helvetica",  size=3.5, hjust = 0)+
  
  geom_line(arrow = arrow(length=unit(0,"cm"), ends="last", type = "closed"),
            linewidth = 0.7, color = "grey", alpha = 0.7)+
  geom_line(data = ecology_data.AMRd.g,
            arrow = arrow(length=unit(0,"cm"), ends="last", type = "closed"),
            linewidth = 0.7, color = "grey", alpha = 0.7)+

  geom_text(data = ecology_data.RMRd.g,
            mapping = aes(label =  n_data_n_species, x = 1.5, color=interaction(test_category3, MR)),
            family = "Helvetica", size=3.5, hjust = 1)+
  geom_text(data = ecology_data.AMRd.g,
            mapping = aes(label =  n_data_n_species, x = 1.55, color=interaction(test_category3, MR)),
            family = "Helvetica", size=3.5, hjust = 0)+

  annotate("text", x = 0.135, y = 19, color = "black", label = expression(italic(b)[RMR]~italic(b)[MMR]),
            family = "Helvetica", size=3.5, hjust = 0, parse = T)+
  annotate("text", x = 1.05, y = 19, color = "black", label = expression(samples~sizes~(data*(species))),
            family = "Helvetica", size=3.5, hjust = 0, parse = T)+
  annotate("text", x = 0.45, y = 18, color = "black", label = 'OPTIMAL',
            family = "Helvetica", size=3.5, hjust = 0, parse = T)+
  annotate("text", x = 0.45, y = 17, color = "black", label = 'WARM',
            family = "Helvetica", size=3.5, hjust = 0, parse = T)+
  geom_linerange(data = ecology_data.RMRd.g,
                 aes(xmin = as.numeric(slope.ciL),
                     xmax = as.numeric(slope.ciH)),
                 alpha = 1, linewidth = 0.5)+
  geom_linerange(data = ecology_data.AMRd.g,
                 mapping = aes(xmin = as.numeric(slope.ciL),
                               xmax = as.numeric(slope.ciH)),
                 alpha = 1, linewidth = 0.5)+
  
  geom_vline(xintercept = 1, lty = "dotted", color="grey50")+
  geom_vline(xintercept = 0.75, lty = "dotted", color="grey50")+

  geom_point(pch=21, stroke=0.1, alpha=0.7, color = "black")+
  geom_point(data = ecology_data.AMRd.g,
             pch=21,  stroke=0.1, alpha=0.7, color = "black")+

  scale_color_manual(values=c("black",cols.rmr[1],"black", cols.amr[2]))+
  scale_fill_manual(values=c("black",cols.rmr[1], "black", cols.amr[2]))+
  xlab(expression(Slope~value~(italic(b))))+
  scale_x_continuous(limits = c(0.13,1.8))+
  scale_size_manual(values = c(4, 2))+
  # scale_y_discrete(limits=rev)+
  theme_classic()+
  theme(axis.text.y = element_text( family="Helvetica", color = "black", size = 12,
                                    vjust = 1.2, hjust = 1),
        axis.text.x = element_text(family="Helvetica",  color = "black", size = 12),
        axis.title.y = element_blank(),
        legend.position = "none",
        axis.line.y=element_line(colour = 'black',size=0.5),
        axis.line.x=element_line(colour = 'black',size=0.5),
        axis.ticks.y=element_line(size=0.5),
        axis.ticks.x=element_line(size=0),
        text=element_text(size=12,  family="Helvetica"), 
        plot.margin = margin(2, 1, 1, 1, "cm"))+
  coord_cartesian(clip = "off")

# rmr.WARM.ER.good


fas.WARM.ER.good<-ggplot(data=ecology_data.FASd.g,
                         aes(x=slope,
                             y=factor(ecol_temp_cat, level = c(as.character(ecol_temp_cat))),
                             color=test_category3,
                             fill=test_category3, 
                             group= ecology_subgroup))+
  scale_y_discrete(limits=rev,
    labels = c("","Any salinity", 
               "", "Marine",
               "","Reef-associated",
              "",  "Demersal",
              "", "Benthopelagic", 
               "","Tropical",
               "","Temperate",
               "","Subtropical",
               "","Fusiform"))+
  geom_text(aes(label = n_data_n_species,  x = 0.2), family = "Helvetica", size=3.5, hjust=0)+
  geom_text(aes(label = slope, x = -0.2, family = "Helvetica"), size=3.5, hjust=1)+
  geom_line(arrow = arrow(length=unit(0,"cm"), ends="first", type = "closed"),
            linewidth = 0.7, color = "grey", alpha = 0.7)+
  geom_linerange(data = ecology_data.FASd.g,
                 aes(xmin = as.numeric(slope.ciL) ,
                     xmax = as.numeric(slope.ciH)),
                 alpha = 1)+
  annotate("text", x = -0.27, y = 19, color = "black", label = expression(italic(b)[FAS]),
            family = "Helvetica", size=3.5, hjust = 0, parse = T)+
  annotate("text", x = 0.023, y = 19, color = "black", label = expression(samples~sizes~(data*(species))),
            family = "Helvetica", size=3.5, hjust = 0, parse = T)+
  geom_vline(xintercept = 0, lty = "dotted", color="grey50")+
  geom_linerange(aes(xmin = slope.ciL, xmax = slope.ciH), alpha = 1)+
  geom_point(pch=21, size=3, stroke=0.1, alpha=0.7, color = "black")+
  scale_color_manual(values=c("black", cols.fas[1]))+
  scale_fill_manual(values=c("black", cols.fas[2]))+
  xlab(expression(Slope~value~(italic(b))))+
  scale_x_continuous(limits = c(-0.3,0.4))+
  # scale_y_discrete(limits=rev)+
  theme_classic()+
  theme(
    # axis.text.y = element_text(face = "italic", color = "black", size = 15),
        axis.text.x = element_text(family="Helvetica",  color = "black", size = 12),
        axis.title.y = element_blank(),
        legend.position = "none",
        axis.line.y=element_line(colour = 'black',size=0.5),
        axis.line.x=element_line(colour = 'black',size=0.5),
        axis.ticks.y=element_line(size=0.5),
        axis.text.y = element_blank(),
        axis.ticks.x=element_line(size=0),
        text=element_text(size=12,  family="Helvetica"),
        plot.margin = margin(2, 1, 1, 0, "cm"))+
  coord_cartesian(clip = "off")
# fas.WARM.ER.good

cowplot::plot_grid(rmr.WARM.ER.good,
                   fas.WARM.ER.good,
                  nrow = 1, 
                  labels = "AUTO", 
                  rel_widths = c(1, 0.5),
                  label_x = c(0.22, -00.1),
                  label_y = c(0.89, 0.89)) 
ggsave(filename = paste("./Figures/Figure4.png", sep=""),
         width = 9.5, height = 6)


## Ecologies: Marginal estimates, all MR performances ------
# 
# ggplot()+
#   geom_pointrange(data = ecol.emmeans, 
#                   mapping = aes(y=EcolGroup, x = exp(emmean),
#                                 xmin = exp(lower.CL), xmax = exp(upper.CL), 
#                                 fill = var, shape = test.categ), 
#                   size = 1, alpha = 0.5)+
#   # scale_fill_manual(values =c(cols.as[1],cols.fas[1], cols.amr[1], cols.rmr[1]))+
#   xlim(0, 20)+
#   theme_pubr()+
#   scale_shape_manual(c(21, 22, 23, 24))+
#   facet_wrap(var~., nrow=4)

## Species ------

sp1_amr<-ggplot(data=species.amr, aes(y=slope, x=fct_reorder(species, (slope.ciL-slope.ciH)), color = test_category3))+
  geom_hline(yintercept = 1, colour="black", lty="solid")+
  geom_hline(yintercept = 0.75, colour="black", lty="solid")+
  ylim(min(species.amr$slope.ciL)-0.1, max(species.amr$slope.ciH)+0.1)+
  geom_linerange(aes(ymin = slope, ymax = slope.ciH, color = test_category3),
                 size=0.6, alpha=1)+
  geom_linerange(aes(ymin = slope, ymax = slope.ciL, color = test_category3),
                 size=0.6, alpha=1)+
  geom_point(pch=1, color="black", size=1)+
  geom_text(aes(label = n, y = -3.4, color = test_category3), hjust = "inward",  family = "Helvetica", size=3)+
  ylab(expression(Slope~value~(italic(b))))+
  scale_color_manual(values = c("black", cols.amr[2]))+
  facet_grid(.~test_category3)+
  coord_flip()+
  theme_pubr()+
  theme(axis.text.y = element_text(face = "italic", color = "black", size = 10),
        axis.text.x = element_text(size = 13),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(5.5, 5.8, 5.5, 5.5), "pt"),
        text=element_text(size=16,  family="Helvetica"))
ggsave(filename = paste("./Figures/suppl_speciesSlopesMMR.png", sep=""),
       plot=sp1_amr, width = 9, height = 9, units = "in")


sp1_rmr<-ggplot(data=species.rmr, aes(y=slope, x=fct_reorder(species, (slope.ciL-slope.ciH)), color = test_category3))+
  geom_hline(yintercept = 1, colour="black", lty="solid")+
  geom_hline(yintercept = 0.75, colour="black", lty="solid")+  ylim(min(species.rmr$slope.ciL)-0.1, max(species.rmr$slope.ciH)+0.1)+
  geom_linerange(aes(ymin = slope, ymax = slope.ciH, color = test_category3),
                 size=0.6, alpha=1)+
  geom_linerange(aes(ymin = slope, ymax = slope.ciL, color = test_category3),
                 size=0.6, alpha=1)+
  geom_point(pch=1, color="black", size=1)+
  geom_text(aes(label = n, y = -4, color = test_category3), hjust = "inward",  family = "Helvetica", size=3)+
  ylab(expression(Slope~value~(italic(b))))+
  scale_color_manual(values = c("black", cols.rmr[2]))+
  facet_grid(.~test_category3)+
  coord_flip()+
  theme_pubr()+
  theme(axis.text.y = element_text(face = "italic", color = "black", size = 10),
        axis.text.x = element_text(size = 13),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(5.5, 5.8, 5.5, 5.5), "pt"),
        text=element_text(size=16,  family="Helvetica"))
ggsave(filename = paste("./Figures/suppl_speciesSlopesRMR.png", sep=""),
       plot=sp1_rmr, width = 9, height = 10, units = "in")


sp1_fas<-ggplot(data=species.fas, aes(y=slope, x=fct_reorder(species, (slope.ciL-slope.ciH)), color = test_category3))+
  geom_hline(yintercept = 0, colour="black", lty="solid")+
  ylim(min(species.fas$slope.ciL)-0.1, max(species.fas$slope.ciH)+0.1)+
  geom_linerange(aes(ymin = slope, ymax = slope.ciH, color = test_category3),
                 size=0.6, alpha=1)+
  geom_linerange(aes(ymin = slope, ymax = slope.ciL, color = test_category3),
                 size=0.6, alpha=1)+
  geom_point(pch=1, color="black", size=1)+
  geom_text(aes(label = n, y = -5, color = test_category3), hjust = "inward",  family = "Helvetica", size=3)+
  ylab(expression(Slope~value~(italic(b))))+
  scale_color_manual(values = c("black", cols.fas[2]))+
  facet_grid(.~test_category3)+
  coord_flip()+
  theme_pubr()+
  theme(axis.text.y = element_text(face = "italic", color = "black", size = 10),
        axis.text.x = element_text(size = 13),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(5.5, 5.8, 5.5, 5.5), "pt"),
        text=element_text(size=16,  family="Helvetica"))
ggsave(filename = paste("./Figures/suppl_speciesSlopesFAS.png", sep=""),
       plot=sp1_fas, width = 9, height = 8.5, units = "in")


 ## SUPPL: histograms of slopes -----
 ### by ecolgies =------
annotationER<-annotate(geom = "text",
                     x= 1.3, y = 4.6, 
                     label = paste(round(mean(ecology_data.AMRd_er$slope),2), " (sd ", round(sd(ecology_data.AMRd_er$slope),2) ,")", sep =""), family = "Helvetica", size=4, color="black")
annotationW<-annotate(geom = "text", 
                      x= 1.3, y = 4,
                      label = paste(round(mean(ecology_data.AMRd_w$slope),2), " (sd ", round(sd(ecology_data.AMRd_w$slope),2) ,")", sep =""),
                      family = "Helvetica", size=4,  color=cols.amr[2])

plot_hist2_amrlm<-ggplot(ecology_data.AMRd, aes(x=slope, fill = test_category3, 
                                                color = test_category3)) +
  geom_vline(xintercept = mean(ecology_data.AMRd_er$slope), color="black",
             lty="solid", lwd=1)+
  geom_vline(xintercept = mean(ecology_data.AMRd_w$slope), color=cols.amr[2],
             lty="solid", lwd=1)+
  geom_histogram(show.legend = FALSE, bins=30, alpha = 0.4, linewidth = 0.)+
  scale_color_manual(values = c("black", cols.amr[2]))+
  scale_fill_manual(values = c("black", cols.amr[2]))+
  ylim(0,5)+
  scale_x_continuous(limits=c(0.5, 1.5))+
  facet_wrap(.~test_category3)+
  at_panel(annotationER, PANEL == 1)+
  at_panel(annotationW, PANEL == 1)

ggformat(plot_hist2_amrlm, x_title = expression(Scaling~slope~(italic(b))), y_title = "Frequency", print = F)
plot_hist2_amrlm <- plot_hist2_amrlm + theme(axis.title.x = element_text(size=12,  family="Helvetica"),
                                             axis.text.x = element_text(size=12,  family="Helvetica"),
                                             axis.title.y = element_text(size=12,  family="Helvetica"),
                                             axis.text.y = element_text(size=12,  family="Helvetica"))

annotationER<-annotate(geom = "text",
                     x= 1.3, y = 4.6, 
                     label = paste(round(mean(ecology_data.RMRd_er$slope),2), " (sd ", round(sd(ecology_data.RMRd_er$slope),2) ,")", sep =""), family = "Helvetica", size=4, color="black")
annotationW<-annotate(geom = "text", 
                      x= 1.3, y = 4,
                      label = paste(round(mean(ecology_data.RMRd_w$slope),2), " (sd ", round(sd(ecology_data.RMRd_w$slope),2) ,")", sep =""),
                      family = "Helvetica", size=4,  color=cols.rmr[2])

plot_hist2_rmrlm<-ggplot(ecology_data.RMRd, aes(x=slope, fill = test_category3, 
                                                color = test_category3)) +
  geom_vline(xintercept = mean(ecology_data.RMRd_er$slope), color="black",
             lty="solid", lwd=1)+
  geom_vline(xintercept = mean(ecology_data.RMRd_w$slope), color=cols.rmr[2],
             lty="solid", lwd=1)+
  geom_histogram(show.legend = FALSE, bins=30, alpha = 0.4, linewidth = 0.)+
  scale_color_manual(values = c("black", cols.rmr[2]))+
  scale_fill_manual(values = c("black", cols.rmr[2]))+
  ylim(0,5)+
  scale_x_continuous(limits=c(0.5, 1.5))+
  facet_wrap(.~test_category3)+
  at_panel(annotationER, PANEL == 1)+
  at_panel(annotationW, PANEL == 1)

ggformat(plot_hist2_rmrlm, x_title = expression(Scaling~slope~(italic(b))), y_title = "Frequency", print = F)
plot_hist2_rmrlm <- plot_hist2_rmrlm + theme(axis.title.x = element_text(size=12,  family="Helvetica"),
                                             axis.text.x = element_text(size=12,  family="Helvetica"),
                                             axis.title.y = element_text(size=12,  family="Helvetica"),
                                             axis.text.y = element_text(size=12,  family="Helvetica"))

annotationER<-annotate(geom = "text",
                     x= -0.5, y = 7.9, 
                     label = paste(round(mean(ecology_data.FASd_er$slope),2), " (sd ", round(sd(ecology_data.FASd_er$slope),2) ,")", sep =""), family = "Helvetica", size=4, color="black")
annotationW<-annotate(geom = "text", 
                      x= -0.5, y = 7,
                      label = paste(round(mean(ecology_data.FASd_w$slope),2), " (sd ", round(sd(ecology_data.FASd_w$slope),2) ,")", sep =""),
                      family = "Helvetica", size=4,  color=cols.fas[2])

plot_hist2_faslm<-ggplot(ecology_data.FASd, aes(x=slope, fill = test_category3, 
                                                color = test_category3)) +
  geom_vline(xintercept = mean(ecology_data.FASd_er$slope), color="black",
             lty="solid", lwd=1)+
  geom_vline(xintercept = mean(ecology_data.FASd_w$slope), color=cols.fas[2],
             lty="solid", lwd=1)+
  geom_histogram(show.legend = FALSE, bins=30, alpha = 0.4, linewidth = 0.)+
  scale_color_manual(values = c("black", cols.fas[2]))+
  scale_fill_manual(values = c("black", cols.fas[2]))+
  ylim(0,8.5)+
  # scale_x_continuous(limits=c(-0.6, 0.1))+
  facet_wrap(.~test_category3)+
  at_panel(annotationER, PANEL == 1)+
  at_panel(annotationW, PANEL == 1)

ggformat(plot_hist2_faslm, x_title = expression(Scaling~slope~(italic(b))), y_title = "Frequency", print = F)
plot_hist2_faslm <- plot_hist2_faslm + theme(axis.title.x = element_text(size=12,  family="Helvetica"),
                                             axis.text.x = element_text(size=12,  family="Helvetica"),
                                             axis.title.y = element_text(size=12,  family="Helvetica"),
                                             axis.text.y = element_text(size=12,  family="Helvetica"))


ecolhist.plot<-cowplot:::plot_grid(plot_hist2_amrlm, plot_hist2_rmrlm,plot_hist2_faslm,
                                   align = "hv",
                                   axis = "l",
                                   nrow = 3,
                                   ncol = 1, 
                              labels = "AUTO") 
save_plot(filename = paste("./Figures/Suppl_ecol_hist.png", sep=""),
          ecolhist.plot, base_width = 10, base_height = 8, units = "in")



# (SUPPL) Hist of slopes for the species ---
# ### by species =------
annotationER<-annotate(geom = "text",
                     x= 2, y = 11.5,
                     label = paste(round(mean(species.amr.er$slope),2), " (", round(sd(species.amr.er$slope),2) ,")", sep =""), family = "Helvetica", size=4, color="black")
annotationW<-annotate(geom = "text", 
                      x= 2, y = 10.1,
                      label = paste(round(mean(species.amr.w$slope),2), " (", round(sd(species.amr.w$slope),2) ,")", sep =""),
                      family = "Helvetica", size=4,  color=cols.amr[2])

plot_hist3_amrlm<-ggplot(species.amr, aes(x=slope, color = test_category3, fill = test_category3)) +
  geom_vline(xintercept = mean(species.amr.er$slope), color="black", lty="solid", lwd=1)+
  geom_vline(xintercept = mean(species.amr.w$slope), color=cols.amr[2], lty="solid", lwd=1)+
  geom_histogram(show.legend = FALSE, bins=50, alpha = 0.5, linewidth = 0)+
  scale_color_manual(values = c("black", cols.amr[2]))+
  scale_fill_manual(values = c("black", cols.amr[2]))+
  ylim(0,12)+
  # scale_x_continuous(limits=c(-1, 3))+
  facet_wrap(.~test_category3)+
  at_panel(annotationER, PANEL == 1)+
  at_panel(annotationW, PANEL == 1)

 ggformat(plot_hist3_amrlm, x_title = expression(Scaling~slope~(italic(b))), y_title = "Frequency", print = FALSE)
plot_hist3_amrlm <- plot_hist3_amrlm + theme(axis.title.x = element_text(size=12,  family="Helvetica"),
                                             axis.text.x = element_text(size=12,  family="Helvetica"),
                                             axis.title.y = element_text(size=12,  family="Helvetica"),
                                             axis.text.y = element_text(size=12,  family="Helvetica"))


annotationER<-annotate(geom = "text",
                     x= -0.5, y = 11.5,
                     label = paste(round(mean(species.rmr.er$slope),2), " (", round(sd(species.rmr.er$slope),2) ,")", sep =""), family = "Helvetica", size=4, color="black")
annotationW<-annotate(geom = "text", 
                      x= -0.5, y = 10.1,
                      label = paste(round(mean(species.rmr.w$slope),2), " (", round(sd(species.rmr.w$slope),2) ,")", sep =""),
                      family = "Helvetica", size=4,  color=cols.rmr[2])

plot_hist3_rmrlm<-ggplot(species.rmr, aes(x=slope, color = test_category3, fill = test_category3)) +
  geom_vline(xintercept = mean(species.rmr.er$slope), color="black", lty="solid", lwd=1)+
  geom_vline(xintercept = mean(species.rmr.w$slope), color=cols.rmr[2], lty="solid", lwd=1)+
  geom_histogram(show.legend = FALSE, bins=50, alpha = 0.5, linewidth = 0)+
  scale_color_manual(values = c("black", cols.rmr[2]))+
  scale_fill_manual(values = c("black", cols.rmr[2]))+
  ylim(0,12)+
  # scale_x_continuous(limits=c(-1, 3))+
  facet_wrap(.~test_category3)+
  at_panel(annotationER, PANEL == 1)+
  at_panel(annotationW, PANEL == 1)

 ggformat(plot_hist3_rmrlm, x_title = expression(Scaling~slope~(italic(b))), y_title = "Frequency", print = FALSE)
plot_hist3_rmrlm <- plot_hist3_rmrlm + theme(axis.title.x = element_text(size=12,  family="Helvetica"),
                                             axis.text.x = element_text(size=12,  family="Helvetica"),
                                             axis.title.y = element_text(size=12,  family="Helvetica"),
                                             axis.text.y = element_text(size=12,  family="Helvetica"))


annotationER<-annotate(geom = "text",
                     x= 2, y = 11.5,
                     label = paste(round(mean(species.fas.er$slope),2), " (", round(sd(species.fas.er$slope),2) ,")", sep =""), family = "Helvetica", size=4, color="black")
annotationW<-annotate(geom = "text", 
                      x= 2, y = 10.1,
                      label = paste(round(mean(species.fas.w$slope),2), " (", round(sd(species.fas.w$slope),2) ,")", sep =""),
                      family = "Helvetica", size=4,  color=cols.fas[2])

plot_hist3_faslm<-ggplot(species.fas, aes(x=slope, color = test_category3, fill = test_category3)) +
  geom_vline(xintercept = mean(species.fas.er$slope), color="black", lty="solid", lwd=1)+
  geom_vline(xintercept = mean(species.fas.w$slope), color=cols.fas[2], lty="solid", lwd=1)+
  geom_histogram(show.legend = FALSE, bins=50, alpha = 0.5, linewidth = 0)+
  scale_color_manual(values = c("black", cols.fas[2]))+
  scale_fill_manual(values = c("black", cols.fas[2]))+
  ylim(0,12)+
  # scale_x_continuous(limits=c(-1.5, 3))+
  facet_wrap(.~test_category3)+
  at_panel(annotationER, PANEL == 1)+
  at_panel(annotationW, PANEL == 1)

 ggformat(plot_hist3_faslm, x_title = expression(Scaling~slope~(italic(b))), y_title = "Frequency", print = FALSE)
plot_hist3_faslm <- plot_hist3_faslm + theme(axis.title.x = element_text(size=12,  family="Helvetica"),
                                             axis.text.x = element_text(size=12,  family="Helvetica"),
                                             axis.title.y = element_text(size=12,  family="Helvetica"),
                                             axis.text.y = element_text(size=12,  family="Helvetica"))


spechist.plot<-cowplot:::plot_grid(plot_hist3_amrlm, plot_hist3_rmrlm,plot_hist3_faslm,
                              align = "hv",
                              axis = "l",
                              nrow = 3,
                              ncol = 1, 
                              labels = "AUTO") 
save_plot(filename = paste("./Figures/Suppl_species_hist.png", sep=""),
      spechist.plot, base_width = 10, base_height = 8, units = "in")


