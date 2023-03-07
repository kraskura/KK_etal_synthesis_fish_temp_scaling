
library(emmeans)
library(cowplot)
library(ggplot2)
library(mgcv)
library(myTAI)
library(lme4)
library(tidyverse)
library(reshape2)
library(ggformat2)
library(here)

# Sys.Date()
data.imports<-"2023-01-23"
getModelParameters<-"NO" 
# ^ this is computationally heavy task, to get saved data: "NO" or 2
# to get compute now: "YES" or 1 

source("./R/phylo_mixed model.R")

# 1. Summary tables: by Species 

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
  AMR_sum<-AMR_sum %>% mutate_if(is.numeric, ~round(., 3))
  
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
 
  RMR_sum<-RMR_sum %>% mutate_if(is.numeric, ~round(., 3))
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
  FAS_sum<-FAS_sum %>% mutate_if(is.numeric, ~round(., 3))

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

fasSpeciesLM <- summarise_speciesLM.FAS(lmdata.FAS, summarylm.FAS)
amrSpeciesLM <- summarise_speciesLM.AMR(lmdata.AMR, summarylm.AMR)
rmrSpeciesLM <- summarise_speciesLM.RMR(lmdata.RMR, summarylm.RMR)

write.csv(file = "./Data_exports/Species/species_FAS.csv", fasSpeciesLM, row.names=FALSE)
write.csv(file = "./Data_exports/Species/species_RMR.csv", rmrSpeciesLM, row.names=FALSE)
write.csv(file = "./Data_exports/Species/species_MMR.csv", amrSpeciesLM, row.names=FALSE)




