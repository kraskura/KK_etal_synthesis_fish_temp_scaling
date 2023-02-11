
# libraries 
library(cowplot)
library(ggplot2)
library(mgcv)
library(myTAI)
library(reshape2)
library(tidyverse)
library(forcats)
library(rfishbase)
library(ggformat2)

# cols.as<-c("#265F73", "#00929A", "#00C5A3")
# cols.fas<-c("#395200", "#89A000", "yellow")
# cols.rmr<-c("#C70039", "#FF6D7C", "#FFA3AC")
# cols.amr<-c("#00749F","#00A8D6", "#9CE9FF")
# cols<-c("#004EBF","#f45905","#004EBF","#f45905", "#8EA8FF","#FFC6A5", "#00B498", "#007E66")# AMR -rmr- AMR dark - rmr dark - light - as - fas
# 

    # # Sys.Date()
    # data.imports<-"2023-01-23"
    # getModelParameters<-"NO" 
    # # ^ this is computationally heavy task, to get saved data: "NO" or 2
    # # to get compute now: "YES" or 1 
    # 
    # setwd("/Users/kristakraskura/Github_repositories/KK_etal_synthesis_fish_temp_scaling/")
    # 
    # source("./R/get_scaling_data_temp.R")
    # source("./R/phylo_mixed model.R")
    # source("/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Codes/MR-fish-metadata/deltaIC_BICdelta_ggformat.R")
    # 
    # # function to order model selection based on the lowest BIC score
    # BICdelta<-function(BICtable){
    #   BIC.t <- BICtable [order(BICtable$BIC), ]
    #   BIC.t$delta <- round(abs(BIC.t$BIC[1] -  BIC.t$BIC), 5)
    #   return( BIC.t)
    # }
    # 
    # data.list<-get_scaling_data_temp(data.amr = "/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Data/MR-fish-metadata-data/Fish_AMR_temp_dataset_mar2022.csv",
    #                                  data.rmr = "/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Data/MR-fish-metadata-data/Fish_RMR_temp_dataset_mar2022.csv",
    #                                  ecology.data = "/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Data/MR-fish-metadata-data/Kraskura_species_ecologies_mar2022.csv", 
    #                                  onlyTop.above = TRUE,
    #                                  exp_rmr = round(fixef(rmr_mod_ER)[2],3),
    #                                  exp_amr = round(fixef(amr_mod_ER)[2],3),
    #                                  exp_as = round(fixef(as_mod_ER)[2],3),
    #                                  exp_rmr_warm = round(fixef(rmr_mod_W)[2],3),
    #                                  exp_amr_warm = round(fixef(amr_mod_W)[2],3),
    #                                  exp_as_warm = round(fixef(as_mod_W)[2],3))
    # 
    # data.amrAC<-data.frame(data.list[1]) 
    # data.rmrAC<-data.frame(data.list[2])
    # data.amrAM<-data.frame(data.list[3])
    # data.rmrAM<-data.frame(data.list[4])
    # data.amrER<-data.frame(data.list[5])
    # data.rmrER<-data.frame(data.list[6])
    # 
    # data.asAC<-data.frame(data.list[7]) 
    # data.fasAC<-data.frame(data.list[8])
    # data.asAM<-data.frame(data.list[9])
    # data.fasAM<-data.frame(data.list[10])
    # data.asER<-data.frame(data.list[11])
    # data.fasER<-data.frame(data.list[12])
    # 
    # data.amr<-data.frame(data.list[13]) 
    # data.rmr<-data.frame(data.list[14]) 
    # dataMR<-data.frame(data.list[15])
    # data.as<-data.frame(data.list[16])
    # data.fas<-data.frame(data.list[17])
    # 
    # data.amr$test_category3<-"ecol_relev"
    # data.amr$test_category3[data.amr$test_category=="acute"] <- "warm"
    # data.amr$test_category3[data.amr$test_category=="acclim"] <- "warm"
    # 
    # data.rmr$test_category3<-"ecol_relev"
    # data.rmr$test_category3[data.rmr$test_category=="acute"] <- "warm"
    # data.rmr$test_category3[data.rmr$test_category=="acclim"] <- "warm"
    # 
    # data.fas$test_category3<-"ecol_relev"
    # data.fas$test_category3[data.fas$test_category=="acute"] <- "warm"
    # data.fas$test_category3[data.fas$test_category=="acclim"] <- "warm"
    # 
    # k<-(8.62*10^(-5)) # Boltzmann's constant
    # E<-0.63 # activation energy MTE
    # 
    # data.amr$test_category2<-"acclim"
    # data.amr$test_category2[data.amr$test_category=="acute"] <- "acute"
    # data.amr$test_category3<-"ecol_relev"
    # data.amr$test_category3[data.amr$test_category=="acute"] <- "warm"
    # data.amr$test_category3[data.amr$test_category=="acclim"] <- "warm"
    # data.rmr$test_category2<-"acclim"
    # data.rmr$test_category2[data.rmr$test_category=="acute"] <- "acute"
    # data.rmr$test_category3<-"ecol_relev"
    # data.rmr$test_category3[data.rmr$test_category=="acute"] <- "warm"
    # data.rmr$test_category3[data.rmr$test_category=="acclim"] <- "warm"
    # data.fas$test_category2<-"acclim"
    # data.fas$test_category2[data.fas$test_category=="acute"] <- "acute"
    # data.fas$test_category3<-"ecol_relev"
    # data.fas$test_category3[data.fas$test_category=="acute"] <- "warm"
    # data.fas$test_category3[data.fas$test_category=="acclim"] <- "warm"
    # data.as$test_category2<-"acclim"
    # data.as$test_category2[data.as$test_category=="acute"] <- "acute"
    # data.as$test_category3<-"ecol_relev"
    # data.as$test_category3[data.as$test_category=="acute"] <- "warm"
    # data.as$test_category3[data.as$test_category=="acclim"] <- "warm"
    # 
    # # these datasets are > Topt temperatures
    # data.amr.test<-rbind(data.amrAC, data.amrAM)
    # data.rmr.test<-rbind(data.rmrAC, data.rmrAM)
    # data.fas.test<-data.amr.test[c(!is.na(data.amr.test$FAS) & is.finite(data.amr.test$FAS)) , ]
    # data.as.test<-data.amr.test[c(!is.na(data.amr.test$lnAS) & is.finite(data.amr.test$lnAS)) , ]
    # 
    # data.amrER$tempTestK1000_inC<-((1000/data.amrER$tempTestK1000))-273.15

# sanity checks
which(is.na(data.amr$DemersPelag))
which(is.na(data.amr$BodyShapeI))
which(is.na(data.amr$Climate))
which(is.na(data.amr$salintyComb))
data.amr[which(is.na(data.amr$salintyComb)),]

# 1. SIMPLE LINEAR MODELS -- Summary tables: by Ecology 

summarise_ecologiesLM<-function(scalingGood=FALSE, NAMEecol, lmdata.AMRW, summarylm.AMR, lmdata.RMRW, summarylm.RMR, lmdata.FASW, summarylm.FAS){
  
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
  sum<-sum %>% mutate_if(is.numeric, ~round(., 3))
  
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
# *************************
# BodyShapeI


lmdata.AMR <- data.amr %>%
  dplyr:::group_by(BodyShapeI, test_category3) %>%
  group_modify(~ broom::tidy(lm(lnAMR ~ lnBWg + tempTest, data = .x), conf.int = TRUE)) %>% 
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



# 

## DATA curating and creation for plotting ------

ecology_data<-rbind(BodyShapeILM, ClimateLM, DemersPelagLM, salintyCombLM)
ecology_data.AMR<-ecology_data[ecology_data$MR=="AMR", ]
ecology_data.AMR$ecol_temp_cat<-factor(ecology_data.AMR$ecol_temp_cat, level = c(as.character(ecology_data.AMR$ecol_temp_cat)))

ecology_data.RMR<-ecology_data[ecology_data$MR=="RMR", ]
ecology_data.RMR$ecol_temp_cat<-factor(ecology_data.RMR$ecol_temp_cat, level = c(as.character(ecology_data.RMR$ecol_temp_cat)))

ecology_data.FAS<-ecology_data[ecology_data$MR=="FAS", ]
ecology_data.FAS$ecol_temp_cat<-factor(ecology_data.FAS$ecol_temp_cat, level = c(as.character(ecology_data.FAS$ecol_temp_cat)))
 
# saving 
write.csv(file = "./Data_exports/Ecologies/ecologies_FAS.csv", ecology_data.FAS, row.names=FALSE)
write.csv(file = "./Data_exports/Ecologies/ecologies_RMR.csv", ecology_data.RMR, row.names=FALSE)
write.csv(file = "./Data_exports/Ecologies/ecologies_MMR.csv", ecology_data.AMR, row.names=FALSE)



