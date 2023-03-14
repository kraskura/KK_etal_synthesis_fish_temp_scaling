

k<-(8.62*10^(-5)) # Boltzmann's constant
E<-0.63 # activation energy MTE

library(reshape2)
# ******************************************************************************************************************************************************
# ******************************************************************************************************************************************************
# GLOBAL DATA (with mass specific values) ---------

data.list<-get_data_temp(data.amr = "./Data/Fish_AMR_temp_dataset_mar2022.csv",
                                 data.rmr = "./Data/Fish_RMR_temp_dataset_mar2022.csv",
                                 ecology.data = "./Data/Kraskura_species_ecologies_mar2022.csv", 
                                 onlyTop.above = TRUE,
                                 calc_mass_specific = TRUE,
                                 exp_rmr = RMR_slope,
                                 exp_amr = round(AMR.slopes$lnBWg.trend[3],3), # at 20ºC
                                 exp_as = AS_slope,
                                 exp_rmr_warm = RMR_slope_w,
                                 exp_amr_warm = MMR_slope_w,
                                 exp_as_warm = AS_slope_w)

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

# acclimated and warm are > Topt temperatures
data.amr.test<-rbind(data.amrAC, data.amrAM)
data.rmr.test<-rbind(data.rmrAC, data.rmrAM)
data.fas.test<-data.amr.test[c(!is.na(data.amr.test$FAS) & is.finite(data.amr.test$FAS)) , ]
data.as.test<-data.amr.test[c(!is.na(data.amr.test$lnAS) & is.finite(data.amr.test$lnAS)) , ]

data.amrER$tempTestK1000_inC<-((1000/data.amrER$tempTestK1000))-273.15

# ******************************************************************************************************************************************************
# ******************************************************************************************************************************************************
# ECOLOGY - SPECIFIC ----
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

## Data curating and creation for plotting

ecology_data<-rbind(BodyShapeILM, ClimateLM, DemersPelagLM, salintyCombLM)
ecology_data.AMR<-ecology_data[ecology_data$MR=="AMR", ]
ecology_data.AMR$ecol_temp_cat<-factor(ecology_data.AMR$ecol_temp_cat, level = c(as.character(ecology_data.AMR$ecol_temp_cat)))

ecology_data.RMR<-ecology_data[ecology_data$MR=="RMR", ]
ecology_data.RMR$ecol_temp_cat<-factor(ecology_data.RMR$ecol_temp_cat, level = c(as.character(ecology_data.RMR$ecol_temp_cat)))

ecology_data.FAS<-ecology_data[ecology_data$MR=="FAS", ]
ecology_data.FAS$ecol_temp_cat<-factor(ecology_data.FAS$ecol_temp_cat, level = c(as.character(ecology_data.FAS$ecol_temp_cat)))

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


# ******************************************************************************************************************************************************
# ******************************************************************************************************************************************************
# SPECIES - SPECIFIC -----------------------
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

species.fas <- summarise_speciesLM.FAS(lmdata.FAS, summarylm.FAS)
species.rmr <- summarise_speciesLM.RMR(lmdata.RMR, summarylm.RMR)
species.amr <- summarise_speciesLM.AMR(lmdata.AMR, summarylm.AMR)

# Data curating for plotting
species.fas.w<-species.fas[species.fas$test_category3=="warm", ]
species.fas.er<-species.fas[species.fas$test_category3=="ecol_relev", ]
species.rmr.w<-species.rmr[species.rmr$test_category3=="warm", ]
species.rmr.er<-species.rmr[species.rmr$test_category3=="ecol_relev", ]
species.amr.w<-species.amr[species.amr$test_category3=="warm", ]
species.amr.er<-species.amr[species.amr$test_category3=="ecol_relev", ]

# saving 
write.csv(file = "./Data_exports/Species/species_FAS.csv", species.fas, row.names=FALSE)
write.csv(file = "./Data_exports/Species/species_RMR.csv", species.rmr, row.names=FALSE)
write.csv(file = "./Data_exports/Species/species_MMR.csv", species.amr, row.names=FALSE)


# Figures -------------
 
data.amr$tempTestK1<-1/data.amr$tempTestK
data.amr$tempTestK1000<-1000/data.amr$tempTestK
    

# Ecology specific groupings (sufficient for scaling ):  -------
ecol_rmr_sum1<-ggplot(ecology_data.RMRd.g, aes(color = test_category3)) +
  geom_segment(aes(x = log(size_min), xend = log(size_max), y = intercept + slope*log(size_min), yend = intercept + slope*log(size_max)), show.legend = FALSE )+
  ylim(x = -7,12 )+
  xlim(x = -7, 12)+
  scale_color_manual(values = c("black", cols.rmr[2]))
  # facet_grid(.~test_category3)
ggformat(ecol_rmr_sum1, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*RMR~(mg~O[2]~h^-1)), print = F)
ggsave(filename = paste("./Figures/FigVAR_ecolRMR_v2_", Sys.Date(), ".png", sep=""),
       plot=ecol_rmr_sum1, width = 4, height = 4, units = "in")

ecol_amr_sum1<-ggplot(ecology_data.AMRd.g, aes(color = test_category3)) +
  geom_segment(aes(x = log(size_min), xend = log(size_max), y = intercept + slope*log(size_min), yend = intercept + slope*log(size_max)) , show.legend = FALSE )+
  ylim(x = -7, 12 )+
  xlim(x = -7, 12)+
  scale_color_manual(values = c("black", cols.amr[2]))
  # facet_grid(.~test_category3)
ggformat(ecol_amr_sum1, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*MMR~(mg~O[2]~h^-1)), print = F)
ggsave(filename = paste("./Figures/FigVAR_ecolAMR_v2_", Sys.Date(), ".png", sep=""),
       plot=ecol_amr_sum1, width = 4, height = 4, units = "in")

# Species specific groupings: 
species_amr_sum1<-ggplot(species.amr, aes(color = test_category3)) +
  geom_segment(aes(x = log(size_min), xend = log(size_max), y = intercept + slope*log(size_min), yend = intercept + slope*log(size_max)) , show.legend = FALSE )+
  ylim(x = -7, 12 )+
  xlim(x = -7, 12)+
  scale_color_manual(values = c("black", cols.amr[2]))
  # facet_grid(.~test_category3)
ggformat(species_amr_sum1, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*MMR~(mg~O[2]~h^-1)), print = F)
ggsave(filename = paste("./Figures/FigVAR_speciesAMR_v2_", Sys.Date(), ".png", sep=""),
       plot=species_amr_sum1, width = 4, height = 4, units = "in")

species_rmr_sum1<-ggplot(species.rmr, aes(color = test_category3)) +
  geom_segment(aes(x = log(size_min), xend = log(size_max), y = intercept + slope*log(size_min), yend = intercept + slope*log(size_max)) , show.legend = FALSE )+
  ylim(x = -7, 12 )+
  xlim(x = -7, 12)+
  scale_color_manual(values = c("black", cols.rmr[2]))
  # facet_grid(.~test_category3)
ggformat(species_rmr_sum1, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*RMR~(mg~O[2]~h^-1)), print = F )
ggsave(filename = paste("./Figures/FigVAR_speciesRMR_v2_", Sys.Date(), ".png", sep=""),
       plot=species_rmr_sum1, width = 4, height = 4, units = "in")


  
# Ecologies -------
plot_hist2_amrlm<-ggplot(ecology_data.AMRd, aes(x=slope)) +
  geom_vline(xintercept = mean(ecology_data.AMRd_er$slope), color="black", lty="solid", lwd=1)+
  geom_vline(xintercept = mean(ecology_data.AMRd_w$slope), color="red", lty="solid", lwd=1)+
  geom_text(mapping = aes(x= 1.3, y = 4.6), label = paste(round(mean(ecology_data.AMRd_er$slope),2), " (", round(sd(ecology_data.AMRd_er$slope),2) ,")", sep =""), family = "Arial", size=6, fontface="bold", color="black")+
  geom_text(mapping = aes(x= 1.3, y = 4), label = paste(round(mean(ecology_data.AMRd_w$slope),2), " (", round(sd(ecology_data.AMRd_w$slope),2) ,")", sep =""), family = "Arial", size=6, fontface="bold", color="red")+
  geom_histogram(show.legend = FALSE, bins=30, fill = cols.amr[1] , alpha = 1)+
  scale_color_manual(values = c("black", "red"))+
  ylim(0,5)+
  scale_x_continuous(limits=c(0.5, 1.5))+
  facet_wrap(.~test_category3)
ggformat(plot_hist2_amrlm, x_title = expression(Scaling~slope~(italic(b))), y_title = "Frequency", print = FALSE)
plot_hist2_amrlm <- plot_hist2_amrlm + theme(axis.title.x = element_text(size=25,  family="Arial"),
                                             axis.text.x = element_text(size=21,  family="Arial"),
                                             axis.title.y = element_text(size=25,  family="Arial"),
                                             axis.text.y = element_text(size=21,  family="Arial"))

plot_hist2_rmrlm<-ggplot(ecology_data.RMRd, aes(x=slope)) +
  geom_vline(xintercept = mean(ecology_data.RMRd_er$slope), color="black", lty="solid", lwd=1)+
  geom_vline(xintercept = mean(ecology_data.RMRd_w$slope), color="red", lty="solid", lwd=1)+
  geom_text(mapping = aes(x= 1.3, y = 4.6), label = paste(round(mean(ecology_data.RMRd_er$slope),2), " (", round(sd(ecology_data.RMRd_er$slope),2) ,")", sep =""), family = "Arial", size=6, fontface="bold", color="black")+
  geom_text(mapping = aes(x= 1.3, y = 4), label = paste(round(mean(ecology_data.RMRd_w$slope),2), " (", round(sd(ecology_data.RMRd_w$slope),2) ,")", sep =""), family = "Arial", size=6, fontface="bold", color="red")+
  geom_histogram(show.legend = FALSE, bins=30, fill = cols.rmr[1] , alpha = 1)+
  scale_color_manual(values = c("black", "red"))+
  ylim(0,5)+
  scale_x_continuous(limits=c(0.4, 1.5))+
  facet_wrap(.~test_category3)
ggformat(plot_hist2_rmrlm, x_title = expression(Scaling~slope~(italic(b))), y_title = "Frequency", print = FALSE)
plot_hist2_rmrlm <- plot_hist2_rmrlm + theme(axis.title.x = element_text(size=25,  family="Arial"),
                                             axis.text.x = element_text(size=21,  family="Arial"),
                                             axis.title.y = element_text(size=25,  family="Arial"),
                                             axis.text.y = element_text(size=21,  family="Arial"))


plot_hist2_faslm<-ggplot(ecology_data.FASd, aes(x=slope)) +
  geom_vline(xintercept = mean(ecology_data.FASd_er$slope), color="black", lty="solid", lwd=1)+
  geom_vline(xintercept = 0, color="grey30", lty="dashed", lwd=0.5)+
  geom_vline(xintercept = mean(ecology_data.FASd_w$slope), color="red", lty="solid", lwd=1)+
  geom_text(mapping = aes(x= -0.5, y = 5.5), label = paste(round(mean(ecology_data.FASd_er$slope),2), " (", round(sd(ecology_data.FASd_er$slope),2) ,")", sep =""), family = "Arial", size=6, fontface="bold", color="black")+
  geom_text(mapping = aes(x= -0.5, y = 4.5), label = paste(round(mean(ecology_data.FASd_w$slope),2), " (", round(sd(ecology_data.FASd_w$slope),2) ,")", sep =""), family = "Arial", size=6, fontface="bold", color="red")+
  geom_histogram(show.legend = FALSE, bins=40, fill = cols.fas[1], alpha = 1)+
  scale_color_manual(values = c("black", "red"))+
  # ylim(0,12)+
  # scale_x_continuous(limits=c(-0.1, 1.7))+
  facet_wrap(.~test_category3)
ggformat(plot_hist2_faslm, x_title = expression(Scaling~slope~(italic(b))), y_title = "Frequency", print = FALSE)
plot_hist2_faslm <- plot_hist2_faslm + theme(axis.title.x = element_text(size=25,  family="Arial"),
                                             axis.text.x = element_text(size=21,  family="Arial"),
                                             axis.title.y = element_text(size=25,  family="Arial"),
                                             axis.text.y = element_text(size=21,  family="Arial"))


ecolhist.plot<-cowplot:::plot_grid(plot_hist2_amrlm, plot_hist2_rmrlm,plot_hist2_faslm,
                                   align = "hv",
                                   axis = "l",
                                   nrow = 3,
                                   ncol = 1) 
save_plot(filename = paste("./Figures/Suppl_ecol_hist", Sys.Date(), ".png", sep=""),
          ecolhist.plot, base_width = 8, base_height = 12, units = "in")


rmr.WARM.ER.good<-ggplot(data=ecology_data.RMRd.g,
                         aes(x=slope,
                             y=factor(ecol_temp_cat, level = c(as.character(ecol_temp_cat))),
                             color=interaction(test_category3, MR),
                             fill=interaction(test_category3, MR), 
                             group= ecology_subgroup))+
  scale_y_discrete(limits=rev,
    labels = c("",
               "Any salinity", 
               "",
               "Marine",
               "",
               "Reef-associated",
              "",
               "Demersal",
              "",
               "Benthopelagic", 
               "",
               "Tropical",
               "",
               "Temperate",
               "",
               "Suptropical",
               "",
               "Fusiform"))+
  geom_text(aes(label = slope,  x = 0.25),
            family = "Arial", fontface ="bold", size=3.5, color = cols.rmr[1], hjust = 1)+
  geom_text(data = ecology_data.AMRd.g, aes(label = slope,  x = 0.3),
            family = "Arial", fontface ="bold", size=3.5, color = cols.amr[1], hjust = 0)+
  
  geom_text(data = ecology_data.RMRd.g,
            mapping = aes(label =  n_data_n_species, x = 1.5),
            family = "Arial", size=3.5, color = cols.rmr[1], hjust = 1)+
  geom_text(data = ecology_data.AMRd.g,
            mapping = aes(label =  n_data_n_species, x = 1.55),
            family = "Arial", size=3.5, color = cols.amr[1], hjust = 0)+
  
  geom_linerange(data = ecology_data.RMRd.g,
                 aes(xmin = as.numeric(slope.ciL),
                     xmax = as.numeric(slope.ciH)), alpha = 0.4)+
  geom_vline(xintercept = 1, lty = "dotted", color="grey")+
  geom_linerange(aes(xmin = slope.ciL, xmax = slope.ciH),
                 alpha = 0.4)+
  geom_linerange(data = ecology_data.AMRd.g,
                 mapping = aes(xmin = as.numeric(slope.ciL), xmax = as.numeric(slope.ciH)),
                 alpha = 0.4)+
  geom_point(pch=21, size=4, stroke=0.6, alpha=1)+
  geom_point(data = ecology_data.AMRd.g,
             pch=21, size=4, stroke=0.6, alpha=1)+
  geom_line(arrow = arrow(length=unit(0.20,"cm"), ends="last", type = "closed"),
            size = 0.7, color = cols.rmr[1])+
  geom_line(data = ecology_data.AMRd.g,
            arrow = arrow(length=unit(0.20,"cm"), ends="last", type = "closed"),
            size = 0.7, color = cols.amr[1])+
  scale_color_manual(values=c("black", "black",cols.amr[1], cols.rmr[1]))+
  scale_fill_manual(values=c(cols.amr[1],cols.rmr[1],cols.amr[2], cols.rmr[2]))+
  xlab(expression(Slope~value~(italic(b))))+
  scale_x_continuous(limits = c(0.1,1.8))+
  # scale_y_discrete(limits=rev)+
  theme_classic()+
  theme(axis.text.y = element_text(face = "italic", color = "black", size = 15),
        axis.text.x = element_text( color = "black", size = 15),
        axis.title.y = element_blank(),
        legend.position = "none",
        axis.line.y=element_line(colour = 'black',size=0.5),
        axis.line.x=element_line(colour = 'black',size=0.5),
        axis.ticks.y=element_line(size=0.5),
        axis.ticks.x=element_line(size=0),
        text=element_text(size=20,  family="Arial"))
rmr.WARM.ER.good


fas.WARM.ER.good<-ggplot(data=ecology_data.FASd.g,
                         aes(x=slope,
                             y=factor(ecol_temp_cat, level = c(as.character(ecol_temp_cat))),
                             color=test_category3,
                             fill=test_category3, 
                             group= ecology_subgroup))+
  scale_y_discrete(limits=rev,
    labels = c("",
               "Any salinity", 
               "",
               "Marine",
               "",
               "Reef-associated",
              "",
               "Demersal",
              "",
               "Benthopelagic", 
               "",
               "Tropical",
               "",
               "Temperate",
               "",
               "Suptropical",
               "",
               "Fusiform"))+
  geom_text(aes(label = n_data_n_species,  x = 0.2), family = "Arial", size=3.5, hjust=0)+
  geom_text(aes(label = slope, x = -0.2, family = "Arial"), size=3.5, hjust=1, fontface="bold")+
  geom_linerange(data = ecology_data.FASd.g,
                 aes(xmin = as.numeric(slope.ciL) , xmax = as.numeric(slope.ciH)), alpha = 0.4)+
  geom_vline(xintercept = 0, lty = "dotted", color="grey")+
  geom_linerange(aes(xmin = slope.ciL, xmax = slope.ciH), alpha = 0.4)+
  geom_point(pch=21, size=4, stroke=0.6, alpha=1)+
  geom_line(arrow = arrow(length=unit(0.20,"cm"), ends="first", type = "closed"), size = 0.7)+
  scale_color_manual(values=c("black", cols.fas[1]))+
  scale_fill_manual(values=c("black", cols.fas[2]))+
  xlab(expression(Slope~value~(italic(b))))+
  scale_x_continuous(limits = c(-0.3,0.4))+
  # scale_y_discrete(limits=rev)+
  theme_classic()+
  theme(
    # axis.text.y = element_text(face = "italic", color = "black", size = 15),
        axis.text.x = element_text( color = "black", size = 15),
        axis.title.y = element_blank(),
        legend.position = "none",
        axis.line.y=element_line(colour = 'black',size=0.5),
        axis.line.x=element_line(colour = 'black',size=0.5),
        axis.ticks.y=element_line(size=0.5),
        axis.text.y = element_blank(),
        axis.ticks.x=element_line(size=0),
        text=element_text(size=20,  family="Arial"))
fas.WARM.ER.good

cowplot::plot_grid(rmr.WARM.ER.good,
                   fas.WARM.ER.good,
                  nrow = 1, 
                  labels = "AUTO", 
                  rel_widths = c(1, 0.5),
                  label_x = c(0.435, 0.25),
                  label_y = c(0.98, 0.98)) %>% 
  ggsave(filename = paste("./Figures/Fig_ecology1_ScalingSuited_FAS", Sys.Date(), ".png", sep=""),
         width = 9.5, height = 6)


# Violin plots: AS and FAS ecology size independent -------------

# Demersal Pelagic
ecolAS1<-ggplot(data=data.as, aes(y=mass_specas, x=DemersPelag, size=BW_g, fill=test_category3, group = interaction(test_category3, DemersPelag), color=tempTest))+
  geom_point(position=position_jitterdodge(jitter.width = 0.3, seed = 1), pch=19, alpha=0.6, show.legend = FALSE)+
  geom_violin(show.legend = FALSE, alpha=0.4)+
  scale_fill_manual(values = c("black", "red"))+
  scale_color_gradient2(low="#004490", high="#C34264", mid = "#FFB654", midpoint = 15)+
  scale_size(breaks=c(10, 1000,75000), range=c(1,8))+
  scale_x_discrete(labels=c("Bentho- \n pelagic", "Demersal", "Pelagic", "Reef- \n associated"))+
  ylim(0,5)+
  labs(fill = expression(degree*C), color=  expression(degree*C), size = "g")+
  annotate("text", label = "Lifestyle", x= 3.25, y=4.8, hjust=-0.5, size=5)+
  annotate("text", label = "*** (ecol subgr)",  x= 0.6, y=4.8, hjust=-0.5, size=3.5)
ggformat(ecolAS1, y_title = bquote("AS" ~ (mgO[2] ~ g^-1 ~ h^-1)), x_title = element_blank() , print=FALSE)
ecolAS1 <- ecolAS1 + theme(
  axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5,  size = 15, family = "Arial"), 
  axis.text.y = element_text( color = "black", size = 15), 
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  legend.position = "right",
  legend.direction = "vertical",
  legend.justification='center',
  legend.margin=margin(0,0,0,0),
  legend.box.margin=margin(0,0,6,0),
  plot.margin = margin(5.5, 5.5, 5.35, 5.5))
# ecolAS1

ecolFAS1<-ggplot(data=data.fas, aes(y=FAS, x=DemersPelag, size=BW_g, fill=test_category3, group = interaction(test_category3, DemersPelag), color=tempTest))+
  geom_point(position=position_jitterdodge(jitter.width = 0.3, seed = 1), pch=19, alpha=0.6, show.legend = FALSE)+
  geom_violin(show.legend = FALSE, alpha=0.4)+
  scale_fill_manual(values = c("black", "red"))+
  scale_color_gradient2(low="#004490", high="#C34264", mid = "#FFB654", midpoint = 15)+
  scale_size(breaks=c(10, 1000,75000), range=c(1,8))+
  scale_x_discrete(labels=c("Bentho- \n pelagic", "Demersal", "Pelagic", "Reef- \n associated"))+ 
  annotate("text", x = 3.5, y = 16, label = "Lifestyle", hjust=0, family="Arial", size=5)+
  ylim(0,16.5)+
  labs(fill = expression(degree*C), color=  expression(degree*C), size = "g")
ggformat(ecolFAS1, y_title = expression(FAS~(MMR/RMR)), x_title = element_blank() , print=FALSE)
ecolFAS1 <- ecolFAS1 + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5,  size = 15, family = "Arial"), 
                             axis.text.y = element_text( color = "black", size = 15), 
                             text = element_text( color = "black", size = 15),
                             legend.position = "right",
                             legend.direction = "vertical",
                             legend.justification='center',
                             legend.margin=margin(0,0,0,0),
                             legend.box.margin=margin(0,0,6,0),
                             plot.margin = margin(5.5, 5.5, 5.35, 5.5))
# ecolFAS1

# Body shape
ecolAS2<-ggplot(data=data.as, aes(y=mass_specas, x=BodyShapeI, size=BW_g, fill=test_category3, group = interaction(test_category3, BodyShapeI), color=tempTest))+
  geom_point(position=position_jitterdodge(jitter.width = 0.3, seed = 1), pch=19, alpha=0.6, show.legend = TRUE)+
  geom_violin(show.legend = FALSE, alpha=0.4)+
  scale_fill_manual(values = c("black", "red"), guide = "none" )+
  scale_color_gradient2(low="#004490", high="#C34264", mid = "#FFB654", midpoint = 15)+# , guide = "none"
  scale_size(breaks=c(10, 1000,75000), range=c(1,8))+
  scale_x_discrete(labels=c("Elongated", "Fusiform", "Short- \n deep"))+
  labs(fill = expression(degree*C), color=  expression(degree*C), size = "g")+
  ylim(0,5)+
  annotate("text", label = "Morphology", x= 2.8, y=4.8, hjust=0.5, size=5)
ggformat(ecolAS2, y_title =  bquote("Aerobic scope" ~ (mgO[2] ~ g^- 1 ~ h^-1)), x_title = element_blank() , print=FALSE)
ecolAS2 <- ecolAS2 + theme(
  axis.title.x = element_blank(), 
  axis.title.y = element_blank(), 
  axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5,  size = 15, family = "Arial"),
  axis.text.y = element_text( color = "black", size = 15, family="Arial"),
  legend.position = "none",
  plot.margin = margin(unit(c(5.5, 100,5.5,5.5), "points")))
# ecolAS2


# ecolAS2<- ecolAS2<-theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5, family="Arial"), 
#                          axis.text.y = element_text(face = "italic", color = "black", size = 10, family="Arial"),
#                           text=element_text(size=15,  family="Arial"), 
#                          legend.position = "none")

ecolFAS2 <- ggplot(data=data.fas, aes(y=FAS, x=BodyShapeI, size=BW_g, fill=test_category3, group = interaction(test_category3, BodyShapeI), color=tempTest))+
  geom_point(position=position_jitterdodge(jitter.width = 0.3, seed = 1), pch=19, alpha=0.6, show.legend = FALSE)+
  geom_violin(show.legend = FALSE, alpha=0.4)+
  scale_fill_manual(values = c("black", "red"))+
  scale_color_gradient2(low="#004490", high="#C34264", mid = "#FFB654", midpoint = 15)+
  scale_x_discrete(labels=c("Elongated", "Fusiform", "Short-\n deep"))+
  ylim(0,16.5)+
  annotate("text", label = "Morphology", x= 2.8, y=16, hjust=0.5, size=5)
ggformat(ecolFAS2, y_title = expression(FAS~(MMR/RMR)), x_title = element_blank() , print=FALSE)
ecolFAS2 <- ecolFAS2 + theme( axis.title.x = element_blank(), 
                              axis.title.y = element_blank(), 
                              axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5,  size = 15, family = "Arial"),
                              axis.text.y = element_text( color = "black", size = 15, family="Arial"),
                              text = element_text( color = "black", size = 15),
                              legend.position = "top",
                              legend.direction = "horizontal",
                              legend.justification='left',
                              plot.margin = margin(unit(c(5.5, 100,5.5, 5.5), "points")))
# ecolFAS2

# Climate
ecolAS3 <-ggplot(data=data.as, aes(y=mass_specas, x=Climate, size=BW_g, fill=test_category3, group = interaction(test_category3, Climate), color=tempTest))+
  geom_point(position=position_jitterdodge(jitter.width = 0.3, seed = 1), pch=19, alpha=0.6, show.legend = FALSE)+
  geom_violin(show.legend = FALSE, alpha=0.4, color = "black")+
  scale_fill_manual(values = c("black", "red"))+
  scale_color_gradient2(low="#004490", high="#C34264", mid = "#FFB654", midpoint = 15, guide="none")+
  labs(fill = expression(degree*C), color=  expression(degree*C), size = "g")+
  ylim(0,5)+
  annotate("text", label = "Climate", x= 3.8, y=4.8, hjust=0.3, size=5)
  # annotate("text", label = "***",  x= 2.5,y=4.8, hjust=-0.5, size=5)
ggformat(ecolAS3, y_title =  bquote("AS" ~ (mgO[2] ~ g^- 1 ~ h^-1)), x_title = element_blank() , print=FALSE)
ecolAS3 <- ecolAS3 + theme(
  axis.title.x = element_blank(), 
  axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5,  size = 15, family = "Arial"),
  axis.text.y = element_text( color = "black", size = 15, family="Arial"),
  text = element_text( color = "black", size = 15),
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  legend.position = "top",
  legend.direction = "horizontal",
  legend.justification='left')
# ecolAS3


ecolFAS3<- ggplot(data=data.fas, aes(y=FAS, x=Climate, size=BW_g, fill=test_category3, group = interaction(test_category3, Climate), color=tempTest))+
  geom_point(position=position_jitterdodge(jitter.width = 0.3, seed = 1), pch=19, alpha=0.6, show.legend = FALSE)+
  ylim(0,16.5)+
  geom_violin(show.legend = FALSE, alpha=0.4, color = "black")+
  scale_fill_manual(values = c("black", "red"))+
  scale_color_gradient2(low="#004490", high="#C34264", mid = "#FFB654", midpoint = 15)+
  annotate("text", label = "Climate", x= 3.8, y=16, hjust=0.3, size=5)+
  annotate("text", label = "*** (ecol subgr)",  x= 0.6, y=16.1, hjust=-0.5, size=3.5)+
  annotate("text", label = "*** (temp categ)",  x= 0.57, y=14.9, hjust=-0.5, size=3.5)
ggformat(ecolFAS3, y_title = expression(FAS~(MMR/RMR)), x_title = element_blank() , print=FALSE)
ecolFAS3 <- ecolFAS3 + theme(
  axis.title.x = element_blank(),
  axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5,  size = 15, family = "Arial"), 
  axis.text.y = element_text( color = "black", size = 15, family="Arial"),
  text = element_text( color = "black", size = 15),
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  legend.position = "top",
  legend.direction = "horizontal",
  legend.justification='left')
# ecolFAS3

# Salinity
ecolAS4 <- ggplot(data=data.as, aes(y=mass_specas, x=salintyComb, size=BW_g, fill=test_category3, group = interaction(test_category3, salintyComb), color=tempTest))+
  geom_point(position=position_jitterdodge(jitter.width = 0.3, seed = 1), pch=19, alpha=0.6, show.legend = FALSE)+
  scale_fill_manual(values = c("black", "red"))+
  scale_color_gradient2(low="#004490", high="#C34264", mid = "#FFB654", midpoint = 15)+
  scale_size(breaks=c(10, 1000,75000), range=c(1,8), guide="none")+
  geom_violin(show.legend = FALSE, alpha=0.4)+
  labs(fill = expression(degree*C), color=  expression(degree*C), size = "g")+
  scale_x_discrete(labels=c("Freshwater" , "Freshwater-\n brackish", "Marine", "Marine-\n brackish" ,"All salinities"))+
  ylim(0,5)+
  annotate("text", label = "Salinity",  x= 4.2 ,y=4.8, hjust=-0.5, size=5)
ggformat(ecolAS4, y_title =  bquote("Aerobic scope" ~ (mgO[2] ~ g^- 1 ~ h^-1)), x_title = element_blank() , print=FALSE)
ecolAS4 <- ecolAS4 + theme(
  axis.title.x = element_blank(), 
  axis.title.y = element_blank(), 
  axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5,  size = 15, family = "Arial"),
  axis.text.y = element_text( color = "black", size = 15, family="Arial"),
  text = element_text( color = "black", size = 15), 
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  legend.position = "top",
  legend.direction = "horizontal",
  legend.justification='left')
# ecolAS4


ecolFAS4 <- ggplot(data=data.fas, aes(y=FAS, x=salintyComb, size=BW_g, color=tempTest, fill=test_category3, group = interaction(test_category3, salintyComb)))+
  geom_point(position=position_jitterdodge(jitter.width = 0.3, seed = 1), pch=19, alpha=0.6, show.legend = FALSE)+
  geom_violin(show.legend = FALSE, alpha=0.4)+
  scale_fill_manual(values = c("black", "red"))+
  scale_color_gradient2(low="#004490", high="#C34264", mid = "#FFB654", midpoint = 15)+
  scale_x_discrete(labels=c("Freshwater" , "Freshwater-\n brackish", "Marine", "Marine-\n brackish" ,"All salinities"))+
  # annotate("text", label = "***",  x= 2.5,y=16, hjust=-0.5, size=5)+
  annotate("text", label = "Salinity",  x= 4.3 ,y=16, hjust=-0.5, size=5)+
  ylim(0,16.5)
ggformat(ecolFAS4, y_title = expression(FAS~(MMR/RMR)), x_title = element_blank() , print=FALSE)
ecolFAS4 <- ecolFAS4 + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5,  size = 15, family = "Arial"), 
                             axis.title.x = element_blank(), 
                             axis.title.y = element_blank(),
                             text = element_text( color = "black", size = 15),
                             plot.title = element_text(face = "bold", size=15, hjust = 0.5),
                             legend.position = "none",
                             legend.direction = "horizontal",
                             legend.justification='left')
# ecolFAS4

plot_grid(ecolAS1, ecolAS2,
          ecolAS3, ecolAS4, 
          align = "h",
          labels = c('A', 'B', 'C', "D"),
          label_size = 15,
          nrow = 2,
          ncol=2,
          label_x = c(0.145, 0.08),
          label_y = c(0.89, 0.89))%>% 
  ggsave(filename = "./Figures/Figure_AS_mar2022-NOlegend.png", width = 10, height = 8)


plot_grid(ecolFAS1, ecolFAS2,
          ecolFAS3, ecolFAS4, 
          align = "h",
          labels = c('A', "B", 'C', "D"),
          label_size = 15,
          nrow = 2,
          ncol=2,
          label_x = c(0.16, 0.08),
          label_y = c(0.89, 0.89)) %>% 
  ggsave(filename = "./Figures/Figure_FAS_mar2022-NOlegend.png", width = 10, height = 8)


# SPECIES ------

# PLOTS:
plot_hist3_amrlm<-ggplot(species.amr, aes(x=slope)) +
  # geom_vline(xintercept = mean(d_bind.rmr$MLEslope[d_bind.rmr$parameter=="lnBWg"]), colour="black", lty=2, lwd=1)+
  geom_vline(xintercept = mean(species.amr.er$slope), color="black", lty="solid", lwd=1)+
  geom_text(mapping = aes(x= 2, y = 11.5), label = paste(round(mean(species.amr.er$slope),2), " (", round(sd(species.amr.er$slope),2) ,")", sep =""), family = "Arial", size=6, fontface="bold", color="black")+
  geom_text(mapping = aes(x= 2, y = 10.1), label = paste(round(mean(species.amr.w$slope),2), " (", round(sd(species.amr.w$slope),2) ,")", sep =""), family = "Arial", size=6, fontface="bold", color="red")+
  geom_vline(xintercept = mean(species.amr.w$slope), color="red", lty="solid", lwd=1)+
  geom_histogram(show.legend = FALSE, bins=50, fill = cols.amr[1], alpha = 1 )+
  scale_color_manual(values = c("black", "red"))+
  ylim(0,12)+
  scale_x_continuous(limits=c(-1, 3))+
  facet_wrap(.~test_category3)
ggformat(plot_hist3_amrlm, x_title = expression(Scaling~slope~(italic(b))), y_title = "Frequency", print = FALSE)
plot_hist3_amrlm <- plot_hist3_amrlm + theme(axis.title.x = element_text(size=25,  family="Arial"),
                                             axis.text.x = element_text(size=21,  family="Arial"),
                                             axis.title.y = element_text(size=25,  family="Arial"),
                                             axis.text.y = element_text(size=21,  family="Arial"))


plot_hist3_rmrlm<-ggplot(species.rmr, aes(x=slope)) +
  # geom_vline(xintercept = mean(d_bind.rmr$MLEslope[d_bind.rmr$parameter=="lnBWg"]), colour="black", lty=2, lwd=1)+
  geom_vline(xintercept = mean(species.rmr.er$slope), color="black", lty="solid", lwd=1)+
  geom_vline(xintercept = mean(species.rmr.w$slope), color="red", lty="solid", lwd=1)+
  geom_text(mapping = aes(x= 1.7, y = 11.5), label = paste(round(mean(species.rmr.er$slope),2), " (", round(sd(species.rmr.er$slope),2) ,")", sep =""), family = "Arial", size=6, fontface="bold", color="black")+
  geom_text(mapping = aes(x= 1.7, y = 10.1), label = paste(round(mean(species.rmr.w$slope),2), " (", round(sd(species.rmr.w$slope),2) ,")", sep =""), family = "Arial", size=6, fontface="bold", color="red")+
  geom_histogram(show.legend = FALSE, bins=50, fill = cols.rmr[1] , alpha = 1 )+
  scale_color_manual(values = c("black", "red"))+
  ylim(0,12)+
  scale_x_continuous(limits=c(-0.5, 2.5))+
  facet_wrap(.~test_category3)
ggformat(plot_hist3_rmrlm, x_title = expression(Scaling~slope~(italic(b))), y_title = "Frequency", print = FALSE)
plot_hist3_rmrlm <- plot_hist3_rmrlm + theme(axis.title.x = element_text(size=25,  family="Arial"),
                                             axis.text.x = element_text(size=21,  family="Arial"),
                                             axis.title.y = element_text(size=25,  family="Arial"),
                                             axis.text.y = element_text(size=21,  family="Arial"))


plot_hist3_faslm<-ggplot(species.fas, aes(x=slope)) +
  # geom_vline(xintercept = mean(d_bind.rmr$MLEslope[d_bind.rmr$parameter=="lnBWg"]), colour="black", lty=2, lwd=1)+
  geom_vline(xintercept = mean(species.fas.er$slope), color="black", lty="solid", lwd=1)+
  geom_vline(xintercept = 0, color="grey30", lty="dashed", lwd=0.5)+
  geom_vline(xintercept = mean(species.fas.w$slope), color="red", lty="solid", lwd=1)+
  geom_text(mapping = aes(x= 2, y = 11.5), label = paste(round(mean(species.fas.er$slope),2), " (", round(sd(species.fas.er$slope),2) ,")", sep =""), family = "Arial", size=6, fontface="bold", color="black")+
  geom_text(mapping = aes(x= 2, y = 10.1), label = paste(round(mean(species.fas.w$slope),2), " (", round(sd(species.fas.w$slope),2) ,")", sep =""), family = "Arial", size=6, fontface="bold", color="red")+
  geom_histogram(show.legend = FALSE, bins=40, fill = cols.fas[1], alpha = 1)+
  scale_color_manual(values = c("black", "red"))+
  ylim(0,12)+
  scale_x_continuous(limits=c(-2, 3.5))+
  facet_wrap(.~test_category3)
ggformat(plot_hist3_faslm, x_title = expression(Scaling~slope~(italic(b))), y_title = "Frequency", print = FALSE)
plot_hist3_faslm <- plot_hist3_faslm + theme(axis.title.x = element_text(size=25,  family="Arial"),
                                             axis.text.x = element_text(size=21,  family="Arial"),
                                             axis.title.y = element_text(size=25,  family="Arial"),
                                             axis.text.y = element_text(size=21,  family="Arial"))

spechist.plot<-cowplot:::plot_grid(plot_hist3_amrlm, plot_hist3_rmrlm,plot_hist3_faslm,
                              align = "hv",
                              axis = "l",
                              nrow = 3,
                              ncol = 1) 
save_plot(filename = paste("./Figures/Suppl_species_hist", Sys.Date(), ".png", sep=""),
      spechist.plot, base_width = 8, base_height = 12, units = "in")


sp1_amr<-ggplot(data=species.amr, aes(y=slope, x=fct_reorder(species, (slope.ciL-slope.ciH)), color = test_category3))+
  geom_hline(yintercept = 0.81, colour="black", lty="solid")+
  ylim(min(species.amr$slope.ciL)-0.1, max(species.amr$slope.ciH)+0.1)+
  geom_linerange(aes(ymin = slope, ymax = slope.ciH, color = test_category3), size=1.2, alpha=1)+
  geom_linerange(aes(ymin = slope, ymax = slope.ciL, color = test_category3), size=1.2, alpha=1)+
  geom_point(pch=1, color="black", size=1)+
  geom_text(aes(label = n, y = -3.4, color = test_category3), hjust = "inward",  family = "Arial", size=3)+
  ylab(expression(Slope~value~(italic(b))))+
  scale_color_manual(values = c("black", "red"))+
  facet_grid(.~test_category3)+
  coord_flip()+
  theme_pubr()+
  theme(axis.text.y = element_text(face = "italic", color = "black", size = 10),
        axis.text.x = element_text(size = 13),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(5.5, 5.8, 5.5, 5.5), "pt"),
        text=element_text(size=20,  family="Arial"))
ggsave(filename = paste("./Figures/FigSUP_sp_slopes_amr_", Sys.Date(), ".png", sep=""),
       plot=sp1_amr, width = 6, height = 8.8, units = "in")


sp1_rmr<-ggplot(data=species.rmr, aes(y=slope, x=fct_reorder(species, (slope.ciL-slope.ciH)), color = test_category3))+
  geom_hline(yintercept = 0.81, colour="black", lty="solid")+
  ylim(min(species.rmr$slope.ciL)-0.1, max(species.rmr$slope.ciH)+0.1)+
  geom_linerange(aes(ymin = slope, ymax = slope.ciH, color = test_category3), size=1.2, alpha=1)+
  geom_linerange(aes(ymin = slope, ymax = slope.ciL, color = test_category3), size=1.2, alpha=1)+
  geom_point(pch=1, color="black", size=1)+
  geom_text(aes(label = n, y = -4, color = test_category3), hjust = "inward",  family = "Arial", size=3)+
  ylab(expression(Slope~value~(italic(b))))+
  scale_color_manual(values = c("black", "red"))+
  facet_grid(.~test_category3)+
  coord_flip()+
  theme_pubr()+
  theme(axis.text.y = element_text(face = "italic", color = "black", size = 10),
        axis.text.x = element_text(size = 13),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(5.5, 5.8, 5.5, 5.5), "pt"),
        text=element_text(size=20,  family="Arial"))
ggsave(filename = paste("./Figures/FigSUP_sp_slopes_rmr_", Sys.Date(), ".png", sep=""),
       plot=sp1_rmr, width = 8, height = 10, units = "in")


sp1_fas<-ggplot(data=species.fas, aes(y=slope, x=fct_reorder(species, (slope.ciL-slope.ciH)), color = test_category3))+
  geom_hline(yintercept = 0, colour="black", lty="solid")+
  ylim(min(species.fas$slope.ciL)-0.1, max(species.fas$slope.ciH)+0.1)+
  geom_linerange(aes(ymin = slope, ymax = slope.ciH, color = test_category3), size=1.2, alpha=1)+
  geom_linerange(aes(ymin = slope, ymax = slope.ciL, color = test_category3), size=1.2, alpha=1)+
  geom_point(pch=1, color="black", size=1)+
  geom_text(aes(label = n, y = -5, color = test_category3), hjust = "inward",  family = "Arial", size=3)+
  ylab(expression(Slope~value~(italic(b))))+
  scale_color_manual(values = c("black", "red"))+
  facet_grid(.~test_category3)+
  coord_flip()+
  theme_pubr()+
  theme(axis.text.y = element_text(face = "italic", color = "black", size = 10),
        axis.text.x = element_text(size = 13),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(5.5, 5.8, 5.5, 5.5), "pt"),
        text=element_text(size=20,  family="Arial"))
ggsave(filename = paste("./Figures/FigSUP_sp_slopes_fas_", Sys.Date(), ".png", sep=""),
       plot=sp1_fas, width = 6, height = 8, units = "in")


  
