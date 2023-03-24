

k<-(8.62*10^(-5)) # Boltzmann's constant
E<-0.63 # activation energy MTE
cols.as<-c("#265F73", "#007E66", "#00C5A3")
cols.fas<-c("#395200", "#89A000", "yellow")
cols.rmr<-c("#C70039", "#FF6D7C", "#FFA3AC")
cols.amr<-c("#00749F","#00A8D6", "#9CE9FF")
cols<-c("#00749F","#C70039","#00A8D6","#FF6D7C", "#9CE9FF","#FFA3AC", "#00C5A3", "#265F73")# AMR -rmr- AMR dark - rmr dark - light - as - fas
# ******************************************************************************************************************************************************


library(reshape2)
library(here)
library(tidyverse)
library(dplyr)
library(emmeans)
options(dplyr.summarise.inform = FALSE)

# for emmeans library
# Ecologically Relevant conditions: 
emm_options(pbkrtest.limit = 7000)
emm_options(lmerTest.limit = 7000)

# ******************************************************************************************************************************************************
# ******************************************************************************************************************************************************
# GLOBAL DATA (with mass specific values) ---------

# 1) for phylo mixed models
source("./R/phylo_mixed_model.R")
# 2) for non-phylo mixed models 
# source("./R/nonPhylo_mixed_models.R")

data.amrER$tempTestK1000_inC<-((1000/data.amrER$tempTestK1000))-273.15
data.amr$tempTestK1<-1/data.amr$tempTestK
data.amr$tempTestK1000<-1000/data.amr$tempTestK

# ******************************************************************************************************************************************************
# ******************************************************************************************************************************************************
# Global models with ecology --------
# standardized metric for contrasts
# tempTest20C<-1000/(((20+273.15)))

ecol.model.update<-function(ecol.model.null, 
                            ecol.data.subset, 
                            temp.test.category, 
                            mr.type,
                            ecol.predictor, 
                            data.BIC = NULL,
                            data.ANOVA = NULL,
                            data.EMMEANS = NULL, 
                            data.CONTRASTS = NULL, 
                            data.PARAMS = NULL,
                            ref.tempTest20C, 
                            ref.lnBWg, 
                            phylo = TRUE){
  
  # if(phylo){
  #   foldername.phylo<-"Phylo"
  # }else{
  #   foldername.phylo<-"nonPhylo"
  # }

  for (i in 1:length(ecol.predictor)){
      
    if(phylo){
      foldername.phylo<-"Phylo"
      ecol.model<-Almer(update.formula(formula(ecol.model.null), paste(' ~ . + ', ecol.predictor[i], collapse = "")), data=ecol.data.subset, REML = F)
    }else{
      foldername.phylo<-"nonPhylo"
      ecol.model<-lmer(update.formula(formula(ecol.model.null), paste(' ~ . + ', ecol.predictor[i], collapse = "")), data=ecol.data.subset, REML = F)
    }
    # print(ecol.model)# ecol.model<-Almer(lnRMR ~ lnBWg + tempTestK1000 + DemersPelag + (1 | species) + (1 | species:trial), data=data.rmr.test, REML = F)
    
    BIC.table<-BICdelta(BIC(ecol.model.null, ecol.model))
    BIC.table$models<-rownames(BIC.table)
    BIC.table$var<-mr.type
    BIC.table$test.categ<-temp.test.category
    BIC.table$ecol<-ecol.predictor[i]
    
    params.table<-as.data.frame(fixef(ecol.model.null))
    ecol.params.table<-as.data.frame(fixef(ecol.model))
    params.table<-as.data.frame(cbind(params.table[1:3,1], ecol.params.table[1:3,1]))
    colnames(params.table)<-c("fixef.model.null", "fixef.model.ecol")
    params.table$estimate<-rownames(params.table)

    params.table$var<-mr.type
    params.table$test.categ<-temp.test.category
    params.table$ecol<-ecol.predictor[i]
    
    if(is.null(data.BIC)){
      data.BIC<-as.data.frame(BIC.table)
    }else{
      data.BIC<-rbind(data.BIC, BIC.table)
    }
    
    if(is.null(data.PARAMS)){
      data.PARAMS<-as.data.frame(params.table)
    }else{
      data.PARAMS<-rbind(data.PARAMS, params.table)
    }
   
    anova.comp<-as.data.frame(car::Anova(ecol.model))
    anova.comp$predictor<-rownames(anova.comp)
    anova.comp$var<-mr.type
    anova.comp$test.categ<-temp.test.category
    anova.comp$ecol<-ecol.predictor[i]
    
    if(is.null(data.ANOVA)){
      data.ANOVA<-as.data.frame(anova.comp)
    }else{
      data.ANOVA<-rbind(data.ANOVA, anova.comp)
    }
      
    if(ecol.predictor[i] == "DemersPelag"){
      if(mr.type == "FAS"){
        ecol.ref.grid<-(ref_grid(ecol.model, at = list(tempTest = ref.tempTest20C, lnBWg = ref.lnBWg), data = ecol.data.subset))
      }else{
        ecol.ref.grid<-(ref_grid(ecol.model, at = list(tempTestK1000 = ref.tempTest20C, lnBWg = ref.lnBWg)))
      }
      emm.emm<-as.data.frame(emmeans(ecol.ref.grid, pairwise ~ DemersPelag)$emmeans)
      colnames(emm.emm)[1]<-c("EcolGroup")
      emm.cont<-as.data.frame(emmeans(ecol.ref.grid, pairwise ~ DemersPelag)$contrasts)
    }
    if(ecol.predictor[i] == "BodyShapeI"){
      if(mr.type == "FAS"){
        ecol.ref.grid<-(ref_grid(ecol.model, at = list(tempTest = ref.tempTest20C, lnBWg = ref.lnBWg), data = ecol.data.subset))
      }else{
        ecol.ref.grid<-(ref_grid(ecol.model, at = list(tempTestK1000 = ref.tempTest20C, lnBWg = ref.lnBWg)))
      }
      emm.emm<-as.data.frame(emmeans(ecol.ref.grid, pairwise ~ BodyShapeI)$emmeans)
      colnames(emm.emm)[1]<-c("EcolGroup")
      emm.cont<-as.data.frame(emmeans(ecol.ref.grid, pairwise ~ BodyShapeI)$contrasts)
    }
    if(ecol.predictor[i] == "Climate"){
      if(mr.type == "FAS"){
        ecol.ref.grid<-(ref_grid(ecol.model, at = list(tempTest = ref.tempTest20C, lnBWg = ref.lnBWg), data = ecol.data.subset))
      }else{
        ecol.ref.grid<-(ref_grid(ecol.model, at = list(tempTestK1000 = ref.tempTest20C, lnBWg = ref.lnBWg)))
      }
      emm.emm<-as.data.frame(emmeans(ecol.ref.grid, pairwise ~ Climate)$emmeans)
      colnames(emm.emm)[1]<-c("EcolGroup")
      emm.cont<-as.data.frame(emmeans(ecol.ref.grid, pairwise ~ Climate)$contrasts)
    }
    if(ecol.predictor[i] == "salintyComb"){
      if(mr.type == "FAS"){
        ecol.ref.grid<-(ref_grid(ecol.model, at = list(tempTest = ref.tempTest20C, lnBWg = ref.lnBWg), data = ecol.data.subset))
      }else{
        ecol.ref.grid<-(ref_grid(ecol.model, at = list(tempTestK1000 = ref.tempTest20C, lnBWg = ref.lnBWg)))
      }
      emm.emm<-as.data.frame(emmeans(ecol.ref.grid, pairwise ~ salintyComb)$emmeans)
      colnames(emm.emm)[1]<-c("EcolGroup")
      emm.cont<-as.data.frame(emmeans(ecol.ref.grid, pairwise ~ salintyComb)$contrasts)
    }
    
     emm.emm$var<-mr.type
     emm.emm$test.categ<-temp.test.category
     emm.emm$ecol<-ecol.predictor[i]
     emm.cont$var<-mr.type
     emm.cont$test.categ<-temp.test.category
     emm.cont$ecol<-ecol.predictor[i]
      
    if(is.null(data.EMMEANS)){
      data.EMMEANS<-as.data.frame(emm.emm)
    }else{
      data.EMMEANS<-rbind(data.EMMEANS, emm.emm)
    }
    
    if(is.null(data.CONTRASTS)){
      data.CONTRASTS<-as.data.frame(emm.cont)
    }else{
      data.CONTRASTS<-rbind(data.CONTRASTS, emm.cont)
    }
  }
  
  return(list(data.BIC , data.ANOVA, data.EMMEANS, data.CONTRASTS, data.PARAMS))
}



# depending on local R setting this may give some errors about displaying all the data. 
output.RMR.ER <- ecol.model.update(
                  ecol.model.null = rmr_mod_ER, 
                  ecol.data.subset = data.rmrER,
                  temp.test.category = "opt",
                  mr.type = "RMR",
                  ecol.predictor = c("DemersPelag", "BodyShapeI", "Climate", "salintyComb"),
                  ref.tempTest20C = c(1000/(((20+273.15)))), 
                  ref.lnBWg = c(log(1)))

output.RMR.W <- ecol.model.update(
                  ecol.model.null = rmr_mod_W, 
                  ecol.data.subset = data.rmr.test,
                  temp.test.category = "warm",
                  mr.type = "RMR",
                  ecol.predictor = c("DemersPelag", "BodyShapeI", "Climate", "salintyComb"),
                  ref.tempTest20C = c(1000/(((20+273.15)))), 
                  ref.lnBWg = c(log(1)))

output.MMR.ER <- ecol.model.update(
                  ecol.model.null = amr_mod_ER, 
                  ecol.data.subset = data.amrER,
                  temp.test.category = "opt",
                  mr.type = "MMR",
                  ecol.predictor = c("DemersPelag", "BodyShapeI", "Climate", "salintyComb"),
                  ref.tempTest20C = c(1000/(((20+273.15)))), 
                  ref.lnBWg = c(log(1)))

output.MMR.W <- ecol.model.update(
                  ecol.model.null = amr_mod_W, 
                  ecol.data.subset = data.amr.test,
                  temp.test.category = "warm",
                  mr.type = "MMR",
                  ecol.predictor = c("DemersPelag", "BodyShapeI", "Climate", "salintyComb"),
                  ref.tempTest20C = c(1000/(((20+273.15)))), 
                  ref.lnBWg = c(log(1)))

output.FAS.ER <- ecol.model.update(
                  ecol.model.null = fas_mod_ER, 
                  ecol.data.subset = data.fasER,
                  temp.test.category = "opt",
                  mr.type = "FAS",
                  ecol.predictor = c("DemersPelag", "BodyShapeI", "Climate", "salintyComb"), 
                  ref.lnBWg = c(log(1)),
                  ref.tempTest20C = c(20))

output.FAS.W <- ecol.model.update(
                  ecol.model.null = fas_mod_W, 
                  ecol.data.subset = data.fas.test,
                  temp.test.category = "warm",
                  mr.type = "FAS",
                  ecol.predictor = c("DemersPelag", "BodyShapeI", "Climate", "salintyComb"),
                  ref.lnBWg = c(log(1)),
                  ref.tempTest20C = c(20))

output.AS.ER <- ecol.model.update(
                  ecol.model.null = as_mod_ER, 
                  ecol.data.subset = data.asER,
                  temp.test.category = "opt",
                  mr.type = "AS",
                  ecol.predictor = c("DemersPelag", "BodyShapeI", "Climate", "salintyComb"),
                  ref.tempTest20C = c(1000/(((20+273.15)))), 
                  ref.lnBWg = c(log(1)))

output.AS.W <- ecol.model.update(
                  ecol.model.null = as_mod_W, 
                  ecol.data.subset = data.as.test,
                  temp.test.category = "warm",
                  mr.type = "AS",
                  ecol.predictor = c("DemersPelag", "BodyShapeI", "Climate", "salintyComb"),
                  ref.tempTest20C = c(1000/(((20+273.15)))), 
                  ref.lnBWg = c(log(1)))

ecol.bic<<-rbind(output.AS.ER[[1]], output.AS.W[[1]],
      output.FAS.ER[[1]], output.FAS.W[[1]],
      output.RMR.ER[[1]], output.RMR.W[[1]],
      output.MMR.ER[[1]], output.MMR.W[[1]])

# best models 
ecol.bic[ecol.bic$delta==0,]

ecol.anova<-rbind(output.AS.ER[[2]], output.AS.W[[2]],
      output.FAS.ER[[2]], output.FAS.W[[2]],
      output.RMR.ER[[2]], output.RMR.W[[2]],
      output.MMR.ER[[2]], output.MMR.W[[2]])

# significant ones. considering alpha 0.05 level 
ecol.anova[!c(grepl("tempTest", ecol.anova$predictor) | grepl("lnBWg", ecol.anova$predictor)) & 
             ecol.anova$`Pr(>Chisq)` < 0.05,  ]

ecol.emmeans<-rbind(output.AS.ER[[3]], output.AS.W[[3]],
      output.FAS.ER[[3]], output.FAS.W[[3]],
      output.RMR.ER[[3]], output.RMR.W[[3]],
      output.MMR.ER[[3]], output.MMR.W[[3]])

ecol.posthoc.comp<-rbind(output.AS.ER[[4]], output.AS.W[[4]],
      output.FAS.ER[[4]], output.FAS.W[[4]],
      output.RMR.ER[[4]], output.RMR.W[[4]],
      output.MMR.ER[[4]], output.MMR.W[[4]])

ecol.model.params<-rbind(output.AS.ER[[5]], output.AS.W[[5]],
      output.FAS.ER[[5]], output.FAS.W[[5]],
      output.RMR.ER[[5]], output.RMR.W[[5]],
      output.MMR.ER[[5]], output.MMR.W[[5]])

write.csv(file = "./Data_exports/Ecologies/ecologies_data_emmeans.csv", ecol.emmeans, row.names=FALSE)
write.csv(file = "./Data_exports/Ecologies/ecologies_data_anova.csv", ecol.anova, row.names=FALSE)
write.csv(file = "./Data_exports/Ecologies/ecologies_data_posthoc.csv", ecol.posthoc.comp, row.names=FALSE)
write.csv(file = "./Data_exports/Ecologies/ecologies_data_bic.csv", ecol.bic, row.names=FALSE)
write.csv(file = "./Data_exports/Ecologies/ecologies_data_modelParams.csv", ecol.model.params, row.names=FALSE)


# ******************************************************************************************************************************************************
# ******************************************************************************************************************************************************
# ECOLOGY - SPECIFIC (estimate scaling) ----
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


# ******************************************************************************************************************************************************
# ******************************************************************************************************************************************************
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

# ******************************************************************************************************************************************************
# ******************************************************************************************************************************************************
# FIGURES ------------

## Ecologies: MMR and RMR trends (warm, optimal)  -------
ecol_amr_sum1<-ggplot(ecology_data.AMRd.g, aes(color = test_category3)) +
  geom_segment(aes(x = log(size_min), xend = log(size_max), y = intercept + slope*log(size_min), yend = intercept + slope*log(size_max)) , show.legend = FALSE )+
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
  geom_segment(aes(x = log(size_min), xend = log(size_max), y = intercept + slope*log(size_min), yend = intercept + slope*log(size_max)), show.legend = FALSE )+
  ylim(x = -7,12 )+
  xlim(x = -7, 12)+
  scale_color_manual(values = c("black", cols.rmr[2]))
ggformat(ecol_rmr_sum1, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*RMR~(mg~O[2]~h^-1)), print = F)
ecol_rmr_sum1<-ecol_rmr_sum1 +
  theme(plot.margin = margin(unit(c(-2.5,0,0,0), "cm")))


# Species specific groupings: 
species_amr_sum1<-ggplot(species.amr, aes(color = test_category3)) +
  geom_segment(aes(x = log(size_min), xend = log(size_max), y = intercept + slope*log(size_min), yend = intercept + slope*log(size_max)) , show.legend = FALSE )+
  ylim(x = -7, 12 )+
  xlim(x = -7, 12)+
  scale_color_manual(values = c("black", cols.amr[2]))
ggformat(species_amr_sum1, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*MMR~(mg~O[2]~h^-1)), print = F)
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
  geom_segment(aes(x = log(size_min), xend = log(size_max), y = intercept + slope*log(size_min), yend = intercept + slope*log(size_max)) , show.legend = FALSE )+
  ylim(x = -7, 12 )+
  xlim(x = -7, 12)+
  scale_color_manual(values = c("black", cols.rmr[2]))
ggformat(species_rmr_sum1, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*RMR~(mg~O[2]~h^-1)), print = F )
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

cowplot:::plot_grid(MRmodel_plot2, grid2, rel_widths = c(1, 1), 
                    labels = c("A", ""),
                    label_x = c(0.05, 0.01),
                    label_y = c(0.95, 0.9)) %>% 
  ggsave(filename = paste("./Figures/Fig5_Scaling_phylo.png", sep=""), width = 10, height = 5, units = "in")

  

 ## Ecologies: FAS trends (warm, optimal) -------
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
                  label_x = c(0.03, 0.25),
                  label_y = c(0.98, 0.98)) %>% 
  ggsave(filename = paste("./Figures/Fig4_ecology1_ScalingSuited_FAS.png", sep=""),
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
# 
#   


  


## Ecologies: violin lots, AS and FAS ecology size independent -------------
# Demersal Pelagic
ecolAS1<-ggplot(data=data.as, aes(y=mass_specas, fill=DemersPelag, x=test_category3,
                                  group = interaction(test_category3, DemersPelag)))+
  geom_violin(data=data.rmr, aes(y=mass_specrmr, fill=DemersPelag, x=test_category3,
                                group = interaction(test_category3, DemersPelag)), 
              alpha=1, size=0.2)+
  geom_violin(data=data.amr, aes(y=mass_specamr, fill=DemersPelag, x=test_category3,
                                  group = interaction(test_category3, DemersPelag)), 
              alpha=0.1, size=0.2)+  # geom_boxplot( outlier.shape = NA) +
  geom_pointrange(data.as, mapping=aes(x=test_category3, y = mass_specas,
                                       fill=DemersPelag, 
                                       group = interaction(test_category3, DemersPelag)),
                  stat = "summary",
                  position=position_dodge(width = 0.9), color = "black", pch=23)+
  scale_fill_viridis_d(option = "A")+
  scale_x_discrete(labels=c("Optimal", "Warm"))+
  ylim(0,8)+
  # geom_segment(aes(x = 0.8, y = 7.8, xend = 1.1, yend = 7.8), size = 0.1)+
  # geom_segment(aes(x = 1.9, y = 7.8, xend = 2.2, yend = 7.8), size = 0.1)+
  annotate(geom = "text", y = 7.9, x = 0.95, label = "(RMR, MMR) BIC *", size = 3)+
  annotate(geom = "text", y = 7.4, x = 0.95, label = "(AS) BIC ns", size = 3)+
  annotate(geom = "text", y = 7.9, x = 2.05, label = "(all) BIC ns", size = 3)+
  geom_vline(xintercept = 1.5, color = "grey", linetype = "dashed")
ggformat(ecolAS1, y_title = bquote("MR" ~ (mgO[2] ~ g^-1 ~ h^-1)), x_title = element_blank() , print=T, size_text = 11)
ecolAS1 <- ecolAS1 + theme(
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  legend.position = "none",
  legend.title = element_blank(),
  legend.direction = "vertical",
  legend.justification='center',
  legend.margin=margin(0,0,0,0),
  legend.box.margin=margin(0,0,6,0),
  plot.margin = margin(5.5, 5.5, 5.35, 5.5))
# ecolAS1

ecolAS2<-ggplot(data=data.as, aes(y=mass_specas, fill=BodyShapeI, x=test_category3,
                                  group = interaction(test_category3, BodyShapeI)))+
  geom_violin(data=data.rmr[!c(data.rmr$BodyShapeI == "dorsoventrflattened" |
                              data.rmr$BodyShapeI == "eel-like"),], aes(y=mass_specrmr, fill=BodyShapeI, x=test_category3,
                                group = interaction(test_category3, BodyShapeI)),
              alpha=1, size=0.2)+
  geom_violin(data=data.amr, aes(y=mass_specamr, fill=BodyShapeI, x=test_category3,
                                  group = interaction(test_category3, BodyShapeI)), 
              alpha=0.1, size=0.2)+ 
  geom_pointrange(data.as, mapping=aes(x=test_category3, y = mass_specas,
                                       fill=BodyShapeI, 
                                       group = interaction(test_category3, BodyShapeI)),
                  stat = "summary", fun.ymin = min,fun.ymax = max,fun.y = mean, 
                  position=position_dodge(width = 0.9), color = "black", pch=23)+
  scale_fill_viridis_d(option = "E")+
  scale_x_discrete(labels=c("Optimal", "Warm"))+
  ylim(0,8)+
  # geom_segment(aes(x = 0.8, y = 7.8, xend = 1.1, yend = 7.8), size = 0.1)+
  # geom_segment(aes(x = 1.9, y = 7.8, xend = 2.2, yend = 7.8), size = 0.1)+
  annotate(geom = "text", y = 7.9, x = 0.95, label = "(RMR) BIC *", size = 3)+
  annotate(geom = "text", y = 7.4, x = 0.95, label = "(AS, MMR) BIC ns", size = 3)+
  annotate(geom = "text", y = 7.9, x = 2.05, label = "(all) BIC ns", size = 3)+
  geom_vline(xintercept = 1.5, color = "grey", linetype = "dashed")
ggformat(ecolAS2, y_title = bquote("MR" ~ (mgO[2] ~ g^-1 ~ h^-1)), x_title = element_blank() , print=T, size_text = 11)
ecolAS2 <- ecolAS2 + theme(
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  legend.position = "none",
  legend.title = element_blank(),
  legend.direction = "vertical",
  legend.justification='center',
  legend.margin=margin(0,0,0,0),
  legend.box.margin=margin(0,0,6,0),
  plot.margin = margin(5.5, 5.5, 5.35, 5.5))


ecolAS3<-ggplot(data=data.as, aes(y=mass_specas, fill=salintyComb, x=test_category3,
                                  group = interaction(test_category3, salintyComb)))+
  geom_violin(data=data.rmr, aes(y=mass_specrmr, fill=salintyComb, x=test_category3,
                                group = interaction(test_category3, salintyComb)), 
              alpha=1, size=0.2)+
  geom_violin(data=data.amr, aes(y=mass_specamr, fill=salintyComb, x=test_category3,
                                  group = interaction(test_category3, salintyComb)), 
              alpha=0.1, size=0.2)+  
  geom_pointrange(data.as, mapping=aes(x=test_category3, y = mass_specas,
                                       fill=salintyComb, 
                                       group = interaction(test_category3, salintyComb)),
                  stat = "summary", fun.ymin = min,fun.ymax = max,fun.y = mean, 
                  position=position_dodge(width = 0.9), color = "black", pch=23)+
  scale_fill_viridis_d(option = "D", direction = -1)+
  scale_x_discrete(labels=c("Optimal", "Warm"))+
  ylim(0,8)+
  # geom_segment(aes(x = 0.8, y = 7.8, xend = 1.1, yend = 7.8), size = 0.1)+
  # geom_segment(aes(x = 1.9, y = 7.8, xend = 2.2, yend = 7.8), size = 0.1)+
  annotate(geom = "text", y = 7.9, x = 0.95, label = "(RMR) BIC *", size = 3)+
  annotate(geom = "text", y = 7.4, x = 0.95, label = "(AS, MMR) BIC ns", size = 3)+
  annotate(geom = "text", y = 7.9, x = 2.05, label = "(all) BIC ns", size = 3)+
  geom_vline(xintercept = 1.5, color = "grey", linetype = "dashed")
ggformat(ecolAS3, y_title = bquote("MR" ~ (mgO[2] ~ g^-1 ~ h^-1)), x_title = element_blank() , print=T, size_text = 11)
ecolAS3 <- ecolAS3 + theme(
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  legend.position = "none",
  legend.title = element_blank(),
  legend.direction = "vertical",
  legend.justification='center',
  legend.margin=margin(0,0,0,0),
  legend.box.margin=margin(0,0,6,0),
  plot.margin = margin(5.5, 5.5, 5.35, 5.5))

ecolAS4<-ggplot(data=data.as, aes(y=mass_specas, fill=Climate, x=test_category3,
                                  group = interaction(test_category3, Climate)))+
  geom_violin(data=data.rmr, aes(y=mass_specrmr, fill=Climate, x=test_category3,
                                group = interaction(test_category3, Climate)), 
              alpha=1, size=0.2)+
  geom_violin(data=data.amr, aes(y=mass_specamr, fill=Climate, x=test_category3,
                                  group = interaction(test_category3, Climate)), 
              alpha=0.1, size=0.2)+  # geom_boxplot( outlier.shape = NA) +
  geom_pointrange(data.as, mapping=aes(x=test_category3, y = mass_specas,
                                       fill=Climate, 
                                       group = interaction(test_category3, Climate)),
                  stat = "summary", fun.ymin = min,fun.ymax = max,fun.y = mean, 
                  position=position_dodge(width = 0.9), color = "black", pch=23)+
  scale_fill_viridis_d(option = "C", direction = 1)+
  scale_x_discrete(labels=c("Optimal", "Warm"))+
  ylim(0,8)+
  # geom_segment(aes(x = 0.8, y = 7.8, xend = 1.1, yend = 7.8), size = 0.1)+
  # geom_segment(aes(x = 1.9, y = 7.8, xend = 2.2, yend = 7.8), size = 0.1)+
  annotate(geom = "text", y = 7.9, x = 0.95, label = "(RMR) BIC *", size = 3)+
  annotate(geom = "text", y = 7.4, x = 0.95, label = "(AS, MMR) BIC ns", size = 3)+
  annotate(geom = "text", y = 7.9, x = 2.05, label = "(all) BIC ns", size = 3)+
  geom_vline(xintercept = 1.5, color = "grey", linetype = "dashed")
ggformat(ecolAS4, y_title = bquote("MR" ~ (mgO[2] ~ g^-1 ~ h^-1)), x_title = element_blank() , print=T, size_text = 11)
ecolAS4 <- ecolAS4 + theme(
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  legend.position = "none",
  legend.title = element_blank(),
  legend.direction = "vertical",
  legend.justification='center',
  legend.margin=margin(0,0,0,0),
  legend.box.margin=margin(0,0,6,0),
  plot.margin = margin(5.5, 5.5, 5.35, 5.5))

# Demersal Pelagic
ecolFAS1<-ggplot(data=data.fas, aes(y=FAS, fill=DemersPelag, x=test_category3,
                                  group = interaction(test_category3, DemersPelag)))+
  geom_violin(alpha = 1, outlier.size = 0.2, size = 0.2)+
  scale_fill_viridis_d(option = "A", labels=c('Benthopelagic', 'Demersal',"Pelagic", "Reef" ), "Lifestyle")+
  scale_x_discrete(labels=c("Optimal", "Warm"))+
  ylim(0,20)+
  # geom_segment(aes(x = 0.8, y = 19, xend = 1.1, yend = 19), size = 0.1)+
  # geom_segment(aes(x = 1.9, y = 19, xend = 2.2, yend = 19), size = 0.1)+
  annotate(geom = "text", y = 19.5, x = 0.95, label = "BIC ns", size = 3)+
  annotate(geom = "text", y = 19.5, x = 2.05, label = "BIC ns", size = 3)+
  geom_vline(xintercept = 1.5, color = "grey", linetype = "dashed")
ggformat(ecolFAS1, y_title = bquote("FAS" ~ (mgO[2] ~ g^-1 ~ h^-1)), x_title = element_blank() , print=T, size_text = 11)
ecolFAS1 <- ecolFAS1 + theme(
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  legend.position = "right",
  legend.direction = "vertical",
  legend.justification='center',
  legend.margin=margin(0,0,0,0),
  legend.box.margin=margin(0,0,6,0),
  plot.margin = margin(5.5, 5.5, 5.35, 5.5),
  legend.spacing.y = unit(0.5, 'cm'),
  legend.text = element_text(margin = margin(l = 0, unit = c("cm"))),
  legend.key.height= unit(0.5, 'cm'),
  legend.key.width= unit(0.5, 'cm'))
# ecolFAS1

ecolFAS2<-ggplot(data=data.fas, aes(FAS, fill=BodyShapeI, x=test_category3,
                                  group = interaction(test_category3, BodyShapeI)))+
  geom_violin(alpha = 1, outlier.size = 0.2, size = 0.2)+
  scale_fill_viridis_d(option = "E", labels=c('Elongated', 'Fusiform',"Short, Deep" ), "Body Shape")+
  scale_x_discrete(labels=c("Optimal", "Warm"))+
  ylim(0,20)+
  # geom_segment(aes(x = 0.8, y = 19, xend = 1.1, yend = 19), size = 0.1)+
  # geom_segment(aes(x = 1.9, y = 19, xend = 2.2, yend = 19), size = 0.1)+
  annotate(geom = "text", y = 19.5, x = 0.95, label = "BIC ns", size = 3)+
  annotate(geom = "text", y = 19.5, x = 2.05, label = "BIC ns", size = 3)+
  geom_vline(xintercept = 1.5, color = "grey", linetype = "dashed")
ggformat(ecolFAS2, y_title = bquote("FAS" ~ (mgO[2] ~ g^-1 ~ h^-1)), x_title = element_blank() , print=T, size_text = 11)
ecolFAS2 <- ecolFAS2 + theme( 
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  legend.position = "right",
  legend.direction = "vertical",
  legend.justification='center',
  legend.margin=margin(0,0,0,0),
  legend.box.margin=margin(0,0,6,0),
  plot.margin = margin(5.5, 5.5, 5.35, 5.5),
  legend.spacing.y = unit(0.5, 'cm'),
  legend.text = element_text(margin = margin(l = 0, unit = c("cm"))),
  legend.key.height= unit(0.5, 'cm'),
  legend.key.width= unit(0.5, 'cm'))


ecolFAS3<-ggplot(data=data.fas, aes(FAS, fill=salintyComb, x=test_category3,
                                  group = interaction(test_category3, salintyComb)))+
  geom_violin(alpha = 1, outlier.size = 0.2, size = 0.2)+
  scale_fill_viridis_d(option = "D", direction = -1,
                       labels=c('Freshw.', 'Freshw.-sal',"Marine", "Marine.-fresh", "All" ), "Salinity")+
  scale_x_discrete(labels=c("Optimal", "Warm"))+
  ylim(0,20)+
  # geom_segment(aes(x = 0.8, y = 19, xend = 1.1, yend = 19), size = 0.1)+
  # geom_segment(aes(x = 1.9, y = 19, xend = 2.2, yend = 19), size = 0.1)+
  annotate(geom = "text", y = 19.5, x = 0.95, label = "BIC ns", size = 3)+
  annotate(geom = "text", y = 19.5, x = 2.05, label = "BIC ns", size = 3)+
  geom_vline(xintercept = 1.5, color = "grey", linetype = "dashed")
ggformat(ecolFAS3, y_title = bquote("FAS" ~ (mgO[2] ~ g^-1 ~ h^-1)), x_title = element_blank() , print=T, size_text = 11)
ecolFAS3 <- ecolFAS3 + theme(
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  legend.position = "right",
  legend.direction = "vertical",
  legend.justification='center',
  legend.margin=margin(0,0,0,0),
  legend.box.margin=margin(0,0,6,0),
  plot.margin = margin(5.5, 5.5, 5.35, 5.5),
  legend.spacing.y = unit(0.5, 'cm'),
  legend.text = element_text(margin = margin(l = 0, unit = c("cm"))),
  legend.key.height= unit(0.5, 'cm'),
  legend.key.width= unit(0.5, 'cm'))

ecolFAS4<-ggplot(data=data.fas, aes(FAS, fill=Climate, x=test_category3,
                                  group = interaction(test_category3, Climate)))+
  geom_violin(alpha = 1, outlier.size = 0.2, size = 0.2)+
  scale_fill_viridis_d(option = "C", direction = 1, "Climate")+
  scale_x_discrete(labels=c("Optimal", "Warm"))+
  ylim(0,20)+
  # geom_segment(aes(x = 0.8, y = 19, xend = 1.1, yend = 19), size = 0.1)+
  # geom_segment(aes(x = 1.9, y = 19, xend = 2.2, yend = 19), size = 0.1)+
  annotate(geom = "text", y = 19.5, x = 0.95, label = "BIC *", size = 3)+
  annotate(geom = "text", y = 19.5, x = 2.05, label = "BIC ns", size = 3)+
  geom_vline(xintercept = 1.5, color = "grey", linetype = "dashed")
ggformat(ecolFAS4, y_title = bquote("FAS" ~ (mgO[2] ~ g^-1 ~ h^-1)), x_title = element_blank() , print=T, size_text = 11)
ecolFAS4 <- ecolFAS4 + theme(
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  legend.position = "right",
  legend.direction = "vertical",
  legend.justification='center',
  legend.margin=margin(0,0,0,0),
  legend.box.margin=margin(0,0,6,0),
  plot.margin = margin(5.5, 5.5, 5.35, 5.5),
  legend.spacing.y = unit(0.5, 'cm'),
  legend.text = element_text(margin = margin(l = 0, unit = c("cm"))),
  legend.key.height= unit(0.5, 'cm'),
  legend.key.width= unit(0.5, 'cm'))

plot_grid(ecolAS1,ecolFAS1,
          ecolAS2, ecolFAS2,
          ecolAS3, ecolFAS3, 
          ecolAS4, ecolFAS4,
          align = "h",
          labels = c("AUTO"),
          label_size = 15,
          nrow = 4,
          ncol = 2,
          label_x = c(0.02, 0.02),
          label_y = c(0.98, 0.98), 
          rel_widths = c(1, 1.35)) %>% 
  ggsave(filename = "./Figures/Fig5_AS_FAS_violin_mar2023.png", width = 7.5, height = 10)


## Species ------

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
ggsave(filename = paste("./Figures/FigSUP_sp_slopes_amr.png", sep=""),
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
ggsave(filename = paste("./Figures/FigSUP_sp_slopes_rmr.png", sep=""),
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
ggsave(filename = paste("./Figures/FigSUP_sp_slopes_fas.png", sep=""),
       plot=sp1_fas, width = 6, height = 8, units = "in")


 ### Misc, not used -----
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
save_plot(filename = paste("./Figures/Suppl_ecol_hist.png", sep=""),
          ecolhist.plot, base_width = 8, base_height = 12, units = "in")



# species hist ---
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
save_plot(filename = paste("./Figures/Suppl_species_hist.png", sep=""),
      spechist.plot, base_width = 8, base_height = 12, units = "in")


