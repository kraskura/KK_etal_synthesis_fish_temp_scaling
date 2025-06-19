

options(dplyr.summarise.inform = FALSE)
# *********************************************************************************

k<-(8.62*10^(-5)) # Boltzmann's constant
E<-0.63 # activation energy MTE


# for emmeans library
# Ecologically Relevant conditions: 
emm_options(pbkrtest.limit = 7000)
emm_options(lmerTest.limit = 7000)

# ************************************************************************************
# ************************************************************************************* 
# GLOBAL DATA (with mass specific values) ---------

# 1) for phylo mixed models
source(here("R", "phylo_mixed_model.R")) # run first if not already done
source(here("R/colors_themes.R"))
# 2) for non-phylo mixed models 
# source("./R/nonPhylo_mixed_models.R")

# data.amrER$tempTestK1000_inC<-((1000/data.amrER$tempTestK1000))-273.15
# data.amr$tempTestK1<-1/data.amr$tempTestK
# data.amr$tempTestK1000<-1000/data.amr$tempTestK

scaling.params<-read.csv(here("./Data_exports/Phylo/scaling_parameters.csv"))

# slopes available from  function call or saved new vars
RMR_slope<-round(as.numeric(scaling.params[scaling.params$performance == "RMR" & 
                            scaling.params$temp_categ == "er", "lnBWg"]),2)
AMR_slope5<-round(as.numeric(scaling.params[scaling.params$performance == "MMR" & 
                            scaling.params$temp_categ == "er" & 
                            scaling.params$tempTest == "5"  , "lnBWg"]),2)
AMR_slope15<-round(as.numeric(scaling.params[scaling.params$performance == "MMR" & 
                            scaling.params$temp_categ == "er" & 
                            scaling.params$tempTest == "15"  , "lnBWg"]),2)
AMR_slope25<-round(as.numeric(scaling.params[scaling.params$performance == "MMR" & 
                            scaling.params$temp_categ == "er" & 
                            scaling.params$tempTest == "25"  , "lnBWg"]),2)
AMR_slope35<-round(as.numeric(scaling.params[scaling.params$performance == "MMR" & 
                            scaling.params$temp_categ == "er" & 
                            scaling.params$tempTest == "35"  , "lnBWg"]),2)
FAS_slope<-round(as.numeric(scaling.params[scaling.params$performance == "FAS" & 
                            scaling.params$temp_categ == "er", "lnBWg"]),2)
AS_slope<-round(as.numeric(scaling.params[scaling.params$performance == "AS" & 
                            scaling.params$temp_categ == "er", "lnBWg"]),2)
FAS_slope_w<-round(as.numeric(scaling.params[scaling.params$performance == "FAS" & 
                            scaling.params$temp_categ == "warm", "lnBWg"]),2)
AS_slope_w<-round(as.numeric(scaling.params[scaling.params$performance == "AS" & 
                            scaling.params$temp_categ == "warm", "lnBWg"]),2)
RMR_slope_w<-round(as.numeric(scaling.params[scaling.params$performance == "RMR" & 
                            scaling.params$temp_categ == "warm", "lnBWg"]),2)
AMR_slope_w<-round(as.numeric(scaling.params[scaling.params$performance == "MMR" & 
                            scaling.params$temp_categ == "warm", "lnBWg"]),2)


# ************************************************************************************
# ************************************************************************************* 
# Reassign the data frames with scaling parameters ----
# # reset datasets with mass specific values using scaling coefficients from the models.
data.list<-get_data_temp(data.amr = "./Data/Fish_AMR_temp_dataset_mar2022.csv",
                         data.rmr = "./Data/Fish_RMR_temp_dataset_mar2022.csv",
                         ecology.data = "./Data/Kraskura_species_ecologies_mar2022.csv",
                         onlyTop.above = TRUE,
                         calc_mass_specific = TRUE,
                         exp_rmr = RMR_slope,
                         exp_amr = AMR_slope25, # at 25ºC
                         exp_as = AS_slope,
                         exp_rmr_warm = RMR_slope_w,
                         exp_amr_warm = AMR_slope_w,
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


# ********************************************************************************************************************
# ********************************************************************************************************************
# Global models with ecology --------
# standardized metric for contrasts
# tempTest<-1000/(((20+273.15)))
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
                            ref.tempTest, 
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
      ecol.model<-Almer(update.formula(formula(ecol.model.null), paste(' ~ . + ', ecol.predictor[i], collapse = "")),
                        data=ecol.data.subset,
                        REML = F)
    }else{
      foldername.phylo<-"nonPhylo"
      ecol.model<-lmer(update.formula(formula(ecol.model.null), paste(' ~ . + ', ecol.predictor[i], collapse = "")),
                       data=ecol.data.subset,
                       REML = F)
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
   
    anova.comp<-as.data.frame(car::Anova(ecol.model), )
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
        ecol.ref.grid<-(ref_grid(ecol.model,
                                 at = list(tempTest = ref.tempTest, lnBWg = ref.lnBWg),
                                 data = ecol.data.subset))

      emm.emm<-as.data.frame(emmeans(ecol.ref.grid, pairwise ~ DemersPelag)$emmeans)
      colnames(emm.emm)[1]<-c("EcolGroup")
      emm.cont<-as.data.frame(emmeans(ecol.ref.grid, pairwise ~ DemersPelag)$contrasts)
    }
    
    if(ecol.predictor[i] == "BodyShapeI"){
        ecol.ref.grid<-(ref_grid(ecol.model,
                                 at = list(tempTest = ref.tempTest, lnBWg = ref.lnBWg),
                                 data = ecol.data.subset))
     
      emm.emm<-as.data.frame(emmeans(ecol.ref.grid, pairwise ~ BodyShapeI)$emmeans)
      colnames(emm.emm)[1]<-c("EcolGroup")
      emm.cont<-as.data.frame(emmeans(ecol.ref.grid, pairwise ~ BodyShapeI)$contrasts)
    }
    if(ecol.predictor[i] == "Climate"){
      ecol.ref.grid<-(ref_grid(ecol.model,
                                 at = list(tempTest = ref.tempTest, lnBWg = ref.lnBWg),
                                 data = ecol.data.subset))

      emm.emm<-as.data.frame(emmeans(ecol.ref.grid, pairwise ~ Climate)$emmeans)
      colnames(emm.emm)[1]<-c("EcolGroup")
      emm.cont<-as.data.frame(emmeans(ecol.ref.grid, pairwise ~ Climate)$contrasts)
    }
    
    if(ecol.predictor[i] == "salintyComb"){
        ecol.ref.grid<-(ref_grid(ecol.model,
                                 at = list(tempTest = ref.tempTest, lnBWg = ref.lnBWg),
                                 data = ecol.data.subset))

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
                  data.BIC = NULL,
                  data.ANOVA = NULL,
                  data.EMMEANS = NULL, 
                  data.CONTRASTS = NULL, 
                  data.PARAMS = NULL,
                  ecol.model.null = rmr_mod_ER, 
                  ecol.data.subset = data.rmrER,
                  temp.test.category = "opt",
                  mr.type = "RMR",
                  ecol.predictor = c("DemersPelag", "BodyShapeI", "Climate", "salintyComb"),
                  ref.tempTest = 25, 
                  ref.lnBWg = c(log(1)))

# causes error sometimes, not other times?
output.RMR.W <- ecol.model.update(
                  data.BIC = NULL,
                  data.ANOVA = NULL,
                  data.EMMEANS = NULL, 
                  data.CONTRASTS = NULL, 
                  data.PARAMS = NULL,
                  ecol.model.null = rmr_mod_W, 
                  ecol.data.subset = data.rmr.test,
                  temp.test.category = "warm",
                  mr.type = "RMR",
                  ecol.predictor = c("DemersPelag", "BodyShapeI", "Climate", "salintyComb"),
                  ref.tempTest = 25, 
                  ref.lnBWg = c(log(1)))

output.MMR.ER <- ecol.model.update(
                  data.BIC = NULL,
                  data.ANOVA = NULL,
                  data.EMMEANS = NULL, 
                  data.CONTRASTS = NULL, 
                  data.PARAMS = NULL,
                  ecol.model.null = amr_mod_ER, 
                  ecol.data.subset = data.amrER,
                  temp.test.category = "opt",
                  mr.type = "MMR",
                  ecol.predictor = c("DemersPelag", "BodyShapeI", "Climate", "salintyComb"),
                  ref.tempTest = 25, 
                  ref.lnBWg = c(log(1)))

output.MMR.W <- ecol.model.update(
                  data.BIC = NULL,
                  data.ANOVA = NULL,
                  data.EMMEANS = NULL, 
                  data.CONTRASTS = NULL, 
                  data.PARAMS = NULL,
                  ecol.model.null = amr_mod_W, 
                  ecol.data.subset = data.amr.test,
                  temp.test.category = "warm",
                  mr.type = "MMR",
                  ecol.predictor = c("DemersPelag", "BodyShapeI", "Climate", "salintyComb"),
                  ref.tempTest = 25, 
                  ref.lnBWg = c(log(1)))

output.FAS.ER <- ecol.model.update(
                  data.BIC = NULL,
                  data.ANOVA = NULL,
                  data.EMMEANS = NULL, 
                  data.CONTRASTS = NULL, 
                  data.PARAMS = NULL,
                  ecol.model.null = fas_mod_ER, 
                  ecol.data.subset = data.fasER,
                  temp.test.category = "opt",
                  mr.type = "FAS",
                  ecol.predictor = c("DemersPelag", "BodyShapeI", "Climate", "salintyComb"), 
                  ref.lnBWg = c(log(1)),
                  ref.tempTest = 25)

output.FAS.W <- ecol.model.update(
                  data.BIC = NULL,
                  data.ANOVA = NULL,
                  data.EMMEANS = NULL, 
                  data.CONTRASTS = NULL, 
                  data.PARAMS = NULL,
                  ecol.model.null = fas_mod_W, 
                  ecol.data.subset = data.fas.test,
                  temp.test.category = "warm",
                  mr.type = "FAS",
                  ecol.predictor = c("DemersPelag", "BodyShapeI", "Climate", "salintyComb"),
                  ref.lnBWg = c(log(1)),
                  ref.tempTest = 25)

output.AS.ER <- ecol.model.update(
                  data.BIC = NULL,
                  data.ANOVA = NULL,
                  data.EMMEANS = NULL, 
                  data.CONTRASTS = NULL, 
                  data.PARAMS = NULL,
                  ecol.model.null = as_mod_ER, 
                  ecol.data.subset = data.asER,
                  temp.test.category = "opt",
                  mr.type = "AS",
                  ecol.predictor = c("DemersPelag", "BodyShapeI", "Climate", "salintyComb"),
                  ref.tempTest = 25, 
                  ref.lnBWg = c(log(1)))

output.AS.W <- ecol.model.update(
                  data.BIC = NULL,
                  data.ANOVA = NULL,
                  data.EMMEANS = NULL, 
                  data.CONTRASTS = NULL, 
                  data.PARAMS = NULL,
                  ecol.model.null = as_mod_W, 
                  ecol.data.subset = data.as.test,
                  temp.test.category = "warm",
                  mr.type = "AS",
                  ecol.predictor = c("DemersPelag", "BodyShapeI", "Climate", "salintyComb"),
                  ref.tempTest = 25, 
                  ref.lnBWg = c(log(1)))


# best models -----
ecol.bic<<-rbind(output.AS.ER[[1]], output.AS.W[[1]],
      output.FAS.ER[[1]], output.FAS.W[[1]],
      output.RMR.ER[[1]], output.RMR.W[[1]],
      output.MMR.ER[[1]], output.MMR.W[[1]])

ecol.bic[ecol.bic$delta==0,]

#models where ecology IMPROVED the fit 
# ecol.model11      10  726.36883     0      ecol.model  AS       warm  BodyShapeI
# ecol.model22      10 1575.08535     0      ecol.model FAS        opt     Climate
# ecol.model7       10 2185.39445     0      ecol.model RMR        opt DemersPelag
# ecol.model14      11 2184.36040     0      ecol.model RMR        opt  BodyShapeI
# ecol.model24      10 2193.66343     0      ecol.model RMR        opt     Climate
# ecol.model34      11 2203.20611     0      ecol.model RMR        opt salintyComb
# ecol.model9       11  227.08662     0      ecol.model MMR        opt DemersPelag


# significant ones. considering alpha 0.05 level  -----
ecol.anova<-rbind(output.AS.ER[[2]], output.AS.W[[2]],
      output.FAS.ER[[2]], output.FAS.W[[2]],
      output.RMR.ER[[2]], output.RMR.W[[2]],
      output.MMR.ER[[2]], output.MMR.W[[2]])

# ecol.anova[!c(grepl("tempTest", ecol.anova$predictor) |
#                 grepl("lnBWg", ecol.anova$predictor)) &
#              ecol.anova$`Pr(>Chisq)` < 0.05,  ]

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

# save outputs -------
write.csv(file = here("./Data_exports/Ecologies/ecologies_data_emmeans.csv"), ecol.emmeans, row.names=FALSE)
write.csv(file = here("./Data_exports/Ecologies/ecologies_data_anova.csv"), ecol.anova, row.names=FALSE)
write.csv(file = here("./Data_exports/Ecologies/ecologies_data_posthoc.csv"), ecol.posthoc.comp, row.names=FALSE)
write.csv(file = here("./Data_exports/Ecologies/ecologies_data_bic.csv"), ecol.bic, row.names=FALSE)
write.csv(file = here("./Data_exports/Ecologies/ecologies_data_modelParams.csv"), ecol.model.params, row.names=FALSE)


## Ecologies: violin lots, AS and FAS ecology size independent -------------
### Demersal Pelagic MMR, RMR, AS ------
data.as$DemersPelag_plot<-paste(data.as$DemersPelag, data.as$test_category3, sep = "-")
data.amr$DemersPelag_plot<-paste(data.amr$DemersPelag, data.amr$test_category3, sep = "-")
data.rmr$DemersPelag_plot<-paste(data.rmr$DemersPelag, data.rmr$test_category3, sep = "-")
data.fas$DemersPelag_plot<-paste(data.fas$DemersPelag, data.fas$test_category3, sep = "-")

data.fas$DemersPelag_plot<-factor(data.fas$DemersPelag_plot, c(
                                  "benthopelagic-ecol_relev",
                                  "demersal-ecol_relev" ,
                                  "pelagic-ecol_relev", 
                                  "reef-associated-ecol_relev",
                                  "benthopelagic-warm",
                                  "demersal-warm" ,
                                  "pelagic-warm" ,
                                  "reef-associated-warm"
                                  ))

data.rmr$DemersPelag_plot<-factor(data.rmr$DemersPelag_plot, c(
                                  "benthopelagic-ecol_relev",
                                  "demersal-ecol_relev" ,
                                  "pelagic-ecol_relev", 
                                  "reef-associated-ecol_relev",
                                  "benthopelagic-warm",
                                  "demersal-warm" ,
                                  "pelagic-warm" ,
                                  "reef-associated-warm"
                                  ))

data.as$DemersPelag_plot<-factor(data.as$DemersPelag_plot, c(
                                  "benthopelagic-ecol_relev",
                                  "demersal-ecol_relev" ,
                                  "pelagic-ecol_relev", 
                                  "reef-associated-ecol_relev",
                                  "benthopelagic-warm",
                                  "demersal-warm" ,
                                  "pelagic-warm" ,
                                  "reef-associated-warm"
                                  ))

data.amr$DemersPelag_plot<-factor(data.amr$DemersPelag_plot, c(
                                    "benthopelagic-ecol_relev",
                                  "demersal-ecol_relev" ,
                                  "pelagic-ecol_relev", 
                                  "reef-associated-ecol_relev",
                                  "benthopelagic-warm",
                                  "demersal-warm" ,
                                  "pelagic-warm" ,
                                  "reef-associated-warm"
                                  ))

ecolAS1<-ggplot(data=data.as, aes(y=mass_specas, fill=DemersPelag_plot,
                                  x=DemersPelag_plot,
                                  group = DemersPelag_plot))+
  geom_violin(data=data.rmr, aes(y=mass_specrmr, fill=DemersPelag_plot,
                                 x=DemersPelag_plot,
                                group = DemersPelag_plot), 
              alpha=1, size=0.2)+
  geom_violin(data=data.amr, aes(y=mass_specamr, fill=DemersPelag_plot,
                                 x=DemersPelag_plot,
                                 group = DemersPelag_plot), 
              alpha=0.3, size=0.2)+  
  geom_pointrange(data.as, mapping=aes(x=DemersPelag_plot, y = mass_specas,
                                       fill=DemersPelag_plot, 
                                       group = DemersPelag_plot),
                  stat = "summary", color = "black", pch=23)+
  scale_fill_manual(values = c("#781c6d", 
                               "#bc3754",
                               "#ed6925", 
                               "#fcffa4", 
                               "#781c6d", 
                               "#bc3754",
                               "white", 
                               "#fcffa4"))+
  scale_x_discrete(labels=c( "benthopelagic-ecol_relev" = "Bentho- \n pelagic",
                                  "demersal-ecol_relev" = "Demersal" ,
                                  "pelagic-ecol_relev" = "Pelagic",
                                  "reef-associated-ecol_relev" = "Reef- \n associac.",
                                  "benthopelagic-warm" = "Bentho- \n pelagic",
                                  "demersal-warm" = "Demersal" ,
                                  "reef-associated-warm" = "Reef- \n associac.",
                                  "pelagic-warm" = "Pelagic"))+
  ylim(0,8)+
  geom_segment(aes(x = 7.02, y = 6, xend = 7.2, yend = 6.2), size = 0.1)+
  geom_segment(aes(x = 7.02, y = 3.5, xend = 7.2, yend = 4.4), size = 0.1)+
  geom_segment(aes(x = 7.02, y = 1.2, xend = 7.2, yend = 3.5), size = 0.1)+
  annotate(geom = "text", y = 6.25, x = 7.21, label = "MMR", size = 3, hjust = 0)+
  annotate(geom = "text", y = 4.45, x = 7.21, label = "AS", size = 3, hjust = 0)+
  annotate(geom = "text", y = 3.55, x = 7.21, label = "RMR", size = 3, hjust = 0)+
  annotate(geom = "text", y = 7.9, x = 1, label = "(RMR, MMR) BIC *", size = 3, hjust = 0)+
  annotate(geom = "text", y = 7.3, x = 1, label = "(AS) BIC ns", size = 3, hjust = 0)+
  annotate(geom = "text", y = 7.9, x = 5.0, label = "(all) BIC ns", size = 3, hjust = 0)+
  annotate(geom = "text", y = 7.55, x = 4.3, label = "OPTIMAL\n TEMPERATURES",
           hjust = 1, fontface = "bold", size = 3)+
  annotate(geom = "text", y = 7.9, x = 8.2, label = "WARM",size = 3, hjust = 1, fontface = "bold")+
  geom_vline(xintercept = 4.5, color = "grey", linetype = "dashed")
ggformat(ecolAS1, y_title = bquote("MR" ~ (mgO[2] ~ g^-1 ~ h^-1)), x_title = element_blank() , print=F, size_text = 11)
ecolAS1 <- ecolAS1 + theme(
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  # axis.text.x = element_text(angle = 45, size = 8, vjust = 0.8),
  axis.text.x = element_blank(),
  legend.position = "none",
  legend.title = element_blank(),
  legend.direction = "vertical",
  legend.justification='center',
  legend.margin=margin(0,0,0,0),
  legend.box.margin=margin(0,0,6,0),
  plot.margin = margin(0, 5, -6, 5))
ecolAS1

### Body shape MMR, RMR, AS ------
data.as$BodyShapeI_plot<-paste(data.as$BodyShapeI, data.as$test_category3, sep = "-")
data.amr$BodyShapeI_plot<-paste(data.amr$BodyShapeI, data.amr$test_category3, sep = "-")
data.rmr$BodyShapeI_plot<-paste(data.rmr$BodyShapeI, data.rmr$test_category3, sep = "-")
data.fas$BodyShapeI_plot<-paste(data.fas$BodyShapeI, data.fas$test_category3, sep = "-")

data.as$BodyShapeI_plot<-factor(data.as$BodyShapeI_plot, c(
                            "elongated-ecol_relev",
                            "fusiform-ecol_relev",
                            "short/deep-ecol_relev",
                            "elongated-warm",
                            "fusiform-warm",
                            "short/deep-warm") )
data.rmr$BodyShapeI_plot<-factor(data.rmr$BodyShapeI_plot, c(
                            "dorsoventrflattened-ecol_relev",
                            "eel-like-ecol_relev",
                            "elongated-ecol_relev",
                            "fusiform-ecol_relev",
                            "short/deep-ecol_relev",
                            "elongated-warm",
                            "fusiform-warm",
                            "short/deep-warm") )
data.amr$BodyShapeI_plot<-factor(data.amr$BodyShapeI_plot, c(                            
                            "elongated-ecol_relev",
                            "fusiform-ecol_relev",
                            "short/deep-ecol_relev",
                            "elongated-warm",
                            "fusiform-warm",
                            "short/deep-warm") )
data.fas$BodyShapeI_plot<-factor(data.fas$BodyShapeI_plot, c(                            
                            "elongated-ecol_relev",
                            "fusiform-ecol_relev",
                            "short/deep-ecol_relev",
                            "elongated-warm",
                            "fusiform-warm",
                            "short/deep-warm") )

ecolAS2<-ggplot(data=data.as, aes(y=mass_specas, fill=BodyShapeI_plot,
                                  x=BodyShapeI_plot))+
  geom_violin(data=data.rmr, aes(y=mass_specrmr, fill=BodyShapeI_plot,
                                 x=BodyShapeI_plot,),
              alpha=1, size=0.2, drop = FALSE)+
  geom_violin(data=data.amr, aes(y=mass_specamr, fill=BodyShapeI_plot,
                                 x=BodyShapeI_plot),
              alpha=0.3, size=0.2)+ 
  geom_pointrange(data.as, mapping=aes(x=BodyShapeI_plot,
                                       y = mass_specas,
                                       fill=BodyShapeI_plot),
                  stat = "summary", fun.ymin = min,fun.ymax = max,mfun.y = mean, 
                  # position=position_dodge2(width = c(0.9,1)), 
                  color = "black", pch=23)+
  scale_fill_manual(values = c("white", "white",
                               "#fde725", "#5ec962", "#440154",
                               "white", "#5ec962", "white"))+
  scale_x_discrete(labels=c("dorsoventrflattened-ecol_relev" = "Dorsoventr. \n flat",
                            "eel-like-ecol_relev" = "Eel-like",
                            "elongated-ecol_relev" = "Elongated",
                            "fusiform-ecol_relev" = "Fusiform",
                            "short/deep-ecol_relev" = "Short/Deep",
                            "elongated-warm" = "Elongated",
                            "fusiform-warm" = "Fusiform",
                            "short/deep-warm" = "Short/Deep"))+
  ylim(0,8)+
  annotate(geom = "text", y = 7.9, x = 1, label = "(RMR, RMR) BIC *", size = 3, hjust =0)+
  annotate(geom = "text", y = 7.3, x = 1, label = "(AS) BIC ns", size = 3, hjust =0)+
  annotate(geom = "text", y = 7.9, x = 6, label = "(AS) BIC *", size = 3, hjust =0)+
  annotate(geom = "text", y = 7.3, x = 6, label = "(RMR, MMR) BIC ns", size = 3, hjust =0)+
  geom_vline(xintercept = 5.5, color = "grey", linetype = "dashed")
ggformat(ecolAS2, y_title = bquote("MR" ~ (mgO[2] ~ g^-1 ~ h^-1)), x_title = element_blank() , print=F, size_text = 11)
ecolAS2 <- ecolAS2 + theme(
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  # axis.text.x = element_text(angle = 45, size = 8, vjust = 0.8),
  axis.text.x = element_blank(),
  legend.position = "none",
  legend.title = element_blank(),
  legend.direction = "vertical",
  legend.justification='center',
  legend.margin=margin(0,0,0,0),
  legend.box.margin=margin(0,0,6,0),
  plot.margin = margin(0, 5, -6, 5))
ecolAS2

### Salinity  MMR, RMR, AS ------
data.as$salintyComb_plot<-paste(data.as$salintyComb, data.as$test_category3, sep = "-")
data.amr$salintyComb_plot<-paste(data.amr$salintyComb, data.amr$test_category3, sep = "-")
data.rmr$salintyComb_plot<-paste(data.rmr$salintyComb, data.rmr$test_category3, sep = "-")
data.fas$salintyComb_plot<-paste(data.fas$salintyComb, data.fas$test_category3, sep = "-")

data.rmr$salintyComb_plot<-factor(data.rmr$salintyComb_plot, c(
      "Freshwater-ecol_relev",
      "Freshwater; brackish-ecol_relev",
      "Marine-ecol_relev",
      "Marine; brackish-ecol_relev",
      "Marine; freshwater; brackish-ecol_relev",
      "Freshwater-warm",
      "Freshwater; brackish-warm",
      "Marine-warm",
      "Marine; brackish-warm",
      "Marine; freshwater; brackish-warm"
))
data.fas$salintyComb_plot<-factor(data.fas$salintyComb_plot, c(
      "Freshwater-ecol_relev",
      "Freshwater; brackish-ecol_relev",
      "Marine-ecol_relev",
      "Marine; brackish-ecol_relev",
      "Marine; freshwater; brackish-ecol_relev",
      "Freshwater-warm",
      "Freshwater; brackish-warm",
      "Marine-warm",
      "Marine; brackish-warm",
      "Marine; freshwater; brackish-warm"
))
data.amr$salintyComb_plot<-factor(data.amr$salintyComb_plot, c(
      "Freshwater-ecol_relev",
      "Freshwater; brackish-ecol_relev",
      "Marine-ecol_relev",
      "Marine; brackish-ecol_relev",
      "Marine; freshwater; brackish-ecol_relev",
      "Freshwater-warm",
      "Freshwater; brackish-warm",
      "Marine-warm",
      "Marine; brackish-warm",
      "Marine; freshwater; brackish-warm"
))
data.as$salintyComb_plot<-factor(data.as$salintyComb_plot, c(
      "Freshwater-ecol_relev",
      "Freshwater; brackish-ecol_relev",
      "Marine-ecol_relev",
      "Marine; brackish-ecol_relev",
      "Marine; freshwater; brackish-ecol_relev",
      "Freshwater-warm",
      "Freshwater; brackish-warm",
      "Marine-warm",
      "Marine; brackish-warm",
      "Marine; freshwater; brackish-warm"
))

ecolAS3<-ggplot(data=data.as, aes(y=mass_specas, fill=salintyComb_plot, x=salintyComb_plot))+
  geom_violin(data=data.rmr, aes(y=mass_specrmr, fill=salintyComb_plot, x=salintyComb_plot), 
              alpha=1, size=0.2)+
  geom_violin(data=data.amr, aes(y=mass_specamr, fill=salintyComb_plot, x=salintyComb_plot), 
              alpha=0.3, size=0.2)+  
  geom_pointrange(data.as, mapping=aes(x=salintyComb_plot, y = mass_specas,
                                       fill=salintyComb_plot),
                  stat = "summary", fun.ymin = min,fun.ymax = max,fun.y = mean, 
                  color = "black", pch=23)+
  scale_fill_manual(values = c("#fcfdbf", "#fe9f6d", "#de4968",
                               "white", "#3b0f70",
                               "white", "white", "#de4968",
                               "white", "#3b0f70"))+
  scale_x_discrete(labels=c( 
      "Freshwater-ecol_relev" = "Freshwater",
      "Freshwater; brackish-ecol_relev" = "Freshwater/ \n brackish",
      "Marine-ecol_relev" = "Marine",
      "Marine; brackish-ecol_relev" = "Marine/ \n brackish",
      "Marine; freshwater; brackish-ecol_relev" = "All \n salinities",
      "Freshwater-warm" = "Freshwater",
      "Freshwater; brackish-warm" ="Freshwater/ \n brackish",
      "Marine-warm" = "Marine",
      "Marine; brackish-warm" = "Marine/ \n brackish",
      "Marine; freshwater; brackish-warm" = "All \n salinities"))+
  ylim(0,8)+
  annotate(geom = "text", y = 7.9, x = 1, label = "(RMR, MMR) BIC *", size = 3, hjust =0)+
  annotate(geom = "text", y = 7.3, x = 1, label = "(AS) BIC ns", size = 3, hjust =0)+
  annotate(geom = "text", y = 7.9, x = 6.05, label = "(all) BIC ns", size = 3, hjust= 0)+
  geom_vline(xintercept = 5.5, color = "grey", linetype = "dashed")
ggformat(ecolAS3, y_title = bquote("MR" ~ (mgO[2] ~ g^-1 ~ h^-1)), x_title = element_blank() , print=F, size_text = 11)
ecolAS3 <- ecolAS3 + theme(
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  # axis.text.x = element_text(angle = 45, size = 8, vjust = 0.8),
  axis.text.x = element_blank(),
  legend.position = "none",
  legend.title = element_blank(),
  legend.direction = "vertical",
  legend.justification='center',
  legend.margin=margin(0,0,0,0),
  legend.box.margin=margin(0,0,6,0),
  plot.margin = margin(0, 5, -6, 5))
ecolAS3

### Climate MMR, RMR, AS ------
data.as$Climate_plot<-paste(data.as$Climate, data.as$test_category3, sep = "-")
data.amr$Climate_plot<-paste(data.amr$Climate, data.amr$test_category3, sep = "-")
data.rmr$Climate_plot<-paste(data.rmr$Climate, data.rmr$test_category3, sep = "-")
data.fas$Climate_plot<-paste(data.fas$Climate, data.fas$test_category3, sep = "-")

data.fas$Climate_plot<-factor(data.fas$Climate_plot, c(
                                "Polar-ecol_relev",
                                "Temperate-ecol_relev",
                                "Subtropical-ecol_relev",
                                "Tropical-ecol_relev",
                                "Polar-warm",
                                "Temperate-warm",
                                "Subtropical-warm",
                                "Tropical-warm"
                                ))
data.rmr$Climate_plot<-factor(data.rmr$Climate_plot, c(
                                "Polar-ecol_relev",
                                "Temperate-ecol_relev",
                                "Subtropical-ecol_relev",
                                "Tropical-ecol_relev",
                                "Polar-warm",
                                "Temperate-warm",
                                "Subtropical-warm",
                                "Tropical-warm"
                                ))
data.amr$Climate_plot<-factor(data.amr$Climate_plot, c(
                                "Polar-ecol_relev",
                                "Temperate-ecol_relev",
                                "Subtropical-ecol_relev",
                                "Tropical-ecol_relev",
                                "Polar-warm",
                                "Temperate-warm",
                                "Subtropical-warm",
                                "Tropical-warm"
                                ))
data.as$Climate_plot<-factor(data.as$Climate_plot, c(
                                "Polar-ecol_relev",
                                "Temperate-ecol_relev",
                                "Subtropical-ecol_relev",
                                "Tropical-ecol_relev",
                                "Polar-warm",
                                "Temperate-warm",
                                "Subtropical-warm",
                                "Tropical-warm"
                                ))

ecolAS4<-ggplot(data=data.as, aes(y=mass_specas, fill=Climate_plot, x=Climate_plot))+
  geom_violin(data=data.rmr, aes(y=mass_specrmr, fill=Climate_plot, x=Climate_plot), 
              alpha=1, size=0.2)+
  geom_violin(data=data.amr, aes(y=mass_specamr, fill=Climate_plot, x=Climate_plot), 
              alpha=0.3, size=0.2)+  # geom_boxplot( outlier.shape = NA) +
  geom_pointrange(data.as, mapping=aes(x=Climate_plot, y = mass_specas,
                                       fill=Climate_plot),
                  stat = "summary", fun.ymin = min,fun.ymax = max,fun.y = mean, color = "black", pch=23)+
  scale_fill_manual(values = c("white", "#fca636","#f0f921", "#b12a90",
                               "white", "#fca636","#f0f921", "#b12a90"))+
  scale_x_discrete(labels=c("Polar-ecol_relev" = "Polar",
                                "Temperate-ecol_relev" = "Temperate",
                                "Subtropical-ecol_relev" = "Subtropical",
                                "Tropical-ecol_relev" = "Tropical",
                                "Polar-warm" = "Polar",
                                "Temperate-warm" = "Temperate",
                                "Subtropical-warm" = "Subtropical",
                                "Tropical-warm" = "Tropical"))+
  ylim(0,8)+
  annotate(geom = "text", y = 7.9, x = 1, label = "(RMR, MMR) BIC *", size = 3, hjust =0)+
  annotate(geom = "text", y = 7.3, x = 1, label = "(AS) BIC ns", size = 3, hjust =0)+
  annotate(geom = "text", y = 7.9, x = 5.0, label = "(all) BIC ns", size = 3, hjust = 0)+
  geom_vline(xintercept = 4.5, color = "grey", linetype = "dashed")
ggformat(ecolAS4, y_title = bquote("MR" ~ (mgO[2] ~ g^-1 ~ h^-1)), x_title = element_blank() , print=F, size_text = 11)
ecolAS4 <- ecolAS4 + theme(
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  # axis.text.x = element_text(angle = 45, size = 8, vjust = 0.8),
  axis.text.x = element_blank(),
  legend.position = "none",
  legend.title = element_blank(),
  legend.direction = "vertical",
  legend.justification='center',
  legend.margin=margin(0,0,0,0),
  legend.box.margin=margin(0,0,6,0),
  plot.margin = margin(0, 5, -6, 5))
ecolAS4

### Demersal Pelagic FAS -----------
ecolFAS1<-ggplot(data=data.fas, aes(y=FAS, fill=DemersPelag_plot, x=DemersPelag_plot))+
  geom_violin(alpha = 1, outlier.size = 0.2, size = 0.2)+
  scale_fill_manual(values = c("#781c6d", 
                               "#bc3754",
                               "#ed6925", 
                               "#fcffa4", 
                               "#781c6d", 
                               "#bc3754",
                               "white", 
                               "#fcffa4"))+ 
  scale_x_discrete(labels=c( "benthopelagic-ecol_relev" = "Bentho- \n pelagic",
                                  "demersal-ecol_relev" = "Demersal" ,
                                  "pelagic-ecol_relev" = "Pelagic",
                                  "reef-associated-ecol_relev" = "Reef- \n associac.",
                                  "benthopelagic-warm" = "Bentho- \n pelagic",
                                  "demersal-warm" = "Demersal" ,
                                  "pelagic-warm" = "Pelagic",
                                  "reef-associated-warm" = "Reef- \n associac."
                                  ))+
  ylim(0,20)+
  annotate(geom = "text", y = 19.5, x = 1, label = "BIC ns", size = 3, hjust = 0)+
  annotate(geom = "text", y = 19.5, x = 5.0, label = "BIC ns", size = 3, hjust = 0)+
  geom_vline(xintercept = 4.5, color = "grey", linetype = "dashed")
ggformat(ecolFAS1, y_title = bquote("FAS"), x_title = element_blank() , print=T, size_text = 11)
ecolFAS1 <- ecolFAS1 + theme(
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  legend.position = "none",
  legend.direction = "vertical",
  legend.justification='center',
  axis.text.x = element_text(angle = 45, size = 10, vjust = 0.8),
  legend.margin=margin(0,0,0,0),
  legend.box.margin=margin(0,0,6,0),
  plot.margin = margin(-12, 5, -5, 5),
  legend.spacing.y = unit(0.5, 'cm'),
  legend.text = element_text(margin = margin(l = 0.1, unit = c("cm"))),
  legend.key.height= unit(0.5, 'cm'),
  legend.key.width= unit(0.5, 'cm'))
ecolFAS1

### Body shape FAS -----------
# add empty layers for plotting x-axis. 
data.fas$BodyShapeI_plot<-factor(data.fas$BodyShapeI_plot, c(
  "dorsoventrflattened-ecol_relev" ,
  "eel-like-ecol_relev",
  "elongated-ecol_relev",
  "fusiform-ecol_relev",
  "short/deep-ecol_relev",
  "elongated-warm",
  "fusiform-warm",
  "short/deep-warm" ))
  
ecolFAS2<-ggplot(data=data.fas, aes(FAS, fill=BodyShapeI_plot, x=BodyShapeI_plot))+
  geom_violin(alpha = 1, outlier.size = 0.2, size = 0.2, drop = FALSE)+
  scale_fill_manual(values = c("#fde725", "#5ec962", "#440154",
                               "white", "#5ec962", "white"))+
  scale_x_discrete(labels=c("dorsoventrflattened-ecol_relev" = "Dorsoventr. \n flat",
                            "eel-like-ecol_relev" = "Eel-like",
                            "elongated-ecol_relev" = "Elongated",
                            "fusiform-ecol_relev" = "Fusiform",
                            "short/deep-ecol_relev" = "Short/Deep",
                            "elongated-warm" = "Elongated",
                            "fusiform-warm" = "Fusiform",
                            "short/deep-warm" = "Short/Deep"), drop = FALSE)+
                   # expand = expansion(mult = c(0.5, 0)))+
  ylim(0,20)+
  annotate(geom = "text", y = 19.5, x = 1, label = "BIC ns", size = 3, hjust = 0)+
  annotate(geom = "text", y = 19.5, x = 6.05, label = "BIC ns", size = 3, hjust = 0)+
  geom_vline(xintercept = 5.5, color = "grey", linetype = "dashed")
ggformat(ecolFAS2, y_title = bquote("FAS"), x_title = element_blank() , print=F, size_text = 11)
ecolFAS2 <- ecolFAS2 + theme( 
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  legend.position = "none",
  axis.text.x = element_text(angle = 45, size = 10, vjust = 0.8),
  legend.direction = "vertical",
  legend.justification='center',
  legend.margin=margin(0,0,0,0),
  legend.box.margin=margin(0,0,6,0),
  plot.margin = margin(-12, 5, -5, 5),
  legend.spacing.y = unit(0.5, 'cm'),
  legend.text = element_text(margin = margin(l = 0.1, unit = c("cm"))),
  legend.key.height= unit(0.5, 'cm'),
  legend.key.width= unit(0.5, 'cm'))
ecolFAS2

### Salinity violin FAS ------
ecolFAS3<-ggplot(data=data.fas, aes(FAS, fill=salintyComb_plot, x=salintyComb_plot))+
  geom_violin(alpha = 1, outlier.size = 0.2, size = 0.2)+
  scale_fill_manual(values = c("#fcfdbf", "#fe9f6d", "#de4968",
                               "white", "#3b0f70",
                               "white", "white", "#de4968",
                               "white", "#3b0f70"))+
  scale_x_discrete(labels=c( 
      "Freshwater-ecol_relev" = "Freshwater",
      "Freshwater; brackish-ecol_relev" = "Freshwater/ \n brackish",
      "Marine-ecol_relev" = "Marine",
      "Marine; brackish-ecol_relev" = "Marine/ \n brackish",
      "Marine; freshwater; brackish-ecol_relev" = "All \n salinities",
      "Freshwater-warm" = "Freshwater",
      "Freshwater; brackish-warm" ="Freshwater/ \n brackish",
      "Marine-warm" = "Marine",
      "Marine; brackish-warm" = "Marine/ \n brackish",
      "Marine; freshwater; brackish-warm" = "All \n salinities"))+
  ylim(0,20)+
  annotate(geom = "text", y = 19.5, x = 1, label = "BIC ns", size = 3, hjust = 0)+
  annotate(geom = "text", y = 19.5, x = 6.05, label = "BIC ns", size = 3, hjust = 0)+
  geom_vline(xintercept = 5.5, color = "grey", linetype = "dashed")
ggformat(ecolFAS3, y_title = bquote("FAS"), x_title = element_blank() , print=F, size_text = 11)
ecolFAS3 <- ecolFAS3 + theme(
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  legend.position = "none",
  legend.direction = "vertical",
  legend.justification='center',
  axis.text.x = element_text(angle = 45, size = 10, vjust = 0.8),
  legend.margin=margin(0,0,0,0),
  legend.box.margin=margin(0,0,6,0),
  plot.margin = margin(-12, 5, -5, 5),
  legend.spacing.y = unit(0.5, 'cm'),
  legend.text = element_text(margin = margin(l = 0.1, unit = c("cm"))),
  legend.key.height= unit(0.5, 'cm'),
  legend.key.width= unit(0.5, 'cm'))
ecolFAS3

### Climate violin FAS ----------
ecolFAS4<-ggplot(data=data.fas, aes(FAS, fill=Climate_plot, x=Climate_plot))+
  geom_violin(alpha = 1, outlier.size = 0.2, size = 0.2)+
  scale_fill_manual(values = c("white", "#fca636","#f0f921", "#b12a90",
                               "white", "#fca636","#f0f921", "#b12a90"))+
  scale_x_discrete(labels=c("Polar-ecol_relev" = "Polar",
                                "Temperate-ecol_relev" = "Temperate",
                                "Subtropical-ecol_relev" = "Subtropical",
                                "Tropical-ecol_relev" = "Tropical",
                                "Polar-warm" = "Polar",
                                "Temperate-warm" = "Temperate",
                                "Subtropical-warm" = "Subtropical",
                                "Tropical-warm" = "Tropical"))+
  ylim(0,20)+
  annotate(geom = "text", y = 19.5, x = 1, label = "BIC *", size = 3, hjust = 0)+
  annotate(geom = "text", y = 19.5, x = 5, label = "BIC ns", size = 3, hjust = 0)+
  geom_vline(xintercept = 4.5, color = "grey", linetype = "dashed")
ggformat(ecolFAS4, y_title = bquote("FAS"), x_title = element_blank() , print=F, size_text = 11)
ecolFAS4 <- ecolFAS4 + theme(
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  legend.position = "none",
  legend.direction = "vertical",
  legend.justification='center',
  axis.text.x = element_text(angle = 45, size = 10, vjust = 0.8),
  legend.margin=margin(0,0,0,0),
  legend.box.margin=margin(0,0,6,0),
  plot.margin = margin(-12, 5, -5, 5),
  legend.spacing.y = unit(0.5, 'cm'),
  legend.text = element_text(margin = margin(l = 0.1, unit = c("cm"))),
  legend.key.height= unit(0.5, 'cm'),
  legend.key.width= unit(0.5, 'cm'))

#### save violin plots -------
 plot_grid(ecolAS1, ecolAS2,
          ecolFAS1, ecolFAS2,
          ecolAS3, ecolAS4, 
          ecolFAS3, ecolFAS4,
          align = "v",
          labels = c("AUTO"),
          label_size = 12,
          nrow = 4,
          ncol = 2,
          label_x = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
          label_y = c(0.88, 0.88, 0.92, 0.92,0.88, 0.88, 0.92, 0.92 )) %>% 
ggsave(filename = "./Figures/Figure5.png", width = 12, height = 9)
 

 # TPC: FAS, AS, MR, mass-independent values ----------
AMRmodel_plot1_t<-ggplot(data=data.amrER, aes(x=tempTest, y=mass_specamr)) +
  geom_point(alpha=0.9,  size=1, pch=1, color="grey75")+
  geom_point(data=data.amr.test, aes(x=tempTest, y=mass_specamr),
             alpha=0.9,  size=1, pch=21, show.legend = FALSE, stroke =0.2,
             fill = cols.amr[3], color = cols.amr[1])+
  scale_color_gradient( low = "grey", high = "black")
ggformat(AMRmodel_plot1_t, x_title=expression("Temperature ºC"),
         y_title=expression(MMR~(mg~O[2]~h^-1~g^-1)), size_text = 12,print = T)

RMRmodel_plot1_t<-ggplot(data=data.rmrER, aes(x=tempTest, y=mass_specrmr)) +
  geom_point(alpha=0.9,  size=1, pch=1, color="grey75")+
  geom_point(data=data.rmr.test, aes(x=tempTest, y=mass_specrmr),
             alpha=0.9,  size=1, pch=21, show.legend = FALSE, stroke =0.2,
             fill = cols.rmr[3], color = cols.rmr[1])+
  scale_color_gradient( low = "grey", high = "black")
ggformat(RMRmodel_plot1_t, x_title=expression("Temperature ºC"),
         y_title=expression(RMR~(mg~O[2]~h^-1~g^-1)), size_text = 12,print = T)

# FAS! 
FASmodel_plot1_t<-ggplot(data=data.fasER, aes(x=tempTest, y=FAS)) +
  geom_point(alpha=0.9,  size=1, pch=1, color="grey75")+
  geom_point(data=data.fas.test, aes(x=tempTest, y=FAS),
             alpha=0.9,  size=1, pch=21, show.legend = FALSE, stroke =0.2,
             fill = cols.fas[3], color = cols.fas[1])+
  scale_color_gradient( low = "grey", high = "black")+
  ylim(0,20)
ggformat(FASmodel_plot1_t, x_title=expression("Temperature ºC"),
         y_title=expression(FAS~(MMR/RMR)), size_text = 12,print = T)

# data.fasER[data.fasER$FAS>20, ]
  
scaling_t<-cowplot:::plot_grid(AMRmodel_plot1_t, RMRmodel_plot1_t, FASmodel_plot1_t,
                              align = "hv",
                              axis = "l",
                              nrow = 1,
                              ncol = 3,
                              labels = "AUTO",
                              label_size = 12)

ggsave(filename = paste("./Figures/Suppl_temp_mass_spec.png", sep=""),
       plot=scaling_t, width = 9, height = 3, units = "in")


# Boxplots: FAS, AS, MR, mass-independent values ----------
fas_boxplot<-ggplot(data=data.fas, aes(y=FAS, x = test_category,  fill=test_category3))+
  geom_boxplot(show.legend = F)+
  scale_fill_manual(values=cols.fas)+
  ylim(0,20)+
  scale_x_discrete(labels=c("acclim" = "Acclimated \n warm", "acute" = "Acute \n warm", "ecol_relev" = "Optimal"))
ggformat(fas_boxplot, y_title = expression(FAS~(MMR/RMR)), x_title = "", print = F)
fas_boxplot<-fas_boxplot+theme(legend.position = "none")

rmr_boxplot<-ggplot(data=data.rmr, aes(y=mass_specrmr, x = test_category, fill=test_category))+
  geom_boxplot(show.legend = F)+
  scale_fill_manual(values=cols.rmr)+
  scale_x_discrete(labels=c("acclim" = "Acclimated \n warm", "acute" = "Acute \n warm", "ecol_relev" = "Optimal"))
ggformat(rmr_boxplot, y_title = expression(RMR~(mgO[2]~g^-1~h^-1)), x_title = "", print = F)
rmr_boxplot<-rmr_boxplot+theme(legend.position = "none")

amr_boxplot<-ggplot(data=data.amr, aes(y=mass_specamr, x = test_category, fill=test_category))+
  geom_boxplot(show.legend = F)+
  scale_fill_manual(values=cols.amr)+
  scale_x_discrete(labels=c("acclim" = "Acclimated \n warm", "acute" = "Acute \n warm", "ecol_relev" = "Optimal"))
ggformat(amr_boxplot, y_title = expression(MMR~(mgO[2]~g^-1~h^-1)), x_title = "", print = F)
amr_boxplot<-amr_boxplot+theme(legend.position = "none")

fas_boxplot2<-ggplot(data=data.fas, aes(y=FAS, x = test_category3,  fill=test_category3))+
  geom_boxplot(show.legend = F)+
  scale_fill_manual(values=cols.fas)+
  ylim(0,20)+
  scale_x_discrete(labels=c("warm" = "Warm", "ecol_relev" = "Optimal"))
ggformat(fas_boxplot2, y_title = expression(FAS~(MMR/RMR)), x_title = "", print = F)
fas_boxplot2<-fas_boxplot2+theme(legend.position = "none")

rmr_boxplot2<-ggplot(data=data.rmr, aes(y=mass_specrmr, x = test_category3, fill=test_category3))+
  geom_boxplot(show.legend = F)+
  scale_fill_manual(values=cols.rmr)+
  scale_x_discrete(labels=c("warm" = "Warm", "ecol_relev" = "Optimal"))
ggformat(rmr_boxplot2, y_title = expression(RMR~(mgO[2]~g^-1~h^-1)), x_title = "", print = F)
rmr_boxplot2<-rmr_boxplot2+theme(legend.position = "none")

amr_boxplot2<-ggplot(data=data.amr, aes(y=mass_specamr, x = test_category3, fill=test_category3))+
  geom_boxplot(show.legend = F)+
  scale_fill_manual(values=cols.amr)+
  scale_x_discrete(labels=c("warm" = "Warm", "ecol_relev" = "Optimal"))
ggformat(amr_boxplot2, y_title = expression(MMR~(mgO[2]~g^-1~h^-1)), x_title = "", print = F)
amr_boxplot2<-amr_boxplot2+theme(legend.position = "none")


cowplot:::plot_grid(amr_boxplot,rmr_boxplot,fas_boxplot,
                    amr_boxplot2,rmr_boxplot2,fas_boxplot2,
          align = "hv",
          axis = "l",
          nrow = 2,
          labels = "AUTO",
          ncol = 3) %>%
ggsave(filename = "./Figures/Supl_fig1_boxplots_Phylo.png", width = 12.5, height =10)



