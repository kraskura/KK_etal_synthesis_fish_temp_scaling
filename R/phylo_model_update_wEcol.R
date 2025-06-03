

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
# Demersal Pelagic MMR, RMR, AS ------
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
  geom_segment(aes(x = 2.12, y = 6, xend = 2.2, yend = 6.2), size = 0.1)+
  geom_segment(aes(x = 2.12, y = 3.5, xend = 2.2, yend = 4.4), size = 0.1)+
  geom_segment(aes(x = 2.12, y = 1.2, xend = 2.2, yend = 3.5), size = 0.1)+
  annotate(geom = "text", y = 6.25, x = 2.21, label = "MMR", size = 3, hjust =0)+
  annotate(geom = "text", y = 4.45, x = 2.21, label = "AS", size = 3, hjust =0)+
  annotate(geom = "text", y = 3.55, x = 2.21, label = "RMR", size = 3, hjust =0)+
  
  annotate(geom = "text", y = 7.9, x = 0.65, label = "(RMR, MMR) BIC *", size = 3, hjust =0)+
  annotate(geom = "text", y = 7.3, x = 0.65, label = "(AS) BIC ns", size = 3, hjust =0)+
  annotate(geom = "text", y = 7.9, x = 2.05, label = "(all) BIC ns", size = 3)+
  geom_vline(xintercept = 1.5, color = "grey", linetype = "dashed")
ggformat(ecolAS1, y_title = bquote("MR" ~ (mgO[2] ~ g^-1 ~ h^-1)), x_title = element_blank() , print=F, size_text = 11)

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

# Body shape MMR, RMR, AS ------
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
  annotate(geom = "text", y = 7.9, x = 0.65, label = "(RMR, RMR) BIC *", size = 3, hjust =0)+
  annotate(geom = "text", y = 7.3, x = 0.65, label = "(AS) BIC ns", size = 3, hjust =0)+
  annotate(geom = "text", y = 7.9, x = 2.05, label = "(AS) BIC *", size = 3)+
  annotate(geom = "text", y = 7.3, x = 2.05, label = "(RMR, MMR) BIC ns", size = 3)+
  geom_vline(xintercept = 1.5, color = "grey", linetype = "dashed")
ggformat(ecolAS2, y_title = bquote("MR" ~ (mgO[2] ~ g^-1 ~ h^-1)), x_title = element_blank() , print=F, size_text = 11)
ecolAS2 <- ecolAS2 + theme(
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  legend.position = "none",
  legend.title = element_blank(),
  legend.direction = "vertical",
  legend.justification='center',
  legend.margin=margin(0,0,0,0),
  legend.box.margin=margin(0,0,6,0),
  plot.margin = margin(5.5, 5.5, 5.35, 5.5))

# Salinity  MMR, RMR, AS ------
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
  annotate(geom = "text", y = 7.9, x = 0.65, label = "(RMR, MMR) BIC *", size = 3, hjust =0)+
  annotate(geom = "text", y = 7.3, x = 0.65, label = "(AS) BIC ns", size = 3, hjust =0)+
  annotate(geom = "text", y = 7.9, x = 2.05, label = "(all) BIC ns", size = 3)+
  geom_vline(xintercept = 1.5, color = "grey", linetype = "dashed")
ggformat(ecolAS3, y_title = bquote("MR" ~ (mgO[2] ~ g^-1 ~ h^-1)), x_title = element_blank() , print=F, size_text = 11)

ecolAS3 <- ecolAS3 + theme(
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  legend.position = "none",
  legend.title = element_blank(),
  legend.direction = "vertical",
  legend.justification='center',
  legend.margin=margin(0,0,0,0),
  legend.box.margin=margin(0,0,6,0),
  plot.margin = margin(5.5, 5.5, 5.35, 5.5))

# Climate MMR, RMR, AS ------
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
  annotate(geom = "text", y = 7.9, x = 0.65, label = "(RMR, MMR) BIC *", size = 3, hjust =0)+
  annotate(geom = "text", y = 7.3, x = 0.65, label = "(AS) BIC ns", size = 3, hjust =0)+
  annotate(geom = "text", y = 7.9, x = 2.05, label = "(all) BIC ns", size = 3)+
  geom_vline(xintercept = 1.5, color = "grey", linetype = "dashed")
ggformat(ecolAS4, y_title = bquote("MR" ~ (mgO[2] ~ g^-1 ~ h^-1)), x_title = element_blank() , print=F, size_text = 11)
ecolAS4 <- ecolAS4 + theme(
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  legend.position = "none",
  legend.title = element_blank(),
  legend.direction = "vertical",
  legend.justification='center',
  legend.margin=margin(0,0,0,0),
  legend.box.margin=margin(0,0,6,0),
  plot.margin = margin(5.5, 5.5, 5.35, 5.5))

# Demersal Pelagic FAS -----------
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
ggformat(ecolFAS1, y_title = bquote("FAS"), x_title = element_blank() , print=T, size_text = 11)
ecolFAS1 <- ecolFAS1 + theme(
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  legend.position = "right",
  legend.direction = "vertical",
  legend.justification='center',
  legend.margin=margin(0,0,0,0),
  legend.box.margin=margin(0,0,6,0),
  plot.margin = margin(5.5, 5.5, 5.35, 5.5),
  legend.spacing.y = unit(0.5, 'cm'),
  legend.text = element_text(margin = margin(l = 0.1, unit = c("cm"))),
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
ggformat(ecolFAS2, y_title = bquote("FAS"), x_title = element_blank() , print=T, size_text = 11)
ecolFAS2 <- ecolFAS2 + theme( 
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  legend.position = "right",
  legend.direction = "vertical",
  legend.justification='center',
  legend.margin=margin(0,0,0,0),
  legend.box.margin=margin(0,0,6,0),
  plot.margin = margin(5.5, 5.5, 5.35, 5.5),
  legend.spacing.y = unit(0.5, 'cm'),
  legend.text = element_text(margin = margin(l = 0.1, unit = c("cm"))),
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
ggformat(ecolFAS3, y_title = bquote("FAS"), x_title = element_blank() , print=T, size_text = 11)
ecolFAS3 <- ecolFAS3 + theme(
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  legend.position = "right",
  legend.direction = "vertical",
  legend.justification='center',
  legend.margin=margin(0,0,0,0),
  legend.box.margin=margin(0,0,6,0),
  plot.margin = margin(5.5, 5.5, 5.35, 5.5),
  legend.spacing.y = unit(0.5, 'cm'),
  legend.text = element_text(margin = margin(l = 0.1, unit = c("cm"))),
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
ggformat(ecolFAS4, y_title = bquote("FAS"), x_title = element_blank() , print=T, size_text = 11)
ecolFAS4 <- ecolFAS4 + theme(
  plot.title = element_text(face = "bold", size=15, hjust = 0.5),
  legend.position = "right",
  legend.direction = "vertical",
  legend.justification='center',
  legend.margin=margin(0,0,0,0),
  legend.box.margin=margin(0,0,6,0),
  plot.margin = margin(5.5, 5.5, 5.35, 5.5),
  legend.spacing.y = unit(0.5, 'cm'),
  legend.text = element_text(margin = margin(l = 0.1, unit = c("cm"))),
  legend.key.height= unit(0.5, 'cm'),
  legend.key.width= unit(0.5, 'cm'))

plot_grid(ecolAS1,ecolFAS1,
          ecolAS2, ecolFAS2,
          ecolAS3, ecolFAS3, 
          ecolAS4, ecolFAS4,
          align = "h",
          labels = c("AUTO"),
          label_size = 12,
          nrow = 4,
          ncol = 2,
          label_x = c(0.163, 0.123),
          label_y = c(0.85, 0.85), 
          rel_widths = c(1, 1.35)) %>% 
ggsave(filename = "./Figures/Figure4.png", width = 7.5, height = 10)
   


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



