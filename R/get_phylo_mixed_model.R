# Author: Krista Kraskura
# Description: 
#   - phylogentic models 
#   - make life stage specific main text figures
#   - analyze data for all lifestages together
# ********************************************
# ********************************************
library(here)

# *************************************************
# SUPPORTING FUNCTIONS
# function to order model selection based on the lowest BIC score
BICdelta<-function(BICtable){
  BIC.t <- BICtable [order(BICtable$BIC), ]
  BIC.t$delta <- round(abs(BIC.t$BIC[1] -  BIC.t$BIC), 5)
  return(BIC.t)
}

# making sure the correct model is extracted
check_best_model <- function(model, datatable) {
  # Extract name from lmer model
  model_name <- deparse(substitute(model))
  # Extract first rowname from datatable
  best_model_name <- rownames(datatable)[1]
  # Compare and write message
  if (model_name == best_model_name) {
    message("Correct best models identified: ", best_model_name)
  } else {
    warning("Mismatch! Model name '", model_name, "' does not match best model '", best_model_name, "'")
  }
}

# **************************************************
# **************************************************
# **************************************************

# lifestage = "juvenile"
# lifestage = NULL

get_phylo_mixed_models<-function(lifestage = NULL){ #default to all data
# Phylo models -----------------
  if (!is.null(lifestage)) {
    #  filter data to only have  specific lifestage
    data.rmrER<-data.rmrER[data.rmrER$lifestage == lifestage & !is.na(data.rmrER$lifestage == lifestage), ]
    data.rmr.test<-data.rmr.test[data.rmr.test$lifestage == lifestage & !is.na(data.rmr.test$lifestage == lifestage), ]
    data.amrER<-data.amrER[data.amrER$lifestage == lifestage & !is.na(data.amrER$lifestage == lifestage), ]
    data.amr.test<-data.amr.test[data.amr.test$lifestage == lifestage & !is.na(data.amr.test$lifestage == lifestage), ]
    data.asER<-data.asER[data.asER$lifestage == lifestage & !is.na(data.asER$lifestage == lifestage), ]
    data.as.test<-data.as.test[data.as.test$lifestage == lifestage & !is.na(data.as.test$lifestage == lifestage), ]
    data.fasER<-data.fasER[data.fasER$lifestage == lifestage & !is.na(data.fasER$lifestage == lifestage), ]
    data.fas.test<-data.fas.test[data.fas.test$lifestage == lifestage & !is.na(data.fas.test$lifestage == lifestage), ]
  }

  # get phylo trees make sure data is clean:
  data.rmrER<-droplevels(data.rmrER)
  data.rmr.test<-droplevels(data.rmr.test)
  data.amrER<-droplevels(data.amrER)
  data.amr.test<-droplevels(data.amr.test)
  data.asER<-droplevels(data.asER)
  data.as.test<-droplevels(data.as.test)
  data.fasER<-droplevels(data.fasER)
  data.fas.test<-droplevels(data.fas.test)


  # this calls custom in function in 'get_data_phylo_matrix.R'. 
  # IMPORTANT: ensure to receive a message: "All species names are identified and mathced with phylo data N species" that marks that all species have been matched with original dataset. 
  if(!is.null(lifestage)){
    if (lifestage == "adult") {
    message("Running adult data")
    # rmr
    get_phylo_matrix(species.list = unique(levels(data.rmrER$species)),
                     matrix.name = "A.adult",
                     dataset.ID = "RMR-optimal-adult")
    get_phylo_matrix(species.list = unique(levels(data.rmr.test$species)),
                     matrix.name = "A.adult.rmr.w",
                     dataset.ID = "RMR-warm-adult")
    # amr 
    get_phylo_matrix(species.list = unique(levels(data.amrER$species)),
                     matrix.name = "A.adult.mmr.er",
                     dataset.ID = "MMR-optimal-adult")
    get_phylo_matrix(species.list = unique(levels(data.amr.test$species)),
                     matrix.name = "A.adult.mmr.w",
                     dataset.ID = "MMR-warm-adult")
    # aas
    get_phylo_matrix(species.list = unique(levels(data.asER$species)),
                     matrix.name = "A.adult.aas.er",
                     dataset.ID = "AAS-optimal-adult")
    get_phylo_matrix(species.list = unique(levels(data.as.test$species)),
                     matrix.name = "A.adult.aas.w",
                     dataset.ID = "AAS-warm-adult")
    # fas
    get_phylo_matrix(species.list = unique(levels(data.fasER$species)),
                     matrix.name = "A.adult.fas.er",
                     dataset.ID = "FAS-optimal-adult")
    get_phylo_matrix(species.list = unique(levels(data.fas.test$species)),
                     matrix.name = "A.adult.fas.w",
                     dataset.ID = "FAS-warm-adult")
    
    
    ## RMR optimal -----------------
    Phylo_RMR_model2 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A.adult))
    Phylo_RMR_model2int <- Almer(lnRMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A.adult))
    
    Phylo_RMR_model4 <- Almer(lnRMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A.adult))
    Phylo_RMR_model4int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A.adult))
    
    Phylo_RMR_model5 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A.adult))
    Phylo_RMR_model5int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A.adult))
    
    
    # No singular fits 
    ### BIC rmr optimal --------------
    RMR_BIC<-BICdelta(BIC(
    
                    Phylo_RMR_model2, Phylo_RMR_model2int, #Phylo_RMR_model2.POLY, Phylo_RMR_model2intPOLY,
                    Phylo_RMR_model4, Phylo_RMR_model4int,# Phylo_RMR_model4.POLY, Phylo_RMR_model4intPOLY,
                    Phylo_RMR_model5, Phylo_RMR_model5int# Phylo_RMR_model5.POLY, Phylo_RMR_model5intPOLY
    ))
    
    ## MMR optimal -----------------
    Phylo_MMR_model2 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.adult.mmr.er))
    Phylo_MMR_model2int <- Almer(lnAMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.adult.mmr.er))
    
    Phylo_MMR_model4 <- Almer(lnAMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.adult.mmr.er))
    Phylo_MMR_model4int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.adult.mmr.er))
    
    Phylo_MMR_model5 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.adult.mmr.er))
    Phylo_MMR_model5int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.adult.mmr.er))
    
    ### BIC --------------
    MMR_BIC<-BICdelta(BIC(
                    Phylo_MMR_model2, Phylo_MMR_model2int, #Phylo_MMR_model2.POLY, Phylo_MMR_model2intPOLY,
                    Phylo_MMR_model4, Phylo_MMR_model4int, #Phylo_MMR_model4.POLY, Phylo_MMR_model4intPOLY,
                    Phylo_MMR_model5, Phylo_MMR_model5int))#, Phylo_MMR_model5.POLY, Phylo_MMR_model5intPOLY))
    
    ## AAS optimal ------------------
    #
    Phylo_AS_model2 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.adult.aas.er))
    Phylo_AS_model2int <- Almer(lnAS ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.adult.aas.er))
    
    Phylo_AS_model4 <- Almer(lnAS ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.adult.aas.er))
    Phylo_AS_model4int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.adult.aas.er))
    
    Phylo_AS_model5 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.adult.aas.er))
    Phylo_AS_model5int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.adult.aas.er))
    
    
    ### BIC --------
    AS_BIC<-BICdelta(BIC(
    
                    Phylo_AS_model2, Phylo_AS_model2int,# Phylo_AS_model2.POLY, Phylo_AS_model2intPOLY,
                    Phylo_AS_model4, Phylo_AS_model4int, #Phylo_AS_model4.POLY, Phylo_AS_model4intPOLY,
                    Phylo_AS_model5, Phylo_AS_model5int))#, #Phylo_AS_model5.POLY, Phylo_AS_model5intPOLY))
    
    
    ## FAS optimal ----------------
    # 
    Phylo_FAS_model2 <- Almer(log(FAS) ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.adult.fas.er))
    Phylo_FAS_model2int <- Almer(log(FAS) ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.adult.fas.er))
    
    Phylo_FAS_model4 <- Almer(log(FAS) ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.adult.fas.er))
    Phylo_FAS_model4int <- Almer(log(FAS) ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.adult.fas.er))
    
    Phylo_FAS_model5 <- Almer(log(FAS) ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.adult.fas.er))
    Phylo_FAS_model5int <- Almer(log(FAS) ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.adult.fas.er))
    
    
    ### BIC -------------
    FAS_BIC<-BICdelta(BIC(
                    Phylo_FAS_model2, Phylo_FAS_model2int,# Phylo_FAS_model2.POLY, Phylo_FAS_model2intPOLY,
                    Phylo_FAS_model4, Phylo_FAS_model4int,# Phylo_FAS_model4.POLY, Phylo_FAS_model4intPOLY,
                    Phylo_FAS_model5, Phylo_FAS_model5int))#, Phylo_FAS_model5.POLY, Phylo_FAS_model5intPOLY))
    
    
    ## RMR warm ------------------
    # 
    Phylo_RMR_W_model2 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.adult.rmr.w))
    Phylo_RMR_W_model2int <- Almer(lnRMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.adult.rmr.w))
    
    Phylo_RMR_W_model4 <- Almer(lnRMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.adult.rmr.w))
    Phylo_RMR_W_model4int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.adult.rmr.w))
    
    Phylo_RMR_W_model5 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.adult.rmr.w))
    Phylo_RMR_W_model5int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.adult.rmr.w))
    
    
    
    ## BIC -------------
    RMR_W_BIC<-BICdelta(BIC(
                    Phylo_RMR_W_model2, Phylo_RMR_W_model2int,# Phylo_RMR_W_model2.POLY, Phylo_RMR_W_model2intPOLY,
                    Phylo_RMR_W_model4, Phylo_RMR_W_model4int,# Phylo_RMR_W_model4.POLY, Phylo_RMR_W_model4intPOLY,
                    Phylo_RMR_W_model5, Phylo_RMR_W_model5int))#, Phylo_RMR_W_model5.POLY, Phylo_RMR_W_model5intPOLY))
    
    # AMR / warm temps -------------
    Phylo_MMR_W_model2 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.adult.mmr.w))
    Phylo_MMR_W_model2int <- Almer(lnAMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.adult.mmr.w))
    
    Phylo_MMR_W_model4 <- Almer(lnAMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.adult.mmr.w))
    Phylo_MMR_W_model4int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.adult.mmr.w))
    
    Phylo_MMR_W_model5 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.adult.mmr.w))
    Phylo_MMR_W_model5int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.adult.mmr.w))
    
    
    ## BIC -------------
    MMR_W_BIC<-BICdelta(BIC(
                    Phylo_MMR_W_model2, Phylo_MMR_W_model2int, #Phylo_MMR_W_model2.POLY, Phylo_MMR_W_model2intPOLY,
                    Phylo_MMR_W_model4, Phylo_MMR_W_model4int, #Phylo_MMR_W_model4.POLY, Phylo_MMR_W_model4intPOLY,
                    Phylo_MMR_W_model5, Phylo_MMR_W_model5int))#, Phylo_MMR_W_model5.POLY, Phylo_MMR_W_model5intPOLY))
    
    
    
    # AS / warm temps --------------
    Phylo_AS_W_model2 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.adult.aas.w))
    Phylo_AS_W_model2int <- Almer(lnAS ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.adult.aas.w))
    
    Phylo_AS_W_model4 <- Almer(lnAS ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.adult.aas.w))
    Phylo_AS_W_model4int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.adult.aas.w))
    
    Phylo_AS_W_model5 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.adult.aas.w))
    Phylo_AS_W_model5int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.adult.aas.w))
    
    ## BIC -------------
    AS_W_BIC<-BICdelta(BIC(
                    Phylo_AS_W_model2, Phylo_AS_W_model2int,# Phylo_AS_W_model2.POLY, Phylo_AS_W_model2intPOLY,
                    Phylo_AS_W_model4, Phylo_AS_W_model4int,# Phylo_AS_W_model4.POLY, Phylo_AS_W_model4intPOLY,
                    Phylo_AS_W_model5, Phylo_AS_W_model5int))#, Phylo_AS_W_model5.POLY, Phylo_AS_W_model5intPOLY))
    
    # FAS / warm temps --------------
    Phylo_FAS_W_model2 <- Almer(lnFAS ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.adult.fas.w))
    Phylo_FAS_W_model2int <- Almer(lnFAS ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.adult.fas.w))
    
    Phylo_FAS_W_model4 <- Almer(lnFAS ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.adult.fas.w))
    Phylo_FAS_W_model4int <- Almer(lnFAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.adult.fas.w))
    
    Phylo_FAS_W_model5 <- Almer(lnFAS ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.adult.fas.w))
    Phylo_FAS_W_model5int <- Almer(lnFAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.adult.fas.w))
    
    
    ## BIC -------------
    FAS_W_BIC<-BICdelta(BIC(
                    Phylo_FAS_W_model2, Phylo_FAS_W_model2int, #Phylo_FAS_W_model2.POLY, Phylo_FAS_W_model2intPOLY,
                    Phylo_FAS_W_model4, Phylo_FAS_W_model4int, #Phylo_FAS_W_model4.POLY, Phylo_FAS_W_model4intPOLY,
                    Phylo_FAS_W_model5, Phylo_FAS_W_model5int)) #Phylo_FAS_W_model5.POLY, Phylo_FAS_W_model5intPOLY))
    
    
    ## BIC results ------
    # RMR_BIC # Phylo_RMR_model4 (dBIC 0)
    # MMR_BIC # Phylo_MMR_model4int (dBIC 0)
    # AS_BIC  # Phylo_AS_model4 (dBIC 0)
    # FAS_BIC # Phylo_FAS_model2 (dBIC 0)
    # 
    # RMR_W_BIC # Phylo_RMR_W_model5 (dBIC 0)
    # MMR_W_BIC # Phylo_MMR_W_model5 (dBIC 0)
    # AS_W_BIC  # Phylo_AS_W_model5int (dBIC 0)
    # FAS_W_BIC # Phylo_FAS_W_model5 (dBIC 0)
    
    # jan 2026 best 
    # ecol relevant
    rmr_mod_ER<- Phylo_RMR_model4
    amr_mod_ER<-Phylo_MMR_model4int
    as_mod_ER<- Phylo_AS_model4
    fas_mod_ER<- Phylo_FAS_model2
  
    # warm
    rmr_mod_W<-Phylo_RMR_W_model5
    amr_mod_W<-Phylo_MMR_W_model5
    as_mod_W<- Phylo_AS_W_model5int
    fas_mod_W<-Phylo_FAS_W_model5
    
    check_best_model(Phylo_RMR_model4, RMR_BIC)
    check_best_model(Phylo_MMR_model4int, MMR_BIC)
    check_best_model(Phylo_AS_model4, AS_BIC)
    check_best_model(Phylo_FAS_model2, FAS_BIC)
    
    check_best_model(Phylo_RMR_W_model5, RMR_W_BIC)
    check_best_model(Phylo_MMR_W_model5, MMR_W_BIC)
    check_best_model(Phylo_AS_W_model5int, AS_W_BIC)
    check_best_model(Phylo_FAS_W_model5, FAS_W_BIC)
    
    write.csv(file = here("Data_exports/BICs/RMR_BIC_adult.csv"), x = RMR_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/MMR_BIC_adult.csv"), x = MMR_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/AS_BIC_adult.csv"), x = AS_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/FAS_BIC_adult.csv"), x = FAS_BIC, row.names = T)
    
    write.csv(file = here("Data_exports/BICs/RMR_W_BIC_adult.csv"), x = RMR_W_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/MMR_W_BIC_adult.csv"), x = MMR_W_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/AS_W_BIC_adult.csv"), x = AS_W_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/FAS_W_BIC_adult.csv"), x = FAS_W_BIC, row.names = T)
    
    }
    
    if(lifestage == "juvenile"){
    message("Running juvenile data")
    # rmr
    get_phylo_matrix(species.list = unique(levels(data.rmrER$species)),
                     matrix.name = "A.juvenile",
                     dataset.ID = "RMR-optimal-juvenile")
    get_phylo_matrix(species.list = unique(levels(data.rmr.test$species)),
                     matrix.name = "A.juvenile.rmr.w",
                     dataset.ID = "RMR-warm-juvenile")
    # amr 
    get_phylo_matrix(species.list = unique(levels(data.amrER$species)),
                     matrix.name = "A.juvenile.mmr.er",
                     dataset.ID = "MMR-optimal-juvenile")
    get_phylo_matrix(species.list = unique(levels(data.amr.test$species)),
                     matrix.name = "A.juvenile.mmr.w",
                     dataset.ID = "MMR-warm-juvenile")
    # aas
    get_phylo_matrix(species.list = unique(levels(data.asER$species)),
                     matrix.name = "A.juvenile.aas.er",
                     dataset.ID = "AAS-optimal-juvenile")
    get_phylo_matrix(species.list = unique(levels(data.as.test$species)),
                     matrix.name = "A.juvenile.aas.w",
                     dataset.ID = "AAS-warm-juvenile")
    # fas
    get_phylo_matrix(species.list = unique(levels(data.fasER$species)),
                     matrix.name = "A.juvenile.fas.er",
                     dataset.ID = "FAS-optimal-juvenile")
    get_phylo_matrix(species.list = unique(levels(data.fas.test$species)),
                     matrix.name = "A.juvenile.fas.w",
                     dataset.ID = "FAS-warm-juvenile")
    
    
    ## RMR optimal -----------------
    Phylo_RMR_model2 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A.juvenile))
    Phylo_RMR_model2int <- Almer(lnRMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A.juvenile))
    
    Phylo_RMR_model4 <- Almer(lnRMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A.juvenile))
    Phylo_RMR_model4int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A.juvenile))
    
    Phylo_RMR_model5 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A.juvenile))
    Phylo_RMR_model5int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A.juvenile))
    
    
    # No singular fits 
    ### BIC rmr optimal --------------
    RMR_BIC<-BICdelta(BIC(
    
                    Phylo_RMR_model2, Phylo_RMR_model2int, #Phylo_RMR_model2.POLY, Phylo_RMR_model2intPOLY,
                    Phylo_RMR_model4, Phylo_RMR_model4int,# Phylo_RMR_model4.POLY, Phylo_RMR_model4intPOLY,
                    Phylo_RMR_model5, Phylo_RMR_model5int# Phylo_RMR_model5.POLY, Phylo_RMR_model5intPOLY
    ))
    
    ## MMR optimal -----------------
    Phylo_MMR_model2 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.juvenile.mmr.er))
    Phylo_MMR_model2int <- Almer(lnAMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.juvenile.mmr.er))
    
    Phylo_MMR_model4 <- Almer(lnAMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.juvenile.mmr.er))
    Phylo_MMR_model4int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.juvenile.mmr.er))
    
    Phylo_MMR_model5 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.juvenile.mmr.er))
    Phylo_MMR_model5int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.juvenile.mmr.er))
    
    ### BIC --------------
    MMR_BIC<-BICdelta(BIC(
                    Phylo_MMR_model2, Phylo_MMR_model2int, #Phylo_MMR_model2.POLY, Phylo_MMR_model2intPOLY,
                    Phylo_MMR_model4, Phylo_MMR_model4int, #Phylo_MMR_model4.POLY, Phylo_MMR_model4intPOLY,
                    Phylo_MMR_model5, Phylo_MMR_model5int))#, Phylo_MMR_model5.POLY, Phylo_MMR_model5intPOLY))
    
    ## AAS optimal ------------------
    #
    Phylo_AS_model2 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.juvenile.aas.er))
    Phylo_AS_model2int <- Almer(lnAS ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.juvenile.aas.er))
    
    Phylo_AS_model4 <- Almer(lnAS ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.juvenile.aas.er))
    Phylo_AS_model4int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.juvenile.aas.er))
    
    Phylo_AS_model5 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.juvenile.aas.er))
    Phylo_AS_model5int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.juvenile.aas.er))
    
    
    ### BIC --------
    AS_BIC<-BICdelta(BIC(
    
                    Phylo_AS_model2, Phylo_AS_model2int,# Phylo_AS_model2.POLY, Phylo_AS_model2intPOLY,
                    Phylo_AS_model4, Phylo_AS_model4int, #Phylo_AS_model4.POLY, Phylo_AS_model4intPOLY,
                    Phylo_AS_model5, Phylo_AS_model5int))#, #Phylo_AS_model5.POLY, Phylo_AS_model5intPOLY))
    
    
    ## FAS optimal ----------------
    # 
    Phylo_FAS_model2 <- Almer(log(FAS) ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.juvenile.fas.er))
    Phylo_FAS_model2int <- Almer(log(FAS) ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.juvenile.fas.er))
    
    Phylo_FAS_model4 <- Almer(log(FAS) ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.juvenile.fas.er))
    Phylo_FAS_model4int <- Almer(log(FAS) ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.juvenile.fas.er))
    
    Phylo_FAS_model5 <- Almer(log(FAS) ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.juvenile.fas.er))
    Phylo_FAS_model5int <- Almer(log(FAS) ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.juvenile.fas.er))
    
    
    ### BIC -------------
    FAS_BIC<-BICdelta(BIC(
                    Phylo_FAS_model2, Phylo_FAS_model2int,# Phylo_FAS_model2.POLY, Phylo_FAS_model2intPOLY,
                    Phylo_FAS_model4, Phylo_FAS_model4int,# Phylo_FAS_model4.POLY, Phylo_FAS_model4intPOLY,
                    Phylo_FAS_model5, Phylo_FAS_model5int))#, Phylo_FAS_model5.POLY, Phylo_FAS_model5intPOLY))
    
    
    ## RMR warm ------------------
    # 
    Phylo_RMR_W_model2 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.juvenile.rmr.w))
    Phylo_RMR_W_model2int <- Almer(lnRMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.juvenile.rmr.w))
    
    Phylo_RMR_W_model4 <- Almer(lnRMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.juvenile.rmr.w))
    Phylo_RMR_W_model4int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.juvenile.rmr.w))
    
    Phylo_RMR_W_model5 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.juvenile.rmr.w))
    Phylo_RMR_W_model5int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.juvenile.rmr.w))
    
    
    
    ## BIC -------------
    RMR_W_BIC<-BICdelta(BIC(
                    Phylo_RMR_W_model2, Phylo_RMR_W_model2int,# Phylo_RMR_W_model2.POLY, Phylo_RMR_W_model2intPOLY,
                    Phylo_RMR_W_model4, Phylo_RMR_W_model4int,# Phylo_RMR_W_model4.POLY, Phylo_RMR_W_model4intPOLY,
                    Phylo_RMR_W_model5, Phylo_RMR_W_model5int))#, Phylo_RMR_W_model5.POLY, Phylo_RMR_W_model5intPOLY))
    
    # AMR / warm temps -------------
    Phylo_MMR_W_model2 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.juvenile.mmr.w))
    Phylo_MMR_W_model2int <- Almer(lnAMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.juvenile.mmr.w))
    
    Phylo_MMR_W_model4 <- Almer(lnAMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.juvenile.mmr.w))
    Phylo_MMR_W_model4int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.juvenile.mmr.w))
    
    Phylo_MMR_W_model5 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.juvenile.mmr.w))
    Phylo_MMR_W_model5int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.juvenile.mmr.w))
    
    
    ## BIC -------------
    MMR_W_BIC<-BICdelta(BIC(
                    Phylo_MMR_W_model2, Phylo_MMR_W_model2int, #Phylo_MMR_W_model2.POLY, Phylo_MMR_W_model2intPOLY,
                    Phylo_MMR_W_model4, Phylo_MMR_W_model4int, #Phylo_MMR_W_model4.POLY, Phylo_MMR_W_model4intPOLY,
                    Phylo_MMR_W_model5, Phylo_MMR_W_model5int))#, Phylo_MMR_W_model5.POLY, Phylo_MMR_W_model5intPOLY))
    
    
    
    # AS / warm temps --------------
    Phylo_AS_W_model2 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.juvenile.aas.w))
    Phylo_AS_W_model2int <- Almer(lnAS ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.juvenile.aas.w))
    
    Phylo_AS_W_model4 <- Almer(lnAS ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.juvenile.aas.w))
    Phylo_AS_W_model4int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.juvenile.aas.w))
    
    Phylo_AS_W_model5 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.juvenile.aas.w))
    Phylo_AS_W_model5int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.juvenile.aas.w))
    
    ## BIC -------------
    AS_W_BIC<-BICdelta(BIC(
                    Phylo_AS_W_model2, Phylo_AS_W_model2int,# Phylo_AS_W_model2.POLY, Phylo_AS_W_model2intPOLY,
                    Phylo_AS_W_model4, Phylo_AS_W_model4int,# Phylo_AS_W_model4.POLY, Phylo_AS_W_model4intPOLY,
                    Phylo_AS_W_model5, Phylo_AS_W_model5int))#, Phylo_AS_W_model5.POLY, Phylo_AS_W_model5intPOLY))
    
    # FAS / warm temps --------------
    Phylo_FAS_W_model2 <- Almer(lnFAS ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.juvenile.fas.w))
    Phylo_FAS_W_model2int <- Almer(lnFAS ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.juvenile.fas.w))
    
    Phylo_FAS_W_model4 <- Almer(lnFAS ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.juvenile.fas.w))
    Phylo_FAS_W_model4int <- Almer(lnFAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.juvenile.fas.w))
    
    Phylo_FAS_W_model5 <- Almer(lnFAS ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.juvenile.fas.w))
    Phylo_FAS_W_model5int <- Almer(lnFAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.juvenile.fas.w))
    
    
    ## BIC -------------
    FAS_W_BIC<-BICdelta(BIC(
                    Phylo_FAS_W_model2, Phylo_FAS_W_model2int, #Phylo_FAS_W_model2.POLY, Phylo_FAS_W_model2intPOLY,
                    Phylo_FAS_W_model4, Phylo_FAS_W_model4int, #Phylo_FAS_W_model4.POLY, Phylo_FAS_W_model4intPOLY,
                    Phylo_FAS_W_model5, Phylo_FAS_W_model5int)) #Phylo_FAS_W_model5.POLY, Phylo_FAS_W_model5intPOLY))
    
    
    ## BIC results ------
    # RMR_BIC # 
    # MMR_BIC #
    # AS_BIC  # 
    # FAS_BIC # 
    # 
    # RMR_W_BIC # 
    # MMR_W_BIC # 
    # AS_W_BIC  # 
    # FAS_W_BIC #
    check_best_model(Phylo_RMR_model5, RMR_BIC)
    check_best_model(Phylo_MMR_model4int, MMR_BIC)
    check_best_model(Phylo_AS_model2int, AS_BIC)
    check_best_model(Phylo_FAS_model2int, FAS_BIC)
    
    check_best_model(Phylo_RMR_W_model2, RMR_W_BIC)
    check_best_model(Phylo_MMR_W_model4, MMR_W_BIC)
    check_best_model(Phylo_AS_W_model2, AS_W_BIC)
    check_best_model(Phylo_FAS_W_model2, FAS_W_BIC)
    # jan 2026 best
    # ecol relevant
    rmr_mod_ER<- Phylo_RMR_model5
    amr_mod_ER<-Phylo_MMR_model4int
    as_mod_ER<-Phylo_AS_model2int
    fas_mod_ER<-Phylo_FAS_model2int
  
    # warm
    rmr_mod_W<-Phylo_RMR_W_model2
    amr_mod_W<-Phylo_MMR_W_model4
    as_mod_W<-Phylo_AS_W_model2
    fas_mod_W<-Phylo_FAS_W_model2
    
    write.csv(file = here("Data_exports/BICs/RMR_BIC_juvenile.csv"), x = RMR_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/MMR_BIC_juvenile.csv"), x = MMR_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/AS_BIC_juvenile.csv"), x = AS_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/FAS_BIC_juvenile.csv"), x = FAS_BIC, row.names = T)
    
    write.csv(file = here("Data_exports/BICs/RMR_W_BIC_juvenile.csv"), x = RMR_W_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/MMR_W_BIC_juvenile.csv"), x = MMR_W_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/AS_W_BIC_juvenile.csv"), x = AS_W_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/FAS_W_BIC_juvenile.csv"), x = FAS_W_BIC, row.names = T)
    
    } 
    
  } else { # full dataset
    
    # rmr
    get_phylo_matrix(species.list = unique(levels(data.rmrER$species)),
                     matrix.name = "A",
                     dataset.ID = "RMR-optimal")
    get_phylo_matrix(species.list = unique(levels(data.rmr.test$species)),
                     matrix.name = "A.rmr.w",
                     dataset.ID = "RMR-warm")
    # amr 
    get_phylo_matrix(species.list = unique(levels(data.amrER$species)),
                     matrix.name = "A.mmr.er",
                     dataset.ID = "MMR-optimal")
    get_phylo_matrix(species.list = unique(levels(data.amr.test$species)),
                     matrix.name = "A.mmr.w",
                     dataset.ID = "MMR-warm")
    # aas
    get_phylo_matrix(species.list = unique(levels(data.asER$species)),
                     matrix.name = "A.aas.er",
                     dataset.ID = "AAS-optimal")
    get_phylo_matrix(species.list = unique(levels(data.as.test$species)),
                     matrix.name = "A.aas.w",
                     dataset.ID = "AAS-warm")
    # fas
    get_phylo_matrix(species.list = unique(levels(data.fasER$species)),
                     matrix.name = "A.fas.er",
                     dataset.ID = "FAS-optimal")
    get_phylo_matrix(species.list = unique(levels(data.fas.test$species)),
                     matrix.name = "A.fas.w",
                     dataset.ID = "FAS-warm")
  
  
    ## RMR optimal -----------------
    Phylo_RMR_model2 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
    Phylo_RMR_model2int <- Almer(lnRMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))

    Phylo_RMR_model4 <- Almer(lnRMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
    Phylo_RMR_model4int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))

    Phylo_RMR_model5 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
    Phylo_RMR_model5int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))

  
    # No singular fits 
    ### BIC rmr optimal --------------
    RMR_BIC<-BICdelta(BIC(
  
                    Phylo_RMR_model2, Phylo_RMR_model2int, #Phylo_RMR_model2.POLY, Phylo_RMR_model2intPOLY,
                    Phylo_RMR_model4, Phylo_RMR_model4int,# Phylo_RMR_model4.POLY, Phylo_RMR_model4intPOLY,
                    Phylo_RMR_model5, Phylo_RMR_model5int# Phylo_RMR_model5.POLY, Phylo_RMR_model5intPOLY
  ))
    
    ## MMR optimal -----------------
    Phylo_MMR_model2 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
    Phylo_MMR_model2int <- Almer(lnAMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))

    Phylo_MMR_model4 <- Almer(lnAMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
    Phylo_MMR_model4int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  
    Phylo_MMR_model5 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
    Phylo_MMR_model5int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))

    ### BIC --------------
    MMR_BIC<-BICdelta(BIC(
                    Phylo_MMR_model2, Phylo_MMR_model2int, #Phylo_MMR_model2.POLY, Phylo_MMR_model2intPOLY,
                    Phylo_MMR_model4, Phylo_MMR_model4int, #Phylo_MMR_model4.POLY, Phylo_MMR_model4intPOLY,
                    Phylo_MMR_model5, Phylo_MMR_model5int))#, Phylo_MMR_model5.POLY, Phylo_MMR_model5intPOLY))
    
    ## AAS optimal ------------------
    #
    Phylo_AS_model2 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
    Phylo_AS_model2int <- Almer(lnAS ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))

    Phylo_AS_model4 <- Almer(lnAS ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
    Phylo_AS_model4int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
 
    Phylo_AS_model5 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
    Phylo_AS_model5int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))

  
    ### BIC --------
    AS_BIC<-BICdelta(BIC(
  
                    Phylo_AS_model2, Phylo_AS_model2int,# Phylo_AS_model2.POLY, Phylo_AS_model2intPOLY,
                    Phylo_AS_model4, Phylo_AS_model4int, #Phylo_AS_model4.POLY, Phylo_AS_model4intPOLY,
                    Phylo_AS_model5, Phylo_AS_model5int))#, #Phylo_AS_model5.POLY, Phylo_AS_model5intPOLY))
    
    
    ## FAS optimal ----------------
    # 
    Phylo_FAS_model2 <- Almer(log(FAS) ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
    Phylo_FAS_model2int <- Almer(log(FAS) ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
 
    Phylo_FAS_model4 <- Almer(log(FAS) ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
    Phylo_FAS_model4int <- Almer(log(FAS) ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))

    Phylo_FAS_model5 <- Almer(log(FAS) ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
    Phylo_FAS_model5int <- Almer(log(FAS) ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))

  
    ### BIC -------------
    FAS_BIC<-BICdelta(BIC(
                    Phylo_FAS_model2, Phylo_FAS_model2int,# Phylo_FAS_model2.POLY, Phylo_FAS_model2intPOLY,
                    Phylo_FAS_model4, Phylo_FAS_model4int,# Phylo_FAS_model4.POLY, Phylo_FAS_model4intPOLY,
                    Phylo_FAS_model5, Phylo_FAS_model5int))#, Phylo_FAS_model5.POLY, Phylo_FAS_model5intPOLY))
    
    
    ## RMR warm ------------------
    # 
    Phylo_RMR_W_model2 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
    Phylo_RMR_W_model2int <- Almer(lnRMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))

    Phylo_RMR_W_model4 <- Almer(lnRMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
    Phylo_RMR_W_model4int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))

    Phylo_RMR_W_model5 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
    Phylo_RMR_W_model5int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))

  
  
    ## BIC -------------
    RMR_W_BIC<-BICdelta(BIC(
                    Phylo_RMR_W_model2, Phylo_RMR_W_model2int,# Phylo_RMR_W_model2.POLY, Phylo_RMR_W_model2intPOLY,
                    Phylo_RMR_W_model4, Phylo_RMR_W_model4int,# Phylo_RMR_W_model4.POLY, Phylo_RMR_W_model4intPOLY,
                    Phylo_RMR_W_model5, Phylo_RMR_W_model5int))#, Phylo_RMR_W_model5.POLY, Phylo_RMR_W_model5intPOLY))
    
    # AMR / warm temps -------------
    Phylo_MMR_W_model2 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
    Phylo_MMR_W_model2int <- Almer(lnAMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))

    Phylo_MMR_W_model4 <- Almer(lnAMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
    Phylo_MMR_W_model4int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))

    Phylo_MMR_W_model5 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
    Phylo_MMR_W_model5int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))

   
    ## BIC -------------
    MMR_W_BIC<-BICdelta(BIC(
                    Phylo_MMR_W_model2, Phylo_MMR_W_model2int, #Phylo_MMR_W_model2.POLY, Phylo_MMR_W_model2intPOLY,
                    Phylo_MMR_W_model4, Phylo_MMR_W_model4int, #Phylo_MMR_W_model4.POLY, Phylo_MMR_W_model4intPOLY,
                    Phylo_MMR_W_model5, Phylo_MMR_W_model5int))#, Phylo_MMR_W_model5.POLY, Phylo_MMR_W_model5intPOLY))
    
    
    
    # AS / warm temps --------------
    Phylo_AS_W_model2 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
    Phylo_AS_W_model2int <- Almer(lnAS ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
 
    Phylo_AS_W_model4 <- Almer(lnAS ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
    Phylo_AS_W_model4int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))

    Phylo_AS_W_model5 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
    Phylo_AS_W_model5int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))

    ## BIC -------------
    AS_W_BIC<-BICdelta(BIC(
                    Phylo_AS_W_model2, Phylo_AS_W_model2int,# Phylo_AS_W_model2.POLY, Phylo_AS_W_model2intPOLY,
                    Phylo_AS_W_model4, Phylo_AS_W_model4int,# Phylo_AS_W_model4.POLY, Phylo_AS_W_model4intPOLY,
                    Phylo_AS_W_model5, Phylo_AS_W_model5int))#, Phylo_AS_W_model5.POLY, Phylo_AS_W_model5intPOLY))
    
    # FAS / warm temps --------------
    Phylo_FAS_W_model2 <- Almer(lnFAS ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
    Phylo_FAS_W_model2int <- Almer(lnFAS ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))

    Phylo_FAS_W_model4 <- Almer(lnFAS ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
    Phylo_FAS_W_model4int <- Almer(lnFAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))

    Phylo_FAS_W_model5 <- Almer(lnFAS ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
    Phylo_FAS_W_model5int <- Almer(lnFAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))

    
    ## BIC -------------
    FAS_W_BIC<-BICdelta(BIC(
                    Phylo_FAS_W_model2, Phylo_FAS_W_model2int, #Phylo_FAS_W_model2.POLY, Phylo_FAS_W_model2intPOLY,
                    Phylo_FAS_W_model4, Phylo_FAS_W_model4int, #Phylo_FAS_W_model4.POLY, Phylo_FAS_W_model4intPOLY,
                    Phylo_FAS_W_model5, Phylo_FAS_W_model5int)) #Phylo_FAS_W_model5.POLY, Phylo_FAS_W_model5intPOLY))
    
    
    ## BIC results ------
    # RMR_BIC #
    # MMR_BIC #
    # AS_BIC #
    # FAS_BIC #
    # 
    # RMR_W_BIC #
    # MMR_W_BIC #
    # AS_W_BIC  #
    # FAS_W_BIC #
    #
    # model checks: 
    check_best_model(Phylo_RMR_model4int, RMR_BIC)
    check_best_model(Phylo_MMR_model4, MMR_BIC)
    check_best_model(Phylo_AS_model5int, AS_BIC)
    check_best_model(Phylo_FAS_model2int, FAS_BIC)
    
    check_best_model(Phylo_RMR_W_model2, RMR_W_BIC)
    check_best_model(Phylo_MMR_W_model2int, MMR_W_BIC)
    check_best_model(Phylo_AS_W_model2int, AS_W_BIC)
    check_best_model(Phylo_FAS_W_model2, FAS_W_BIC)
    
    # jan 2026 best 
    # ecol relevant
    rmr_mod_ER<- Phylo_RMR_model4int
    amr_mod_ER<-Phylo_MMR_model4
    as_mod_ER<- Phylo_AS_model5int
    fas_mod_ER<- Phylo_FAS_model2int
  
    # warm
    rmr_mod_W<-Phylo_RMR_W_model2
    amr_mod_W<-Phylo_MMR_W_model2int
    as_mod_W<- Phylo_AS_W_model2int
    fas_mod_W<-Phylo_FAS_W_model2
    

    write.csv(file = here("Data_exports/BICs/RMR_BIC.csv"), x = RMR_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/MMR_BIC.csv"), x = MMR_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/AS_BIC.csv"), x = AS_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/FAS_BIC.csv"), x = FAS_BIC, row.names = T)
  
    write.csv(file = here("Data_exports/BICs/RMR_W_BIC.csv"), x = RMR_W_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/MMR_W_BIC.csv"), x = MMR_W_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/AS_W_BIC.csv"), x = AS_W_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/FAS_W_BIC.csv"), x = FAS_W_BIC, row.names = T)
    
  }
    

  # *****************************************************************
  # *****************************************************************
  # un-comment to see best model summaries or view output from R markdown
  # # Ecol relev
  # summary(rmr_mod_ER) # one scaling 0.823, temp 0.063
  # summary(amr_mod_ER) # (poly) temp dependent scaling (poly)
  # summary(fas_mod_ER) # (linear) temp dependent scaling
  # summary(as_mod_ER) # (poly) temp dependent scaling
  # 
  # # # warm
  # summary(rmr_mod_W)# one scaling 0.89, poly temp 
  # summary(amr_mod_W) # one scaling 0.78, poly temp 
  # summary(fas_mod_W) # one scaling -0.07, poly temp 
  # summary(as_mod_W) # one scaling 0.79, poly temp 
  
  # # # Ecol relev
  # plot(rmr_mod_ER) # residuals
  # plot(amr_mod_ER) # residuals
  # plot(fas_mod_ER) # residuals
  # plot(as_mod_ER) # residuals
  # # # warm
  # plot(rmr_mod_W) # residuals
  # plot(amr_mod_W) # residuals
  # plot(fas_mod_W) # residuals
  # plot(as_mod_W) # residuals
  # # 
  # hist(resid(amr_mod_ER), breaks = 50) # residuals
  # hist(resid(rmr_mod_ER), breaks = 50) # residuals
  # hist(resid(fas_mod_ER), breaks = 50) # residuals
  # hist(resid(as_mod_ER), breaks = 50) # residuals
  # # 
  # hist(resid(amr_mod_W), breaks = 50) # residuals
  # hist(resid(rmr_mod_W), breaks = 50) # residuals
  # hist(resid(fas_mod_W), breaks = 50) # residuals
  
  print(car::Anova(rmr_mod_ER))
  print(car::Anova(amr_mod_ER))
  print(car::Anova(as_mod_ER))
  print(car::Anova(fas_mod_ER))

  print(car::Anova(rmr_mod_W))
  print(car::Anova(amr_mod_W))
  print(car::Anova(as_mod_W))
  print(car::Anova(fas_mod_W))
  
  # hist(resid(as_mod_W), breaks = 50) # little skew not too bad, only 5 measurements, all reasonable biologically

  # *****************************************************************************
  # *****************************************************************************
  # Model scaling parameters and CIs ----------
  # custom function to obtain model parameters and recalculated mass specific MR using model estimate scaling slopes
  # this also expands the dataset to get mass independent values of metabolic rates
  # function source script available at 'mixed_model_outputs.R'
  # tic() # time the function run
  model_out<-get_model_outputs(
                best.model.rmr.er = rmr_mod_ER,
                best.model.amr.er= amr_mod_ER,
                best.model.as.er = as_mod_ER,
                best.model.fas.er = fas_mod_ER,
                best.model.rmr.w = rmr_mod_W,
                best.model.amr.w = amr_mod_W,
                best.model.as.w = as_mod_W,
                best.model.fas.w = fas_mod_W,
                data.rmr.test = data.rmr.test,
                data.rmrER = data.rmrER,
                data.amr.test = data.amr.test,
                data.amrER = data.amrER,                        
                data.as.test = data.as.test,
                data.asER = data.asER, 
                data.fas.test = data.fas.test,
                data.fasER = data.fasER, 
                name_extension = paste(lifestage,"models", sep="")) #
  # toc()

  sum_CItable<-data.frame(model_out[[1]])
  colnames(sum_CItable)<-c("ci5", "ci95","var_repeat", "MR", "temp_cat")
  sum_data<-data.frame(model_out[[2]])
  scaling_params<-data.frame(model_out[[3]])

  # *****************************************************************************
  # *****************************************************************************

  # Figures -------
  set.seed(51423)
  
  # General scaling plots ------
  scaling_params_amr_w<-scaling_params %>% 
    filter(performance == "MMR",
           temp_categ =="warm",
           Temperature == 20) # one slope
  scaling_params_amr_er<-scaling_params %>% 
    filter(performance == "MMR",
           temp_categ =="ecol_relev",
           Temperature == 20) # change with temp
  
  AMRmodel_plot1<- ggplot(data=data.amrER, aes(x=lnBWg, y=lnAMR)) +
    geom_point(alpha=0.9,  size=1, pch=1, color="grey75")+
    geom_point(data=data.amr.test, aes(x=lnBWg, y=lnAMR),
               alpha=0.9,  size=1, pch=21, show.legend = FALSE, stroke =0.2,
               fill = cols.amr[3], color = cols.amr[3])+
    geom_segment(aes(x = 0, xend = 9,
                   y = scaling_params_amr_w$Intercept + scaling_params_amr_w$Slope * 0,
                   yend = scaling_params_amr_w$Intercept + scaling_params_amr_w$Slope * 9),
                color = cols.amr[1])+
    geom_segment(aes(x = -2.5, xend = 9,
                   y = scaling_params_amr_er$Intercept + scaling_params_amr_er$Slope * -2.5,
                   yend = scaling_params_amr_er$Intercept + scaling_params_amr_er$Slope * 9),
                color = "black")+
    annotate("text",  x = -5.2, y = 11.5,
             label = bquote(Optimal:~italic(b)[MMR] == .(round(scaling_params_amr_er$Slope,2))), 
                            # "\u00b1" ~ .(round(scaling_params_amr_er$SE_slope,2))),
             size=4, hjust=0, family="Helvetica", color = "black")+
    annotate("text",  x = -5.2, y = 9.8,
             label = bquote(Warm:~italic(b)[MMR] == .(round(scaling_params_amr_w$Slope,2))),
             size=4, hjust=0, family="Helvetica", color = cols.amr[1])+
    scale_color_gradient( low = "grey", high = "black")+
    ylim(x = -6.3, 12)+
    xlim(x = -6.3, 12)+
    annotate("text", label = paste("n = ", nrow(data.amrER), sep=""),
             x = -5.2, y = 8.5, size=3, hjust=0, family="Helvetica", color = "black")+
    annotate("text", label = paste("n = ", nrow(data.amr.test), sep=""),
             x = -5.2, y = 7.6, size=3, hjust=0, family="Helvetica", color = cols.amr[1])+
    annotate("text", label = lifestage,
             x = 5, y = -5, size=5, hjust=0, family="Helvetica", color = "black")
  
    if(is.null(lifestage)){
      AMRmodel_plot1<-AMRmodel_plot1+ annotate("text",  x = 7.5, y = 9.8,
             label = bquote("\u2193" ~ "with" ~ degree*C),
             size=4, hjust=0,  color = cols.amr[1])
    }else{
      if(lifestage == "juvenile"){
        AMRmodel_plot1<-AMRmodel_plot1+ annotate("text",  x = 7.5, y = 11.5,
           label = bquote("\u2193" ~ "with" ~ degree*C),
           size=4, hjust=0,  color = "black")
      }
      if(lifestage == "adult"){
        AMRmodel_plot1<-AMRmodel_plot1+ annotate("text",  x = 7.5, y = 11.5,
           label = bquote("\u2191" ~ "with" ~ degree*C),
           size=4, hjust=0,  color = "black")
      }
    }
  
  ggformat(AMRmodel_plot1, x_title=expression(italic(ln)*Body~weight~(g)),
           y_title=expression(italic(ln)*MMR~(mg~O[2]~h^-1)), size_text = 12, print = T)
  
  
  # RMR 
  scaling_params_rmr_w<-scaling_params %>%
    filter(performance == "RMR",
           temp_categ =="warm",
           Temperature == 20) # one slope
  scaling_params_rmr_er<-scaling_params %>% 
    filter(performance == "RMR",
           temp_categ =="ecol_relev",
           Temperature == 20) # one slope
  
  RMRmodel_plot1<- ggplot(data=data.rmrER, aes(x=lnBWg, y=lnRMR)) +
    geom_point(alpha=0.9,  size=1, pch=1, color="grey75")+
    geom_point(data=data.rmr.test, aes(x=lnBWg, y=lnRMR),
               alpha=0.9,  size=1, pch=21, show.legend = FALSE, stroke =0.2,
               fill = cols.rmr[2], color = cols.rmr[2])+
    geom_segment(aes(x = 0, xend = 9,
                   y = scaling_params_rmr_w$Intercept + scaling_params_rmr_w$Slope * 0,
                   yend = scaling_params_rmr_w$Intercept + scaling_params_rmr_w$Slope * 9),
                color = cols.rmr[1])+
    geom_segment(aes(x = -2.5, xend = 9,
                   y = scaling_params_rmr_er$Intercept + scaling_params_rmr_er$Slope * -2.5,
                   yend = scaling_params_rmr_er$Intercept + scaling_params_rmr_er$Slope * 9),
                color = "black")+
    annotate("text",  x = -5.2, y = 11.5,
             label = bquote(Optimal:~italic(b)[RMR] == .(round(scaling_params_rmr_er$Slope,2))),
             size=4, hjust=0, family="Helvetica", color = "black")+
    annotate("text",  x = -5.2, y = 9.8,
             label = bquote(Warm:~italic(b)[RMR] == .(round(scaling_params_rmr_w$Slope,2))),
             size=4, hjust=0, family="Helvetica", color = cols.rmr[1])+
    scale_color_gradient( low = "grey", high = "black")+
    ylim(x = -6.3, 12)+
    xlim(x = -6.3, 12)+
    annotate("text", label = paste("n = ", nrow(data.rmrER), sep=""),
             x = -5.2, y = 8.5, size=3, hjust=0, family="Helvetica", color = "black")+
    annotate("text", label = paste("n = ", nrow(data.rmr.test), sep=""),
             x = -5.2, y = 7.6, size=3, hjust=0, family="Helvetica", color = cols.rmr[1])
  
    if(is.null(lifestage)){
      RMRmodel_plot1<-RMRmodel_plot1+ annotate("text",  x = 7.5, y = 11.5,
             label = bquote("\u2191" ~ "with" ~ degree*C),
             size=4, hjust=0,  color = "black")
    }else{

    }
    # annotate("text", label = "ADULT FISH",
             # x = 5, y = -5, size=5, hjust=0, family="Helvetica", color = "black")
  ggformat(RMRmodel_plot1, x_title=expression(italic(ln)*Body~weight~(g)),
           y_title=expression(italic(ln)*RMR~(mg~O[2]~h^-1)), size_text = 12, print = T)
  
  
  # AS 
  scaling_params_as_w<-scaling_params %>%
    filter(performance == "AS",
           temp_categ =="warm",
           Temperature == 20) # change with temp POLY
  scaling_params_as_er<-scaling_params %>% 
    filter(performance == "AS",
           temp_categ =="ecol_relev",
           Temperature == 20) # 
  
  ASmodel_plot1<- ggplot(data=data.asER, aes(x=lnBWg, y=lnAS)) +
    geom_point(alpha=0.9,  size=1, pch=1, color="grey75")+
    geom_point(data=data.as.test, aes(x=lnBWg, y=lnAS),
               alpha=0.9,  size=1, pch=21, show.legend = FALSE, stroke =0.2,
               fill = cols.as[3], color = cols.as[3])+
    geom_segment(aes(x = 0, xend = 9,
                   y = scaling_params_as_w$Intercept + scaling_params_as_w$Slope * 0,
                   yend = scaling_params_as_w$Intercept + scaling_params_as_w$Slope * 9),
                color = cols.as[1])+
    geom_segment(aes(x = -2.5, xend = 9,
                   y = scaling_params_as_er$Intercept + scaling_params_as_er$Slope * -2.5,
                   yend = scaling_params_as_er$Intercept + scaling_params_as_er$Slope * 9),
                color = "black")+
    annotate("text",  x = -5.2, y = 11.5,
             label = bquote(Optimal:~italic(b)[AS] == .(round(scaling_params_as_er$Slope,2))),
             size=4, hjust=0, family="Helvetica", color = "black")+
    annotate("text",  x = -5.2, y = 9.8,
             label = bquote(Warm:~italic(b)[AS] == .(round(scaling_params_as_w$Slope,2))),
             size=4, hjust=0, family="Helvetica", color = cols.as[1])+
    scale_color_gradient( low = "grey", high = "black")+
    ylim(x = -6.3, 12)+
    xlim(x = -6.3, 12)+
    annotate("text", label = paste("n = ", nrow(data.asER), sep=""),
             x = -5.2, y = 8.5, size=3, hjust=0, family="Helvetica", color = "black")+
    annotate("text", label = paste("n = ", nrow(data.as.test), sep=""),
             x = -5.2, y = 7.6, size=3, hjust=0, family="Helvetica", color = cols.as[1])
    if(is.null(lifestage)){
      ASmodel_plot1<-ASmodel_plot1+ annotate("text",  x = 7.5, y = 11.5,
             label = bquote("\u2193" ~ "with" ~ degree*C),
             size=4, hjust=0,  color = "black")
      ASmodel_plot1<-ASmodel_plot1+ annotate("text",  x = 7, y = 9.8,
             label = bquote("\u2193" ~ "with" ~ degree*C),
             size=4, hjust=0,  color = cols.as[1])
    }else{
      if(lifestage == "juvenile"){
        ASmodel_plot1<-ASmodel_plot1+ annotate("text",  x = 7, y = 9.8,
             label = bquote("\u2193" ~ "with" ~ degree*C),
             size=4, hjust=0,  color = cols.as[1])
      }
      if(lifestage == "adult"){
        ASmodel_plot1<-ASmodel_plot1+ annotate("text",  x = 7, y = 9.8,
             label = bquote("\u2191" ~ "with" ~ degree*C),
             size=4, hjust=0,  color = cols.as[1])
      }

    }
    # annotate("text", label = "ADULT FISH",
    #          x = 5, y = -5, size=5, hjust=0, family="Helvetica", color = "black")
  ggformat(ASmodel_plot1, x_title=expression(italic(ln)*Body~weight~(g)),
           y_title=expression(italic(ln)*AS~(mg~O[2]~h^-1)), size_text = 12, print = T)
  
  # FAS
  scaling_params_fas_w<-scaling_params %>%
    filter(performance == "FAS",
           temp_categ =="warm",
           Temperature == 20) # change with temp POLY
  scaling_params_fas_er<-scaling_params %>% 
    filter(performance == "FAS",
           temp_categ =="ecol_relev",
           Temperature == 20) # 
  
  FASmodel_plot1<- ggplot(data=data.fasER, aes(x=lnBWg, y=lnFAS)) +
    geom_point(alpha=0.9,  size=1, pch=1, color="grey75")+
    geom_point(data=data.fas.test, aes(x=lnBWg, y=lnFAS),
               alpha=0.9,  size=1, pch=21, show.legend = FALSE, stroke =0.2,
               fill = cols.fas[3], color = cols.fas[3])+
    geom_segment(aes(x = 0, xend = 9,
                   y = scaling_params_fas_w$Intercept + scaling_params_fas_w$Slope * 0,
                   yend = scaling_params_fas_w$Intercept + scaling_params_fas_w$Slope * 9),
                color = cols.fas[1])+
    geom_segment(aes(x = -2.5, xend = 9,
                   y = scaling_params_fas_er$Intercept + scaling_params_fas_er$Slope * -2.5,
                   yend = scaling_params_fas_er$Intercept + scaling_params_fas_er$Slope * 9),
                color = "black")+
    annotate("text",  x = -5.2, y = 11.5,
             label = bquote(Optimal:~italic(b)[FAS] == .(round(scaling_params_fas_er$Slope,2))),
             size=4, hjust=0, family="Helvetica", color = "black")+
    annotate("text",  x = -5.2, y = 9.8,
             label = bquote(Warm:~italic(b)[FAS] == .(round(scaling_params_fas_w$Slope,2))),
             size=4, hjust=0, family="Helvetica", color = cols.fas[1])+
    scale_color_gradient( low = "grey", high = "black")+
    ylim(x = -6.3, 12)+
    xlim(x = -6.3, 12)+
    annotate("text", label = paste("n = ", nrow(data.fasER), sep=""),
             x = -5.2, y = 8.5, size=3, hjust=0, family="Helvetica", color = "black")+
    annotate("text", label = paste("n = ", nrow(data.fas.test), sep=""),
             x = -5.2, y = 7.6, size=3, hjust=0, family="Helvetica", color = cols.fas[1])
    if(is.null(lifestage)){
      FASmodel_plot1<-FASmodel_plot1+ annotate("text",  x = 7.5, y = 11.5,
             label = bquote("\u2193" ~ "with" ~ degree*C),
             size=4, hjust=0,  color = "black")
    }else{
      if(lifestage=="juvenile"){
        FASmodel_plot1<-FASmodel_plot1+ annotate("text",  x = 7.5, y = 11.5,
             label = bquote("\u2193" ~ "with" ~ degree*C),
             size=4, hjust=0,  color = "black")
      }
    }
  
  ggformat(FASmodel_plot1, x_title=expression(italic(ln)*Body~weight~(g)),
           y_title=expression(italic(ln)*FAS~(mg~O[2]~h^-1)), size_text = 12, print = T)
  
  
  scaling<-cowplot:::plot_grid(AMRmodel_plot1, RMRmodel_plot1,
                              ASmodel_plot1, FASmodel_plot1,
                              align = "hv",
                              axis = "l",
                              nrow = 2,
                              ncol = 2,
                              labels = "AUTO",
                              label_x = c(0.18, 0.18),
                              label_y = c(0.895, 0.895),
                              label_size = 12)
  # scaling
  ggsave(filename = paste("./Figures/Figure3", lifestage,".png", sep=""),
         plot=scaling, width = 6.8, height = 6.8, units = "in")
  
  
  # Global outcomes by scaled MR ~ temprature  --------
  AMRmodel_plot1_t<-ggplot()+
    geom_point(data=data.amrER,
               aes(x=tempTest, y=AMR/(BW_g^scaling_params_amr_er$Slope), size = BW_g),
               alpha=0.9, pch=1, size = 1,color="grey75")+
    geom_point(data=data.amr.test,
               mapping = aes(x=tempTest, y=AMR/(BW_g^scaling_params_amr_w$Slope), size = BW_g),
               alpha=0.9, pch=21, size = 1,show.legend = FALSE, stroke =0.2,
               fill = cols.amr[3], color = cols.amr[1])+
    ylim(0,8)+
    xlim(0,40)+
    annotate("text", label = lifestage,
             x = 0, y = 6.5, size=5, hjust=0, family="Helvetica", color = "black")
  ggformat(AMRmodel_plot1_t, x_title="Temperature ºC",
           y_title=bquote(MMR~(mg~O[2]~h^-1~g^italic(b))),
           size_text = 12, print = T)
  
  
  RMRmodel_plot1_t<-ggplot()+
    geom_point(data=data.rmrER,
               aes(x=tempTest, y=RMR/(BW_g^scaling_params_rmr_er$Slope), size = BW_g),
               alpha=0.9, pch=1,size = 1, color="grey75")+
    geom_point(data=data.rmr.test,
               mapping = aes(x=tempTest, y=RMR/(BW_g^scaling_params_rmr_w$Slope), size = BW_g),
               alpha=0.9, pch=21,size = 1, show.legend = FALSE, stroke =0.2,
               fill = cols.rmr[3], color = cols.rmr[1])+
    ylim(0,1.5)+
    xlim(-1,40)
  ggformat(RMRmodel_plot1_t, x_title="Temperature ºC",
           y_title=bquote(RMR~(mg~O[2]~h^-1~g^italic(b))),
           size_text = 12, print = T)
  
  
  ASmodel_plot1_t<-ggplot()+
    geom_point(data=data.asER,
               aes(x=tempTest, y=AS/(BW_g^scaling_params_as_er$Slope), size = BW_g),
               alpha=0.9, pch=1, size = 1,color="grey75")+
    geom_point(data=data.as.test,
               mapping = aes(x=tempTest, y=AS/(BW_g^scaling_params_as_w$Slope), size = BW_g),
               alpha=0.9, pch=21, size = 1,show.legend = FALSE, stroke =0.2,
               fill = cols.as[3], color = cols.as[1])+
    ylim(0,3)+
    xlim(0,35)
  ggformat(ASmodel_plot1_t, x_title="Temperature ºC",
           y_title=bquote(AS~(mg~O[2]~h^-1~g^italic(b))),
           size_text = 12, print = T)
  
  FASmodel_plot1_t<-ggplot()+
    geom_point(data=data.fasER,
               aes(x=tempTest, y=FAS/(BW_g^scaling_params_fas_er$Slope), size = BW_g),
               alpha=0.9, size = 1, pch=1, color="grey75")+
    geom_point(data=data.fas.test,
               mapping = aes(x=tempTest, y=FAS/(BW_g^scaling_params_fas_w$Slope), size = BW_g),
               alpha=0.9, pch=21,size = 1, show.legend = FALSE, stroke =0.2,
               fill = cols.fas[3], color = cols.fas[1])+
    ylim(0,30)+
    xlim(0,35)
  ggformat(FASmodel_plot1_t, x_title="Temperature ºC",
           y_title=bquote(FAS~(mg~O[2]~h^-1~g^italic(b))),
           size_text = 12, print = T)
  
  scaling_t<-cowplot:::plot_grid(AMRmodel_plot1_t, RMRmodel_plot1_t,
                               ASmodel_plot1_t, FASmodel_plot1_t,
                                align = "hv",
                                axis = "l",
                                nrow = 2,
                                ncol = 2,
                                labels = "AUTO",
                               label_x = c(0.2, 0.2),
                               label_y = c(0.895, 0.895),
                                label_size = 12)
  # scaling
  ggsave(filename = paste("./Figures/Figure3-temp", lifestage,".png", sep=""),
         plot=scaling_t, width = 6.8, height = 6.8, units = "in")
  
  
  # **************************************************************
  # **************************************************************
  # ## residuals correlate with any methodological metrics? --------
  # type of MMR protocol 
  if(is.null(lifestage)){
    data.amrER$resid<-resid(amr_mod_ER)
    data.rmrER$resid<-resid(rmr_mod_ER)
    data.amr.test$resid<-resid(amr_mod_W)
    data.rmr.test$resid<-resid(rmr_mod_W)
    
    # MMR
    MMR_method_er<-ggplot(data = data.amrER ,aes( x = MMR_method, y = resid))+
      geom_point(size = 1, position = position_jitter(width = 0.2), 
                 alpha = 0.3, color  = "grey")+
      geom_boxplot(alpha = 0)
    ggformat(MMR_method_er, x_title="MMR method",
             y_title=bquote(residual~MMR),
             size_text = 12, print = F, title = "Ecol Relev")
    
    # warm temps
    MMR_method_w<-ggplot(data = data.amr.test ,aes( x = MMR_method, y = resid))+
      geom_point(size = 1, position = position_jitter(width = 0.2), 
                 alpha = 0.3, color  = cols.amr[3])+
      geom_boxplot(alpha = 0)
    ggformat(MMR_method_w, x_title="MMR method",
             y_title=bquote(residual~MMR),
             size_text = 12, print = F, title = "Warm")
    
    # hours starved
    # RMR 
    RMR_method_starved<-ggplot(data = data.rmrER, aes(x = h_starved, y = resid))+
      geom_smooth(color = "black", method = "loess", se = F)+
      geom_smooth(data = data.rmr.test, method = "loess", se = F, 
                  mapping = aes(x = h_starved, y = resid),
                  color  = cols.rmr[3], fill  = cols.rmr[3])+
      geom_point(size = 1, alpha = 1, color  = "grey")+
      geom_point(data = data.rmr.test, mapping = aes(x = h_starved, y = resid),
                 size = 1, alpha = 1, color  = cols.rmr[3])+
      geom_hline(yintercept = 0)

    ggformat(RMR_method_starved, x_title="Hours starved prior test",
             y_title=bquote(residual~RMR),
             size_text = 12, print = F, title = "All fish together")
    
    # MMR 
    MMR_method_starved<-ggplot(data = data.amrER, aes(x = h_starved, y = resid))+
      geom_smooth(color = "black", method = "loess", se = F)+
      geom_smooth(data = data.amr.test, method = "loess", se = F, 
                  mapping = aes(x = h_starved, y = resid),
                  color  = cols.amr[3], fill  = cols.amr[3])+
      geom_point(size = 1, alpha = 0.3, color  = "grey")+
      geom_point(data = data.amr.test, mapping = aes(x = h_starved, y = resid),
                 size = 1, alpha = 1, color  = cols.amr[3])+
      geom_hline(yintercept = 0)
    ggformat(MMR_method_starved, x_title="Hours starved prior test",
             y_title=bquote(residual~MMR),
             size_text = 12, print = F, title = "All fish together")
    
    # hours in respo 
    # RMR 
    RMR_method_in_respo<-ggplot(data = data.rmrER, aes(x = h_in_respo, y = resid))+
      geom_smooth(color = "black", method = "loess", se = F)+
      geom_smooth(data = data.rmr.test, method = "loess", se = F, 
                  mapping = aes(x = h_in_respo, y = resid),
                  color  = cols.rmr[3], fill  = cols.rmr[3])+
      geom_point(size = 1, alpha = 0.3, color  = "grey")+
      geom_point(data = data.rmr.test, mapping = aes(x = h_in_respo, y = resid),
                 size = 1, alpha = 1, color  = cols.rmr[3])+
      geom_hline(yintercept = 0)
    ggformat(RMR_method_in_respo, x_title="Hours in respo chamber",
             y_title=bquote(residual~RMR),
             size_text = 12, print = F, title = "All fish together")
    
    # MMR 
    MMR_method_in_respo<-ggplot(data = data.amrER, aes(x = h_in_respo, y = resid))+
      geom_smooth(color = "black", method = "loess", se = F)+
      geom_smooth(data = data.amr.test, method = "loess", se = F, 
                  mapping = aes(x = h_in_respo, y = resid),
                  color  = cols.amr[3], fill  = cols.amr[3])+
      geom_point(size = 1, alpha = 0.3, color  = "grey")+
      geom_point(data = data.amr.test, mapping = aes(x = h_in_respo, y = resid),
                 size = 1, alpha = 1, color  = cols.amr[3])+
      geom_hline(yintercept = 0)
    ggformat(MMR_method_in_respo, x_title="Hours in respo chamber",
             y_title=bquote(residual~MMR),
             size_text = 12, print = F, title = "All fish together")
    
    # resp to size ratio
    # RMR 
    RMR_method_resp_vol<-ggplot(data = data.rmrER, aes(x = (chamber_vol_L*1000)/BW_g, y = resid))+
      geom_smooth(color = "black", method = "loess", se = F)+
      geom_smooth(data = data.rmr.test, method = "loess", se = F, 
                  mapping = aes(x = (chamber_vol_L*1000)/BW_g, y = resid),
                  color  = cols.rmr[3], fill  = cols.rmr[3])+
      geom_point(size = 1, alpha = 0.3, color  = "grey")+
      geom_point(data = data.rmr.test, mapping = aes(x = (chamber_vol_L*1000)/BW_g, y = resid),
                 size = 1, alpha = 1, color  = cols.rmr[3])+
      geom_hline(yintercept = 0)
    ggformat(RMR_method_resp_vol, x_title="Respo chamber to body mass ratio",
             y_title=bquote(residual~RMR),
             size_text = 12, print = F, title = "All fish together")
    
    # MMR 
    MMR_method_resp_vol<-ggplot(data = data.amrER, aes(x = (chamber_vol_L*1000)/BW_g, y = resid))+
      geom_smooth(color = "black", method = "loess", se = F)+
      geom_smooth(data = data.amr.test, method = "loess", se = F, 
                  mapping = aes(x = (chamber_vol_L*1000)/BW_g, y = resid),
                  color  = cols.amr[3], fill  = cols.amr[3])+
      geom_point(size = 1, alpha = 0.3, color  = "grey")+
      geom_point(data = data.amr.test, mapping = aes(x = (chamber_vol_L*1000)/BW_g, y = resid),
                 size = 1, alpha = 1, color  = cols.amr[3])+
      geom_hline(yintercept = 0)
    ggformat(MMR_method_resp_vol, x_title="Respo chamber to body mass ratio",
             y_title=bquote(residual~MMR),
             size_text = 12, print = F, title = "All fish together")
    
    
    # save methods plots -----
    scaling_t<-cowplot:::plot_grid(MMR_method_in_respo, RMR_method_in_respo,
                                 MMR_method_starved, RMR_method_starved,
                                 MMR_method_resp_vol, RMR_method_resp_vol,
                                 MMR_method_er,MMR_method_w,
                                  align = "hv",
                                  axis = "l",
                                  nrow = 4,
                                  ncol = 2,
                                  labels = "AUTO",
                                 label_x = c(0.2, 0.2),
                                 label_y = c(0.895, 0.895),
                                  label_size = 12)
    # scaling
    ggsave(filename = paste("./Figures/FigureSUPPL_Review-methods.png", sep=""),
           plot=scaling_t, width = 10, height = 13, units = "in")
  }
  
  return(list(
   rmr_mod_ER,
   amr_mod_ER,
   as_mod_ER,
   fas_mod_ER,
   
   rmr_mod_W,
   amr_mod_W,
   as_mod_W,
   fas_mod_W
  ))
}

# Function to update global models with ecology groupings
# standardized metric for contrasts
ecol_model_update<-function(ecol.model.null, 
                            ecol.data.subset, 
                            temp.test.category, 
                            mr.type,
                            ecol.predictor, 
                            data.BIC = NULL,
                            data.ANOVA = NULL,
                            data.EMMEANS = NULL, 
                            data.CONTRASTS = NULL, 
                            data.PARAMS = NULL,
                            ref.tempTest, # for estimate means reference grid
                            ref.lnBWg # for estimate means reference grid
                            ){

  
    for (i in 1:length(ecol.predictor)){
        
      ecol.model<-Almer(update.formula(formula(ecol.model.null), paste(' ~ . + ', ecol.predictor[i], collapse = "")), data=ecol.data.subset, REML = F)

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




