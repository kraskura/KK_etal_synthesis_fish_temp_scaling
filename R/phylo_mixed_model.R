# last run: april 27 2025


# **************************************************
# run full model comparison analysis ( run = TRUE), or figures, best models, only (run = FALSE)?
# this is set for full execution of code (sourcing code); see message below
# if run is FALSE the last run's 90% CI are used; not estimated new (take a few)
# run = FALSE
run = TRUE

if(run){
  message("running all models and BIC comparison: this will take a while")
}else{
  message("NOT re-running all models and BIC comparisons, only known best-fit models: quicker")
}

# **************************************************

## Source the data and functions------------
source(here("R", "get_data_temp.R")) # dataset warngling (gets the right species names from fishbase, calculates AS, FAS, and seperates everything in temperature categories)
source(here("R", "mixed_model_outputs.R")) # function for extracting 90% CI and making prediction dataframes 
source(here("R", "get_data_phylo_matrix.R")) # phylogenetics work, extracts phylo matrix, saves it, makes plots and saves them, 
source(here("R", "colors_themes.R")) # common color schemes

# function to order model selection based on the lowest BIC score
BICdelta<-function(BICtable){
  BIC.t <- BICtable [order(BICtable$BIC), ]
  BIC.t$delta <- round(abs(BIC.t$BIC[1] -  BIC.t$BIC), 5)
  return(BIC.t)
}

# Data for models ------
data.list<-get_data_temp(data.amr = here("Data", "Fish_AMR_temp_dataset_mar2022.csv"),
                         data.rmr = here("Data","Fish_RMR_temp_dataset_mar2022.csv"),
                         ecology.data = here("Data", "Kraskura_species_ecologies_mar2022.csv"),
                         onlyTop.above = TRUE, save.FishBase.species.data = F,
                         calc_mass_specific = FALSE)

# data.list.NO.SALMON<-get_data_temp(data.amr = here("Data", "Fish_AMR_temp_dataset_mar2022.csv"),
#                                  data.rmr = here("Data","Fish_RMR_temp_dataset_mar2022.csv"),
#                                  ecology.data = here("Data", "Kraskura_species_ecologies_mar2022.csv"),
#                                  onlyTop.above = TRUE, save.FishBase.species.data = F,
#                                  calc_mass_specific = FALSE)

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

# the warm ones
data.amr.test<-rbind(data.amrAC, data.amrAM)
data.rmr.test<-rbind(data.rmrAC, data.rmrAM)
data.fas.test<-rbind(data.fasAC, data.fasAM)
data.as.test<-rbind(data.asAC, data.asAM)
data.fas.test<-data.fas.test[c(!is.na(data.fas.test$FAS) & is.finite(data.fas.test$FAS)) , ]
data.as.test<-data.as.test[c(!is.na(data.as.test$lnAS) & is.finite(data.as.test$lnAS)) , ]

# test categories
data.fas$test_category3 <- "warm"
data.fas[data.fas$test_category == "ecol_relev", "test_category3"] <- "optimal"
data.as$test_category3 <- "warm"
data.as[data.as$test_category == "ecol_relev", "test_category3"] <- "optimal"
data.amr$test_category3 <- "warm"
data.amr[data.amr$test_category == "ecol_relev", "test_category3"] <- "optimal"
data.rmr$test_category3 <- "warm"
data.rmr[data.rmr$test_category == "ecol_relev", "test_category3"] <- "optimal"
 
# Phylo models -----------------
# get phylo trees:
data.rmrER<-droplevels(data.rmrER)
data.rmr.test<-droplevels(data.rmr.test)
data.amrER<-droplevels(data.amrER)
data.amr.test<-droplevels(data.amr.test)
data.asER<-droplevels(data.asER)
data.as.test<-droplevels(data.as.test)
data.fasER<-droplevels(data.fasER)
data.fas.test<-droplevels(data.fas.test)

get_phylo_matrix(species.list = unique(levels(data.rmrER$species)), matrix.name = "A", tree.name = "tree.kk", dataset.ID = "RMR optimal")
get_phylo_matrix(species.list = unique(levels(data.rmr.test$species)), matrix.name = "A.rmr.w", tree.name = "tr.rmr.w", dataset.ID = "RMR warm")

get_phylo_matrix(species.list = unique(levels(data.amrER$species)), matrix.name = "A.mmr.er", tree.name = "tr.mmr.er", dataset.ID = "MMR optimal")
get_phylo_matrix(species.list = unique(levels(data.amr.test$species)), matrix.name = "A.mmr.w", tree.name = "tr.mmr.w", dataset.ID = "MMR warm")

get_phylo_matrix(species.list = unique(levels(data.asER$species)), matrix.name = "A.aas.er", tree.name = "tr.aas.er", dataset.ID = "AAS optimal")
get_phylo_matrix(species.list = unique(levels(data.as.test$species)), matrix.name = "A.aas.w", tree.name = "tr.aas.w", dataset.ID = "AAS warm")

get_phylo_matrix(species.list = unique(levels(data.fasER$species)), matrix.name = "A.fas.er", tree.name = "tr.fas.er", dataset.ID = "FAS optimal")
get_phylo_matrix(species.list = unique(levels(data.fas.test$species)), matrix.name = "A.fas.w", tree.name = "tree.fas.w", dataset.ID = "FAS warm")

# runnign take like 5 or so minutes. 
if (run){
  # RMR optimal -----------------
  Phylo_RMR_model0 <- Almer(lnRMR ~ lnBWg + tempTest + (1|species) , data=data.rmrER, REML=FALSE,A = list(species = A))
  Phylo_RMR_model0.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model0intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model0int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) , data=data.rmrER, REML=FALSE, A = list(species = A))
  
  Phylo_RMR_model1 <- Almer(lnRMR ~ lnBWg + tempTest + (1|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model1.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model1intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model1int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  
  Phylo_RMR_model2 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model2int <- Almer(lnRMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model2.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model2intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE)  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  
  Phylo_RMR_model4 <- Almer(lnRMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model4int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model4.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model4intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  
  Phylo_RMR_model5 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model5int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model5.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model5intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  
  Phylo_RMR_model6 <- Almer(lnRMR ~ lnBWg + tempTest + (lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model6int <- Almer(lnRMR ~ lnBWg * tempTest+ (lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model6.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model6intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  
  

  
  # No singular fits 
  ### BIC rmr optimal --------------
  RMR_BIC<-BICdelta(BIC(Phylo_RMR_model0, Phylo_RMR_model0int, Phylo_RMR_model0.POLY, Phylo_RMR_model0intPOLY,
                  Phylo_RMR_model1, Phylo_RMR_model1int, Phylo_RMR_model1.POLY, Phylo_RMR_model1intPOLY,
                  Phylo_RMR_model2, Phylo_RMR_model2int, Phylo_RMR_model2.POLY, Phylo_RMR_model2intPOLY,
                  Phylo_RMR_model4, Phylo_RMR_model4int, Phylo_RMR_model4.POLY, Phylo_RMR_model4intPOLY,
                  Phylo_RMR_model5, Phylo_RMR_model5int, Phylo_RMR_model5.POLY, Phylo_RMR_model5intPOLY,
                  Phylo_RMR_model6, Phylo_RMR_model6int, Phylo_RMR_model6.POLY, Phylo_RMR_model6intPOLY))
  
  ## MMR optimal -----------------
  Phylo_MMR_model0 <- Almer(lnAMR ~ lnBWg + tempTest + (1|species) , data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model0.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model0intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model0int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) , data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  
  Phylo_MMR_model1 <- Almer(lnAMR ~ lnBWg + tempTest + (1|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model1.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model1intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model1int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  
  Phylo_MMR_model2 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model2int <- Almer(lnAMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model2.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model2intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE)  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  
  Phylo_MMR_model4 <- Almer(lnAMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model4int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model4.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model4intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  
  Phylo_MMR_model5 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model5int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model5.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model5intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  
  Phylo_MMR_model6 <- Almer(lnAMR ~ lnBWg + tempTest + (lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model6int <- Almer(lnAMR ~ lnBWg * tempTest+ (lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model6.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model6intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  
  ### BIC --------------
  MMR_BIC<-BICdelta(BIC(Phylo_MMR_model0, Phylo_MMR_model0int, Phylo_MMR_model0.POLY, Phylo_MMR_model0intPOLY,
                  Phylo_MMR_model1, Phylo_MMR_model1int, Phylo_MMR_model1.POLY, Phylo_MMR_model1intPOLY,
                  Phylo_MMR_model2, Phylo_MMR_model2int, Phylo_MMR_model2.POLY, Phylo_MMR_model2intPOLY,
                  Phylo_MMR_model4, Phylo_MMR_model4int, Phylo_MMR_model4.POLY, Phylo_MMR_model4intPOLY,
                  Phylo_MMR_model5, Phylo_MMR_model5int, Phylo_MMR_model5.POLY, Phylo_MMR_model5intPOLY,
                  Phylo_MMR_model6, Phylo_MMR_model6int, Phylo_MMR_model6.POLY, Phylo_MMR_model6intPOLY))
  
  ## AAS optimal ------------------
  Phylo_AS_model0 <- Almer(lnAS ~ lnBWg + tempTest + (1|species) , data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model0.POLY <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model0intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model0int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) , data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  
  Phylo_AS_model1 <- Almer(lnAS ~ lnBWg + tempTest + (1|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model1.POLY <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model1intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model1int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  
  Phylo_AS_model2 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model2int <- Almer(lnAS ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model2.POLY <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model2intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE)  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  
  Phylo_AS_model4 <- Almer(lnAS ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model4int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model4.POLY <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model4intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  
  Phylo_AS_model5 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model5int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model5.POLY <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model5intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  
  Phylo_AS_model6 <- Almer(lnAS ~ lnBWg + tempTest + (lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model6int <- Almer(lnAS ~ lnBWg * tempTest+ (lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model6.POLY <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model6intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  
  ### BIC --------
  AS_BIC<-BICdelta(BIC(Phylo_AS_model0, Phylo_AS_model0int, Phylo_AS_model0.POLY, Phylo_AS_model0intPOLY,
                  Phylo_AS_model1, Phylo_AS_model1int, Phylo_AS_model1.POLY, Phylo_AS_model1intPOLY,
                  Phylo_AS_model2, Phylo_AS_model2int, Phylo_AS_model2.POLY, Phylo_AS_model2intPOLY,
                  Phylo_AS_model4, Phylo_AS_model4int, Phylo_AS_model4.POLY, Phylo_AS_model4intPOLY,
                  Phylo_AS_model5, Phylo_AS_model5int, Phylo_AS_model5.POLY, Phylo_AS_model5intPOLY,
                  Phylo_AS_model6, Phylo_AS_model6int, Phylo_AS_model6.POLY, Phylo_AS_model6intPOLY))
  
  
  ## FAS optimal ----------------
  Phylo_FAS_model0 <- Almer(log(FAS) ~ lnBWg + tempTest + (1|species) , data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model0.POLY <- Almer(log(FAS) ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model0intPOLY <- Almer(log(FAS) ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model0int <- Almer(log(FAS) ~ lnBWg * tempTest + (1|species) , data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  
  Phylo_FAS_model1 <- Almer(log(FAS) ~ lnBWg + tempTest + (1|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model1.POLY <- Almer(log(FAS) ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model1intPOLY <- Almer(log(FAS) ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model1int <- Almer(log(FAS) ~ lnBWg * tempTest + (1|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  
  Phylo_FAS_model2 <- Almer(log(FAS) ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model2int <- Almer(log(FAS) ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model2.POLY <- Almer(log(FAS) ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model2intPOLY <- Almer(log(FAS) ~ lnBWg * poly(tempTest,2, raw = TRUE)  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  
  Phylo_FAS_model4 <- Almer(log(FAS) ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model4int <- Almer(log(FAS) ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model4.POLY <- Almer(log(FAS) ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model4intPOLY <- Almer(log(FAS) ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  
  Phylo_FAS_model5 <- Almer(log(FAS) ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model5int <- Almer(log(FAS) ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model5.POLY <- Almer(log(FAS) ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model5intPOLY <- Almer(log(FAS) ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  
  Phylo_FAS_model6 <- Almer(lnFAS ~ lnBWg + tempTest + (lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model6int <- Almer(lnFAS ~ lnBWg * tempTest+ (lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model6.POLY <- Almer(lnFAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model6intPOLY <- Almer(lnFAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  
  ### BIC -------------
  FAS_BIC<-BICdelta(BIC(Phylo_FAS_model0, Phylo_FAS_model0int, Phylo_FAS_model0.POLY, Phylo_FAS_model0intPOLY,
                  Phylo_FAS_model1, Phylo_FAS_model1int, Phylo_FAS_model1.POLY, Phylo_FAS_model1intPOLY,
                  Phylo_FAS_model2, Phylo_FAS_model2int, Phylo_FAS_model2.POLY, Phylo_FAS_model2intPOLY,
                  Phylo_FAS_model4, Phylo_FAS_model4int, Phylo_FAS_model4.POLY, Phylo_FAS_model4intPOLY,
                  Phylo_FAS_model5, Phylo_FAS_model5int, Phylo_FAS_model5.POLY, Phylo_FAS_model5intPOLY,
                  Phylo_FAS_model6, Phylo_FAS_model6int, Phylo_FAS_model6.POLY, Phylo_FAS_model6intPOLY))
  
  
  # RMR warm ------------------
  Phylo_RMR_W_model0 <- Almer(lnRMR ~ lnBWg + tempTest + (1|species) , data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model0.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model0intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model0int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) , data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  
  Phylo_RMR_W_model1 <- Almer(lnRMR ~ lnBWg + tempTest + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model1.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model1intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model1int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  
  Phylo_RMR_W_model2 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model2int <- Almer(lnRMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model2.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model2intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE)  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  
  Phylo_RMR_W_model4 <- Almer(lnRMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model4int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model4.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model4intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  
  Phylo_RMR_W_model5 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model5int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model5.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model5intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  
  Phylo_RMR_W_model6 <- Almer(lnRMR ~ lnBWg + tempTest + (lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model6int <- Almer(lnRMR ~ lnBWg * tempTest+ (lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model6.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model6intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  
  
  ## BIC -------------
  RMR_W_BIC<-BICdelta(BIC(Phylo_RMR_W_model0, Phylo_RMR_W_model0int, Phylo_RMR_W_model0.POLY, Phylo_RMR_W_model0intPOLY,
                  Phylo_RMR_W_model1, Phylo_RMR_W_model1int, Phylo_RMR_W_model1.POLY, Phylo_RMR_W_model1intPOLY,
                  Phylo_RMR_W_model2, Phylo_RMR_W_model2int, Phylo_RMR_W_model2.POLY, Phylo_RMR_W_model2intPOLY,
                  Phylo_RMR_W_model4, Phylo_RMR_W_model4int, Phylo_RMR_W_model4.POLY, Phylo_RMR_W_model4intPOLY,
                  Phylo_RMR_W_model5, Phylo_RMR_W_model5int, Phylo_RMR_W_model5.POLY, Phylo_RMR_W_model5intPOLY,
                  Phylo_RMR_W_model6, Phylo_RMR_W_model6int, Phylo_RMR_W_model6.POLY, Phylo_RMR_W_model6intPOLY))
  
  # AMR / warm temps --------------
  Phylo_MMR_W_model0 <- Almer(lnAMR ~ lnBWg + tempTest + (1|species) , data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model0.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model0intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model0int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) , data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  
  Phylo_MMR_W_model1 <- Almer(lnAMR ~ lnBWg + tempTest + (1|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model1.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model1intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model1int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  
  Phylo_MMR_W_model2 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model2int <- Almer(lnAMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model2.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model2intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE)  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  
  Phylo_MMR_W_model4 <- Almer(lnAMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model4int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model4.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model4intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  
  Phylo_MMR_W_model5 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model5int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model5.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model5intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  
  Phylo_MMR_W_model6 <- Almer(lnAMR ~ lnBWg + tempTest + (lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model6int <- Almer(lnAMR ~ lnBWg * tempTest+ (lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model6.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model6intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  
  ## BIC -------------
  MMR_W_BIC<-BICdelta(BIC(Phylo_MMR_W_model0, Phylo_MMR_W_model0int, Phylo_MMR_W_model0.POLY, Phylo_MMR_W_model0intPOLY,
                  Phylo_MMR_W_model1, Phylo_MMR_W_model1int, Phylo_MMR_W_model1.POLY, Phylo_MMR_W_model1intPOLY,
                  Phylo_MMR_W_model2, Phylo_MMR_W_model2int, Phylo_MMR_W_model2.POLY, Phylo_MMR_W_model2intPOLY,
                  Phylo_MMR_W_model4, Phylo_MMR_W_model4int, Phylo_MMR_W_model4.POLY, Phylo_MMR_W_model4intPOLY,
                  Phylo_MMR_W_model5, Phylo_MMR_W_model5int, Phylo_MMR_W_model5.POLY, Phylo_MMR_W_model5intPOLY,
                  Phylo_MMR_W_model6, Phylo_MMR_W_model6int, Phylo_MMR_W_model6.POLY, Phylo_MMR_W_model6intPOLY))
  
  
  
  # AS / warm temps --------------
  Phylo_AS_W_model0 <- Almer(lnAS ~ lnBWg + tempTest + (1|species) , data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model0.POLY <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model0intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model0int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) , data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  
  Phylo_AS_W_model1 <- Almer(lnAS ~ lnBWg + tempTest + (1|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model1.POLY <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model1intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model1int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  
  Phylo_AS_W_model2 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model2int <- Almer(lnAS ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model2.POLY <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model2intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE)  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  
  Phylo_AS_W_model4 <- Almer(lnAS ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model4int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model4.POLY <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model4intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  
  Phylo_AS_W_model5 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model5int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model5.POLY <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model5intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  
  Phylo_AS_W_model6 <- Almer(lnAS ~ lnBWg + tempTest + (lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model6int <- Almer(lnAS ~ lnBWg * tempTest+ (lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model6.POLY <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model6intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  
  ## BIC -------------
  AS_W_BIC<-BICdelta(BIC(Phylo_AS_W_model0, Phylo_AS_W_model0int, Phylo_AS_W_model0.POLY, Phylo_AS_W_model0intPOLY,
                  Phylo_AS_W_model1, Phylo_AS_W_model1int, Phylo_AS_W_model1.POLY, Phylo_AS_W_model1intPOLY,
                  Phylo_AS_W_model2, Phylo_AS_W_model2int, Phylo_AS_W_model2.POLY, Phylo_AS_W_model2intPOLY,
                  Phylo_AS_W_model4, Phylo_AS_W_model4int, Phylo_AS_W_model4.POLY, Phylo_AS_W_model4intPOLY,
                  Phylo_AS_W_model5, Phylo_AS_W_model5int, Phylo_AS_W_model5.POLY, Phylo_AS_W_model5intPOLY,
                  Phylo_AS_W_model6, Phylo_AS_W_model6int, Phylo_AS_W_model6.POLY, Phylo_AS_W_model6intPOLY))
  
  
  # FAS / warm temps --------------
  Phylo_FAS_W_model0 <- Almer(lnFAS ~ lnBWg + tempTest + (1|species) , data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model0.POLY <- Almer(lnFAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model0intPOLY <- Almer(lnFAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model0int <- Almer(lnFAS ~ lnBWg * tempTest + (1|species) , data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  
  Phylo_FAS_W_model1 <- Almer(lnFAS ~ lnBWg + tempTest + (1|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model1.POLY <- Almer(lnFAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model1intPOLY <- Almer(lnFAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model1int <- Almer(lnFAS ~ lnBWg * tempTest + (1|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  
  Phylo_FAS_W_model2 <- Almer(lnFAS ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model2int <- Almer(lnFAS ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model2.POLY <- Almer(lnFAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model2intPOLY <- Almer(lnFAS ~ lnBWg * poly(tempTest,2, raw = TRUE)  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  
  Phylo_FAS_W_model4 <- Almer(lnFAS ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model4int <- Almer(lnFAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model4.POLY <- Almer(lnFAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model4intPOLY <- Almer(lnFAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  
  Phylo_FAS_W_model5 <- Almer(lnFAS ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model5int <- Almer(lnFAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model5.POLY <- Almer(lnFAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model5intPOLY <- Almer(lnFAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  
  Phylo_FAS_W_model6 <- Almer(lnFAS ~ lnBWg + tempTest + (lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model6int <- Almer(lnFAS ~ lnBWg * tempTest+ (lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model6.POLY <- Almer(lnFAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model6intPOLY <- Almer(lnFAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  
  ## BIC -------------
  FAS_W_BIC<-BICdelta(BIC(Phylo_FAS_W_model0, Phylo_FAS_W_model0int, Phylo_FAS_W_model0.POLY, Phylo_FAS_W_model0intPOLY,
                  Phylo_FAS_W_model1, Phylo_FAS_W_model1int, Phylo_FAS_W_model1.POLY, Phylo_FAS_W_model1intPOLY,
                  Phylo_FAS_W_model2, Phylo_FAS_W_model2int, Phylo_FAS_W_model2.POLY, Phylo_FAS_W_model2intPOLY,
                  Phylo_FAS_W_model4, Phylo_FAS_W_model4int, Phylo_FAS_W_model4.POLY, Phylo_FAS_W_model4intPOLY,
                  Phylo_FAS_W_model5, Phylo_FAS_W_model5int, Phylo_FAS_W_model5.POLY, Phylo_FAS_W_model5intPOLY,
                  Phylo_FAS_W_model6, Phylo_FAS_W_model6int, Phylo_FAS_W_model6.POLY, Phylo_FAS_W_model6intPOLY))
  
  ## BIC results ------
  RMR_BIC # Phylo_RMR_model5
  MMR_BIC # Phylo_MMR_model4int
  AS_BIC  # Phylo_AS_model4
  FAS_BIC # Phylo_FAS_model4
  
  RMR_W_BIC # Phylo_RMR_W_model1
  MMR_W_BIC # Phylo_MMR_W_model4
  AS_W_BIC  # PPhylo_AS_W_model2.POLY
  FAS_W_BIC # Phylo_FAS_W_model2.POLY
  
  write.csv(file = here("Data_exports/RMR_BIC.csv"), x = RMR_BIC, row.names = T)
  write.csv(file = here("Data_exports/MMR_BIC.csv"), x = MMR_BIC, row.names = T)
  write.csv(file = here("Data_exports/AS_BIC.csv"), x = AS_BIC, row.names = T)
  write.csv(file = here("Data_exports/FAS_BIC.csv"), x = FAS_BIC, row.names = T)

  write.csv(file = here("Data_exports/RMR_W_BIC.csv"), x = RMR_W_BIC, row.names = T)
  write.csv(file = here("Data_exports/MMR_W_BIC.csv"), x = MMR_W_BIC, row.names = T)
  write.csv(file = here("Data_exports/AS_W_BIC.csv"), x = AS_W_BIC, row.names = T)
  write.csv(file = here("Data_exports/FAS_W_BIC.csv"), x = FAS_W_BIC, row.names = T)

  # Best models  --------
  rmr_mod_ER<-Phylo_RMR_model5 # without the interaction
  # rmr_mod_ER<-Phylo_RMR_model4int # with the interaction when zebrafish study excluded # 415
  amr_mod_ER<-Phylo_MMR_model4int
  as_mod_ER<-Phylo_AS_model4
  fas_mod_ER<-Phylo_FAS_model4

  rmr_mod_W<-Phylo_RMR_W_model1
  amr_mod_W<-Phylo_MMR_W_model4
  as_mod_W<-Phylo_AS_W_model2.POLY
  fas_mod_W<-Phylo_FAS_W_model2.POLY
  
}else{
  
  # Phylo_RMR_model5 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_MMR_model4int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # running best models only 
    # Phylo_AS_model4 <- Almer(lnAS ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
    # Phylo_FAS_model4 <- Almer(log(FAS) ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
    
  # Phylo_RMR_W_model1 <- Almer(lnRMR ~ lnBWg + tempTest + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_MMR_W_model4 <- Almer(lnAMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))   
  # Phylo_AS_W_model2.POLY <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_FAS_W_model2.POLY <- Almer(lnFAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w)) 
  
  # with interaction
  # rmr_mod_ER <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))

  # without the interaction.
  rmr_mod_ER <-  Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  amr_mod_ER <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  as_mod_ER <- Almer(lnAS ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  fas_mod_ER <- Almer(log(FAS) ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  
  rmr_mod_W <- Almer(lnRMR ~ lnBWg + tempTest + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  amr_mod_W <- Almer(lnAMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))   
  as_mod_W <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  fas_mod_W <- Almer(lnFAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w)) 
  
}


# un-comment to see best model summaries or view output from R markdown
# # Ecol relev
# summary(rmr_mod_ER)
# summary(amr_mod_ER)
# summary(fas_mod_ER)
# summary(as_mod_ER)
# 
# # warm
# summary(rmr_mod_W)
# summary(amr_mod_W)
# summary(fas_mod_W)
# summary(as_mod_W)
# 
# # Ecol relev
# plot(rmr_mod_ER)
# plot(amr_mod_ER)
# plot(fas_mod_ER)
# plot(as_mod_ER)
# 
# # warm
# plot(rmr_mod_W)
# plot(amr_mod_W)
# plot(fas_mod_W)
# plot(as_mod_W)
# 
# hist(resid(amr_mod_ER))
# hist(resid(rmr_mod_ER))
# hist(resid(fas_mod_ER))
# hist(resid(as_mod_ER))
# 
# hist(resid(amr_mod_W))
# hist(resid(rmr_mod_W))
# hist(resid(fas_mod_W))
# hist(resid(as_mod_W)) # little skew not too bad, only 5 measurements, all reasonable biologically

# data.as.test$resid<-resid(as_mod_W)
# plot(resid(as_mod_W))
# # data.as.test[which(data.as.test$resid < -1.5),] # could be considered outlier data, n = 5
# data.fasER[which(data.fasER$FAS > 20),] # considered outliers, measurement extremes for zebrafish; n =14


# *****************************************************************************
# *****************************************************************************

# Model scaling parameters and CIs ----------
# custom function to obtain model parameters and recalculated mass specific MR using model estimate scaling slopes
# this also expands the dataset to get mass independent values of metabolic rates
# function source script available at 'mixed_model_outputs.R'
model_out<-model_outputs(phylo = TRUE, 
              best.model.rmr.er = rmr_mod_ER,
              best.model.amr.er= amr_mod_ER,
              best.model.as.er = as_mod_ER,
              best.model.fas.er = fas_mod_ER,
              best.model.rmr.w = rmr_mod_W,
              best.model.amr.w = amr_mod_W,
              best.model.as.w = as_mod_W,
              best.model.fas.w = fas_mod_W,
              estimate.CI = run) # takes a few hours; outputs saved 

read.csv("./Data_exports/Phylo/Table_CIsummary.csv")
colnames(sum_CItable)<-c("var", "ci5", "ci95","var_repeat", "MR", "temp_cat") # if from read.csv

# use this 
sum_CItable<-data.frame(model_out[[1]])
colnames(sum_CItable)<-c("ci5", "ci95","var_repeat", "MR", "temp_cat")

scaling.params<-data.frame(model_out[[2]])
sum_data<-data.frame(model_out[[3]])


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
AMR_slope31<-round(as.numeric(scaling.params[scaling.params$performance == "MMR" & 
                            scaling.params$temp_categ == "er" & 
                            scaling.params$tempTest == "31"  , "lnBWg"]),2)
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

data.plotRMR_ER<-read.csv(file = here("Data_exports/Phylo/dataPred_RMR_er.csv"))
data.plotAMRint_ER<-read.csv(file = here("./Data_exports/Phylo/dataPred_AMR_er.csv"))
data.plotAS_ER<-read.csv(file = here("./Data_exports/Phylo/dataPred_AS_er.csv"))
data.plotFAS_ER<-read.csv(file = here("./Data_exports/Phylo/dataPred_FAS_er.csv"))

data.plotRMR_warm<-read.csv(here("./Data_exports/Phylo/dataPred_RMR_warm.csv"))
data.plotAMR_warm<-read.csv(file = here("./Data_exports/Phylo/dataPred_AMR_warm.csv"))
data.plotAS_warm<-read.csv(file = here("./Data_exports/Phylo/dataPred_AS_warm.csv"))
data.plotFAS_warm<-read.csv(file = here("./Data_exports/Phylo/dataPred_FAS_warm.csv"))

# *****************************************************************************
# *****************************************************************************




# ******************************************************************************
# Breakpoint? (NOT USED) ------------
# library(segmented)
# 
# 
# seg.amrER <-lm(log(mass_specamr) ~ tempTestK1000, data=data.amrER)
# seg.amr.test <-lm(log(mass_specamr) ~ tempTestK1000, data=data.amr.test)
# 
# seg.rmrER <-lm(log(mass_specrmr) ~ tempTestK1000, data=data.rmrER)
# seg.rmr.test <-lm(log(mass_specrmr) ~ tempTestK1000, data=data.rmr.test)
# 
# seg.asER <-lm(log(mass_specas) ~ tempTestK1000, data=data.asER)
# seg.as.test <-lm(log(mass_specas) ~ tempTestK1000, data=data.as.test)
# 
# 
# seg.amrER.r <- segmented(seg.amrER , ~tempTestK1000)
# seg.amr.test.r <- segmented(seg.amr.test , ~tempTestK1000)
# 
# plot(y = log(data.amrER$mass_specamr), x = data.amrER$tempTestK1000)
# plot(seg.amrER.r, add = TRUE)
# plot(y = log(data.amr.test$mass_specamr), x = data.amr.test$tempTestK1000)
# plot(seg.amr.test.r, add = TRUE)
# 


# ******************************************************************************
# ******************************************************************************
# models within SIC 7  -------
# 
# models_SIC<-list(Phylo_RMR_model5, 
#               Phylo_RMR_model4,
#               Phylo_RMR_model4int,
#               Phylo_MMR_model4int,
#               Phylo_MMR_model4,
#               Phylo_AS_model4,
#               Phylo_AS_model4int,
#               Phylo_FAS_model4,
#               Phylo_FAS_model4int,
#               Phylo_FAS_model5,
#               Phylo_RMR_W_model1,
#               Phylo_RMR_W_model2,
#               Phylo_RMR_W_model1int,
#               Phylo_RMR_W_model1.POLY,
#               Phylo_RMR_W_model4,
#               Phylo_RMR_W_model2.POLY,
#               Phylo_RMR_W_model5,
#               Phylo_MMR_W_model4,
#               Phylo_MMR_W_model5,
#               Phylo_MMR_W_model5int,
#               Phylo_MMR_W_model4int,
#               Phylo_MMR_W_model1int,
#               Phylo_MMR_W_model1,
#               Phylo_MMR_W_model2int,
#               Phylo_MMR_W_model4.POLY,
#               Phylo_MMR_W_model5.POLY,
#               Phylo_AS_W_model2.POLY,
#               Phylo_AS_W_model2int,
#               Phylo_FAS_W_model2.POLY,
#               Phylo_FAS_W_model1,
#               Phylo_FAS_W_model1intPOLY,
#               Phylo_FAS_W_model2intPOLY)
# f <- function(x) {
#     d <- substitute(x)
#     n <- sapply(d[-1],deparse)
#     return(n)
#  }
# models_SIC_names<-f(c(Phylo_RMR_model5, 
#               Phylo_RMR_model4,
#               Phylo_RMR_model4int,
#               Phylo_MMR_model4int,
#               Phylo_MMR_model4,
#               Phylo_AS_model4,
#               Phylo_AS_model4int,
#               Phylo_FAS_model4,
#               Phylo_FAS_model4int,
#               Phylo_FAS_model5,
#               Phylo_RMR_W_model1,
#               Phylo_RMR_W_model2,
#               Phylo_RMR_W_model1int,
#               Phylo_RMR_W_model1.POLY,
#               Phylo_RMR_W_model4,
#               Phylo_RMR_W_model2.POLY,
#               Phylo_RMR_W_model5,
#               Phylo_MMR_W_model4,
#               Phylo_MMR_W_model5,
#               Phylo_MMR_W_model5int,
#               Phylo_MMR_W_model4int,
#               Phylo_MMR_W_model1int,
#               Phylo_MMR_W_model1,
#               Phylo_MMR_W_model2int,
#               Phylo_MMR_W_model4.POLY,
#               Phylo_MMR_W_model5.POLY,
#               Phylo_AS_W_model2.POLY,
#               Phylo_AS_W_model2int,
#               Phylo_FAS_W_model2.POLY,
#               Phylo_FAS_W_model1,
#               Phylo_FAS_W_model1intPOLY,
#               Phylo_FAS_W_model2intPOLY))
# 
# new_df <- data.frame(matrix(ncol = 6, nrow = 0))
# # loop to get all neccessary data for model comparisons 
# for(i in 1:length(models_SIC)){
#   print(i)
#   Almer_model<-models_SIC[[i]]
#   
#   model_name<-capture.output(summary(Almer_model))[2]
#   BIC_val<-BIC(Almer_model)
#   slope_val<-fixef(Almer_model)["lnBWg"]
#   interc_val<-fixef(Almer_model)["(Intercept)"]
#   temp_val<-fixef(Almer_model)["tempTest"]
#   model_type<-models_SIC_names[[i]]
#   
#   colnames(new_df)<- c("metabolism metric",	"model",	"BIC",
#                        "scaling slope",	"intercept coef",	"temperature coef")
#   new_vals<-as.data.frame(t(c(model_type, model_name,
#                      BIC_val, slope_val,
#                      interc_val, temp_val)))
#   colnames(new_vals)<- c("metabolism metric",	"model",	"BIC",
#                        "scaling slope",	"intercept coef",	"temperature coef")
#   
#   new_df<-rbind(new_df, new_vals)
# 
# }
# write.csv(file = here("Data_exports/model_comparison.csv"), x = new_df, row.names = F)


# ******************************************************************************
# Figures -------


set.seed(51423)
# MMR and AMR used interchangeably throughout 

# predict vs real 
# ggplot(data.amrER, aes(predict(amr_mod_ER), lnAMR, color = tempTest))+
#   geom_point()
# ggplot(data.rmrER, aes(predict(rmr_mod_ER), lnRMR, color = tempTest))+
#   geom_point()
# ggplot(data.asER, aes(predict(as_mod_ER), lnAS, color = tempTest))+
#   geom_point()
# ggplot(data.fasER, aes(predict(fas_mod_ER), lnFAS, color = tempTest))+
#   geom_point()
# plot(resid(fas_mod_ER))

# predict vs real 
# ggplot(data.amr.test, aes(predict(amr_mod_W), lnAMR, color = tempTest))+
#   geom_point()
# ggplot(data.rmr.test, aes(predict(rmr_mod_W), lnRMR, color = tempTest))+
#   geom_point()
# ggplot(data.as.test, aes(predict(as_mod_W), lnAS, color = tempTest))+
#   geom_point()
# ggplot(data.fas.test, aes(predict(fas_mod_W), lnFAS, color = tempTest))+
#   geom_point()
# plot(resid(fas_mod_W))


# General scaling plots ------
AMRmodel_plot1<-ggplot(data=data.amrER, aes(x=lnBWg, y=lnAMR)) +
  geom_point(alpha=0.9,  size=1, pch=1, color="grey75")+
  geom_line(data=data.plotAMRint_ER[round(data.plotAMRint_ER$tempTest,2)== 5 |
                                      round(data.plotAMRint_ER$tempTest,2)== 15|
                                      round(data.plotAMRint_ER$tempTest,2)== 25|
                                      round(data.plotAMRint_ER$tempTest,2)== 31,],
            aes(y = model_predFE, x=lnBWg,  group=tempTest, color = tempTest),
            linewidth=0.5, lty=1,alpha=0.8, show.legend=FALSE) +
  geom_point(data=data.amr.test, aes(x=lnBWg, y=lnAMR),
             alpha=0.9,  size=1, pch=21, show.legend = FALSE, stroke =0.2,
             fill = cols.amr[3], color = cols.amr[1])+
  geom_line(data=data.plotAMR_warm[round(data.plotAMR_warm$tempTest,2)==25,],
            aes(y = model_predFE, x=lnBWg,  group=tempTest), color=cols.amr[1],
            linewidth=1, lty=1, show.legend=FALSE) +
  annotate("text",  x = -5.2, y = 11.5,
           label = bquote(Optimal:~italic(b)[MMR] == change~with~degree*C),
           size=4, hjust=0, family="Helvetica", color = "black")+
  annotate("text",  x = -5.2, y = 9.8,
           label = bquote(Warm:~italic(b)[MMR] == .(AMR_slope_w)),
           size=4, hjust=0, family="Helvetica", color = cols.amr[1])+
  annotate("text",  x = 2.5, y = -3,
           label = bquote(~5*degree*C:~italic(b)[MMR] == .(AMR_slope5)),
           size=3.5, hjust=0, family="Helvetica", color = "black")+
  annotate("text",  x = 2.5, y = -4,
           label = bquote(15*degree*C:~italic(b)[MMR] == .(AMR_slope15)),
           size=3.5, hjust=0, family="Helvetica", color = "black")+
  annotate("text",  x = 2.5, y = -5,
           label = bquote(25*degree*C:~italic(b)[MMR] == .(AMR_slope25)),
           size=3.5, hjust=0, family="Helvetica", color = "black")+
  annotate("text",  x = 2.5, y = -6,
           label = bquote(31*degree*C:~italic(b)[MMR] == .(AMR_slope31)),
           size=3.5, hjust=0, family="Helvetica", color = "black")+
  # scale_fill_gradient( low = cols.amr[5], high = cols.amr[1])+
  scale_color_gradient( low = "grey", high = "black")+
  ylim(x = -6.3, 12)+
  xlim(x = -6.3, 12)+
  annotate("text", label = paste("n = ", nrow(data.amrER), sep=""),
           x = -5.2, y = 8.5, size=3, hjust=0, family="Helvetica", color = "black")+
  annotate("text", label = paste("n = ", nrow(data.amr.test), sep=""),
           x = -5.2, y = 7.6, size=3, hjust=0, family="Helvetica", color = cols.amr[1])
ggformat(AMRmodel_plot1, x_title=expression(italic(ln)*Body~weight~(g)), y_title=expression(italic(ln)*MMR~(mg~O[2]~h^-1)), size_text = 12, print = F)

RMRmodel_plot1<-ggplot(data=data.rmrER, aes(x=lnBWg, y=lnRMR)) +
  geom_point(alpha=0.9, size=1, pch=1, color="grey75")+
  geom_line(data=data.plotRMR_ER[round(data.plotRMR_ER$tempTest,2)==25,],
            mapping=aes(y = model_predFE, x=lnBWg,  group=tempTest, color= tempTest),
            color="black", linewidth=0.5, lty=1, alpha=0.8, show.legend=FALSE) +
  geom_point(data=data.rmr.test, aes(x=lnBWg, y=lnRMR),
             alpha=0.9,  size=1, pch=21, show.legend = FALSE, stroke =0.2,
             fill = cols.rmr[3], color = cols.rmr[1])+
  geom_line(data=data.plotRMR_warm[round(data.plotRMR_warm$tempTest,2)==25,],
            aes(y = model_predFE, x=lnBWg,  group=tempTest),
            color=cols.rmr[1], linewidth=1, lty=1, show.legend=FALSE) +
  annotate("text",  x = -5.2, y = 11.5, label = bquote(Optimal:~italic(b)[RMR] == .(RMR_slope)),
           size=4, hjust=0, family="Helvetica", color = "black")+
  annotate("text",  x = -5.2, y = 9.8, label = bquote(Warm:~italic(b)[RMR] == .(RMR_slope_w)),
           size=4, hjust=0, family="Helvetica", color = cols.rmr[1])+
  scale_fill_gradient( low = cols.rmr[5], high = cols.rmr[1])+
  scale_color_gradient( low = cols.rmr[5], high = cols.rmr[1])+
  ylim(x = -6.3, 12)+
  xlim(x = -6.3, 12)+
  annotate("text", label = paste("n = ", nrow(data.rmrER), sep=""),
           x = -5.2, y = 8.5, size=3, hjust=0, family="Helvetica", color = "black")+
  annotate("text", label = paste("n = ", nrow(data.rmr.test), sep=""),
           x = -5.2, y = 7.6, size=3, hjust=0, family="Helvetica", color = cols.rmr[1])
ggformat(RMRmodel_plot1, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*RMR~(mg~O[2]~h^-1)),size_text = 12, print = F)

# FAS! 
FASmodel_plot1<-ggplot(data=data.fasER, aes(x=lnBWg, y=log(FAS))) +
  geom_point(alpha=0.9, size=1, pch=1, color="grey75")+
  geom_line(data=data.plotFAS_ER[round(data.plotFAS_ER$tempTest,2)==25,],
            aes(y = model_predFE, x=lnBWg,  group=tempTest),
            color="black", linewidth=0.5, lty=1, alpha=0.8, show.legend=FALSE) +
  geom_point(data=data.fas.test, aes(x=lnBWg, y=log(FAS), fill=tempTest, color = tempTest),
             alpha=0.9, size=1, pch=21, show.legend = FALSE, stroke =0.2,
             fill = cols.fas[3], color = cols.fas[1])+
  geom_line(data=data.plotFAS_warm[round(data.plotFAS_warm$tempTest,2)==25,],
            aes(y = model_predFE, x=lnBWg,  group=tempTest), color=cols.fas[1], linewidth=1, lty=1, show.legend=FALSE) +
  annotate("text",  x = 2, y = 3.9, label = bquote(Optimal:~italic(b)[FAS] == .(FAS_slope)),
           size=4, hjust=0, family="Helvetica", color = "black")+
  annotate("text",  x = 2, y = 3.5, label = bquote(Warm:~italic(b)[FAS] == .(FAS_slope_w)),
           size=4, hjust=0, family="Helvetica", color = cols.fas[1])+
  annotate("text", label = paste("n = ", nrow(data.fasER), sep=""),
           x = 4, y = 3.1, size=3, hjust=0, family="Helvetica", color = "black")+
  annotate("text", label = paste("n = ", nrow(data.fas.test), sep=""),
           x = 4, y = 2.85, size=3, hjust=0, family="Helvetica", color = cols.fas[1])+
  scale_fill_gradient(high = cols.fas[1], low = cols.fas[5])+
  scale_color_gradient(high = cols.fas[1], low = cols.fas[5])+
  ylim(0,4)
ggformat(FASmodel_plot1, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*FAS),
         size_text = 12, print = FALSE)

# AS 
ASmodel_plot1<-ggplot(data=data.asER, aes(x=lnBWg, y=lnAS)) +
  geom_point(alpha=0.9, size=1, pch=1, color="grey75")+
  geom_line(data=data.plotAS_ER[round(data.plotAS_ER$tempTest,2)==25,], 
            mapping=aes(y = model_predFE, x=lnBWg,  group=tempTest, color= tempTest),
            color="black", linewidth=0.5, lty=1, alpha=0.8, show.legend=FALSE) +
  geom_point(data=data.as.test, aes(x=lnBWg, y=lnAS),
             alpha=0.9, size=1, pch=21, show.legend = FALSE, stroke =0.2,
             fill = cols.as[3], color = cols.as[1])+
  geom_line(data=data.plotAS_warm[round(data.plotAS_warm$tempTest,1)==25,],
            aes(y = model_predFE, x=lnBWg,  group=tempTest),
            color=cols.as[1], linewidth=1, lty=1, show.legend=FALSE) +
  annotate("text",  x = -5.2, y = 11.5, label = bquote(Optimal:~italic(b)[AS] == .(AS_slope)),
           size=4, hjust=0, family="Helvetica", color = "black")+
  annotate("text",  x = -5.2, y = 9.8, label = bquote(Warm:~italic(b)[AS] == .(AS_slope_w)),
           size=4, hjust=0, family="Helvetica", color = cols.as[1])+
  scale_fill_gradient( low = cols.as[5], high = cols.as[1])+
  scale_color_gradient( low = cols.as[5], high = cols.as[1])+
  ylim(x = -6.3, 12)+
  xlim(x = -6.3, 12)+
  annotate("text", label = paste("n = ", nrow(data.asER), sep=""),
           x = -5.2, y = 8.5, size=3, hjust=0, family="Helvetica", color = "black")+
  annotate("text", label = paste("n = ", nrow(data.as.test), sep=""),
           x = -5.2, y = 7.6, size=3, hjust=0, family="Helvetica", color = cols.as[1])
ggformat(ASmodel_plot1, x_title=expression(italic(ln)*Body~mass~(g)),
         y_title=expression(italic(ln)*AS~(mg~O[2]~h^-1)), size_text = 12, print = F)


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
ggsave(filename = paste("./Figures/Figure3.png", sep=""),
       plot=scaling, width = 6.8, height = 6.8, units = "in")



# Global outcomes by temprature  --------

AMRmodel_plot1_t<-ggplot(data=data.amrER, aes(x=tempTest, y=lnAMR)) +
  geom_point(alpha=0.9, size=1, pch=1, color="grey75")+
  geom_line(data=data.plotAMRint_ER[
                                round(data.plotAMRint_ER$lnBWg,2)==-2.54 |
                                round(data.plotAMRint_ER$lnBWg,2)== 2.46 | # ~ 10
                                round(data.plotAMRint_ER$lnBWg,2)== 4.46 | # ~ 100; 119.104 grams
                                round(data.plotAMRint_ER$lnBWg,2)== 6.96
                                ,], # ~ 1000
            mapping=aes(y = model_predFE, x=tempTest,  group=lnBWg, color= lnBWg),
            # color="black",
            linewidth=0.5, lty=1, alpha=0.8, show.legend=FALSE) +
  geom_point(data=data.amr.test, aes(x=tempTest, y=lnAMR),
             alpha=0.9, size=1, pch=21, show.legend = FALSE, stroke =0.2,
             fill = cols.amr[3], color = cols.amr[1])+
  geom_line(data=data.plotAMR_warm[round(data.plotAMR_warm$lnBWg,2)== 4.88,],
            aes(y = model_predFE, x=tempTest,  group=lnBWg),
            color=cols.amr[1], linewidth=1, lty=1, show.legend=FALSE) +
  scale_fill_gradient( low = cols.amr[5], high = cols.amr[1])+
  scale_color_gradient( low = "grey", high = "black")+
  ylim(x = -6.3, 12)+
  xlim(x = 0,40)
ggformat(AMRmodel_plot1_t, x_title="Temperature ºC",
         y_title=expression(italic(ln)*MMR~(mg~O[2]~h^-1)), size_text = 12, print = F)


# RMR 
RMRmodel_plot1_t<-ggplot(data=data.rmrER, aes(x=tempTest, y=lnRMR)) +
  geom_point(alpha=0.9, size=1, pch=1, color="grey75")+
  geom_line(data=data.plotRMR_ER[round(data.plotRMR_ER$lnBWg,2)== 4.78 ,], # ~ 1000
            mapping=aes(y = model_predFE, x=tempTest,  group=lnBWg, color= tempTest),
            color="black", linewidth=0.5, lty=1, alpha=0.8, show.legend=FALSE) +
  geom_point(data=data.rmr.test, aes(x=tempTest, y=lnRMR),
             alpha=0.9, size=1, pch=21, show.legend = FALSE, stroke =0.2,
             fill = cols.rmr[3], color = cols.rmr[1])+
  geom_line(data=data.plotRMR_warm[round(data.plotRMR_warm$lnBWg,2)== 4.88,],
            aes(y = model_predFE, x=tempTest,  group=lnBWg),
            color=cols.rmr[1], linewidth=1, lty=1, show.legend=FALSE) +
  scale_fill_gradient( low = cols.rmr[5], high = cols.rmr[1])+
  scale_color_gradient( low = cols.rmr[5], high = cols.rmr[1])+
  ylim(x = -6.3, 12)+
  xlim(x = 0,40)
ggformat(RMRmodel_plot1_t, x_title="Temperature ºC",
         y_title=expression(italic(ln)*RMR~(mg~O[2]~h^-1)), size_text = 12, print = F)



# AS 
ASmodel_plot1_t<-ggplot(data=data.asER, aes(x=tempTest, y=lnAS)) +
  geom_point(alpha=0.9, size=1, pch=1, color="grey75")+
  geom_line(data=data.plotAS_ER[round(data.plotAS_ER$lnBWg,2)== 4.78 ,], # ~ 1000
            mapping=aes(y = model_predFE, x=tempTest,  group=lnBWg, color= tempTest),
            color="black", linewidth=0.5, lty=1, alpha=0.8, show.legend=FALSE) +
  geom_point(data=data.as.test, aes(x=tempTest, y=lnAS),
             alpha=0.9, size=1, pch=21, show.legend = FALSE, stroke =0.2,
             fill = cols.as[3], color = cols.as[1])+
  geom_line(data=data.plotAS_warm[round(data.plotAS_warm$lnBWg,2)== 4.88,],
            aes(y = model_predFE, x=tempTest,  group=lnBWg),
            color=cols.as[1], linewidth=1, lty=1, show.legend=FALSE) +
  scale_fill_gradient( low = cols.as[5], high = cols.as[1])+
  scale_color_gradient( low = cols.as[5], high = cols.as[1])+
  ylim(x = -6.3, 12)+
  xlim(x = 0,40)+
  geom_point(data=data.as.test[grepl(x = data.as.test$Common_name, pattern = "salmon"),],
             aes(x=tempTest, y=lnAS),
             alpha=1, size=2, pch=21, show.legend = FALSE, stroke =0.2,
             fill = "red", color = "black")  
ggformat(ASmodel_plot1_t, x_title="Temperature ºC",
         y_title=expression(italic(ln)*AS~(mg~O[2]~h^-1)), size_text = 12, print = F)


# FAS 
FASmodel_plot1_t<-ggplot(data=data.fasER, aes(x=tempTest, y=lnFAS)) +
  geom_point(alpha=0.9, size=1, pch=1, color="grey75")+
  geom_line(data=data.plotFAS_ER[round(data.plotFAS_ER$lnBWg,2)== 4.78 ,], # ~ 1000
            mapping=aes(y = model_predFE, x=tempTest,  group=lnBWg, color= tempTest),
            color="black", linewidth=0.5, lty=1, alpha=0.8, show.legend=FALSE) +
  geom_point(data=data.fas.test, aes(x=tempTest, y=lnFAS),
             alpha=0.9, size=1, pch=21, show.legend = FALSE, stroke =0.2,
             fill = cols.fas[3], color = cols.fas[1])+
  geom_line(data=data.plotFAS_warm[round(data.plotFAS_warm$lnBWg,2)== 4.88,],
            aes(y = model_predFE, x=tempTest,  group=lnBWg),
            color=cols.fas[1], linewidth=1, lty=1, show.legend=FALSE) +
  scale_fill_gradient( low = cols.fas[5], high = cols.fas[1])+
  scale_color_gradient( low = cols.fas[5], high = cols.fas[1])+
  ylim(x = -1, 6)+
  xlim(x = 0,40)+
  geom_point(data=data.fas.test[grepl(x = data.fas.test$Common_name, pattern = "salmon"),],
             aes(x=tempTest, y=lnFAS),
             alpha=1, size=2, pch=21, show.legend = FALSE, stroke =0.2,
             fill = "red", color = "black")  
ggformat(FASmodel_plot1_t, x_title="Temperature ºC",
         y_title=expression(italic(ln)*FAS~(mg~O[2]~h^-1)), size_text = 12, print = F)


scaling_t<-cowplot:::plot_grid(AMRmodel_plot1_t, RMRmodel_plot1_t,
                             ASmodel_plot1_t, FASmodel_plot1_t,
                              align = "hv",
                              axis = "l",
                              nrow = 2,
                              ncol = 2,
                              labels = "AUTO",
                             label_x = c(0.18, 0.18),
                             label_y = c(0.895, 0.895),
                              label_size = 12)
# scaling
ggsave(filename = paste("./Figures/Suppl_Figure_t.png", sep=""),
       plot=scaling_t, width = 6.8, height = 6.8, units = "in")

# Both Activation energies together: (NOT USED) -------
# model_predFE << is in lnMR units 
# MMR
# AMRmodel_plot3.E_ALL<-ggplot(data.amrER[data.amrER$lnBWg==0,]) +
#   geom_point(data.amrER, mapping=aes(x=tempTest, y=log(mass_specamr), fill= tempTest, size = lnBWg), color="grey50", fill="grey70", alpha=1, pch=21, show.legend = FALSE)+
#   geom_point(data.amr.test, mapping=aes(x=tempTest, y=log(mass_specamr), fill= tempTest, size = lnBWg), color="black",  alpha=1, pch=21,  show.legend = FALSE)+
#   scale_fill_gradient( low = cols.amr[3], high = cols.amr[1])+
#   scale_color_gradient( low = "grey80", high = "grey0")+
#   ylim(-5.7,4)+
#   annotate("text",  x = 3.27, y = 3.5, label = bquote(italic(E)[MMR] == change~with~mass),size=5, hjust=0, family="Helvetica", color = "black")+
#   annotate("text",  x = 3.27, y = 2.5, label = bquote(italic(E)[MMR] == .(round(MMR_E_W_eV,3))),size=5, hjust=0, family="Helvetica", color = cols.amr[1])+
#   annotate("text",  x = 3.27, y = -3.2, label = bquote(1*g~italic(E)== .(round(MMR_E_ER$MMR_E_ER_eV[1],3))),size=4, hjust=0, family="Helvetica", color = "black")+
#   annotate("text",  x = 3.27, y = -4, label = bquote(10*g~italic(E) == .(round(MMR_E_ER$MMR_E_ER_eV[2],3))),size=4, hjust=0, family="Helvetica", color = "black")+
#   annotate("text",  x = 3.27, y = -4.8, label = bquote(100*g~italic(E) == .(round(MMR_E_ER$MMR_E_ER_eV[3],3))),size=4, hjust=0, family="Helvetica", color = "black")+
#   annotate("text",  x = 3.27, y = -5.6, label = bquote(1000*g~italic(E) == .(round(MMR_E_ER$MMR_E_ER_eV[4],3))),size=4, hjust=0, family="Helvetica", color = "black")+
#   scale_x_continuous(limits = c(3.25, 3.661), sec.axis = sec_axis(~ ((1000/.))-273.15, name = expression(Temperature~degree*C), breaks = c(32, 25, 17, 10, 3 )))+
#   geom_line(data = data.plotAMR_warm[which(round(data.plotAMR_warm$lnBWg, 1) == round(log(1.4), 1)),],
#             aes(y = log((exp(model_predFE)/exp(lnBWg))), x=tempTest, group=lnBWg) ,color=cols.amr[2], linewidth=1, lty=1, show.legend=FALSE)+
#   geom_line(data = data.plotAMRint_ER[exp(data.plotAMRint_ER$lnBWg) <= 1000, ], aes(y = log((exp(model_predFE)/exp(lnBWg))), x=tempTest, group=lnBWg, color = lnBWg), linewidth=0.5, lty=1, show.legend=FALSE)
# ggformat(AMRmodel_plot3.E_ALL, x_title= expression(Temperature^-1~(1000/K)), y_title=expression(italic(ln)*MMR~(mg~O[2]~h^-1~g^-1)), print = F)
# 
# RMRmodel_plot3.E_ALL<-ggplot(data.rmrER[data.rmrER$lnBWg==0,]) +
#   geom_point(data.rmrER, mapping=aes(x=tempTest, y=log(mass_specrmr), fill= tempTest, size = lnBWg), color="grey50", fill="grey70", alpha=1, pch=21)+
#   geom_point(data.rmr.test, mapping=aes(x=tempTest, y=log(mass_specrmr), fill= tempTest, size = lnBWg), color = "black", alpha=1, pch=21, show.legend = FALSE)+
#   scale_fill_gradient( low = cols.rmr[3], high = cols.rmr[1])+
#   annotate("text",  x = 3.27, y = 3.5, label = bquote(italic(E)[RMR] == .(round(RMR_E_ER_eV,3))),size=5, hjust=0, family="Helvetica", color = "black")+
#   annotate("text",  x = 3.27, y = 2.5, label = bquote(italic(E)[RMR] == .(round(RMR_E_W_eV,3))),size=5, hjust=0, family="Helvetica", color = cols.rmr[1])+
#   ylim(-5.7,4)+
#   scale_size_continuous(name="Mass (g)",
#                       breaks=c(log(0.1), log(10), log(1000)),
#                       labels=c("0.1", "10", "1000"))+
#   scale_x_continuous(limits = c(3.25, 3.661), sec.axis = sec_axis(~ ((1000/.))-273.15, name = expression(Temperature~degree*C), breaks = c(32, 25, 17, 10, 3 )))+
#   geom_line(data = data.plotRMR_ER[which(round(data.plotRMR_ER$lnBWg, 1) == round(log(1.35), 1)),], aes(y = log((exp(model_predFE)/exp(lnBWg))), x=tempTest, group=lnBWg ) ,color="black", linewidth=1, lty=1, show.legend=FALSE)+
#   geom_line(data = data.plotRMR_warm[which(round(data.plotRMR_warm$lnBWg, 1) == round(log(1.35), 1)),], aes(y = log((exp(model_predFE)/exp(lnBWg))), x=tempTest, group=lnBWg), color = cols.rmr[2], linewidth=1, lty=1, show.legend=FALSE)
# ggformat(RMRmodel_plot3.E_ALL, x_title= expression(Temperature^-1~(1000/K)), y_title=expression(italic(ln)*RMR~(mg~O[2]~h^-1~g^-1)), print = FALSE)
# RMRmodel_plot3.E_ALL <- RMRmodel_plot3.E_ALL + theme(legend.position = c(0.87, 0.82))
# 
# ASmodel_plot3.E_ALL<-ggplot(data.asER[data.asER$lnBWg==0,]) +
#   geom_point(data.asER, mapping=aes(x=tempTest, y=log(mass_specas), fill= tempTest, size = lnBWg), color="grey50", fill="grey70", alpha=1, pch=21, show.legend = FALSE)+
#   geom_point(data.as[!c(data.as$test_category=="ecol_relev"),], mapping=aes(x=tempTest, y=log(mass_specas), fill= tempTest, size = lnBWg), color = "black", alpha=1, pch=21, show.legend = FALSE)+
#   scale_fill_gradient( low = cols.as[3], high = cols.as[1])+
#   annotate("text",  x = 3.27, y = 3.5, label = bquote(italic(E)[AS] == .(round(AS_E_ER_eV,3))),size=5, hjust=0, family="Helvetica", color = "black")+
#   annotate("text",  x = 3.27, y = 2.5, label = bquote(italic(E)[AS] == .(round(AS_E_W_eV,3))),size=5, hjust=0, family="Helvetica", color = cols.as[1])+
#   ylim(-5.7,4)+
#   scale_x_continuous(limits = c(3.25, 3.661), sec.axis = sec_axis(~ ((1000/.))-273.15, name = expression(Temperature~degree*C), breaks = c(32, 25, 17, 10, 3 )))+
#   geom_line(data = data.plotAS_ER[which(round(data.plotAS_ER$lnBWg, 1) == round(log(1.35), 1)),], aes(y = log((exp(model_predFE)/exp(lnBWg))), x=tempTest, group=lnBWg ) ,color="black", linewidth=1, lty=1, show.legend=FALSE)+
#   geom_line(data = data.plotAS_warm[which(round(data.plotAS_warm$lnBWg, 1) == round(log(1.4), 1)),], aes(y = log((exp(model_predFE)/exp(lnBWg))), x=tempTest, group=lnBWg), color = cols.as[2], linewidth=1, lty=1, show.legend=FALSE)
# ggformat(ASmodel_plot3.E_ALL, x_title= expression(Temperature^-1~(1000/K)), y_title=expression(italic(ln)*AS~(mg~O[2]~h^-1~g^-1)), print = F)
# 
# Arh.plot<-cowplot:::plot_grid(AMRmodel_plot3.E_ALL, RMRmodel_plot3.E_ALL, ASmodel_plot3.E_ALL, 
#           align = "hv",
#           axis = "l",
#           nrow = 1,
#           ncol = 3,
#           labels = "AUTO",
#           label_y = 0.95,
#           label_size = 17)
# ggsave(filename = paste("./Figures/Figure5.png", sep=""),
#        plot=Arh.plot, width = 13, height = 5, units = "in")


# Final plots with conceptual figure ---------
# All fish together:
MRmodel_plot_inset<-
  ggplot(data = sum_CItable[sum_CItable$var_repeat == "lnBWg" &  sum_CItable$MR == "RMR" &
                              sum_CItable$temp_cat == "ER",])+
  geom_linerange(aes(ymin = ci5, ymax = ci95, x = MR))+
  geom_point(mapping = aes(x = "MMR", y = AMR_slope31), size = 2, pch = 21,
             color = "black", stroke = 0.5,fill = "#5d7d7c")+
  geom_point(mapping = aes(x = "MMR", y = AMR_slope5), size = 2, pch = 21,
             color = "black", stroke = 0.5,fill = "#c4c7cc")+
  geom_point(mapping = aes(x = "MMR", y = AMR_slope15), size = 2, pch = 21,
             color = "black", stroke = 0.5,fill = "#a2adb5")+
  geom_point(mapping = aes(x = "MMR", y = AMR_slope25), size = 2, pch = 21,
             color = "black", stroke = 0.5, fill = "#7e959b")+
  geom_point(mapping = aes(x = "RMR", y = RMR_slope), size = 2)+
  annotate(geom = "text", x = 1.25, y = 0.7, label = "5", size = 2, angle = 45, hjust =0)+
  annotate(geom = "text", x = 1.25, y = 0.75, label = "15", size = 2, angle = 45, hjust =0)+
  annotate(geom = "text", x = 1.25, y = 0.80, label = "25", size = 2, angle = 45, hjust =0)+
  annotate(geom = "text", x = 1.25, y = 0.84, label = "31ºC", size = 2, angle = 45, hjust =0)+
  theme_classic()+
  coord_flip()+
  ylim(0.5, 1)+
  xlab("")+
  ylab("Slope value")+
  geom_hline(yintercept = c(0.75, 1), linetype = "dashed", linewidth = 0.1)+
  theme(axis.text = element_text(size = 8, family = "Helvetica"),
        axis.title.x = element_text(size = 9, family = "Helvetica"),
      panel.background = element_rect(fill = "transparent",
                                      colour = NA_character_), # necessary to avoid drawing panel outline
      panel.grid.major = element_blank(), # get rid of major grid
      panel.grid.minor = element_blank(), # get rid of minor grid
      plot.background = element_rect(fill = "transparent",
                                     colour = NA_character_), # necessary to avoid drawing plot outline
      legend.background = element_rect(fill = "transparent"),
      legend.box.background = element_rect(fill = "transparent"),
      legend.key = element_rect(fill = "transparent"))


MRmodel_plot2<-
  ggplot(data=data.rmrER, aes(x=lnBWg, y=lnRMR)) +
  geom_line(data=data.plotRMR_ER[round(data.plotRMR_ER$tempTest,2)==25,],
            aes(y = model_predFE, x=lnBWg,  group=tempTest, color= tempTest),
            color="black", linewidth=0.7, lty=1, show.legend=FALSE) +
  geom_line(data=data.plotAMRint_ER[round(data.plotAMRint_ER$tempTest,2)==5 | 
                                        round(data.plotAMRint_ER$tempTest,2)==15 |
                                        round(data.plotAMRint_ER$tempTest,2)==25 |
                                        round(data.plotAMRint_ER$tempTest,2)==31,],
            aes(y = model_predFE, x=lnBWg,  group=factor(tempTest), color= factor(tempTest)),
            size=0.7, lty=1, show.legend=FALSE) +
  scale_color_manual(values = c("#c4c7cc", "#a2adb5", "#7e959b","#5d7d7c"))+
  ylim(x = -6.3, 12)+
  xlim(x = -6.3, 12)+
  annotate("text",  x = -5.8, y = 11.8, label = "B", fontface =2,
           size=4.5, hjust=1, family="Helvetica", color = "black")+
  # annotate("text",  x = 9, y = 8.4, label = expression(paste("MMR")),
  #          size=5, hjust=0, family="Helvetica", color = "black", angle = 42)+
  # annotate("text",  x = 9, y = 4.1, label = expression(paste("RMR")),
  #          size=5, hjust=0, family="Helvetica", color = "black", angle = 42)+
  annotate("text",  x = -1, y = -5.9, label = expression(paste("OPTIMAL TEMPERATURES")),
           size=3, hjust=0, family="Helvetica", color = "black", angle = 0)+
  xlab(expression(italic(ln)*Body~mass~(g)))+
  ylab(expression(italic(ln)*MR~(mg~O[2]~h^-1)))+
  theme(axis.text.y=element_text(size=12, colour= 'black'),
		axis.text.x=element_text(size=12, colour= 'black'),
		axis.line.y=element_line(colour = 'black',linewidth=0.5),
		axis.line.x=element_line(colour = 'black',linewidth=0.5),
		axis.ticks.y=element_line(size=0.5),
		panel.background = element_blank(),
		axis.ticks.x.bottom = element_line(linewidth=0.5, colour = "black"),
	  axis.title.y=element_text(size=12),
		axis.title.x=element_text(size=12),
		panel.border = element_rect(linetype = "solid",fill=NA, colour = "black"))+
  inset_element(MRmodel_plot_inset, 0.01, 0.6, 0.7, 1)
MRmodel_plot2
  
# warm final plots 
# # All fish together:
MRmodel_plot_inset_w<-
  ggplot(data = sum_CItable[sum_CItable$var_repeat == "lnBWg" &
                              c(sum_CItable$MR == "RMR" | sum_CItable$MR == "MMR") &
                              sum_CItable$temp_cat == "W",])+
  geom_hline(yintercept = c(0.75, 1), linetype = "dashed", linewidth = 0.1)+
  geom_linerange(aes(ymin = ci5, ymax = ci95, x = MR, color = MR), show.legend = F)+
  geom_point(mapping = aes(x = "MMR", y = AMR_slope_w), size = 2, color = cols.amr[3])+
  geom_point(mapping = aes(x = "RMR", y = RMR_slope_w), size = 2, color = cols.rmr[3])+
  theme_classic()+
  coord_flip()+
  scale_color_manual(values = c(cols.amr[2], cols.rmr[2]))+
  ylim(0.5, 1)+
  xlab("")+
  ylab("Slope value")+
  theme(axis.text = element_text(size = 8, family = "Helvetica"),
        axis.title.x = element_text(size = 9, family = "Helvetica"),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_), # necessary to avoid drawing plot outline
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"))

MRmodel_plot2_w<-
  ggplot(data=data.rmr.test, aes(x=lnBWg, y=lnRMR)) +
  geom_line(data=data.plotRMR_warm[round(data.plotRMR_warm$tempTest,2)==25,],
            aes(y = model_predFE, x=lnBWg,  group=tempTest, color= tempTest),
            linewidth=0.7, lty=1, show.legend=FALSE, color = cols.rmr[2]) +
  geom_line(data=data.plotAMR_warm[round(data.plotAMR_warm$tempTest,2)==25,],
            aes(y = model_predFE, x=lnBWg,  group=tempTest),
            linewidth=0.7, lty=1, show.legend=FALSE, color = cols.amr[2]) +
  scale_color_gradient( low = "grey", high = "black")+
  ylim(x = -6.3, 12)+
  xlim(x = -6.3, 12)+
  annotate("text",  x = -5.8, y = 11.8, label = "D", fontface =2,
           size=4.5, hjust=1, family="Helvetica", color = "black")+
  # annotate("text",  x = 6.2, y = 6, label = expression(paste("MMR")),
  #          size=5, hjust=0, family="Helvetica", angle = 30, color = cols.amr[1])+
  # annotate("text",  x = 6.7, y = 3.1, label = expression(paste("RMR")),
  #          size=5, hjust=0, family="Helvetica", angle = 42, color = cols.rmr[1])+
  annotate("text",  x = -1, y = -5.9, label = expression(paste("WARM TEMPERATURES")),
           size=3, hjust=0, family="Helvetica", color = "black", angle = 0)+
  xlab(expression(italic(ln)*Body~mass~(g)))+
  ylab(expression(italic(ln)*MR~(mg~O[2]~h^-1)))+
  theme(axis.text.y=element_text(size=12, colour= 'black'),
		axis.text.x=element_text(size=12, colour= 'black'),
		axis.line.y=element_line(colour = 'black',linewidth=0.5),
		axis.line.x=element_line(colour = 'black',linewidth=0.5),
		axis.ticks.y=element_line(size=0.5),
		panel.background = element_blank(),
		axis.ticks.x.bottom = element_line(linewidth=0.5, colour = "black"),
	  axis.title.y=element_text(size=12),
		axis.title.x=element_text(size=12),
		panel.border = element_rect(linetype = "solid",fill=NA, colour = "black"))+
  inset_element(MRmodel_plot_inset_w, 0.01, 0.6, 0.7, 1)
# MRmodel_plot2_w

# conceptual figures 
data_sim<-read.csv(here("./Data/Conceptual_fig_data.csv"))

psim<-ggplot(data_sim[data_sim$simID==3,], aes(x=log(BW_kg), y=log(pred.mmr.mgO2min) ))+
  geom_abline(slope = data_sim[data_sim$simID==3,"slope.MMR"][1],
              intercept = data_sim[data_sim$simID==3,"int.MMR"][1],
              color = "#7e959b", linewidth=0.7)+ # MMR
  geom_abline(slope = data_sim[data_sim$simID==3,"slope.SMR"][1],
              intercept = data_sim[data_sim$simID==3,"int.SMR"][1],
              color = "black", linewidth=0.7)+ # MMR
  lims(x = c(-4.5, 4.5), y = c(-4.5, 4.5))+
  annotate("text",  x = -4.3, y = 4.4, label = "A", fontface =2,
           size=4.5, hjust=1, family="Helvetica", color = "black")+
  annotate("text",  x = 0.8, y = 3.5, label = expression(paste("MMR")),
       size=5, hjust=0, family="Helvetica", angle = 42, color = "#7e959b")+
  annotate("text",  x = 1, y = 0.9, label = expression(paste("RMR")),
           size=5, hjust=0, family="Helvetica", angle = 33, color = "black")+
  annotate("text",  x = -2, y = -4.15, label = expression(paste("OPTIMAL TEMPERATURES")),
           size=3, hjust=0, family="Helvetica", color = "black", angle = 0)+
  annotate("text",  x = -3.8, y = 4.3, label = expression(paste("HYPOTHESIS")),
           size=3, hjust=0, family="Helvetica", color = "black", angle = 0)+
  annotate("text",  x = -3.8, y = 3.6, label = expression(italic(b)[MMR]~">"~italic(b)[RMR]),
           size=3, hjust=0, family="Helvetica", color = "black", angle = 0)+
  annotate("text",  x = 2, y = -0.8,
           label = expression(italic(b)[MMR]~"="~1.00),
           size=3, hjust=0, family="Helvetica", color = "black", angle = 0)+
  annotate("text",  x = 2, y = -1.5,
           label = expression(italic(b)[RMR]~"="~0.89),
           size=3, hjust=0, family="Helvetica", color = "black", angle = 0)+
  annotate("text",  x = 2, y = -2.2,
           label = expression(italic(b)[AS]~"="~1.05),
           size=3, hjust=0, family="Helvetica", color = "black", angle = 0)+
  annotate("text",  x = 2, y = -2.9,
           label = expression(italic(b)[FAS]~"="~0.11),
           size=3, hjust=0, family="Helvetica", color = "black", angle = 0)+
  xlab(expression(italic(ln)*Body~mass~(g)))+
  ylab(expression(italic(ln)*MR~(mg~O[2]~h^-1)))+
  theme(axis.text.y=element_text(size=12, colour= 'black'),
		axis.text.x=element_text(size=12, colour= 'black'),
		axis.line.y=element_line(colour = 'black',size=0.5),
		axis.line.x=element_line(colour = 'black',size=0.5),
		axis.ticks.y=element_line(size=0.5),
		panel.background = element_blank(),
		axis.ticks.x.bottom = element_line(size=0.5, colour = "black"),
	  axis.title.y=element_text(size=12),
		axis.title.x=element_text(size=12),
		panel.border = element_rect(linetype = "solid",fill=NA, colour = "black"))
# psim
    
# conceptual WArm 
psim_w<-ggplot(data_sim[data_sim$simID==1,], aes(x=log(BW_kg), y=log(pred.mmr.mgO2min) ))+
  geom_abline(slope = data_sim[data_sim$simID==1,"slope.MMR"][1],
              intercept = data_sim[data_sim$simID==1,"int.MMR"][1],
              color = cols.amr[2], linewidth=0.7)+ # MMR
  geom_abline(slope = data_sim[data_sim$simID==1,"slope.SMR"][1],
              intercept = data_sim[data_sim$simID==1,"int.SMR"][1],
              color = cols.rmr[2], linewidth=0.7)+ # MMR
  lims(x = c(-4.5, 4.5), y = c(-4.5, 4.5))+
  annotate("text",  x = -4.3, y = 4.4, label = "C", fontface =2,
           size=4.5, hjust=1, family="Helvetica", color = "black")+
  annotate("text",  x = 0.8, y = 3.3, label = expression(paste("MMR")),
       size=5, hjust=0, family="Helvetica", angle = 32, color = cols.amr[2])+
  annotate("text",  x = 1, y = 0.9, label = expression(paste("RMR")),
           size=5, hjust=0, family="Helvetica", angle = 47, color = cols.rmr[2])+
  annotate("text",  x = -2, y = -4.15, label = expression(paste("WARM TEMPERATURES")),
           size=3, hjust=0, family="Helvetica", color = "black", angle = 0)+
  annotate("text",  x = -3.8, y = 4.3, label = expression(paste("HYPOTHESIS")),
           size=3, hjust=0, family="Helvetica", color = "black", angle = 0)+
  annotate("text",  x = -3.8, y = 3.6, label = expression(italic(b)[MMR]~"<"~italic(b)[RMR]),
           size=3, hjust=0, family="Helvetica", color = "black", angle = 0)+
  annotate("text",  x = 2, y = -0.8,
           label = expression(italic(b)[MMR]~"="~0.75),
           size=3, hjust=0, family="Helvetica", color = "black", angle = 0)+
  annotate("text",  x = 2, y = -1.5,
           label = expression(italic(b)[RMR]~"="~0.95),
           size=3, hjust=0, family="Helvetica", color = "black", angle = 0)+
  annotate("text",  x = 2, y = -2.2,
           label = expression(italic(b)[AS]~"="~0.69),
           size=3, hjust=0, family="Helvetica", color = "black", angle = 0)+
  annotate("text",  x = 2, y = -2.9,
           label = expression(italic(b)[FAS]~"="~-0.14),
           size=3, hjust=0, family="Helvetica", color = "black", angle = 0)+
  xlab(expression(italic(ln)*Body~mass~(g)))+
  ylab(expression(italic(ln)*MR~(mg~O[2]~h^-1)))+
  theme(axis.text.y=element_text(size=12, colour= 'black'),
		axis.text.x=element_text(size=12, colour= 'black'),
		axis.line.y=element_line(colour = 'black',size=0.5),
		axis.line.x=element_line(colour = 'black',size=0.5),
		axis.ticks.y=element_line(size=0.5),
		panel.background = element_blank(),
		axis.ticks.x.bottom = element_line(size=0.5, colour = "black"),
	  axis.title.y=element_text(size=12),
		axis.title.x=element_text(size=12),
		panel.border = element_rect(linetype = "solid",fill=NA, colour = "black"))
# psim_w
# patchwork library 
top<-  psim + MRmodel_plot2+ plot_layout(axis_titles = "collect") 
# top
bottom <- psim_w + MRmodel_plot2_w +  plot_layout(axis_titles = "collect")
  
# plot_layout(widths = c(2, 1), heights = unit(c(5, 1), c('cm', 'null')))

# save the plots 
mainplot<-
  cowplot::plot_grid(top, bottom,
                   align = "hv", nrow = 2)

ggsave(filename = paste("./Figures/FigureMAIN.png", sep=""),
       plot=mainplot, width = 7, height = 7, units = "in")

# --- misc -----
# fas_temp<-ggplot(data=data.fas, aes(y=FAS, x = tempTest))+
#   geom_point(show.legend = F)+
#   ylim(0,20)+
#   facet_grid(DemersPelag~test_category3)+
#   geom_smooth(method = "lm")
# ggformat(fas_temp, y_title = expression(FAS~(MMR/RMR)), x_title = "Temperature ºC", print = F)
# fas_temp<-fas_temp+theme(legend.position = "none")
# fas_temp

# Not used in manuscript -------------
# Overall all fish together, with booted CI:
# MRmodel_plot1<-ggplot() +
#   geom_ribbon(data=data.plotRMR_ER, mapping = aes(y = model_predFE, x=lnBWg, ymin=CI_2.5, ymax=CI_97.5, group = tempTestK1000), linetype=2, alpha=0.1, fill = "grey30")+
#   geom_ribbon(data.plotAMRint_ER, mapping = aes(y = model_predFE, x=lnBWg, ymin=CI_2.5, ymax=CI_97.5, group = tempTestK1000), linetype=2, alpha=0.1, fill = "grey30")+
#   geom_line(data=data.plotAMRint_ER[round(data.plotAMRint_ER$tempTestK1000,2)==3.49,], mapping = aes(y = model_predFE, x=lnBWg,  group=tempTestK1000), color="black", linewidth=0.7, lty=1, show.legend=FALSE) +
#   geom_line(data=data.plotRMR_ER[round(data.plotRMR_ER$tempTestK1000,2)==3.49,], aes(y = model_predFE, x=lnBWg,  group=tempTestK1000, color= tempTestK1000), color="black", linewidth=0.7, lty=1, show.legend=FALSE) +
#   scale_y_continuous(limits = c(-7, 12), breaks = seq(-7,12, 2))+
#   scale_x_continuous(limits = c(-7, 12), breaks = seq(-7,12, 2))
# ggformat(MRmodel_plot1, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*MR~(mg~O[2]~h^-1)), print = F)
# 
# MRmodel_plotW<-ggplot() +
#   geom_ribbon(data.plotAMR_warm, mapping = aes(y = model_predFE, x=lnBWg, ymin=CI_2.5, ymax=CI_97.5, group = tempTestK1000), linetype=2, alpha=0.1, fill = cols.amr[2])+
#   geom_ribbon(data=data.plotRMR_warm, mapping = aes(y = model_predFE, x=lnBWg, ymin=CI_2.5, ymax=CI_97.5, group = tempTestK1000), linetype=2, alpha=0.1, fill = cols.rmr[2])+
#   geom_line(data=data.plotAMR_warm[round(data.plotAMR_warm$tempTestK1000,2)==3.49,], mapping = aes(y = model_predFE, x=lnBWg,  group=tempTestK1000), color=cols.amr[1], size=0.7, lty=1, show.legend=FALSE) +
#   geom_line(data=data.plotRMR_warm[round(data.plotRMR_warm$tempTestK1000,2)==3.49,], aes(y = model_predFE, x=lnBWg,  group=tempTestK1000, color= tempTestK1000), color=cols.rmr[1], size=0.7, lty=1, show.legend=FALSE) +
#   scale_y_continuous(limits = c(-7, 12), breaks = seq(-7,12, 2))+
#   scale_x_continuous(limits = c(-7, 12), breaks = seq(-7,12, 2))
# ggformat(MRmodel_plotW, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*MR~(mg~O[2]~h^-1)), print = F)