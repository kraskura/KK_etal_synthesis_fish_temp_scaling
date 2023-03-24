
library(lme4)
library(emmeans)
library(car)

library(evolvability) # Almer function
library(ape)
library(rotl)
library(ggtree)
library(Matrix)

library(pryr)
library(TDbook)
library(ggimage)

library(ggformatKK) # from github
library(weathermetrics)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(forcats)

library(here)

# **************************************************
# run full model comparison analysis ( run = TRUE), or figures, best models, only (run = FALSE)?

run = FALSE
# run = TRUE

# **************************************************

if(run){
  message("running all models and BIC comparison: this will take a while")
}else{
  message("not re-running all models and BIC comparison: using models run before")
}

## Source the data ------------
source("./R/get_data_temp.R") 
source("./R/mixed_model_outputs.R")

# function to order model selection based on the lowest BIC score
BICdelta<-function(BICtable){
  BIC.t <- BICtable [order(BICtable$BIC), ]
  BIC.t$delta <- round(abs(BIC.t$BIC[1] -  BIC.t$BIC), 5)
  return( BIC.t)
}

# function to get phylogenetic relatedness matrix for each data subset
get_phylo_matrix<-function(species.list, matrix.name, tree.name, dataset.ID){

  taxon_search <- tnrs_match_names(names=species.list, context_name="Vertebrates") 
  ott_in_tree <- ott_id(taxon_search)[is_in_tree(ott_id(taxon_search))]
  tr <- tol_induced_subtree(ott_ids = ott_in_tree)
  
  tr$tip.label <- strip_ott_ids(tr$tip.label, remove_underscores = TRUE)
  labels <- as.data.frame(tr$tip.label)

  labels$`tr$tip.label`[which(labels$`tr$tip.label` == "Oncorhynchus mykiss (species in domain Eukaryota)")]<-"Oncorhynchus mykiss"
  labels$`tr$tip.label`[which(labels$`tr$tip.label` == "Gadus morhua (species in domain Eukaryota)")]<-"Gadus morhua"
  labels$`tr$tip.label`[which(labels$`tr$tip.label` == "Leiocassis longirostris")]<-"Tachysurus dumerili" 
  labels$`tr$tip.label`[which(labels$`tr$tip.label` == "Tachysurus vachellii")]<-"Pseudobagrus vachellii" 
  labels$`tr$tip.label`[which(labels$`tr$tip.label` == "Rhinogobius similis")]<- "Rhinogobius giurinus" 
  tr$tip.label <- labels$`tr$tip.label`
  
  if(all(species.list %in% tr$tip.label)){
    # message(paste(dataset.ID, ": All species names are identified and mathced with phylo data", sep =""))
    message(paste(dataset.ID,": All species names are identified and mathced with phylo data \n",  "N species:", length(species.list), sep = ""))
  }
  
  tr2<-compute.brlen(tr)

  A <- ape::vcv.phylo(tr2)
  tree <- compute.brlen(tr2)
  cor <- vcv(tree, cor = T)
  
  # if((dim(A)[1] * dim(A)[2] == sum(!is.na(A))) && all(colnames(A) == rownames(A))){
  # 
  #   message("The phylo matrix has no missing values, and column/row names has matching species IDs") 
  #   message(paste("Size of the matrix ", dim(A)[1], " X ", dim(A)[2], sep =""))
  #  # The column names of A must be the species identifier. and the same
  # }  
  
  A <- Matrix::Matrix(ape::vcv(tree), sparse = TRUE)

  assign(matrix.name, value = A, envir = .GlobalEnv)
  assign(tree.name, value = tree, envir = .GlobalEnv)

}

# 1. Data for models ------
data.list<-get_data_temp(data.amr = "./Data/Fish_AMR_temp_dataset_mar2022.csv",
                                 data.rmr = "./Data//Fish_RMR_temp_dataset_mar2022.csv",
                                 ecology.data = "./Data/Kraskura_species_ecologies_mar2022.csv",
                                 onlyTop.above = TRUE, save.FishBase.species.data = F,
                                 calc_mass_specific = FALSE)

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
data.fas.test<-data.amr.test[c(!is.na(data.amr.test$FAS) & is.finite(data.amr.test$FAS)) , ]
data.as.test<-data.amr.test[c(!is.na(data.amr.test$lnAS) & is.finite(data.amr.test$lnAS)) , ]

# get model data set specific phylogentics model matrixes -------
data.rmrER<-droplevels(data.rmrER)
data.rmr.test<-droplevels(data.rmr.test)
data.amrER<-droplevels(data.amrER)
data.amr.test<-droplevels(data.amr.test)
data.asER<-droplevels(data.asER)
data.as.test<-droplevels(data.as.test)
data.fasER<-droplevels(data.fasER)
data.fas.test<-droplevels(data.fas.test)

# get model specific matrices
get_phylo_matrix(species.list = unique(levels(data.rmrER$species)), matrix.name = "A", tree.name = "tree.kk", dataset.ID = "RMR optimal")
get_phylo_matrix(species.list = unique(levels(data.rmr.test$species)), matrix.name = "A.rmr.w", tree.name = "tr.rmr.w", dataset.ID = "RMR warm")

get_phylo_matrix(species.list = unique(levels(data.amrER$species)), matrix.name = "A.mmr.er", tree.name = "tr.mmr.er", dataset.ID = "MMR optimal")
get_phylo_matrix(species.list = unique(levels(data.amr.test$species)), matrix.name = "A.mmr.w", tree.name = "tr.mmr.w", dataset.ID = "MMR warm")

get_phylo_matrix(species.list = unique(levels(data.asER$species)), matrix.name = "A.aas.er", tree.name = "tr.aas.er", dataset.ID = "AAS optimal")
get_phylo_matrix(species.list = unique(levels(data.as.test$species)), matrix.name = "A.aas.w", tree.name = "tr.aas.w", dataset.ID = "AAS warm")

get_phylo_matrix(species.list = unique(levels(data.fasER$species)), matrix.name = "A.fas.er", tree.name = "tr.fas.er", dataset.ID = "FAS optimal")
get_phylo_matrix(species.list = unique(levels(data.fas.test$species)), matrix.name = "A.fas.w", tree.name = "tree.fas.w", dataset.ID = "FAS warm")

if (run){
  # phylo models ---------
  ## RMR: ecologically relevant -----------------
  Phylo_RMR_model0 <- Almer(lnRMR ~ lnBWg + tempTestK1000 + (1|species) , data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model0.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model0intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (1|species), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model0int <- Almer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) , data=data.rmrER, REML=FALSE, A = list(species = A))
  
  Phylo_RMR_model1 <- Almer(lnRMR ~ lnBWg + tempTestK1000 + (1|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model1.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model1intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model1int <- Almer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  
  Phylo_RMR_model2 <- Almer(lnRMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model2int <- Almer(lnRMR ~ lnBWg * tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model2.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model2intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTestK1000,2)  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  
  Phylo_RMR_model4 <- Almer(lnRMR ~ lnBWg + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model4int <- Almer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model4.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model4intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  
  Phylo_RMR_model5 <- Almer(lnRMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model5int <- Almer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model5.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model5intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  
  Phylo_RMR_model6 <- Almer(lnRMR ~ lnBWg + tempTestK1000 + (lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model6int <- Almer(lnRMR ~ lnBWg * tempTestK1000+ (lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model6.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model6intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  
  ### BIC --------------
  RMR_BIC<-BICdelta(BIC(Phylo_RMR_model0, Phylo_RMR_model0int, Phylo_RMR_model0.POLY, Phylo_RMR_model0intPOLY,
                  Phylo_RMR_model1, Phylo_RMR_model1int, Phylo_RMR_model1.POLY, Phylo_RMR_model1intPOLY,
                  Phylo_RMR_model2, Phylo_RMR_model2int, Phylo_RMR_model2.POLY, Phylo_RMR_model2intPOLY,
                  Phylo_RMR_model4, Phylo_RMR_model4int, Phylo_RMR_model4.POLY, Phylo_RMR_model4intPOLY,
                  Phylo_RMR_model5, Phylo_RMR_model5int, Phylo_RMR_model5.POLY, Phylo_RMR_model5intPOLY,
                  Phylo_RMR_model6, Phylo_RMR_model6int, Phylo_RMR_model6.POLY, Phylo_RMR_model6intPOLY))
  
  ## MMR ecologically relevant -----------------
  Phylo_MMR_model0 <- Almer(lnAMR ~ lnBWg + tempTestK1000 + (1|species) , data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model0.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model0intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (1|species), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model0int <- Almer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) , data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  
  Phylo_MMR_model1 <- Almer(lnAMR ~ lnBWg + tempTestK1000 + (1|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model1.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model1intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model1int <- Almer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  
  Phylo_MMR_model2 <- Almer(lnAMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model2int <- Almer(lnAMR ~ lnBWg * tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model2.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model2intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTestK1000,2)  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  
  Phylo_MMR_model4 <- Almer(lnAMR ~ lnBWg + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model4int <- Almer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model4.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model4intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  
  Phylo_MMR_model5 <- Almer(lnAMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model5int <- Almer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model5.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model5intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  
  Phylo_MMR_model6 <- Almer(lnAMR ~ lnBWg + tempTestK1000 + (lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model6int <- Almer(lnAMR ~ lnBWg * tempTestK1000+ (lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model6.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model6intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  
  ### BIC --------------
  MMR_BIC<-BICdelta(BIC(Phylo_MMR_model0, Phylo_MMR_model0int, Phylo_MMR_model0.POLY, Phylo_MMR_model0intPOLY,
                  Phylo_MMR_model1, Phylo_MMR_model1int, Phylo_MMR_model1.POLY, Phylo_MMR_model1intPOLY,
                  Phylo_MMR_model2, Phylo_MMR_model2int, Phylo_MMR_model2.POLY, Phylo_MMR_model2intPOLY,
                  Phylo_MMR_model4, Phylo_MMR_model4int, Phylo_MMR_model4.POLY, Phylo_MMR_model4intPOLY,
                  Phylo_MMR_model5, Phylo_MMR_model5int, Phylo_MMR_model5.POLY, Phylo_MMR_model5intPOLY,
                  Phylo_MMR_model6, Phylo_MMR_model6int, Phylo_MMR_model6.POLY, Phylo_MMR_model6intPOLY))
  
  ## AAS ecologically relevant ------------------
  Phylo_AS_model0 <- Almer(lnAS ~ lnBWg + tempTestK1000 + (1|species) , data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model0.POLY <- Almer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model0intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTestK1000,2) + (1|species), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model0int <- Almer(lnAS ~ lnBWg * tempTestK1000 + (1|species) , data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  
  Phylo_AS_model1 <- Almer(lnAS ~ lnBWg + tempTestK1000 + (1|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model1.POLY <- Almer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model1intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model1int <- Almer(lnAS ~ lnBWg * tempTestK1000 + (1|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  
  Phylo_AS_model2 <- Almer(lnAS ~ lnBWg + tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model2int <- Almer(lnAS ~ lnBWg * tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model2.POLY <- Almer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model2intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTestK1000,2)  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  
  Phylo_AS_model4 <- Almer(lnAS ~ lnBWg + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model4int <- Almer(lnAS ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model4.POLY <- Almer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model4intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  
  Phylo_AS_model5 <- Almer(lnAS ~ lnBWg + tempTestK1000  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model5int <- Almer(lnAS ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model5.POLY <- Almer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model5intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  
  Phylo_AS_model6 <- Almer(lnAS ~ lnBWg + tempTestK1000 + (lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model6int <- Almer(lnAS ~ lnBWg * tempTestK1000+ (lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model6.POLY <- Almer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model6intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  
  ## BIC --------
  AS_BIC<-BICdelta(BIC(Phylo_AS_model0, Phylo_AS_model0int, Phylo_AS_model0.POLY, Phylo_AS_model0intPOLY,
                  Phylo_AS_model1, Phylo_AS_model1int, Phylo_AS_model1.POLY, Phylo_AS_model1intPOLY,
                  Phylo_AS_model2, Phylo_AS_model2int, Phylo_AS_model2.POLY, Phylo_AS_model2intPOLY,
                  Phylo_AS_model4, Phylo_AS_model4int, Phylo_AS_model4.POLY, Phylo_AS_model4intPOLY,
                  Phylo_AS_model5, Phylo_AS_model5int, Phylo_AS_model5.POLY, Phylo_AS_model5intPOLY,
                  Phylo_AS_model6, Phylo_AS_model6int, Phylo_AS_model6.POLY, Phylo_AS_model6intPOLY))
  
  
  ## FAS ecologically relevant ----------------
  Phylo_FAS_model0 <- Almer(log(FAS) ~ lnBWg + tempTest + (1|species) , data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model0.POLY <- Almer(log(FAS) ~ lnBWg + poly(tempTest,2) + (1|species), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model0intPOLY <- Almer(log(FAS) ~ lnBWg * poly(tempTest,2) + (1|species), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model0int <- Almer(log(FAS) ~ lnBWg * tempTest + (1|species) , data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  
  Phylo_FAS_model1 <- Almer(log(FAS) ~ lnBWg + tempTest + (1|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model1.POLY <- Almer(log(FAS) ~ lnBWg + poly(tempTest,2) + (1|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model1intPOLY <- Almer(log(FAS) ~ lnBWg * poly(tempTest,2) + (1|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model1int <- Almer(log(FAS) ~ lnBWg * tempTest + (1|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  
  Phylo_FAS_model2 <- Almer(log(FAS) ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model2int <- Almer(log(FAS) ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model2.POLY <- Almer(log(FAS) ~ lnBWg + poly(tempTest,2) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model2intPOLY <- Almer(log(FAS) ~ lnBWg * poly(tempTest,2)  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  
  Phylo_FAS_model4 <- Almer(log(FAS) ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model4int <- Almer(log(FAS) ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model4.POLY <- Almer(log(FAS) ~ lnBWg + poly(tempTest,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model4intPOLY <- Almer(log(FAS) ~ lnBWg * poly(tempTest,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  
  Phylo_FAS_model5 <- Almer(log(FAS) ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model5int <- Almer(log(FAS) ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model5.POLY <- Almer(log(FAS) ~ lnBWg + poly(tempTest,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model5intPOLY <- Almer(log(FAS) ~ lnBWg * poly(tempTest,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  
  Phylo_FAS_model6 <- Almer(lnFAS ~ lnBWg + tempTestK1000 + (lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model6int <- Almer(lnFAS ~ lnBWg * tempTestK1000+ (lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model6.POLY <- Almer(lnFAS ~ lnBWg + poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model6intPOLY <- Almer(lnFAS ~ lnBWg * poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  
  ## BIC -------------
  FAS_BIC<-BICdelta(BIC(Phylo_FAS_model0, Phylo_FAS_model0int, Phylo_FAS_model0.POLY, Phylo_FAS_model0intPOLY,
                  Phylo_FAS_model1, Phylo_FAS_model1int, Phylo_FAS_model1.POLY, Phylo_FAS_model1intPOLY,
                  Phylo_FAS_model2, Phylo_FAS_model2int, Phylo_FAS_model2.POLY, Phylo_FAS_model2intPOLY,
                  Phylo_FAS_model4, Phylo_FAS_model4int, Phylo_FAS_model4.POLY, Phylo_FAS_model4intPOLY,
                  Phylo_FAS_model5, Phylo_FAS_model5int, Phylo_FAS_model5.POLY, Phylo_FAS_model5intPOLY,
                  Phylo_FAS_model6, Phylo_FAS_model6int, Phylo_FAS_model6.POLY, Phylo_FAS_model6intPOLY))
  
  
  # RMR warm ------------------
  Phylo_RMR_W_model0 <- Almer(lnRMR ~ lnBWg + tempTestK1000 + (1|species) , data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model0.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model0intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (1|species), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model0int <- Almer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) , data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  
  Phylo_RMR_W_model1 <- Almer(lnRMR ~ lnBWg + tempTestK1000 + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model1.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model1intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model1int <- Almer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  
  Phylo_RMR_W_model2 <- Almer(lnRMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model2int <- Almer(lnRMR ~ lnBWg * tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model2.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model2intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTestK1000,2)  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  
  Phylo_RMR_W_model4 <- Almer(lnRMR ~ lnBWg + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model4int <- Almer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model4.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model4intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  
  Phylo_RMR_W_model5 <- Almer(lnRMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model5int <- Almer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model5.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model5intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  
  Phylo_RMR_W_model6 <- Almer(lnRMR ~ lnBWg + tempTestK1000 + (lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model6int <- Almer(lnRMR ~ lnBWg * tempTestK1000+ (lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model6.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model6intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  
  
  ## BIC -------------
  RMR_W_BIC<-BICdelta(BIC(Phylo_RMR_W_model0, Phylo_RMR_W_model0int, Phylo_RMR_W_model0.POLY, Phylo_RMR_W_model0intPOLY,
                  Phylo_RMR_W_model1, Phylo_RMR_W_model1int, Phylo_RMR_W_model1.POLY, Phylo_RMR_W_model1intPOLY,
                  Phylo_RMR_W_model2, Phylo_RMR_W_model2int, Phylo_RMR_W_model2.POLY, Phylo_RMR_W_model2intPOLY,
                  Phylo_RMR_W_model4, Phylo_RMR_W_model4int, Phylo_RMR_W_model4.POLY, Phylo_RMR_W_model4intPOLY,
                  Phylo_RMR_W_model5, Phylo_RMR_W_model5int, Phylo_RMR_W_model5.POLY, Phylo_RMR_W_model5intPOLY,
                  Phylo_RMR_W_model6, Phylo_RMR_W_model6int, Phylo_RMR_W_model6.POLY, Phylo_RMR_W_model6intPOLY))
  
  # AMR / warm temps --------------
  Phylo_MMR_W_model0 <- Almer(lnAMR ~ lnBWg + tempTestK1000 + (1|species) , data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model0.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model0intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (1|species), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model0int <- Almer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) , data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  
  Phylo_MMR_W_model1 <- Almer(lnAMR ~ lnBWg + tempTestK1000 + (1|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model1.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model1intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model1int <- Almer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  
  Phylo_MMR_W_model2 <- Almer(lnAMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model2int <- Almer(lnAMR ~ lnBWg * tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model2.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model2intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTestK1000,2)  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  
  Phylo_MMR_W_model4 <- Almer(lnAMR ~ lnBWg + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model4int <- Almer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model4.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model4intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  
  Phylo_MMR_W_model5 <- Almer(lnAMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model5int <- Almer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model5.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model5intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  
  Phylo_MMR_W_model6 <- Almer(lnAMR ~ lnBWg + tempTestK1000 + (lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model6int <- Almer(lnAMR ~ lnBWg * tempTestK1000+ (lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model6.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model6intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  
  ## BIC -------------
  MMR_W_BIC<-BICdelta(BIC(Phylo_MMR_W_model0, Phylo_MMR_W_model0int, Phylo_MMR_W_model0.POLY, Phylo_MMR_W_model0intPOLY,
                  Phylo_MMR_W_model1, Phylo_MMR_W_model1int, Phylo_MMR_W_model1.POLY, Phylo_MMR_W_model1intPOLY,
                  Phylo_MMR_W_model2, Phylo_MMR_W_model2int, Phylo_MMR_W_model2.POLY, Phylo_MMR_W_model2intPOLY,
                  Phylo_MMR_W_model4, Phylo_MMR_W_model4int, Phylo_MMR_W_model4.POLY, Phylo_MMR_W_model4intPOLY,
                  Phylo_MMR_W_model5, Phylo_MMR_W_model5int, Phylo_MMR_W_model5.POLY, Phylo_MMR_W_model5intPOLY,
                  Phylo_MMR_W_model6, Phylo_MMR_W_model6int, Phylo_MMR_W_model6.POLY, Phylo_MMR_W_model6intPOLY))
  
  
  
  # AS / warm temps --------------
  Phylo_AS_W_model0 <- Almer(lnAS ~ lnBWg + tempTestK1000 + (1|species) , data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model0.POLY <- Almer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model0intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTestK1000,2) + (1|species), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model0int <- Almer(lnAS ~ lnBWg * tempTestK1000 + (1|species) , data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  
  Phylo_AS_W_model1 <- Almer(lnAS ~ lnBWg + tempTestK1000 + (1|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model1.POLY <- Almer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model1intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model1int <- Almer(lnAS ~ lnBWg * tempTestK1000 + (1|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  
  Phylo_AS_W_model2 <- Almer(lnAS ~ lnBWg + tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model2int <- Almer(lnAS ~ lnBWg * tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model2.POLY <- Almer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model2intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTestK1000,2)  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  
  Phylo_AS_W_model4 <- Almer(lnAS ~ lnBWg + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model4int <- Almer(lnAS ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model4.POLY <- Almer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model4intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  
  Phylo_AS_W_model5 <- Almer(lnAS ~ lnBWg + tempTestK1000  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model5int <- Almer(lnAS ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model5.POLY <- Almer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model5intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  
  Phylo_AS_W_model6 <- Almer(lnAS ~ lnBWg + tempTestK1000 + (lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model6int <- Almer(lnAS ~ lnBWg * tempTestK1000+ (lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model6.POLY <- Almer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model6intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  
  ## BIC -------------
  AS_W_BIC<-BICdelta(BIC(Phylo_AS_W_model0, Phylo_AS_W_model0int, Phylo_AS_W_model0.POLY, Phylo_AS_W_model0intPOLY,
                  Phylo_AS_W_model1, Phylo_AS_W_model1int, Phylo_AS_W_model1.POLY, Phylo_AS_W_model1intPOLY,
                  Phylo_AS_W_model2, Phylo_AS_W_model2int, Phylo_AS_W_model2.POLY, Phylo_AS_W_model2intPOLY,
                  Phylo_AS_W_model4, Phylo_AS_W_model4int, Phylo_AS_W_model4.POLY, Phylo_AS_W_model4intPOLY,
                  Phylo_AS_W_model5, Phylo_AS_W_model5int, Phylo_AS_W_model5.POLY, Phylo_AS_W_model5intPOLY,
                  Phylo_AS_W_model6, Phylo_AS_W_model6int, Phylo_AS_W_model6.POLY, Phylo_AS_W_model6intPOLY))
  
  
  # FAS / warm temps --------------
  Phylo_FAS_W_model0 <- Almer(lnFAS ~ lnBWg + tempTest + (1|species) , data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model0.POLY <- Almer(lnFAS ~ lnBWg + poly(tempTest,2) + (1|species), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model0intPOLY <- Almer(lnFAS ~ lnBWg * poly(tempTest,2) + (1|species), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model0int <- Almer(lnFAS ~ lnBWg * tempTest + (1|species) , data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  
  Phylo_FAS_W_model1 <- Almer(lnFAS ~ lnBWg + tempTest + (1|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model1.POLY <- Almer(lnFAS ~ lnBWg + poly(tempTest,2) + (1|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model1intPOLY <- Almer(lnFAS ~ lnBWg * poly(tempTest,2) + (1|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model1int <- Almer(lnFAS ~ lnBWg * tempTest + (1|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  
  Phylo_FAS_W_model2 <- Almer(lnFAS ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model2int <- Almer(lnFAS ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model2.POLY <- Almer(lnFAS ~ lnBWg + poly(tempTest,2) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model2intPOLY <- Almer(lnFAS ~ lnBWg * poly(tempTest,2)  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  
  Phylo_FAS_W_model4 <- Almer(lnFAS ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model4int <- Almer(lnFAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model4.POLY <- Almer(lnFAS ~ lnBWg + poly(tempTest,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model4intPOLY <- Almer(lnFAS ~ lnBWg * poly(tempTest,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  
  Phylo_FAS_W_model5 <- Almer(lnFAS ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model5int <- Almer(lnFAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model5.POLY <- Almer(lnFAS ~ lnBWg + poly(tempTest,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model5intPOLY <- Almer(lnFAS ~ lnBWg * poly(tempTest,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  
  Phylo_FAS_W_model6 <- Almer(lnFAS ~ lnBWg + tempTestK1000 + (lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model6int <- Almer(lnFAS ~ lnBWg * tempTestK1000+ (lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model6.POLY <- Almer(lnFAS ~ lnBWg + poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model6intPOLY <- Almer(lnFAS ~ lnBWg * poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  
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
  AS_W_BIC  # Phylo_AS_W_model4
  FAS_W_BIC # Phylo_FAS_W_model2.POLY
  
  # Phylo_AS_W_model1intPOLY  9 737.4351  0.00000
  # Phylo_AS_W_model1int      7 738.2306  0.79549
  # Phylo_AS_W_model4         7 738.3042  0.86917 # << use this
  
  # Best models  --------
  rmr_mod_ER<-Phylo_RMR_model5
  amr_mod_ER<-Phylo_MMR_model4int
  as_mod_ER<-Phylo_AS_model4
  fas_mod_ER<-Phylo_FAS_model4
  
  rmr_mod_W<-Phylo_RMR_W_model1
  amr_mod_W<-Phylo_MMR_W_model4
  as_mod_W<-Phylo_AS_W_model4
  fas_mod_W<-Phylo_FAS_W_model2.POLY
  
}else{
  # running best models only 
  rmr_mod_ER <- Almer(lnRMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  amr_mod_ER <- Almer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  as_mod_ER <- Almer(lnAS ~ lnBWg + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  fas_mod_ER <- Almer(log(FAS) ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  
  rmr_mod_W <- Almer(lnRMR ~ lnBWg + tempTestK1000 + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  amr_mod_W <- Almer(lnAMR ~ lnBWg + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  as_mod_W <- Almer(lnAS ~ lnBWg + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  fas_mod_W <- Almer(lnFAS ~ lnBWg + poly(tempTest,2) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
}

# un-comment to see - or run R markdown html summary script

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

data.as.test$resid<-resid(as_mod_W)
data.as.test[which(data.as.test$resid < -1.5),] # could be considered outlier data
data.fasER[which(data.fasER$FAS > 20),] # considered outliers, measurement extremes for zebrafish

# *****************************************************************************
# *****************************************************************************

# Model scaling parameters and CIs ----------
# custom function to obtain model parameters and recalculated mass specific MR using model estimate scaling slopes
# this also expands the dataset to get mass independent values of metabolic rates
model_outputs(phylo = TRUE, 
              best.model.rmr.er = rmr_mod_ER,
              best.model.amr.er= amr_mod_ER,
              best.model.as.er = as_mod_ER,
              best.model.fas.er = fas_mod_ER,
              best.model.rmr.w = rmr_mod_W,
              best.model.amr.w = amr_mod_W,
              best.model.as.w = as_mod_W,
              best.model.fas.w = fas_mod_W)

# *****************************************************************************
# *****************************************************************************



# ******************************************************************************
# Figures -------
cols.as<-c("#265F73", "#007E66", "#00C5A3")
cols.fas<-c("#395200", "#89A000", "yellow")
cols.rmr<-c("#C70039", "#FF6D7C", "#FFA3AC")
cols.amr<-c("#00749F","#00A8D6", "#9CE9FF")
cols<-c("#00749F","#C70039","#00A8D6","#FF6D7C", "#9CE9FF","#FFA3AC", "#00C5A3", "#265F73")# AMR -rmr- AMR dark - rmr dark - light - as - fas

set.seed(51423)
# MMR and AMR used interchangeably throughout 

# General scaling plots ------
AMRmodel_plot1<-ggplot(data=data.amrER, aes(x=lnBWg, y=lnAMR)) +
  geom_point(alpha=0.9,  size=2, pch=19, color="grey70")+
  geom_line(data=data.plotAMRint_ER,
            aes(y = model_predFE, x=lnBWg,  group=tempTestK1000_inC, color = tempTestK1000_inC),
            linewidth=0.3, lty=1,alpha=0.8, show.legend=FALSE) +
  geom_point(alpha=0.9, size=2, pch=21, color="grey50",fill="grey70" )+
  geom_point(data=data.amr.test, aes(x=lnBWg, y=lnAMR, fill=tempTest),
             alpha=0.9,  size=2, pch=21, show.legend = FALSE)+
  geom_line(data=data.plotAMR_warm[round(data.plotAMR_warm$tempTestK1000,2)==3.39,],
            aes(y = model_predFE, x=lnBWg,  group=tempTestK1000), color="#002E53", linewidth=1, lty=1, show.legend=FALSE) +
  annotate("text",  x = -5.5, y = 11.5, label = bquote(Optimal:~italic(b)[MMR] == change~with~degree*C),size=5, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = -5.5, y = 9.8, label = bquote(Warm:~italic(b)[MMR] == .(MMR_slope_w)),size=5, hjust=0, family="Arial", color = cols.amr[1])+
  annotate("text",  x = 3.0, y = -2, label = bquote(~0*degree*C:~italic(b)[MMR] == .(round(AMR.slopes$lnBWg.trend[1],3))),size=4, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = 3.0, y = -3, label = bquote(10*degree*C:~italic(b)[MMR] == .(round(AMR.slopes$lnBWg.trend[2],3))),size=4, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = 3.0, y = -4, label = bquote(20*degree*C:~italic(b)[MMR] == .(round(AMR.slopes$lnBWg.trend[3],3))),size=4, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = 3.0, y = -5, label = bquote(30*degree*C:~italic(b)[MMR] == .(round(AMR.slopes$lnBWg.trend[4],3))),size=4, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = 3.0, y = -6, label = bquote(40*degree*C:~italic(b)[MMR] == .(round(AMR.slopes$lnBWg.trend[5],3))),size=4, hjust=0, family="Arial", color = "black")+
  scale_fill_gradient( low = cols.amr[3], high = cols.amr[1])+
  scale_color_gradient( low = "grey80", high = "grey0")+
  ylim(x = -6.5, 12)+
  xlim(x = -6.5, 12)+
  annotate("text", label = paste("n = ", nrow(data.amrER), sep=""),   x = -5.5, y = 8.5, size=3, hjust=0, family="Arial", color = "black")+
  annotate("text", label = paste("n = ", nrow(data.amr.test), sep=""),  x = -5.5, y = 7.6, size=3, hjust=0, family="Arial", color = cols.amr[1])
ggformat(AMRmodel_plot1, x_title=expression(italic(ln)*Body~weight~(g)), y_title=expression(italic(ln)*MMR~(mg~O[2]~h^-1)), print = F)

RMRmodel_plot1<-ggplot(data=data.rmrER, aes(x=lnBWg, y=lnRMR)) +
  geom_point(alpha=0.9,  size=2, pch=21, color="grey50",fill="grey70" )+
  geom_line(data=data.plotRMR_ER[round(data.plotRMR_ER$tempTestK1000,2)==3.39,], mapping=aes(y = model_predFE, x=lnBWg,  group=tempTestK1000, color= tempTestK1000), color="black", linewidth=1, lty=1, show.legend=FALSE) +
  geom_point(data=data.rmr.test, aes(x=lnBWg, y=lnRMR, fill=tempTest), alpha=0.9,  size=2, pch=21, show.legend = FALSE)+
  geom_line(data=data.plotRMR_warm[round(data.plotRMR_warm$tempTestK1000,2)==3.39,],
            aes(y = model_predFE, x=lnBWg,  group=tempTestK1000), color="#7C0003", linewidth=1, lty=1, show.legend=FALSE) +
  annotate("text",  x = -5.5, y = 11.5, label = bquote(Optimal:~italic(b)[RMR] == .(RMR_slope)),size=5, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = -5.5, y = 9.8, label = bquote(Warm:~italic(b)[RMR] == .(RMR_slope_w)),size=5, hjust=0, family="Arial", color = cols.rmr[1])+
  scale_fill_gradient( low = cols.rmr[3], high = cols.rmr[1])+
  scale_color_gradient( low = cols.rmr[3], high = cols.rmr[1])+
  ylim(x = -6.5, 12)+
  xlim(x = -6.5, 12)+
  annotate("text", label = paste("n = ", nrow(data.rmrER), sep=""),   x = -5.5, y = 8.5, size=3, hjust=0, family="Arial", color = "black")+
  annotate("text", label = paste("n = ", nrow(data.rmr.test), sep=""),  x = -5.5, y = 7.6, size=3, hjust=0, family="Arial", color = cols.rmr[1])
ggformat(RMRmodel_plot1, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*RMR~(mg~O[2]~h^-1)), print = FALSE)

# FAS! 
FASmodel_plot1<-ggplot(data=data.fasER, aes(x=lnBWg, y=log(FAS))) +
  geom_point(alpha=0.9,  size=2, pch=21, color="grey50",fill="grey70" )+
  geom_line(data=data.plotFAS_ER, aes(y = model_predFE, x=lnBWg,  group=tempTest), color="black", linewidth=1, lty=1, show.legend=FALSE) +
  geom_point(data=data.fas.test, aes(x=lnBWg, y=log(FAS), fill=tempTest), alpha=0.9,  size= 2, pch=21, show.legend=FALSE)+
  geom_line(data=data.plotFAS_warm, aes(y = model_predFE, x=lnBWg,  group=tempTest), color="#475500", linewidth=1, lty=1, show.legend=FALSE) +
  annotate("text",  x = 4, y = 3.9, label = bquote(Optimal:~italic(b)[FAS] == .(FAS_slope)),size=5, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = 4, y = 3.5, label = bquote(Warm:~italic(b)[FAS] == .(FAS_slope_w)),size=5, hjust=0, family="Arial", color = "#475500")+
  annotate("text", label = paste("n = ", nrow(data.fasER), sep=""),   x = 4, y = 3.1, size=3, hjust=0, family="Arial", color = "black")+
  annotate("text", label = paste("n = ", nrow(data.fas.test), sep=""),  x = 4, y = 2.85, size=3, hjust=0, family="Arial", color = "#475500")+
  scale_fill_gradient(high = "yellow", low = "black")+
  ylim(0,4)
ggformat(FASmodel_plot1, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*FAS), print = FALSE)

# AS 
ASmodel_plot1<-ggplot(data=data.asER, aes(x=lnBWg, y=lnAS)) +
  geom_point(alpha=0.9,  size=2, pch=21, color="grey50",fill="grey70" )+
  geom_line(data=data.plotAS_ER[round(data.plotAS_ER$tempTestK1000,2)==3.39,], mapping=aes(y = model_predFE, x=lnBWg,  group=tempTestK1000, color= tempTestK1000), color="black", linewidth=1, lty=1, show.legend=FALSE) +
  geom_point(data=data.as.test, aes(x=lnBWg, y=lnAS, fill=tempTest), alpha=0.9,  size=2, pch=21, show.legend = FALSE)+
  geom_line(data=data.plotAS_warm[round(data.plotAS_warm$tempTestK1000,1)==round(C20inTempTestK1000,1),],
            aes(y = model_predFE, x=lnBWg,  group=tempTestK1000), color="#00785B", linewidth=1, lty=1, show.legend=FALSE) +
  annotate("text",  x = -5.5, y = 11.5, label = bquote(Optimal:~italic(b)[AS] == .(AS_slope)),size=5, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = -5.5, y = 9.8, label = bquote(Warm:~italic(b)[AS] == .(AS_slope_w)),size=5, hjust=0, family="Arial", color = cols.as[1])+
  scale_fill_gradient( low = cols.as[3], high = cols.as[1])+
  scale_color_gradient( low = cols.as[3], high = cols.as[1])+
  ylim(x = -6.5, 12)+
  xlim(x = -6.5, 12)+
  annotate("text", label = paste("n = ", nrow(data.asER), sep=""),   x = -5.5, y = 8.5, size=3, hjust=0, family="Arial", color = "black")+
  annotate("text", label = paste("n = ", nrow(data.as.test), sep=""),  x = -5.5, y = 7.6, size=3, hjust=0, family="Arial", color = cols.as[1])
ggformat(ASmodel_plot1, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*AS~(mg~O[2]~h^-1)), print = F)


scaling<-cowplot:::plot_grid(AMRmodel_plot1, RMRmodel_plot1,
                             ASmodel_plot1, FASmodel_plot1,
                              align = "hv",
                              axis = "l",
                              nrow = 2,
                              ncol = 2,
                              labels = "AUTO",
                              label_size = 17)
ggsave(filename = paste("./Figures/Fig2_Scaling_phylo_.png", sep=""),
       plot=scaling, width = 8.5, height = 8.5, units = "in")

# Both Activation energies together: -------
# model_predFE << is in lnMR units 
# MMR
AMRmodel_plot3.E_ALL<-ggplot(data.amrER[data.amrER$lnBWg==0,]) +
  geom_point(data.amrER, mapping=aes(x=tempTestK1000, y=log(mass_specamr), fill= tempTest, size = lnBWg), color="grey50", fill="grey70", alpha=1, pch=21, show.legend = FALSE)+
  geom_point(data.amr.test, mapping=aes(x=tempTestK1000, y=log(mass_specamr), fill= tempTest, size = lnBWg), color="black",  alpha=1, pch=21,  show.legend = FALSE)+
  scale_fill_gradient( low = cols.amr[3], high = cols.amr[1])+
  scale_color_gradient( low = "grey80", high = "grey0")+
  ylim(-5.7,4)+
  annotate("text",  x = 3.27, y = 3.5, label = bquote(italic(E)[MMR] == change~with~mass),size=5, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = 3.27, y = 2.5, label = bquote(italic(E)[MMR] == .(round(MMR_E_W_eV,3))),size=5, hjust=0, family="Arial", color = cols.amr[1])+
  annotate("text",  x = 3.27, y = -3.2, label = bquote(1*g~italic(E)== .(round(MMR_E_ER$MMR_E_ER_eV[1],3))),size=4, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = 3.27, y = -4, label = bquote(10*g~italic(E) == .(round(MMR_E_ER$MMR_E_ER_eV[2],3))),size=4, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = 3.27, y = -4.8, label = bquote(100*g~italic(E) == .(round(MMR_E_ER$MMR_E_ER_eV[3],3))),size=4, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = 3.27, y = -5.6, label = bquote(1000*g~italic(E) == .(round(MMR_E_ER$MMR_E_ER_eV[4],3))),size=4, hjust=0, family="Arial", color = "black")+
  scale_x_continuous(limits = c(3.25, 3.661), sec.axis = sec_axis(~ ((1000/.))-273.15, name = expression(Temperature~degree*C), breaks = c(32, 25, 17, 10, 3 )))+
  geom_line(data = data.plotAMR_warm[which(round(data.plotAMR_warm$lnBWg, 1) == round(log(1.4), 1)),],
            aes(y = log((exp(model_predFE)/exp(lnBWg))), x=tempTestK1000, group=lnBWg) ,color=cols.amr[2], linewidth=1, lty=1, show.legend=FALSE)+
  geom_line(data = data.plotAMRint_ER[exp(data.plotAMRint_ER$lnBWg) <= 1000, ], aes(y = log((exp(model_predFE)/exp(lnBWg))), x=tempTestK1000, group=lnBWg, color = lnBWg), linewidth=0.5, lty=1, show.legend=FALSE)
ggformat(AMRmodel_plot3.E_ALL, x_title= expression(Temperature^-1~(1000/K)), y_title=expression(italic(ln)*MMR~(mg~O[2]~h^-1~g^-1)), print = F)

RMRmodel_plot3.E_ALL<-ggplot(data.rmrER[data.rmrER$lnBWg==0,]) +
  geom_point(data.rmrER, mapping=aes(x=tempTestK1000, y=log(mass_specrmr), fill= tempTest, size = lnBWg), color="grey50", fill="grey70", alpha=1, pch=21)+
  geom_point(data.rmr.test, mapping=aes(x=tempTestK1000, y=log(mass_specrmr), fill= tempTest, size = lnBWg), color = "black", alpha=1, pch=21, show.legend = FALSE)+
  scale_fill_gradient( low = cols.rmr[3], high = cols.rmr[1])+
  annotate("text",  x = 3.27, y = 3.5, label = bquote(italic(E)[RMR] == .(round(RMR_E_ER_eV,3))),size=5, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = 3.27, y = 2.5, label = bquote(italic(E)[RMR] == .(round(RMR_E_W_eV,3))),size=5, hjust=0, family="Arial", color = cols.rmr[1])+
  ylim(-5.7,4)+
  scale_size_continuous(name="Mass (g)",
                      breaks=c(log(0.1), log(10), log(1000)),
                      labels=c("0.1", "10", "1000"))+
  scale_x_continuous(limits = c(3.25, 3.661), sec.axis = sec_axis(~ ((1000/.))-273.15, name = expression(Temperature~degree*C), breaks = c(32, 25, 17, 10, 3 )))+
  geom_line(data = data.plotRMR_ER[which(round(data.plotRMR_ER$lnBWg, 1) == round(log(1.35), 1)),], aes(y = log((exp(model_predFE)/exp(lnBWg))), x=tempTestK1000, group=lnBWg ) ,color="black", linewidth=1, lty=1, show.legend=FALSE)+
  geom_line(data = data.plotRMR_warm[which(round(data.plotRMR_warm$lnBWg, 1) == round(log(1.35), 1)),], aes(y = log((exp(model_predFE)/exp(lnBWg))), x=tempTestK1000, group=lnBWg), color = cols.rmr[2], linewidth=1, lty=1, show.legend=FALSE)
ggformat(RMRmodel_plot3.E_ALL, x_title= expression(Temperature^-1~(1000/K)), y_title=expression(italic(ln)*RMR~(mg~O[2]~h^-1~g^-1)), print = FALSE)
RMRmodel_plot3.E_ALL <- RMRmodel_plot3.E_ALL + theme(legend.position = c(0.87, 0.82))

ASmodel_plot3.E_ALL<-ggplot(data.asER[data.asER$lnBWg==0,]) +
  geom_point(data.asER, mapping=aes(x=tempTestK1000, y=log(mass_specas), fill= tempTest, size = lnBWg), color="grey50", fill="grey70", alpha=1, pch=21, show.legend = FALSE)+
  geom_point(data.as[!c(data.as$test_category=="ecol_relev"),], mapping=aes(x=tempTestK1000, y=log(mass_specas), fill= tempTest, size = lnBWg), color = "black", alpha=1, pch=21, show.legend = FALSE)+
  scale_fill_gradient( low = cols.as[3], high = cols.as[1])+
  annotate("text",  x = 3.27, y = 3.5, label = bquote(italic(E)[AS] == .(round(AS_E_ER_eV,3))),size=5, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = 3.27, y = 2.5, label = bquote(italic(E)[AS] == .(round(AS_E_W_eV,3))),size=5, hjust=0, family="Arial", color = cols.as[1])+
  ylim(-5.7,4)+
  scale_x_continuous(limits = c(3.25, 3.661), sec.axis = sec_axis(~ ((1000/.))-273.15, name = expression(Temperature~degree*C), breaks = c(32, 25, 17, 10, 3 )))+
  geom_line(data = data.plotAS_ER[which(round(data.plotAS_ER$lnBWg, 1) == round(log(1.35), 1)),], aes(y = log((exp(model_predFE)/exp(lnBWg))), x=tempTestK1000, group=lnBWg ) ,color="black", linewidth=1, lty=1, show.legend=FALSE)+
  geom_line(data = data.plotAS_warm[which(round(data.plotAS_warm$lnBWg, 1) == round(log(1.4), 1)),], aes(y = log((exp(model_predFE)/exp(lnBWg))), x=tempTestK1000, group=lnBWg), color = cols.as[2], linewidth=1, lty=1, show.legend=FALSE)
ggformat(ASmodel_plot3.E_ALL, x_title= expression(Temperature^-1~(1000/K)), y_title=expression(italic(ln)*AS~(mg~O[2]~h^-1~g^-1)), print = F)

Arh.plot<-cowplot:::plot_grid(AMRmodel_plot3.E_ALL, RMRmodel_plot3.E_ALL, ASmodel_plot3.E_ALL, 
          align = "hv",
          axis = "l",
          nrow = 1,
          ncol = 3,
          labels = "AUTO",
          label_y = 0.95,
          label_size = 17)
ggsave(filename = paste("./Figures/Fig3_ArrheniusFigMMR-RMR-AS_Phylo.png", sep=""),
       plot=Arh.plot, width = 13, height = 5, units = "in")


# Both MMR and RMR together - AS punchline plots ---------
# All fish together:
MRmodel_plot2<-ggplot(data=data.rmrER, aes(x=lnBWg, y=lnRMR)) +
  geom_line(data=data.plotRMR_ER[round(data.plotRMR_ER$tempTestK1000,2)==3.39,], aes(y = model_predFE, x=lnBWg,  group=tempTestK1000, color= tempTestK1000), color="black", size=0.7, lty=1, show.legend=FALSE) +
  geom_line(data=data.plotAMRint_ER[round(data.plotAMRint_ER$tempTestK1000,2)==3.39,], aes(y = model_predFE, x=lnBWg,  group=tempTestK1000, color= tempTestK1000), color="black", size=0.7, lty=1, show.legend=FALSE) +
  geom_line(data=data.plotRMR_warm[round(data.plotRMR_warm$tempTestK1000,2)==3.39,], aes(y = model_predFE, x=lnBWg,  group=tempTestK1000), color=cols.rmr[2], size=1, lty=1, show.legend=FALSE) +
  geom_line(data=data.plotAMR_warm[round(data.plotAMR_warm$tempTestK1000,2)==3.39,], aes(y = model_predFE, x=lnBWg,  group=tempTestK1000), color=cols.amr[2], size=1, lty=1, show.legend=FALSE) +
  ylim(x = -6.5, 12)+
  xlim(x = -6.5, 12)+
  annotate("text",  x = -0, y = 11.5, label = bquote(italic(b)[MMR] ~"("* .(MMR_slope_w) *")"),size=5, hjust=0, family="Arial", color = cols.amr[1])+
  annotate("text",  x = -6.5, y = 11.5, label = bquote(italic(b)[RMR] ~"("* .(RMR_slope_w) *")"),size=5, hjust=0, family="Arial", color = cols.rmr[1])+
  annotate("text",  x = -0.9, y = 11.5, label = expression(paste("">"")),size=5, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = -0, y = 9.8 , label = bquote(italic(b)[MMR~20*degree*C] ~"("* .(round(AMR.slopes$lnBWg.trend[3],3)) *")"),size=5, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = -6.5, y = 9.8, label = bquote(italic(b)[RMR] ~"("* .(RMR_slope) *")"),size=5, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = -0.9, y = 9.8, label = expression(paste("" %~~% "")),size=5, hjust=0, family="Arial", color = "black")+
  scale_fill_gradient( low = cols.amr[3], high = cols.amr[1])+
  scale_color_gradient( low = "grey80", high = "grey0")+
  annotate("segment", x = 5.2, xend = 5.2, y = 6, yend = 4.5,
           colour = cols.amr[1], size = 1, arrow = arrow(length = unit(0.3,"cm"), type="closed",))+
  annotate("segment", x = -4, xend = -4, yend = -4, y = -5.2,
           colour = cols.rmr[1], size = 1, arrow = arrow(length = unit(0.3,"cm"), type="closed",))
ggformat(MRmodel_plot2, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*MR~(mg~O[2]~h^-1)), print = T)
MRmodel_plot2<<-MRmodel_plot2

# ggsave(filename = paste("./Figures/Fig5-final_MMR-RMR_Phylo.png", sep=""),
#        plot=MRmodel_plot2, width = 4, height = 4, units = "in")


# Boxplots: FAS, AS, MR, mass-independent models ----------
fas_boxplot<-ggplot(data=data.fas, aes(y=FAS, x = test_category,  fill=test_category3))+
  geom_boxplot(show.legend = F)+
  scale_fill_manual(values=cols.fas)+
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

cowplot:::plot_grid(amr_boxplot,rmr_boxplot,fas_boxplot,
          align = "hv",
          axis = "l",
          nrow = 1,
          labels = "AUTO",
          ncol = 3) %>%
ggsave(filename = "./Figures/Supl_fig1_boxplots_Phylo.png", width = 12.5, height =4)

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
# 
#   