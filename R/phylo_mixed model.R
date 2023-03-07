

library(lme4)
library(evolvability) # Almer function
library(ape)
library(rotl)
library(Matrix)
library(emmeans)
library(pryr)
library(ggtree)
library(TDbook)
library(ggimage)
library(phylobase)
library(ggformatKK) # from github

setwd("/Users/kristakraskura/Github_repositories/KK_etal_synthesis_fish_temp_scaling/")
## Source the data ------------
source("./R/get_scaling_data_temp.R")
source("./R/get_phylo_matrix.R")

# function to order model selection based on the lowest BIC score
BICdelta<-function(BICtable){
  BIC.t <- BICtable [order(BICtable$BIC), ]
  BIC.t$delta <- round(abs(BIC.t$BIC[1] -  BIC.t$BIC), 5)
  return( BIC.t)
}

# 1. DATA FOR MODELS ------
data.list<-get_scaling_data_temp(data.amr = "/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Data/MR-fish-metadata-data/Fish_AMR_temp_dataset_mar2022.csv",
                                 data.rmr = "/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Data/MR-fish-metadata-data/Fish_RMR_temp_dataset_mar2022.csv",
                                 ecology.data = "/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Data/MR-fish-metadata-data/Kraskura_species_ecologies_mar2022.csv",
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

get_phylo_matrix(species.list = unique(levels(data.rmrER$species)), matrix.name = "A", tree.name = "tree.kk", dataset.ID = "RMR optimal")
get_phylo_matrix(species.list = unique(levels(data.rmr.test$species)), matrix.name = "A.rmr.w", tree.name = "tr.rmr.w", dataset.ID = "RMR warm")

get_phylo_matrix(species.list = unique(levels(data.amrER$species)), matrix.name = "A.mmr.er", tree.name = "tr.mmr.er", dataset.ID = "MMR optimal")
get_phylo_matrix(species.list = unique(levels(data.amr.test$species)), matrix.name = "A.mmr.w", tree.name = "tr.mmr.w", dataset.ID = "MMR warm")

get_phylo_matrix(species.list = unique(levels(data.asER$species)), matrix.name = "A.aas.er", tree.name = "tr.aas.er", dataset.ID = "AAS optimal")
get_phylo_matrix(species.list = unique(levels(data.as.test$species)), matrix.name = "A.aas.w", tree.name = "tr.aas.w", dataset.ID = "AAS warm")

get_phylo_matrix(species.list = unique(levels(data.fasER$species)), matrix.name = "A.fas.er", tree.name = "tr.fas.er", dataset.ID = "FAS optimal")
get_phylo_matrix(species.list = unique(levels(data.fas.test$species)), matrix.name = "A.fas.w", tree.name = "tree.fas.w", dataset.ID = "FAS warm")

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

# tests herehere
Phylo_RMR_model6 <- Almer(lnRMR ~ lnBWg + tempTestK1000 + (lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
Phylo_RMR_model6int <- Almer(lnRMR ~ lnBWg * tempTestK1000+ (lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
Phylo_RMR_model6.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
Phylo_RMR_model6intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))

### BIC --------------
RMR_BIC<-BICdelta(BIC(Phylo_RMR_model0, Phylo_RMR_model0int, Phylo_RMR_model0.POLY, Phylo_RMR_model0intPOLY,
                Phylo_RMR_model1, Phylo_RMR_model1int, Phylo_RMR_model1.POLY, Phylo_RMR_model1intPOLY,
                Phylo_RMR_model2, Phylo_RMR_model2int, Phylo_RMR_model2.POLY, Phylo_RMR_model2intPOLY,
                Phylo_RMR_model4, Phylo_RMR_model4int, Phylo_RMR_model4.POLY, Phylo_RMR_model4intPOLY,
                Phylo_RMR_model5, Phylo_RMR_model5int, Phylo_RMR_model5.POLY, Phylo_RMR_model5intPOLY))


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

### BIC --------------
MMR_BIC<-BICdelta(BIC(Phylo_MMR_model0, Phylo_MMR_model0int, Phylo_MMR_model0.POLY, Phylo_MMR_model0intPOLY,
                Phylo_MMR_model1, Phylo_MMR_model1int, Phylo_MMR_model1.POLY, Phylo_MMR_model1intPOLY,
                Phylo_MMR_model2, Phylo_MMR_model2int, Phylo_MMR_model2.POLY, Phylo_MMR_model2intPOLY,
                Phylo_MMR_model4, Phylo_MMR_model4int, Phylo_MMR_model4.POLY, Phylo_MMR_model4intPOLY,
                Phylo_MMR_model5, Phylo_MMR_model5int, Phylo_MMR_model5.POLY, Phylo_MMR_model5intPOLY))

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

## BIC --------
AS_BIC<-BICdelta(BIC(Phylo_AS_model0, Phylo_AS_model0int, Phylo_AS_model0.POLY, Phylo_AS_model0intPOLY,
                Phylo_AS_model1, Phylo_AS_model1int, Phylo_AS_model1.POLY, Phylo_AS_model1intPOLY,
                Phylo_AS_model2, Phylo_AS_model2int, Phylo_AS_model2.POLY, Phylo_AS_model2intPOLY,
                Phylo_AS_model4, Phylo_AS_model4int, Phylo_AS_model4.POLY, Phylo_AS_model4intPOLY,
                Phylo_AS_model5, Phylo_AS_model5int, Phylo_AS_model5.POLY, Phylo_AS_model5intPOLY))


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

## BIC -------------
FAS_BIC<-BICdelta(BIC(Phylo_FAS_model0, Phylo_FAS_model0int, Phylo_FAS_model0.POLY, Phylo_FAS_model0intPOLY,
                Phylo_FAS_model1, Phylo_FAS_model1int, Phylo_FAS_model1.POLY, Phylo_FAS_model1intPOLY,
                Phylo_FAS_model2, Phylo_FAS_model2int, Phylo_FAS_model2.POLY, Phylo_FAS_model2intPOLY,
                Phylo_FAS_model4, Phylo_FAS_model4int, Phylo_FAS_model4.POLY, Phylo_FAS_model4intPOLY,
                Phylo_FAS_model5, Phylo_FAS_model5int, Phylo_FAS_model5.POLY, Phylo_FAS_model5intPOLY))


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

## BIC -------------
RMR_W_BIC<-BICdelta(BIC(Phylo_RMR_W_model0, Phylo_RMR_W_model0int, Phylo_RMR_W_model0.POLY, Phylo_RMR_W_model0intPOLY,
                Phylo_RMR_W_model1, Phylo_RMR_W_model1int, Phylo_RMR_W_model1.POLY, Phylo_RMR_W_model1intPOLY,
                Phylo_RMR_W_model2, Phylo_RMR_W_model2int, Phylo_RMR_W_model2.POLY, Phylo_RMR_W_model2intPOLY,
                Phylo_RMR_W_model4, Phylo_RMR_W_model4int, Phylo_RMR_W_model4.POLY, Phylo_RMR_W_model4intPOLY,
                Phylo_RMR_W_model5, Phylo_RMR_W_model5int, Phylo_RMR_W_model5.POLY, Phylo_RMR_W_model5intPOLY))

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

## BIC -------------
MMR_W_BIC<-BICdelta(BIC(Phylo_MMR_W_model0, Phylo_MMR_W_model0int, Phylo_MMR_W_model0.POLY, Phylo_MMR_W_model0intPOLY,
                Phylo_MMR_W_model1, Phylo_MMR_W_model1int, Phylo_MMR_W_model1.POLY, Phylo_MMR_W_model1intPOLY,
                Phylo_MMR_W_model2, Phylo_MMR_W_model2int, Phylo_MMR_W_model2.POLY, Phylo_MMR_W_model2intPOLY,
                Phylo_MMR_W_model4, Phylo_MMR_W_model4int, Phylo_MMR_W_model4.POLY, Phylo_MMR_W_model4intPOLY,
                Phylo_MMR_W_model5, Phylo_MMR_W_model5int, Phylo_MMR_W_model5.POLY, Phylo_MMR_W_model5intPOLY))


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

## BIC -------------
AS_W_BIC<-BICdelta(BIC(Phylo_AS_W_model0, Phylo_AS_W_model0int, Phylo_AS_W_model0.POLY, Phylo_AS_W_model0intPOLY,
                Phylo_AS_W_model1, Phylo_AS_W_model1int, Phylo_AS_W_model1.POLY, Phylo_AS_W_model1intPOLY,
                Phylo_AS_W_model2, Phylo_AS_W_model2int, Phylo_AS_W_model2.POLY, Phylo_AS_W_model2intPOLY,
                Phylo_AS_W_model4, Phylo_AS_W_model4int, Phylo_AS_W_model4.POLY, Phylo_AS_W_model4intPOLY,
                Phylo_AS_W_model5, Phylo_AS_W_model5int, Phylo_AS_W_model5.POLY, Phylo_AS_W_model5intPOLY))


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

## BIC -------------
FAS_W_BIC<-BICdelta(BIC(Phylo_FAS_W_model0, Phylo_FAS_W_model0int, Phylo_FAS_W_model0.POLY, Phylo_FAS_W_model0intPOLY,
                Phylo_FAS_W_model1, Phylo_FAS_W_model1int, Phylo_FAS_W_model1.POLY, Phylo_FAS_W_model1intPOLY,
                Phylo_FAS_W_model2, Phylo_FAS_W_model2int, Phylo_FAS_W_model2.POLY, Phylo_FAS_W_model2intPOLY,
                Phylo_FAS_W_model4, Phylo_FAS_W_model4int, Phylo_FAS_W_model4.POLY, Phylo_FAS_W_model4intPOLY,
                Phylo_FAS_W_model5, Phylo_FAS_W_model5int, Phylo_FAS_W_model5.POLY, Phylo_FAS_W_model5intPOLY))

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

# Ecol relev
summary(rmr_mod_ER)
summary(amr_mod_ER)
summary(fas_mod_ER)
summary(as_mod_ER)

# warm
summary(rmr_mod_W)
summary(amr_mod_W)
summary(fas_mod_W)
summary(as_mod_W)

# Ecol relev
plot(rmr_mod_ER)
plot(amr_mod_ER)
plot(fas_mod_ER)
plot(as_mod_ER)

# warm
plot(rmr_mod_W)
plot(amr_mod_W)
plot(fas_mod_W)
plot(as_mod_W)

hist(resid(amr_mod_ER))
hist(resid(rmr_mod_ER))
hist(resid(fas_mod_ER))
hist(resid(as_mod_ER))

hist(resid(amr_mod_W))
hist(resid(rmr_mod_W))
hist(resid(fas_mod_W))
hist(resid(as_mod_W)) # little skew not too bad, only 5 measurements, all reasonable biologically

data.as.test$resid<-resid(as_mod_W)
data.as.test[which(data.as.test$resid < -1.5),] # could be considered outlier data

# Model scaling parameters and CIs ----------

if (getModelCIs){

  # RMR optimal - ecologically relevant grid; plot/predict full range (All data coverage for reference)
  data.plotRMR_ER<-expand.grid(tempTestK1000 = seq(C40inTempTestK1000, C0inTempTestK1000 , 0.05), lnBWg = seq(sum_rmr_ALL$min_bwLOG,sum_rmr_ALL$max_bwLOG, 0.5))
  data.plotRMR_ER$model_predFE <- predict(rmr_mod_ER, re.form = NA, newdata = data.plotRMR_ER)
  conf.int_RMR_ER <- confint(rmr_mod_ER, method = "boot", FUN=function(x)predict(x, data.plotRMR_ER, re.form=NA), nsim=10)
  data.plotRMR_ER$CI_2.5<-conf.int_RMR_ER[,1]
  data.plotRMR_ER$CI_97.5<-conf.int_RMR_ER[,2]

  # FAS optimal - ecologically relevant grid; plot/predict full range (All data coverage for reference)
  data.plotFAS_ER<-expand.grid(tempTest = 20, lnBWg = seq(sum_fas_ER$min_bwLOG,sum_fas_ER$max_bwLOG, 0.5))
  data.plotFAS_ER$model_predFE <- predict(fas_mod_ER, re.form = NA, newdata = data.plotFAS_ER)
  conf.int_FAS_ER <- confint(fas_mod_ER, method = "boot", FUN=function(x)predict(x, data.plotFAS_ER, re.form=NA), nsim=10)
  data.plotFAS_ER$CI_2.5<-conf.int_FAS_ER[,1]
  data.plotFAS_ER$CI_97.5<-conf.int_FAS_ER[,2]

  # AS optimal - ecologically relevant grid; plot/predict full range (All data coverage for reference)
  data.plotAS_ER<-expand.grid(tempTestK1000 = seq(C40inTempTestK1000, C0inTempTestK1000 , 0.05), lnBWg = seq(sum_as_ALL$min_bwLOG,sum_as_ALL$max_bwLOG, 0.5))
  data.plotAS_ER$model_predFE <- predict(as_mod_ER, re.form = NA, newdata = data.plotAS_ER)
  conf.int_AS_ER <- confint(as_mod_ER, method = "boot", FUN=function(x)predict(x, data.plotAS_ER, re.form=NA), nsim=10)
  data.plotAS_ER$CI_2.5<-conf.int_AS_ER[,1]
  data.plotAS_ER$CI_97.5<-conf.int_AS_ER[,2]

  # MMR or AMR  optimal - interaction model; plot/predict full range (All data coverage for reference)
  data.plotAMRint_ER <- expand.grid(tempTestK1000 = seq(C40inTempTestK1000,C0inTempTestK1000, 0.05), lnBWg= seq(sum_amr_ALL$min_bwLOG,sum_amr_ALL$max_bwLOG, 0.5))
  data.plotAMRint_ER$model_predFE <- predict(amr_mod_ER, re.form = NA, newdata = data.plotAMRint_ER)
  conf.int_AMRint_ER <- confint(as_mod_ER, method = "boot", FUN=function(x)predict(x, data.plotAMRint_ER, re.form=NA), nsim=10)
  data.plotAMRint_ER$CI_2.5<-conf.int_AMRint_ER[,1]
  data.plotAMRint_ER$CI_97.5<-conf.int_AMRint_ER[,2]
  data.plotAMRint_ER$tempTestK1000_inC <- ((1000/data.plotAMRint_ER$tempTestK1000))-275.15

  # All temps for the Arrhenius plot:
  # RMR - warm temp grid; plot/predict full range (All data coverage for reference)
  data.plotRMR_warm<-expand.grid(tempTestK1000 = seq(C40inTempTestK1000, C0inTempTestK1000 , 0.05), lnBWg = seq(sum_rmr_ALL$min_bwLOG,sum_rmr_ALL$max_bwLOG, 0.5))
  data.plotRMR_warm$model_predFE <- predict(rmr_mod_W, re.form = NA, newdata = data.plotRMR_warm)
  conf.int_RMR_warm <- confint(rmr_mod_W, method = "boot", FUN=function(x)predict(x, data.plotRMR_warm, re.form=NA), nsim=10)
  data.plotRMR_warm$CI_2.5<-conf.int_RMR_warm[,1]
  data.plotRMR_warm$CI_97.5<-conf.int_RMR_warm[,2]

  # MMR - warm temp plot grid; plot/predict full range (All data coverage for reference)
  data.plotAMR_warm<-expand.grid(tempTestK1000 = seq(C40inTempTestK1000, C0inTempTestK1000 , 0.05), lnBWg = seq(sum_amr_ALL$min_bwLOG,sum_amr_ALL$max_bwLOG, 0.5))
  data.plotAMR_warm$model_predFE <- predict(amr_mod_W, re.form = NA, newdata = data.plotAMR_warm)
  conf.int_AMR_warm <- confint(amr_mod_W, method = "boot", FUN=function(x)predict(x, data.plotAMR_warm, re.form=NA), nsim=10)
  data.plotAMR_warm$CI_2.5<-conf.int_AMR_warm[,1]
  data.plotAMR_warm$CI_97.5<-conf.int_AMR_warm[,2]

  # FAS - warm temp plot grid; plot/predict full range (All data coverage for reference)
  data.plotFAS_warm <- expand.grid(lnBWg = seq(sum_fas_ALL$min_bwLOG,sum_fas_ALL$max_bwLOG, 0.5), tempTest = 20)
  data.plotFAS_warm$model_predFE <- predict(fas_mod_W, re.form = NA, newdata = data.plotFAS_warm)
  conf.int_FAS_warm <- confint(fas_mod_W, method = "boot", FUN=function(x)predict(x, data.plotFAS_warm, re.form=NA), nsim=10)
  data.plotFAS_warm$CI_2.5<-conf.int_FAS_warm[,1]
  data.plotFAS_warm$CI_97.5<-conf.int_FAS_warm[,2]

  # AS - warm temp plot grid; plot/predict full range (All data coverage for reference)
  data.plotAS_warm <- expand.grid(lnBWg = seq(sum_as_ALL$min_bwLOG,sum_as_ALL$max_bwLOG, 0.5), tempTestK1000 = seq(C40inTempTestK1000, C0inTempTestK1000 , 0.1))
  data.plotAS_warm$model_predFE <- predict(as_mod_W, re.form = NA, newdata = data.plotAS_warm)
  conf.int_AS_warm <- confint(as_mod_W, method = "boot", FUN=function(x)predict(x, data.plotAS_warm, re.form=NA), nsim=10)
  data.plotAS_warm$CI_2.5<-conf.int_AS_warm[,1]
  data.plotAS_warm$CI_97.5<-conf.int_AS_warm[,2]

  # CIs for all model parameters
  CI.amr.ER<-as.data.frame(confint.merMod(amr_mod_ER, level = 0.90))
  CI.amr.ER$var<-rownames(CI.amr.ER)
  CI.amr.ER$MR<-"MMR"
  CI.amr.ER$temp_cat<-"ER"

  CI.rmr.ER<-as.data.frame(confint.merMod(rmr_mod_ER, level = 0.90))
  CI.rmr.ER$var<-rownames(CI.rmr.ER)
  CI.rmr.ER$MR<-"RMR"
  CI.rmr.ER$temp_cat<-"ER"

  CI.fas.ER<-as.data.frame(confint.merMod(fas_mod_ER, level = 0.90))
  CI.fas.ER$var<-rownames(CI.fas.ER)
  CI.fas.ER$MR<-"FAS"
  CI.fas.ER$temp_cat<-"ER"

  CI.as.ER<-as.data.frame(confint.merMod(as_mod_ER, level = 0.90))
  CI.as.ER$var<-rownames(CI.as.ER)
  CI.as.ER$MR<-"AS"
  CI.as.ER$temp_cat<-"ER"

  CI.amr.W<-as.data.frame(confint.merMod(amr_mod_W, level = 0.90))
  CI.amr.W$var<-rownames(CI.amr.W)
  CI.amr.W$MR<-"MMR"
  CI.amr.W$temp_cat<-"W"

  CI.rmr.W<-as.data.frame(confint.merMod(rmr_mod_W, level = 0.90))
  CI.rmr.W$var<-rownames(CI.rmr.W)
  CI.rmr.W$MR<-"RMR"
  CI.rmr.W$temp_cat<-"W"

  CI.fas.W<-as.data.frame(confint.merMod(fas_mod_W, level = 0.90))
  CI.fas.W$var<-rownames(CI.fas.W)
  CI.fas.W$MR<-"FAS"
  CI.fas.W$temp_cat<-"W"

  CI.as.W<-as.data.frame(confint.merMod(as_mod_W, level = 0.90))
  CI.as.W$var<-rownames(CI.as.W)
  CI.as.W$MR<-"AS"
  CI.as.W$temp_cat<-"W"

  sum_CItable<-rbind(CI.amr.ER, CI.rmr.ER, CI.fas.ER, CI.as.ER, # ecol relev
                     CI.amr.W, CI.rmr.W, CI.fas.W, CI.as.W) # warm

  # Saving all export files
  write.csv(file = "./Data_exports/Phylo/Table_CIsummary.csv", sum_CItable, row.names=TRUE)

  write.csv(file = paste("./Data_exports/Phylo/dataPred_RMR_er_", Sys.Date(), ".csv", sep=""),  data.plotRMR_ER, row.names = F)
  write.csv(file = paste("./Data_exports/Phylo/dataPred_AMR_er_", Sys.Date(), ".csv", sep=""),  data.plotAMRint_ER, row.names = F)
  write.csv(file = paste("./Data_exports/Phylo/dataPred_AS_er_", Sys.Date(), ".csv", sep=""),  data.plotAS_ER, row.names = F)
  write.csv(file = paste("./Data_exports/Phylo/dataPred_FAS_er_", Sys.Date(), ".csv", sep=""),  data.plotFAS_ER, row.names = F)

  write.csv(file = paste("./Data_exports/Phylo/dataPred_RMR_warm_", Sys.Date(), ".csv", sep=""),  data.plotRMR_warm, row.names = F)
  write.csv(file = paste("./Data_exports/Phylo/dataPred_AMR_warm_", Sys.Date(), ".csv", sep=""),  data.plotAMR_warm, row.names = F)
  write.csv(file = paste("./Data_exports/Phylo/dataPred_AS_warm_", Sys.Date(), ".csv", sep=""),  data.plotAS_warm, row.names = F)
  write.csv(file = paste("./Data_exports/Phylo/dataPred_FAS_warm_", Sys.Date(), ".csv", sep=""),  data.plotFAS_warm, row.names = F)


}else{
 
  k<-(8.62*10^(-5)) # Boltzmann's constant

  # ECOL RELEV
  RMR_slope<-round(fixef(rmr_mod_ER)[2],3)
  RMR_int<-round(fixef(rmr_mod_ER)[1],3)
  FAS_slope<-round(fixef(fas_mod_ER)[2],3)
  FAS_int<-round(fixef(fas_mod_ER)[1],3)
  AS_slope<-round(fixef(as_mod_ER)[2],3)
  AS_int<-round(fixef(as_mod_ER)[1],3)

  # WARM
  MMR_slope_w<-round(fixef(amr_mod_W)[2],3)
  MMR_int_w<-round(fixef(amr_mod_W)[1],3)
  RMR_slope_w<-round(fixef(rmr_mod_W)[2],3)
  RMR_int_w<-round(fixef(rmr_mod_W)[1],3)
  FAS_slope_w<-round(fixef(fas_mod_W)[2],3)
  FAS_int_w<-round(fixef(fas_mod_W)[1],3)
  AS_slope_w<-round(fixef(as_mod_W)[2],3)
  AS_int_w<-round(fixef(as_mod_W)[1],3)

  C20inTempTestK1000<-1000/(((20+273.15)))
  C10inTempTestK1000<-1000/((10+273.15))
  C0inTempTestK1000<-1000/((0+273.15))
  C30inTempTestK1000<-1000/((30+273.15))
  C40inTempTestK1000<-1000/((40+273.15))

  # predicted slopes for MMR interaction
  AMR.slopes<-(emtrends(amr_mod_ER, pairwise ~ tempTestK1000,
                        pbkrtest.limit = 4000, var="lnBWg",
                        at=list(tempTestK1000=c(C0inTempTestK1000,
                                                C10inTempTestK1000,
                                                C20inTempTestK1000,
                                                C30inTempTestK1000,
                                                C40inTempTestK1000)))) # 0, 10, 20, 30, 40 C
  AMR.slopes<-as.data.frame(AMR.slopes$emtrends)
  AMR.slopes$tempTestK1000_inC<-((1000/AMR.slopes$tempTestK1000))-273.15
  AMR.slopes$performance<-"MMR"
  AMR.slopes$temp_categ<-"er"
  AMR.slopes$`(Intercept)` <- NA
  AMR.slopes$tempTestK1000 <- NA
  AMR.slopes$`E(ev)`<-NA
  # fixed effects activation energies:
  MMR_E_ER<-(emtrends(amr_mod_ER, pairwise ~ lnBWg , pbkrtest.limit = 4000,
                      var = "tempTestK1000",
                      at = list(lnBWg = c(log(1),log(10), log(100),log(1000)))))
  MMR_E_ER<-as.data.frame(MMR_E_ER$emtrends)
  MMR_E_ER$performance<-"MMR"
  MMR_E_ER$temp_categ<-"er"
  MMR_E_ER$lnBWg.trend<-NA
  MMR_E_ER$`(Intercept)` <- NA
  MMR_E_ER$tempTestK1000_inC <- NA

  RMR_E_ER<-round(fixef(rmr_mod_ER)[3],3)
  AS_E_ER<-round(fixef(as_mod_ER)[3],3)
  MMR_E_W<-round(fixef(amr_mod_ER)[3],3)
  RMR_E_W<-round(fixef(rmr_mod_W)[3],3)
  AS_E_W<-round(fixef(as_mod_W)[3],3)

  # activation energies:
  # -1*(E*k)*1000
  MMR_E_ER$MMR_E_ER_eV<--1*(MMR_E_ER$tempTestK1000.trend)*k*1000
  RMR_E_ER_eV<--1*(RMR_E_ER)*k*1000
  AS_E_ER_eV<--1*(AS_E_ER)*k*1000

  MMR_E_W_eV<--1*(MMR_E_W)*k*1000
  RMR_E_W_eV<--1*(RMR_E_W)*k*1000
  AS_E_W_eV<--1*(AS_E_W)*k*1000

  # dataframe with slopes and intercepts
  scaling.params<-as.data.frame(matrix(ncol = 7, nrow = 0))
  colnames(scaling.params)<-c("performance", "temp_categ",
                              "lnBWg",
                              "(Intercept)",
                              "tempTestK1000",
                              "E(ev)", "tempTestK1000_inC")

  RMR_er_row<-as.data.frame(t(c("performance" = "RMR",
                  "temp_categ" = "er",
                  RMR_slope, RMR_int, RMR_E_ER,
                  "E(ev)" = RMR_E_ER_eV,
                  "tempTestK1000_inC" = NA)))
  RMR_W_row<-as.data.frame(t(c("performance" = "RMR",
                  "temp_categ" = "warm",
                  RMR_slope_w, RMR_int_w, RMR_E_W,
                  "E(ev)" = RMR_E_W_eV,
                  "tempTestK1000_inC" = NA)))
  FAS_er_row<-as.data.frame(t(c("performance" = "FAS",
                  "temp_categ" = "er",
                  FAS_slope, FAS_int,
                  "tempTestK1000" = NA,
                  "E(ev)" = NA,
                  "tempTestK1000_inC" = NA)))
  FAS_W_row<-as.data.frame(t(c("performance" = "FAS",
                  "temp_categ" = "warm",
                  FAS_slope_w, FAS_int_w,
                  "tempTestK1000" = NA,
                  "E(ev)" = NA,
                  "tempTestK1000_inC" = NA)))
  AS_er_row<-as.data.frame(t(c("performance" = "AS",
                  "temp_categ" = "er",
                  AS_slope, AS_int, AS_E_ER,
                  "E(ev)" = AS_E_ER_eV,
                  "tempTestK1000_inC" = NA)))
  AS_W_row<-as.data.frame(t(c("performance" = "AS",
                  "temp_categ" = "warm",
                  AS_slope_w, AS_int_w, AS_E_W,
                  "E(ev)" = AS_E_W_eV,
                  "tempTestK1000_inC" = NA)))

  MMR_W_row<-as.data.frame(t(c("performance" = "MMR",
                  "temp_categ" = "warm",
                  MMR_slope_w, MMR_int_w, MMR_E_W,
                  "E(ev)" = MMR_E_W_eV,
                  "tempTestK1000_inC" = NA)))

  MMR_er_row<-AMR.slopes[, c("performance", "temp_categ", "lnBWg.trend",
                              "(Intercept)", "tempTestK1000","E(ev)","tempTestK1000_inC")]
  MMR_er_row2<-MMR_E_ER[, c("performance", "temp_categ", "lnBWg.trend",
                              "(Intercept)", "tempTestK1000.trend","MMR_E_ER_eV","tempTestK1000_inC")]
  colnames(MMR_er_row)<-colnames(scaling.params)
  colnames(MMR_er_row2)<-colnames(scaling.params)
  colnames(MMR_W_row)<-colnames(scaling.params)
  colnames(RMR_er_row)<-colnames(scaling.params)
  colnames(RMR_W_row)<-colnames(scaling.params)
  colnames(AS_er_row)<-colnames(scaling.params)
  colnames(AS_W_row)<-colnames(scaling.params)
  colnames(FAS_er_row)<-colnames(scaling.params)
  colnames(FAS_W_row)<-colnames(scaling.params)

  scaling.params<-rbind(MMR_er_row, MMR_er_row2, MMR_W_row, 
                        RMR_er_row, RMR_W_row,
                        AS_er_row, AS_W_row, 
                        FAS_er_row, FAS_W_row)
  
  # MMR_W_ER$MMR_E_ER_eV<--1*(MMR_E_ER$tempTestK1000.trend)*k*1000
  # RMR_W_ER_eV<--1*(RMR_W_ER_eV)*k*1000

  # Figure parameters:
  summarise.dataframes<- function(data.summarise, export.name){
    dataframe.new <- data.summarise %>%
      dplyr:::summarise(min_bw = min(BW_g), max_bw = max(BW_g), mean_bw = mean(BW_g),
                min_bwLOG = min(lnBWg), max_bwLOG = max(lnBWg), mean_bwLOG = mean(lnBWg),
                min_temp = min(tempTest), max_temp = max(tempTest), mean_temp = mean(tempTest),
                min_tempArrh = min(tempTestK1000), max_tempArrh = max(tempTestK1000), mean_tempArrh = mean(tempTestK1000),
                n = length(BW_g)) %>%
      as.data.frame()

    assign(export.name, dataframe.new)
  }

  sum_amr_ER<-summarise.dataframes(data.amrER, "sum_amr_ER")
  sum_rmr_ER<-summarise.dataframes(data.rmrER, "sum_rmr_ER")
  sum_fas_ER<-summarise.dataframes(data.fasER, "sum_fas_ER")
  sum_as_ER<-summarise.dataframes(data.asER, "sum_as_ER")
  sum_amr_warm<-summarise.dataframes(data.amr.test, "sum_amr_warm")
  sum_rmr_warm<-summarise.dataframes(data.rmr.test, "sum_rmr_warm")
  sum_fas_warm<-summarise.dataframes(data.fas.test, "sum_fas_warm")
  sum_as_warm<-summarise.dataframes(data.as.test, "sum_as_warm")
  sum_amr_ALL<-summarise.dataframes(data.amr, "sum_amr_ALL")
  sum_rmr_ALL<-summarise.dataframes(data.rmr, "sum_rmr_ALL")
  sum_fas_ALL<-summarise.dataframes(data.fas, "sum_fas_ALL")
  sum_as_ALL<-summarise.dataframes(data.as, "sum_as_ALL")
  
  # row-bind in one dataframe 
  sum_data<-rbind(sum_amr_ER, sum_amr_warm, 
                  sum_rmr_ER, sum_rmr_warm,
                  sum_as_ER, sum_as_warm,
                  sum_fas_ER, sum_fas_warm)
  
  sum_data$MR<-c("MMR", "MMR", "RMR", "RMR", "AS", "AS", "FAS", "FAS")
  sum_data$Temp<-c("er", "warm","er", "warm","er", "warm","er", "warm")

  data.amr$tempTestK1<-1/data.amr$tempTestK
  data.amr$tempTestK1000<-1000/data.amr$tempTestK
  # 20 deg
  # -1*(RMR_E*k)*1000
  # tempTestK1000

  C20inTempTestK1000<-1000/(((20+273.15)))
  C10inTempTestK1000<-1000/((10+273.15))
  C0inTempTestK1000<-1000/((0+273.15))
  C30inTempTestK1000<-1000/((30+273.15))
  C40inTempTestK1000<-1000/((40+273.15))

  # 40 = (1000 / (C40inTempTestK1000))- 273.15
  # (1000 / (C40inTempTestK1000))- 273.15
  # ((1/C20inTempTestkT)/k)-273.15
  

  message("Data to plot predicted data and scaling are imported from './Data_exports/'")
  
  sum_CItable<-read.csv(file = "./Data_exports/Phylo/Table_CIsummary.csv")

  data.plotRMR_ER<-read.csv(file = paste("./Data_exports/Phylo/dataPred_RMR_er_", date.imports, ".csv", sep=""))
  data.plotAMRint_ER<-read.csv(file = paste("./Data_exports/Phylo/dataPred_AMR_er_", date.imports, ".csv", sep=""))
  data.plotAS_ER<-read.csv(file = paste("./Data_exports/Phylo/dataPred_AS_er_", date.imports, ".csv", sep=""))
  data.plotFAS_ER<-read.csv(file = paste("./Data_exports/Phylo/dataPred_FAS_er_", date.imports, ".csv", sep=""))

  data.plotRMR_warm<-read.csv(paste("./Data_exports/Phylo/dataPred_RMR_warm_", date.imports, ".csv", sep=""))
  data.plotAMR_warm<-read.csv(file = paste("./Data_exports/Phylo/dataPred_AMR_warm_", date.imports, ".csv", sep=""))
  data.plotAS_warm<-read.csv(file = paste("./Data_exports/Phylo/dataPred_AS_warm_", date.imports, ".csv", sep=""))
  data.plotFAS_warm<-read.csv(file = paste("./Data_exports/Phylo/dataPred_FAS_warm_", date.imports, ".csv", sep=""))

}




