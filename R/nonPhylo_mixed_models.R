
library(lme4)
library(emmeans)
library(pryr)
library(ggformat2) # from github
library(here)

## Source the data ------------
source(here("R", "get_data_temp.R"))
source(here("R", "mixed_model_outputs.R"))

# function to order model selection based on the lowest BIC score
BICdelta<-function(BICtable){
  BIC.t <- BICtable [order(BICtable$BIC), ]
  BIC.t$delta <- round(abs(BIC.t$BIC[1] -  BIC.t$BIC), 5)
  return( BIC.t)
}

# 1. DATA FOR MODELS ------
data.list<-get_data_temp(data.amr = "./Data/Fish_AMR_temp_dataset_mar2022.csv",
                                 data.rmr = "./Data/Fish_RMR_temp_dataset_mar2022.csv",
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
data.fas.test<-rbind(data.fasAC, data.fasAM)
data.as.test<-rbind(data.asAC, data.asAM)
data.fas.test<-data.fas.test[c(!is.na(data.fas.test$FAS) & is.finite(data.fas.test$FAS)) , ]
data.as.test<-data.as.test[c(!is.na(data.as.test$lnAS) & is.finite(data.as.test$lnAS)) , ]

# get model data set specific phylogentics model matrixes -------
data.rmrER<-droplevels(data.rmrER)
data.rmr.test<-droplevels(data.rmr.test)
data.amrER<-droplevels(data.amrER)
data.amr.test<-droplevels(data.amr.test)
data.asER<-droplevels(data.asER)
data.as.test<-droplevels(data.as.test)
data.fasER<-droplevels(data.fasER)
data.fas.test<-droplevels(data.fas.test)

## RMR: ecologically relevant -----------------
RMR_model0 <- lmer(lnRMR ~ lnBWg + tempTestK1000 + (1|species) , data=data.rmrER, REML=FALSE)
RMR_model0.POLY <- lmer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species), data=data.rmrER, REML=FALSE)
RMR_model0intPOLY <- lmer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (1|species), data=data.rmrER, REML=FALSE)
RMR_model0int <- lmer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) , data=data.rmrER, REML=FALSE)

RMR_model1 <- lmer(lnRMR ~ lnBWg + tempTestK1000 + (1|species) + (1|species:trial), data=data.rmrER, REML=FALSE)
RMR_model1.POLY <- lmer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.rmrER, REML=FALSE)
RMR_model1intPOLY <- lmer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.rmrER, REML=FALSE)
RMR_model1int <- lmer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) + (1|species:trial), data=data.rmrER, REML=FALSE)

RMR_model2 <- lmer(lnRMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.rmrER, REML=FALSE)
RMR_model2int <- lmer(lnRMR ~ lnBWg * tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.rmrER, REML=FALSE)
RMR_model2.POLY <- lmer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.rmrER, REML=FALSE)
# ^^^ singular fit 
RMR_model2intPOLY <- lmer(lnRMR ~ lnBWg * poly(tempTestK1000,2)  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.rmrER, REML=FALSE)

RMR_model4 <- lmer(lnRMR ~ lnBWg + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE)
RMR_model4int <- lmer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE)
RMR_model4.POLY <- lmer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE)
RMR_model4intPOLY <- lmer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE)

RMR_model5 <- lmer(lnRMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE)
RMR_model5int <- lmer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE)
RMR_model5.POLY <- lmer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE)
RMR_model5intPOLY <- lmer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE)

RMR_model6 <- lmer(lnRMR ~ lnBWg + tempTestK1000 + (lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE)
RMR_model6int <- lmer(lnRMR ~ lnBWg * tempTestK1000+ (lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE)
RMR_model6.POLY <- lmer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE)
RMR_model6intPOLY <- lmer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE)

### BIC --------------
RMR_BIC.nP<-BICdelta(BIC(RMR_model0, RMR_model0int, RMR_model0.POLY, RMR_model0intPOLY,
                RMR_model1, RMR_model1int, RMR_model1.POLY, RMR_model1intPOLY,
                RMR_model2, RMR_model2int, RMR_model2.POLY, RMR_model2intPOLY,
                RMR_model4, RMR_model4int, RMR_model4.POLY, RMR_model4intPOLY,
                RMR_model5, RMR_model5int, RMR_model5.POLY, RMR_model5intPOLY,
                RMR_model6, RMR_model6int, RMR_model6.POLY, RMR_model6intPOLY))

## MMR ecologically relevant -----------------
MMR_model0 <- lmer(lnAMR ~ lnBWg + tempTestK1000 + (1|species) , data=data.amrER, REML=FALSE)
MMR_model0.POLY <- lmer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species), data=data.amrER, REML=FALSE)
MMR_model0intPOLY <- lmer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (1|species), data=data.amrER, REML=FALSE)
MMR_model0int <- lmer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) , data=data.amrER, REML=FALSE)

MMR_model1 <- lmer(lnAMR ~ lnBWg + tempTestK1000 + (1|species) + (1|species:trial), data=data.amrER, REML=FALSE)
MMR_model1.POLY <- lmer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.amrER, REML=FALSE)
MMR_model1intPOLY <- lmer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.amrER, REML=FALSE)
MMR_model1int <- lmer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) + (1|species:trial), data=data.amrER, REML=FALSE)

MMR_model2 <- lmer(lnAMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.amrER, REML=FALSE)
MMR_model2int <- lmer(lnAMR ~ lnBWg * tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.amrER, REML=FALSE)
MMR_model2.POLY <- lmer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.amrER, REML=FALSE)
MMR_model2intPOLY <- lmer(lnAMR ~ lnBWg * poly(tempTestK1000,2)  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.amrER, REML=FALSE)

MMR_model4 <- lmer(lnAMR ~ lnBWg + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE)
MMR_model4int <- lmer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE)
MMR_model4.POLY <- lmer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE)
MMR_model4intPOLY <- lmer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE)

MMR_model5 <- lmer(lnAMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE)
MMR_model5int <- lmer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE)
MMR_model5.POLY <- lmer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE)
MMR_model5intPOLY <- lmer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE)

MMR_model6 <- lmer(lnAMR ~ lnBWg + tempTestK1000 + (lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE)
MMR_model6int <- lmer(lnAMR ~ lnBWg * tempTestK1000+ (lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE)
MMR_model6.POLY <- lmer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE)
MMR_model6intPOLY <- lmer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE)

### BIC --------------
MMR_BIC.nP<-BICdelta(BIC(MMR_model0, MMR_model0int, MMR_model0.POLY, MMR_model0intPOLY,
                MMR_model1, MMR_model1int, MMR_model1.POLY, MMR_model1intPOLY,
                MMR_model2, MMR_model2int, MMR_model2.POLY, MMR_model2intPOLY,
                MMR_model4, MMR_model4int, MMR_model4.POLY, MMR_model4intPOLY,
                MMR_model5, MMR_model5int, MMR_model5.POLY, MMR_model5intPOLY,
                MMR_model6, MMR_model6int, MMR_model6.POLY, MMR_model6intPOLY))

## AAS ecologically relevant ------------------
AS_model0 <- lmer(lnAS ~ lnBWg + tempTestK1000 + (1|species) , data=data.asER, REML=FALSE)
AS_model0.POLY <- lmer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species), data=data.asER, REML=FALSE)
AS_model0intPOLY <- lmer(lnAS ~ lnBWg * poly(tempTestK1000,2) + (1|species), data=data.asER, REML=FALSE)
AS_model0int <- lmer(lnAS ~ lnBWg * tempTestK1000 + (1|species) , data=data.asER, REML=FALSE)

AS_model1 <- lmer(lnAS ~ lnBWg + tempTestK1000 + (1|species) + (1|species:trial), data=data.asER, REML=FALSE)
AS_model1.POLY <- lmer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.asER, REML=FALSE)
AS_model1intPOLY <- lmer(lnAS ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.asER, REML=FALSE)
AS_model1int <- lmer(lnAS ~ lnBWg * tempTestK1000 + (1|species) + (1|species:trial), data=data.asER, REML=FALSE)

AS_model2 <- lmer(lnAS ~ lnBWg + tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.asER, REML=FALSE)
# ^^^ convergence issues: Model failed to converge with max|grad| = 0.00350204 (tol = 0.002, component 1)
AS_model2int <- lmer(lnAS ~ lnBWg * tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.asER, REML=FALSE)
AS_model2.POLY <- lmer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.asER, REML=FALSE)
# ^^^ convergence issues: Model failed to converge: degenerate  Hessian with 1 negative eigenvalues
AS_model2intPOLY <- lmer(lnAS ~ lnBWg * poly(tempTestK1000,2)  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.asER, REML=FALSE)

AS_model4 <- lmer(lnAS ~ lnBWg + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE)
AS_model4int <- lmer(lnAS ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE)
AS_model4.POLY <- lmer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE)
AS_model4intPOLY <- lmer(lnAS ~ lnBWg * poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE)

AS_model5 <- lmer(lnAS ~ lnBWg + tempTestK1000  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE)
AS_model5int <- lmer(lnAS ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE)
AS_model5.POLY <- lmer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE)
AS_model5intPOLY <- lmer(lnAS ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE)

AS_model6 <- lmer(lnAS ~ lnBWg + tempTestK1000 + (lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE)
AS_model6int <- lmer(lnAS ~ lnBWg * tempTestK1000+ (lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE)
AS_model6.POLY <- lmer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE)
AS_model6intPOLY <- lmer(lnAS ~ lnBWg * poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE)

## BIC --------
AS_BIC.nP<-BICdelta(BIC(AS_model0, AS_model0int, AS_model0.POLY, AS_model0intPOLY,
                AS_model1, AS_model1int, AS_model1.POLY, AS_model1intPOLY,
                AS_model2, AS_model2int, AS_model2.POLY, AS_model2intPOLY,
                AS_model4, AS_model4int, AS_model4.POLY, AS_model4intPOLY,
                AS_model5, AS_model5int, AS_model5.POLY, AS_model5intPOLY,
                AS_model6, AS_model6int, AS_model6.POLY, AS_model6intPOLY))


## FAS ecologically relevant ----------------
FAS_model0 <- lmer(log(FAS) ~ lnBWg + tempTest + (1|species) , data=data.fasER, REML=FALSE)
FAS_model0.POLY <- lmer(log(FAS) ~ lnBWg + poly(tempTest,2) + (1|species), data=data.fasER, REML=FALSE)
FAS_model0intPOLY <- lmer(log(FAS) ~ lnBWg * poly(tempTest,2) + (1|species), data=data.fasER, REML=FALSE)
FAS_model0int <- lmer(log(FAS) ~ lnBWg * tempTest + (1|species) , data=data.fasER, REML=FALSE)

FAS_model1 <- lmer(log(FAS) ~ lnBWg + tempTest + (1|species) + (1|species:trial), data=data.fasER, REML=FALSE)
FAS_model1.POLY <- lmer(log(FAS) ~ lnBWg + poly(tempTest,2) + (1|species) + (1|species:trial), data=data.fasER, REML=FALSE)
FAS_model1intPOLY <- lmer(log(FAS) ~ lnBWg * poly(tempTest,2) + (1|species) + (1|species:trial), data=data.fasER, REML=FALSE)
# ^^^ convergence issues: Model failed to converge with max|grad| = 0.00551275 (tol = 0.002, component 1
FAS_model1int <- lmer(log(FAS) ~ lnBWg * tempTest + (1|species) + (1|species:trial), data=data.fasER, REML=FALSE)

FAS_model2 <- lmer(log(FAS) ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE)
# ^^^ convergence issues: Model failed to converge with max|grad| =
FAS_model2int <- lmer(log(FAS) ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE)
# ^^^ convergence issues: Model failed to converge with max|grad| =
FAS_model2.POLY <- lmer(log(FAS) ~ lnBWg + poly(tempTest,2) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE)
FAS_model2intPOLY <- lmer(log(FAS) ~ lnBWg * poly(tempTest,2)  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE)

FAS_model4 <- lmer(log(FAS) ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE)
FAS_model4int <- lmer(log(FAS) ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE)
FAS_model4.POLY <- lmer(log(FAS) ~ lnBWg + poly(tempTest,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE)
FAS_model4intPOLY <- lmer(log(FAS) ~ lnBWg * poly(tempTest,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE)

FAS_model5 <- lmer(log(FAS) ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE)
FAS_model5int <- lmer(log(FAS) ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE)
FAS_model5.POLY <- lmer(log(FAS) ~ lnBWg + poly(tempTest,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE)
FAS_model5intPOLY <- lmer(log(FAS) ~ lnBWg * poly(tempTest,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE)

FAS_model6 <- lmer(lnFAS ~ lnBWg + tempTestK1000 + (lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE)
FAS_model6int <- lmer(lnFAS ~ lnBWg * tempTestK1000+ (lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE)
FAS_model6.POLY <- lmer(lnFAS ~ lnBWg + poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE)
FAS_model6intPOLY <- lmer(lnFAS ~ lnBWg * poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE)
# ^^^ singular fit 

## BIC -------------
FAS_BIC.nP<-BICdelta(BIC(FAS_model0, FAS_model0int, FAS_model0.POLY, FAS_model0intPOLY,
                FAS_model1, FAS_model1int, FAS_model1.POLY, FAS_model1intPOLY,
                FAS_model2, FAS_model2int, FAS_model2.POLY, FAS_model2intPOLY,
                FAS_model4, FAS_model4int, FAS_model4.POLY, FAS_model4intPOLY,
                FAS_model5, FAS_model5int, FAS_model5.POLY, FAS_model5intPOLY,
                FAS_model6, FAS_model6int, FAS_model6.POLY, FAS_model6intPOLY))


# RMR warm ------------------
RMR_W_model0 <- lmer(lnRMR ~ lnBWg + tempTestK1000 + (1|species) , data=data.rmr.test, REML=FALSE)
RMR_W_model0.POLY <- lmer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species), data=data.rmr.test, REML=FALSE)
RMR_W_model0intPOLY <- lmer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (1|species), data=data.rmr.test, REML=FALSE)
RMR_W_model0int <- lmer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) , data=data.rmr.test, REML=FALSE)

RMR_W_model1 <- lmer(lnRMR ~ lnBWg + tempTestK1000 + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)
RMR_W_model1.POLY <- lmer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)
RMR_W_model1intPOLY <- lmer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)
RMR_W_model1int <- lmer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)

RMR_W_model2 <- lmer(lnRMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)
RMR_W_model2int <- lmer(lnRMR ~ lnBWg * tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)
# ^^^ singular fit 
RMR_W_model2.POLY <- lmer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)
# ^^^ singular fit 
RMR_W_model2intPOLY <- lmer(lnRMR ~ lnBWg * poly(tempTestK1000,2)  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)

RMR_W_model4 <- lmer(lnRMR ~ lnBWg + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)
RMR_W_model4int <- lmer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)
# ^^^ singular fit 
RMR_W_model4.POLY <- lmer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)
RMR_W_model4intPOLY <- lmer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)
# ^^^ singular fit 

RMR_W_model5 <- lmer(lnRMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE)
# ^^^ singular fit 
RMR_W_model5int <- lmer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE)
# ^^^ singular fit 
RMR_W_model5.POLY <- lmer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE)
# ^^^ singular fit 
RMR_W_model5intPOLY <- lmer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE)
# ^^^ singular fit 

RMR_W_model6 <- lmer(lnRMR ~ lnBWg + tempTestK1000 + (lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)
# ^^^ singular fit 
RMR_W_model6int <- lmer(lnRMR ~ lnBWg * tempTestK1000+ (lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)
# ^^^ convergence issues: Model failed to converge with max|grad|
RMR_W_model6.POLY <- lmer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)
# ^^^ singular fit 
RMR_W_model6intPOLY <- lmer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)
# ^^^ singular fit 


## BIC -------------
RMR_W_BIC.nP<-BICdelta(BIC(RMR_W_model0, RMR_W_model0int, RMR_W_model0.POLY, RMR_W_model0intPOLY,
                RMR_W_model1, RMR_W_model1int, RMR_W_model1.POLY, RMR_W_model1intPOLY,
                RMR_W_model2, RMR_W_model2int, RMR_W_model2.POLY, RMR_W_model2intPOLY,
                RMR_W_model4, RMR_W_model4int, RMR_W_model4.POLY, RMR_W_model4intPOLY,
                RMR_W_model5, RMR_W_model5int, RMR_W_model5.POLY, RMR_W_model5intPOLY,
                RMR_W_model6, RMR_W_model6int, RMR_W_model6.POLY, RMR_W_model6intPOLY))

# AMR / warm temps --------------
MMR_W_model0 <- lmer(lnAMR ~ lnBWg + tempTestK1000 + (1|species) , data=data.amr.test, REML=FALSE)
MMR_W_model0.POLY <- lmer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species), data=data.amr.test, REML=FALSE)
MMR_W_model0intPOLY <- lmer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (1|species), data=data.amr.test, REML=FALSE)
MMR_W_model0int <- lmer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) , data=data.amr.test, REML=FALSE)

MMR_W_model1 <- lmer(lnAMR ~ lnBWg + tempTestK1000 + (1|species) + (1|species:trial), data=data.amr.test, REML=FALSE)
MMR_W_model1.POLY <- lmer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.amr.test, REML=FALSE)
MMR_W_model1intPOLY <- lmer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.amr.test, REML=FALSE)
MMR_W_model1int <- lmer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) + (1|species:trial), data=data.amr.test, REML=FALSE)

MMR_W_model2 <- lmer(lnAMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.amr.test, REML=FALSE)
MMR_W_model2int <- lmer(lnAMR ~ lnBWg * tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.amr.test, REML=FALSE)
# ^^^ convergence issues: Model failed to converge: degenerate  Hessian with 1 negative eigenvalues
MMR_W_model2.POLY <- lmer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.amr.test, REML=FALSE)
MMR_W_model2intPOLY <- lmer(lnAMR ~ lnBWg * poly(tempTestK1000,2)  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.amr.test, REML=FALSE)
# ^^^ convergence issues: Model failed to converge: degenerate  Hessian with 1 negative eigenvalues

MMR_W_model4 <- lmer(lnAMR ~ lnBWg + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE)
MMR_W_model4int <- lmer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE)
MMR_W_model4.POLY <- lmer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE)
MMR_W_model4intPOLY <- lmer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE)

MMR_W_model5 <- lmer(lnAMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE)
MMR_W_model5int <- lmer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE)
MMR_W_model5.POLY <- lmer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE)
MMR_W_model5intPOLY <- lmer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE)

MMR_W_model6 <- lmer(lnAMR ~ lnBWg + tempTestK1000 + (lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE)
MMR_W_model6int <- lmer(lnAMR ~ lnBWg * tempTestK1000+ (lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE)
MMR_W_model6.POLY <- lmer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE)
MMR_W_model6intPOLY <- lmer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE)

## BIC -------------
MMR_W_BIC.nP<-BICdelta(BIC(MMR_W_model0, MMR_W_model0int, MMR_W_model0.POLY, MMR_W_model0intPOLY,
                MMR_W_model1, MMR_W_model1int, MMR_W_model1.POLY, MMR_W_model1intPOLY,
                MMR_W_model2, MMR_W_model2int, MMR_W_model2.POLY, MMR_W_model2intPOLY,
                MMR_W_model4, MMR_W_model4int, MMR_W_model4.POLY, MMR_W_model4intPOLY,
                MMR_W_model5, MMR_W_model5int, MMR_W_model5.POLY, MMR_W_model5intPOLY,
                MMR_W_model6, MMR_W_model6int, MMR_W_model6.POLY, MMR_W_model6intPOLY))



# AS / warm temps --------------
AS_W_model0 <- lmer(lnAS ~ lnBWg + tempTestK1000 + (1|species) , data=data.as.test, REML=FALSE)
AS_W_model0.POLY <- lmer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species), data=data.as.test, REML=FALSE)
AS_W_model0intPOLY <- lmer(lnAS ~ lnBWg * poly(tempTestK1000,2) + (1|species), data=data.as.test, REML=FALSE)
AS_W_model0int <- lmer(lnAS ~ lnBWg * tempTestK1000 + (1|species) , data=data.as.test, REML=FALSE)

AS_W_model1 <- lmer(lnAS ~ lnBWg + tempTestK1000 + (1|species) + (1|species:trial), data=data.as.test, REML=FALSE)
AS_W_model1.POLY <- lmer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.as.test, REML=FALSE)
AS_W_model1intPOLY <- lmer(lnAS ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.as.test, REML=FALSE)
AS_W_model1int <- lmer(lnAS ~ lnBWg * tempTestK1000 + (1|species) + (1|species:trial), data=data.as.test, REML=FALSE)

AS_W_model2 <- lmer(lnAS ~ lnBWg + tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.as.test, REML=FALSE)
AS_W_model2int <- lmer(lnAS ~ lnBWg * tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.as.test, REML=FALSE)
# ^^^ convergence issues: boundary (singular) fit: see help('isSingular')
AS_W_model2.POLY <- lmer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.as.test, REML=FALSE)
AS_W_model2intPOLY <- lmer(lnAS ~ lnBWg * poly(tempTestK1000,2)  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.as.test, REML=FALSE)

AS_W_model4 <- lmer(lnAS ~ lnBWg + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE)
AS_W_model4int <- lmer(lnAS ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE)
# ^^^ singular fit 
AS_W_model4.POLY <- lmer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE)
AS_W_model4intPOLY <- lmer(lnAS ~ lnBWg * poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE)
# ^^^ singular fit 

AS_W_model5 <- lmer(lnAS ~ lnBWg + tempTestK1000  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE)
AS_W_model5int <- lmer(lnAS ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE)
# ^^^ convergence issues: Model failed to converge with max|grad|
AS_W_model5.POLY <- lmer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE)
AS_W_model5intPOLY <- lmer(lnAS ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE)

AS_W_model6 <- lmer(lnAS ~ lnBWg + tempTestK1000 + (lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE)
# ^^^ singular fit 
AS_W_model6int <- lmer(lnAS ~ lnBWg * tempTestK1000+ (lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE)
# ^^^ singular fit 
AS_W_model6.POLY <- lmer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE)
# ^^^ singular fit 
AS_W_model6intPOLY <- lmer(lnAS ~ lnBWg * poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE)
# ^^^ singular fit 

## BIC -------------
AS_W_BIC.nP<-BICdelta(BIC(AS_W_model0, AS_W_model0int, AS_W_model0.POLY, AS_W_model0intPOLY,
                AS_W_model1, AS_W_model1int, AS_W_model1.POLY, AS_W_model1intPOLY,
                AS_W_model2, AS_W_model2int, AS_W_model2.POLY, AS_W_model2intPOLY,
                AS_W_model4, AS_W_model4int, AS_W_model4.POLY, AS_W_model4intPOLY,
                AS_W_model5, AS_W_model5int, AS_W_model5.POLY, AS_W_model5intPOLY,
                AS_W_model6, AS_W_model6int, AS_W_model6.POLY, AS_W_model6intPOLY))


# FAS / warm temps --------------
FAS_W_model0 <- lmer(lnFAS ~ lnBWg + tempTest + (1|species) , data=data.fas.test, REML=FALSE)
FAS_W_model0.POLY <- lmer(lnFAS ~ lnBWg + poly(tempTest,2) + (1|species), data=data.fas.test, REML=FALSE)
FAS_W_model0intPOLY <- lmer(lnFAS ~ lnBWg * poly(tempTest,2) + (1|species), data=data.fas.test, REML=FALSE)
FAS_W_model0int <- lmer(lnFAS ~ lnBWg * tempTest + (1|species) , data=data.fas.test, REML=FALSE)

FAS_W_model1 <- lmer(lnFAS ~ lnBWg + tempTest + (1|species) + (1|species:trial), data=data.fas.test, REML=FALSE)
FAS_W_model1.POLY <- lmer(lnFAS ~ lnBWg + poly(tempTest,2) + (1|species) + (1|species:trial), data=data.fas.test, REML=FALSE)
FAS_W_model1intPOLY <- lmer(lnFAS ~ lnBWg * poly(tempTest,2) + (1|species) + (1|species:trial), data=data.fas.test, REML=FALSE)
FAS_W_model1int <- lmer(lnFAS ~ lnBWg * tempTest + (1|species) + (1|species:trial), data=data.fas.test, REML=FALSE)

FAS_W_model2 <- lmer(lnFAS ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE)
FAS_W_model2int <- lmer(lnFAS ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE)
FAS_W_model2.POLY <- lmer(lnFAS ~ lnBWg + poly(tempTest,2) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE)
# ^^^ convergence issues: Model failed to converge with max|grad|
FAS_W_model2intPOLY <- lmer(lnFAS ~ lnBWg * poly(tempTest,2)  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE)

FAS_W_model4 <- lmer(lnFAS ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE)
# ^^^ singular fit
FAS_W_model4int <- lmer(lnFAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE)
FAS_W_model4.POLY <- lmer(lnFAS ~ lnBWg + poly(tempTest,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE)
# ^^^ singular fit 
FAS_W_model4intPOLY <- lmer(lnFAS ~ lnBWg * poly(tempTest,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE)

FAS_W_model5 <- lmer(lnFAS ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE)
# ^^^ singular fit 
FAS_W_model5int <- lmer(lnFAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE)
# ^^^ singular fit 
FAS_W_model5.POLY <- lmer(lnFAS ~ lnBWg + poly(tempTest,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE)
# ^^^ singular fit 
FAS_W_model5intPOLY <- lmer(lnFAS ~ lnBWg * poly(tempTest,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE)
# ^^^ singular fit 


FAS_W_model6 <- lmer(lnFAS ~ lnBWg + tempTestK1000 + (lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE)
# ^^^ singular fit 
FAS_W_model6int <- lmer(lnFAS ~ lnBWg * tempTestK1000+ (lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE)
FAS_W_model6.POLY <- lmer(lnFAS ~ lnBWg + poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE)
# ^^^ singular fit 
FAS_W_model6intPOLY <- lmer(lnFAS ~ lnBWg * poly(tempTestK1000,2) + (lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE)

## BIC -------------
FAS_W_BIC.nP<-BICdelta(BIC(FAS_W_model0, FAS_W_model0int, FAS_W_model0.POLY, FAS_W_model0intPOLY,
                FAS_W_model1, FAS_W_model1int, FAS_W_model1.POLY, FAS_W_model1intPOLY,
                FAS_W_model2, FAS_W_model2int, FAS_W_model2.POLY, FAS_W_model2intPOLY,
                FAS_W_model4, FAS_W_model4int, FAS_W_model4.POLY, FAS_W_model4intPOLY,
                FAS_W_model5, FAS_W_model5int, FAS_W_model5.POLY, FAS_W_model5intPOLY,
                FAS_W_model6, FAS_W_model6int, FAS_W_model6.POLY, FAS_W_model6intPOLY))

## BIC results ------
RMR_BIC.nP # RMR_model5 
MMR_BIC.nP # MMR_model4int
AS_BIC.nP  # AS_model4
FAS_BIC.nP # FAS_model5 # different????

RMR_W_BIC.nP # RMR_W_model1
MMR_W_BIC.nP # MMR_W_model4
AS_W_BIC.nP  # AS_W_model1 # different????
FAS_W_BIC.nP # FAS_W_model2.POLY


# Best models  --------
rmr_mod_nonPhyloER<-RMR_model5
amr_mod_nonPhyloER<-MMR_model4int
as_mod_nonPhyloER<-AS_model4
fas_mod_nonPhyloER<-FAS_model5 # different??? than Phylo?

rmr_mod_nonPhyloW<-RMR_W_model1
amr_mod_nonPhyloW<-MMR_W_model4
as_mod_nonPhyloW<-AS_W_model1 # different???? than Phylo?
fas_mod_nonPhyloW<-FAS_W_model2.POLY

# un-comment to see - or run R markdown html summary script
# # Ecol relev
# summary(rmr_mod_nonPhyloER)
# summary(amr_mod_nonPhyloER)
# summary(fas_mod_nonPhyloER)
# summary(as_mod_nonPhyloER)
# 
# # warm
# summary(rmr_mod_nonPhyloW)
# summary(amr_mod_nonPhyloW)
# summary(fas_mod_nonPhyloW)
# summary(as_mod_nonPhyloW)
# 
# # Ecol relev
# plot(rmr_mod_nonPhyloER)
# plot(amr_mod_nonPhyloER)
# plot(fas_mod_nonPhyloER)
# plot(as_mod_nonPhyloER)
# 
# # warm
# plot(rmr_mod_nonPhyloW)
# plot(amr_mod_nonPhyloW)
# plot(fas_mod_nonPhyloW)
# plot(as_mod_nonPhyloW)
# 
# hist(resid(amr_mod_nonPhyloER))
# hist(resid(rmr_mod_nonPhyloER))
# hist(resid(fas_mod_nonPhyloER))
# hist(resid(as_mod_nonPhyloER))
# 
# hist(resid(amr_mod_nonPhyloW))
# hist(resid(rmr_mod_nonPhyloW))
# hist(resid(fas_mod_nonPhyloW))
# hist(resid(as_mod_nonPhyloW)) # little skew not too bad, only 5 measurements, all reasonable biologically

data.as.test$resid<-resid(as_mod_nonPhyloW)
data.as.test[which(data.as.test$resid < -1.5),] # could be considered outlier data

model_outputs(phylo = FALSE, 
              best.model.rmr.er = rmr_mod_nonPhyloER,
              best.model.amr.er= amr_mod_nonPhyloER,
              best.model.as.er= as_mod_nonPhyloER,
              best.model.fas.er= fas_mod_nonPhyloER,
              best.model.rmr.w= rmr_mod_nonPhyloW,
              best.model.amr.w= amr_mod_nonPhyloW,
              best.model.as.w= as_mod_nonPhyloW,
              best.model.fas.w= fas_mod_nonPhyloW)



# Can Run this part seperately: Figures -------
cols.as<<-c("#265F73", "#007E66", "#00C5A3")
cols.fas<<-c("#395200", "#89A000", "yellow")
cols.rmr<<-c("#C70039", "#FF6D7C", "#FFA3AC")
cols.amr<<-c("#00749F","#00A8D6", "#9CE9FF")
cols<-c("#00749F","#C70039","#00A8D6","#FF6D7C", "#9CE9FF","#FFA3AC", "#00C5A3", "#265F73")# AMR -rmr- AMR dark - rmr dark - light - as - fas

set.seed(51423)
# MMR and AMR used interchangeably throughout 

# uncomment and run if the libraries were not installed at the top 
# library(ggplot2)
# library(ggpubr)
# library(cowplot)
# library(forcats)
# library(here)


  # General scaling plots ------
  AMRmodel_plot1<-ggplot(data=data.amrER, aes(x=lnBWg, y=lnAMR)) +
    geom_point(alpha=0.9,  size=2, pch=19, color="grey70")+
    geom_line(data=data.plotAMRint_ER,
              aes(y = model_predFE, x=lnBWg,  group=tempTestK1000_inC, color = tempTestK1000_inC),
              linewidth=0.3, lty=1,alpha=0.8, show.legend=FALSE) +
    geom_point(alpha=0.9,  size=2, pch=21, color="grey50",fill="grey70" )+
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
    annotate("text",  x = -5.5, y = 3.9, label = bquote(Optimal:~italic(b)[FAS] == .(FAS_slope)),size=5, hjust=0, family="Arial", color = "black")+
    annotate("text",  x = -5.5, y = 3.5, label = bquote(Warm:~italic(b)[FAS] == .(FAS_slope_w)),size=5, hjust=0, family="Arial", color = "#475500")+
    annotate("text", label = paste("n = ", nrow(data.fasER), sep=""),   x = -5.5, y = 3.1, size=3, hjust=0, family="Arial", color = "black")+
    annotate("text", label = paste("n = ", nrow(data.fas.test), sep=""),  x = -5.5, y = 2.85, size=3, hjust=0, family="Arial", color = "#475500")+
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
  ggsave(filename = paste("./Figures/Fig1_Scaling_aLL_", Sys.Date(), ".png", sep=""),
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
  ggformat(ASmodel_plot3.E_ALL, x_title= expression(Temperature^-1~(1000/K)), y_title=expression(italic(ln)*AS~(mg~O[2]~h^-1~g^-1)), print = T)
  
  Arh.plot<-cowplot:::plot_grid(AMRmodel_plot3.E_ALL, RMRmodel_plot3.E_ALL, ASmodel_plot3.E_ALL, 
            align = "hv",
            axis = "l",
            nrow = 1,
            ncol = 3,
            label_size = 17)
  ggsave(filename = paste("./Figures/Fig_ArrheniusFigMMR-RMR-AS_", Sys.Date(), ".png", sep=""),
         plot=Arh.plot, width = 13, height = 5, units = "in")
  
  
  # Both MMR and RMR together - AS punchline plots ---------
  # Overall all fish together:
  MRmodel_plot1<-ggplot() +
    geom_ribbon(data=data.plotRMR_ER, mapping = aes(y = model_predFE, x=lnBWg, ymin=CI_2.5, ymax=CI_97.5, group = tempTestK1000), linetype=2, alpha=0.1, fill = "grey30")+
    geom_ribbon(data.plotAMRint_ER, mapping = aes(y = model_predFE, x=lnBWg, ymin=CI_2.5, ymax=CI_97.5, group = tempTestK1000), linetype=2, alpha=0.1, fill = "grey30")+
    geom_line(data=data.plotAMRint_ER[round(data.plotAMRint_ER$tempTestK1000,2)==3.49,], mapping = aes(y = model_predFE, x=lnBWg,  group=tempTestK1000), color="black", size=0.7, lty=1, show.legend=FALSE) +
    geom_line(data=data.plotRMR_ER[round(data.plotRMR_ER$tempTestK1000,2)==3.49,], aes(y = model_predFE, x=lnBWg,  group=tempTestK1000, color= tempTestK1000), color="black", size=0.7, lty=1, show.legend=FALSE) +
    scale_y_continuous(limits = c(-7, 12), breaks = seq(-7,12, 2))+
    scale_x_continuous(limits = c(-7, 12), breaks = seq(-7,12, 2))
  ggformat(MRmodel_plot1, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*MR~(mg~O[2]~h^-1)), print = T)
  
  MRmodel_plotW<-ggplot() +
    geom_ribbon(data.plotAMR_warm, mapping = aes(y = model_predFE, x=lnBWg, ymin=CI_2.5, ymax=CI_97.5, group = tempTestK1000), linetype=2, alpha=0.1, fill = cols.amr[2])+
    geom_ribbon(data=data.plotRMR_warm, mapping = aes(y = model_predFE, x=lnBWg, ymin=CI_2.5, ymax=CI_97.5, group = tempTestK1000), linetype=2, alpha=0.1, fill = cols.rmr[2])+
    geom_line(data=data.plotAMR_warm[round(data.plotAMR_warm$tempTestK1000,2)==3.49,], mapping = aes(y = model_predFE, x=lnBWg,  group=tempTestK1000), color=cols.amr[1], size=0.7, lty=1, show.legend=FALSE) +
    geom_line(data=data.plotRMR_warm[round(data.plotRMR_warm$tempTestK1000,2)==3.49,], aes(y = model_predFE, x=lnBWg,  group=tempTestK1000, color= tempTestK1000), color=cols.rmr[1], size=0.7, lty=1, show.legend=FALSE) +
    scale_y_continuous(limits = c(-7, 12), breaks = seq(-7,12, 2))+
    scale_x_continuous(limits = c(-7, 12), breaks = seq(-7,12, 2))
  ggformat(MRmodel_plotW, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*MR~(mg~O[2]~h^-1)), pr = T)
  
  
  FinalMR<-cowplot:::plot_grid(MRmodel_plot1, MRmodel_plotW,
            align = "hv",
            axis = "l",
            nrow = 1,
            ncol = 2,
            label_size = 17)
  ggsave(filename = paste("./Figures/Fig_Final_", Sys.Date(), ".png", sep=""),
         plot=FinalMR, width = 7.5, height = 3.5, units = "in")
  
  
  # Overall all fish together:
  MRmodel_plot2<-ggplot(data=data.rmrER, aes(x=lnBWg, y=lnRMR)) +
    geom_line(data=data.plotRMR_ER[round(data.plotRMR_ER$tempTestK1000,2)==3.39,], aes(y = model_predFE, x=lnBWg,  group=tempTestK1000, color= tempTestK1000), color="black", size=0.7, lty=1, show.legend=FALSE) +
    geom_line(data=data.plotAMRint_ER[round(data.plotAMRint_ER$tempTestK1000,2)==3.39,], aes(y = model_predFE, x=lnBWg,  group=tempTestK1000, color= tempTestK1000), color="black", size=0.7, lty=1, show.legend=FALSE) +
    geom_line(data=data.plotRMR_warm[round(data.plotRMR_warm$tempTestK1000,2)==3.39,], aes(y = model_predFE, x=lnBWg,  group=tempTestK1000), color=cols.rmr[2], size=1, lty=1, show.legend=FALSE) +
    geom_line(data=data.plotAMR_warm[round(data.plotAMR_warm$tempTestK1000,2)==3.39,], aes(y = model_predFE, x=lnBWg,  group=tempTestK1000), color=cols.amr[2], size=1, lty=1, show.legend=FALSE) +
     ylim(x = -5, 12)+
     xlim(x = -5, 12)
  ggformat(MRmodel_plot2, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*MR~(mg~O[2]~h^-1)), print = F)
  
  ggsave(filename = paste("./Figures/Fig5-final_MMR-RMR_", Sys.Date(), ".png", sep=""),
         plot=MRmodel_plot2, width = 4, height = 4, units = "in")
  
  

# Boxplots: FAS, AS, MR, mass-independent models ----------
fas_boxplot<-ggplot(data=data.fas, aes(y=FAS, x = test_category, label=species, fill=test_category3))+
  geom_boxplot(show.legend = F)+
  scale_fill_manual(values=cols.fas)+
  scale_x_discrete(labels=c("acclim" = "Acclimated \n warm", "acute" = "Acute \n warm", "ecol_relev" = "Optimal"))
ggformat(fas_boxplot, y_title = expression(FAS~(MMR/RMR)), x_title = "", print = F)
fas_boxplot<-fas_boxplot+theme(legend.position = "none")

rmr_boxplot<-ggplot(data=data.rmr, aes(y=mass_specrmr, x = test_category, label=species, fill=test_category))+
  geom_boxplot(show.legend = F)+
  scale_fill_manual(values=cols.rmr)+
  scale_x_discrete(labels=c("acclim" = "Acclimated \n warm", "acute" = "Acute \n warm", "ecol_relev" = "Optimal"))
ggformat(rmr_boxplot, y_title = expression(RMR~(mgO[2]~g^-1~h^-1)), x_title = "", print = F)
rmr_boxplot<-rmr_boxplot+theme(legend.position = "none")

amr_boxplot<-ggplot(data=data.amr, aes(y=mass_specamr, x = test_category, label=species, fill=test_category))+
  geom_boxplot(show.legend = F)+
  scale_fill_manual(values=cols.amr)+
  scale_x_discrete(labels=c("acclim" = "Acclimated \n warm", "acute" = "Acute \n warm", "ecol_relev" = "Optimal"))
ggformat(amr_boxplot, y_title = expression(MMR~(mgO[2]~g^-1~h^-1)), x_title = "", print = F)
amr_boxplot<-amr_boxplot+theme(legend.position = "none")

cowplot:::plot_grid(amr_boxplot,rmr_boxplot,fas_boxplot,
          align = "hv",
          axis = "l",
          nrow = 1,
          ncol = 3) %>%
ggsave(filename = "./Figures/FigSUP_boxplots_mar2022.png", width = 12.5, height =4)



# COMBO warm fish and ecol relev ********************
rmr_boxplot<-ggplot(data=data.rmr, aes(y=mass_specrmr, x = test_category3, label=species, fill=test_category3))+
  geom_boxplot(show.legend = F)+
  scale_fill_manual(values=c("grey50", cols.rmr[3]))+
  scale_x_discrete(labels=c("warm" = "Warm", "ecol_relev" = "Optimal"))
ggformat(rmr_boxplot, y_title = expression(RMR~(mgO[2]~g^-1~h^-1)), x_title = "", print = F)
rmr_boxplot<-rmr_boxplot+theme(legend.position = "none")
ggsave(filename = paste("./Figures/Fig_box_rmr1_", Sys.Date(), ".png", sep=""),
       plot=rmr_boxplot, width = 4, height = 4, units = "in")

amr_boxplot<-ggplot(data=data.amr, aes(y=mass_specamr, x = test_category3, label=species, fill=test_category3))+
  geom_boxplot(show.legend = F)+
  scale_fill_manual(values=c("grey50", cols.amr[3]))+
  scale_x_discrete(labels=c("warm" = "Warm", "ecol_relev" = "Optimal"))
ggformat(amr_boxplot, y_title = expression(MMR~(mgO[2]~g^-1~h^-1)), x_title = "", print = F)
amr_boxplot<-amr_boxplot+theme(legend.position = "none")
ggsave(filename = paste("./Figures/Fig_box_amr1_", Sys.Date(), ".png", sep=""),
       plot=amr_boxplot, width = 4, height = 4, units = "in")

fas_boxplot<-ggplot(data=data.fas, aes(y=FAS, x = test_category3, label=species, fill=test_category3))+
  geom_boxplot(show.legend = FALSE)+
  scale_fill_manual(values=c("grey50", cols.fas[3]))+
  scale_color_manual(values=cols.fas)+
  scale_x_discrete(labels=c("warm" = "Warm", "ecol_relev" = "Optimal"))
ggformat(fas_boxplot, y_title = "FAS (MMR / RMR)", x_title = "", print = F)
fas_boxplot<-fas_boxplot+theme(legend.position = "none")
ggsave(filename = paste("./Figures/Fig_box_fas1_", Sys.Date(), ".png", sep=""),
       plot=fas_boxplot, width = 4, height = 4, units = "in")




