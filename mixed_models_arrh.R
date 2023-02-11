

# libraries 
library(emmeans)
library(cowplot)
library(ggplot2)
library(mgcv)
library(myTAI)

# source the data --------
source("/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Codes/MR-fish-metadata/deltaIC_BICdelta_ggformat.R")

source("/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Codes/Metadata-temperature/JEB/get_scaling_data_temp.R")

data.list<-get_scaling_data_temp(data.amr = "/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Data/MR-fish-metadata-data/Fish_AMR_temp_dataset_mar2022.csv",
                                 data.rmr = "/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Data/MR-fish-metadata-data/Fish_RMR_temp_dataset_mar2022.csv",
                                 ecology.data = "/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Data/MR-fish-metadata-data/Kraskura_species_ecologies_mar2022.csv", 
                                 onlyTop.above = TRUE, 
                                 exp_rmr = 0.81, exp_amr_warm = 0.81, 
                                 exp_rmr_warm = 0.84, exp_amr = 0.82, 
                                 exp_as = 0.82, exp_as_warm = 0.78)

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


# Data where FAS is greater than 12. why?
data.fas12<-data.fas[(data.fas$FAS>12),]
# write.csv(file= "/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Data/MR-fish-metadata-data/dataFAS_largerThan12_mar22-2022_ALLtempCATEG.csv", data.fas12, row.names=FALSE)

unique(data.fas12$study_ID)
mean(data.fas$FAS[data.fas$study_ID==4])
mean(data.fas$FAS[data.fas$study_ID==6])
mean(data.fas$FAS[data.fas$study_ID==15])
mean(data.fas$FAS[data.fas$study_ID==20])
mean(data.fas$FAS[data.fas$study_ID==61])
mean(data.fas$FAS[data.fas$study_ID==410])
mean(data.fas$FAS[data.fas$study_ID==405])

length(data.fas$FAS[data.fas$study_ID==4])
length(data.fas$FAS[data.fas$study_ID==6])
length(data.fas$FAS[data.fas$study_ID==15])
length(data.fas$FAS[data.fas$study_ID==20])
length(data.fas$FAS[data.fas$study_ID==61])
length(data.fas$FAS[data.fas$study_ID==410])
length(data.fas$FAS[data.fas$study_ID==405])

k<-(8.62*10^(-5)) # Boltzmann's constant - units eV/K
# Boltz<-1.38 * 10^23; units J/K
E<-0.63 # activation energy MTE
# cols <- c("#d55e00", "#0072b2", "#009e73")

cols.as<-c("#265F73", "#00929A", "#00C5A3")
cols.rmr<-c("#C70039", "#FF6D7C", "#FFA3AC")
cols.amr<-c("#00749F","#00A8D6", "#9CE9FF")
cols.fas<-c("#C94F00", "#F37121", "#FFD784")

# the warm ones 
data.amr.test<-rbind(data.amrAC, data.amrAM)
data.rmr.test<-rbind(data.rmrAC, data.rmrAM)
data.fas.test<-data.amr.test[c(!is.na(data.amr.test$FAS) & is.finite(data.amr.test$FAS)) , ]
data.as.test<-data.amr.test[c(!is.na(data.amr.test$lnAS) & is.finite(data.amr.test$lnAS)) , ]

nrow(data.amr.test)+nrow(data.amrER)
nrow(data.rmr.test)+nrow(data.rmrER)

dataMR$test_category2<-"acclim"
dataMR$test_category2[dataMR$test_category=="acute"] <- "acute"

dataMR$test_category3<-"ecol_relev"
dataMR$test_category3[dataMR$test_category=="acute"] <- "warm"
dataMR$test_category3[dataMR$test_category=="acclim"] <- "warm"


# *****************
# Multiple regression model (Downs et al)---------
# ln(MR) = ln(a) + b*ln(BW) + c(1000/T) 
# c = -1E/1000k 

# all out-commented models are singular fits

# RMR / ecologically relevant temps --------------
#
#           RMR
# **************************
RMR_MultReg_ER_model0 <- lmer(lnRMR ~ lnBWg + tempTestK1000 + (1|species) , data=data.rmrER, REML=FALSE)
RMR_MultReg_ER_model0.POLY <- lmer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species), data=data.rmrER, REML=FALSE)
RMR_MultReg_ER_model0intPOLY <- lmer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (1|species), data=data.rmrER, REML=FALSE)
RMR_MultReg_ER_model0int <- lmer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) , data=data.rmrER, REML=FALSE)

RMR_MultReg_ER_model1 <- lmer(lnRMR ~ lnBWg + tempTestK1000 + (1|species) + (1|species:trial), data=data.rmrER, REML=FALSE)
RMR_MultReg_ER_model1.POLY <- lmer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.rmrER, REML=FALSE)
RMR_MultReg_ER_model1intPOLY <- lmer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.rmrER, REML=FALSE)
RMR_MultReg_ER_model1int <- lmer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) + (1|species:trial), data=data.rmrER, REML=FALSE)

RMR_MultReg_ER_model2 <- lmer(lnRMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.rmrER, REML=FALSE)
RMR_MultReg_ER_model2int <- lmer(lnRMR ~ lnBWg * tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.rmrER, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# RMR_MultReg_ER_model2.POLY <- lmer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.rmrER, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
RMR_MultReg_ER_model2intPOLY <- lmer(lnRMR ~ lnBWg * poly(tempTestK1000,2)  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.rmrER, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

# not enough data for this - take out from all models;  singular fits 
# RMR_MultReg_ER_model3 <- lmer(lnRMR ~ lnBWg + tempTestK1000 + (1|species) +(0 + tempTestK1000|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# RMR_MultReg_ER_model3int <- lmer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + tempTestK1000|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# RMR_MultReg_ER_model3.POLY <- lmer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (0 + tempTestK1000|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# RMR_MultReg_ER_model3intPOLY <- lmer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (0 + tempTestK1000|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

RMR_MultReg_ER_model4 <- lmer(lnRMR ~ lnBWg + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE) 
RMR_MultReg_ER_model4int <- lmer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE)
RMR_MultReg_ER_model4.POLY <- lmer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE)
RMR_MultReg_ER_model4intPOLY <- lmer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE)

RMR_MultReg_ER_model5 <- lmer(lnRMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE)
RMR_MultReg_ER_model5int <- lmer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE)
RMR_MultReg_ER_model5.POLY <- lmer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
RMR_MultReg_ER_model5intPOLY <- lmer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))


# AMR / ecologically relevant temps --------------
#
#           AMR
# **************************
AMR_MultReg_ER_model0 <- lmer(lnAMR ~ lnBWg + tempTestK1000 + (1|species) , data=data.amrER, REML=FALSE)
AMR_MultReg_ER_model0.POLY <- lmer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species), data=data.amrER, REML=FALSE)
AMR_MultReg_ER_model0intPOLY <- lmer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (1|species), data=data.amrER, REML=FALSE)
AMR_MultReg_ER_model0int <- lmer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) , data=data.amrER, REML=FALSE)

AMR_MultReg_ER_model1 <- lmer(lnAMR ~ lnBWg + tempTestK1000 + (1|species) + (1|species:trial), data=data.amrER, REML=FALSE)
AMR_MultReg_ER_model1.POLY <- lmer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.amrER, REML=FALSE)
AMR_MultReg_ER_model1intPOLY <- lmer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.amrER, REML=FALSE)
AMR_MultReg_ER_model1int <- lmer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) + (1|species:trial), data=data.amrER, REML=FALSE)

AMR_MultReg_ER_model2 <- lmer(lnAMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.amrER, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
AMR_MultReg_ER_model2int <- lmer(lnAMR ~ lnBWg * tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.amrER, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
AMR_MultReg_ER_model2.POLY <- lmer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.amrER, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
AMR_MultReg_ER_model2intPOLY <- lmer(lnAMR ~ lnBWg * poly(tempTestK1000,2)  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.amrER, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

# AMR_MultReg_ER_model3 <- lmer(lnAMR ~ lnBWg + tempTestK1000 + (1|species) +(0 + tempTestK1000|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# AMR_MultReg_ER_model3int <- lmer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + tempTestK1000|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# AMR_MultReg_ER_model3.POLY <- lmer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (0 + tempTestK1000|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# AMR_MultReg_ER_model3intPOLY <- lmer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (0 + tempTestK1000|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

AMR_MultReg_ER_model4 <- lmer(lnAMR ~ lnBWg + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE) 
AMR_MultReg_ER_model4int <- lmer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE)
AMR_MultReg_ER_model4.POLY <- lmer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE)
AMR_MultReg_ER_model4intPOLY <- lmer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE)

AMR_MultReg_ER_model5 <- lmer(lnAMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE)
AMR_MultReg_ER_model5int <- lmer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE)
AMR_MultReg_ER_model5.POLY <- lmer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
AMR_MultReg_ER_model5intPOLY <- lmer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

# AS / ecologically relevant temps --------------
#
#           AS
# **************************
AS_MultReg_ER_model0 <- lmer(lnAS ~ lnBWg + tempTestK1000 + (1|species) , data=data.asER, REML=FALSE)
AS_MultReg_ER_model0.POLY <- lmer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species), data=data.asER, REML=FALSE)
AS_MultReg_ER_model0intPOLY <- lmer(lnAS ~ lnBWg * poly(tempTestK1000,2) + (1|species), data=data.asER, REML=FALSE)
AS_MultReg_ER_model0int <- lmer(lnAS ~ lnBWg * tempTestK1000 + (1|species) , data=data.asER, REML=FALSE)

AS_MultReg_ER_model1 <- lmer(lnAS ~ lnBWg + tempTestK1000 + (1|species) + (1|species:trial), data=data.asER, REML=FALSE)
AS_MultReg_ER_model1.POLY <- lmer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.asER, REML=FALSE)
AS_MultReg_ER_model1intPOLY <- lmer(lnAS ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.asER, REML=FALSE)
AS_MultReg_ER_model1int <- lmer(lnAS ~ lnBWg * tempTestK1000 + (1|species) + (1|species:trial), data=data.asER, REML=FALSE)

AS_MultReg_ER_model2 <- lmer(lnAS ~ lnBWg + tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.asER, REML=FALSE)
AS_MultReg_ER_model2int <- lmer(lnAS ~ lnBWg * tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.asER, REML=FALSE)
AS_MultReg_ER_model2.POLY <- lmer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.asER, REML=FALSE)
AS_MultReg_ER_model2intPOLY <- lmer(lnAS ~ lnBWg * poly(tempTestK1000,2)  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.asER, REML=FALSE)

AS_MultReg_ER_model4 <- lmer(lnAS ~ lnBWg + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE)
AS_MultReg_ER_model4int <- lmer(lnAS ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE)
AS_MultReg_ER_model4.POLY <- lmer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE)
AS_MultReg_ER_model4intPOLY <- lmer(lnAS ~ lnBWg * poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE)

AS_MultReg_ER_model5 <- lmer(lnAS ~ lnBWg + tempTestK1000  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE)
AS_MultReg_ER_model5int <- lmer(lnAS ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE)
AS_MultReg_ER_model5.POLY <- lmer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
AS_MultReg_ER_model5intPOLY <- lmer(lnAS ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))


# # FAS / ecologically relevant temps --------------
# #
# #           lnFAS
# # **************************
FAS_MultReg_ER_model0 <- lmer(log(FAS) ~ lnBWg + tempTest + (1|species) , data=data.fasER, REML=FALSE)
FAS_MultReg_ER_model0.POLY <- lmer(log(FAS) ~ lnBWg + poly(tempTest,2) + (1|species), data=data.fasER, REML=FALSE)
FAS_MultReg_ER_model0intPOLY <- lmer(log(FAS) ~ lnBWg * poly(tempTest,2) + (1|species), data=data.fasER, REML=FALSE)
FAS_MultReg_ER_model0int <- lmer(log(FAS) ~ lnBWg * tempTest + (1|species) , data=data.fasER, REML=FALSE)

FAS_MultReg_ER_model1 <- lmer(log(FAS) ~ lnBWg + tempTest + (1|species) + (1|species:trial), data=data.fasER, REML=FALSE)
FAS_MultReg_ER_model1.POLY <- lmer(log(FAS) ~ lnBWg + poly(tempTest,2) + (1|species) + (1|species:trial), data=data.fasER, REML=FALSE)
FAS_MultReg_ER_model1intPOLY <- lmer(log(FAS) ~ lnBWg * poly(tempTest,2) + (1|species) + (1|species:trial), data=data.fasER, REML=FALSE)
FAS_MultReg_ER_model1int <- lmer(log(FAS) ~ lnBWg * tempTest + (1|species) + (1|species:trial), data=data.fasER, REML=FALSE)

FAS_MultReg_ER_model2 <- lmer(log(FAS) ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE,control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
FAS_MultReg_ER_model2int <- lmer(log(FAS) ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE)
FAS_MultReg_ER_model2.POLY <- lmer(log(FAS) ~ lnBWg + poly(tempTest,2) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE,control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)) )
FAS_MultReg_ER_model2intPOLY <- lmer(log(FAS) ~ lnBWg * poly(tempTest,2)  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE)

FAS_MultReg_ER_model4 <- lmer(log(FAS) ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE)
FAS_MultReg_ER_model4int <- lmer(log(FAS) ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE)
FAS_MultReg_ER_model4.POLY <- lmer(log(FAS) ~ lnBWg + poly(tempTest,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE)
FAS_MultReg_ER_model4intPOLY <- lmer(log(FAS) ~ lnBWg * poly(tempTest,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE)

FAS_MultReg_ER_model5 <- lmer(log(FAS) ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE)
FAS_MultReg_ER_model5int <- lmer(log(FAS) ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE)
FAS_MultReg_ER_model5.POLY <- lmer(log(FAS) ~ lnBWg + poly(tempTest,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
FAS_MultReg_ER_model5intPOLY <- lmer(log(FAS) ~ lnBWg * poly(tempTest,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))



# Model compariso: BIC ------------------

RMR_ER_BIC<-BIC(RMR_MultReg_ER_model0, RMR_MultReg_ER_model0int, RMR_MultReg_ER_model0.POLY, RMR_MultReg_ER_model0intPOLY,
                RMR_MultReg_ER_model1, RMR_MultReg_ER_model1int, RMR_MultReg_ER_model1.POLY, RMR_MultReg_ER_model1intPOLY,
                RMR_MultReg_ER_model2, RMR_MultReg_ER_model2int, RMR_MultReg_ER_model2intPOLY,
                # RMR_MultReg_ER_model3, RMR_MultReg_ER_model3int, RMR_MultReg_ER_model3.POLY, AMR_MultReg_ER_model3intPOLY, RMR_MultReg_ER_model2.POLY, # << singular fits, exclude here 
                RMR_MultReg_ER_model4, RMR_MultReg_ER_model4int, RMR_MultReg_ER_model4.POLY, RMR_MultReg_ER_model4intPOLY,
                RMR_MultReg_ER_model5, RMR_MultReg_ER_model5int, RMR_MultReg_ER_model5.POLY, RMR_MultReg_ER_model5intPOLY)


AMR_ER_BIC<-BIC(AMR_MultReg_ER_model0, AMR_MultReg_ER_model0int, AMR_MultReg_ER_model0.POLY, AMR_MultReg_ER_model0intPOLY,
                AMR_MultReg_ER_model1, AMR_MultReg_ER_model1int, AMR_MultReg_ER_model1.POLY, AMR_MultReg_ER_model1intPOLY,
                AMR_MultReg_ER_model2, AMR_MultReg_ER_model2int, AMR_MultReg_ER_model2.POLY, AMR_MultReg_ER_model2intPOLY,
                # AMR_MultReg_ER_model3, AMR_MultReg_ER_model3int, AMR_MultReg_ER_model3.POLY, AMR_MultReg_ER_model3intPOLY, # singular fits 
                AMR_MultReg_ER_model4, AMR_MultReg_ER_model4int, AMR_MultReg_ER_model4.POLY, AMR_MultReg_ER_model4intPOLY,
                AMR_MultReg_ER_model5, AMR_MultReg_ER_model5int, AMR_MultReg_ER_model5.POLY, AMR_MultReg_ER_model5intPOLY)

FAS_ER_BIC <- BIC(FAS_MultReg_ER_model0, FAS_MultReg_ER_model0.POLY, FAS_MultReg_ER_model0int, FAS_MultReg_ER_model0intPOLY,
                  FAS_MultReg_ER_model1, FAS_MultReg_ER_model1.POLY, FAS_MultReg_ER_model1int, FAS_MultReg_ER_model1intPOLY,
                  FAS_MultReg_ER_model2, FAS_MultReg_ER_model2.POLY, FAS_MultReg_ER_model2int, FAS_MultReg_ER_model2intPOLY,
                  FAS_MultReg_ER_model4, FAS_MultReg_ER_model4.POLY, FAS_MultReg_ER_model4int, FAS_MultReg_ER_model4intPOLY,
                  FAS_MultReg_ER_model5, FAS_MultReg_ER_model5.POLY, FAS_MultReg_ER_model5int, FAS_MultReg_ER_model5intPOLY)

AS_ER_BIC <- BIC(AS_MultReg_ER_model0, AS_MultReg_ER_model0.POLY, AS_MultReg_ER_model0int, AS_MultReg_ER_model0intPOLY,
                 AS_MultReg_ER_model1, AS_MultReg_ER_model1.POLY, AS_MultReg_ER_model1int, AS_MultReg_ER_model1intPOLY,
                 AS_MultReg_ER_model2, AS_MultReg_ER_model2.POLY,  AS_MultReg_ER_model2intPOLY,AS_MultReg_ER_model2int,
                 AS_MultReg_ER_model4, AS_MultReg_ER_model4.POLY, AS_MultReg_ER_model4int, AS_MultReg_ER_model4intPOLY,
                 AS_MultReg_ER_model5, AS_MultReg_ER_model5.POLY, AS_MultReg_ER_model5int, AS_MultReg_ER_model5intPOLY)

BICdelta(RMR_ER_BIC) # RMR_MultReg_ER_model5
BICdelta(AMR_ER_BIC) # AMR_MultReg_ER_model4int
BICdelta(FAS_ER_BIC) # FAS_MultReg_ER_model4
BICdelta(AS_ER_BIC) # AS_MultReg_ER_model4


# RMR / warm temps --------------
RMR_MultReg_W_model0 <- lmer(lnRMR ~ lnBWg + tempTestK1000 + (1|species) , data=data.rmr.test, REML=FALSE)
RMR_MultReg_W_model0.POLY <- lmer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species), data=data.rmr.test, REML=FALSE)
RMR_MultReg_W_model0intPOLY <- lmer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (1|species), data=data.rmr.test, REML=FALSE)
RMR_MultReg_W_model0int <- lmer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) , data=data.rmr.test, REML=FALSE)

RMR_MultReg_W_model1 <- lmer(lnRMR ~ lnBWg + tempTestK1000 + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)
RMR_MultReg_W_model1.POLY <- lmer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)
RMR_MultReg_W_model1intPOLY <- lmer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)
RMR_MultReg_W_model1int <- lmer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)

RMR_MultReg_W_model2 <- lmer(lnRMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.rmr.test, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
RMR_MultReg_W_model2int <- lmer(lnRMR ~ lnBWg * tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# RMR_MultReg_W_model2.POLY <- lmer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.rmr.test, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# RMR_MultReg_W_model2intPOLY <- lmer(lnRMR ~ lnBWg * poly(tempTestK1000,2)  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.rmr.test, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

RMR_MultReg_W_model4 <- lmer(lnRMR ~ lnBWg + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE) 
# RMR_MultReg_W_model4int <- lmer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)
RMR_MultReg_W_model4.POLY <- lmer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)
# RMR_MultReg_W_model4intPOLY <- lmer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)

# RMR_MultReg_W_model5 <- lmer(lnRMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE)
# RMR_MultReg_W_model5int <- lmer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE)
# RMR_MultReg_W_model5.POLY <- lmer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# RMR_MultReg_W_model5intPOLY <- lmer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))


# AMR / warm temps --------------

AMR_MultReg_W_model0 <- lmer(lnAMR ~ lnBWg + tempTestK1000 + (1|species) , data=data.amr.test, REML=FALSE)
AMR_MultReg_W_model0.POLY <- lmer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species), data=data.amr.test, REML=FALSE)
AMR_MultReg_W_model0intPOLY <- lmer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (1|species), data=data.amr.test, REML=FALSE)
AMR_MultReg_W_model0int <- lmer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) , data=data.amr.test, REML=FALSE)

AMR_MultReg_W_model1 <- lmer(lnAMR ~ lnBWg + tempTestK1000 + (1|species) + (1|species:trial), data=data.amr.test, REML=FALSE)
AMR_MultReg_W_model1.POLY <- lmer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.amr.test, REML=FALSE)
AMR_MultReg_W_model1intPOLY <- lmer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.amr.test, REML=FALSE)
AMR_MultReg_W_model1int <- lmer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) + (1|species:trial), data=data.amr.test, REML=FALSE)

# AMR_MultReg_W_model2 <- lmer(lnAMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.amr.test, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
AMR_MultReg_W_model2int <- lmer(lnAMR ~ lnBWg * tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.amr.test, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# AMR_MultReg_W_model2.POLY <- lmer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.amr.test, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
AMR_MultReg_W_model2intPOLY <- lmer(lnAMR ~ lnBWg * poly(tempTestK1000,2)  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.amr.test, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

AMR_MultReg_W_model4 <- lmer(lnAMR ~ lnBWg + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE) 
AMR_MultReg_W_model4int <- lmer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE)
AMR_MultReg_W_model4.POLY <- lmer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE)
AMR_MultReg_W_model4intPOLY <- lmer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE)

AMR_MultReg_W_model5 <- lmer(lnAMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE)
AMR_MultReg_W_model5int <- lmer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE)
AMR_MultReg_W_model5.POLY <- lmer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
AMR_MultReg_W_model5intPOLY <- lmer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))


# FAS / warm temps --------------

FAS_MultReg_W_model0 <- lmer(lnFAS ~ lnBWg + tempTest + (1|species) , data=data.fas.test, REML=FALSE)
FAS_MultReg_W_model0.POLY <- lmer(lnFAS ~ lnBWg + poly(tempTest,2) + (1|species), data=data.fas.test, REML=FALSE)
FAS_MultReg_W_model0intPOLY <- lmer(lnFAS ~ lnBWg * poly(tempTest,2) + (1|species), data=data.fas.test, REML=FALSE)
FAS_MultReg_W_model0int <- lmer(lnFAS ~ lnBWg * tempTest + (1|species) , data=data.fas.test, REML=FALSE)

FAS_MultReg_W_model1 <- lmer(lnFAS ~ lnBWg + tempTest + (1|species) + (1|species:trial), data=data.fas.test, REML=FALSE)
FAS_MultReg_W_model1.POLY <- lmer(lnFAS ~ lnBWg + poly(tempTest,2) + (1|species) + (1|species:trial), data=data.fas.test, REML=FALSE)
FAS_MultReg_W_model1intPOLY <- lmer(lnFAS ~ lnBWg * poly(tempTest,2) + (1|species) + (1|species:trial), data=data.fas.test, REML=FALSE)
FAS_MultReg_W_model1int <- lmer(lnFAS ~ lnBWg * tempTest + (1|species) + (1|species:trial), data=data.fas.test, REML=FALSE)

FAS_MultReg_W_model2 <- lmer(lnFAS ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
FAS_MultReg_W_model2int <- lmer(lnFAS ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
FAS_MultReg_W_model2.POLY <- lmer(lnFAS ~ lnBWg + poly(tempTest,2) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
FAS_MultReg_W_model2intPOLY <- lmer(lnFAS ~ lnBWg * poly(tempTest,2)  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

FAS_MultReg_W_model4 <- lmer(lnFAS ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 
FAS_MultReg_W_model4int <- lmer(lnFAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE)
FAS_MultReg_W_model4.POLY <- lmer(lnFAS ~ lnBWg + poly(tempTest,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
FAS_MultReg_W_model4intPOLY <- lmer(lnFAS ~ lnBWg * poly(tempTest,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE)

# FAS_MultReg_W_model5 <- lmer(lnFAS ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# FAS_MultReg_W_model5int <- lmer(lnFAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# FAS_MultReg_W_model5.POLY <- lmer(lnFAS ~ lnBWg + poly(tempTest,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# FAS_MultReg_W_model5intPOLY <- lmer(lnFAS ~ lnBWg * poly(tempTest,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))


# AS / warm temps --------------

AS_MultReg_W_model0 <- lmer(lnAS ~ lnBWg + tempTestK1000 + (1|species) , data=data.as.test, REML=FALSE)
AS_MultReg_W_model0.POLY <- lmer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species), data=data.as.test, REML=FALSE)
AS_MultReg_W_model0intPOLY <- lmer(lnAS ~ lnBWg * poly(tempTestK1000,2) + (1|species), data=data.as.test, REML=FALSE)
AS_MultReg_W_model0int <- lmer(lnAS ~ lnBWg * tempTestK1000 + (1|species) , data=data.as.test, REML=FALSE)

AS_MultReg_W_model1 <- lmer(lnAS ~ lnBWg + tempTestK1000 + (1|species) + (1|species:trial), data=data.as.test, REML=FALSE)
AS_MultReg_W_model1.POLY <- lmer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.as.test, REML=FALSE)
AS_MultReg_W_model1intPOLY <- lmer(lnAS ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.as.test, REML=FALSE)
AS_MultReg_W_model1int <- lmer(lnAS ~ lnBWg * tempTestK1000 + (1|species) + (1|species:trial), data=data.as.test, REML=FALSE)

AS_MultReg_W_model2 <- lmer(lnAS ~ lnBWg + tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.as.test, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# AS_MultReg_W_model2int <- lmer(lnAS ~ lnBWg * tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.as.test, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
AS_MultReg_W_model2.POLY <- lmer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.as.test, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
AS_MultReg_W_model2intPOLY <- lmer(lnAS ~ lnBWg * poly(tempTestK1000,2)  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.as.test, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

AS_MultReg_W_model4 <- lmer(lnAS ~ lnBWg + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 
# AS_MultReg_W_model4int <- lmer(lnAS ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE)
AS_MultReg_W_model4.POLY <- lmer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# AS_MultReg_W_model4intPOLY <- lmer(lnAS ~ lnBWg * poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE)

AS_MultReg_W_model5 <- lmer(lnAS ~ lnBWg + tempTestK1000  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
AS_MultReg_W_model5int <- lmer(lnAS ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
AS_MultReg_W_model5.POLY <- lmer(lnAS ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
AS_MultReg_W_model5intPOLY <- lmer(lnAS ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))


# BIC warm fish model comparison: --------------
BICdelta(BIC(RMR_MultReg_W_model0, RMR_MultReg_W_model0.POLY, RMR_MultReg_W_model0intPOLY,RMR_MultReg_W_model0int, 
             RMR_MultReg_W_model1, RMR_MultReg_W_model1.POLY, RMR_MultReg_W_model1intPOLY,RMR_MultReg_W_model1int,
             RMR_MultReg_W_model2, RMR_MultReg_W_model2int,
             RMR_MultReg_W_model4, RMR_MultReg_W_model4.POLY
             # RMR_MultReg_W_model5, RMR_MultReg_W_model5.POLY, RMR_MultReg_W_model5intPOLY,RMR_MultReg_W_model5int, # singular fits 
             # RMR_MultReg_W_model4intPOLY,RMR_MultReg_W_model4int, RMR_MultReg_W_model2.POLY, RMR_MultReg_W_model2intPOLY # singular fits 
))

BICdelta(BIC(AMR_MultReg_W_model0, AMR_MultReg_W_model0.POLY, AMR_MultReg_W_model0intPOLY,AMR_MultReg_W_model0int, 
             AMR_MultReg_W_model1, AMR_MultReg_W_model1.POLY, AMR_MultReg_W_model1intPOLY,AMR_MultReg_W_model1int,
             AMR_MultReg_W_model2intPOLY,AMR_MultReg_W_model2int,
             # AMR_MultReg_W_model2,AMR_MultReg_W_model2.POLY, # singular fits 
             AMR_MultReg_W_model4, AMR_MultReg_W_model4.POLY, AMR_MultReg_W_model4intPOLY,AMR_MultReg_W_model4int,
             AMR_MultReg_W_model5, AMR_MultReg_W_model5.POLY, AMR_MultReg_W_model5intPOLY,AMR_MultReg_W_model5int))

BICdelta(BIC(FAS_MultReg_W_model0, FAS_MultReg_W_model0.POLY, FAS_MultReg_W_model0intPOLY,FAS_MultReg_W_model0int, 
             FAS_MultReg_W_model1, FAS_MultReg_W_model1.POLY, FAS_MultReg_W_model1intPOLY,FAS_MultReg_W_model1int,
             FAS_MultReg_W_model4, FAS_MultReg_W_model4.POLY, FAS_MultReg_W_model4intPOLY,FAS_MultReg_W_model4int,
             FAS_MultReg_W_model2, FAS_MultReg_W_model2.POLY, FAS_MultReg_W_model2intPOLY,FAS_MultReg_W_model2int)) # all set 5 = singular fits 

BICdelta(BIC(AS_MultReg_W_model0, AS_MultReg_W_model0.POLY, AS_MultReg_W_model0intPOLY,AS_MultReg_W_model0int, 
             AS_MultReg_W_model1, AS_MultReg_W_model1.POLY, AS_MultReg_W_model1intPOLY,AS_MultReg_W_model1int,
             AS_MultReg_W_model4, AS_MultReg_W_model4.POLY,
             AS_MultReg_W_model2, AS_MultReg_W_model2.POLY, AS_MultReg_W_model2intPOLY,
             # AS_MultReg_W_model2int, AS_MultReg_W_model4int, AS_MultReg_W_model4intPOLY, # singular fits
             AS_MultReg_W_model5, AS_MultReg_W_model5.POLY, AS_MultReg_W_model5intPOLY,AS_MultReg_W_model5int))

## Best model summaries, CI, residual plots ------------
# changes updates march 22, 2022
RMR_MultReg_ER_model5 <- lmer(lnRMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE)
AMR_MultReg_ER_model4int <- lmer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE)
FAS_MultReg_ER_model4 <- lmer(log(FAS) ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE)
AS_MultReg_ER_model4 <- lmer(lnAS ~ lnBWg + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE)

RMR_MultReg_W_model1 <- lmer(lnRMR ~ lnBWg + tempTestK1000 + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)
AMR_MultReg_W_model4 <- lmer(lnAMR ~ lnBWg + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE) 
FAS_MultReg_W_model2.POLY <- lmer(lnFAS ~ lnBWg + poly(tempTest,2) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
AS_MultReg_W_model1 <- lmer(lnAS ~ lnBWg + tempTestK1000 + (1|species) + (1|species:trial), data=data.as.test, REML=FALSE)

rmr_mod_ER<-RMR_MultReg_ER_model5
amr_mod_ER<-AMR_MultReg_ER_model4int
fas_mod_ER<-FAS_MultReg_ER_model4
as_mod_ER<-AS_MultReg_ER_model4

rmr_mod_W<-RMR_MultReg_W_model1
amr_mod_W<-AMR_MultReg_W_model4
fas_mod_W<-FAS_MultReg_W_model2.POLY
as_mod_W<-AS_MultReg_W_model1

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
data.as.test[which(data.as.test$resid < -1.5),]

# CIs
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

sum_CItable<-rbind(CI.amr.ER, CI.rmr.ER, CI.fas.ER, CI.as.ER, 
                   CI.amr.W, CI.rmr.W, CI.fas.W, CI.as.W)
write.csv(file = "/Users/kristakraskura/Desktop/BOX/UCSB/Research/Metabolic_scaling/ms-AMR-RMR-Temperature/thesis/data files/tableX_CIsummary.csv", sum_CItable, row.names=FALSE)





# NOT USED models ------

# 1. AMR and RMR: full datasets, no temp test categories, Arrhenius rships----------------
#           RMR
# **************************
RMR_MultReg_ALL_model0 <- lmer(lnRMR ~ lnBWg + tempTestK1000 + (1|species) , data=data.rmr, REML=FALSE)
RMR_MultReg_ALL_model0.POLY <- lmer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species), data=data.rmr, REML=FALSE)
RMR_MultReg_ALL_model0intPOLY <- lmer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (1|species), data=data.rmr, REML=FALSE)
RMR_MultReg_ALL_model0int <- lmer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) , data=data.rmr, REML=FALSE)

RMR_MultReg_ALL_model1 <- lmer(lnRMR ~ lnBWg + tempTestK1000 + (1|species) + (1|species:trial), data=data.rmr, REML=FALSE)
RMR_MultReg_ALL_model1.POLY <- lmer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.rmr, REML=FALSE)
RMR_MultReg_ALL_model1intPOLY <- lmer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.rmr, REML=FALSE)
RMR_MultReg_ALL_model1int <- lmer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) + (1|species:trial), data=data.rmr, REML=FALSE)

RMR_MultReg_ALL_model2 <- lmer(lnRMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.rmr, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
RMR_MultReg_ALL_model2int <- lmer(lnRMR ~ lnBWg * tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.rmr, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
RMR_MultReg_ALL_model2.POLY <- lmer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.rmr, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
RMR_MultReg_ALL_model2intPOLY <- lmer(lnRMR ~ lnBWg * poly(tempTestK1000,2)  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.rmr, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

RMR_MultReg_ALL_model4 <- lmer(lnRMR ~ lnBWg + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr, REML=FALSE) 
RMR_MultReg_ALL_model4int <- lmer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr, REML=FALSE)
RMR_MultReg_ALL_model4.POLY <- lmer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr, REML=FALSE)
RMR_MultReg_ALL_model4intPOLY <- lmer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr, REML=FALSE)

RMR_MultReg_ALL_model5 <- lmer(lnRMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr, REML=FALSE)
RMR_MultReg_ALL_model5int <- lmer(lnRMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr, REML=FALSE)
RMR_MultReg_ALL_model5.POLY <- lmer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
RMR_MultReg_ALL_model5intPOLY <- lmer(lnRMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

#           AMR
# **************************
AMR_MultReg_ALL_model0 <- lmer(lnAMR ~ lnBWg + tempTestK1000 + (1|species) , data=data.amr, REML=FALSE)
AMR_MultReg_ALL_model0.POLY <- lmer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species), data=data.amr, REML=FALSE)
AMR_MultReg_ALL_model0intPOLY <- lmer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (1|species), data=data.amr, REML=FALSE)
AMR_MultReg_ALL_model0int <- lmer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) , data=data.amr, REML=FALSE)

AMR_MultReg_ALL_model1 <- lmer(lnAMR ~ lnBWg + tempTestK1000 + (1|species) + (1|species:trial), data=data.amr, REML=FALSE)
AMR_MultReg_ALL_model1.POLY <- lmer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.amr, REML=FALSE)
AMR_MultReg_ALL_model1intPOLY <- lmer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (1|species:trial), data=data.amr, REML=FALSE)
AMR_MultReg_ALL_model1int <- lmer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) + (1|species:trial), data=data.amr, REML=FALSE)

AMR_MultReg_ALL_model2 <- lmer(lnAMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.amr, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
AMR_MultReg_ALL_model2int <- lmer(lnAMR ~ lnBWg * tempTestK1000  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.amr, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# AMR_MultReg_ALL_model2.POLY <- lmer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.amr, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# AMR_MultReg_ALL_model2intPOLY <- lmer(lnAMR ~ lnBWg * poly(tempTestK1000,2)  + (1|species) +(0 + tempTestK1000|species) + (1|species:trial), data=data.amr, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

AMR_MultReg_ALL_model4 <- lmer(lnAMR ~ lnBWg + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr, REML=FALSE) 
AMR_MultReg_ALL_model4int <- lmer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr, REML=FALSE)
AMR_MultReg_ALL_model4.POLY <- lmer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr, REML=FALSE)
AMR_MultReg_ALL_model4intPOLY <- lmer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr, REML=FALSE)

AMR_MultReg_ALL_model5 <- lmer(lnAMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amr, REML=FALSE)
AMR_MultReg_ALL_model5int <- lmer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amr, REML=FALSE)
AMR_MultReg_ALL_model5.POLY <- lmer(lnAMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amr, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
AMR_MultReg_ALL_model5intPOLY <- lmer(lnAMR ~ lnBWg * poly(tempTestK1000,2) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amr, REML=FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

RMR_ALL_BIC<-BIC(RMR_MultReg_ALL_model0, RMR_MultReg_ALL_model0int, RMR_MultReg_ALL_model0.POLY, RMR_MultReg_ALL_model0intPOLY,
                 RMR_MultReg_ALL_model1, RMR_MultReg_ALL_model1int, RMR_MultReg_ALL_model1.POLY, RMR_MultReg_ALL_model1intPOLY,
                 RMR_MultReg_ALL_model2, RMR_MultReg_ALL_model2int,                           
                 RMR_MultReg_ALL_model4, RMR_MultReg_ALL_model4int, RMR_MultReg_ALL_model4.POLY, RMR_MultReg_ALL_model4intPOLY,
                 RMR_MultReg_ALL_model5, RMR_MultReg_ALL_model5int, RMR_MultReg_ALL_model5.POLY, RMR_MultReg_ALL_model5intPOLY)


AMR_ALL_BIC<-BIC(AMR_MultReg_ALL_model0, AMR_MultReg_ALL_model0int, AMR_MultReg_ALL_model0.POLY, AMR_MultReg_ALL_model0intPOLY,
                 AMR_MultReg_ALL_model1, AMR_MultReg_ALL_model1int, AMR_MultReg_ALL_model1.POLY, AMR_MultReg_ALL_model1intPOLY,
                 AMR_MultReg_ALL_model4, AMR_MultReg_ALL_model4int, AMR_MultReg_ALL_model4.POLY, AMR_MultReg_ALL_model4intPOLY,
                 AMR_MultReg_ALL_model5, AMR_MultReg_ALL_model5int, AMR_MultReg_ALL_model5.POLY, AMR_MultReg_ALL_model5intPOLY)

BICdelta(RMR_ALL_BIC) # 
BICdelta(AMR_ALL_BIC) # 

# 2. warm temp category - NON-Arrhenius relationships -------

AMR_model1     <- lmer(lnAMR ~ lnBWg + (1|species) + (1|species:trial), data=data.amr.test, REML=FALSE)
AMR_model12     <- lmer(lnAMR ~ lnBWg + tempTest + (1|species) + (1|species:trial), data=data.amr.test, REML=FALSE)
AMR_model12int     <- lmer(lnAMR ~ lnBWg * tempTest + (1|species) + (1|species:trial), data=data.amr.test, REML=FALSE)

AMR_model2     <- lmer(lnAMR ~ lnBWg + (tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
AMR_model22     <- lmer(lnAMR ~ lnBWg + tempTest + (tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
AMR_model22int     <- lmer(lnAMR ~ lnBWg * tempTest + (tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

AMR_model4     <- lmer(lnAMR ~ lnBWg  + (1|species) + (0+lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE)
AMR_model42     <- lmer(lnAMR ~ lnBWg + tempTest + (1|species) + (0+lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE)
AMR_model42int     <- lmer(lnAMR ~ lnBWg * tempTest + (1|species) + (0+lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE)

# AMR_model5     <- lmer(lnAMR ~ lnBWg + (lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE)
# AMR_model52     <- lmer(lnAMR ~ lnBWg + tempTest + (lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE)
# AMR_model52int     <- lmer(lnAMR ~ lnBWg * tempTest + (lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE)

# RMR
RMR_model1     <- lmer(lnRMR ~ lnBWg + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)
RMR_model12     <- lmer(lnRMR ~ lnBWg + tempTest + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)
RMR_model12int     <- lmer(lnRMR ~ lnBWg * tempTest + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)
# 
RMR_model2     <- lmer(lnRMR ~ lnBWg + (tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
RMR_model22     <- lmer(lnRMR ~ lnBWg + tempTest + (tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)
RMR_model22int     <- lmer(lnRMR ~ lnBWg * tempTest + (tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

RMR_model4     <- lmer(lnRMR ~ lnBWg  + (1|species) + (0+lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)
RMR_model42     <- lmer(lnRMR ~ lnBWg + tempTest + (1|species) + (0+lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)
RMR_model42int     <- lmer(lnRMR ~ lnBWg * tempTest + (1|species) + (0+lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE) 

# RMR_model5     <- lmer(lnRMR ~ lnBWg + (lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)
# RMR_model52     <- lmer(lnRMR ~ lnBWg + tempTest + (lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 
# RMR_model52int     <- lmer(lnRMR ~ lnBWg * tempTest + (lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 

bic.amr.warm<-BIC(AMR_model1,AMR_model12,AMR_model12int,
                  AMR_model2,AMR_model22,
                  AMR_model4,AMR_model42,AMR_model42int
                  # AMR_model5,AMR_model52,AMR_model52int
)

bic.rmr.warm<-BIC(RMR_model1,RMR_model12,RMR_model12int,
                  RMR_model2,RMR_model22,RMR_model22int,
                  RMR_model4,RMR_model42,RMR_model42int
                  # RMR_model5,RMR_model52,RMR_model52int
)

BICdelta(bic.amr.warm)
BICdelta(bic.rmr.warm)

# summary(AMR_model1) # << best 
summary(AMR_model2) # << best june 27, 2021, still holds best mar 2022
summary(RMR_model22) # << best  june 27, 2021, still holds best mar 2022

FASmodel1     <- glmer(FAS ~ lnBWg + (1|species) , data=data.fas.test, family=Gamma(link=log))
FASmodel2     <- glmer(FAS ~ lnBWg + tempTest + (1|species) , data=data.fas.test, family=Gamma(link=log))
FASmodel3     <- lmer(log(FAS) ~ lnBWg + (1|species) , data=data.fas.test, REML=FALSE)
FASmodel5     <- lmer(log(FAS) ~ lnBWg + tempTest + (1|species) , data=data.fas.test, REML=FALSE)
FASmodel5int  <- lmer(log(FAS) ~ lnBWg * tempTest + (1|species) , data=data.fas.test, REML = FALSE)

# gets too complex, singular fit
# FASmodel4   <- lmer(log(FAS) ~ lnBWg + tempTest + (lnBWg|species) + (1|species:trial) , data=data.fas.test, REML=FALSE)

FASmodel6     <- lmer(log(FAS) ~ lnBWg + tempTest + (1|species) + (1|species:trial) , data=data.fas.test, REML = FALSE)
FASmodel6int  <- lmer(log(FAS) ~ lnBWg * tempTest + (1|species) + (1|species:trial) , data=data.fas.test, REML = FALSE)


FAS.table<-BIC(FASmodel1, FASmodel2, FASmodel3,
               FASmodel5,FASmodel5int,
               FASmodel6,FASmodel6int)

BICdelta(FAS.table) # best FASmodel6 mar 22, 2022
summary(FASmodel6)






