








# Obtaining parameters: Scaling_models_multReg_categories.R data analyses --------

# # best models:
# RMR_MultReg_ER_model5 <- lmer(lnRMR ~ lnBWg + tempTestK1000  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE)
# AMR_MultReg_ER_model4int <- lmer(lnAMR ~ lnBWg * tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE)
# FAS_MultReg_ER_model4 <- lmer(log(FAS) ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE)
# AS_MultReg_ER_model4 <- lmer(lnAS ~ lnBWg + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE)
# 
# RMR_MultReg_W_model1 <- lmer(lnRMR ~ lnBWg + tempTestK1000 + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE)
# AMR_MultReg_W_model4 <- lmer(lnAMR ~ lnBWg + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE) 
# FAS_MultReg_W_model2.POLY <- lmer(lnFAS ~ lnBWg + poly(tempTest,2) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE,  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
# AS_MultReg_W_model1 <- lmer(lnAS ~ lnBWg + tempTestK1000 + (1|species) + (1|species:trial), data=data.as.test, REML=FALSE)
# 
# rmr_mod_ER<-RMR_MultReg_ER_model5
# amr_mod_ER<-AMR_MultReg_ER_model4int
# fas_mod_ER<-FAS_MultReg_ER_model4
# as_mod_ER<-AS_MultReg_ER_model4
# 
# rmr_mod_W<-RMR_MultReg_W_model1
# amr_mod_W<-AMR_MultReg_W_model4
# fas_mod_W<-FAS_MultReg_W_model2.POLY
# as_mod_W<-AS_MultReg_W_model1

# ECOL RELEV  
# # MMR_slope<-round(fixef(AMR_MultReg_ER_model4)[2],3)
# RMR_slope<-round(fixef(rmr_mod_ER)[2],3)
# # MMR_int<-round(fixef(AMR_MultReg_ER_model4)[1],3)
# RMR_int<-round(fixef(rmr_mod_ER)[1],3)
# # MMR_E<-round(fixef(AMR_MultReg_ER_model4)[3],3)
# RMR_E<-round(fixef(rmr_mod_ER)[3],3)
# FAS_slope<-round(fixef(fas_mod_ER)[2],3)
# FAS_int<-round(fixef(fas_mod_ER)[1],3)
# AS_slope<-round(fixef(as_mod_ER)[2],3)
# AS_int<-round(fixef(as_mod_ER)[1],3)


# activation energies: (Downs et al)
# slope = (-1*MR_E)/(1000*k)
# (-1*MR_E) = (1000*k) * slope
# (MR_E) = -1 * (1000*k) * slope


# # WARM 
# MMR_slope_w<-round(fixef(amr_mod_W)[2],3)
# RMR_slope_w<-round(fixef(rmr_mod_W)[2],3)
# MMR_int_w<-round(fixef(amr_mod_W)[1],3)
# RMR_int_w<-round(fixef(rmr_mod_W)[1],3)
# FAS_slope_w<-round(fixef(fas_mod_W)[2],3)
# FAS_int_w<-round(fixef(fas_mod_W)[1],3)
# AS_slope_w<-round(fixef(as_mod_W)[2],3)
# AS_int_w<-round(fixef(as_mod_W)[1],3)
# 
# # fixed effects activation energies:
# MMR_E_ER<-(emtrends(amr_mod_ER, pairwise ~ lnBWg , pbkrtest.limit = 4000, var="tempTestK1000", at=list(lnBWg=c(log(1),log(10), log(100),log(1000)))))
# MMR_E_ER<-as.data.frame(MMR_E_ER$emtrends)
# RMR_E_ER<-round(fixef(rmr_mod_ER)[3],3)
# AS_E_ER<-round(fixef(as_mod_ER)[3],3)
# 
# MMR_E_W<-round(fixef(amr_mod_ER)[3],3)
# RMR_E_W<-round(fixef(rmr_mod_W)[3],3)
# AS_E_W<-round(fixef(as_mod_W)[3],3)
# 
# # activation energies:
# # -1*(E*k)*1000
# MMR_E_ER$MMR_E_ER_eV<--1*(MMR_E_ER$tempTestK1000.trend)*k*1000
# RMR_E_ER_eV<--1*(RMR_E_ER)*k*1000
# AS_E_ER_eV<--1*(AS_E_ER)*k*1000
# 
# MMR_E_W_eV<--1*(MMR_E_W)*k*1000
# RMR_E_W_eV<--1*(RMR_E_W)*k*1000
# AS_E_W_eV<--1*(AS_E_W)*k*1000
# 
# # MMR_W_ER$MMR_E_ER_eV<--1*(MMR_E_ER$tempTestK1000.trend)*k*1000
# # RMR_W_ER_eV<--1*(RMR_W_ER_eV)*k*1000
# 
# # Obtaining predicted slopes: Scaling_models_multReg_categories.R data analyses ----------
# 
# # Figure parameters: 
# summarise.dataframes<- function(data.summarise, export.name){  
#   dataframe.new <- data.summarise %>% 
#     summarise(min_bw = min(BW_g), max_bw = max(BW_g), mean_bw = mean(BW_g), 
#               min_bwLOG = min(lnBWg), max_bwLOG = max(lnBWg), mean_bwLOG = mean(lnBWg),
#               min_temp = min(tempTest), max_temp = max(tempTest), mean_temp = mean(tempTest),
#               min_tempArrh = min(tempTestK1000), max_tempArrh = max(tempTestK1000), mean_tempArrh = mean(tempTestK1000),
#               n = length(BW_g)) %>% 
#     as.data.frame()
#   
#   assign(export.name, dataframe.new)
# }
# 
# sum_amr_ER<-summarise.dataframes(data.amrER, "sum_amr_ER")
# sum_rmr_ER<-summarise.dataframes(data.rmrER, "sum_rmr_ER")
# sum_fas_ER<-summarise.dataframes(data.fasER, "sum_fas_ER")
# sum_as_ER<-summarise.dataframes(data.asER, "sum_as_ER")
# sum_amr_warm<-summarise.dataframes(data.amr.test, "sum_amr_warm")
# sum_rmr_warm<-summarise.dataframes(data.rmr.test, "sum_rmr_warm")
# sum_fas_warm<-summarise.dataframes(data.fas.test, "sum_fas_warm")
# sum_as_warm<-summarise.dataframes(data.as.test, "sum_as_warm")
# sum_amr_ALL<-summarise.dataframes(data.amr, "sum_amr_ALL")
# sum_rmr_ALL<-summarise.dataframes(data.rmr, "sum_rmr_ALL")
# sum_fas_ALL<-summarise.dataframes(data.fas, "sum_fas_ALL")
# sum_as_ALL<-summarise.dataframes(data.as, "sum_as_ALL")
# 
# data.amr$tempTestK1<-1/data.amr$tempTestK
# data.amr$tempTestK1000<-1000/data.amr$tempTestK
# # 20 deg
# # -1*(RMR_E*k)*1000
# # tempTestK1000
# 
# C20inTempTestK1000<-1000/(((20+273.15)))
# C10inTempTestK1000<-1000/((10+273.15))
# C0inTempTestK1000<-1000/((0+273.15))
# C30inTempTestK1000<-1000/((30+273.15))
# C40inTempTestK1000<-1000/((40+273.15))
# 
# # 40 = (1000 / (C40inTempTestK1000))- 273.15
# (1000 / (C40inTempTestK1000))- 273.15 
# # ((1/C20inTempTestkT)/k)-273.15
# 
# # RMR optimal - ecologically relevant grid; plot/predict full range (All data coverage for reference)
# size_temp_grid<-expand.grid(tempTestK1000 = seq(C40inTempTestK1000, C0inTempTestK1000 , 0.1), lnBWg = seq(sum_rmr_ALL$min_bwLOG,sum_rmr_ALL$max_bwLOG, 0.5))
# data.plotRMR_ER <- data.rmrER[!duplicated(data.rmrER[,c("species","trial","test_category")]),c("species","trial", "test_category")]
# data.plotRMR_ER <- merge(size_temp_grid,data.plotRMR_ER)
# data.plotRMR_ER$model_predFE <- predict(rmr_mod_ER, re.form = NA, newdata = data.plotRMR_ER)
# 
# # FAS optimal - ecologically relevant grid; plot/predict full range (All data coverage for reference)
# size_temp_grid<-expand.grid(tempTest = 20, lnBWg = seq(sum_fas_ER$min_bwLOG,sum_fas_ER$max_bwLOG, 0.5))
# data.plotFAS_ER <- data.fasER[!duplicated(data.fasER[,c("species","trial","test_category")]),c("species","trial", "test_category")]
# data.plotFAS_ER <- merge(size_temp_grid, data.plotFAS_ER)
# data.plotFAS_ER$model_predFE <- predict(fas_mod_ER, re.form = NA, newdata = data.plotFAS_ER)
# 
# # AS optimal - ecologically relevant grid; plot/predict full range (All data coverage for reference)
# size_temp_grid<-expand.grid(tempTestK1000 = seq(C40inTempTestK1000, C0inTempTestK1000 , 0.1), lnBWg = seq(sum_as_ALL$min_bwLOG,sum_as_ALL$max_bwLOG, 0.5))
# data.plotAS_ER <- data.asER[!duplicated(data.asER[,c("species","trial","test_category")]),c("species","trial", "test_category")]
# data.plotAS_ER <- merge(size_temp_grid, data.plotAS_ER)
# data.plotAS_ER$model_predFE <- predict(as_mod_ER, re.form = NA, newdata = data.plotAS_ER)
# 
# # MMR optimal - interaction model; plot/predict full range (All data coverage for reference)
# size_temp_grid<-expand.grid(tempTestK1000 = seq(C40inTempTestK1000,C0inTempTestK1000, 0.05), lnBWg= seq(sum_amr_ALL$min_bwLOG,sum_amr_ALL$max_bwLOG, 0.5))
# data.plotAMRint_ER <- data.amrER[!duplicated(data.amrER[,c("species","trial","test_category")]),c("species","trial", "test_category")]
# data.plotAMRint_ER<-merge(size_temp_grid,data.plotAMRint_ER)
# data.plotAMRint_ER$model_predFE <- predict(amr_mod_ER, re.form = NA, newdata = data.plotAMRint_ER)
# data.plotAMRint_ER$tempTestK1000_inC<-((1000/data.plotAMRint_ER$tempTestK1000))-275.15
# 
# # predicted slopes for MMR interaction
# # emmip(amr_mod_ER, lnAMR ~ lnBWg | tempTestK1000, mult.name = c(32, 45), cov.reduce = FALSE)
# AMR.slopes<-(emtrends(amr_mod_ER, pairwise ~ tempTestK1000,
#                       pbkrtest.limit = 4000, var="lnBWg",
#                       at=list(tempTestK1000=c(C0inTempTestK1000,
#                                               C10inTempTestK1000,
#                                               C20inTempTestK1000,
#                                               C30inTempTestK1000,
#                                               C40inTempTestK1000)))) # 0, 10, 20, 30, 40 C 
# AMR.slopes<-as.data.frame(AMR.slopes$emtrends)
# AMR.slopes$tempTestK1000_inC<-((1000/AMR.slopes$tempTestK1000))-273.15
# data.amrER$tempTestK1000_inC<-((1000/data.amrER$tempTestK1000))-273.15
# 
# # All temps for the Arrhenius plot:
# 
# # # july 18, 2021
# # RMR_MultReg_ALL_model4.POLY <- lmer(lnRMR ~ lnBWg + poly(tempTestK1000,2) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr, REML=FALSE)
# # AMR_MultReg_ALL_model4 <- lmer(lnAMR ~ lnBWg + tempTestK1000 + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr, REML=FALSE) 
# # 
# # size_temp_grid<-expand.grid(tempTestK1000 = seq(C40inTempTestK1000, C0inTempTestK1000 , 0.1),
# #                             lnBWg= seq(sum_amr_ALL$min_bwLOG,sum_amr_ALL$max_bwLOG, 0.5))
# # data.plotAMR_ALL <- data.amr[!duplicated(data.amr[,c("species","trial","test_category")]),c("species","trial", "test_category")]
# # data.plotAMR_ALL<-merge(size_temp_grid,data.plotAMR_ALL)
# # data.plotAMR_ALL$model_predFE <- predict(amr_mod_ER, re.form = NA, newdata = data.plotAMR_ALL)
# # 
# # size_temp_grid<-expand.grid(tempTestK1000 = seq(C40inTempTestK1000, C0inTempTestK1000 , 0.1),
# #                             lnBWg = seq(sum_rmr_ALL$min_bwLOG,sum_rmr_ALL$max_bwLOG, 0.5))
# # data.plotRMR_ALL <- data.rmr[!duplicated(data.rmr[,c("species","trial","test_category")]),c("species","trial", "test_category")]
# # data.plotRMR_ALL<-merge(size_temp_grid,data.plotRMR_ALL)
# # data.plotRMR_ALL$model_predFE <- predict(RMR_MultReg_ALL_model4.POLY, re.form = NA, newdata = data.plotRMR_ALL) 
# # 
# # data.plotAMR_ALL$tempTestK1000_inC<-((1000/data.plotAMR_ALL$tempTestK1000))-273.15
# # data.plotRMR_ALL$tempTestK1000_inC<-((1000/data.plotRMR_ALL$tempTestK1000))-273.15
# 
# 
# # # fixed effects ALL TEMPS FOR Act Energ
# # data.plotAMR_ALL$model_predFE <- predict(AMR_MultReg_ALL_model4, re.form = NA, newdata = data.plotAMR_ALL) 
# # MMR_slope_ALL<-round(fixef(AMR_MultReg_ALL_model4)[2],2)
# # MMR_int_ALL<-round(fixef(AMR_MultReg_ALL_model4)[1],2)
# # RMR_slope_ALL<-round(fixef(RMR_MultReg_ALL_model4.POLY)[2],2)
# # RMR_int_ALL<-round(fixef(RMR_MultReg_ALL_model4.POLY)[1],2)
# 
# # RMR - warm temp grid; plot/predict full range (All data coverage for reference)
# size_temp_grid<-expand.grid(tempTestK1000 = seq(C40inTempTestK1000, C0inTempTestK1000 , 0.1), lnBWg = seq(sum_rmr_ALL$min_bwLOG,sum_rmr_ALL$max_bwLOG, 0.5))
# data.plotRMR_warm <- data.rmr.test[!duplicated(data.rmr.test[,c("species","trial","test_category")]),c("species","trial", "test_category")]
# data.plotRMR_warm<-merge(size_temp_grid,data.plotRMR_warm)
# data.plotRMR_warm$model_predFE <- predict(rmr_mod_W, re.form = NA, newdata = data.plotRMR_warm)
# 
# # size_temp_grid <- expand.grid(lnBWg = seq(sum_rmr_ER$min_bwLOG,sum_rmr_ER$max_bwLOG, 0.5), tempTestK1000 = seq())
# # data.plotRMR_warm <- data.rmr.test[!duplicated(data.rmr.test[,c("species","trial")]),c("species","trial")]
# # data.plotRMR_warm <- merge(size_temp_grid,data.plotRMR_warm)
# # data.plotRMR_warm$model_predFE <- predict(rmr_mod_W, re.form = NA, newdata = data.plotRMR_warm)
# 
# # size_temp_grid <- expand.grid(lnBWg = seq(sum_amr_ER$min_bwLOG,sum_amr_ER$max_bwLOG, 0.5), tempTest = 20)
# # data.plotAMR_warm <- data.amr.test[!duplicated(data.amr.test[,c("species","trial")]),c("species","trial")]
# # data.plotAMR_warm <- merge(size_temp_grid,data.plotAMR_warm)
# # data.plotAMR_warm$model_predFE <- predict(amr_mod_W, re.form = NA, newdata = data.plotAMR_warm)
# 
# # MMR - warm temp plot grid; plot/predict full range (All data coverage for reference)
# size_temp_grid<-expand.grid(tempTestK1000 = seq(C40inTempTestK1000, C0inTempTestK1000 , 0.1), lnBWg = seq(sum_amr_ALL$min_bwLOG,sum_amr_ALL$max_bwLOG, 0.5))
# data.plotAMR_warm <- data.amr.test[!duplicated(data.amr.test[,c("species","trial","test_category")]),c("species","trial", "test_category")]
# data.plotAMR_warm<-merge(size_temp_grid,data.plotAMR_warm)
# data.plotAMR_warm$model_predFE <- predict(amr_mod_W, re.form = NA, newdata = data.plotAMR_warm)
# 
# # FAS - warm temp plot grid; plot/predict full range (All data coverage for reference)
# size_temp_grid <- expand.grid(lnBWg = seq(sum_fas_ALL$min_bwLOG,sum_fas_ALL$max_bwLOG, 0.5), tempTest = 20)
# data.plotFAS_warm <- data.fas.test[!duplicated(data.fas.test[,c("species","trial")]),c("species","trial")]
# data.plotFAS_warm <- merge(size_temp_grid,data.plotFAS_warm)
# data.plotFAS_warm$model_predFE <- predict(fas_mod_W, re.form = NA, newdata = data.plotFAS_warm)
# 
# # AS - warm temp plot grid; plot/predict full range (All data coverage for reference)
# size_temp_grid <- expand.grid(lnBWg = seq(sum_as_ALL$min_bwLOG,sum_as_ALL$max_bwLOG, 0.5), tempTestK1000 = seq(C40inTempTestK1000, C0inTempTestK1000 , 0.1))
# data.plotAS_warm <- data.as.test[!duplicated(data.as.test[,c("species","trial")]),c("species","trial")]
# data.plotAS_warm <- merge(size_temp_grid, data.plotAS_warm)
# data.plotAS_warm$model_predFE <- predict(as_mod_W, re.form = NA, newdata = data.plotAS_warm)
# 
# 
# 
# # Figures: ecol relev scaling -----------
# 
# AMRmodel_plot1.ecol_relevINTERACTION<-ggplot(data.plotAMRint_ER) +
#   geom_line(data=data.plotAMRint_ER, aes(y = model_predFE, x=lnBWg,  group=tempTestK1000_inC, color= tempTestK1000_inC), size=0.5, lty=1, show.legend=TRUE) +
#   geom_point(data=data.amr[data.amr$test_category=="ecol_relev",], aes(x=lnBWg, y=lnAMR, fill= tempTestK1000), color="black", alpha=1, pch=21, size=2, show.legend = FALSE)+
#   ylim(-6.5, 12)+
#   xlim(-6.5, 12)+
#   scale_color_gradient( low = cols.amr[3], high = cols.amr[1])+
#   scale_fill_gradient( low = cols.amr[3], high = cols.amr[1])
# ggformat(AMRmodel_plot1.ecol_relevINTERACTION, x_title=expression(italic(ln)*Body~weight~(g)), y_title=expression(italic(ln)*MMR~(mg~O[2]~h^-1)), print = FALSE)
# 
# RMRmodel_plot1.ecol_relev<-ggplot(data.plotRMR_ER[data.plotRMR_ER$tempTestK1000==mean(data.plotRMR_ER$tempTestK1000),]) +
#   geom_point(data=data.rmr[data.rmr$test_category=="ecol_relev",], aes(x=lnBWg, y=lnRMR, fill= tempTest ),  alpha=1, pch=21, size=2, show.legend = FALSE)+
#   geom_line(data=data.plotRMR_ER[round(data.plotRMR_ER$tempTestK1000,1)==round(C20inTempTestK1000,1),], aes(y = model_predFE, x=lnBWg,  group=tempTestK1000, color= tempTestK1000), color="black", size=1, lty=1, show.legend=FALSE) +
#   ylim(-6.5, 12)+
#   xlim(-6.5, 12)+
#   scale_color_gradient( high = cols.rmr[3], low = cols.rmr[1])+
#   scale_fill_gradient( low = cols.rmr[3], high = cols.rmr[1])
# ggformat(RMRmodel_plot1.ecol_relev, x_title=expression(italic(ln)*Body~weight~(g)), y_title=expression(italic(ln)*RMR~(mg~O[2]~h^-1)), print = FALSE)
# 
# FASmodel_plot1.ecol_relev<-ggplot(data=data.fasER, aes(x=lnBWg, y=log(FAS))) +
#   geom_point(aes(x=lnBWg, y=log(FAS), fill= tempTest ),  alpha=1, pch=21, size=2, show.legend = FALSE)+
#   geom_line(data=data.plotFAS_ER[data.plotFAS_ER$tempTest==20,], aes(y = model_predFE, x=lnBWg,  group=tempTest), color="black", size=1, lty=1, show.legend=FALSE) +
#   scale_fill_gradient(high = "yellow", low = "black")+
#   ylim(0,4)
# ggformat(FASmodel_plot1.ecol_relev, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*FAS), print = FALSE)
# 
# ASmodel_plot1.ecol_relev <- ggplot(data=data.asER, aes(x=lnBWg, y=log(AS))) +
#   geom_point(aes(x=lnBWg, y=log(AS), fill= tempTest ),  alpha=1, pch=21, size=2, show.legend = FALSE)+
#   scale_fill_gradient(low = "orange", high = cols.as[3])+
#   geom_line(data=data.plotAS_ER[round(data.plotAS_ER$tempTestK1000,1)==round(C20inTempTestK1000,1),], aes(y = model_predFE, x=lnBWg,  group=tempTestK1000, color= tempTestK1000), color="black", size=1, lty=1, show.legend=FALSE) +
#   ylim(-6.5, 12)+
#   xlim(-6.5, 12)
# ggformat(ASmodel_plot1.ecol_relev, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*AS), print = FALSE)
# 
# # ggsave(filename = paste("./Figures/Fig3_fas1_",Sys.Date(), ".png", sep=""),
# #        plot=FASmodel_plot1.ecol_relev, width = 4, height = 4, units = "in")
# # 
# 
# 
# # library(LMERConvenienceFunctions)
# # library(fields)
# # 
# # plotLMER3d.fnc(model = rmr_mod_ER, intr="lnBWg", pred = "tempTestK1000", plot.type = "persp",
# #                shift = 0, scale = 1, cex = 0, contourstepsize = 0,
# #                fun = NA, color = "heat", theta = 45, phi = 25,
# #                legend.args = NULL, rug = FALSE)
# # plotLMER3d.fnc(model = rmr_mod_ER, intr="lnBWg", pred = "tempTestK1000", plot.type = "contour",
# #                shift = 0, scale = 0.5, cex = 1, contourstepsize = 0.4,
# #                fun = NA, color = "heat", 
# #                legend.args = NULL, rug = FALSE)
# 

# Figures: ecologies ----------------








# Figures: Main Scaling:  -----
data.amrAC$tt_indicVar<-1
data.amrAM$tt_indicVar<-0
data.rmrAC$tt_indicVar<-1
data.rmrAM$tt_indicVar<-0

# RMR predictions
# create data.plotRMR (data frame) that is trial-level specific 

# FAS predictions
# # create data.plotFAS (data frame) that is trial-level specific 
# data.plotFAS <- data.fas.test[,c("species","trial", "lnBWg")]
# data.plotFAS$sp_t <- as.factor(data.plotFAS$species:data.plotFAS$trial)	
# # create dummy data sets for the mean temp, 
# data.plotFAS$tempTest<-20 # 20.02848
# # adding predicted values to a dataframe 
# # random effects
# data.plotFAS$model_predRE <- predict(FASmodel5, data.plotFAS) 
# # fixed effects
# data.plotFAS$model_predFE <- predict(FASmodel5, re.form = NA, newdata = data.plotFAS) 

AMRmodel_plot1<-ggplot(data=data.amrER, aes(x=lnBWg, y=lnAMR)) +
  geom_point(alpha=0.9,  size=2, pch=19, color="grey70")+
  geom_line(data=data.plotAMRint_ER, aes(y = model_predFE, x=lnBWg,  group=tempTestK1000_inC, color = tempTestK1000_inC), size=0.3, lty=1,alpha=0.8, show.legend=FALSE) +
  geom_point(alpha=0.9,  size=2, pch=21, color="grey50",fill="grey70" )+
  geom_point(data=data.amr.test, aes(x=lnBWg, y=lnAMR, fill=tempTest), alpha=0.9,  size=2, pch=21, show.legend = FALSE)+
  geom_line(data=data.plotAMR_warm[round(data.plotAMR_warm$tempTestK1000,1)==round(C20inTempTestK1000,1),],
            aes(y = model_predFE, x=lnBWg,  group=tempTestK1000), color="#002E53", size=1, lty=1, show.legend=FALSE) +
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

ggformat(AMRmodel_plot1, x_title=expression(italic(ln)*Body~weight~(g)), y_title=expression(italic(ln)*MMR~(mg~O[2]~h^-1)), print = T)


ggsave(filename = paste("./Figures/Fig3_amr3_", Sys.Date(), ".png", sep=""),
       plot=AMRmodel_plot1, width = 4, height = 4, units = "in")
# ggsave(filename = paste("./Figures/Fig3_amr2_", Sys.Date(), ".png", sep=""),
#        plot=AMRmodel_plot1, width = 4, height = 4, units = "in")
# 

RMRmodel_plot1<-ggplot(data=data.rmrER, aes(x=lnBWg, y=lnRMR)) +
  geom_point(alpha=0.9,  size=2, pch=21, color="grey50",fill="grey70" )+
  geom_line(data=data.plotRMR_ER[round(data.plotRMR_ER$tempTestK1000,1)==round(C20inTempTestK1000,1),], mapping=aes(y = model_predFE, x=lnBWg,  group=tempTestK1000, color= tempTestK1000), color="black", size=1, lty=1, show.legend=FALSE) +
  geom_point(data=data.rmr.test, aes(x=lnBWg, y=lnRMR, fill=tempTest), alpha=0.9,  size=2, pch=21, show.legend = FALSE)+
  geom_line(data=data.plotRMR_warm[round(data.plotRMR_warm$tempTestK1000,1)==round(C20inTempTestK1000,1),],
            aes(y = model_predFE, x=lnBWg,  group=tempTestK1000), color="#7C0003", size=1, lty=1, show.legend=FALSE) +
  annotate("text",  x = -5.5, y = 11.5, label = bquote(Optimal:~italic(b)[RMR] == .(RMR_slope)),size=5, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = -5.5, y = 9.8, label = bquote(Warm:~italic(b)[RMR] == .(RMR_slope_w)),size=5, hjust=0, family="Arial", color = cols.rmr[1])+
  scale_fill_gradient( low = cols.rmr[3], high = cols.rmr[1])+
  scale_color_gradient( low = cols.rmr[3], high = cols.rmr[1])+
  ylim(x = -6.5, 12)+
  xlim(x = -6.5, 12)+
  annotate("text", label = paste("n = ", nrow(data.rmrER), sep=""),   x = -5.5, y = 8.5, size=3, hjust=0, family="Arial", color = "black")+
  annotate("text", label = paste("n = ", nrow(data.rmr.test), sep=""),  x = -5.5, y = 7.6, size=3, hjust=0, family="Arial", color = cols.rmr[1])
ggformat(RMRmodel_plot1, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*RMR~(mg~O[2]~h^-1)), print = FALSE)

ggsave(filename = paste("./Figures/Fig3_rmr3_", Sys.Date(), ".png", sep=""),
       plot=RMRmodel_plot1, width = 4, height = 4, units = "in")
# ggsave(filename = paste("./Figures/Fig3_rmr2_", Sys.Date(), ".png", sep=""),
#        plot=RMRmodel_plot1, width = 4, height = 4, units = "in")

# FAS! 
FASmodel_plot1<-ggplot(data=data.fasER, aes(x=lnBWg, y=log(FAS))) +
  geom_point(alpha=0.9,  size=2, pch=21, color="grey50",fill="grey70" )+
  geom_line(data=data.plotFAS_ER, aes(y = model_predFE, x=lnBWg,  group=tempTest), color="black", size=1, lty=1, show.legend=FALSE) +
  geom_point(data=data.fas.test, aes(x=lnBWg, y=log(FAS), fill=tempTest), alpha=0.9,  size= 2, pch=21, show.legend=FALSE)+
  geom_line(data=data.plotFAS_warm, aes(y = model_predFE, x=lnBWg,  group=tempTest), color="#475500", size=1, lty=1, show.legend=FALSE) +
  annotate("text",  x = -5.5, y = 3.9, label = bquote(Optimal:~italic(b)[FAS] == .(FAS_slope)),size=5, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = -5.5, y = 3.5, label = bquote(Warm:~italic(b)[FAS] == .(FAS_slope_w)),size=5, hjust=0, family="Arial", color = "#475500")+
  annotate("text", label = paste("n = ", nrow(data.fasER), sep=""),   x = -5.5, y = 3.1, size=3, hjust=0, family="Arial", color = "black")+
  annotate("text", label = paste("n = ", nrow(data.fas.test), sep=""),  x = -5.5, y = 2.85, size=3, hjust=0, family="Arial", color = "#475500")+
  scale_fill_gradient(high = "yellow", low = "black")+
  ylim(0,4)
ggformat(FASmodel_plot1, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*FAS), print = FALSE)
ggsave(filename = paste("./Figures/Fig3_fas3_", Sys.Date(), ".png", sep=""),
       plot=FASmodel_plot1, width = 4, height = 4, units = "in")
# ggsave(filename = paste("./Figures/Fig3_fas2_", Sys.Date(), ".png", sep=""),
#        plot=FASmodel_plot1, width = 4, height = 4, units = "in")

# AS 

ASmodel_plot1<-ggplot(data=data.asER, aes(x=lnBWg, y=lnAS)) +
  geom_point(alpha=0.9,  size=2, pch=21, color="grey50",fill="grey70" )+
  geom_line(data=data.plotAS_ER[round(data.plotAS_ER$tempTestK1000,1)==round(C20inTempTestK1000,1),], mapping=aes(y = model_predFE, x=lnBWg,  group=tempTestK1000, color= tempTestK1000), color="black", size=1, lty=1, show.legend=FALSE) +
  geom_point(data=data.as.test, aes(x=lnBWg, y=lnAS, fill=tempTest), alpha=0.9,  size=2, pch=21, show.legend = FALSE)+
  geom_line(data=data.plotAS_warm[round(data.plotAS_warm$tempTestK1000,1)==round(C20inTempTestK1000,1),],
            aes(y = model_predFE, x=lnBWg,  group=tempTestK1000), color="#00785B", size=1, lty=1, show.legend=FALSE) +
  annotate("text",  x = -5.5, y = 11.5, label = bquote(Optimal:~italic(b)[AS] == .(AS_slope)),size=5, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = -5.5, y = 9.8, label = bquote(Warm:~italic(b)[AS] == .(AS_slope_w)),size=5, hjust=0, family="Arial", color = cols.as[1])+
  scale_fill_gradient( low = cols.as[3], high = cols.as[1])+
  scale_color_gradient( low = cols.as[3], high = cols.as[1])+
  ylim(x = -6.5, 12)+
  xlim(x = -6.5, 12)+
  annotate("text", label = paste("n = ", nrow(data.asER), sep=""),   x = -5.5, y = 8.5, size=3, hjust=0, family="Arial", color = "black")+
  annotate("text", label = paste("n = ", nrow(data.as.test), sep=""),  x = -5.5, y = 7.6, size=3, hjust=0, family="Arial", color = cols.as[1])
ggformat(ASmodel_plot1, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*AS~(mg~O[2]~h^-1)), print = F)


# ggsave(filename = paste("./Figures/Fig3_fas2_", Sys.Date(), ".png", sep=""),
#        plot=FASmodel_plot1, width = 4, height = 4, units = "in")


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





# Both Activation energies together:
# MMR
AMRmodel_plot3.E_ALL<-ggplot(data.amrER[data.amrER$lnBWg==0,]) +
  geom_point(data.amrER, mapping=aes(x=tempTestK1000, y=log(mass_specamr), fill= tempTest, size = lnBWg), color="grey50", fill="grey70", alpha=1, pch=21, show.legend = FALSE)+
  geom_point(data.amr.test, mapping=aes(x=tempTestK1000, y=log(mass_specamr), fill= tempTest, size = lnBWg), color="black",  alpha=1, pch=21,  show.legend = FALSE)+
  scale_fill_gradient( low = cols.amr[3], high = cols.amr[1])+
  scale_color_gradient( low = cols.amr[3], high = cols.amr[1])+
  ylim(-5.7,4)+
  annotate("text",  x = 3.27, y = 3.5, label = bquote(italic(E)[MMR] == change~with~mass),size=5, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = 3.27, y = 2.5, label = bquote(italic(E)[MMR] == .(round(MMR_E_W_eV,3))),size=5, hjust=0, family="Arial", color = "#002E53")+
  annotate("text",  x = 3.5, y = -3.2, label = bquote(1*g~italic(E)== .(round(MMR_E_ER$MMR_E_ER_eV[1],3))),size=4, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = 3.5, y = -4, label = bquote(10*g~italic(E) == .(round(MMR_E_ER$MMR_E_ER_eV[2],3))),size=4, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = 3.5, y = -4.8, label = bquote(100*g~italic(E) == .(round(MMR_E_ER$MMR_E_ER_eV[3],3))),size=4, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = 3.5, y = -5.6, label = bquote(1000*g~italic(E) == .(round(MMR_E_ER$MMR_E_ER_eV[4],3))),size=4, hjust=0, family="Arial", color = "black")+
  scale_x_continuous(limits = c(3.25, 3.661), sec.axis = sec_axis(~ ((1000/.))-273.15, name = expression(Temperature~degree*C), breaks = c(32, 25, 17, 10, 3 )))+
  geom_line(data = data.plotAMR_warm[which(round(data.plotAMR_warm$lnBWg, 1) == round(log(1.5), 1)),], aes(y = model_predFE, x=tempTestK1000, group=lnBWg ) ,color="#002E53", size=1, lty=1, show.legend=FALSE)+
  geom_line(data = data.plotAMRint_ER[which(round(data.plotAMRint_ER$lnBWg, 1) == round(log(1.4),1)),], aes(y = model_predFE, x=tempTestK1000, group=lnBWg ) ,color="black", size=1, lty=1, show.legend=FALSE)
ggformat(AMRmodel_plot3.E_ALL, x_title= expression(Temperature^-1~(1000/T)), y_title=expression(italic(ln)*MMR~(mg~O[2]~h^-1~g^-1)), print = FALSE)

RMRmodel_plot3.E_ALL<-ggplot(data.rmrER[data.rmrER$lnBWg==0,]) +
  geom_point(data.rmrER, mapping=aes(x=tempTestK1000, y=log(mass_specrmr), fill= tempTest, size = lnBWg), color="grey50", fill="grey70", alpha=1, pch=21, show.legend = FALSE)+
  geom_point(data.rmr.test, mapping=aes(x=tempTestK1000, y=log(mass_specrmr), fill= tempTest, size = lnBWg), color = "black", alpha=1, pch=21, show.legend = FALSE)+
  scale_fill_gradient( low = cols.rmr[3], high = cols.rmr[1])+
  scale_color_gradient( low = cols.rmr[3], high = cols.rmr[1])+
  annotate("text",  x = 3.27, y = 3.5, label = bquote(italic(E)[RMR] == .(round(RMR_E_ER_eV,3))),size=5, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = 3.27, y = 2.5, label = bquote(italic(E)[RMR] == .(round(RMR_E_W_eV,3))),size=5, hjust=0, family="Arial", color = cols.rmr[1])+
  ylim(-5.7,4)+
  scale_x_continuous(limits = c(3.25, 3.661), sec.axis = sec_axis(~ ((1000/.))-273.15, name = expression(Temperature~degree*C), breaks = c(32, 25, 17, 10, 3 )))+
  geom_line(data = data.plotRMR_ER[which(round(data.plotRMR_ER$lnBWg, 1) == round(log(1.35), 1)),], aes(y = model_predFE, x=tempTestK1000, group=lnBWg ) ,color="black", size=1, lty=1, show.legend=FALSE)+
  geom_line(data = data.plotRMR_warm[which(round(data.plotRMR_warm$lnBWg, 1) == round(log(1.35), 1)),], aes(y = model_predFE, x=tempTestK1000, group=lnBWg), color = "#7C0003", size=1, lty=1, show.legend=FALSE)
ggformat(RMRmodel_plot3.E_ALL, x_title= expression(Temperature^-1~(1000/T)), y_title=expression(italic(ln)*RMR~(mg~O[2]~h^-1~g^-1)), print = FALSE)

ASmodel_plot3.E_ALL<-ggplot(data.asER[data.asER$lnBWg==0,]) +
  geom_point(data.asER, mapping=aes(x=tempTestK1000, y=log(mass_specas), fill= tempTest, size = lnBWg), color="grey50", fill="grey70", alpha=1, pch=21, show.legend = FALSE)+
  geom_point(data.as[!c(data.as$test_category=="ecol_relev"),], mapping=aes(x=tempTestK1000, y=log(mass_specas), fill= tempTest, size = lnBWg), color = "black", alpha=1, pch=21, show.legend = FALSE)+
  scale_fill_gradient( low = cols.as[3], high = cols.as[1])+
  scale_color_gradient( low = cols.as[3], high = cols.as[1])+
  annotate("text",  x = 3.27, y = 3.5, label = bquote(italic(E)[AS] == .(round(AS_E_ER_eV,3))),size=5, hjust=0, family="Arial", color = "black")+
  annotate("text",  x = 3.27, y = 2.5, label = bquote(italic(E)[AS] == .(round(AS_E_W_eV,3))),size=5, hjust=0, family="Arial", color = cols.as[1])+
  ylim(-5.7,4)+
  scale_x_continuous(limits = c(3.25, 3.661), sec.axis = sec_axis(~ ((1000/.))-273.15, name = expression(Temperature~degree*C), breaks = c(32, 25, 17, 10, 3 )))+
  geom_line(data = data.plotAS_ER[which(round(data.plotAS_ER$lnBWg, 1) == round(log(1.35), 1)),], aes(y = model_predFE, x=tempTestK1000, group=lnBWg ) ,color="black", size=1, lty=1, show.legend=FALSE)+
  geom_line(data = data.plotAS_warm[which(round(data.plotAS_warm$lnBWg, 1) == round(log(1.4), 1)),], aes(y = model_predFE, x=tempTestK1000, group=lnBWg), color = "#00785B", size=1, lty=1, show.legend=FALSE)
ggformat(ASmodel_plot3.E_ALL, x_title= expression(Temperature^-1~(1000/T)), y_title=expression(italic(ln)*as~(mg~O[2]~h^-1~g^-1)), print = T)


Arh.plot<-cowplot:::plot_grid(RMRmodel_plot3.E_ALL, AMRmodel_plot3.E_ALL,
          align = "hv",
          axis = "l",
          nrow = 1,
          ncol = 2,
          label_size = 17)
ggsave(filename = paste("./Figures/Figxx_ArrheniusFigMMR-RMR_", Sys.Date(), ".png", sep=""),
       plot=Arh.plot, width = 8.5, height = 4.5, units = "in")

Arh.plotamr<-cowplot:::plot_grid( AMRmodel_plot3.E_ALL,
                              align = "hv",
                              axis = "l",
                              nrow = 1,
                              ncol = 1,
                              label_size = 17)
ggsave(filename = paste("./Figures/Figxx_ArrheniusFigMMR_", Sys.Date(), ".png", sep=""),
       plot=Arh.plotamr, width = 4, height = 4.2, units = "in")

Arh.plotrmr<-cowplot:::plot_grid( RMRmodel_plot3.E_ALL,
                                  align = "hv",
                                  axis = "l",
                                  nrow = 1,
                                  ncol = 1,
                                  label_size = 17)
ggsave(filename = paste("./Figures/Figxx_ArrheniusFigRMR_", Sys.Date(), ".png", sep=""),
       plot=Arh.plotrmr, width = 4, height = 4.2, units = "in")














# Both MMR and RMR together - AS punchline plots ---------
polygog.x <- c( data.plotRMR_warm[which(data.plotRMR_warm$lnBWg==min(data.plotRMR_warm$lnBWg))[1], "lnBWg"],
                data.plotAMR_warm[which(data.plotAMR_warm$lnBWg==min(data.plotAMR_warm$lnBWg))[1], "lnBWg"],
                data.plotAMR_warm[which(data.plotAMR_warm$lnBWg==max(data.plotAMR_warm$lnBWg))[1], "lnBWg"],
                data.plotRMR_warm[which(data.plotRMR_warm$lnBWg==max(data.plotRMR_warm$lnBWg))[1], "lnBWg"]
)

polygog.y <- c(data.plotRMR_warm[which(data.plotRMR_warm$lnBWg==min(data.plotRMR_warm$lnBWg))[1], "model_predFE"],
               data.plotAMR_warm[which(data.plotAMR_warm$lnBWg==min(data.plotAMR_warm$lnBWg))[1], "model_predFE"],
               data.plotAMR_warm[which(data.plotAMR_warm$lnBWg==max(data.plotAMR_warm$lnBWg))[1], "model_predFE"],
               data.plotRMR_warm[which(data.plotRMR_warm$lnBWg==max(data.plotRMR_warm$lnBWg))[1], "model_predFE"]
)

polygon.plot<-data.frame(x=polygog.x, y=polygog.y)

# recreate the plot, to capture only body size range we have for warm:
size_temp_grid<-expand.grid(tempTestK1000 = seq(C40inTempTestK1000, C0inTempTestK1000 , 0.1), lnBWg = seq(sum_rmr_warm$min_bwLOG,sum_rmr_warm$max_bwLOG, 0.5))
data.plotRMR_warm <- data.rmr.test[!duplicated(data.rmr.test[,c("species","trial","test_category")]),c("species","trial", "test_category")]
data.plotRMR_warm<-merge(size_temp_grid,data.plotRMR_warm)
data.plotRMR_warm$model_predFE <- predict(rmr_mod_W, re.form = NA, newdata = data.plotRMR_warm)
conf.int_RMR_warm <- confint(rmr_mod_W, method = "boot", FUN=function(x)predict(x, data.plotRMR_warm, re.form=NA), nsim=200)
data.plotRMR_warm$CI_2.5<-conf.int_RMR_warm[,1]
data.plotRMR_warm$CI_97.5<-conf.int_RMR_warm[,2]

# MMR - warm temp plot grid; plot/predict full range (All data coverage for reference)
size_temp_grid_amr<-expand.grid(tempTestK1000 = seq(C40inTempTestK1000, C0inTempTestK1000 , 0.1), lnBWg = seq(sum_amr_warm$min_bwLOG,sum_amr_warm$max_bwLOG, 0.5))
data.plotAMR_warm <- data.amr.test[!duplicated(data.amr.test[,c("species","trial","test_category")]),c("species","trial", "test_category")]
data.plotAMR_warm<-merge(size_temp_grid_amr,data.plotAMR_warm)
data.plotAMR_warm$model_predFE <- predict(amr_mod_W, re.form = NA, newdata = data.plotAMR_warm)



        # ## Run only when needed *************
        # 
        # # last run mar 24 2022
        # 
        # # species_fam_data<-as.data.frame(matrix(nrow=0, ncol=3))
        # # colnames(species_fam_data)<-c("Species", "Species.otl", "Family")
        # # 
        # # for (i in 1:length(unique(dataMR2$species))){
        # # 
        # #   t1<-taxonomy(as.character(unique(dataMR2$Species.otl)[i]), db = "ncbi", output = "classification")
        # #   t1$name[t1$rank=="family"]
        # #   print(c(i, as.character(unique(dataMR2$Species.otl)[i])))
        # #   values<-t(c(as.character(unique(dataMR2$species)[i]), as.character(unique(dataMR2$Species.otl)[i]), t1$name[t1$rank=="family"]))
        # #   colnames(values)<-c("Species", "Species.otl", "Family")
        # #   species_fam_data<-rbind(species_fam_data,values)
        # # 
        # # }
        # 
        # # ***************
        # 
        # # write.csv(species_fam_data, "/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Data/MR-fish-metadata-data/Species_fam_mar242022.csv", row.names = FALSE )
        # species_fam_data<-read.csv("/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Data/MR-fish-metadata-data/Species_fam_mar242022.csv")
        # # ecol_scaling_sum<-read.csv("/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Data/MR-fish-metadata-data/.csv")
        # # get only DemersPelag data
        # # Climate_ecol<-ecol_scaling_sum[c(ecol_scaling_sum$Ecol=="Climate" & (ecol_scaling_sum$MR_type2 == "AMR" | ecol_scaling_sum$MR_type2 == "RMR")),c("Variable", "lnBWg_est",  "MR_type2" )]
        # # colnames(Climate_ecol)<-c("Climate",  "lnBWg_est",  "MR_type2")
        # # 
        # # dataMR2<-left_join(dataMR2, Climate_ecol, by=c("Climate", "MR_type2"))
        # # 
        # # # which are NA?
        # # # view(dataMR2[is.na(dataMR2$lnBWg_est),])
        # # dataMR2[c(is.na(dataMR2$lnBWg_est) & dataMR2$MR_type2 == "AMR"), "lnBWg_est"] <- 0.80
        # # dataMR2[c(is.na(dataMR2$lnBWg_est) &  dataMR2$MR_type2 == "RMR"), "lnBWg_est"] <- 0.84
        # # 
        # # dataMR2<-left_join(dataMR2, species_fam_data, by=c("Species"))
        # # split_fam<-split(dataMR2, dataMR2$Family) # 32 families
        # # # add another category - MR type = ecol_relev and acclim lump together and acute is different
        # # 
        # # for (i in 1:length(split_fam)){
        # # 
        # #   fam_data<-as.data.frame(split_fam[i])
        # #   colnames(fam_data)<-colnames(dataMR2)
        # #   fam_data$MR_mass_corr<-NA
        # #   fam_data$MR_mass_corr[fam_data$MR_type2=="AMR"] <- fam_data$MR[fam_data$MR_type2=="AMR"] /fam_data$BW_g[fam_data$MR_type2=="AMR"] ^ fam_data$lnBWg_est[fam_data$MR_type2=="AMR"][1] # temperate
        # #   fam_data$MR_mass_corr[fam_data$MR_type2=="RMR"] <- fam_data$MR[fam_data$MR_type2=="RMR"] /fam_data$BW_g[fam_data$MR_type2=="RMR"] ^ fam_data$lnBWg_est[fam_data$MR_type2=="RMR"][1] # temperate
        # # 
        # #   # summary(split_fam)
        # #   plot<-ggplot(fam_data, aes(tempTest, MR_mass_corr, fill=test_category, color=test_category, group = interaction(test_category2, MR_type2), shape=MR_type2, label=study_ID))+
        # #     geom_point(size=3,  alpha=0.6)+
        # #     # geom_boxplot(aes(group =  interaction(test_category2, MR_type2)))+
        # #     # geom_text(check_overlap = TRUE)+
        # #     scale_shape_manual(values = c("AMR" = 8, "RMR" = 25))+
        # #     scale_color_manual(values = c("acclim" = "#C70039", "ecol_relev" = "#00929A", "acute" = "#00749F")) +
        # #     scale_fill_manual(values = c("acclim" = "#C70039", "ecol_relev" = "#00929A", "acute" = "#00749F")) +
        # #     facet_grid(.~test_category2)+
        # #     ylab(expression(MR~(mgO[2]~g^-1~h^-1)))+
        # #     xlab(Temperature~degree*C)+
        # #     ggtitle(fam_data$Family[1])+
        # #     theme_classic()+
        # #     geom_smooth(method="gam", formula = y ~ s(x, k = 3, bs = "cs"), color = "black", fill = "grey")
        # #   # stat_smooth(method = "lm", formula = y ~ poly(x, 3), size = 1)
        # #   assign(paste("plot", i, fam_data$Family[1], sep=""), plot)
        # #   ggsave(filename = paste("/Users/kristakraskura/Desktop/BOX/UCSB/Research/Metabolic_scaling/ms-AMR-RMR-Temperature/Figures/Family_taxa/plot", i, fam_data$Family[1],".png", sep=""), plot = plot, height = 4, width = 7, units = "in", dpi=300 )
        # #   print(paste("plot", i, fam_data$Family[1], sep=""))
        # # 
        # #   plot2<-ggplot(fam_data, aes(x=lnBWg, y=lnMR, fill=tempTest, color=tempTest,linetype=MR_type2,  group = interaction(test_category2, MR_type2), shape=MR_type2, label=study_ID))+
        # #     geom_point(size=2,  alpha=0.6)+
        # #     # geom_boxplot(aes(group =  interaction(test_category2, MR_type2)))+
        # #     # geom_text(check_overlap = TRUE)+
        # #     scale_shape_manual(values = c("AMR" = 8, "RMR" = 25))+
        # #     # scale_color_manual(values = cols.amr) +
        # #     # scale_fill_manual(values = cols.amr) +
        # #     facet_grid(.~test_category2)+
        # #     ylab(expression(italic(ln)*MR~(mgO[2]~g^-1~h^-1)))+
        # #     xlab(expression(italic(ln)~Body~mass~(g)))+
        # #     ggtitle(fam_data$Family[1])+
        # #     theme_classic()+
        # #     geom_smooth(method="lm", formula = y ~ x, fill = "black", color="black")
        # #   # stat_smooth(method = "lm", formula = y ~ poly(x, 3), size = 1)
        # #   assign(paste("plot", i, "-2",fam_data$Family[1], sep=""), plot2)
        # #   ggsave(filename = paste("/Users/kristakraskura/Desktop/BOX/UCSB/Research/Metabolic_scaling/ms-AMR-RMR-Temperature/Figures/Family_taxa/plot2", i,"-2", fam_data$Family[1],".png", sep=""), plot = plot2, height = 4, width = 7, units = "in", dpi=300 )
        # #   print(paste("plot", i, "-2", fam_data$Family[1], sep=""))
        # # 
        # # 
        # # }
        # 
        # # # salmonid data
        # colnames(species_fam_data)<-c("species", "Species.otl", "Family")
        # dataMR2<-merge(dataMR2, species_fam_data, by = "Species.otl", all.x = TRUE)
        # salmonids <- dataMR2[dataMR2$Family == "Salmonidae",]
        # cyrinids <- dataMR2[dataMR2$Family == "Cyprinidae",]
        # 
        # salmonids$sizeClass <- "Adult"
        # salmonids[salmonids$BW_g<800,"sizeClass"] <- "Juvies"
        # 
        # cyrinids$sizeClass <- "Adult"
        # cyrinids[cyrinids$BW_g<30,"sizeClass"] <- "Juvies"
        # 
        # # salmonids$MR_mass_corr[salmonids$MR_type2=="AMR"] <- salmonids$MR[salmonids$MR_type2=="AMR"] /salmonids$BW_g[salmonids$MR_type2=="AMR"] ^ salmonids$lnBWg_est[salmonids$MR_type2=="AMR"][1] # temperate
        # # salmonids$MR_mass_corr[salmonids$MR_type2=="RMR"] <- salmonids$MR[salmonids$MR_type2=="RMR"] /salmonids$BW_g[salmonids$MR_type2=="RMR"] ^ salmonids$lnBWg_est[salmonids$MR_type2=="RMR"][1] # temperate
        # 
        # # fit the linear model to this
        # salmonids$sizeClass<-as.factor(salmonids$sizeClass)
        # model<-lm(lnMR ~ lnBWg + sizeClass:tempTestK1_1000 + lnBWg:sizeClass, data = salmonids[salmonids$MR_type2=="RMR",])
        # summary(model)
        # 
        # # interaction between maximum, resting and temperature 
        # # fisheries, mmr and rmr physio 
        # 
        # prediction<-as.data.frame(ggpredict(model, terms = ~ lnBWg*sizeClass + tempTestK1_1000*sizeClass))
        # plot(prediction)
        # 
        # ggplot(salmonids, aes(tempTestK1_1000, MR/BW_g, fill=test_category, group = interaction(MR_type2, sizeClass), shape=MR_type2, linetype = MR_type2, size = BW_g, label=study_ID))+
        #   geom_point( size = 3, color="black", alpha=0.6)+
        #   scale_shape_manual(values = c(21, 25))+
        #   scale_fill_viridis_d(option = "D")+
        #   theme_classic()+
        #   facet_wrap(.~sizeClass)+
        #   stat_smooth(method='lm', color = "black")
        # # formula = y~poly(x,2), 
        # 
        # ggplot(cyrinids, aes(tempTest, MR/BW_g, fill=test_category, group = interaction(MR_type2, sizeClass), shape=MR_type2, size = BW_g, label=study_ID))+
        #   geom_point( size = 3, color="black", alpha=0.6)+
        #   # geom_boxplot(aes(group =  interaction(test_category2, MR_type2)))+
        #   # geom_text(check_overlap = TRUE)+
        #   scale_shape_manual(values = c(21, 25))+
        #   scale_fill_viridis_d(option = "D")+
        #   theme_classic()+
        #   # facet_grid(.~sizeClass, scales = "free")+
        #   stat_smooth(method='lm', formula = y~poly(x,2), color = "black")
        # 
        # 
        # 
        # 
        # 
