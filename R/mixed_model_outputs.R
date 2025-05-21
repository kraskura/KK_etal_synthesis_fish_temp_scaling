# Model scaling parameters and CIs ----------

model_outputs<-function(phylo=TRUE,
                        best.model.rmr.er,
                        best.model.amr.er,
                        best.model.as.er,
                        best.model.fas.er,
                        best.model.rmr.w,
                        best.model.amr.w,
                        best.model.as.w,
                        best.model.fas.w,
                        estimate.CI = TRUE){
  if(phylo){
    foldername.phylo<-"Phylo"
  }else{
    foldername.phylo<-"nonPhylo"
  }
  
  
  if(estimate.CI){
    
    message("Estimating confidence intervals, output saved in './Data_exports/' ")
    message("RMR, AMR, FAS, and AS optimal temp slope estimates over temp range 0 - 31 ºC")

    dir.create(paste("./Data_exports/", foldername.phylo, sep =""), recursive =TRUE)

    # RMR optimal - ecologically relevant grid; plot/predict full range (All data coverage for reference)
    data.plotRMR_ER<-expand.grid(tempTest = seq(0, 31 , 1), lnBWg = seq(-3.2189, 11.7440, 0.5))
    data.plotRMR_ER$model_predFE <- predict(best.model.rmr.er, re.form = NA, newdata = data.plotRMR_ER)
    conf.int_RMR_ER <- confint(best.model.rmr.er, method = "boot",
                               FUN=function(x)predict(x, data.plotRMR_ER, re.form=NA),
                               nsim=10)
    data.plotRMR_ER$CI_2.5<-conf.int_RMR_ER[,1]
    data.plotRMR_ER$CI_97.5<-conf.int_RMR_ER[,2]
  
    # FAS optimal - ecologically relevant grid; plot/predict full range (All data coverage for reference)
    data.plotFAS_ER<-expand.grid(tempTest = seq(0, 31 , 1), lnBWg = seq(-3.219, 11.744, 0.5))
    data.plotFAS_ER$model_predFE <- predict(best.model.fas.er, re.form = NA, newdata = data.plotFAS_ER)
    conf.int_FAS_ER <- confint(best.model.fas.er, method = "boot",
                               FUN=function(x)predict(x, data.plotFAS_ER, re.form=NA),
                               nsim=10)
    data.plotFAS_ER$CI_2.5<-conf.int_FAS_ER[,1]
    data.plotFAS_ER$CI_97.5<-conf.int_FAS_ER[,2]
  
    # AS optimal - ecologically relevant grid; plot/predict full range (All data coverage for reference)
    data.plotAS_ER<-expand.grid(tempTest = seq(0, 31 , 1), lnBWg = seq(-3.219, 11.744, 0.5))
    data.plotAS_ER$model_predFE <- predict(best.model.as.er, re.form = NA, newdata = data.plotAS_ER)
    conf.int_AS_ER <- confint(best.model.as.er, method = "boot",
                              FUN=function(x)predict(x, data.plotAS_ER,
                                                     re.form=NA), nsim=10)
    data.plotAS_ER$CI_2.5<-conf.int_AS_ER[,1]
    data.plotAS_ER$CI_97.5<-conf.int_AS_ER[,2]
  
    # MMR or AMR  optimal - interaction model; plot/predict full range (All data coverage for reference)
    data.plotAMRint_ER <- expand.grid(tempTest = seq(0, 31 , 1), lnBWg= seq(-3.5405,11.9512, 0.5))
    data.plotAMRint_ER$model_predFE <- predict(best.model.amr.er, re.form = NA, newdata = data.plotAMRint_ER)
    conf.int_AMRint_ER <- confint(best.model.amr.er, method = "boot",
                                  FUN=function(x)predict(x, data.plotAMRint_ER, re.form=NA),
                                  nsim=10)
    data.plotAMRint_ER$CI_2.5<-conf.int_AMRint_ER[,1]
    data.plotAMRint_ER$CI_97.5<-conf.int_AMRint_ER[,2]

    # All temps :
    # RMR - warm temp grid; plot/predict full range (All data coverage for reference)
    data.plotRMR_warm<-expand.grid(tempTest = seq(1, 37 ,1), lnBWg = seq(-1.121, 8.756 , 0.5))
    data.plotRMR_warm$model_predFE <- predict(best.model.rmr.w, re.form = NA, newdata = data.plotRMR_warm)
    conf.int_RMR_warm <- confint(best.model.rmr.w, method = "boot", FUN=function(x)predict(x, data.plotRMR_warm, re.form=NA), nsim=10)
    data.plotRMR_warm$CI_2.5<-conf.int_RMR_warm[,1]
    data.plotRMR_warm$CI_97.5<-conf.int_RMR_warm[,2]
  
    # MMR - warm temp plot grid; plot/predict full range (All data coverage for reference)
    data.plotAMR_warm<-expand.grid(tempTest = seq(1, 39 , 0.05), lnBWg = seq(-1.619,8.824, 0.5))
    data.plotAMR_warm$model_predFE <- predict(best.model.amr.w, re.form = NA, newdata = data.plotAMR_warm)
    conf.int_AMR_warm <- confint(best.model.amr.w, method = "boot", FUN=function(x)predict(x, data.plotAMR_warm, re.form=NA), nsim=10)
    data.plotAMR_warm$CI_2.5<-conf.int_AMR_warm[,1]
    data.plotAMR_warm$CI_97.5<-conf.int_AMR_warm[,2]
  
    # FAS - warm temp plot grid; plot/predict full range (All data coverage for reference)
    data.plotFAS_warm <- expand.grid(tempTest = seq(1, 35 , 1), lnBWg = seq(-1.121,8.068, 0.5))
    data.plotFAS_warm$model_predFE <- predict(best.model.fas.w, re.form = NA, newdata = data.plotFAS_warm)
    conf.int_FAS_warm <- confint(best.model.fas.w, method = "boot", FUN=function(x)predict(x, data.plotFAS_warm, re.form=NA), nsim=10)
    data.plotFAS_warm$CI_2.5<-conf.int_FAS_warm[,1]
    data.plotFAS_warm$CI_97.5<-conf.int_FAS_warm[,2]
  
    # AS - warm temp plot grid; plot/predict full range (All data coverage for reference)
    data.plotAS_warm <- expand.grid(tempTest = seq(1, 35 , 1), lnBWg = seq(-1.121,8.068, 0.5))
    data.plotAS_warm$model_predFE <- predict(best.model.as.w, re.form = NA, newdata = data.plotAS_warm)
    conf.int_AS_warm <- confint(best.model.as.w, method = "boot", FUN=function(x)predict(x, data.plotAS_warm, re.form=NA), nsim=10)
    data.plotAS_warm$CI_2.5<-conf.int_AS_warm[,1]
    data.plotAS_warm$CI_97.5<-conf.int_AS_warm[,2]
  
    # CIs for all model parameters ---------
    CI.amr.ER<-as.data.frame(confint.merMod(best.model.amr.er, level = 0.90))
    CI.amr.ER$var<-rownames(CI.amr.ER)
    CI.amr.ER$MR<-"MMR"
    CI.amr.ER$temp_cat<-"ER"
  
    CI.rmr.ER<-as.data.frame(confint.merMod(best.model.rmr.er, level = 0.90))
    CI.rmr.ER$var<-rownames(CI.rmr.ER)
    CI.rmr.ER$MR<-"RMR"
    CI.rmr.ER$temp_cat<-"ER"
  
    CI.fas.ER<-as.data.frame(confint.merMod(best.model.fas.er, level = 0.90))
    CI.fas.ER$var<-rownames(CI.fas.ER)
    CI.fas.ER$MR<-"FAS"
    CI.fas.ER$temp_cat<-"ER"
  
    CI.as.ER<-as.data.frame(confint.merMod(best.model.as.er, level = 0.90))
    CI.as.ER$var<-rownames(CI.as.ER)
    CI.as.ER$MR<-"AS"
    CI.as.ER$temp_cat<-"ER"
  
    CI.amr.W<-as.data.frame(confint.merMod(best.model.amr.w, level = 0.90))
    CI.amr.W$var<-rownames(CI.amr.W)
    CI.amr.W$MR<-"MMR"
    CI.amr.W$temp_cat<-"W"
  
    CI.rmr.W<-as.data.frame(confint.merMod(best.model.rmr.w, level = 0.90))
    CI.rmr.W$var<-rownames(CI.rmr.W)
    CI.rmr.W$MR<-"RMR"
    CI.rmr.W$temp_cat<-"W"
  
    CI.fas.W<-as.data.frame(confint.merMod(best.model.fas.w, level = 0.90))
    CI.fas.W$var<-rownames(CI.fas.W)
    CI.fas.W$MR<-"FAS"
    CI.fas.W$temp_cat<-"W"
  
    CI.as.W<-as.data.frame(confint.merMod(best.model.as.w, level = 0.90))
    CI.as.W$var<-rownames(CI.as.W)
    CI.as.W$MR<-"AS"
    CI.as.W$temp_cat<-"W"
  
    sum_CItable<-rbind(CI.amr.ER, CI.rmr.ER, CI.fas.ER, CI.as.ER, # ecol relev
                       CI.amr.W, CI.rmr.W, CI.fas.W, CI.as.W) # warm
  
    # Saving all confidence interval and predictions data frame  ---------
    write.csv(file = paste("./Data_exports/", foldername.phylo,"/Table_CIsummary.csv", sep = ""),
              sum_CItable, row.names=TRUE)
  
    write.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_RMR_er.csv", sep=""),
              data.plotRMR_ER, row.names = F)
    write.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_AMR_er.csv", sep=""),
              data.plotAMRint_ER, row.names = F)
    write.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_AS_er.csv", sep=""),
              data.plotAS_ER, row.names = F)
    write.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_FAS_er.csv", sep=""),
              data.plotFAS_ER, row.names = F)
  
    write.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_RMR_warm.csv", sep=""),
              data.plotRMR_warm, row.names = F)
    write.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_AMR_warm.csv", sep=""),
              data.plotAMR_warm, row.names = F)
    write.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_AS_warm.csv", sep=""),
              data.plotAS_warm, row.names = F)
    write.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_FAS_warm.csv", sep=""),
              data.plotFAS_warm, row.names = F)


  }else{
    message("NOT estimating CI, last saved data imported from './Data_exports/'")
    
    sum_CItable<-read.csv(file = paste("./Data_exports/", foldername.phylo,"/Table_CIsummary.csv", sep =""))
  
    data.plotRMR_ER<-read.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_RMR_er.csv", sep=""))
    data.plotAMRint_ER<-read.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_AMR_er.csv", sep=""))
    data.plotAS_ER<-read.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_AS_er.csv", sep=""))
    data.plotFAS_ER<-read.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_FAS_er.csv", sep=""))
  
    data.plotRMR_warm<-read.csv(paste("./Data_exports/", foldername.phylo,"/dataPred_RMR_warm.csv", sep=""))
    data.plotAMR_warm<-read.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_AMR_warm.csv", sep=""))
    data.plotAS_warm<-read.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_AS_warm.csv", sep=""))
    data.plotFAS_warm<-read.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_FAS_warm.csv", sep=""))
  }
 
  k<-(8.62*10^(-5)) # Boltzmann's constant

  # ECOL RELEV
  RMR_slope<-round(coef(summary(best.model.rmr.er))[2],3)
  RMR_int<-round(coef(summary(best.model.rmr.er))[1],3)
  FAS_slope<-round(coef(summary(best.model.fas.er))[2],3)
  FAS_int<-round(coef(summary(best.model.fas.er))[1],3)
  AS_slope<-round(coef(summary(best.model.as.er))[2],3)
  AS_int<-round(coef(summary(best.model.as.er))[1],3)

  # WARM
  MMR_slope_w<-round(coef(summary(best.model.amr.w))[2],3)
  MMR_int_w<-round(coef(summary(best.model.amr.w))[1],3)
  RMR_slope_w<-round(coef(summary(best.model.rmr.w))[2],3)
  RMR_int_w<-round(coef(summary(best.model.rmr.w))[1],3)
  FAS_slope_w<-round(coef(summary(best.model.fas.w))[2],3)
  FAS_int_w<-round(coef(summary(best.model.fas.w))[1],3)
  AS_slope_w<-round(coef(summary(best.model.as.w))[2],3)
  AS_int_w<-round(coef(summary(best.model.as.w))[1],3)
  
  # parameter correlations 
  # plot(model.matrix(best.model.fas.w)[, 3], model.matrix(best.model.fas.w)[, 4])
  # plot(model.matrix(best.model.as.w)[, 3], model.matrix(best.model.as.w)[, 4])

  # C20inTempTestK1000<-1000/(((20+273.15)))
  # C10inTempTestK1000<-1000/((10+273.15))
  # C0inTempTestK1000<-1000/((0+273.15))
  # C30inTempTestK1000<-1000/((30+273.15))
  # C40inTempTestK1000<-1000/((40+273.15))

  # predicted slopes for MMR interaction 
  AMR.slopes<-(emtrends(best.model.amr.er, pairwise ~ tempTest,
                        pbkrtest.limit = 4000, var="lnBWg",
                        at=list(tempTest=c(5,
                                          15,
                                          25,
                                          35)))) # 0, 10, 20, 30, 35 C
  
  AMR.slopes<-data.frame(AMR.slopes$emtrends)
  # AMR.slopes<-AMR.slopes
  AMR.slopes$performance<-"MMR"
  AMR.slopes$temp_categ<-"er"
  AMR.slopes$`(Intercept)` <- NA
  
  # AMR.slopes$`E(ev)`<-NA
  # fixed effects activation energies:
  # MMR_E_ER<-(emtrends(best.model.amr.er, pairwise ~ lnBWg , pbkrtest.limit = 4000,
  #                     var = "tempTest",
  #                     at = list(lnBWg = c(log(1),log(10), log(100),log(1000)))))
  # MMR_E_ER<-as.data.frame(MMR_E_ER$emtrends)
  # MMR_E_ER$performance<-"MMR"
  # MMR_E_ER$temp_categ<-"er"
  # MMR_E_ER$lnBWg.trend<-NA
  # MMR_E_ER$`(Intercept)` <- NA
  # MMR_E_ER$tempTest <- NA

  # RMR_E_ER<-round(coef(summary(best.model.rmr.er))[3],3)
  # AS_E_ER<-round(coef(summary(best.model.as.er))[3],3)
  # MMR_E_W<-round(coef(summary(best.model.amr.w))[3],3)
  # RMR_E_W<-round(coef(summary(best.model.rmr.w))[3],3)
  # AS_E_W<-round(coef(summary(best.model.as.w))[3],3)

  # activation energies:
  # -1*(E*k)*1000
  # MMR_E_ER$MMR_E_ER_eV<--1*(MMR_E_ER$tempTestK1000.trend)*k*1000
  # RMR_E_ER_eV<--1*(RMR_E_ER)*k*1000
  # AS_E_ER_eV<--1*(AS_E_ER)*k*1000
  # 
  # MMR_E_W_eV<--1*(MMR_E_W)*k*1000
  # RMR_E_W_eV<--1*(RMR_E_W)*k*1000
  # AS_E_W_eV<--1*(AS_E_W)*k*1000

  # Get dataframe with slopes and intercepts ----- 
  scaling.params<-as.data.frame(matrix(ncol = 5, nrow = 0))
  colnames(scaling.params)<-c("performance",
                              "temp_categ",
                              "lnBWg",
                              "(Intercept)",
                              "tempTest")

  RMR_er_row<-as.data.frame(t(c("performance" = "RMR",
                  "temp_categ" = "er",
                  RMR_slope,
                  RMR_int, 
                  "tempTest" = NA)))
  RMR_W_row<-as.data.frame(t(c("performance" = "RMR",
                  "temp_categ" = "warm",
                  RMR_slope_w,
                  RMR_int_w, 
                  "tempTest" = NA)))
  FAS_er_row<-as.data.frame(t(c("performance" = "FAS",
                  "temp_categ" = "er",
                  FAS_slope,
                  FAS_int,
                  "tempTest" = NA)))
  FAS_W_row<-as.data.frame(t(c("performance" = "FAS",
                  "temp_categ" = "warm",
                  FAS_slope_w,
                  FAS_int_w,
                  "tempTest" = NA)))
  AS_er_row<-as.data.frame(t(c("performance" = "AS",
                  "temp_categ" = "er",
                  AS_slope, 
                  AS_int,
                  "tempTest" = NA)))
  AS_W_row<-as.data.frame(t(c("performance" = "AS",
                  "temp_categ" = "warm",
                  AS_slope_w, 
                  AS_int_w,
                  "tempTestK" = NA)))
  MMR_W_row<-as.data.frame(t(c("performance" = "MMR",
                  "temp_categ" = "warm",
                  MMR_slope_w,
                  MMR_int_w,
                  "tempTestK" = NA)))

  MMR_er_row<-AMR.slopes[, c("performance", "temp_categ", "lnBWg.trend",
                              "(Intercept)", "tempTest")]

  colnames(MMR_er_row)<-colnames(scaling.params)
  colnames(MMR_W_row)<-colnames(scaling.params)
  colnames(RMR_er_row)<-colnames(scaling.params)
  colnames(RMR_W_row)<-colnames(scaling.params)
  colnames(AS_er_row)<-colnames(scaling.params)
  colnames(AS_W_row)<-colnames(scaling.params)
  colnames(FAS_er_row)<-colnames(scaling.params)
  colnames(FAS_W_row)<-colnames(scaling.params)

  scaling.params<-rbind(MMR_er_row, MMR_W_row, 
                        RMR_er_row, RMR_W_row,
                        AS_er_row, AS_W_row, 
                        FAS_er_row, FAS_W_row)
  
  # MMR_W_ER$MMR_E_ER_eV<--1*(MMR_E_ER$tempTestK1000.trend)*k*1000
  # RMR_W_ER_eV<--1*(RMR_W_ER_eV)*k*1000

  # Summary and figure parameters: ---- 
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

  # save files -----
  # save scaling parameters 
  write.csv(file = paste("./Data_exports/", foldername.phylo,"/scaling_parameters.csv", sep=""),
              scaling.params, row.names = F)
  # save summary data for each metric
  write.csv(file = paste("./Data_exports/", foldername.phylo,"/summary_data_all_metrics.csv", sep=""), sum_data, row.names = F)
          
  # 20 deg
  # -1*(RMR_E*k)*1000
  # tempTestK1000
  # 
  # # needed values; Arrh only, not used
  # C20inTempTestK1000<-1000/(((20+273.15)))
  # C10inTempTestK1000<-1000/((10+273.15))
  # C0inTempTestK1000<-1000/((0+273.15))
  # C30inTempTestK1000<-1000/((30+273.15))
  # C40inTempTestK1000<-1000/((40+273.15))

  message("Return: 1) sum_CItable, 2) scaling.params, 3) sum_data")
  return(list(sum_CItable, scaling.params, sum_data))
}
