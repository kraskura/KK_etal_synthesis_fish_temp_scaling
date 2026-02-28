# Model scaling parameters and CIs ----------
# test
# best.model.rmr.er= rmr_mod_ER
#                 best.model.amr.er= amr_mod_ER
#                 best.model.as.er = as_mod_ER
#                 best.model.fas.er = fas_mod_ER
#                 best.model.rmr.w = rmr_mod_W
#                 best.model.amr.w = amr_mod_W
#                 best.model.as.w = as_mod_W
#                 best.model.fas.w = fas_mod_W
#                 data.rmr.test = data.rmr.test
#                 data.rmrER = data.rmrER
#                 data.amr.test = data.amr.test
#                 data.amrER = data.amrER                        
#                 data.as.test = data.as.test
#                 data.asER = data.asER 
#                 data.fas.test = data.fas.test
#                 data.fasER = data.fasER 
#                 name_extension = paste(lifestage,"models", sep="")
              
get_model_outputs<-function(
                        best.model.rmr.er,
                        best.model.amr.er,
                        best.model.as.er,
                        best.model.fas.er,
                        best.model.rmr.w,
                        best.model.amr.w,
                        best.model.as.w,
                        best.model.fas.w,
                        # used to get min max values for datasets
                        data.rmr.test,
                        data.rmrER, 
                        data.amr.test,
                        data.amrER,                        
                        data.as.test,
                        data.asER, 
                        data.fas.test,
                        data.fasER, 
                        name_extension){
  
  dir.create(paste("./Data_exports/", name_extension, sep =""), recursive =TRUE)

  # RMR optimal - ecologically relevant grid; plot/predict full range (All data coverage for reference)
  data.plotRMR_ER<-expand.grid(tempTest = seq(min(data.rmrER$tempTest), max(data.rmrER$tempTest) , 1),
                               lnBWg = seq(min(data.rmrER$lnBWg), max(data.rmrER$lnBWg), 0.5))
  data.plotRMR_ER$model_predFE <- predict(best.model.rmr.er, re.form = NA, newdata = data.plotRMR_ER)
  # conf.int_RMR_ER <- confint(best.model.rmr.er, method = "boot",
  #                            FUN=function(x)predict(x, data.plotRMR_ER, re.form=NA),
  #                            nsim=10)
  # data.plotRMR_ER$CI_2.5<-conf.int_RMR_ER[,1]
  # data.plotRMR_ER$CI_97.5<-conf.int_RMR_ER[,2]

  # FAS optimal - ecologically relevant grid; plot/predict full range (All data coverage for reference)
  data.plotFAS_ER<-expand.grid(tempTest = seq(min(data.fasER$tempTest), max(data.fasER$tempTest) , 1),
                               lnBWg = seq(min(data.fasER$lnBWg), max(data.fasER$lnBWg), 0.5))
  data.plotFAS_ER$model_predFE <- predict(best.model.fas.er, re.form = NA, newdata = data.plotFAS_ER)
  # conf.int_FAS_ER <- confint(best.model.fas.er, method = "boot",
  #                            FUN=function(x)predict(x, data.plotFAS_ER, re.form=NA),
  #                            nsim=10)
  # data.plotFAS_ER$CI_2.5<-conf.int_FAS_ER[,1]
  # data.plotFAS_ER$CI_97.5<-conf.int_FAS_ER[,2]

  # AS optimal - ecologically relevant grid; plot/predict full range (All data coverage for reference)
  data.plotAS_ER<-expand.grid(tempTest = seq(min(data.asER$tempTest), max(data.asER$tempTest) , 1),
                               lnBWg = seq(min(data.asER$lnBWg), max(data.asER$lnBWg), 0.5))
  data.plotAS_ER$model_predFE <- predict(best.model.as.er, re.form = NA, newdata = data.plotAS_ER)
  # conf.int_AS_ER <- confint(best.model.as.er, method = "boot",
  #                           FUN=function(x)predict(x, data.plotAS_ER,
  #                                                  re.form=NA), nsim=10)
  # data.plotAS_ER$CI_2.5<-conf.int_AS_ER[,1]
  # data.plotAS_ER$CI_97.5<-conf.int_AS_ER[,2]

  # MMR or AMR  optimal - interaction model; plot/predict full range (All data coverage for reference)
  data.plotAMRint_ER <- expand.grid(tempTest = seq(min(data.amrER$tempTest), max(data.amrER$tempTest) , 1),
                               lnBWg = seq(min(data.amrER$lnBWg), max(data.amrER$lnBWg), 0.5))
  data.plotAMRint_ER$model_predFE <- predict(best.model.amr.er, re.form = NA, newdata = data.plotAMRint_ER)
  # conf.int_AMRint_ER <- confint(best.model.amr.er, method = "boot",
  #                               FUN=function(x)predict(x, data.plotAMRint_ER, re.form=NA),
  #                               nsim=10)
  # data.plotAMRint_ER$CI_2.5<-conf.int_AMRint_ER[,1]
  # data.plotAMRint_ER$CI_97.5<-conf.int_AMRint_ER[,2]

  # All temps :
  # RMR - warm temp grid; plot/predict full range (All data coverage for reference)
  data.plotRMR_warm<-expand.grid(tempTest = seq(min(data.rmr.test$tempTest), max(data.rmr.test$tempTest) , 1),
                               lnBWg = seq(min(data.rmr.test$lnBWg), max(data.rmr.test$lnBWg), 0.5))
  data.plotRMR_warm$model_predFE <- predict(best.model.rmr.w, re.form = NA, newdata = data.plotRMR_warm)
  # conf.int_RMR_warm <- confint(best.model.rmr.w, method = "boot", FUN=function(x)predict(x, data.plotRMR_warm, re.form=NA), nsim=10)
  # data.plotRMR_warm$CI_2.5<-conf.int_RMR_warm[,1]
  # data.plotRMR_warm$CI_97.5<-conf.int_RMR_warm[,2]

  # MMR - warm temp plot grid; plot/predict full range (All data coverage for reference)
  data.plotAMR_warm<-expand.grid(tempTest = seq(min(data.amr.test$tempTest), max(data.amr.test$tempTest) , 1),
                               lnBWg = seq(min(data.amr.test$lnBWg), max(data.amr.test$lnBWg), 0.5))
  data.plotAMR_warm$model_predFE <- predict(best.model.amr.w, re.form = NA, newdata = data.plotAMR_warm)
  # conf.int_AMR_warm <- confint(best.model.amr.w, method = "boot", FUN=function(x)predict(x, data.plotAMR_warm, re.form=NA), nsim=10)
  # data.plotAMR_warm$CI_2.5<-conf.int_AMR_warm[,1]
  # data.plotAMR_warm$CI_97.5<-conf.int_AMR_warm[,2]

  # FAS - warm temp plot grid; plot/predict full range (All data coverage for reference)
  data.plotFAS_warm <- expand.grid(tempTest = seq(min(data.fas.test$tempTest), max(data.fas.test$tempTest) , 1),
                               lnBWg = seq(min(data.fas.test$lnBWg), max(data.fas.test$lnBWg), 0.5))
  data.plotFAS_warm$model_predFE <- predict(best.model.fas.w, re.form = NA, newdata = data.plotFAS_warm)
  # conf.int_FAS_warm <- confint(best.model.fas.w, method = "boot", FUN=function(x)predict(x, data.plotFAS_warm, re.form=NA), nsim=10)
  # data.plotFAS_warm$CI_2.5<-conf.int_FAS_warm[,1]
  # data.plotFAS_warm$CI_97.5<-conf.int_FAS_warm[,2]

  # AS - warm temp plot grid; plot/predict full range (All data coverage for reference)
  data.plotAS_warm <-expand.grid(tempTest = seq(min(data.as.test$tempTest), max(data.as.test$tempTest) , 1),
                               lnBWg = seq(min(data.as.test$lnBWg), max(data.as.test$lnBWg), 0.5))
  data.plotAS_warm$model_predFE <- predict(best.model.as.w, re.form = NA, newdata = data.plotAS_warm)
  # conf.int_AS_warm <- confint(best.model.as.w, method = "boot", FUN=function(x)predict(x, data.plotAS_warm, re.form=NA), nsim=10)
  # data.plotAS_warm$CI_2.5<-conf.int_AS_warm[,1]
  # data.plotAS_warm$CI_97.5<-conf.int_AS_warm[,2]

  # CIs for all model parameters ---------
  CI.amr.ER<-as.data.frame(confint.merMod(best.model.amr.er, level = 0.90, method = "Wald"))
  CI.amr.ER$var<-rownames(CI.amr.ER)
  CI.amr.ER$MR<-"MMR"
  CI.amr.ER$temp_cat<-"ER"

  CI.rmr.ER<-as.data.frame(confint.merMod(best.model.rmr.er, level = 0.90, method = "Wald"))
  CI.rmr.ER$var<-rownames(CI.rmr.ER)
  CI.rmr.ER$MR<-"RMR"
  CI.rmr.ER$temp_cat<-"ER"

  CI.fas.ER<-as.data.frame(confint.merMod(best.model.fas.er, level = 0.90, method = "Wald"))
  CI.fas.ER$var<-rownames(CI.fas.ER)
  CI.fas.ER$MR<-"FAS"
  CI.fas.ER$temp_cat<-"ER"

  CI.as.ER<-as.data.frame(confint.merMod(best.model.as.er, level = 0.90, method = "Wald"))
  CI.as.ER$var<-rownames(CI.as.ER)
  CI.as.ER$MR<-"AS"
  CI.as.ER$temp_cat<-"ER"

  CI.amr.W<-as.data.frame(confint.merMod(best.model.amr.w, level = 0.90, method = "Wald"))
  CI.amr.W$var<-rownames(CI.amr.W)
  CI.amr.W$MR<-"MMR"
  CI.amr.W$temp_cat<-"W"

  CI.rmr.W<-as.data.frame(confint.merMod(best.model.rmr.w, level = 0.90, method = "Wald"))
  CI.rmr.W$var<-rownames(CI.rmr.W)
  CI.rmr.W$MR<-"RMR"
  CI.rmr.W$temp_cat<-"W"

  CI.fas.W<-as.data.frame(confint.merMod(best.model.fas.w, level = 0.90, method = "Wald"))
  CI.fas.W$var<-rownames(CI.fas.W)
  CI.fas.W$MR<-"FAS"
  CI.fas.W$temp_cat<-"W"

  CI.as.W<-as.data.frame(confint.merMod(best.model.as.w, level = 0.90, method = "Wald"))
  CI.as.W$var<-rownames(CI.as.W)
  CI.as.W$MR<-"AS"
  CI.as.W$temp_cat<-"W"

  sum_CItable<-rbind(CI.amr.ER, CI.rmr.ER, CI.fas.ER, CI.as.ER, # ecol relev
                     CI.amr.W, CI.rmr.W, CI.fas.W, CI.as.W) # warm

  # Saving all confidence interval and predictions data frame  ---------
  write.csv(file = paste("./Data_exports/", name_extension,"/Table_CIsummary.csv", sep = ""),
            sum_CItable, row.names=TRUE)

  write.csv(file = paste("./Data_exports/", name_extension,"/dataPred_RMR_er.csv", sep=""),
            data.plotRMR_ER, row.names = F)
  write.csv(file = paste("./Data_exports/", name_extension,"/dataPred_AMR_er.csv", sep=""),
            data.plotAMRint_ER, row.names = F)
  write.csv(file = paste("./Data_exports/", name_extension,"/dataPred_AS_er.csv", sep=""),
            data.plotAS_ER, row.names = F)
  write.csv(file = paste("./Data_exports/", name_extension,"/dataPred_FAS_er.csv", sep=""),
            data.plotFAS_ER, row.names = F)

  write.csv(file = paste("./Data_exports/", name_extension,"/dataPred_RMR_warm.csv", sep=""),
            data.plotRMR_warm, row.names = F)
  write.csv(file = paste("./Data_exports/", name_extension,"/dataPred_AMR_warm.csv", sep=""),
            data.plotAMR_warm, row.names = F)
  write.csv(file = paste("./Data_exports/", name_extension,"/dataPred_AS_warm.csv", sep=""),
            data.plotAS_warm, row.names = F)
  write.csv(file = paste("./Data_exports/", name_extension,"/dataPred_FAS_warm.csv", sep=""),
            data.plotFAS_warm, row.names = F)

  k<-(8.62*10^(-5)) # Boltzmann's constant

  # Summary and figure parameters: ---- 
  summarise.dataframes<- function(data.summarise, export.name){
    dataframe.new <- data.summarise %>%
      dplyr:::summarise(min_bw = min(BW_g), max_bw = max(BW_g), mean_bw = mean(BW_g),
                min_bwLOG = min(lnBWg), max_bwLOG = max(lnBWg), mean_bwLOG = mean(lnBWg),
                min_temp = min(tempTest), max_temp = max(tempTest), mean_temp = mean(tempTest),
                #min_tempArrh = min(tempTestK1000), max_tempArrh = max(tempTestK1000), mean_tempArrh = mean(tempTestK1000),
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

  # scaling parameters ------------------
  # if(grepl(pattern = "lnBWg *", formula(best.model.rmr.er)[2])){ 
  # slope predicitons
  est_scaling_param<-function(best_model,
                              temp_cat,
                              perform_cat){
    # estimate slopes
    est.slopes<-emtrends(best_model,
                         pairwise ~ tempTest,
                         var="lnBWg",
                         at=list(tempTest=seq(5, 30, 1)),
                         mode = "asymptotic") 
    # Intercepts (predicted value at lnBWg = 4.60517; 100 grams)
    est.intercepts <- emmeans(best_model, pairwise ~ tempTest,
                            pbkrtest.limit = 4000,
                            at = list(lnBWg = log(1), tempTest = seq(5, 30, 1)),
                            mode = "asymptotic")
    # Extract summaries
    slopes.df <- as.data.frame(est.slopes$emtrends)[, c("tempTest", "lnBWg.trend", "SE")]
    intercepts.df <- as.data.frame(summary(est.intercepts$emmeans))[, c("tempTest", "emmean", "SE")]
    # Combine in one table
    est.table <- merge(intercepts.df, slopes.df, by = "tempTest")
    colnames(est.table) <- c("Temperature", "Intercept", "SE_int", "Slope", "SE_slope")
    est.table$temp_categ = temp_cat
    est.table$performance = perform_cat
    
    return(est.table)
  }
  
       
  message("Estimate slopes to common size 1 g ")
  # estimate 
  RMR_er_row<-est_scaling_param(best.model.rmr.er, "ecol_relev", "RMR")
  MMR_er_row<-est_scaling_param(best.model.amr.er, "ecol_relev", "MMR")
  AS_er_row<-est_scaling_param(best.model.as.er, "ecol_relev", "AS")
  FAS_er_row<-est_scaling_param(best.model.fas.er, "ecol_relev", "FAS")
  
  RMR_w_row<-est_scaling_param(best.model.rmr.w, "warm", "RMR")
  MMR_w_row<-est_scaling_param(best.model.amr.w, "warm", "MMR")
  AS_w_row<-est_scaling_param(best.model.as.w, "warm", "AS")
  FAS_w_row<-est_scaling_param(best.model.fas.w, "warm", "FAS")
  
  scaling_params<-rbind(MMR_er_row, MMR_w_row, 
                        RMR_er_row, RMR_w_row,
                        AS_er_row, AS_w_row, 
                        FAS_er_row, FAS_w_row)
      
  # save files -----
  # save scaling parameters 
  write.csv(file = paste("./Data_exports/", name_extension,"/scaling_parameters.csv", sep=""),
              scaling_params, row.names = F)
  # save summary data for each metric
  write.csv(file = paste("./Data_exports/", name_extension,"/summary_data_all_metrics.csv", sep=""),
            sum_data, row.names = F)
          

  message("Return: 1) sum_CItable, 2) sum_data, 3) scaling parameters")
  return(list(sum_CItable, sum_data, scaling_params))
}
