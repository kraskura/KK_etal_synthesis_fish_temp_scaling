# Model scaling parameters and CIs ----------

model_outputs<-function(phylo=TRUE,
                        best.model.rmr.er,
                        best.model.amr.er,
                        best.model.as.er,
                        best.model.fas.er,
                        best.model.rmr.w,
                        best.model.amr.w,
                        best.model.as.w,
                        best.model.fas.w){
  if(phylo){
    foldername.phylo<-"Phylo"
  }else{
    foldername.phylo<-"nonPhylo"
  }
  
  if(!file.exists(paste("./Data_exports/", foldername.phylo, "/Table_CIsummary.csv", sep=""))){
    dir.create(paste("./Data_exports/", foldername.phylo, sep =""), recursive =TRUE)

    # RMR optimal - ecologically relevant grid; plot/predict full range (All data coverage for reference)
    data.plotRMR_ER<-expand.grid(tempTestK1000 = seq(C40inTempTestK1000, C0inTempTestK1000 , 0.05), lnBWg = seq(sum_rmr_ALL$min_bwLOG,sum_rmr_ALL$max_bwLOG, 0.5))
    data.plotRMR_ER$model_predFE <- predict(best.model.rmr.er, re.form = NA, newdata = data.plotRMR_ER)
    conf.int_RMR_ER <- confint(best.model.rmr.er, method = "boot", FUN=function(x)predict(x, data.plotRMR_ER, re.form=NA), nsim=10)
    data.plotRMR_ER$CI_2.5<-conf.int_RMR_ER[,1]
    data.plotRMR_ER$CI_97.5<-conf.int_RMR_ER[,2]
  
    # FAS optimal - ecologically relevant grid; plot/predict full range (All data coverage for reference)
    data.plotFAS_ER<-expand.grid(tempTest = 20, lnBWg = seq(sum_fas_ER$min_bwLOG,sum_fas_ER$max_bwLOG, 0.5))
    data.plotFAS_ER$model_predFE <- predict(best.model.fas.er, re.form = NA, newdata = data.plotFAS_ER)
    conf.int_FAS_ER <- confint(best.model.fas.er, method = "boot", FUN=function(x)predict(x, data.plotFAS_ER, re.form=NA), nsim=10)
    data.plotFAS_ER$CI_2.5<-conf.int_FAS_ER[,1]
    data.plotFAS_ER$CI_97.5<-conf.int_FAS_ER[,2]
  
    # AS optimal - ecologically relevant grid; plot/predict full range (All data coverage for reference)
    data.plotAS_ER<-expand.grid(tempTestK1000 = seq(C40inTempTestK1000, C0inTempTestK1000 , 0.05), lnBWg = seq(sum_as_ALL$min_bwLOG,sum_as_ALL$max_bwLOG, 0.5))
    data.plotAS_ER$model_predFE <- predict(best.model.as.er, re.form = NA, newdata = data.plotAS_ER)
    conf.int_AS_ER <- confint(best.model.as.er, method = "boot", FUN=function(x)predict(x, data.plotAS_ER, re.form=NA), nsim=10)
    data.plotAS_ER$CI_2.5<-conf.int_AS_ER[,1]
    data.plotAS_ER$CI_97.5<-conf.int_AS_ER[,2]
  
    # MMR or AMR  optimal - interaction model; plot/predict full range (All data coverage for reference)
    data.plotAMRint_ER <- expand.grid(tempTestK1000 = seq(C40inTempTestK1000,C0inTempTestK1000, 0.05), lnBWg= seq(sum_amr_ALL$min_bwLOG,sum_amr_ALL$max_bwLOG, 0.5))
    data.plotAMRint_ER$model_predFE <- predict(best.model.amr.er, re.form = NA, newdata = data.plotAMRint_ER)
    conf.int_AMRint_ER <- confint(best.model.amr.er, method = "boot", FUN=function(x)predict(x, data.plotAMRint_ER, re.form=NA), nsim=10)
    data.plotAMRint_ER$CI_2.5<-conf.int_AMRint_ER[,1]
    data.plotAMRint_ER$CI_97.5<-conf.int_AMRint_ER[,2]
    data.plotAMRint_ER$tempTestK1000_inC <- ((1000/data.plotAMRint_ER$tempTestK1000))-275.15
  
    # All temps for the Arrhenius plot:
    # RMR - warm temp grid; plot/predict full range (All data coverage for reference)
    data.plotRMR_warm<-expand.grid(tempTestK1000 = seq(C40inTempTestK1000, C0inTempTestK1000 , 0.05), lnBWg = seq(sum_rmr_ALL$min_bwLOG,sum_rmr_ALL$max_bwLOG, 0.5))
    data.plotRMR_warm$model_predFE <- predict(best.model.rmr.w, re.form = NA, newdata = data.plotRMR_warm)
    conf.int_RMR_warm <- confint(best.model.rmr.w, method = "boot", FUN=function(x)predict(x, data.plotRMR_warm, re.form=NA), nsim=10)
    data.plotRMR_warm$CI_2.5<-conf.int_RMR_warm[,1]
    data.plotRMR_warm$CI_97.5<-conf.int_RMR_warm[,2]
  
    # MMR - warm temp plot grid; plot/predict full range (All data coverage for reference)
    data.plotAMR_warm<-expand.grid(tempTestK1000 = seq(C40inTempTestK1000, C0inTempTestK1000 , 0.05), lnBWg = seq(sum_amr_ALL$min_bwLOG,sum_amr_ALL$max_bwLOG, 0.5))
    data.plotAMR_warm$model_predFE <- predict(best.model.amr.w, re.form = NA, newdata = data.plotAMR_warm)
    conf.int_AMR_warm <- confint(best.model.amr.w, method = "boot", FUN=function(x)predict(x, data.plotAMR_warm, re.form=NA), nsim=10)
    data.plotAMR_warm$CI_2.5<-conf.int_AMR_warm[,1]
    data.plotAMR_warm$CI_97.5<-conf.int_AMR_warm[,2]
  
    # FAS - warm temp plot grid; plot/predict full range (All data coverage for reference)
    data.plotFAS_warm <- expand.grid(lnBWg = seq(sum_fas_ALL$min_bwLOG,sum_fas_ALL$max_bwLOG, 0.5), tempTest = 20)
    data.plotFAS_warm$model_predFE <- predict(best.model.fas.w, re.form = NA, newdata = data.plotFAS_warm)
    conf.int_FAS_warm <- confint(best.model.fas.w, method = "boot", FUN=function(x)predict(x, data.plotFAS_warm, re.form=NA), nsim=10)
    data.plotFAS_warm$CI_2.5<-conf.int_FAS_warm[,1]
    data.plotFAS_warm$CI_97.5<-conf.int_FAS_warm[,2]
  
    # AS - warm temp plot grid; plot/predict full range (All data coverage for reference)
    data.plotAS_warm <- expand.grid(lnBWg = seq(sum_as_ALL$min_bwLOG,sum_as_ALL$max_bwLOG, 0.5), tempTestK1000 = seq(C40inTempTestK1000, C0inTempTestK1000 , 0.1))
    data.plotAS_warm$model_predFE <- predict(best.model.as.w, re.form = NA, newdata = data.plotAS_warm)
    conf.int_AS_warm <- confint(best.model.as.w, method = "boot", FUN=function(x)predict(x, data.plotAS_warm, re.form=NA), nsim=10)
    data.plotAS_warm$CI_2.5<-conf.int_AS_warm[,1]
    data.plotAS_warm$CI_97.5<-conf.int_AS_warm[,2]
  
    # CIs for all model parameters
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
  
    sum_CItable<<-rbind(CI.amr.ER, CI.rmr.ER, CI.fas.ER, CI.as.ER, # ecol relev
                       CI.amr.W, CI.rmr.W, CI.fas.W, CI.as.W) # warm
  
    # Saving all export files
    write.csv(file = paste("./Data_exports/", foldername.phylo,"/Table_CIsummary.csv", sep = ""), sum_CItable, row.names=TRUE)
  
    write.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_RMR_er.csv", sep=""),  data.plotRMR_ER, row.names = F)
    write.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_AMR_er.csv", sep=""),  data.plotAMRint_ER, row.names = F)
    write.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_AS_er.csv", sep=""),  data.plotAS_ER, row.names = F)
    write.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_FAS_er.csv", sep=""),  data.plotFAS_ER, row.names = F)
  
    write.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_RMR_warm.csv", sep=""),  data.plotRMR_warm, row.names = F)
    write.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_AMR_warm.csv", sep=""),  data.plotAMR_warm, row.names = F)
    write.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_AS_warm.csv", sep=""),  data.plotAS_warm, row.names = F)
    write.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_FAS_warm.csv", sep=""),  data.plotFAS_warm, row.names = F)


  }else{
    message("Data imported from './Data_exports/'")
    
    sum_CItable<<-read.csv(file = paste("./Data_exports/", foldername.phylo,"/Table_CIsummary.csv", sep =""))
  
    data.plotRMR_ER<<-read.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_RMR_er.csv", sep=""))
    data.plotAMRint_ER<<-read.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_AMR_er.csv", sep=""))
    data.plotAS_ER<<-read.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_AS_er.csv", sep=""))
    data.plotFAS_ER<<-read.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_FAS_er.csv", sep=""))
  
    data.plotRMR_warm<<-read.csv(paste("./Data_exports/", foldername.phylo,"/dataPred_RMR_warm.csv", sep=""))
    data.plotAMR_warm<<-read.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_AMR_warm.csv", sep=""))
    data.plotAS_warm<<-read.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_AS_warm.csv", sep=""))
    data.plotFAS_warm<<-read.csv(file = paste("./Data_exports/", foldername.phylo,"/dataPred_FAS_warm.csv", sep=""))
  }
 
  k<-(8.62*10^(-5)) # Boltzmann's constant

  # ECOL RELEV
  RMR_slope<<-round(fixef(best.model.rmr.er)[2],3)
  RMR_int<<-round(fixef(best.model.rmr.er)[1],3)
  FAS_slope<<-round(fixef(best.model.fas.er)[2],3)
  FAS_int<<-round(fixef(best.model.fas.er)[1],3)
  AS_slope<<-round(fixef(best.model.as.er)[2],3)
  AS_int<<-round(fixef(best.model.as.er)[1],3)

  # WARM
  MMR_slope_w<<-round(fixef(best.model.amr.w)[2],3)
  MMR_int_w<<-round(fixef(best.model.amr.w)[1],3)
  RMR_slope_w<<-round(fixef(best.model.rmr.w)[2],3)
  RMR_int_w<<-round(fixef(best.model.rmr.w)[1],3)
  FAS_slope_w<<-round(fixef(best.model.fas.w)[2],3)
  FAS_int_w<<-round(fixef(best.model.fas.w)[1],3)
  AS_slope_w<<-round(fixef(best.model.as.w)[2],3)
  AS_int_w<<-round(fixef(best.model.as.w)[1],3)

  C20inTempTestK1000<<-1000/(((20+273.15)))
  C10inTempTestK1000<<-1000/((10+273.15))
  C0inTempTestK1000<<-1000/((0+273.15))
  C30inTempTestK1000<<-1000/((30+273.15))
  C40inTempTestK1000<<-1000/((40+273.15))

  # predicted slopes for MMR interaction
  AMR.slopes<<-(emtrends(best.model.amr.er, pairwise ~ tempTestK1000,
                        pbkrtest.limit = 4000, var="lnBWg",
                        at=list(tempTestK1000=c(C0inTempTestK1000,
                                                C10inTempTestK1000,
                                                C20inTempTestK1000,
                                                C30inTempTestK1000,
                                                C40inTempTestK1000)))) # 0, 10, 20, 30, 40 C
  AMR.slopes<<-as.data.frame(AMR.slopes$emtrends)
  AMR.slopes$tempTestK1000_inC<<-((1000/AMR.slopes$tempTestK1000))-273.15
  AMR.slopes$performance<<-"MMR"
  AMR.slopes$temp_categ<<-"er"
  AMR.slopes$`(Intercept)` <<- NA
  AMR.slopes$tempTestK1000 <<- NA
  AMR.slopes$`E(ev)`<<-NA
  # fixed effects activation energies:
  MMR_E_ER<<-(emtrends(best.model.amr.er, pairwise ~ lnBWg , pbkrtest.limit = 4000,
                      var = "tempTestK1000",
                      at = list(lnBWg = c(log(1),log(10), log(100),log(1000)))))
  MMR_E_ER<<-as.data.frame(MMR_E_ER$emtrends)
  MMR_E_ER$performance<<-"MMR"
  MMR_E_ER$temp_categ<<-"er"
  MMR_E_ER$lnBWg.trend<<-NA
  MMR_E_ER$`(Intercept)` <<- NA
  MMR_E_ER$tempTestK1000_inC <<- NA

  RMR_E_ER<<-round(fixef(best.model.rmr.er)[3],3)
  AS_E_ER<<-round(fixef(best.model.as.er)[3],3)
  MMR_E_W<<-round(fixef(best.model.amr.w)[3],3)
  RMR_E_W<<-round(fixef(best.model.rmr.w)[3],3)
  AS_E_W<<-round(fixef(best.model.as.w)[3],3)

  # activation energies:
  # -1*(E*k)*1000
  MMR_E_ER$MMR_E_ER_eV<<--1*(MMR_E_ER$tempTestK1000.trend)*k*1000
  RMR_E_ER_eV<<--1*(RMR_E_ER)*k*1000
  AS_E_ER_eV<<--1*(AS_E_ER)*k*1000

  MMR_E_W_eV<<--1*(MMR_E_W)*k*1000
  RMR_E_W_eV<<--1*(RMR_E_W)*k*1000
  AS_E_W_eV<<--1*(AS_E_W)*k*1000

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

  scaling.params<<-rbind(MMR_er_row, MMR_er_row2, MMR_W_row, 
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
  sum_data<<-rbind(sum_amr_ER, sum_amr_warm, 
                  sum_rmr_ER, sum_rmr_warm,
                  sum_as_ER, sum_as_warm,
                  sum_fas_ER, sum_fas_warm)
  
  sum_data$MR<<-c("MMR", "MMR", "RMR", "RMR", "AS", "AS", "FAS", "FAS")
  sum_data$Temp<<-c("er", "warm","er", "warm","er", "warm","er", "warm")

  # 20 deg
  # -1*(RMR_E*k)*1000
  # tempTestK1000
  
  # needed values
  C20inTempTestK1000<-1000/(((20+273.15)))
  C10inTempTestK1000<-1000/((10+273.15))
  C0inTempTestK1000<-1000/((0+273.15))
  C30inTempTestK1000<-1000/((30+273.15))
  C40inTempTestK1000<-1000/((40+273.15))

  # reset datasets with mass specific values:
  data.list<-get_data_temp(data.amr = "./Data/Fish_AMR_temp_dataset_mar2022.csv",
                                   data.rmr = "./Data/Fish_RMR_temp_dataset_mar2022.csv",
                                   ecology.data = "./Data/Kraskura_species_ecologies_mar2022.csv", 
                                   onlyTop.above = TRUE,
                                   calc_mass_specific = TRUE,
                                   exp_rmr = RMR_slope,
                                   exp_amr = round(AMR.slopes$lnBWg.trend[3],3), # at 20ºC
                                   exp_as = AS_slope,
                                   exp_rmr_warm = RMR_slope_w,
                                   exp_amr_warm = MMR_slope_w,
                                   exp_as_warm = AS_slope_w)
  
  data.amrAC<<-data.frame(data.list[1]) 
  data.rmrAC<<-data.frame(data.list[2])
  data.amrAM<<-data.frame(data.list[3])
  data.rmrAM<<-data.frame(data.list[4])
  data.amrER<<-data.frame(data.list[5])
  data.rmrER<<-data.frame(data.list[6])
  
  data.asAC<<-data.frame(data.list[7]) 
  data.fasAC<<-data.frame(data.list[8])
  data.asAM<<-data.frame(data.list[9])
  data.fasAM<<-data.frame(data.list[10])
  data.asER<<-data.frame(data.list[11])
  data.fasER<<-data.frame(data.list[12])
  
  data.amr<<-data.frame(data.list[13]) 
  data.rmr<<-data.frame(data.list[14]) 
  dataMR<<-data.frame(data.list[15])
  data.as<<-data.frame(data.list[16])
  data.fas<<-data.frame(data.list[17])
  
  data.amr$test_category2<<-"acclim"
  data.amr$test_category2[data.amr$test_category=="acute"] <<- "acute"
  data.amr$test_category3<<-"ecol_relev"
  data.amr$test_category3[data.amr$test_category=="acute"] <<- "warm"
  data.amr$test_category3[data.amr$test_category=="acclim"] <<- "warm"
  data.rmr$test_category2<<-"acclim"
  data.rmr$test_category2[data.rmr$test_category=="acute"] <<- "acute"
  data.rmr$test_category3<<-"ecol_relev"
  data.rmr$test_category3[data.rmr$test_category=="acute"] <<- "warm"
  data.rmr$test_category3[data.rmr$test_category=="acclim"] <<- "warm"
  data.fas$test_category2<<-"acclim"
  data.fas$test_category2[data.fas$test_category=="acute"] <<- "acute"
  data.fas$test_category3<<-"ecol_relev"
  data.fas$test_category3[data.fas$test_category=="acute"] <<- "warm"
  data.fas$test_category3[data.fas$test_category=="acclim"] <<- "warm"
  data.as$test_category2<<-"acclim"
  data.as$test_category2[data.as$test_category=="acute"] <<- "acute"
  data.as$test_category3<<-"ecol_relev"
  data.as$test_category3[data.as$test_category=="acute"] <<- "warm"
  data.as$test_category3[data.as$test_category=="acclim"] <<- "warm"
  
  # acclimated and warm are > Topt temperatures
  data.amr.test<<-rbind(data.amrAC, data.amrAM)
  data.rmr.test<<-rbind(data.rmrAC, data.rmrAM)
  data.fas.test<<-data.amr.test[c(!is.na(data.amr.test$FAS) & is.finite(data.amr.test$FAS)) , ]
  data.as.test<<-data.amr.test[c(!is.na(data.amr.test$lnAS) & is.finite(data.amr.test$lnAS)) , ]
  
  data.amrER$tempTestK1000_inC<<-((1000/data.amrER$tempTestK1000))-273.15
  
}
