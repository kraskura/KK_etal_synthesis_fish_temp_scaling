

get_data_temp <- function(data.amr, data.rmr,
                                  ecology.data,
                                  onlyTop.above = TRUE,
                                  save.FishBase.species.data = FALSE, 
                                  calc_mass_specific = FALSE,
                                  exp_rmr=NULL,
                                  exp_amr=NULL,
                                  exp_as=NULL,
                                  exp_rmr_warm=NULL,
                                  exp_amr_warm=NULL,
                                  exp_as_warm=NULL){
  
  # 1. import the latest AMR file , oct 4 2020 data
  data.amr<-read.csv(data.amr)
  names(data.amr)<-c( "tempAccl","TempAcclDays", "test_category", "tempTest", "fish_ID", "species", "Common_name" , "BW_g", "AMR", "RMR","study_ID", "trial", "trial_ID") 
  
  # 2. import the latest RMR file , oct 4 2020 note 
  data.rmr<-read.csv(data.rmr, header=TRUE, stringsAsFactors=FALSE, fileEncoding="latin1")
  names(data.rmr)<-c( "tempAccl","TempAcclDays", "test_category", "tempTest", "fish_ID",  "species", "Common_name" , "BW_g", "RMR", "study_ID", "trial", "trial_ID") 
  
  # 3. import ecology data
  ecology.data<-read.csv(ecology.data)
  
  # curing AMR *****************************
  data.amr$species<-as.character(data.amr$species)
  # Species names that are different between authors specs and what is on fishbase, here:
  data.amr[data.amr$species ==  "Ophiocephalus argus", "species"] <- "Channa argus"  # source: https://nas.er.usgs.gov/queries/factsheet.aspx?speciesid=2265
  data.amr[data.amr$species ==  "Aristichthys nobilis" , "species"]<- "Hypophthalmichthys nobilis"  # https://nas.er.usgs.gov/queries/FactSheet.aspx?speciesID=551, http://www.iucngisd.org/gisd/species.php?sc=773
  data.amr[data.amr$species ==  "Ctenopharyngodon idellus" , "species"]<- "Ctenopharyngodon idella" # grass carp fishbase.org.
  data.amr[data.amr$species ==  "Leiocassis longirostris" , "species"]<- "Tachysurus dumerili" # https://www-ncbi-nlm-nih-gov.proxy.library.ucsb.edu:9443/Taxonomy/Browser/wwwtax.cgi?id=175787, http://www.catalogueoflife.org/col/details/species/id/af3c089360f1c76ee1a6b204ae4315bc
  data.amr[data.amr$species ==  "Mystus macropterus" , "species"]<- "Hemibagrus macropterus" # https://www.fishbase.se/Nomenclature/SynonymsList.php?ID=52976&SynCode=128451&GenusName=Hemibagrus&SpeciesName=macropterus, https://www-ncbi-nlm-nih-gov.proxy.library.ucsb.edu:9443/Taxonomy/Browser/wwwtax.cgi
  data.amr[data.amr$species == "Onychostoma sima" , "species"]<- "Onychostoma simum" # https://www.fishbase.se/Nomenclature/SynonymsList.php?ID=54868&SynCode=152530&GenusName=Onychostoma&SpeciesName=simum,https://www-ncbi-nlm-nih-gov.proxy.library.ucsb.edu:9443/Taxonomy/Browser/wwwtax.cgi?id=369674&lvl=0
  data.amr[data.amr$species == "Pelteobagrus vachelli" , "species"]<- "Pseudobagrus vachellii" # https://www.fishbase.se/Nomenclature/SynonymsList.php?ID=50950&SynCode=165984&GenusName=Pseudobagrus&SpeciesName=vachellii
  data.amr[data.amr$species == "Sarcocheilichys parvus" , "species"]<- "Sarcocheilichthys parvus" # spelling
  data.amr[data.amr$species == "Oncorhynchus kisutch " , "species"]<- "Oncorhynchus kisutch" # spelling
  data.amr[data.amr$species == "Acrossocheilus monticolus" , "species"]<- "Acrossocheilus monticola" # https://www-ncbi-nlm-nih-gov.proxy.library.ucsb.edu:9443/Taxonomy/Browser/wwwtax.cgi?id=356813
  data.amr[data.amr$species == "Archocentrus nigrofasciatus" , "species"]<- "Amatitlania nigrofasciata" #https://www.fishbase.se/summary/Archocentrus-nigrofasciatus.html, https://nas.er.usgs.gov/queries/factsheet.aspx?SpeciesID=447
  # species("Sinibrama taeniatus") # fishbase is finding this 
  
  # new additions 2022:
  # shorthorn sculpin 
  # zebrafish 
  # Arctic charr
  # Opaleye
  # mahi mahi
  # California killifish

    # curing AMR *****************************
  
  cols.numeric<-c(1,4,8,9,10)
  data.amr[, cols.numeric]<-sapply(data.amr[, cols.numeric], as.numeric)
  
  data.amr$species<-factor(data.amr$species)
  data.amr$lnRMR<-log(data.amr$RMR)
  data.amr$lnAMR<-log(data.amr$AMR)
  data.amr$lnBWg<-log(data.amr$BW_g)
  
  data.amr$trial_ID<-as.factor(data.amr$trial_ID)
  data.amr$study_ID<-as.factor(data.amr$study_ID)
  data.amr$fish_ID<-as.factor(data.amr$fish_ID)
  data.amr$test_category<-factor(data.amr$test_category)
  data.amr$trial<-as.factor(data.amr$trial)
  data.amr$species<-factor(data.amr$species)
  
  # adding Factorial Scope and Aerobic scope to ony AMR dataframe
  data.amr$FAS<-data.amr$AMR/data.amr$RMR
  data.amr$AS<-data.amr$AMR-data.amr$RMR
  data.amr$lnAS<-log(data.amr$AS)
  data.amr$lnFAS<-log(data.amr$FAS)
  data.amr$tempTestK<-celsius.to.kelvin(data.amr$tempTest, round = 2)
  
  # dummy variable for temp category (k-1 variables so 2) 
  # d1(acute): ecol_relev = 0, acute=1, acclim =0 
  # d2(acclim): ecol_relev = 0, acute=0, acclim =1
  
  data.amr$d1<-NA
  data.amr$d2<-NA
  
  data.amr$d1[data.amr$test_category=="ecol_relev"]<-0
  data.amr$d2[data.amr$test_category=="ecol_relev"]<-0
  data.amr$d1[data.amr$test_category=="acute"]<-1
  data.amr$d2[data.amr$test_category=="acute"]<-0
  data.amr$d1[data.amr$test_category=="acclim"]<-0
  data.amr$d2[data.amr$test_category=="acclim"]<-1
  
  # curing RMR ******************************
  data.rmr$species<-as.character(data.rmr$species)
  # Species names that are different between authors specs and what is on fishbase, here:
  data.rmr[data.rmr$species ==  "Ophiocephalus argus", "species"] <- "Channa argus"  # source: https://nas.er.usgs.gov/queries/factsheet.aspx?speciesid=2265
  data.rmr[data.rmr$species ==  "Aristichthys nobilis" , "species"]<- "Hypophthalmichthys nobilis"  # https://nas.er.usgs.gov/queries/FactSheet.aspx?speciesID=551, http://www.iucngisd.org/gisd/species.php?sc=773
  data.rmr[data.rmr$species ==  "Ctenopharyngodon idellus" , "species"]<- "Ctenopharyngodon idella" # grass carp fishbase.org.
  data.rmr[data.rmr$species ==  "Leiocassis longirostris" , "species"]<- "Tachysurus dumerili" # https://www-ncbi-nlm-nih-gov.proxy.library.ucsb.edu:9443/Taxonomy/Browser/wwwtax.cgi?id=175787, http://www.catalogueoflife.org/col/details/species/id/af3c089360f1c76ee1a6b204ae4315bc
  data.rmr[data.rmr$species ==  "Mystus macropterus" , "species"]<- "Hemibagrus macropterus" # https://www.fishbase.se/Nomenclature/SynonymsList.php?ID=52976&SynCode=128451&GenusName=Hemibagrus&SpeciesName=macropterus, https://www-ncbi-nlm-nih-gov.proxy.library.ucsb.edu:9443/Taxonomy/Browser/wwwtax.cgi
  data.rmr[data.rmr$species == "Onychostoma sima" , "species"]<- "Onychostoma simum" # https://www.fishbase.se/Nomenclature/SynonymsList.php?ID=54868&SynCode=152530&GenusName=Onychostoma&SpeciesName=simum,https://www-ncbi-nlm-nih-gov.proxy.library.ucsb.edu:9443/Taxonomy/Browser/wwwtax.cgi?id=369674&lvl=0
  data.rmr[data.rmr$species == "Pelteobagrus vachelli" , "species"]<- "Pseudobagrus vachellii" # https://www.fishbase.se/Nomenclature/SynonymsList.php?ID=50950&SynCode=165984&GenusName=Pseudobagrus&SpeciesName=vachellii
  data.rmr[data.rmr$species == "Sarcocheilichys parvus" , "species"]<- "Sarcocheilichthys parvus" # spelling
  data.rmr[data.rmr$species == "Oncorhynchus kisutch " , "species"]<- "Oncorhynchus kisutch" # spelling
  data.rmr[data.rmr$species == "Acrossocheilus monticolus" , "species"]<- "Acrossocheilus monticola" # https://www-ncbi-nlm-nih-gov.proxy.library.ucsb.edu:9443/Taxonomy/Browser/wwwtax.cgi?id=356813
  data.rmr[data.rmr$species == "Archocentrus nigrofasciatus" , "species"]<- "Amatitlania nigrofasciata" #https://www.fishbase.se/summary/Archocentrus-nigrofasciatus.html, https://nas.er.usgs.gov/queries/factsheet.aspx?SpeciesID=447
  # species("Sinibrama taeniatus") # fishbase is finding this 
  
  # data.rmr[which(is.na(data.rmr[,2])),]
  cols.numeric<-c(1,2,4,8,9)
  
  data.rmr[, cols.numeric]<-sapply(data.rmr[, cols.numeric], as.numeric) 
  message(" NAs present in days acclimated at temp acclimation, when going 'as.numeric'")
  
  data.rmr$log10RMR<-log10(data.rmr$RMR)
  data.rmr$log10BWg<-log10(data.rmr$BW_g)
  data.rmr$lnRMR<-log(data.rmr$RMR)
  data.rmr$lnBWg<-log(data.rmr$BW_g)
  
  data.rmr$trial_ID<-as.factor(data.rmr$trial_ID)
  data.rmr$study_ID<-as.factor(data.rmr$study_ID)
  data.rmr$fish_ID<-as.factor(data.rmr$fish_ID)
  data.rmr$test_category<-factor(data.rmr$test_category)
  data.rmr$trial<-as.factor(data.rmr$trial)
  data.rmr$species<-factor(data.rmr$species)
  data.rmr$tempTestK<-celsius.to.kelvin(data.rmr$tempTest, round = 2)
  
  
  # add dummy variable for temp category (k-1 variables so 2) 
  # d1(acute): ecol_relev = 0, acute=1, acclim =0 
  # d2(acclim): ecol_relev = 0, acute=0, acclim =1
  
  data.rmr$d1<-NA
  data.rmr$d2<-NA
  
  data.rmr$d1[data.rmr$test_category=="ecol_relev"]<-0
  data.rmr$d2[data.rmr$test_category=="ecol_relev"]<-0
  data.rmr$d1[data.rmr$test_category=="acute"]<-1
  data.rmr$d2[data.rmr$test_category=="acute"]<-0
  data.rmr$d1[data.rmr$test_category=="acclim"]<-0
  data.rmr$d2[data.rmr$test_category=="acclim"]<-1
  
  
  if(onlyTop.above){
    # take out fish that are acclimated to Top min ranges.
    # "Gasterosteus aculeatus" acclimated at 10 for warm Populations: "POPMyvW", "POPGTS", "POPAshnW"
    data.rmr<-data.rmr[!c((data.rmr$test_category=="acclim") & data.rmr$species=="Gasterosteus aculeatus" & (data.rmr$trial_ID == "POPMyvW" | data.rmr$trial_ID == "POPGTS" | data.rmr$trial_ID == "POPAshnW") & !is.na(data.rmr$trial_ID)),]
    
    # salmon O. nerka populations that are below 14 c
    data.rmr<-data.rmr[!c((data.rmr$test_category=="acclim") & data.rmr$species=="Oncorhynchus nerka"  & !is.na(data.rmr$trial_ID) & data.rmr$tempTest < 14),]
    
    # "Gasterosteus aculeatus" acclimated at 10 for warm Populations: "POPMyvW", "POPGTS", "POPAshnW"
    data.amr<-data.amr[!c((data.amr$test_category=="acclim") & data.amr$species=="Gasterosteus aculeatus" & (data.amr$trial_ID == "POPMyvW" | data.amr$trial_ID == "POPGTS" | data.amr$trial_ID == "POPAshnW") & !is.na(data.amr$trial_ID)),]
    
    # salmon O. nerka populations that are below 14 c
    data.amr<-data.amr[!c((data.amr$test_category=="acclim") & data.amr$species=="Oncorhynchus nerka"  & !is.na(data.amr$trial_ID) & data.amr$tempTest < 14),]
    
  }
  
  # Formatting/organizing/combining AMR and RMR datasets --------------
  # all(is.numeric(data.amr$AMR))
  # all(!is.na(data.amr$AMR)) # making sure all AMR are numeric values, no NA's
  # 
  # getting the needed data from RMR 
  data.rmr.combine<-data.rmr[, c( "test_category", "tempTest", "fish_ID",  "species", "Common_name" , "BW_g", "RMR","study_ID", "trial", "tempTestK")] 
  colnames(data.rmr.combine)<-c( "test_category", "tempTest", "fish_ID",  "species", "Common_name" , "BW_g", "MR","study_ID", "trial", "tempTestK")
  data.rmr.combine$MR_type<-as.numeric(0)
  data.rmr.combine$MR_type2<-"RMR"
  
  # getting the needed data from AMR 
  data.amr.combine<-data.amr[, c( "test_category", "tempTest", "fish_ID",  "species", "Common_name" , "BW_g", "AMR","study_ID", "trial", "tempTestK")] 
  colnames(data.amr.combine)<-c( "test_category", "tempTest", "fish_ID",  "species", "Common_name" , "BW_g", "MR","study_ID", "trial", "tempTestK")
  data.amr.combine$MR_type<-as.numeric(1)
  data.amr.combine$MR_type2<-"AMR"
  
  nrow(data.amr.combine) # 
  nrow(data.rmr.combine) # 
  
  dataMR<-rbind(data.amr.combine, data.rmr.combine) # combining RMR and AMR in long format 
  dataMR$lnMR<-log(dataMR$MR)
  dataMR$lnBWg<-log(dataMR$BW_g)
  
  category_names <- list('acclim'="Acclimated", 'acute'="Acute", 'ecol_relev'="Ecologically relevant")
  category_labeller <- function(variable,value){
    return(category_names[value])
  }
  
  
  k<-(8.62*10^(-5)) # Boltzmann's constant
  # kBoltz<-1.38 * 10^23 
  E<-0.63 # activation energy MTE
  data.amr$tempTestK1<-1/data.amr$tempTestK
  data.amr$tempTestK1000<-1000/data.amr$tempTestK
  
  data.rmr$tempTestK1<-1/data.rmr$tempTestK
  data.rmr$tempTestK1000<-1000/data.rmr$tempTestK
  
  data.amr$tempTestK1_1000 <- data.amr$tempTestK1*1000
  data.rmr$tempTestK1_1000 <- data.rmr$tempTestK1*1000
  dataMR$tempTestK1<-1/dataMR$tempTestK
  dataMR$tempTestK1_1000 <- dataMR$tempTestK1*1000
  
  data.amr$tempTestkT<-1/(data.amr$tempTestK*k)
  data.rmr$tempTestkT<-1/(data.rmr$tempTestK*k)
  dataMR$tempTestkT<-1/(dataMR$tempTestK*k)
  
  
  
  # adding body shape from R FishBase online: 
  if(save.FishBase.species.data){
    
    colnames(dataMR)[4]<-"Species"
    mySpecies<-as.character(unique(dataMR$Species))
    mySpecies
    species_data<-as.data.frame(species(mySpecies))
    
    species_data$BodyShapeI <- as.factor(gsub(' / ', '-', as.character(species_data$BodyShapeI)))
    
    names(species_data)[names(species_data) == "Species"] <- "species"
    names(dataMR)[names(dataMR) == "Species"] <- "species"
    
    filename<-paste("./Data_exports/Species/Species_data_FishBase_", Sys.Date(), ".csv", sep="")
    write.csv(x = species_data, file = filename, row.names = FALSE)
    
  }

  
  # potentially useful variables taken out: i) length female out, ii)  Weight, iii) "WeightFemale, iv) PriceCateg,
  # ***************************
  # merging data files
  # 1. Body shap, demersPelag, Vulnerability Price categories, used for aquaculture
  names(dataMR)[names(dataMR) == "Species"] <- "species"
  names(ecology.data)<-c("species", "DemersPelag", "BodyShapeI", "Climate", "salintyComb")
  dataMR<- merge(x = dataMR, y = ecology.data , by = "species", all.x = TRUE)
  data.amr<- merge(x = data.amr, y = ecology.data, by = "species", all.x = TRUE)
  data.rmr<-merge(x = data.rmr, y = ecology.data, by = "species", all.x = TRUE)
  
  # cleaning 
  data.as<-data.amr[!c(is.na(data.amr$AS) | is.infinite(data.amr$lnAS)), ]
  data.fas<-data.amr[!c(is.na(data.amr$FAS) | is.infinite(data.amr$lnFAS)), ]
  
  data.as<-data.as[!c(is.na(data.as$AS) | is.infinite(data.as$lnAS)), ]
  data.fas<-data.fas[!c(is.na(data.fas$FAS)), ]
  
  # caculate mass specific metabolic rates
  if(calc_mass_specific){
  
    # adjust mass specific MMR and RMR; ecol relevant 
    data.amr$mass_specamr<-(data.amr$AMR/(data.amr$BW_g^exp_amr)) # mass specific; linear  
    data.rmr$mass_specrmr<-data.rmr$RMR/(data.rmr$BW_g^exp_rmr) # mass specific ; linear
    # log transformed
    data.amr$ln_mass_specamr<-log(data.amr$mass_specamr) # 
    data.rmr$ln_mass_specrmr<-log(data.rmr$mass_specrmr) #
    
    # adjust warm: AMR/MMR and RMR
    data.amr$mass_specamr[!c(data.amr$test_category == "ecol_relev")]<-(data.amr$AMR[!c(data.amr$test_category == "ecol_relev")]/(data.amr$BW_g[!c(data.amr$test_category == "ecol_relev")]^exp_amr_warm)) # 
    data.rmr$mass_specrmr[!c(data.rmr$test_category == "ecol_relev")]<-data.rmr$RMR[!c(data.rmr$test_category == "ecol_relev")]/(data.rmr$BW_g[!c(data.rmr$test_category == "ecol_relev")]^exp_rmr_warm) # 
  
    
    # aerobic scope dataset; ecolo relevant 
    data.as$mass_specamr<-(data.as$AMR/(data.as$BW_g^exp_amr)) # 
    data.as$mass_specrmr<-data.as$RMR/(data.as$BW_g^exp_rmr) # 
    data.as$mass_specas<-data.as$AS/(data.as$BW_g^exp_as) # 
    
    # aerobic scope dataset; warm
    data.as$mass_specamr[!c(data.as$test_category == "ecol_relev")] <- (data.as$AMR[!c(data.as$test_category == "ecol_relev")]/(data.as$BW_g[!c(data.as$test_category == "ecol_relev")]^exp_amr_warm)) # 
    data.as$mass_specrmr[!c(data.as$test_category == "ecol_relev")] <- (data.as$RMR[!c(data.as$test_category == "ecol_relev")]/(data.as$BW_g[!c(data.as$test_category == "ecol_relev")]^exp_rmr_warm)) #
    data.as$mass_specas[!c(data.as$test_category == "ecol_relev")] <- (data.as$AS[!c(data.as$test_category == "ecol_relev")]/(data.as$BW_g[!c(data.as$test_category == "ecol_relev")]^exp_as_warm)) # 
    
    
    # add the exponents to the dataset 
    data.amr$exp_amr<-exp_amr
    data.rmr$exp_rmr<-exp_rmr
    data.as$exp_amr<-exp_amr
    data.as$exp_rmr<-exp_rmr
    data.as$exp_as<-exp_as
    data.amr$exp_amr[!c(data.amr$test_category == "ecol_relev")]<-exp_amr_warm
    data.rmr$exp_rmr[!c(data.rmr$test_category == "ecol_relev")]<-exp_rmr_warm
    data.as$exp_amr[!c(data.as$test_category == "ecol_relev")]<-exp_amr_warm
    data.as$exp_rmr[!c(data.as$test_category == "ecol_relev")]<-exp_rmr_warm
    data.as$exp_as[!c(data.as$test_category == "ecol_relev")]<-exp_as_warm # poly term
    
    message(paste("Mass specific measures: 1) AMR (b = ", exp_amr, "); 2) RMR (b = ",exp_rmr,")", sep = ""))
    message(paste("Mass specific measures for WARM: 1) AMR (b = ", exp_amr_warm, "); 2) RMR (b = ",exp_rmr_warm,")", sep = ""))
  
  }else{
    data.amr$mass_specamr<-NA
    data.rmr$mass_specrmr<-NA
    data.amr$ln_mass_specamr<-NA
    data.rmr$ln_mass_specrmr<-NA
    data.as$mass_specamr<-NA
    data.as$mass_specrmr<-NA
    data.as$mass_specas<-NA
   
    data.amr$exp_amr<-exp_amr
    data.rmr$exp_rmr<-exp_rmr
    data.as$exp_amr<-exp_amr
    data.as$exp_rmr<-exp_rmr
    data.as$exp_as<-exp_as
    data.amr$exp_amr[!c(data.amr$test_category == "ecol_relev")]<-NA
    data.rmr$exp_rmr[!c(data.rmr$test_category == "ecol_relev")]<-NA
    data.as$exp_amr[!c(data.as$test_category == "ecol_relev")]<-NA
    data.as$exp_rmr[!c(data.as$test_category == "ecol_relev")]<-NA
    data.as$exp_as[!c(data.as$test_category == "ecol_relev")]<-NA
    
    message("Mass specific values not calculated")

  }

  # ------------- finalising data -----------
  # create data frames from each group:
  data.amrAC<-data.amr[data.amr$test_category=="acute",]
  data.rmrAC<-data.rmr[data.rmr$test_category=="acute",]
  data.amrAM<-data.amr[data.amr$test_category=="acclim",]
  data.rmrAM<-data.rmr[data.rmr$test_category=="acclim",]
  data.amrER<-data.amr[data.amr$test_category=="ecol_relev",]
  data.rmrER<-data.rmr[data.rmr$test_category=="ecol_relev",]
  data.asAC<-data.as[data.as$test_category=="acute",]
  data.fasAC<-data.fas[data.fas$test_category=="acute",]
  data.asAM<-data.as[data.as$test_category=="acclim",]
  data.fasAM<-data.fas[data.fas$test_category=="acclim",]
  data.asER<-data.as[data.as$test_category=="ecol_relev",]
  data.fasER<-data.fas[data.fas$test_category=="ecol_relev",]
  
  data.amrAC$test_category<-factor(data.amrAC$test_category)
  data.rmrAC$test_category<-factor(data.rmrAC$test_category)
  data.amrAM$test_category<-factor(data.amrAM$test_category)
  data.rmrAM$test_category<-factor(data.rmrAM$test_category)
  data.amrER$test_category<-factor(data.amrER$test_category)
  data.rmrER$test_category<-factor(data.rmrER$test_category)
  data.asAC$test_category<-factor(data.asAC$test_category)
  data.fasAC$test_category<-factor(data.fasAC$test_category)
  data.asAM$test_category<-factor(data.asAM$test_category)
  data.fasAM$test_category<-factor(data.fasAM$test_category)
  data.asER$test_category<-factor(data.asER$test_category)
  data.fasER$test_category<-factor(data.fasER$test_category)
  
  
  # return(data.amrMean)
  # return(data.rmrMean)
  data.list<-list(data.amrAC, 
                  data.rmrAC,
                  data.amrAM,
                  data.rmrAM,
                  data.amrER,
                  data.rmrER,
                  
                  data.asAC, 
                  data.fasAC,
                  data.asAM,
                  data.fasAM,
                  data.asER,
                  data.fasER,
                  
                  data.amr, 
                  data.rmr, 
                  dataMR,
                  data.as,
                  data.fas
                  )
  return(data.list)
 
}
