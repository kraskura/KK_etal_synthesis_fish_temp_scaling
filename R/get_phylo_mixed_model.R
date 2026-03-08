# Author: Krista Kraskura
# Description: 
#   - phylogentic models 
#   - make life stage specific main text figures
#   - analyze data for all lifestages together
# ********************************************
# ********************************************
# library(here)


# *************************************************
get_data_temp <- function(data.amr.readin,
                          data.rmr.readin,
                          ecology.data.readin,
                          onlyTop.above = TRUE,
                          save.FishBase.species.data = FALSE, 
                          calc_mass_specific = FALSE,
                          exp_rmr=NULL,
                          exp_amr=NULL,
                          exp_as=NULL,
                          exp_rmr_warm=NULL,
                          exp_amr_warm=NULL,
                          exp_as_warm=NULL, 
                          exclude_Wootton = FALSE){
  
  if(exclude_Wootton){
    message("Wootton et al data on zebrafish exculded")
  }
  
  # 1. import the latest AMR file , oct 4 2020 data
  data.amr<-read.csv(data.amr.readin)
  names(data.amr)<-c( "tempAccl","TempAcclDays", "test_category", "tempTest", "fish_ID", "species", "Common_name" , "BW_g", "AMR", "RMR","study_ID", "trial", "trial_ID", "lifestage", "chamber_vol_L", "h_starved", "MMR_method", "h_in_respo") 
  
  # 2. import the latest RMR file , oct 4 2020 note 
  data.rmr<-read.csv(data.rmr.readin, header=TRUE, stringsAsFactors=FALSE, fileEncoding="latin1")
  names(data.rmr)<-c( "tempAccl","TempAcclDays", "test_category", "tempTest", "fish_ID",  "species", "Common_name" , "BW_g", "RMR", "study_ID", "trial", "trial_ID","lifestage", "chamber_vol_L", "h_starved", "h_in_respo") 
  
  # 3. import ecology data
  ecology.data<-read.csv(ecology.data.readin)
  
  if(exclude_Wootton){
    data.amr<-data.amr %>% 
      filter(study_ID != "415")
    # view(data.amr)
    data.rmr<-data.rmr %>% 
      filter(study_ID != "415")
  }

  
  # ******************************************  
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
  # data.amr[data.amr$species == "Aphanius iberus" , "species"]<- "Apricaphanius iberus" # https://www.fishbase.se/summary/Apricaphanius-iberus # dont change because the 

  # new additions 2022: ******
  # shorthorn sculpin 
  # zebrafish 
  # Arctic charr
  # Opaleye
  # mahi mahi
  # California killifish
  # 
  # new additions 2026: ******
  # spanish toothcarp
  # Arc-eye hawkfish
  # blacktail snapper
  # convict tang
  # green sturgeon
  # golden perch
  # roman seabream
  # westslope cutthroat trout
  # pacific cod
  # barred surfperch
  # emerald rockcod
  # three-striped dwarf cichlid
  # cardinal tetra
  
  # str(data.amr)
  cols.numeric<-c(1,4,8,9,10, 15, 16, 18)
  data.amr[, cols.numeric]<-sapply(data.amr[, cols.numeric], as.numeric)
  
  data.amr$species<-factor(data.amr$species)
  data.amr$lnRMR<-log(data.amr$RMR)
  data.amr$lnAMR<-log(data.amr$AMR)
  data.amr$lnBWg<-log(data.amr$BW_g)
  
  # adding Factorial Scope and Aerobic scope to ony AMR dataframe
  data.amr$FAS<-data.amr$AMR/data.amr$RMR
  data.amr$AS<-data.amr$AMR-data.amr$RMR
  data.amr$lnAS<-log(data.amr$AS)
  data.amr$lnFAS<-log(data.amr$FAS)
  data.amr$tempTestK<-celsius.to.kelvin(data.amr$tempTest, round = 2)
  
  data.amr$trial_ID<-as.factor(data.amr$trial_ID)
  data.amr$study_ID<-as.factor(data.amr$study_ID)
  data.amr$fish_ID<-as.factor(data.amr$fish_ID)
  data.amr$test_category<-factor(data.amr$test_category)
  data.amr$trial<-as.factor(data.amr$trial)
  data.amr$species<-factor(data.amr$species)
  data.amr$MMR_method<-factor(data.amr$MMR_method) 
  data.amr$lifestage<-factor(data.amr$lifestage) 
   
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
  
  
  # ******************************************  
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
  data.rmr[data.rmr$species == "Aphanius iberus" , "species"]<- "Apricaphanius iberus" # https://www.fishbase.se/summary/Apricaphanius-iberus
    
  # data.rmr[which(is.na(data.rmr[,2])),]
  cols.numeric<-c(1,2,4,8,9, 14, 15, 16)
  
  data.rmr[, cols.numeric]<-sapply(data.rmr[, cols.numeric], as.numeric) 
  message(" NAs present in RMR dataset when numeric variables have NA (days acclimated, h in respo, h starved, respo chamber volume")
  
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
  data.rmr$lifestage<-factor(data.rmr$lifestage) 
  
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
  
  # ********************************************************************
  # Take out MMR from study # 415:
  # Wootton, H.F., Morrongiello, J.R., Schmitt, T., Audzijonyte, A., 2022 Smaller adult fish size in warmer water is not explained by elevated metabolism. Ecology Letters. https://doi.org/10.1111/ele.13989
  # Measurements produce FAS > 20; only RMR assumed usable. 
  data.amr<-data.amr[!data.amr$study_ID == 415,]
  message("EXCLUDING Study # 415 MMR measurements; Wootton et al 2022; FAS > 20")
  # ********************************************************************
  
  # dont include data that were below optimal temperature range. 
  if(onlyTop.above){
    # take out fish that are acclimated to Top min ranges.
    # 
  
    # RMR *****
    # "Gasterosteus aculeatus" acclimated at 10 for warm Populations: "POPMyvW", "POPGTS", "POPAshnW"
    data.rmr<-data.rmr[!c((data.rmr$test_category=="acclim") & data.rmr$species=="Gasterosteus aculeatus" & (data.rmr$trial_ID == "POPMyvW" | data.rmr$trial_ID == "POPGTS" | data.rmr$trial_ID == "POPAshnW") & !is.na(data.rmr$trial_ID)),]
    
    # salmon O. nerka populations that are below 14 c
    data.rmr<-data.rmr[!c((data.rmr$test_category=="acclim") & data.rmr$species=="Oncorhynchus nerka"  & !is.na(data.rmr$trial_ID) & data.rmr$tempTest < 14),]
    
    # Bailey et al take out 10 C roman seabream
    data.rmr<-data.rmr[!c(data.rmr$study_ID=="1019" & data.rmr$tempTest == 10),]
    # Zillig et al juvenile chinook
    data.rmr<-data.rmr[!c(data.rmr$study_ID=="1007" & data.rmr$tempAccl < 16 & data.rmr$tempTest < 16),]
    # Zillig et al  sturgeon 
    data.rmr<-data.rmr[!c(data.rmr$study_ID=="1013" & data.rmr$tempAccl < 19 & data.rmr$tempTest < 15),]
    # Kraskura barred surfperch tested at 12C 
    data.rmr<-data.rmr[!c((data.rmr$species=="Amphistichus argenteus") & data.rmr$tempTest < 16),]
    
    # AMR *****
    # "Gasterosteus aculeatus" acclimated at 10 for warm Populations: "POPMyvW", "POPGTS", "POPAshnW"
    data.amr<-data.amr[!c((data.amr$test_category=="acclim") & data.amr$species=="Gasterosteus aculeatus" & (data.amr$trial_ID == "POPMyvW" | data.amr$trial_ID == "POPGTS" | data.amr$trial_ID == "POPAshnW") & !is.na(data.amr$trial_ID)),]
    
    # salmon O. nerka populations that are below 14 c
    data.amr<-data.amr[!c((data.amr$test_category=="acclim") & data.amr$species=="Oncorhynchus nerka"  & !is.na(data.amr$trial_ID) & data.amr$tempTest < 14),]
    
    # Bailey et al take out 10 C roman seabream
    data.amr<-data.amr[!c(data.amr$study_ID=="1019" & data.amr$tempTest == 10),]
    # Zillig et al juvenile chinook
    data.amr<-data.amr[!c(data.amr$study_ID=="1007" & data.amr$tempAccl < 16 & data.amr$tempTest < 16),]
    # Zillig et al  sturgeon 
    data.amr<-data.amr[!c(data.amr$study_ID=="1013" & data.amr$tempAccl < 19 & data.amr$tempTest < 16),]
    # Kraskura barred surfperch tested at 12C 
    data.amr<-data.amr[!c((data.amr$species=="Amphistichus argenteus") & data.amr$tempTest < 16),]

    
  }
  
  # Formatting/organizing/combining AMR and RMR datasets *************************
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
  
  nrow(data.amr.combine) # 5476 feb 2026
  nrow(data.rmr.combine) # 7427 feb 2026
  
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
    
    require(rfishbase)
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
  data.as<-data.amr[!c(is.na(data.amr$AS) | is.infinite(data.amr$lnAS)), ] # no NA for AS
  data.fas<-data.amr[!c(is.na(data.amr$FAS) | is.infinite(data.amr$lnFAS)), ] # no NA for FAS
  
  data.as<-data.as[!c(is.na(data.as$AS) | is.infinite(data.as$lnAS)), ] # no NA fod AS 
  data.fas<-data.fas[!c(is.na(data.fas$FAS)), ] # no NA for FAS 
  
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
    
    message("Mass specific values not calculated; provide exponents to do that")

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

# ***************************************************************
#  ------------------- SUPPORTING FUNCTIONS ----------------------
# ***************************************************************
# function to get phylogenetic relatedness matrix for each data subset
get_phylo_matrix<-function(species.list,
                           matrix.name,
                           dataset.ID,
                           plot = TRUE){
  
  taxon_search <- tnrs_match_names(names=species.list, context_name="Vertebrates") 
  ott_in_tree <- ott_id(taxon_search)[is_in_tree(ott_id(taxon_search))]
  tr <- tol_induced_subtree(ott_ids = ott_in_tree)
  
  tr$tip.label <- strip_ott_ids(tr$tip.label, remove_underscores = TRUE)
  labels <- as.data.frame(tr$tip.label)

  # any unmatching tip labels with our data? change those 
  labels$`tr$tip.label`[which(labels$`tr$tip.label` == "Oncorhynchus mykiss (species in domain Eukaryota)")]<-"Oncorhynchus mykiss"
  labels$`tr$tip.label`[which(labels$`tr$tip.label` == "Gadus morhua (species in domain Eukaryota)")]<-"Gadus morhua"
  labels$`tr$tip.label`[which(labels$`tr$tip.label` == "Leiocassis longirostris")]<-"Tachysurus dumerili" 
  labels$`tr$tip.label`[which(labels$`tr$tip.label` == "Tachysurus vachellii")]<-"Pseudobagrus vachellii" 
  labels$`tr$tip.label`[which(labels$`tr$tip.label` == "Rhinogobius similis")]<- "Rhinogobius giurinus" 
  
  labels$`tr$tip.label`[which(labels$`tr$tip.label` == "Rhinogobius similis")]<- "Rhinogobius giurinus" 
  
  # (sort(labels$`tr$tip.label`) == sort(species.list)) "Amphistichus argenteus" is the same too. 

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
  write.csv(as.data.frame.matrix(A), here(paste("Data_exports/Phylo_matrices/", matrix.name,".csv", sep="")))

  if(plot){
    A.plot<-ggtree(tree) +
    geom_tiplab(as_ylab=TRUE, color='black', size = 12, align = TRUE)
    ggsave(plot = A.plot,
           filename = here(paste("Data_exports/Phylo_plots/",dataset.ID, "_treeplot.png", sep="")),
           width = 7,
           height = 15)
  } 

  assign(matrix.name, value = A, envir = .GlobalEnv)
  # assign(tree.name, value = tree, envir = .GlobalEnv)

}

# ***************************************************************
# # function to order model selection based on the lowest BIC score
BICdelta<-function(BICtable){
  BIC.t <- BICtable [order(BICtable$BIC), ]
  BIC.t$delta <- round(abs(BIC.t$BIC[1] -  BIC.t$BIC), 5)
  return(BIC.t)
}

# **************************************************
# making sure the correct model is extracted; created with help from Cloude Code 
check_best_model <- function(model,
                             datatable) {
  # Extract name from lmer model
  model_name <- deparse(substitute(model))
  # Extract first rowname from datatable
  best_model_name <- rownames(datatable)[1]
  # Compare and write message
  if (model_name == best_model_name) {
    message("Correct best models identified: ", best_model_name)
  } else {
    warning("Mismatch! Model name '", model_name, "' does not match best model '", best_model_name, "'")
  }
}

# **************************************************
# identifies if the model has interaction; created with help from Cloude Code 
has_interaction <- function(model) {
  # Extract only the fixed-effects formula (drops all random effect terms)
  # lme4::fixef formula is obtained by removing the random parts
  fixed_formula <- nobars(formula(model))  # nobars() from lme4 strips (|) terms

  # extract term labels — only fixed effects remain
  fixed_terms <- attr(terms(fixed_formula), "term.labels")

  # print("Fixed-effect terms detected:\n")
  # print(fixed_terms)
  # print("\n")

  # Interactions are represented with ":" after formula expansion
  any(grepl(":", fixed_terms, fixed = TRUE))
}

# **************************************************
# run appropriate Anova
run_anova <- function(model,
                      save.path) {

  interaction_detected <- has_interaction(model)

  if (interaction_detected) {
    cat(">>> Interaction term(s) detected — using Type III ANOVA\n\n")
    anova_result <- car::Anova(model, type = "III")
    anova_type   <- "III"
  } else {
    cat(">>> No interaction terms detected — using Type II ANOVA\n\n")
    anova_result <- car::Anova(model, type = "II")
    anova_type   <- "II"
  }

  # Export formatted text report with Anovas, model summary, etc
  report_path <- save.path
    sink(report_path)
    
    cat("============================================================\n")
    cat("      ANOVA RESULTS REPORT — lmer model (car package)\n")
    cat("============================================================\n\n")
    
    cat("Model formula:\n")
    print(formula(model))
    
    cat(paste0("\nANOVA Type selected automatically: Type ", anova_type, "\n"))

    cat(" **************** ANOVA Table ****************  \n")
    print(anova_result)
    
    cat("\n ****************  Model Summary **************** \n")
    print(summary(model))
    
    sink()
}

# **************************************************
# if(grepl(pattern = "lnBWg *", formula(best.model.rmr.er)[2])){ 
# slope predicitons
est_scaling_param<-function(best_model,
                            temp_cat,
                            perform_cat){
    # estimate slopes
    est.slopes<-emtrends(best_model,
                         pairwise ~ tempTest,
                         var="lnBWg",
                         at=list(tempTest=seq(5, 35, 1)),
                         mode = "asymptotic") 
    # Intercepts (predicted value at 1 gram)
    message("Estimate slopes to common size 1 g ")

    est.intercepts <- emmeans(best_model, pairwise ~ tempTest,
                            pbkrtest.limit = 4000,
                            at = list(lnBWg = log(1), tempTest = seq(5, 35, 1)),
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

# *************************************************
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
                        data.rmr.test.modpar,
                        data.rmrER.modpar, 
                        data.amr.test.modpar,
                        data.amrER.modpar,                        
                        data.as.test.modpar,
                        data.asER.modpar, 
                        data.fas.test.modpar,
                        data.fasER.modpar, 
                        name_extension){
  
  # dir.create(paste("./Data_exports/", name_extension, sep =""), recursive =TRUE)

  # RMR optimal - ecologically relevant grid; plot/predict full range (All data coverage for reference)
  data.plotRMR_ER<-expand.grid(tempTest = seq(min(data.rmrER.modpar$tempTest),
                                              max(data.rmrER.modpar$tempTest) , 1),
                               lnBWg = seq(min(data.rmrER.modpar$lnBWg),
                                           max(data.rmrER.modpar$lnBWg), 0.5))
  data.plotRMR_ER$model_predFE <- predict(best.model.rmr.er, re.form = NA, newdata = data.plotRMR_ER)


  # FAS optimal - ecologically relevant grid; plot/predict full range (All data coverage for reference)
  data.plotFAS_ER<-expand.grid(tempTest = seq(min(data.fasER.modpar$tempTest),
                                              max(data.fasER.modpar$tempTest) , 1),
                               lnBWg = seq(min(data.fasER.modpar$lnBWg), 
                                           max(data.fasER.modpar$lnBWg), 0.5))
  data.plotFAS_ER$model_predFE <- predict(best.model.fas.er, re.form = NA, newdata = data.plotFAS_ER)


  # AS optimal - ecologically relevant grid; plot/predict full range (All data coverage for reference)
  data.plotAS_ER<-expand.grid(tempTest = seq(min(data.asER.modpar$tempTest),
                                             max(data.asER.modpar$tempTest) , 1),
                               lnBWg = seq(min(data.asER.modpar$lnBWg),
                                           max(data.asER.modpar$lnBWg), 0.5))
  data.plotAS_ER$model_predFE <- predict(best.model.as.er, re.form = NA, newdata = data.plotAS_ER)


  # MMR or AMR  optimal - interaction model; plot/predict full range (All data coverage for reference)
  data.plotAMR_ER <- expand.grid(tempTest = seq(min(data.amrER.modpar$tempTest),
                                                   max(data.amrER.modpar$tempTest) , 1),
                               lnBWg = seq(min(data.amrER.modpar$lnBWg),
                                           max(data.amrER.modpar$lnBWg), 0.5))
  data.plotAMR_ER$model_predFE <- predict(best.model.amr.er, re.form = NA, newdata = data.plotAMR_ER)


  # All temps :
  # RMR - warm temp grid; plot/predict full range (All data coverage for reference)
  data.plotRMR_warm<-expand.grid(tempTest = seq(min(data.rmr.test.modpar$tempTest),
                                                max(data.rmr.test.modpar$tempTest) , 1),
                               lnBWg = seq(min(data.rmr.test.modpar$lnBWg),
                                           max(data.rmr.test.modpar$lnBWg), 0.5))
  data.plotRMR_warm$model_predFE <- predict(best.model.rmr.w, re.form = NA, newdata = data.plotRMR_warm)

  # MMR - warm temp plot grid; plot/predict full range (All data coverage for reference)
  data.plotAMR_warm<-expand.grid(tempTest = seq(min(data.amr.test.modpar$tempTest),
                                                max(data.amr.test.modpar$tempTest) , 1),
                               lnBWg = seq(min(data.amr.test.modpar$lnBWg),
                                           max(data.amr.test.modpar$lnBWg), 0.5))
  data.plotAMR_warm$model_predFE <- predict(best.model.amr.w, re.form = NA, newdata = data.plotAMR_warm)


  # FAS - warm temp plot grid; plot/predict full range (All data coverage for reference)
  data.plotFAS_warm <- expand.grid(tempTest = seq(min(data.fas.test.modpar$tempTest),
                                                  max(data.fas.test.modpar$tempTest) , 1),
                               lnBWg = seq(min(data.fas.test.modpar$lnBWg),
                                           max(data.fas.test.modpar$lnBWg), 0.5))
  data.plotFAS_warm$model_predFE <- predict(best.model.fas.w, re.form = NA, newdata = data.plotFAS_warm)

  # AS - warm temp plot grid; plot/predict full range (All data coverage for reference)
  data.plotAS_warm <-expand.grid(tempTest = seq(min(data.as.test.modpar$tempTest),
                                                max(data.as.test.modpar$tempTest) , 1),
                               lnBWg = seq(min(data.as.test.modpar$lnBWg), 
                                           max(data.as.test.modpar$lnBWg), 0.5))
  data.plotAS_warm$model_predFE <- predict(best.model.as.w, re.form = NA, newdata = data.plotAS_warm)

  # CIs for all model parameters ---------
  CI.amr.ER<-as.data.frame(confint.merMod(best.model.amr.er, level = 0.95, method = "Wald"))
  CI.amr.ER$var<-rownames(CI.amr.ER)
  CI.amr.ER$MR<-"MMR"
  CI.amr.ER$temp_cat<-"ER"

  CI.rmr.ER<-as.data.frame(confint.merMod(best.model.rmr.er, level = 0.95, method = "Wald"))
  CI.rmr.ER$var<-rownames(CI.rmr.ER)
  CI.rmr.ER$MR<-"RMR"
  CI.rmr.ER$temp_cat<-"ER"

  CI.fas.ER<-as.data.frame(confint.merMod(best.model.fas.er, level = 0.95, method = "Wald"))
  CI.fas.ER$var<-rownames(CI.fas.ER)
  CI.fas.ER$MR<-"FAS"
  CI.fas.ER$temp_cat<-"ER"

  CI.as.ER<-as.data.frame(confint.merMod(best.model.as.er, level = 0.95, method = "Wald"))
  CI.as.ER$var<-rownames(CI.as.ER)
  CI.as.ER$MR<-"AS"
  CI.as.ER$temp_cat<-"ER"

  CI.amr.W<-as.data.frame(confint.merMod(best.model.amr.w, level = 0.95, method = "Wald"))
  CI.amr.W$var<-rownames(CI.amr.W)
  CI.amr.W$MR<-"MMR"
  CI.amr.W$temp_cat<-"W"

  CI.rmr.W<-as.data.frame(confint.merMod(best.model.rmr.w, level = 0.95, method = "Wald"))
  CI.rmr.W$var<-rownames(CI.rmr.W)
  CI.rmr.W$MR<-"RMR"
  CI.rmr.W$temp_cat<-"W"

  CI.fas.W<-as.data.frame(confint.merMod(best.model.fas.w, level = 0.95, method = "Wald"))
  CI.fas.W$var<-rownames(CI.fas.W)
  CI.fas.W$MR<-"FAS"
  CI.fas.W$temp_cat<-"W"

  CI.as.W<-as.data.frame(confint.merMod(best.model.as.w, level = 0.95, method = "Wald"))
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
            data.plotAMR_ER, row.names = F)
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

  sum_amr_ER<-summarise.dataframes(data.amrER.modpar, "sum_amr_ER")
  sum_rmr_ER<-summarise.dataframes(data.rmrER.modpar, "sum_rmr_ER")
  sum_fas_ER<-summarise.dataframes(data.fasER.modpar, "sum_fas_ER")
  sum_as_ER<-summarise.dataframes(data.asER.modpar, "sum_as_ER")
  sum_amr_warm<-summarise.dataframes(data.amr.test.modpar, "sum_amr_warm")
  sum_rmr_warm<-summarise.dataframes(data.rmr.test.modpar, "sum_rmr_warm")
  sum_fas_warm<-summarise.dataframes(data.fas.test.modpar, "sum_fas_warm")
  sum_as_warm<-summarise.dataframes(data.as.test.modpar, "sum_as_warm")
  
  # row-bind in one dataframe 
  sum_data<-rbind(sum_amr_ER, sum_amr_warm, 
                  sum_rmr_ER, sum_rmr_warm,
                  sum_as_ER, sum_as_warm,
                  sum_fas_ER, sum_fas_warm)
  
  sum_data$MR<-c("MMR", "MMR", "RMR", "RMR", "AS", "AS", "FAS", "FAS")
  sum_data$Temp<-c("er", "warm","er", "warm","er", "warm","er", "warm")

  # scaling parameters ------------------
  
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



# ***************************************************************
#  ------------------- DATA ANALYSIS FUNCTIONS ----------------------
# ***************************************************************
# lifestage = "juvenile"
# lifestage = "adult" 

get_phylo_mixed_models<-function(data.rmr.test.modrun,
                                data.rmrER.modrun, 
                                data.amr.test.modrun,
                                data.amrER.modrun,                        
                                data.as.test.modrun,
                                data.asER.modrun, 
                                data.fas.test.modrun,
                                data.fasER.modrun, 
                                lifestage = NULL){ #default to all data
  
# Phylo models -----------------
  # get correct data frames
  if (!is.null(lifestage)) {
    #  filter data to only have  specific lifestage
    data.rmrER<-data.rmrER.modrun[data.rmrER.modrun$lifestage == lifestage &
                             !is.na(data.rmrER.modrun$lifestage == lifestage), ]
    data.rmr.test<-data.rmr.test.modrun[data.rmr.test.modrun$lifestage == lifestage &
                                   !is.na(data.rmr.test.modrun$lifestage == lifestage), ]
    data.amrER<-data.amrER.modrun[data.amrER.modrun$lifestage == lifestage &
                             !is.na(data.amrER.modrun$lifestage == lifestage), ]
    data.amr.test<-data.amr.test.modrun[data.amr.test.modrun$lifestage == lifestage &
                                   !is.na(data.amr.test.modrun$lifestage == lifestage), ]
    data.asER<-data.asER.modrun[data.asER.modrun$lifestage == lifestage &
                           !is.na(data.asER.modrun$lifestage == lifestage), ]
    data.as.test<-data.as.test.modrun[data.as.test.modrun$lifestage == lifestage &
                                 !is.na(data.as.test.modrun$lifestage == lifestage), ]
    data.fasER<-data.fasER.modrun[data.fasER.modrun$lifestage == lifestage &
                             !is.na(data.fasER.modrun$lifestage == lifestage), ]
    data.fas.test<-data.fas.test.modrun[data.fas.test.modrun$lifestage == lifestage &
                                   !is.na(data.fas.test.modrun$lifestage == lifestage), ]
    
  }else{
    data.rmrER<-data.rmrER.modrun
    data.rmr.test<-data.rmr.test.modrun
    data.amrER<-data.amrER.modrun
    data.amr.test<-data.amr.test.modrun
    data.asER<-data.asER.modrun
    data.as.test<-data.as.test.modrun
    data.fasER<-data.fasER.modrun
    data.fas.test<-data.fas.test.modrun
  }

  # get phylo trees make sure data is clean:
  data.rmrER<-droplevels(data.rmrER)
  data.rmr.test<-droplevels(data.rmr.test)
  data.amrER<-droplevels(data.amrER)
  data.amr.test<-droplevels(data.amr.test)
  data.asER<-droplevels(data.asER)
  data.as.test<-droplevels(data.as.test)
  data.fasER<-droplevels(data.fasER)
  data.fas.test<-droplevels(data.fas.test)


  # this calls custom function 'get_data_phylo_matrix.R'. 
  # IMPORTANT: ensure to receive a message: "All species names are identified and mathced with phylo data N species" that marks that all species have been matched with original dataset. 
  if(!is.null(lifestage)){
    if (lifestage == "adult") {
    message("Running adult data")
    # rmr
    get_phylo_matrix(species.list = unique(levels(data.rmrER$species)),
                     matrix.name = "A.adult",
                     dataset.ID = "RMR-optimal-adult")
    get_phylo_matrix(species.list = unique(levels(data.rmr.test$species)),
                     matrix.name = "A.adult.rmr.w",
                     dataset.ID = "RMR-warm-adult")
    # amr 
    get_phylo_matrix(species.list = unique(levels(data.amrER$species)),
                     matrix.name = "A.adult.mmr.er",
                     dataset.ID = "MMR-optimal-adult")
    get_phylo_matrix(species.list = unique(levels(data.amr.test$species)),
                     matrix.name = "A.adult.mmr.w",
                     dataset.ID = "MMR-warm-adult")
    # aas
    get_phylo_matrix(species.list = unique(levels(data.asER$species)),
                     matrix.name = "A.adult.aas.er",
                     dataset.ID = "AAS-optimal-adult")
    get_phylo_matrix(species.list = unique(levels(data.as.test$species)),
                     matrix.name = "A.adult.aas.w",
                     dataset.ID = "AAS-warm-adult")
    # fas
    get_phylo_matrix(species.list = unique(levels(data.fasER$species)),
                     matrix.name = "A.adult.fas.er",
                     dataset.ID = "FAS-optimal-adult")
    get_phylo_matrix(species.list = unique(levels(data.fas.test$species)),
                     matrix.name = "A.adult.fas.w",
                     dataset.ID = "FAS-warm-adult")
    
    
    ## RMR optimal -----------------
    Phylo_RMR_model2 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A.adult))
    Phylo_RMR_model2int <- Almer(lnRMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A.adult))
    
    Phylo_RMR_model4 <- Almer(lnRMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A.adult))
    Phylo_RMR_model4int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER,, REML=FALSE, A = list(species = A.adult))
    
    Phylo_RMR_model5 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A.adult))
    Phylo_RMR_model5int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A.adult))
    
    
    # No singular fits 
    ### BIC rmr optimal --------------
    RMR_BIC<-BICdelta(BIC(
    
                    Phylo_RMR_model2, Phylo_RMR_model2int, #Phylo_RMR_model2.POLY, Phylo_RMR_model2intPOLY,
                    Phylo_RMR_model4, Phylo_RMR_model4int,# Phylo_RMR_model4.POLY, Phylo_RMR_model4intPOLY,
                    Phylo_RMR_model5, Phylo_RMR_model5int# Phylo_RMR_model5.POLY, Phylo_RMR_model5intPOLY
    ))
    
    ## MMR optimal -----------------
    Phylo_MMR_model2 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.adult.mmr.er))
    Phylo_MMR_model2int <- Almer(lnAMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.adult.mmr.er))
    
    Phylo_MMR_model4 <- Almer(lnAMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.adult.mmr.er))
    Phylo_MMR_model4int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.adult.mmr.er))
    
    Phylo_MMR_model5 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.adult.mmr.er))
    Phylo_MMR_model5int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.adult.mmr.er))
    
    ### BIC --------------
    MMR_BIC<-BICdelta(BIC(
                    Phylo_MMR_model2, Phylo_MMR_model2int, #Phylo_MMR_model2.POLY, Phylo_MMR_model2intPOLY,
                    Phylo_MMR_model4, Phylo_MMR_model4int, #Phylo_MMR_model4.POLY, Phylo_MMR_model4intPOLY,
                    Phylo_MMR_model5, Phylo_MMR_model5int))#, Phylo_MMR_model5.POLY, Phylo_MMR_model5intPOLY))
    
    ## AAS optimal ------------------
    #
    Phylo_AS_model2 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.adult.aas.er))
    Phylo_AS_model2int <- Almer(lnAS ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.adult.aas.er))
    
    Phylo_AS_model4 <- Almer(lnAS ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.adult.aas.er))
    Phylo_AS_model4int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.adult.aas.er))
    
    Phylo_AS_model5 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.adult.aas.er))
    Phylo_AS_model5int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.adult.aas.er))
    
    
    ### BIC --------
    AS_BIC<-BICdelta(BIC(
    
                    Phylo_AS_model2, Phylo_AS_model2int,# Phylo_AS_model2.POLY, Phylo_AS_model2intPOLY,
                    Phylo_AS_model4, Phylo_AS_model4int, #Phylo_AS_model4.POLY, Phylo_AS_model4intPOLY,
                    Phylo_AS_model5, Phylo_AS_model5int))#, #Phylo_AS_model5.POLY, Phylo_AS_model5intPOLY))
    
    
    ## FAS optimal ----------------
    # 
    Phylo_FAS_model2 <- Almer(log(FAS) ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.adult.fas.er))
    Phylo_FAS_model2int <- Almer(log(FAS) ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.adult.fas.er))
    
    Phylo_FAS_model4 <- Almer(log(FAS) ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.adult.fas.er))
    Phylo_FAS_model4int <- Almer(log(FAS) ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.adult.fas.er))
    
    Phylo_FAS_model5 <- Almer(log(FAS) ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.adult.fas.er))
    Phylo_FAS_model5int <- Almer(log(FAS) ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.adult.fas.er))
    
    
    ### BIC -------------
    FAS_BIC<-BICdelta(BIC(
                    Phylo_FAS_model2, Phylo_FAS_model2int,# Phylo_FAS_model2.POLY, Phylo_FAS_model2intPOLY,
                    Phylo_FAS_model4, Phylo_FAS_model4int,# Phylo_FAS_model4.POLY, Phylo_FAS_model4intPOLY,
                    Phylo_FAS_model5, Phylo_FAS_model5int))#, Phylo_FAS_model5.POLY, Phylo_FAS_model5intPOLY))
    
    
    ## RMR warm ------------------
    # 
    Phylo_RMR_W_model2 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.adult.rmr.w))
    Phylo_RMR_W_model2int <- Almer(lnRMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.adult.rmr.w))
    
    Phylo_RMR_W_model4 <- Almer(lnRMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.adult.rmr.w))
    Phylo_RMR_W_model4int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.adult.rmr.w))
    
    Phylo_RMR_W_model5 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.adult.rmr.w))
    Phylo_RMR_W_model5int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.adult.rmr.w))
    
    
    
    ## BIC -------------
    RMR_W_BIC<-BICdelta(BIC(
                    Phylo_RMR_W_model2, Phylo_RMR_W_model2int,# Phylo_RMR_W_model2.POLY, Phylo_RMR_W_model2intPOLY,
                    Phylo_RMR_W_model4, Phylo_RMR_W_model4int,# Phylo_RMR_W_model4.POLY, Phylo_RMR_W_model4intPOLY,
                    Phylo_RMR_W_model5, Phylo_RMR_W_model5int))#, Phylo_RMR_W_model5.POLY, Phylo_RMR_W_model5intPOLY))
    
    # AMR / warm temps -------------
    Phylo_MMR_W_model2 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.adult.mmr.w))
    Phylo_MMR_W_model2int <- Almer(lnAMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.adult.mmr.w))
    
    Phylo_MMR_W_model4 <- Almer(lnAMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.adult.mmr.w))
    Phylo_MMR_W_model4int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.adult.mmr.w))
    
    Phylo_MMR_W_model5 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.adult.mmr.w))
    Phylo_MMR_W_model5int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.adult.mmr.w))
    
    
    ## BIC -------------
    MMR_W_BIC<-BICdelta(BIC(
                    Phylo_MMR_W_model2, Phylo_MMR_W_model2int, #Phylo_MMR_W_model2.POLY, Phylo_MMR_W_model2intPOLY,
                    Phylo_MMR_W_model4, Phylo_MMR_W_model4int, #Phylo_MMR_W_model4.POLY, Phylo_MMR_W_model4intPOLY,
                    Phylo_MMR_W_model5, Phylo_MMR_W_model5int))#, Phylo_MMR_W_model5.POLY, Phylo_MMR_W_model5intPOLY))
    
    
    
    # AS / warm temps --------------
    Phylo_AS_W_model2 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.adult.aas.w))
    Phylo_AS_W_model2int <- Almer(lnAS ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.adult.aas.w))
    
    Phylo_AS_W_model4 <- Almer(lnAS ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.adult.aas.w))
    Phylo_AS_W_model4int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.adult.aas.w))
    
    Phylo_AS_W_model5 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.adult.aas.w))
    Phylo_AS_W_model5int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.adult.aas.w))
    
    ## BIC -------------
    AS_W_BIC<-BICdelta(BIC(
                    Phylo_AS_W_model2, Phylo_AS_W_model2int,# Phylo_AS_W_model2.POLY, Phylo_AS_W_model2intPOLY,
                    Phylo_AS_W_model4, Phylo_AS_W_model4int,# Phylo_AS_W_model4.POLY, Phylo_AS_W_model4intPOLY,
                    Phylo_AS_W_model5, Phylo_AS_W_model5int))#, Phylo_AS_W_model5.POLY, Phylo_AS_W_model5intPOLY))
    
    # FAS / warm temps --------------
    Phylo_FAS_W_model2 <- Almer(lnFAS ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.adult.fas.w))
    Phylo_FAS_W_model2int <- Almer(lnFAS ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.adult.fas.w))
    
    Phylo_FAS_W_model4 <- Almer(lnFAS ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.adult.fas.w))
    Phylo_FAS_W_model4int <- Almer(lnFAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.adult.fas.w))
    
    Phylo_FAS_W_model5 <- Almer(lnFAS ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.adult.fas.w))
    Phylo_FAS_W_model5int <- Almer(lnFAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.adult.fas.w))
    
    
    ## BIC -------------
    FAS_W_BIC<-BICdelta(BIC(
                    Phylo_FAS_W_model2, Phylo_FAS_W_model2int, #Phylo_FAS_W_model2.POLY, Phylo_FAS_W_model2intPOLY,
                    Phylo_FAS_W_model4, Phylo_FAS_W_model4int, #Phylo_FAS_W_model4.POLY, Phylo_FAS_W_model4intPOLY,
                    Phylo_FAS_W_model5, Phylo_FAS_W_model5int)) #Phylo_FAS_W_model5.POLY, Phylo_FAS_W_model5intPOLY))
    
    
    ## BIC results ------
    # RMR_BIC # Phylo_RMR_model4 (dBIC 0)
    # MMR_BIC # Phylo_MMR_model4int (dBIC 0)
    # AS_BIC  # Phylo_AS_model4 (dBIC 0)
    # FAS_BIC # Phylo_FAS_model2 (dBIC 0)
    # 
    # RMR_W_BIC # Phylo_RMR_W_model5 (dBIC 0)
    # MMR_W_BIC # Phylo_MMR_W_model5 (dBIC 0)
    # AS_W_BIC  # Phylo_AS_W_model5int (dBIC 0)
    # FAS_W_BIC # Phylo_FAS_W_model5 (dBIC 0)
    
    # jan 2026 best 
    # ecol relevant
    rmr_mod_ER<- Phylo_RMR_model4
    amr_mod_ER<-Phylo_MMR_model4int
    as_mod_ER<- Phylo_AS_model4
    fas_mod_ER<- Phylo_FAS_model2
  
    # warm
    rmr_mod_W<-Phylo_RMR_W_model5
    amr_mod_W<-Phylo_MMR_W_model5
    as_mod_W<- Phylo_AS_W_model5int
    fas_mod_W<-Phylo_FAS_W_model5
    
    check_best_model(Phylo_RMR_model4, RMR_BIC)
    check_best_model(Phylo_MMR_model4int, MMR_BIC)
    check_best_model(Phylo_AS_model4, AS_BIC)
    check_best_model(Phylo_FAS_model2, FAS_BIC)
    
    check_best_model(Phylo_RMR_W_model5, RMR_W_BIC)
    check_best_model(Phylo_MMR_W_model5, MMR_W_BIC)
    check_best_model(Phylo_AS_W_model5int, AS_W_BIC)
    check_best_model(Phylo_FAS_W_model5, FAS_W_BIC)
    
    write.csv(file = here("Data_exports/BICs/RMR_BIC_adult.csv"), x = RMR_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/MMR_BIC_adult.csv"), x = MMR_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/AS_BIC_adult.csv"), x = AS_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/FAS_BIC_adult.csv"), x = FAS_BIC, row.names = T)
    
    write.csv(file = here("Data_exports/BICs/RMR_W_BIC_adult.csv"), x = RMR_W_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/MMR_W_BIC_adult.csv"), x = MMR_W_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/AS_W_BIC_adult.csv"), x = AS_W_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/FAS_W_BIC_adult.csv"), x = FAS_W_BIC, row.names = T)
    
    }
    
    if(lifestage == "juvenile"){
    message("Running juvenile data")
    # rmr
    get_phylo_matrix(species.list = unique(levels(data.rmrER$species)),
                     matrix.name = "A.juvenile",
                     dataset.ID = "RMR-optimal-juvenile")
    get_phylo_matrix(species.list = unique(levels(data.rmr.test$species)),
                     matrix.name = "A.juvenile.rmr.w",
                     dataset.ID = "RMR-warm-juvenile")
    # amr 
    get_phylo_matrix(species.list = unique(levels(data.amrER$species)),
                     matrix.name = "A.juvenile.mmr.er",
                     dataset.ID = "MMR-optimal-juvenile")
    get_phylo_matrix(species.list = unique(levels(data.amr.test$species)),
                     matrix.name = "A.juvenile.mmr.w",
                     dataset.ID = "MMR-warm-juvenile")
    # aas
    get_phylo_matrix(species.list = unique(levels(data.asER$species)),
                     matrix.name = "A.juvenile.aas.er",
                     dataset.ID = "AAS-optimal-juvenile")
    get_phylo_matrix(species.list = unique(levels(data.as.test$species)),
                     matrix.name = "A.juvenile.aas.w",
                     dataset.ID = "AAS-warm-juvenile")
    # fas
    get_phylo_matrix(species.list = unique(levels(data.fasER$species)),
                     matrix.name = "A.juvenile.fas.er",
                     dataset.ID = "FAS-optimal-juvenile")
    get_phylo_matrix(species.list = unique(levels(data.fas.test$species)),
                     matrix.name = "A.juvenile.fas.w",
                     dataset.ID = "FAS-warm-juvenile")
    
    
    ## RMR optimal -----------------
    Phylo_RMR_model2 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A.juvenile))
    Phylo_RMR_model2int <- Almer(lnRMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A.juvenile))
    
    Phylo_RMR_model4 <- Almer(lnRMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A.juvenile))
    Phylo_RMR_model4int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A.juvenile))
    
    Phylo_RMR_model5 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A.juvenile))
    Phylo_RMR_model5int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A.juvenile))
    
    
    # No singular fits 
    ### BIC rmr optimal --------------
    RMR_BIC<-BICdelta(BIC(
    
                    Phylo_RMR_model2, Phylo_RMR_model2int, #Phylo_RMR_model2.POLY, Phylo_RMR_model2intPOLY,
                    Phylo_RMR_model4, Phylo_RMR_model4int,# Phylo_RMR_model4.POLY, Phylo_RMR_model4intPOLY,
                    Phylo_RMR_model5, Phylo_RMR_model5int# Phylo_RMR_model5.POLY, Phylo_RMR_model5intPOLY
    ))
    
    ## MMR optimal -----------------
    Phylo_MMR_model2 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.juvenile.mmr.er))
    Phylo_MMR_model2int <- Almer(lnAMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.juvenile.mmr.er))
    
    Phylo_MMR_model4 <- Almer(lnAMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.juvenile.mmr.er))
    Phylo_MMR_model4int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.juvenile.mmr.er))
    
    Phylo_MMR_model5 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.juvenile.mmr.er))
    Phylo_MMR_model5int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.juvenile.mmr.er))
    
    ### BIC --------------
    MMR_BIC<-BICdelta(BIC(
                    Phylo_MMR_model2, Phylo_MMR_model2int, #Phylo_MMR_model2.POLY, Phylo_MMR_model2intPOLY,
                    Phylo_MMR_model4, Phylo_MMR_model4int, #Phylo_MMR_model4.POLY, Phylo_MMR_model4intPOLY,
                    Phylo_MMR_model5, Phylo_MMR_model5int))#, Phylo_MMR_model5.POLY, Phylo_MMR_model5intPOLY))
    
    ## AAS optimal ------------------
    #
    Phylo_AS_model2 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.juvenile.aas.er))
    Phylo_AS_model2int <- Almer(lnAS ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.juvenile.aas.er))
    
    Phylo_AS_model4 <- Almer(lnAS ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.juvenile.aas.er))
    Phylo_AS_model4int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.juvenile.aas.er))
    
    Phylo_AS_model5 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.juvenile.aas.er))
    Phylo_AS_model5int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.juvenile.aas.er))
    
    
    ### BIC --------
    AS_BIC<-BICdelta(BIC(
    
                    Phylo_AS_model2, Phylo_AS_model2int,# Phylo_AS_model2.POLY, Phylo_AS_model2intPOLY,
                    Phylo_AS_model4, Phylo_AS_model4int, #Phylo_AS_model4.POLY, Phylo_AS_model4intPOLY,
                    Phylo_AS_model5, Phylo_AS_model5int))#, #Phylo_AS_model5.POLY, Phylo_AS_model5intPOLY))
    
    
    ## FAS optimal ----------------
    # 
    Phylo_FAS_model2 <- Almer(log(FAS) ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.juvenile.fas.er))
    Phylo_FAS_model2int <- Almer(log(FAS) ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.juvenile.fas.er))
    
    Phylo_FAS_model4 <- Almer(log(FAS) ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.juvenile.fas.er))
    Phylo_FAS_model4int <- Almer(log(FAS) ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.juvenile.fas.er))
    
    Phylo_FAS_model5 <- Almer(log(FAS) ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.juvenile.fas.er))
    Phylo_FAS_model5int <- Almer(log(FAS) ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.juvenile.fas.er))
    
    
    ### BIC -------------
    FAS_BIC<-BICdelta(BIC(
                    Phylo_FAS_model2, Phylo_FAS_model2int,# Phylo_FAS_model2.POLY, Phylo_FAS_model2intPOLY,
                    Phylo_FAS_model4, Phylo_FAS_model4int,# Phylo_FAS_model4.POLY, Phylo_FAS_model4intPOLY,
                    Phylo_FAS_model5, Phylo_FAS_model5int))#, Phylo_FAS_model5.POLY, Phylo_FAS_model5intPOLY))
    
    
    ## RMR warm ------------------
    # 
    Phylo_RMR_W_model2 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.juvenile.rmr.w))
    Phylo_RMR_W_model2int <- Almer(lnRMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.juvenile.rmr.w))
    
    Phylo_RMR_W_model4 <- Almer(lnRMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.juvenile.rmr.w))
    Phylo_RMR_W_model4int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.juvenile.rmr.w))
    
    Phylo_RMR_W_model5 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.juvenile.rmr.w))
    Phylo_RMR_W_model5int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.juvenile.rmr.w))
    
    
    
    ## BIC -------------
    RMR_W_BIC<-BICdelta(BIC(
                    Phylo_RMR_W_model2, Phylo_RMR_W_model2int,# Phylo_RMR_W_model2.POLY, Phylo_RMR_W_model2intPOLY,
                    Phylo_RMR_W_model4, Phylo_RMR_W_model4int,# Phylo_RMR_W_model4.POLY, Phylo_RMR_W_model4intPOLY,
                    Phylo_RMR_W_model5, Phylo_RMR_W_model5int))#, Phylo_RMR_W_model5.POLY, Phylo_RMR_W_model5intPOLY))
    
    # AMR / warm temps -------------
    Phylo_MMR_W_model2 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.juvenile.mmr.w))
    Phylo_MMR_W_model2int <- Almer(lnAMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.juvenile.mmr.w))
    
    Phylo_MMR_W_model4 <- Almer(lnAMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.juvenile.mmr.w))
    Phylo_MMR_W_model4int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.juvenile.mmr.w))
    
    Phylo_MMR_W_model5 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.juvenile.mmr.w))
    Phylo_MMR_W_model5int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.juvenile.mmr.w))
    
    
    ## BIC -------------
    MMR_W_BIC<-BICdelta(BIC(
                    Phylo_MMR_W_model2, Phylo_MMR_W_model2int, #Phylo_MMR_W_model2.POLY, Phylo_MMR_W_model2intPOLY,
                    Phylo_MMR_W_model4, Phylo_MMR_W_model4int, #Phylo_MMR_W_model4.POLY, Phylo_MMR_W_model4intPOLY,
                    Phylo_MMR_W_model5, Phylo_MMR_W_model5int))#, Phylo_MMR_W_model5.POLY, Phylo_MMR_W_model5intPOLY))
    
    
    
    # AS / warm temps --------------
    Phylo_AS_W_model2 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.juvenile.aas.w))
    Phylo_AS_W_model2int <- Almer(lnAS ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.juvenile.aas.w))
    
    Phylo_AS_W_model4 <- Almer(lnAS ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.juvenile.aas.w))
    Phylo_AS_W_model4int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.juvenile.aas.w))
    
    Phylo_AS_W_model5 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.juvenile.aas.w))
    Phylo_AS_W_model5int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.juvenile.aas.w))
    
    ## BIC -------------
    AS_W_BIC<-BICdelta(BIC(
                    Phylo_AS_W_model2, Phylo_AS_W_model2int,# Phylo_AS_W_model2.POLY, Phylo_AS_W_model2intPOLY,
                    Phylo_AS_W_model4, Phylo_AS_W_model4int,# Phylo_AS_W_model4.POLY, Phylo_AS_W_model4intPOLY,
                    Phylo_AS_W_model5, Phylo_AS_W_model5int))#, Phylo_AS_W_model5.POLY, Phylo_AS_W_model5intPOLY))
    
    # FAS / warm temps --------------
    Phylo_FAS_W_model2 <- Almer(lnFAS ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.juvenile.fas.w))
    Phylo_FAS_W_model2int <- Almer(lnFAS ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.juvenile.fas.w))
    
    Phylo_FAS_W_model4 <- Almer(lnFAS ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.juvenile.fas.w))
    Phylo_FAS_W_model4int <- Almer(lnFAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.juvenile.fas.w))
    
    Phylo_FAS_W_model5 <- Almer(lnFAS ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.juvenile.fas.w))
    Phylo_FAS_W_model5int <- Almer(lnFAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.juvenile.fas.w))
    
    
    ## BIC -------------
    FAS_W_BIC<-BICdelta(BIC(
                    Phylo_FAS_W_model2, Phylo_FAS_W_model2int, #Phylo_FAS_W_model2.POLY, Phylo_FAS_W_model2intPOLY,
                    Phylo_FAS_W_model4, Phylo_FAS_W_model4int, #Phylo_FAS_W_model4.POLY, Phylo_FAS_W_model4intPOLY,
                    Phylo_FAS_W_model5, Phylo_FAS_W_model5int)) #Phylo_FAS_W_model5.POLY, Phylo_FAS_W_model5intPOLY))
    
    
    ## BIC results ------
    # RMR_BIC#
    # print(any(data.amr$study_ID == 415))
    # MMR_BIC #
    # AS_BIC  # 
    # FAS_BIC # 
    # 
    # RMR_W_BIC # 
    # MMR_W_BIC # 
    # AS_W_BIC  # 
    # FAS_W_BIC #
    check_best_model(Phylo_MMR_model4int, MMR_BIC)
    check_best_model(Phylo_AS_model2int, AS_BIC)
    check_best_model(Phylo_FAS_model2int, FAS_BIC)
    
    check_best_model(Phylo_RMR_W_model2, RMR_W_BIC)
    check_best_model(Phylo_MMR_W_model4, MMR_W_BIC)
    check_best_model(Phylo_AS_W_model2, AS_W_BIC)
    check_best_model(Phylo_FAS_W_model2, FAS_W_BIC)
    # jan 2026 best
    # ecol relevant
    if(any(data.amrER$study_ID == 415)){
      rmr_mod_ER<- Phylo_RMR_model5 # this is when including wootton data!! 
      check_best_model(Phylo_RMR_model5, RMR_BIC)
    }else{
      rmr_mod_ER<- Phylo_RMR_model4 # this is when NOT including wootton data!! 
      check_best_model(Phylo_RMR_model4, RMR_BIC)

    }
    amr_mod_ER<-Phylo_MMR_model4int
    as_mod_ER<-Phylo_AS_model2int
    fas_mod_ER<-Phylo_FAS_model2int
  
    # warm
    rmr_mod_W<-Phylo_RMR_W_model2
    amr_mod_W<-Phylo_MMR_W_model4
    as_mod_W<-Phylo_AS_W_model2
    fas_mod_W<-Phylo_FAS_W_model2
    
    write.csv(file = here("Data_exports/BICs/RMR_BIC_juvenile.csv"), x = RMR_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/MMR_BIC_juvenile.csv"), x = MMR_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/AS_BIC_juvenile.csv"), x = AS_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/FAS_BIC_juvenile.csv"), x = FAS_BIC, row.names = T)
    
    write.csv(file = here("Data_exports/BICs/RMR_W_BIC_juvenile.csv"), x = RMR_W_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/MMR_W_BIC_juvenile.csv"), x = MMR_W_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/AS_W_BIC_juvenile.csv"), x = AS_W_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/FAS_W_BIC_juvenile.csv"), x = FAS_W_BIC, row.names = T)
    
    } 
    
  } else { # full dataset
    
    # rmr
    get_phylo_matrix(species.list = unique(levels(data.rmrER$species)),
                     matrix.name = "A",
                     dataset.ID = "RMR-optimal")
    get_phylo_matrix(species.list = unique(levels(data.rmr.test$species)),
                     matrix.name = "A.rmr.w",
                     dataset.ID = "RMR-warm")
    # amr 
    get_phylo_matrix(species.list = unique(levels(data.amrER$species)),
                     matrix.name = "A.mmr.er",
                     dataset.ID = "MMR-optimal")
    get_phylo_matrix(species.list = unique(levels(data.amr.test$species)),
                     matrix.name = "A.mmr.w",
                     dataset.ID = "MMR-warm")
    # aas
    get_phylo_matrix(species.list = unique(levels(data.asER$species)),
                     matrix.name = "A.aas.er",
                     dataset.ID = "AAS-optimal")
    get_phylo_matrix(species.list = unique(levels(data.as.test$species)),
                     matrix.name = "A.aas.w",
                     dataset.ID = "AAS-warm")
    # fas
    get_phylo_matrix(species.list = unique(levels(data.fasER$species)),
                     matrix.name = "A.fas.er",
                     dataset.ID = "FAS-optimal")
    get_phylo_matrix(species.list = unique(levels(data.fas.test$species)),
                     matrix.name = "A.fas.w",
                     dataset.ID = "FAS-warm")
  
  
    ## RMR optimal -----------------
    Phylo_RMR_model2 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
    Phylo_RMR_model2int <- Almer(lnRMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))

    Phylo_RMR_model4 <- Almer(lnRMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
    Phylo_RMR_model4int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))

    Phylo_RMR_model5 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
    Phylo_RMR_model5int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))

  
    # No singular fits 
    ### BIC rmr optimal --------------
    RMR_BIC<-BICdelta(BIC(
  
                    Phylo_RMR_model2, Phylo_RMR_model2int, #Phylo_RMR_model2.POLY, Phylo_RMR_model2intPOLY,
                    Phylo_RMR_model4, Phylo_RMR_model4int,# Phylo_RMR_model4.POLY, Phylo_RMR_model4intPOLY,
                    Phylo_RMR_model5, Phylo_RMR_model5int# Phylo_RMR_model5.POLY, Phylo_RMR_model5intPOLY
  ))
    
    ## MMR optimal -----------------
    Phylo_MMR_model2 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
    Phylo_MMR_model2int <- Almer(lnAMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))

    Phylo_MMR_model4 <- Almer(lnAMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
    Phylo_MMR_model4int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  
    Phylo_MMR_model5 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
    Phylo_MMR_model5int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))

    ### BIC --------------
    MMR_BIC<-BICdelta(BIC(
                    Phylo_MMR_model2, Phylo_MMR_model2int, #Phylo_MMR_model2.POLY, Phylo_MMR_model2intPOLY,
                    Phylo_MMR_model4, Phylo_MMR_model4int, #Phylo_MMR_model4.POLY, Phylo_MMR_model4intPOLY,
                    Phylo_MMR_model5, Phylo_MMR_model5int))#, Phylo_MMR_model5.POLY, Phylo_MMR_model5intPOLY))
    
    ## AAS optimal ------------------
    #
    Phylo_AS_model2 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
    Phylo_AS_model2int <- Almer(lnAS ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))

    Phylo_AS_model4 <- Almer(lnAS ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
    Phylo_AS_model4int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
 
    Phylo_AS_model5 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
    Phylo_AS_model5int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))

  
    ### BIC --------
    AS_BIC<-BICdelta(BIC(
  
                    Phylo_AS_model2, Phylo_AS_model2int,# Phylo_AS_model2.POLY, Phylo_AS_model2intPOLY,
                    Phylo_AS_model4, Phylo_AS_model4int, #Phylo_AS_model4.POLY, Phylo_AS_model4intPOLY,
                    Phylo_AS_model5, Phylo_AS_model5int))#, #Phylo_AS_model5.POLY, Phylo_AS_model5intPOLY))
    
    
    ## FAS optimal ----------------
    # 
    Phylo_FAS_model2 <- Almer(log(FAS) ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
    Phylo_FAS_model2int <- Almer(log(FAS) ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
 
    Phylo_FAS_model4 <- Almer(log(FAS) ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
    Phylo_FAS_model4int <- Almer(log(FAS) ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))

    Phylo_FAS_model5 <- Almer(log(FAS) ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
    Phylo_FAS_model5int <- Almer(log(FAS) ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))

  
    ### BIC -------------
    FAS_BIC<-BICdelta(BIC(
                    Phylo_FAS_model2, Phylo_FAS_model2int,# Phylo_FAS_model2.POLY, Phylo_FAS_model2intPOLY,
                    Phylo_FAS_model4, Phylo_FAS_model4int,# Phylo_FAS_model4.POLY, Phylo_FAS_model4intPOLY,
                    Phylo_FAS_model5, Phylo_FAS_model5int))#, Phylo_FAS_model5.POLY, Phylo_FAS_model5intPOLY))
    
    
    ## RMR warm ------------------
    # 
    Phylo_RMR_W_model2 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
    Phylo_RMR_W_model2int <- Almer(lnRMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))

    Phylo_RMR_W_model4 <- Almer(lnRMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
    Phylo_RMR_W_model4int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))

    Phylo_RMR_W_model5 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
    Phylo_RMR_W_model5int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))

  
  
    ## BIC -------------
    RMR_W_BIC<-BICdelta(BIC(
                    Phylo_RMR_W_model2, Phylo_RMR_W_model2int,# Phylo_RMR_W_model2.POLY, Phylo_RMR_W_model2intPOLY,
                    Phylo_RMR_W_model4, Phylo_RMR_W_model4int,# Phylo_RMR_W_model4.POLY, Phylo_RMR_W_model4intPOLY,
                    Phylo_RMR_W_model5, Phylo_RMR_W_model5int))#, Phylo_RMR_W_model5.POLY, Phylo_RMR_W_model5intPOLY))
    
    # AMR / warm temps -------------
    Phylo_MMR_W_model2 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
    Phylo_MMR_W_model2int <- Almer(lnAMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))

    Phylo_MMR_W_model4 <- Almer(lnAMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
    Phylo_MMR_W_model4int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))

    Phylo_MMR_W_model5 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
    Phylo_MMR_W_model5int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))

   
    ## BIC -------------
    MMR_W_BIC<-BICdelta(BIC(
                    Phylo_MMR_W_model2, Phylo_MMR_W_model2int, #Phylo_MMR_W_model2.POLY, Phylo_MMR_W_model2intPOLY,
                    Phylo_MMR_W_model4, Phylo_MMR_W_model4int, #Phylo_MMR_W_model4.POLY, Phylo_MMR_W_model4intPOLY,
                    Phylo_MMR_W_model5, Phylo_MMR_W_model5int))#, Phylo_MMR_W_model5.POLY, Phylo_MMR_W_model5intPOLY))
    
    
    
    # AS / warm temps --------------
    Phylo_AS_W_model2 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
    Phylo_AS_W_model2int <- Almer(lnAS ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
 
    Phylo_AS_W_model4 <- Almer(lnAS ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
    Phylo_AS_W_model4int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))

    Phylo_AS_W_model5 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
    Phylo_AS_W_model5int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))

    ## BIC -------------
    AS_W_BIC<-BICdelta(BIC(
                    Phylo_AS_W_model2, Phylo_AS_W_model2int,# Phylo_AS_W_model2.POLY, Phylo_AS_W_model2intPOLY,
                    Phylo_AS_W_model4, Phylo_AS_W_model4int,# Phylo_AS_W_model4.POLY, Phylo_AS_W_model4intPOLY,
                    Phylo_AS_W_model5, Phylo_AS_W_model5int))#, Phylo_AS_W_model5.POLY, Phylo_AS_W_model5intPOLY))
    
    # FAS / warm temps --------------
    Phylo_FAS_W_model2 <- Almer(lnFAS ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
    Phylo_FAS_W_model2int <- Almer(lnFAS ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))

    Phylo_FAS_W_model4 <- Almer(lnFAS ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
    Phylo_FAS_W_model4int <- Almer(lnFAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))

    Phylo_FAS_W_model5 <- Almer(lnFAS ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
    Phylo_FAS_W_model5int <- Almer(lnFAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))

    
    ## BIC -------------
    FAS_W_BIC<-BICdelta(BIC(
                    Phylo_FAS_W_model2, Phylo_FAS_W_model2int, #Phylo_FAS_W_model2.POLY, Phylo_FAS_W_model2intPOLY,
                    Phylo_FAS_W_model4, Phylo_FAS_W_model4int, #Phylo_FAS_W_model4.POLY, Phylo_FAS_W_model4intPOLY,
                    Phylo_FAS_W_model5, Phylo_FAS_W_model5int)) #Phylo_FAS_W_model5.POLY, Phylo_FAS_W_model5intPOLY))
    
    
    ## BIC results ------
    # RMR_BIC #
    # MMR_BIC #
    # AS_BIC #
    # FAS_BIC #
    # 
    # RMR_W_BIC #
    # MMR_W_BIC #
    # AS_W_BIC  #
    # FAS_W_BIC #
    #
    # model checks: 
    check_best_model(Phylo_RMR_model4int, RMR_BIC)
    check_best_model(Phylo_MMR_model4, MMR_BIC)
    check_best_model(Phylo_AS_model5int, AS_BIC)
    check_best_model(Phylo_FAS_model2int, FAS_BIC)
    
    check_best_model(Phylo_RMR_W_model2, RMR_W_BIC)
    check_best_model(Phylo_MMR_W_model2int, MMR_W_BIC)
    check_best_model(Phylo_AS_W_model2int, AS_W_BIC)
    check_best_model(Phylo_FAS_W_model2, FAS_W_BIC)
    
    # jan 2026 best 
    # ecol relevant
    rmr_mod_ER<- Phylo_RMR_model4int
    amr_mod_ER<-Phylo_MMR_model4
    as_mod_ER<- Phylo_AS_model5int
    fas_mod_ER<- Phylo_FAS_model2int
  
    # warm
    rmr_mod_W<-Phylo_RMR_W_model2
    amr_mod_W<-Phylo_MMR_W_model2int
    as_mod_W<- Phylo_AS_W_model2int
    fas_mod_W<-Phylo_FAS_W_model2
    

    write.csv(file = here("Data_exports/BICs/RMR_BIC.csv"), x = RMR_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/MMR_BIC.csv"), x = MMR_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/AS_BIC.csv"), x = AS_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/FAS_BIC.csv"), x = FAS_BIC, row.names = T)
  
    write.csv(file = here("Data_exports/BICs/RMR_W_BIC.csv"), x = RMR_W_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/MMR_W_BIC.csv"), x = MMR_W_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/AS_W_BIC.csv"), x = AS_W_BIC, row.names = T)
    write.csv(file = here("Data_exports/BICs/FAS_W_BIC.csv"), x = FAS_W_BIC, row.names = T)
    
  }
    

  # *****************************************************************
  # *****************************************************************
  # # # Ecol relev
  # plot(rmr_mod_ER) # residuals
  # plot(amr_mod_ER) # residuals
  # plot(fas_mod_ER) # residuals
  # plot(as_mod_ER) # residuals
  # # # warm
  # plot(rmr_mod_W) # residuals
  # plot(amr_mod_W) # residuals
  # plot(fas_mod_W) # residuals
  # plot(as_mod_W) # residuals
  # # 
  # hist(resid(amr_mod_ER), breaks = 50) # residuals
  # hist(resid(rmr_mod_ER), breaks = 50) # residuals
  # hist(resid(fas_mod_ER), breaks = 50) # residuals
  # hist(resid(as_mod_ER), breaks = 50) # residuals
  # # 
  # hist(resid(amr_mod_W), breaks = 50) # residuals
  # hist(resid(rmr_mod_W), breaks = 50) # residuals
  # hist(resid(fas_mod_W), breaks = 50) # residuals
  # hist(resid(as_mod_W), breaks = 50) # little skew not too bad, only 5 measurements, all reasonable biologically
  
  # warm models
  run_anova(rmr_mod_ER, save.path = here(paste("Data_exports/", lifestage, "models/RMR_ER_anova_summary_report.txt", sep = "")))
  run_anova(amr_mod_ER, save.path = here(paste("Data_exports/", lifestage, "models/MMR_ER_anova_summary_report.txt", sep = "")))
  run_anova(as_mod_ER, save.path = here(paste("Data_exports/", lifestage, "models/AS_ER_anova_summary_report.txt", sep = "")))
  run_anova(fas_mod_ER, save.path = here(paste("Data_exports/", lifestage, "models/FAS_ER_anova_summary_report.txt", sep = "")))
  
  # ER models
  run_anova(rmr_mod_W, save.path = here(paste("Data_exports/", lifestage, "models/RMR_W_anova_summary_report.txt", sep = "")))
  run_anova(amr_mod_W, save.path = here(paste("Data_exports/", lifestage, "models/MMR_W_anova_summary_report.txt", sep = "")))
  run_anova(as_mod_W, save.path = here(paste("Data_exports/", lifestage, "models/AS_W_anova_summary_report.txt", sep = "")))
  run_anova(fas_mod_W, save.path = here(paste("Data_exports/", lifestage, "models/FAS_W_anova_summary_report.txt", sep = "")))
  
  # *****************************************************************************
  # *****************************************************************************
  # Model scaling parameters and CIs ----------
  # custom function to obtain model parameters and recalculated mass specific MR using model estimate scaling slopes
  # this also expands the dataset to get mass independent values of metabolic rates
  # function source script available at 'mixed_model_outputs.R'
  # tic() # time the function run
  model_out<-get_model_outputs(
                best.model.rmr.er = rmr_mod_ER,
                best.model.amr.er= amr_mod_ER,
                best.model.as.er = as_mod_ER,
                best.model.fas.er = fas_mod_ER,
                best.model.rmr.w = rmr_mod_W,
                best.model.amr.w = amr_mod_W,
                best.model.as.w = as_mod_W,
                best.model.fas.w = fas_mod_W,
                data.rmr.test.modpar = data.rmr.test,
                data.rmrER.modpar = data.rmrER,
                data.amr.test.modpar = data.amr.test,
                data.amrER.modpar = data.amrER,                        
                data.as.test.modpar = data.as.test,
                data.asER.modpar = data.asER, 
                data.fas.test.modpar = data.fas.test,
                data.fasER.modpar = data.fasER, 
                name_extension = paste(lifestage, "models", sep="")) #
  # toc()

  sum_CItable<-data.frame(model_out[[1]])
  colnames(sum_CItable)<-c("ci5", "ci95","var_repeat", "MR", "temp_cat")
  sum_data<-data.frame(model_out[[2]])
  scaling_params<-data.frame(model_out[[3]])

  # *****************************************************************************
  # *****************************************************************************

  # Figures -------
  set.seed(51423)
  
  # General scaling plots ------
  scaling_params_amr_w<-scaling_params %>% 
    filter(performance == "MMR",
           temp_categ =="warm",
           Temperature == 20) # one slope
  scaling_params_amr_er<-scaling_params %>% 
    filter(performance == "MMR",
           temp_categ =="ecol_relev",
           Temperature == 20) # change with temp
  
  AMRmodel_plot1<- ggplot(data=data.amrER, aes(x=lnBWg, y=lnAMR)) +
    geom_point(alpha=0.9,  size=1, pch=1, color="grey75")+
    geom_point(data=data.amr.test, aes(x=lnBWg, y=lnAMR),
               alpha=0.9,  size=1, pch=21, show.legend = FALSE, stroke =0.2,
               fill = cols.amr[3], color = cols.amr[3])+
    geom_segment(aes(x = min(data.amr.test$lnBWg),
                     xend = max(data.amr.test$lnBWg),
                   y = scaling_params_amr_w$Intercept +
                     scaling_params_amr_w$Slope *
                     min(data.amr.test$lnBWg),
                   yend = scaling_params_amr_w$Intercept +
                     scaling_params_amr_w$Slope *
                     max(data.amr.test$lnBWg)),
                color = cols.amr[1])+
    geom_segment(aes(x = min(data.amrER$lnBWg),
                     xend = max(data.amrER$lnBWg),
                   y = scaling_params_amr_er$Intercept +
                     scaling_params_amr_er$Slope *
                     min(data.amrER$lnBWg),
                   yend = scaling_params_amr_er$Intercept +
                     scaling_params_amr_er$Slope *
                     max(data.amrER$lnBWg)),
                color = "black")+
    annotate("text",  x = -5.2, y = 11.5,
             label = bquote(Optimal:~italic(b)[MMR] == .(round(scaling_params_amr_er$Slope,3))), 
                            # "\u00b1" ~ .(round(scaling_params_amr_er$SE_slope,2))),
             size=4, hjust=0, family="Helvetica", color = "black")+
    annotate("text",  x = -5.2, y = 9.8,
             label = bquote(Warm:~italic(b)[MMR] == .(round(scaling_params_amr_w$Slope,3))),
             size=4, hjust=0, family="Helvetica", color = cols.amr[1])+
    scale_color_gradient( low = "grey", high = "black")+
    ylim(x = -6.3, 12)+
    xlim(x = -6.3, 12)+
    annotate("text", label = paste("n = ", nrow(data.amrER), sep=""),
             x = -5.2, y = 8.5, size=3, hjust=0, family="Helvetica", color = "black")+
    annotate("text", label = paste("n = ", nrow(data.amr.test), sep=""),
             x = -5.2, y = 7.6, size=3, hjust=0, family="Helvetica", color = cols.amr[1])+
    annotate("text", label = lifestage,
             x = 5, y = -5, size=5, hjust=0, family="Helvetica", color = "black")
  
    if(is.null(lifestage)){
      AMRmodel_plot1<-AMRmodel_plot1+ annotate("text",  x = 7.5, y = 9.8,
             label = bquote("\u2193" ~ "with" ~ degree*C),
             size=4, hjust=0,  color = cols.amr[1])
    }else{
      if(lifestage == "juvenile"){
        AMRmodel_plot1<-AMRmodel_plot1+ annotate("text",  x = 7.5, y = 11.5,
           label = bquote("\u2193" ~ "with" ~ degree*C),
           size=4, hjust=0,  color = "black")
      }
      if(lifestage == "adult"){
        AMRmodel_plot1<-AMRmodel_plot1+ annotate("text",  x = 7.5, y = 11.5,
           label = bquote("\u2191" ~ "with" ~ degree*C),
           size=4, hjust=0,  color = "black")
      }
    }
  
  ggformat(AMRmodel_plot1, x_title=expression(italic(ln)*Body~weight~(g)),
           y_title=expression(italic(ln)*MMR~(mg~O[2]~h^-1)), size_text = 12, print = T)
  
  
  # RMR 
  scaling_params_rmr_w<-scaling_params %>%
    filter(performance == "RMR",
           temp_categ =="warm",
           Temperature == 20) # one slope
  scaling_params_rmr_er<-scaling_params %>% 
    filter(performance == "RMR",
           temp_categ =="ecol_relev",
           Temperature == 20) # one slope
  
  RMRmodel_plot1<- ggplot(data=data.rmrER, aes(x=lnBWg, y=lnRMR)) +
    geom_point(alpha=0.9,  size=1, pch=1, color="grey75")+
    geom_point(data=data.rmr.test, aes(x=lnBWg, y=lnRMR),
               alpha=0.9,  size=1, pch=21, show.legend = FALSE, stroke =0.2,
               fill = cols.rmr[2], color = cols.rmr[2])+
    geom_segment(aes(x = min(data.rmr.test$lnBWg),
                     xend = max(data.rmr.test$lnBWg),
                   y = scaling_params_rmr_w$Intercept +
                     scaling_params_rmr_w$Slope *
                     min(data.rmr.test$lnBWg),
                   yend = scaling_params_rmr_w$Intercept +
                     scaling_params_rmr_w$Slope *
                     max(data.rmr.test$lnBWg)),
                color = cols.rmr[1])+
    geom_segment(aes(x = min(data.rmrER$lnBWg),
                     xend = max(data.rmrER$lnBWg),
                   y = scaling_params_rmr_er$Intercept +
                     scaling_params_rmr_er$Slope *
                     min(data.rmrER$lnBWg),
                   yend = scaling_params_rmr_er$Intercept +
                     scaling_params_rmr_er$Slope *
                     max(data.rmrER$lnBWg)),
                color = "black")+
    annotate("text",  x = -5.2, y = 11.5,
             label = bquote(Optimal:~italic(b)[RMR] == .(round(scaling_params_rmr_er$Slope,3))),
             size=4, hjust=0, family="Helvetica", color = "black")+
    annotate("text",  x = -5.2, y = 9.8,
             label = bquote(Warm:~italic(b)[RMR] == .(round(scaling_params_rmr_w$Slope,3))),
             size=4, hjust=0, family="Helvetica", color = cols.rmr[1])+
    scale_color_gradient( low = "grey", high = "black")+
    ylim(x = -6.3, 12)+
    xlim(x = -6.3, 12)+
    annotate("text", label = paste("n = ", nrow(data.rmrER), sep=""),
             x = -5.2, y = 8.5, size=3, hjust=0, family="Helvetica", color = "black")+
    annotate("text", label = paste("n = ", nrow(data.rmr.test), sep=""),
             x = -5.2, y = 7.6, size=3, hjust=0, family="Helvetica", color = cols.rmr[1])
  
    if(is.null(lifestage)){
      RMRmodel_plot1<-RMRmodel_plot1+ annotate("text",  x = 7.5, y = 11.5,
             label = bquote("\u2191" ~ "with" ~ degree*C),
             size=4, hjust=0,  color = "black")
    }else{

    }
    # annotate("text", label = "ADULT FISH",
             # x = 5, y = -5, size=5, hjust=0, family="Helvetica", color = "black")
  ggformat(RMRmodel_plot1, x_title=expression(italic(ln)*Body~weight~(g)),
           y_title=expression(italic(ln)*RMR~(mg~O[2]~h^-1)), size_text = 12, print = T)
  
  
  # AS 
  scaling_params_as_w<-scaling_params %>%
    filter(performance == "AS",
           temp_categ =="warm",
           Temperature == 20) # change with temp POLY
  scaling_params_as_er<-scaling_params %>% 
    filter(performance == "AS",
           temp_categ =="ecol_relev",
           Temperature == 20) # 
  
  ASmodel_plot1<- ggplot(data=data.asER, aes(x=lnBWg, y=lnAS)) +
    geom_point(alpha=0.9,  size=1, pch=1, color="grey75")+
    geom_point(data=data.as.test, aes(x=lnBWg, y=lnAS),
               alpha=0.9,  size=1, pch=21, show.legend = FALSE, stroke =0.2,
               fill = cols.as[3], color = cols.as[3])+
    geom_segment(aes(x = min(data.as.test$lnBWg),
                     xend = max(data.as.test$lnBWg),
                   y = scaling_params_as_w$Intercept +
                     scaling_params_as_w$Slope *
                     min(data.as.test$lnBWg),
                   yend = scaling_params_as_w$Intercept +
                     scaling_params_as_w$Slope *
                     max(data.as.test$lnBWg)),
                color = cols.as[1])+
    geom_segment(aes(x = min(data.asER$lnBWg),
                     xend = max(data.asER$lnBWg),
                   y = scaling_params_as_er$Intercept +
                     scaling_params_as_er$Slope *
                     min(data.asER$lnBWg),
                   yend = scaling_params_as_er$Intercept +
                     scaling_params_as_er$Slope *
                     max(data.asER$lnBWg)),
                color = "black")+
    annotate("text",  x = -5.2, y = 11.5,
             label = bquote(Optimal:~italic(b)[AS] == .(round(scaling_params_as_er$Slope,3))),
             size=4, hjust=0, family="Helvetica", color = "black")+
    annotate("text",  x = -5.2, y = 9.8,
             label = bquote(Warm:~italic(b)[AS] == .(round(scaling_params_as_w$Slope,3))),
             size=4, hjust=0, family="Helvetica", color = cols.as[1])+
    scale_color_gradient( low = "grey", high = "black")+
    ylim(x = -6.3, 12)+
    xlim(x = -6.3, 12)+
    annotate("text", label = paste("n = ", nrow(data.asER), sep=""),
             x = -5.2, y = 8.5, size=3, hjust=0, family="Helvetica", color = "black")+
    annotate("text", label = paste("n = ", nrow(data.as.test), sep=""),
             x = -5.2, y = 7.6, size=3, hjust=0, family="Helvetica", color = cols.as[1])
    if(is.null(lifestage)){
      ASmodel_plot1<-ASmodel_plot1+ annotate("text",  x = 7.5, y = 11.5,
             label = bquote("\u2193" ~ "with" ~ degree*C),
             size=4, hjust=0,  color = "black")
      ASmodel_plot1<-ASmodel_plot1+ annotate("text",  x = 7, y = 9.8,
             label = bquote("\u2193" ~ "with" ~ degree*C),
             size=4, hjust=0,  color = cols.as[1])
    }else{
      if(lifestage == "juvenile"){
        ASmodel_plot1<-ASmodel_plot1+ annotate("text",  x = 7, y = 11.5,
             label = bquote("\u2193" ~ "with" ~ degree*C),
             size=4, hjust=0,  color = "black")
      }
      if(lifestage == "adult"){
        ASmodel_plot1<-ASmodel_plot1+ annotate("text",  x = 7, y = 9.8,
             label = bquote("\u2191" ~ "with" ~ degree*C),
             size=4, hjust=0,  color = cols.as[1])
      }

    }
    # annotate("text", label = "ADULT FISH",
    #          x = 5, y = -5, size=5, hjust=0, family="Helvetica", color = "black")
  ggformat(ASmodel_plot1, x_title=expression(italic(ln)*Body~weight~(g)),
           y_title=expression(italic(ln)*AS~(mg~O[2]~h^-1)), size_text = 12, print = T)
  
  # FAS
  scaling_params_fas_w<-scaling_params %>%
    filter(performance == "FAS",
           temp_categ =="warm",
           Temperature == 20) # change with temp POLY
  scaling_params_fas_er<-scaling_params %>% 
    filter(performance == "FAS",
           temp_categ =="ecol_relev",
           Temperature == 20) # 
  
  FASmodel_plot1<- ggplot(data=data.fasER, aes(x=lnBWg, y=lnFAS)) +
    geom_point(alpha=0.9,  size=1, pch=1, color="grey75")+
    geom_point(data=data.fas.test, aes(x=lnBWg, y=lnFAS),
               alpha=0.9,  size=1, pch=21, show.legend = FALSE, stroke =0.2,
               fill = cols.fas[3], color = cols.fas[3])+
    geom_segment(aes(x = min(data.fas.test$lnBWg),
                     xend = max(data.fas.test$lnBWg),
                   y = scaling_params_fas_w$Intercept +
                     scaling_params_fas_w$Slope *
                     min(data.fas.test$lnBWg),
                   yend = scaling_params_fas_w$Intercept +
                     scaling_params_fas_w$Slope *
                     max(data.fas.test$lnBWg)),
                color = cols.fas[1])+
    geom_segment(aes(x = min(data.fasER$lnBWg),
                     xend = max(data.fasER$lnBWg),
                   y = scaling_params_fas_er$Intercept +
                     scaling_params_fas_er$Slope *
                     min(data.fasER$lnBWg),
                   yend = scaling_params_fas_er$Intercept +
                     scaling_params_fas_er$Slope *
                     max(data.fasER$lnBWg)),
                color = "black")+
    annotate("text",  x = -5.2, y = 11.5,
             label = bquote(Optimal:~italic(b)[FAS] == .(round(scaling_params_fas_er$Slope,3))),
             size=4, hjust=0, family="Helvetica", color = "black")+
    annotate("text",  x = -5.2, y = 9.8,
             label = bquote(Warm:~italic(b)[FAS] == .(round(scaling_params_fas_w$Slope,3))),
             size=4, hjust=0, family="Helvetica", color = cols.fas[1])+
    scale_color_gradient( low = "grey", high = "black")+
    ylim(x = -6.3, 12)+
    xlim(x = -6.3, 12)+
    annotate("text", label = paste("n = ", nrow(data.fasER), sep=""),
             x = -5.2, y = 8.5, size=3, hjust=0, family="Helvetica", color = "black")+
    annotate("text", label = paste("n = ", nrow(data.fas.test), sep=""),
             x = -5.2, y = 7.6, size=3, hjust=0, family="Helvetica", color = cols.fas[1])
    if(is.null(lifestage)){
      FASmodel_plot1<-FASmodel_plot1+ annotate("text",  x = 7.5, y = 11.5,
             label = bquote("\u2193" ~ "with" ~ degree*C),
             size=4, hjust=0,  color = "black")
    }else{
      if(lifestage=="juvenile"){
        FASmodel_plot1<-FASmodel_plot1+ annotate("text",  x = 7.5, y = 11.5,
             label = bquote("\u2193" ~ "with" ~ degree*C),
             size=4, hjust=0,  color = "black")
      }
    }
  
  ggformat(FASmodel_plot1, x_title=expression(italic(ln)*Body~weight~(g)),
           y_title=expression(italic(ln)*FAS~(mg~O[2]~h^-1)), size_text = 12, print = T)
  
  
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
  ggsave(filename = paste("./Figures/Figure3", lifestage,".png", sep=""),
         plot=scaling, width = 6.8, height = 6.8, units = "in")
  
  
  # Global outcomes by scaled MR ~ temprature  --------
  AMRmodel_plot1_t<-ggplot()+
    geom_point(data=data.amrER,
               aes(x=tempTest, y=AMR/(BW_g^scaling_params_amr_er$Slope), size = BW_g),
               alpha=0.9, pch=1, size = 1,color="grey75")+
    geom_point(data=data.amr.test,
               mapping = aes(x=tempTest, y=AMR/(BW_g^scaling_params_amr_w$Slope), size = BW_g),
               alpha=0.9, pch=21, size = 1,show.legend = FALSE, stroke =0.2,
               fill = cols.amr[3], color = cols.amr[1])+
    ylim(0,8)+
    xlim(0,40)+
    annotate("text", label = lifestage,
             x = 0, y = 6.5, size=5, hjust=0, family="Helvetica", color = "black")
  ggformat(AMRmodel_plot1_t, x_title="Temperature ºC",
           y_title=bquote(MMR~(mg~O[2]~h^-1~g^italic(b))),
           size_text = 12, print = T)
  
  
  RMRmodel_plot1_t<-ggplot()+
    geom_point(data=data.rmrER,
               aes(x=tempTest, y=RMR/(BW_g^scaling_params_rmr_er$Slope), size = BW_g),
               alpha=0.9, pch=1,size = 1, color="grey75")+
    geom_point(data=data.rmr.test,
               mapping = aes(x=tempTest, y=RMR/(BW_g^scaling_params_rmr_w$Slope), size = BW_g),
               alpha=0.9, pch=21,size = 1, show.legend = FALSE, stroke =0.2,
               fill = cols.rmr[3], color = cols.rmr[1])+
    ylim(0,1.5)+
    xlim(-1,40)
  ggformat(RMRmodel_plot1_t, x_title="Temperature ºC",
           y_title=bquote(RMR~(mg~O[2]~h^-1~g^italic(b))),
           size_text = 12, print = T)
  
  
  ASmodel_plot1_t<-ggplot()+
    geom_point(data=data.asER,
               aes(x=tempTest, y=AS/(BW_g^scaling_params_as_er$Slope), size = BW_g),
               alpha=0.9, pch=1, size = 1,color="grey75")+
    geom_point(data=data.as.test,
               mapping = aes(x=tempTest, y=AS/(BW_g^scaling_params_as_w$Slope), size = BW_g),
               alpha=0.9, pch=21, size = 1,show.legend = FALSE, stroke =0.2,
               fill = cols.as[3], color = cols.as[1])+
    ylim(0,3)+
    xlim(0,35)
  ggformat(ASmodel_plot1_t, x_title="Temperature ºC",
           y_title=bquote(AS~(mg~O[2]~h^-1~g^italic(b))),
           size_text = 12, print = T)
  
  FASmodel_plot1_t<-ggplot()+
    geom_point(data=data.fasER,
               aes(x=tempTest, y=FAS/(BW_g^scaling_params_fas_er$Slope), size = BW_g),
               alpha=0.9, size = 1, pch=1, color="grey75")+
    geom_point(data=data.fas.test,
               mapping = aes(x=tempTest, y=FAS/(BW_g^scaling_params_fas_w$Slope), size = BW_g),
               alpha=0.9, pch=21,size = 1, show.legend = FALSE, stroke =0.2,
               fill = cols.fas[3], color = cols.fas[1])+
    ylim(0,30)+
    xlim(0,35)
  ggformat(FASmodel_plot1_t, x_title="Temperature ºC",
           y_title=bquote(FAS~(mg~O[2]~h^-1~g^italic(b))),
           size_text = 12, print = T)
  
  scaling_t<-cowplot:::plot_grid(AMRmodel_plot1_t, RMRmodel_plot1_t,
                               ASmodel_plot1_t, FASmodel_plot1_t,
                                align = "hv",
                                axis = "l",
                                nrow = 2,
                                ncol = 2,
                                labels = "AUTO",
                               label_x = c(0.2, 0.2),
                               label_y = c(0.895, 0.895),
                                label_size = 12)
  # scaling
  ggsave(filename = paste("./Figures/Figure3-temp", lifestage,".png", sep=""),
         plot=scaling_t, width = 6.8, height = 6.8, units = "in")
  
  
  # **************************************************************
  # **************************************************************
  # ## residuals correlate with any methodological metrics? --------
  # type of MMR protocol 
  if(is.null(lifestage)){
    data.amrER$resid<-resid(amr_mod_ER)
    data.rmrER$resid<-resid(rmr_mod_ER)
    data.amr.test$resid<-resid(amr_mod_W)
    data.rmr.test$resid<-resid(rmr_mod_W)
    
    # MMR
    MMR_method_er<-ggplot(data = data.amrER ,aes( x = MMR_method, y = resid))+
      geom_point(size = 1, position = position_jitter(width = 0.2), 
                 alpha = 0.3, color  = "grey")+
      geom_boxplot(alpha = 0)
    ggformat(MMR_method_er, x_title="MMR method",
             y_title=bquote(residual~MMR),
             size_text = 12, print = F, title = "Ecol Relev")
    
    # warm temps
    MMR_method_w<-ggplot(data = data.amr.test ,aes( x = MMR_method, y = resid))+
      geom_point(size = 1, position = position_jitter(width = 0.2), 
                 alpha = 0.3, color  = cols.amr[3])+
      geom_boxplot(alpha = 0)
    ggformat(MMR_method_w, x_title="MMR method",
             y_title=bquote(residual~MMR),
             size_text = 12, print = F, title = "Warm")
    
    # hours starved
    # RMR 
    RMR_method_starved<-ggplot(data = data.rmrER, aes(x = h_starved, y = resid))+
      geom_smooth(color = "black", method = "loess", se = F)+
      geom_smooth(data = data.rmr.test, method = "loess", se = F, 
                  mapping = aes(x = h_starved, y = resid),
                  color  = cols.rmr[3], fill  = cols.rmr[3])+
      geom_point(size = 1, alpha = 1, color  = "grey")+
      geom_point(data = data.rmr.test, mapping = aes(x = h_starved, y = resid),
                 size = 1, alpha = 1, color  = cols.rmr[3])+
      geom_hline(yintercept = 0)

    ggformat(RMR_method_starved, x_title="Hours starved prior test",
             y_title=bquote(residual~RMR),
             size_text = 12, print = F, title = "All fish together")
    
    # MMR 
    MMR_method_starved<-ggplot(data = data.amrER, aes(x = h_starved, y = resid))+
      geom_smooth(color = "black", method = "loess", se = F)+
      geom_smooth(data = data.amr.test, method = "loess", se = F, 
                  mapping = aes(x = h_starved, y = resid),
                  color  = cols.amr[3], fill  = cols.amr[3])+
      geom_point(size = 1, alpha = 0.3, color  = "grey")+
      geom_point(data = data.amr.test, mapping = aes(x = h_starved, y = resid),
                 size = 1, alpha = 1, color  = cols.amr[3])+
      geom_hline(yintercept = 0)
    ggformat(MMR_method_starved, x_title="Hours starved prior test",
             y_title=bquote(residual~MMR),
             size_text = 12, print = F, title = "All fish together")
    
    # hours in respo 
    # RMR 
    RMR_method_in_respo<-ggplot(data = data.rmrER, aes(x = h_in_respo, y = resid))+
      geom_smooth(color = "black", method = "loess", se = F)+
      geom_smooth(data = data.rmr.test, method = "loess", se = F, 
                  mapping = aes(x = h_in_respo, y = resid),
                  color  = cols.rmr[3], fill  = cols.rmr[3])+
      geom_point(size = 1, alpha = 0.3, color  = "grey")+
      geom_point(data = data.rmr.test, mapping = aes(x = h_in_respo, y = resid),
                 size = 1, alpha = 1, color  = cols.rmr[3])+
      geom_hline(yintercept = 0)
    ggformat(RMR_method_in_respo, x_title="Hours in respo chamber",
             y_title=bquote(residual~RMR),
             size_text = 12, print = F, title = "All fish together")
    
    # MMR 
    MMR_method_in_respo<-ggplot(data = data.amrER, aes(x = h_in_respo, y = resid))+
      geom_smooth(color = "black", method = "loess", se = F)+
      geom_smooth(data = data.amr.test, method = "loess", se = F, 
                  mapping = aes(x = h_in_respo, y = resid),
                  color  = cols.amr[3], fill  = cols.amr[3])+
      geom_point(size = 1, alpha = 0.3, color  = "grey")+
      geom_point(data = data.amr.test, mapping = aes(x = h_in_respo, y = resid),
                 size = 1, alpha = 1, color  = cols.amr[3])+
      geom_hline(yintercept = 0)
    ggformat(MMR_method_in_respo, x_title="Hours in respo chamber",
             y_title=bquote(residual~MMR),
             size_text = 12, print = F, title = "All fish together")
    
    # resp to size ratio
    # RMR 
    RMR_method_resp_vol<-ggplot(data = data.rmrER, aes(x = (chamber_vol_L*1000)/BW_g, y = resid))+
      geom_smooth(color = "black", method = "loess", se = F)+
      geom_smooth(data = data.rmr.test, method = "loess", se = F, 
                  mapping = aes(x = (chamber_vol_L*1000)/BW_g, y = resid),
                  color  = cols.rmr[3], fill  = cols.rmr[3])+
      geom_point(size = 1, alpha = 0.3, color  = "grey")+
      geom_point(data = data.rmr.test, mapping = aes(x = (chamber_vol_L*1000)/BW_g, y = resid),
                 size = 1, alpha = 1, color  = cols.rmr[3])+
      geom_hline(yintercept = 0)
    ggformat(RMR_method_resp_vol, x_title="Respo chamber to body mass ratio",
             y_title=bquote(residual~RMR),
             size_text = 12, print = F, title = "All fish together")
    
    # MMR 
    MMR_method_resp_vol<-ggplot(data = data.amrER, aes(x = (chamber_vol_L*1000)/BW_g, y = resid))+
      geom_smooth(color = "black", method = "loess", se = F)+
      geom_smooth(data = data.amr.test, method = "loess", se = F, 
                  mapping = aes(x = (chamber_vol_L*1000)/BW_g, y = resid),
                  color  = cols.amr[3], fill  = cols.amr[3])+
      geom_point(size = 1, alpha = 0.3, color  = "grey")+
      geom_point(data = data.amr.test, mapping = aes(x = (chamber_vol_L*1000)/BW_g, y = resid),
                 size = 1, alpha = 1, color  = cols.amr[3])+
      geom_hline(yintercept = 0)
    ggformat(MMR_method_resp_vol, x_title="Respo chamber to body mass ratio",
             y_title=bquote(residual~MMR),
             size_text = 12, print = F, title = "All fish together")
    
    
    # save methods plots -----
    scaling_t<-cowplot:::plot_grid(MMR_method_in_respo, RMR_method_in_respo,
                                 MMR_method_starved, RMR_method_starved,
                                 MMR_method_resp_vol, RMR_method_resp_vol,
                                 MMR_method_er,MMR_method_w,
                                  align = "hv",
                                  axis = "l",
                                  nrow = 4,
                                  ncol = 2,
                                  labels = "AUTO",
                                 label_x = c(0.2, 0.2),
                                 label_y = c(0.895, 0.895),
                                  label_size = 12)
    # scaling
    ggsave(filename = paste("./Figures/FigureSUPPL_Review-methods.png", sep=""),
           plot=scaling_t, width = 10, height = 13, units = "in")
  }
  
  return(list(
   rmr_mod_ER,
   amr_mod_ER,
   as_mod_ER,
   fas_mod_ER,
   
   rmr_mod_W,
   amr_mod_W,
   as_mod_W,
   fas_mod_W
  ))
}

# **************************************************
# Function to update global models with ecology groupings
# standardized metric for contrasts

ecol_model_update<-function(ecol.model.null, 
                            ecol.data.subset, 
                            temp.test.category, 
                            mr.type,
                            ecol.predictor, 
                            data.BIC = NULL,
                            data.ANOVA = NULL,
                            data.EMMEANS = NULL, 
                            data.CONTRASTS = NULL, 
                            data.PARAMS = NULL,
                            ref.tempTest, # for estimate means reference grid
                            ref.lnBWg # for estimate means reference grid
                            ){

  
    for (i in 1:length(ecol.predictor)){
        
      ecol.model<-Almer(update.formula(formula(ecol.model.null), paste(' ~ . + ', ecol.predictor[i], collapse = "")), data=ecol.data.subset, REML = F)

      # print(ecol.model)# ecol.model<-Almer(lnRMR ~ lnBWg + tempTestK1000 + DemersPelag + (1 | species) + (1 | species:trial), data=data.rmr.test, REML = F)
      
      # setup BIC table with all data/model specific factors as columns 
      BIC.table<-BICdelta(BIC(ecol.model.null, ecol.model))
      BIC.table$models<-rownames(BIC.table)
      BIC.table$var<-mr.type
      BIC.table$test.categ<-temp.test.category
      BIC.table$ecol<-ecol.predictor[i]
      
      params.table<-as.data.frame(fixef(ecol.model.null))
      ecol.params.table<-as.data.frame(fixef(ecol.model))
      params.table<-as.data.frame(cbind(params.table[1:3,1], ecol.params.table[1:3,1]))
      colnames(params.table)<-c("fixef.model.null", "fixef.model.ecol")
      params.table$estimate<-rownames(params.table)
  
      params.table$var<-mr.type
      params.table$test.categ<-temp.test.category
      params.table$ecol<-ecol.predictor[i]
      
      # start new data frames if not already there
      if(is.null(data.BIC)){
        data.BIC<-as.data.frame(BIC.table)
      }else{
        data.BIC<-rbind(data.BIC, BIC.table)
      }
      
      if(is.null(data.PARAMS)){
        data.PARAMS<-as.data.frame(params.table)
      }else{
        data.PARAMS<-rbind(data.PARAMS, params.table)
      }
     
      # for interaction model run type III anova, for non-interaction models type II 
      if(has_interaction(ecol.model)){
        anova.comp<-as.data.frame(car::Anova(ecol.model, type = "III"), )
        anova.comp$anova_type<-"III, car package"
      }else{
        anova.comp<-as.data.frame(car::Anova(ecol.model, type = "II"), )
        anova.comp$anova_type<-"II, car package"
      }
   
      anova.comp$predictor<-rownames(anova.comp)
      anova.comp$var<-mr.type
      anova.comp$test.categ<-temp.test.category
      anova.comp$ecol<-ecol.predictor[i]
      
      # if not already start a new data frame
      if(is.null(data.ANOVA)){
        data.ANOVA<-as.data.frame(anova.comp)
      }else{
        data.ANOVA<-rbind(data.ANOVA, anova.comp)
      }
      
      # Get marginal estimated means and post-hoc pairwise comparisons
      if(ecol.predictor[i] == "DemersPelag"){
          ecol.ref.grid<-(ref_grid(ecol.model,
                                   at = list(tempTest = ref.tempTest, lnBWg = ref.lnBWg),
                                   data = ecol.data.subset))
  
        emm.emm<-as.data.frame(emmeans(ecol.ref.grid, pairwise ~ DemersPelag)$emmeans)
        colnames(emm.emm)[1]<-c("EcolGroup")
        emm.cont<-as.data.frame(emmeans(ecol.ref.grid, pairwise ~ DemersPelag)$contrasts)
      }
      
      if(ecol.predictor[i] == "BodyShapeI"){
          ecol.ref.grid<-(ref_grid(ecol.model,
                                   at = list(tempTest = ref.tempTest, lnBWg = ref.lnBWg),
                                   data = ecol.data.subset))
       
        emm.emm<-as.data.frame(emmeans(ecol.ref.grid, pairwise ~ BodyShapeI)$emmeans)
        colnames(emm.emm)[1]<-c("EcolGroup")
        emm.cont<-as.data.frame(emmeans(ecol.ref.grid, pairwise ~ BodyShapeI)$contrasts)
      }
      if(ecol.predictor[i] == "Climate"){
        ecol.ref.grid<-(ref_grid(ecol.model,
                                   at = list(tempTest = ref.tempTest, lnBWg = ref.lnBWg),
                                   data = ecol.data.subset))
  
        emm.emm<-as.data.frame(emmeans(ecol.ref.grid, pairwise ~ Climate)$emmeans)
        colnames(emm.emm)[1]<-c("EcolGroup")
        emm.cont<-as.data.frame(emmeans(ecol.ref.grid, pairwise ~ Climate)$contrasts)
      }
      
      if(ecol.predictor[i] == "salintyComb"){
          ecol.ref.grid<-(ref_grid(ecol.model,
                                   at = list(tempTest = ref.tempTest, lnBWg = ref.lnBWg),
                                   data = ecol.data.subset))
  
        emm.emm<-as.data.frame(emmeans(ecol.ref.grid, pairwise ~ salintyComb)$emmeans)
        colnames(emm.emm)[1]<-c("EcolGroup")
        emm.cont<-as.data.frame(emmeans(ecol.ref.grid, pairwise ~ salintyComb)$contrasts)
      }
      
       emm.emm$var<-mr.type
       emm.emm$test.categ<-temp.test.category
       emm.emm$ecol<-ecol.predictor[i]
       emm.cont$var<-mr.type
       emm.cont$test.categ<-temp.test.category
       emm.cont$ecol<-ecol.predictor[i]
        
      # if not alread existant, start new dataframes
      if(is.null(data.EMMEANS)){
        data.EMMEANS<-as.data.frame(emm.emm)
      }else{
        data.EMMEANS<-rbind(data.EMMEANS, emm.emm)
      }
      
      if(is.null(data.CONTRASTS)){
        data.CONTRASTS<-as.data.frame(emm.cont)
      }else{
        data.CONTRASTS<-rbind(data.CONTRASTS, emm.cont)
      }
    }
    
    return(list(data.BIC , data.ANOVA, data.EMMEANS, data.CONTRASTS, data.PARAMS))
}




