# this calls custom in function in 'get_data_phylo_matrix.R'. 
# IMPORTANT: ensure to receive a message: "All species names are identified and mathced with phylo data N species" that marks that all species have been matched with original dataset. 

# rmr
get_phylo_matrix(species.list = unique(levels(data.rmrER$species)), matrix.name = "A", tree.name = "tree.kk", dataset.ID = "RMR optimal")
get_phylo_matrix(species.list = unique(levels(data.rmr.test$species)), matrix.name = "A.rmr.w", tree.name = "tr.rmr.w", dataset.ID = "RMR warm")
# amr 
get_phylo_matrix(species.list = unique(levels(data.amrER$species)), matrix.name = "A.mmr.er", tree.name = "tr.mmr.er", dataset.ID = "MMR optimal")
get_phylo_matrix(species.list = unique(levels(data.amr.test$species)), matrix.name = "A.mmr.w", tree.name = "tr.mmr.w", dataset.ID = "MMR warm")
# aas
get_phylo_matrix(species.list = unique(levels(data.asER$species)), matrix.name = "A.aas.er", tree.name = "tr.aas.er", dataset.ID = "AAS optimal")
get_phylo_matrix(species.list = unique(levels(data.as.test$species)), matrix.name = "A.aas.w", tree.name = "tr.aas.w", dataset.ID = "AAS warm")
# fas
get_phylo_matrix(species.list = unique(levels(data.fasER$species)), matrix.name = "A.fas.er", tree.name = "tr.fas.er", dataset.ID = "FAS optimal")
get_phylo_matrix(species.list = unique(levels(data.fas.test$species)), matrix.name = "A.fas.w", tree.name = "tree.fas.w", dataset.ID = "FAS warm")

# running take ~ on my personal comupeter

  ## RMR optimal -----------------
  # Phylo_RMR_model0 <- Almer(lnRMR ~ lnBWg + tempTest + (1|species) , data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model0.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model0intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model0int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) , data=data.rmrER, REML=FALSE, A = list(species = A))
  # 
  # Phylo_RMR_model1 <- Almer(lnRMR ~ lnBWg + tempTest + (1|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model1.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model1intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model1int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # 
  Phylo_RMR_model2 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model2int <- Almer(lnRMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model2.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model2intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE)  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # 
  Phylo_RMR_model4 <- Almer(lnRMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model4int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model4.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model4intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # 
  Phylo_RMR_model5 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model5int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model5.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model5intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # 
  # Phylo_RMR_model6 <- Almer(lnRMR ~ lnBWg + tempTest + (lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model6int <- Almer(lnRMR ~ lnBWg * tempTest+ (lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model6.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model6intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # 
  ## RMR optimal JUVENILE INCLUDED-----------------
  # Phylo_RMR_model0.j <- Almer(lnRMR ~ lnBWg + tempTest + lifestage + (1|species) , data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model0.POLY.j <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model0intPOLY.j <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage +  (1|species), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model0int.j <- Almer(lnRMR ~ lnBWg * tempTest + lifestage + (1|species) , data=data.rmrER, REML=FALSE, A = list(species = A))
  # 
  # Phylo_RMR_model1.j <- Almer(lnRMR ~ lnBWg + tempTest + lifestage + (1|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model1.POLY.j <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model1intPOLY.j <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model1int.j <- Almer(lnRMR ~ lnBWg * tempTest + lifestage + (1|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # 
  Phylo_RMR_model2.j <- Almer(lnRMR ~ lnBWg + tempTest  + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model2int.j <- Almer(lnRMR ~ lnBWg * tempTest  + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model2.POLY.j <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model2intPOLY.j <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE)  + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # # 
  Phylo_RMR_model4.j <- Almer(lnRMR ~ lnBWg + tempTest + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model4int.j <- Almer(lnRMR ~ lnBWg * tempTest + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model4.POLY.j <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model4intPOLY.j <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # 
  Phylo_RMR_model5.j <- Almer(lnRMR ~ lnBWg + tempTest  + lifestage + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model5int.j <- Almer(lnRMR ~ lnBWg * tempTest + lifestage + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model5.POLY.j <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model5intPOLY.j <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # # 
  # Phylo_RMR_model6.j <- Almer(lnRMR ~ lnBWg + tempTest + lifestage + (lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model6int.j <- Almer(lnRMR ~ lnBWg * tempTest+ lifestage + (lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model6.POLY.j <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model6intPOLY.j <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  
    ## RMR optimal Life stage X scaling -----------------
  # Phylo_RMR_model0.jint <- Almer(lnRMR ~ lnBWg * lifestage + tempTest + (1|species) , data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model0.POLY.jint <- Almer(lnRMR ~ lnBWg * lifestage+ poly(tempTest,2, raw = TRUE) +  (1|species), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model0intPOLY.jint <- Almer(lnRMR ~ lnBWg * lifestage* poly(tempTest,2, raw = TRUE) +   (1|species), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model0int.jint <- Almer(lnRMR ~ lnBWg * lifestage* tempTest +  (1|species) , data=data.rmrER, REML=FALSE, A = list(species = A))
  # 
  # Phylo_RMR_model1.jint <- Almer(lnRMR ~ lnBWg * lifestage+ tempTest +  (1|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model1.POLY.jint <- Almer(lnRMR ~ lnBWg * lifestage+ poly(tempTest,2, raw = TRUE) +  (1|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model1intPOLY.jint <- Almer(lnRMR ~ lnBWg * lifestage* poly(tempTest,2, raw = TRUE) +  (1|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model1int.jint <- Almer(lnRMR ~ lnBWg * lifestage* tempTest +  (1|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # 
  Phylo_RMR_model2.jint <- Almer(lnRMR ~ lnBWg * lifestage+ tempTest  +  (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model2int.jint <- Almer(lnRMR ~ lnBWg * lifestage* tempTest  +  (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model2.POLY.jint <- Almer(lnRMR ~ lnBWg * lifestage+ poly(tempTest,2, raw = TRUE) +  (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model2intPOLY.jint <- Almer(lnRMR ~ lnBWg * lifestage* poly(tempTest,2, raw = TRUE)  +  (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # 
  Phylo_RMR_model4.jint <- Almer(lnRMR ~ lnBWg * lifestage+ tempTest +  (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model4int.jint <- Almer(lnRMR ~ lnBWg * lifestage* tempTest +  (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model4.POLY.jint <- Almer(lnRMR ~ lnBWg * lifestage+ poly(tempTest,2, raw = TRUE) +  (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model4intPOLY.jint <- Almer(lnRMR ~ lnBWg * lifestage* poly(tempTest,2, raw = TRUE) +  (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # 
  Phylo_RMR_model5.jint <- Almer(lnRMR ~ lnBWg * lifestage+ tempTest  +  (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  Phylo_RMR_model5int.jint <- Almer(lnRMR ~ lnBWg * lifestage* tempTest +  (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model5.POLY.jint <- Almer(lnRMR ~ lnBWg * lifestage+ poly(tempTest,2, raw = TRUE) +  (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model5intPOLY.jint <- Almer(lnRMR ~ lnBWg * lifestage* poly(tempTest,2, raw = TRUE) +  (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # 
  # Phylo_RMR_model6.jint <- Almer(lnRMR ~ lnBWg * lifestage+ tempTest +  (lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model6int.jint <- Almer(lnRMR ~ lnBWg * lifestage* tempTest+  (lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model6.POLY.jint <- Almer(lnRMR ~ lnBWg * lifestage+ poly(tempTest,2, raw = TRUE) +  (lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  # Phylo_RMR_model6intPOLY.jint <- Almer(lnRMR ~ lnBWg * lifestage* poly(tempTest,2, raw = TRUE) +  (lnBWg|species) + (1|species:trial), data=data.rmrER, REML=FALSE, A = list(species = A))
  
  
  # No singular fits 
  ### BIC rmr optimal --------------
  RMR_BIC<-BICdelta(BIC(
                  # Phylo_RMR_model0, Phylo_RMR_model0int, Phylo_RMR_model0.POLY, Phylo_RMR_model0intPOLY,
                  # Phylo_RMR_model1, Phylo_RMR_model1int, Phylo_RMR_model1.POLY, Phylo_RMR_model1intPOLY,
                  Phylo_RMR_model2, Phylo_RMR_model2int,# Phylo_RMR_model2.POLY, Phylo_RMR_model2intPOLY,
                  Phylo_RMR_model4, Phylo_RMR_model4int, #Phylo_RMR_model4.POLY, Phylo_RMR_model4intPOLY,
                  Phylo_RMR_model5, Phylo_RMR_model5int, #Phylo_RMR_model5.POLY, Phylo_RMR_model5intPOLY,
                  # Phylo_RMR_model6, Phylo_RMR_model6int, Phylo_RMR_model6.POLY, Phylo_RMR_model6intPOLY,
                  # lifestage 
                  # Phylo_RMR_model0.j, Phylo_RMR_model0int.j, Phylo_RMR_model0.POLY.j, Phylo_RMR_model0intPOLY.j,
                  # Phylo_RMR_model1.j, Phylo_RMR_model1int.j, Phylo_RMR_model1.POLY.j, Phylo_RMR_model1intPOLY.j,
                  Phylo_RMR_model2.j, Phylo_RMR_model2int.j,# Phylo_RMR_model2.POLY.j, Phylo_RMR_model2intPOLY.j,
                  Phylo_RMR_model4.j, Phylo_RMR_model4int.j, #Phylo_RMR_model4.POLY.j, Phylo_RMR_model4intPOLY.j,
                  Phylo_RMR_model5.j, Phylo_RMR_model5int.j, #Phylo_RMR_model5.POLY.j, Phylo_RMR_model5intPOLY.j,
                  # Phylo_RMR_model6.j, Phylo_RMR_model6int.j, Phylo_RMR_model6.POLY.j, Phylo_RMR_model6intPOLY.j, 
                  # lifestage x scaling
                  # Phylo_RMR_model0.jint, Phylo_RMR_model0int.jint, Phylo_RMR_model0.POLY.jint, Phylo_RMR_model0intPOLY.jint,
                  # Phylo_RMR_model1.jint, Phylo_RMR_model1int.jint, Phylo_RMR_model1.POLY.jint, Phylo_RMR_model1intPOLY.jint,
                  Phylo_RMR_model2.jint, Phylo_RMR_model2int.jint,# Phylo_RMR_model2.POLY.jint, Phylo_RMR_model2intPOLY.jint,
                  Phylo_RMR_model4.jint, Phylo_RMR_model4int.jint, #Phylo_RMR_model4.POLY.jint, Phylo_RMR_model4intPOLY.jint,
                  Phylo_RMR_model5.jint, Phylo_RMR_model5int.jint))#, Phylo_RMR_model5.POLY.jint, Phylo_RMR_model5intPOLY.jint))
                  # Phylo_RMR_model6.jint, Phylo_RMR_model6int.jint, Phylo_RMR_model6.POLY.jint, Phylo_RMR_model6intPOLY.jint))
  
  ## MMR optimal -----------------
  # Phylo_MMR_model0 <- Almer(lnAMR ~ lnBWg + tempTest + (1|species) , data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model0.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model0intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model0int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) , data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # 
  # Phylo_MMR_model1 <- Almer(lnAMR ~ lnBWg + tempTest + (1|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model1.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model1intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model1int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))

  Phylo_MMR_model2 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model2int <- Almer(lnAMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model2.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model2intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE)  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # # 
  Phylo_MMR_model4 <- Almer(lnAMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model4int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model4.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model4intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # 
  Phylo_MMR_model5 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model5int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model5.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model5intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # # 
  # Phylo_MMR_model6 <- Almer(lnAMR ~ lnBWg + tempTest + (lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model6int <- Almer(lnAMR ~ lnBWg * tempTest+ (lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model6.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model6intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  
  ## MMR optimal juveniles -----------------
  #   Phylo_MMR_model0.j<- Almer(lnAMR ~ lnBWg + tempTest + lifestage + (1|species) , data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model0.POLY.j <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model0intPOLY.j <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model0int.j <- Almer(lnAMR ~ lnBWg * tempTest + lifestage + (1|species) , data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # 
  # Phylo_MMR_model1.j <- Almer(lnAMR ~ lnBWg + tempTest + lifestage + (1|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model1.POLY.j <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model1intPOLY.j <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model1int.j <- Almer(lnAMR ~ lnBWg * tempTest + lifestage + (1|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # 
  Phylo_MMR_model2.j <- Almer(lnAMR ~ lnBWg + tempTest  + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model2int.j <- Almer(lnAMR ~ lnBWg * tempTest  + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model2.POLY.j <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model2intPOLY.j <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE)  + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # # 
  Phylo_MMR_model4.j <- Almer(lnAMR ~ lnBWg + tempTest + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model4int.j <- Almer(lnAMR ~ lnBWg * tempTest + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model4.POLY.j <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model4intPOLY.j <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # 
  Phylo_MMR_model5.j <- Almer(lnAMR ~ lnBWg + tempTest  + lifestage + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model5int.j <- Almer(lnAMR ~ lnBWg * tempTest + lifestage + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model5.POLY.j <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model5intPOLY.j <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # # 
  # Phylo_MMR_model6.j <- Almer(lnAMR ~ lnBWg + tempTest + lifestage + (lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model6int.j <- Almer(lnAMR ~ lnBWg * tempTest+ lifestage + (lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model6.POLY.j <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model6intPOLY.j <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  
  # MMR lifestage X scaling ------------------------
  Phylo_MMR_model2.jint <- Almer(lnAMR ~ lnBWg * lifestage  + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model2int.jint <- Almer(lnAMR ~ lnBWg * lifestage  * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model2.POLY.jint <- Almer(lnAMR ~ lnBWg * lifestage  + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model2intPOLY.jint <- Almer(lnAMR ~ lnBWg * lifestage  * poly(tempTest,2, raw = TRUE)  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # 
  Phylo_MMR_model4.jint <- Almer(lnAMR ~ lnBWg * lifestage  + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model4int.jint <- Almer(lnAMR ~ lnBWg * lifestage  * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model4.POLY.jint <- Almer(lnAMR ~ lnBWg * lifestage  + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model4intPOLY.jint <- Almer(lnAMR ~ lnBWg * lifestage  * poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # 
  Phylo_MMR_model5.jint <- Almer(lnAMR ~ lnBWg * lifestage  + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  Phylo_MMR_model5int.jint <- Almer(lnAMR ~ lnBWg * lifestage  * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model5.POLY.jint <- Almer(lnAMR ~ lnBWg * lifestage  + poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))
  # Phylo_MMR_model5intPOLY.jint <- Almer(lnAMR ~ lnBWg * lifestage  * poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amrER, REML=FALSE, A = list(species = A.mmr.er))

  ### BIC --------------
  MMR_BIC<-BICdelta(BIC(
                  # Phylo_MMR_model0, Phylo_MMR_model0int, Phylo_MMR_model0.POLY, Phylo_MMR_model0intPOLY,
                  # Phylo_MMR_model1, Phylo_MMR_model1int, Phylo_MMR_model1.POLY, Phylo_MMR_model1intPOLY,
                  Phylo_MMR_model2, Phylo_MMR_model2int,# Phylo_MMR_model2.POLY, Phylo_MMR_model2intPOLY,
                  Phylo_MMR_model4, Phylo_MMR_model4int, #Phylo_MMR_model4.POLY, Phylo_MMR_model4intPOLY,
                  Phylo_MMR_model5, Phylo_MMR_model5int, #Phylo_MMR_model5.POLY, Phylo_MMR_model5intPOLY,
                  # Phylo_MMR_model6, Phylo_MMR_model6int, Phylo_MMR_model6.POLY, Phylo_MMR_model6intPOLY,
                  # # juveniles
                  # Phylo_MMR_model0.j, Phylo_MMR_model0int.j, Phylo_MMR_model0.POLY.j, Phylo_MMR_model0intPOLY.j,
                  # Phylo_MMR_model1.j, Phylo_MMR_model1int.j, Phylo_MMR_model1.POLY.j, Phylo_MMR_model1intPOLY.j,
                  Phylo_MMR_model2.j, Phylo_MMR_model2int.j,# Phylo_MMR_model2.POLY.j, Phylo_MMR_model2intPOLY.j,
                  Phylo_MMR_model4.j, Phylo_MMR_model4int.j, #Phylo_MMR_model4.POLY.j, Phylo_MMR_model4intPOLY.j,
                  Phylo_MMR_model5.j, Phylo_MMR_model5int.j, #Phylo_MMR_model5.POLY.j, Phylo_MMR_model5intPOLY.j,
                  # Phylo_MMR_model6.j, Phylo_MMR_model6int.j, Phylo_MMR_model6.POLY.j, Phylo_MMR_model6intPOLY.j
                  # juvenile x sacling
                  Phylo_MMR_model2.jint, Phylo_MMR_model2int.jint,# Phylo_MMR_model2.POLY.jint, Phylo_MMR_model2intPOLY.jint,
                  Phylo_MMR_model4.jint, Phylo_MMR_model4int.jint, #Phylo_MMR_model4.POLY.jint, Phylo_MMR_model4intPOLY.jint,
                  Phylo_MMR_model5.jint, Phylo_MMR_model5int.jint))#, Phylo_MMR_model5.POLY.jint, Phylo_MMR_model5intPOLY.jint
                
  
  ## AAS optimal ------------------
  # Phylo_AS_model0 <- Almer(lnAS ~ lnBWg + tempTest + (1|species) , data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model0.POLY <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model0intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model0int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) , data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # 
  # Phylo_AS_model1 <- Almer(lnAS ~ lnBWg + tempTest + (1|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model1.POLY <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model1intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model1int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # 
  Phylo_AS_model2 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model2int <- Almer(lnAS ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model2.POLY <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model2intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE)  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # # 
  Phylo_AS_model4 <- Almer(lnAS ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model4int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model4.POLY <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model4intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # 
  Phylo_AS_model5 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model5int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model5.POLY <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model5intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # # 
  # Phylo_AS_model6 <- Almer(lnAS ~ lnBWg + tempTest + (lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model6int <- Almer(lnAS ~ lnBWg * tempTest+ (lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model6.POLY <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model6intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  
  ## AAS optimal with lifestage ------------------
  # Phylo_AS_model0.j <- Almer(lnAS ~ lnBWg + tempTest + lifestage + (1|species) , data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model0.POLY.j <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model0intPOLY.j <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model0int.j <- Almer(lnAS ~ lnBWg * tempTest + lifestage + (1|species) , data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # 
  # Phylo_AS_model1.j <- Almer(lnAS ~ lnBWg + tempTest + lifestage + (1|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model1.POLY.j <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model1intPOLY.j <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model1int.j <- Almer(lnAS ~ lnBWg * tempTest + lifestage + (1|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # 
  Phylo_AS_model2.j <- Almer(lnAS ~ lnBWg + tempTest  + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model2int.j <- Almer(lnAS ~ lnBWg * tempTest  + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model2.POLY.j <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model2intPOLY.j <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE)  + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # # 
  Phylo_AS_model4.j <- Almer(lnAS ~ lnBWg + tempTest + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model4int.j <- Almer(lnAS ~ lnBWg * tempTest + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model4.POLY.j <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model4intPOLY.j <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # 
  Phylo_AS_model5.j <- Almer(lnAS ~ lnBWg + tempTest  + lifestage + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model5int.j <- Almer(lnAS ~ lnBWg * tempTest + lifestage + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model5.POLY.j <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model5intPOLY.j <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # # 
  # Phylo_AS_model6.j <- Almer(lnAS ~ lnBWg + tempTest + lifestage + (lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model6int.j <- Almer(lnAS ~ lnBWg * tempTest+ lifestage + (lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model6.POLY.j <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model6intPOLY.j <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # 
  # AS lifestage X scaling --------------
  Phylo_AS_model2.jint <- Almer(lnAS ~ lnBWg * lifestage  + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model2int.jint <- Almer(lnAS ~ lnBWg * lifestage  * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model2.POLY.jint <- Almer(lnAS ~ lnBWg * lifestage  + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model2intPOLY.jint <- Almer(lnAS ~ lnBWg * lifestage  * poly(tempTest,2, raw = TRUE)  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # 
  Phylo_AS_model4.jint <- Almer(lnAS ~ lnBWg * lifestage + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model4int.jint <- Almer(lnAS ~ lnBWg * lifestage * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model4.POLY.jint <- Almer(lnAS ~ lnBWg * lifestage + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model4intPOLY.jint <- Almer(lnAS ~ lnBWg * lifestage * poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # 
  Phylo_AS_model5.jint <- Almer(lnAS ~ lnBWg * lifestage + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  Phylo_AS_model5int.jint <- Almer(lnAS ~ lnBWg * lifestage * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model5.POLY.jint <- Almer(lnAS ~ lnBWg * lifestage + poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # Phylo_AS_model5intPOLY.jint <- Almer(lnAS ~ lnBWg * lifestage * poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.asER, REML=FALSE, A = list(species = A.aas.er))
  # 
  ### BIC --------
  AS_BIC<-BICdelta(BIC(#Phylo_AS_model0, Phylo_AS_model0int, Phylo_AS_model0.POLY, Phylo_AS_model0intPOLY,
                  # Phylo_AS_model1, Phylo_AS_model1int, Phylo_AS_model1.POLY, Phylo_AS_model1intPOLY,
                  Phylo_AS_model2, Phylo_AS_model2int, #Phylo_AS_model2.POLY, Phylo_AS_model2intPOLY,
                  Phylo_AS_model4, Phylo_AS_model4int, #Phylo_AS_model4.POLY, Phylo_AS_model4intPOLY,
                  Phylo_AS_model5, Phylo_AS_model5int, #Phylo_AS_model5.POLY, Phylo_AS_model5intPOLY,
                  # Phylo_AS_model6, Phylo_AS_model6int, Phylo_AS_model6.POLY, Phylo_AS_model6intPOLY, 
                  # juveniles 
                  # Phylo_AS_model0.j, Phylo_AS_model0int.j, Phylo_AS_model0.POLY.j, Phylo_AS_model0intPOLY.j,
                  # Phylo_AS_model1.j, Phylo_AS_model1int.j, Phylo_AS_model1.POLY.j, Phylo_AS_model1intPOLY.j,
                  Phylo_AS_model2.j, Phylo_AS_model2int.j,# Phylo_AS_model2.POLY.j, Phylo_AS_model2intPOLY.j,
                  Phylo_AS_model4.j, Phylo_AS_model4int.j, #Phylo_AS_model4.POLY.j, Phylo_AS_model4intPOLY.j,
                  Phylo_AS_model5.j, Phylo_AS_model5int.j,# Phylo_AS_model5.POLY.j, Phylo_AS_model5intPOLY.j,
                  # Phylo_AS_model6.j, Phylo_AS_model6int.j, Phylo_AS_model6.POLY.j, Phylo_AS_model6intPOLY.j
                  Phylo_AS_model2.jint, Phylo_AS_model2int.jint, #Phylo_AS_model2.POLY.jint, Phylo_AS_model2intPOLY.jint,
                  Phylo_AS_model4.jint, Phylo_AS_model4int.jint, #Phylo_AS_model4.POLY.jint, Phylo_AS_model4intPOLY.jint,
                  Phylo_AS_model5.jint, Phylo_AS_model5int.jint))#, Phylo_AS_model5.POLY.jint, Phylo_AS_model5intPOLY.jint))
  
  
  ## FAS optimal ----------------
  # Phylo_FAS_model0 <- Almer(log(FAS) ~ lnBWg + tempTest + (1|species) , data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model0.POLY <- Almer(log(FAS) ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model0intPOLY <- Almer(log(FAS) ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model0int <- Almer(log(FAS) ~ lnBWg * tempTest + (1|species) , data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # 
  # Phylo_FAS_model1 <- Almer(log(FAS) ~ lnBWg + tempTest + (1|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model1.POLY <- Almer(log(FAS) ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model1intPOLY <- Almer(log(FAS) ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model1int <- Almer(log(FAS) ~ lnBWg * tempTest + (1|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # 
  Phylo_FAS_model2 <- Almer(log(FAS) ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model2int <- Almer(log(FAS) ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model2.POLY <- Almer(log(FAS) ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model2intPOLY <- Almer(log(FAS) ~ lnBWg * poly(tempTest,2, raw = TRUE)  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # # 
  Phylo_FAS_model4 <- Almer(log(FAS) ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model4int <- Almer(log(FAS) ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model4.POLY <- Almer(log(FAS) ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model4intPOLY <- Almer(log(FAS) ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # 
  Phylo_FAS_model5 <- Almer(log(FAS) ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model5int <- Almer(log(FAS) ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model5.POLY <- Almer(log(FAS) ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model5intPOLY <- Almer(log(FAS) ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # # 
  # Phylo_FAS_model6 <- Almer(lnFAS ~ lnBWg + tempTest + (lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model6int <- Almer(lnFAS ~ lnBWg * tempTest+ (lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model6.POLY <- Almer(lnFAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model6intPOLY <- Almer(lnFAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  
  ## FAS optimal with lifestage----------------
  # Phylo_FAS_model0.j <- Almer(log(FAS) ~ lnBWg + tempTest + lifestage + (1|species) , data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model0.POLY.j <- Almer(log(FAS) ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model0intPOLY.j <- Almer(log(FAS) ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model0int.j <- Almer(log(FAS) ~ lnBWg * tempTest + lifestage + (1|species) , data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # 
  # Phylo_FAS_model1.j <- Almer(log(FAS) ~ lnBWg + tempTest + lifestage + (1|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model1.POLY.j <- Almer(log(FAS) ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model1intPOLY.j <- Almer(log(FAS) ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model1int.j <- Almer(log(FAS) ~ lnBWg * tempTest + (1|species) + lifestage + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # 
  Phylo_FAS_model2.j <- Almer(log(FAS) ~ lnBWg + tempTest  + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model2int.j <- Almer(log(FAS) ~ lnBWg * tempTest  + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model2.POLY.j <- Almer(log(FAS) ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model2intPOLY.j <- Almer(log(FAS) ~ lnBWg * poly(tempTest,2, raw = TRUE)  + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # # 
  Phylo_FAS_model4.j <- Almer(log(FAS) ~ lnBWg + tempTest + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model4int.j <- Almer(log(FAS) ~ lnBWg * tempTest + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model4.POLY.j <- Almer(log(FAS) ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model4intPOLY.j <- Almer(log(FAS) ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # # 
  Phylo_FAS_model5.j <- Almer(log(FAS) ~ lnBWg + tempTest  + lifestage + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model5int.j <- Almer(log(FAS) ~ lnBWg * tempTest + lifestage + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model5.POLY.j <- Almer(log(FAS) ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model5intPOLY.j <- Almer(log(FAS) ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # # 
  # Phylo_FAS_model6.j <- Almer(lnFAS ~ lnBWg + tempTest + lifestage + (lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model6int.j <- Almer(lnFAS ~ lnBWg * tempTest+ lifestage + (lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model6.POLY.j <- Almer(lnFAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model6intPOLY.j <- Almer(lnFAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # 
  # FAS lifestage X scaling -----------------
  Phylo_FAS_model2.jint <- Almer(log(FAS) ~ lnBWg * lifestage  + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model2int.jint <- Almer(log(FAS) ~ lnBWg * lifestage  * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model2.POLY.jint <- Almer(log(FAS) ~ lnBWg * lifestage  + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model2intPOLY.jint <- Almer(log(FAS) ~ lnBWg * lifestage  * poly(tempTest,2, raw = TRUE)  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # 
  Phylo_FAS_model4.jint <- Almer(log(FAS) ~ lnBWg * lifestage  + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model4int.jint <- Almer(log(FAS) ~ lnBWg * lifestage * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model4.POLY.jint <- Almer(log(FAS) ~ lnBWg * lifestage + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model4intPOLY.jint <- Almer(log(FAS) ~ lnBWg * lifestage * poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # 
  Phylo_FAS_model5.jint <- Almer(log(FAS) ~ lnBWg * lifestage + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  Phylo_FAS_model5int.jint <- Almer(log(FAS) ~ lnBWg * lifestage * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model5.POLY.jint <- Almer(log(FAS) ~ lnBWg * lifestage + poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # Phylo_FAS_model5intPOLY.jint <- Almer(log(FAS) ~ lnBWg * lifestage * poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.fasER, REML=FALSE, A = list(species = A.fas.er))
  # 
  ### BIC -------------
  FAS_BIC<-BICdelta(BIC(#Phylo_FAS_model0, Phylo_FAS_model0int, Phylo_FAS_model0.POLY, Phylo_FAS_model0intPOLY,
                  # Phylo_FAS_model1, Phylo_FAS_model1int, Phylo_FAS_model1.POLY, Phylo_FAS_model1intPOLY,
                  Phylo_FAS_model2, Phylo_FAS_model2int,# Phylo_FAS_model2.POLY, Phylo_FAS_model2intPOLY,
                  Phylo_FAS_model4, Phylo_FAS_model4int, #Phylo_FAS_model4.POLY, Phylo_FAS_model4intPOLY,
                  Phylo_FAS_model5, Phylo_FAS_model5int, #Phylo_FAS_model5.POLY, Phylo_FAS_model5intPOLY,
                  # Phylo_FAS_model6, Phylo_FAS_model6int, Phylo_FAS_model6.POLY, Phylo_FAS_model6intPOLY,
                  # juveniles 
                  # Phylo_FAS_model0.j, Phylo_FAS_model0int.j, Phylo_FAS_model0.POLY.j, Phylo_FAS_model0intPOLY.j,
                  # Phylo_FAS_model1.j, Phylo_FAS_model1int.j, Phylo_FAS_model1.POLY.j, Phylo_FAS_model1intPOLY.j,
                  Phylo_FAS_model2.j, Phylo_FAS_model2int.j,# Phylo_FAS_model2.POLY.j, Phylo_FAS_model2intPOLY.j,
                  Phylo_FAS_model4.j, Phylo_FAS_model4int.j,# Phylo_FAS_model4.POLY.j, Phylo_FAS_model4intPOLY.j,
                  Phylo_FAS_model5.j, Phylo_FAS_model5int.j,# Phylo_FAS_model5.POLY.j, Phylo_FAS_model5intPOLY.j,
                  # Phylo_FAS_model6.j, Phylo_FAS_model6int.j, Phylo_FAS_model6.POLY.j, Phylo_FAS_model6intPOLY.j,
                  Phylo_FAS_model2.jint, Phylo_FAS_model2int.jint, #Phylo_FAS_model2.POLY.jint, Phylo_FAS_model2intPOLY.jint,
                  Phylo_FAS_model4.jint, Phylo_FAS_model4int.jint, #Phylo_FAS_model4.POLY.jint, Phylo_FAS_model4intPOLY.jint,
                  Phylo_FAS_model5.jint, Phylo_FAS_model5int.jint))#, Phylo_FAS_model5.POLY.jint, Phylo_FAS_model5intPOLY.jint))
  
  
  ## RMR warm ------------------
  # Phylo_RMR_W_model0 <- Almer(lnRMR ~ lnBWg + tempTest + (1|species) , data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model0.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model0intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model0int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) , data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # 
  # Phylo_RMR_W_model1 <- Almer(lnRMR ~ lnBWg + tempTest + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model1.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model1intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model1int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # 
  Phylo_RMR_W_model2 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model2int <- Almer(lnRMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model2.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model2intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE)  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # # 
  Phylo_RMR_W_model4 <- Almer(lnRMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model4int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model4.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model4intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # 
  Phylo_RMR_W_model5 <- Almer(lnRMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model5int <- Almer(lnRMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model5.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model5intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # # 
  # Phylo_RMR_W_model6 <- Almer(lnRMR ~ lnBWg + tempTest + (lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model6int <- Almer(lnRMR ~ lnBWg * tempTest+ (lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model6.POLY <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model6intPOLY <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # 
  ## RMR juvenile warm ------------------
  # Phylo_RMR_W_model0.j <- Almer(lnRMR ~ lnBWg + tempTest + lifestage + (1|species) , data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model0.POLY.j <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model0intPOLY.j <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model0int.j <- Almer(lnRMR ~ lnBWg * tempTest + lifestage + (1|species) , data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # 
  # Phylo_RMR_W_model1.j <- Almer(lnRMR ~ lnBWg + tempTest + lifestage + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model1.POLY.j <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model1intPOLY.j <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model1int.j <- Almer(lnRMR ~ lnBWg * tempTest + lifestage + (1|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # 
  Phylo_RMR_W_model2.j <- Almer(lnRMR ~ lnBWg + tempTest  + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model2int.j <- Almer(lnRMR ~ lnBWg * tempTest  + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model2.POLY.j <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model2intPOLY.j <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE)  + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))

  Phylo_RMR_W_model4.j <- Almer(lnRMR ~ lnBWg + tempTest + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model4int.j <- Almer(lnRMR ~ lnBWg * tempTest + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model4.POLY.j <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model4intPOLY.j <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # 
  Phylo_RMR_W_model5.j <- Almer(lnRMR ~ lnBWg + tempTest  + lifestage + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model5int.j <- Almer(lnRMR ~ lnBWg * tempTest + lifestage + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model5.POLY.j <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model5intPOLY.j <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # # 
  # Phylo_RMR_W_model6.j <- Almer(lnRMR ~ lnBWg + tempTest + lifestage + (lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model6int.j <- Almer(lnRMR ~ lnBWg * tempTest+ lifestage + (lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model6.POLY.j <- Almer(lnRMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model6intPOLY.j <- Almer(lnRMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # 
  # 
  # RMR warm lifestage X scaling -----------------
  Phylo_RMR_W_model2.jint <- Almer(lnRMR ~ lnBWg * lifestage  + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model2int.jint <- Almer(lnRMR ~ lnBWg * lifestage  * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model2.POLY.jint <- Almer(lnRMR ~ lnBWg * lifestage  + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model2intPOLY.jint <- Almer(lnRMR ~ lnBWg * lifestage  * poly(tempTest,2, raw = TRUE)  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # 
  Phylo_RMR_W_model4.jint <- Almer(lnRMR ~ lnBWg * lifestage + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model4int.jint <- Almer(lnRMR ~ lnBWg * lifestage * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model4.POLY.jint <- Almer(lnRMR ~ lnBWg * lifestage + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model4intPOLY.jint <- Almer(lnRMR ~ lnBWg * lifestage * poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # 
  Phylo_RMR_W_model5.jint <- Almer(lnRMR ~ lnBWg * lifestage + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  Phylo_RMR_W_model5int.jint <- Almer(lnRMR ~ lnBWg * lifestage * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model5.POLY.jint <- Almer(lnRMR ~ lnBWg * lifestage + poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # Phylo_RMR_W_model5intPOLY.jint <- Almer(lnRMR ~ lnBWg * lifestage * poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.rmr.test, REML=FALSE, A = list(species = A.rmr.w))
  # 
  ## BIC -------------
  RMR_W_BIC<-BICdelta(BIC(#Phylo_RMR_W_model0, Phylo_RMR_W_model0int, Phylo_RMR_W_model0.POLY, Phylo_RMR_W_model0intPOLY,
                  # Phylo_RMR_W_model1, Phylo_RMR_W_model1int, Phylo_RMR_W_model1.POLY, Phylo_RMR_W_model1intPOLY,
                  Phylo_RMR_W_model2, Phylo_RMR_W_model2int, #Phylo_RMR_W_model2.POLY, Phylo_RMR_W_model2intPOLY,
                  Phylo_RMR_W_model4, Phylo_RMR_W_model4int, #Phylo_RMR_W_model4.POLY, Phylo_RMR_W_model4intPOLY,
                  Phylo_RMR_W_model5, Phylo_RMR_W_model5int, #Phylo_RMR_W_model5.POLY, Phylo_RMR_W_model5intPOLY,
                  # Phylo_RMR_W_model6, Phylo_RMR_W_model6int, Phylo_RMR_W_model6.POLY, Phylo_RMR_W_model6intPOLY,
                  # juvenile
                  # Phylo_RMR_W_model0.j,  Phylo_RMR_W_model0int.j,  Phylo_RMR_W_model0.POLY.j,  Phylo_RMR_W_model0intPOLY.j, 
                  # Phylo_RMR_W_model1.j,  Phylo_RMR_W_model1int.j,  Phylo_RMR_W_model1.POLY.j,  Phylo_RMR_W_model1intPOLY.j, 
                  Phylo_RMR_W_model2.j,  Phylo_RMR_W_model2int.j, # Phylo_RMR_W_model2.POLY.j,  Phylo_RMR_W_model2intPOLY.j,
                  Phylo_RMR_W_model4.j,  Phylo_RMR_W_model4int.j, # Phylo_RMR_W_model4.POLY.j,  Phylo_RMR_W_model4intPOLY.j, 
                  Phylo_RMR_W_model5.j,  Phylo_RMR_W_model5int.j, # Phylo_RMR_W_model5.POLY.j,  Phylo_RMR_W_model5intPOLY.j, 
                  Phylo_RMR_W_model2.jint,  Phylo_RMR_W_model2int.jint,  #Phylo_RMR_W_model2.POLY.jint,  Phylo_RMR_W_model2intPOLY.jint,
                  Phylo_RMR_W_model4.jint,  Phylo_RMR_W_model4int.jint,  #Phylo_RMR_W_model4.POLY.jint,  Phylo_RMR_W_model4intPOLY.jint, 
                  Phylo_RMR_W_model5.jint,  Phylo_RMR_W_model5int.jint))#,  Phylo_RMR_W_model5.POLY.jint,  Phylo_RMR_W_model5intPOLY.jint))
                  # Phylo_RMR_W_model6.j,  Phylo_RMR_W_model6int.j,  Phylo_RMR_W_model6.POLY.j,  Phylo_RMR_W_model6intPOLY.j))
  
  # AMR / warm temps --------------
  # Phylo_MMR_W_model0 <- Almer(lnAMR ~ lnBWg + tempTest + (1|species) , data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model0.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model0intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model0int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) , data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # 
  # Phylo_MMR_W_model1 <- Almer(lnAMR ~ lnBWg + tempTest + (1|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model1.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model1intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model1int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))

  Phylo_MMR_W_model2 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model2int <- Almer(lnAMR ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model2.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model2intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE)  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # # 
  Phylo_MMR_W_model4 <- Almer(lnAMR ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model4int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model4.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model4intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # 
  Phylo_MMR_W_model5 <- Almer(lnAMR ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model5int <- Almer(lnAMR ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model5.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model5intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # # 
  # Phylo_MMR_W_model6 <- Almer(lnAMR ~ lnBWg + tempTest + (lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model6int <- Almer(lnAMR ~ lnBWg * tempTest+ (lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model6.POLY <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model6intPOLY <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # 
  # AMR / warm JUVENILE temps --------------
  # Phylo_MMR_W_model0.j <- Almer(lnAMR ~ lnBWg + tempTest + lifestage + (1|species) , data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model0.POLY.j <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model0intPOLY.j <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model0int.j <- Almer(lnAMR ~ lnBWg * tempTest + lifestage + (1|species) , data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # 
  # Phylo_MMR_W_model1.j <- Almer(lnAMR ~ lnBWg + tempTest + lifestage + (1|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model1.POLY.j <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model1intPOLY.j <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model1int.j <- Almer(lnAMR ~ lnBWg * tempTest + lifestage + (1|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # 
  Phylo_MMR_W_model2.j <- Almer(lnAMR ~ lnBWg + tempTest  + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model2int.j <- Almer(lnAMR ~ lnBWg * tempTest  + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model2.POLY.j <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model2intPOLY.j <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE)  + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # 
  Phylo_MMR_W_model4.j <- Almer(lnAMR ~ lnBWg + tempTest + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model4int.j <- Almer(lnAMR ~ lnBWg * tempTest + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model4.POLY.j <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model4intPOLY.j <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # 
  Phylo_MMR_W_model5.j <- Almer(lnAMR ~ lnBWg + tempTest  + lifestage + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model5int.j <- Almer(lnAMR ~ lnBWg * tempTest + lifestage + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model5.POLY.j <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model5intPOLY.j <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # # 
  # Phylo_MMR_W_model6.j <- Almer(lnAMR ~ lnBWg + tempTest + lifestage + (lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model6int.j <- Almer(lnAMR ~ lnBWg * tempTest+ lifestage + (lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model6.POLY.j <- Almer(lnAMR ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model6intPOLY.j <- Almer(lnAMR ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # 
  
  # MMR warm lifestage X scaling ------------
  Phylo_MMR_W_model2.jint <- Almer(lnAMR ~ lnBWg * lifestage  + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model2int.jint <- Almer(lnAMR ~ lnBWg * lifestage  * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model2.POLY.jint <- Almer(lnAMR ~ lnBWg * lifestage  + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model2intPOLY.jint <- Almer(lnAMR ~ lnBWg * lifestage  * poly(tempTest,2, raw = TRUE)  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # 
  Phylo_MMR_W_model4.jint <- Almer(lnAMR ~ lnBWg * lifestage + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model4int.jint <- Almer(lnAMR ~ lnBWg * lifestage * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model4.POLY.jint <- Almer(lnAMR ~ lnBWg * lifestage + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model4intPOLY.jint <- Almer(lnAMR ~ lnBWg * lifestage * poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # 
  Phylo_MMR_W_model5.jint <- Almer(lnAMR ~ lnBWg * lifestage + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  Phylo_MMR_W_model5int.jint <- Almer(lnAMR ~ lnBWg * lifestage * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model5.POLY.jint <- Almer(lnAMR ~ lnBWg * lifestage + poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # Phylo_MMR_W_model5intPOLY.jint <- Almer(lnAMR ~ lnBWg * lifestage * poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.amr.test, REML=FALSE, A = list(species = A.mmr.w))
  # 
  
  ## BIC -------------
  MMR_W_BIC<-BICdelta(BIC(#Phylo_MMR_W_model0, Phylo_MMR_W_model0int, Phylo_MMR_W_model0.POLY, Phylo_MMR_W_model0intPOLY,
                  # Phylo_MMR_W_model1, Phylo_MMR_W_model1int, Phylo_MMR_W_model1.POLY, Phylo_MMR_W_model1intPOLY,
                  Phylo_MMR_W_model2, Phylo_MMR_W_model2int,# Phylo_MMR_W_model2.POLY, Phylo_MMR_W_model2intPOLY,
                  Phylo_MMR_W_model4, Phylo_MMR_W_model4int,# Phylo_MMR_W_model4.POLY, Phylo_MMR_W_model4intPOLY,
                  Phylo_MMR_W_model5, Phylo_MMR_W_model5int,# Phylo_MMR_W_model5.POLY, Phylo_MMR_W_model5intPOLY,
                  # Phylo_MMR_W_model6, Phylo_MMR_W_model6int, Phylo_MMR_W_model6.POLY, Phylo_MMR_W_model6intPOLY,
                  # juvenile
                  # Phylo_MMR_W_model0.j,  Phylo_MMR_W_model0int.j,  Phylo_MMR_W_model0.POLY.j,  Phylo_MMR_W_model0intPOLY.j, 
                  # Phylo_MMR_W_model1.j,  Phylo_MMR_W_model1int.j,  Phylo_MMR_W_model1.POLY.j,  Phylo_MMR_W_model1intPOLY.j, 
                  Phylo_MMR_W_model2.j,  Phylo_MMR_W_model2int.j, # Phylo_MMR_W_model2.POLY.j,  Phylo_MMR_W_model2intPOLY.j,
                  Phylo_MMR_W_model4.j,  Phylo_MMR_W_model4int.j, # Phylo_MMR_W_model4.POLY.j,  Phylo_MMR_W_model4intPOLY.j, 
                  Phylo_MMR_W_model5.j,  Phylo_MMR_W_model5int.j, # Phylo_MMR_W_model5.POLY.j,  Phylo_MMR_W_model5intPOLY.j, 
                  Phylo_MMR_W_model2.jint,  Phylo_MMR_W_model2int.jint,  #Phylo_MMR_W_model2.POLY.jint,  Phylo_MMR_W_model2intPOLY.jint,
                  Phylo_MMR_W_model4.jint,  Phylo_MMR_W_model4int.jint,  #Phylo_MMR_W_model4.POLY.jint,  Phylo_MMR_W_model4intPOLY.jint, 
                  Phylo_MMR_W_model5.jint,  Phylo_MMR_W_model5int.jint))#,  Phylo_MMR_W_model5.POLY.jint,  Phylo_MMR_W_model5intPOLY.jint))
                  # Phylo_MMR_W_model6.j,  Phylo_MMR_W_model6int.j,  Phylo_MMR_W_model6.POLY.j,  Phylo_MMR_W_model6intPOLY.j))
  
  
  
  # AS / warm temps --------------
  # Phylo_AS_W_model0 <- Almer(lnAS ~ lnBWg + tempTest + (1|species) , data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model0.POLY <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model0intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model0int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) , data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # 
  # Phylo_AS_W_model1 <- Almer(lnAS ~ lnBWg + tempTest + (1|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model1.POLY <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model1intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model1int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # 
  Phylo_AS_W_model2 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model2int <- Almer(lnAS ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model2.POLY <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model2intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE)  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # # 
  Phylo_AS_W_model4 <- Almer(lnAS ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model4int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model4.POLY <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model4intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))

  Phylo_AS_W_model5 <- Almer(lnAS ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model5int <- Almer(lnAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model5.POLY <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model5intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # # 
  # Phylo_AS_W_model6 <- Almer(lnAS ~ lnBWg + tempTest + (lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model6int <- Almer(lnAS ~ lnBWg * tempTest+ (lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model6.POLY <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model6intPOLY <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  
  
  # AS / warm JUVENILE temps --------------
  # Phylo_AS_W_model0.j <- Almer(lnAS ~ lnBWg + tempTest + lifestage + (1|species) , data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model0.POLY.j <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model0intPOLY.j <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model0int.j <- Almer(lnAS ~ lnBWg * tempTest + lifestage + (1|species) , data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # 
  # Phylo_AS_W_model1.j <- Almer(lnAS ~ lnBWg + tempTest + lifestage + (1|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model1.POLY.j <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model1intPOLY.j <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) +lifestage +  (1|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model1int.j <- Almer(lnAS ~ lnBWg * tempTest + lifestage + (1|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # 
  Phylo_AS_W_model2.j <- Almer(lnAS ~ lnBWg + tempTest  + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model2int.j <- Almer(lnAS ~ lnBWg * tempTest  + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model2.POLY.j <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model2intPOLY.j <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE)  + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # # 
  Phylo_AS_W_model4.j <- Almer(lnAS ~ lnBWg + tempTest + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model4int.j <- Almer(lnAS ~ lnBWg * tempTest + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model4.POLY.j <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model4intPOLY.j <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # 
  Phylo_AS_W_model5.j <- Almer(lnAS ~ lnBWg + tempTest  + lifestage + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model5int.j <- Almer(lnAS ~ lnBWg * tempTest + lifestage + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model5.POLY.j <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model5intPOLY.j <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # # 
  # Phylo_AS_W_model6.j <- Almer(lnAS ~ lnBWg + tempTest + lifestage + (lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model6int.j <- Almer(lnAS ~ lnBWg * tempTest + lifestage + (lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model6.POLY.j <- Almer(lnAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model6intPOLY.j <- Almer(lnAS ~ lnBWg * poly(tempTest,2, raw = TRUE) +lifestage + (lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  
  
  # AS lifestageX scaling -----------
  Phylo_AS_W_model2.jint <- Almer(lnAS ~ lnBWg * lifestage  + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model2int.jint <- Almer(lnAS ~ lnBWg * lifestage  * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model2.POLY.jint <- Almer(lnAS ~ lnBWg * lifestage  + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model2intPOLY.jint <- Almer(lnAS ~ lnBWg * lifestage  * poly(tempTest,2, raw = TRUE)  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # 
  Phylo_AS_W_model4.jint <- Almer(lnAS ~ lnBWg * lifestage + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model4int.jint <- Almer(lnAS ~ lnBWg * lifestage * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model4.POLY.jint <- Almer(lnAS ~ lnBWg * lifestage + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model4intPOLY.jint <- Almer(lnAS ~ lnBWg * lifestage * poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))

  Phylo_AS_W_model5.jint <- Almer(lnAS ~ lnBWg * lifestage + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  Phylo_AS_W_model5int.jint <- Almer(lnAS ~ lnBWg * lifestage * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model5.POLY.jint <- Almer(lnAS ~ lnBWg * lifestage + poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # Phylo_AS_W_model5intPOLY.jint <- Almer(lnAS ~ lnBWg * lifestage * poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.as.test, REML=FALSE, A = list(species = A.aas.w))
  # 
  
  ## BIC -------------
  AS_W_BIC<-BICdelta(BIC(#Phylo_AS_W_model0, Phylo_AS_W_model0int, Phylo_AS_W_model0.POLY, Phylo_AS_W_model0intPOLY,
                  # Phylo_AS_W_model1, Phylo_AS_W_model1int, Phylo_AS_W_model1.POLY, Phylo_AS_W_model1intPOLY,
                  Phylo_AS_W_model2, Phylo_AS_W_model2int,# Phylo_AS_W_model2.POLY, Phylo_AS_W_model2intPOLY,
                  Phylo_AS_W_model4, Phylo_AS_W_model4int,# Phylo_AS_W_model4.POLY, Phylo_AS_W_model4intPOLY,
                  Phylo_AS_W_model5, Phylo_AS_W_model5int,# Phylo_AS_W_model5.POLY, Phylo_AS_W_model5intPOLY,
                  # Phylo_AS_W_model6, Phylo_AS_W_model6int, Phylo_AS_W_model6.POLY, Phylo_AS_W_model6intPOLY,
                  # juvenile
                  # Phylo_AS_W_model0.j,  Phylo_AS_W_model0int.j,  Phylo_AS_W_model0.POLY.j,  Phylo_AS_W_model0intPOLY.j, 
                  # Phylo_AS_W_model1.j,  Phylo_AS_W_model1int.j,  Phylo_AS_W_model1.POLY.j,  Phylo_AS_W_model1intPOLY.j, 
                  Phylo_AS_W_model2.j,  Phylo_AS_W_model2int.j, # Phylo_AS_W_model2.POLY.j,  Phylo_AS_W_model2intPOLY.j,
                  Phylo_AS_W_model4.j,  Phylo_AS_W_model4int.j,  # Phylo_AS_W_model4.POLY.j,  Phylo_AS_W_model4intPOLY.j, 
                  Phylo_AS_W_model5.j,  Phylo_AS_W_model5int.j, # Phylo_AS_W_model5.POLY.j,  Phylo_AS_W_model5intPOLY.j, 
                  Phylo_AS_W_model2.jint,  Phylo_AS_W_model2int.jint,#  Phylo_AS_W_model2.POLY.jint,  Phylo_AS_W_model2intPOLY.jint,
                  Phylo_AS_W_model4.jint,  Phylo_AS_W_model4int.jint,#  Phylo_AS_W_model4.POLY.jint,  Phylo_AS_W_model4intPOLY.jint, 
                  Phylo_AS_W_model5.jint,  Phylo_AS_W_model5int.jint))#,  Phylo_AS_W_model5.POLY.jint,  Phylo_AS_W_model5intPOLY.jint))                
                  # Phylo_AS_W_model6.j,  Phylo_AS_W_model6int.j,  Phylo_AS_W_model6.POLY.j,  Phylo_AS_W_model6intPOLY.j))
  
  
  # FAS / warm temps --------------
  # Phylo_FAS_W_model0 <- Almer(lnFAS ~ lnBWg + tempTest + (1|species) , data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model0.POLY <- Almer(lnFAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model0intPOLY <- Almer(lnFAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model0int <- Almer(lnFAS ~ lnBWg * tempTest + (1|species) , data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # 
  # Phylo_FAS_W_model1 <- Almer(lnFAS ~ lnBWg + tempTest + (1|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model1.POLY <- Almer(lnFAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model1intPOLY <- Almer(lnFAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model1int <- Almer(lnFAS ~ lnBWg * tempTest + (1|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  
  Phylo_FAS_W_model2 <- Almer(lnFAS ~ lnBWg + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model2int <- Almer(lnFAS ~ lnBWg * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model2.POLY <- Almer(lnFAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model2intPOLY <- Almer(lnFAS ~ lnBWg * poly(tempTest,2, raw = TRUE)  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # # 
  Phylo_FAS_W_model4 <- Almer(lnFAS ~ lnBWg + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model4int <- Almer(lnFAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model4.POLY <- Almer(lnFAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model4intPOLY <- Almer(lnFAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # 
  Phylo_FAS_W_model5 <- Almer(lnFAS ~ lnBWg + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model5int <- Almer(lnFAS ~ lnBWg * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model5.POLY <- Almer(lnFAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model5intPOLY <- Almer(lnFAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # # 
  # Phylo_FAS_W_model6 <- Almer(lnFAS ~ lnBWg + tempTest + (lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model6int <- Almer(lnFAS ~ lnBWg * tempTest+ (lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model6.POLY <- Almer(lnFAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model6intPOLY <- Almer(lnFAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + (lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # 
  # FAS W juvenile --------------------------------
  # Phylo_FAS_W_model0.j <- Almer(lnFAS ~ lnBWg + tempTest + lifestage + (1|species) , data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model0.POLY.j <- Almer(lnFAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model0intPOLY.j <- Almer(lnFAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model0int.j <- Almer(lnFAS ~ lnBWg * tempTest + lifestage + (1|species) , data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # 
  # Phylo_FAS_W_model1.j <- Almer(lnFAS ~ lnBWg + tempTest + lifestage + (1|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model1.POLY.j <- Almer(lnFAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model1intPOLY.j <- Almer(lnFAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model1int.j <- Almer(lnFAS ~ lnBWg * tempTest + lifestage + (1|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # 
  Phylo_FAS_W_model2.j <- Almer(lnFAS ~ lnBWg + tempTest  + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model2int.j <- Almer(lnFAS ~ lnBWg * tempTest  + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model2.POLY.j <- Almer(lnFAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model2intPOLY.j <- Almer(lnFAS ~ lnBWg * poly(tempTest,2, raw = TRUE)  + lifestage + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # # 
  Phylo_FAS_W_model4.j <- Almer(lnFAS ~ lnBWg + tempTest + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model4int.j <- Almer(lnFAS ~ lnBWg * tempTest + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model4.POLY.j <- Almer(lnFAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model4intPOLY.j <- Almer(lnFAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # 
  Phylo_FAS_W_model5.j <- Almer(lnFAS ~ lnBWg + tempTest  + lifestage + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model5int.j <- Almer(lnFAS ~ lnBWg * tempTest + lifestage + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model5.POLY.j <- Almer(lnFAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model5intPOLY.j <- Almer(lnFAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # # 
  # Phylo_FAS_W_model6.j <- Almer(lnFAS ~ lnBWg + tempTest + lifestage + (lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model6int.j <- Almer(lnFAS ~ lnBWg * tempTest + lifestage + (lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model6.POLY.j <- Almer(lnFAS ~ lnBWg + poly(tempTest,2, raw = TRUE) + lifestage + (lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model6intPOLY.j <- Almer(lnFAS ~ lnBWg * poly(tempTest,2, raw = TRUE) + lifestage + (lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  #
  
  # FAS lifestage X scaling --------------
  Phylo_FAS_W_model2.jint <- Almer(lnFAS ~ lnBWg * lifestage  + tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model2int.jint <- Almer(lnFAS ~ lnBWg * lifestage  * tempTest  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model2.POLY.jint <- Almer(lnFAS ~ lnBWg * lifestage  + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model2intPOLY.jint <- Almer(lnFAS ~ lnBWg * lifestage  * poly(tempTest,2, raw = TRUE)  + (1|species) +(0 + tempTest|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # 
  # 
  Phylo_FAS_W_model4.jint <- Almer(lnFAS ~ lnBWg * lifestage + tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model4int.jint <- Almer(lnFAS ~ lnBWg * lifestage * tempTest + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model4.POLY.jint <- Almer(lnFAS ~ lnBWg * lifestage + poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model4intPOLY.jint <- Almer(lnFAS ~ lnBWg * lifestage * poly(tempTest,2, raw = TRUE) + (1|species) +(0 + lnBWg|species) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # 
  Phylo_FAS_W_model5.jint <- Almer(lnFAS ~ lnBWg * lifestage + tempTest  + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  Phylo_FAS_W_model5int.jint <- Almer(lnFAS ~ lnBWg * lifestage * tempTest + (1|species) +(0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model5.POLY.jint <- Almer(lnFAS ~ lnBWg * lifestage + poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # Phylo_FAS_W_model5intPOLY.jint <- Almer(lnFAS ~ lnBWg * lifestage * poly(tempTest,2, raw = TRUE) + (1|species) + (0 + lnBWg|species:trial) + (1|species:trial), data=data.fas.test, REML=FALSE, A = list(species = A.fas.w))
  # 
  
  ## BIC -------------
  FAS_W_BIC<-BICdelta(BIC(#Phylo_FAS_W_model0, Phylo_FAS_W_model0int, Phylo_FAS_W_model0.POLY, Phylo_FAS_W_model0intPOLY,
                  # Phylo_FAS_W_model1, Phylo_FAS_W_model1int, Phylo_FAS_W_model1.POLY, Phylo_FAS_W_model1intPOLY,
                  Phylo_FAS_W_model2, Phylo_FAS_W_model2int, Phylo_FAS_W_model2.POLY, Phylo_FAS_W_model2intPOLY,
                  Phylo_FAS_W_model4, Phylo_FAS_W_model4int, Phylo_FAS_W_model4.POLY, Phylo_FAS_W_model4intPOLY,
                  Phylo_FAS_W_model5, Phylo_FAS_W_model5int, Phylo_FAS_W_model5.POLY, Phylo_FAS_W_model5intPOLY,
                  # Phylo_FAS_W_model6, Phylo_FAS_W_model6int, Phylo_FAS_W_model6.POLY, Phylo_FAS_W_model6intPOLY,
                  # juvenile
                  # Phylo_FAS_W_model0.j,  Phylo_FAS_W_model0int.j,  Phylo_FAS_W_model0.POLY.j,  Phylo_FAS_W_model0intPOLY.j, 
                  # Phylo_FAS_W_model1.j,  Phylo_FAS_W_model1int.j,  Phylo_FAS_W_model1.POLY.j,  Phylo_FAS_W_model1intPOLY.j, 
                  Phylo_FAS_W_model2.j,  Phylo_FAS_W_model2int.j,  Phylo_FAS_W_model2.POLY.j,  Phylo_FAS_W_model2intPOLY.j,
                  Phylo_FAS_W_model4.j,  Phylo_FAS_W_model4int.j,  Phylo_FAS_W_model4.POLY.j,  Phylo_FAS_W_model4intPOLY.j, 
                  Phylo_FAS_W_model5.j,  Phylo_FAS_W_model5int.j,  Phylo_FAS_W_model5.POLY.j,  Phylo_FAS_W_model5intPOLY.j, 
                  Phylo_FAS_W_model2.jint,  Phylo_FAS_W_model2int.jint,  Phylo_FAS_W_model2.POLY.jint,  Phylo_FAS_W_model2intPOLY.jint,
                  Phylo_FAS_W_model4.jint,  Phylo_FAS_W_model4int.jint,  Phylo_FAS_W_model4.POLY.jint,  Phylo_FAS_W_model4intPOLY.jint, 
                  Phylo_FAS_W_model5.jint,  Phylo_FAS_W_model5int.jint,  Phylo_FAS_W_model5.POLY.jint,  Phylo_FAS_W_model5intPOLY.jint))
                  # Phylo_FAS_W_model6.j,  Phylo_FAS_W_model6int.j,  Phylo_FAS_W_model6.POLY.j,  Phylo_FAS_W_model6intPOLY.j))
  
  
  # generally adding lifestage improves the model, but must compare a smaller dataset/trials becuase not all species have lifes stage assignment. 
  
  ## BIC results ------
  ## Feb 2026
  RMR_BIC # Phylo_RMR_model5 (dBIC 0)
  MMR_BIC # Phylo_MMR_model4intPOLY (dBIC 0)
  AS_BIC  # Phylo_AS_model4intPOLY (dBIC 0)
  FAS_BIC # Phylo_FAS_model2int (dBIC 0)
  
  RMR_W_BIC # Phylo_RMR_W_model4.POLY (dBIC 0)
  MMR_W_BIC # Phylo_MMR_W_model4.POLY (dBIC 0)
  AS_W_BIC  # PPhylo_AS_W_model2.POLY (dBIC 0)
  FAS_W_BIC # Phylo_FAS_W_model2.POLY (dBIC 0)
  
  write.csv(file = here("Data_exports/BICs/RMR_BIC.csv"), x = RMR_BIC, row.names = T)
  write.csv(file = here("Data_exports/BICs/MMR_BIC.csv"), x = MMR_BIC, row.names = T)
  write.csv(file = here("Data_exports/BICs/AS_BIC.csv"), x = AS_BIC, row.names = T)
  write.csv(file = here("Data_exports/BICs/FAS_BIC.csv"), x = FAS_BIC, row.names = T)

  write.csv(file = here("Data_exports/BICs/RMR_W_BIC.csv"), x = RMR_W_BIC, row.names = T)
  write.csv(file = here("Data_exports/BICs/MMR_W_BIC.csv"), x = MMR_W_BIC, row.names = T)
  write.csv(file = here("Data_exports/BICs/AS_W_BIC.csv"), x = AS_W_BIC, row.names = T)
  write.csv(file = here("Data_exports/BICs/FAS_W_BIC.csv"), x = FAS_W_BIC, row.names = T)

  # Best models  --------

  # jan 2026 best 
  # ecol relevant
  # rmr_mod_ER<- Phylo_RMR_model5
  # amr_mod_ER<-Phylo_MMR_model4intPOLY
  # as_mod_ER<-Phylo_AS_model4intPOLY
  # fas_mod_ER<-Phylo_FAS_model2int
  # 
  # # warm
  # rmr_mod_W<-Phylo_RMR_W_model4.POLY
  # amr_mod_W<-Phylo_MMR_W_model4.POLY
  # as_mod_W<-Phylo_AS_W_model2.POLY
  # fas_mod_W<-Phylo_FAS_W_model2.POLY
  
    # jan 2026 best 
  # ecol relevant
  rmr_mod_ER<- Phylo_RMR_model5
  amr_mod_ER<-Phylo_MMR_model4intPOLY
  as_mod_ER<-Phylo_AS_model4intPOLY
  fas_mod_ER<-Phylo_FAS_model2int

  # warm
  rmr_mod_W<-Phylo_RMR_W_model4.POLY
  amr_mod_W<-Phylo_MMR_W_model4.POLY
  as_mod_W<-Phylo_AS_W_model2.POLY
  fas_mod_W<-Phylo_FAS_W_model2.POLY

  
  
  # FIGURES
  # 

# MMR and AMR used interchangeably throughout 
# predict vs real 
# ggplot(data.amrER, aes(predict(amr_mod_ER), lnAMR, color = tempTest))+
#   geom_point()
# ggplot(data.rmrER, aes(predict(rmr_mod_ER), lnRMR, color = tempTest))+
#   geom_point()
# ggplot(data.asER, aes(predict(as_mod_ER), lnAS, color = tempTest))+
#   geom_point()
# ggplot(data.fasER, aes(predict(fas_mod_ER), lnFAS, color = tempTest))+
#   geom_point()
# plot(resid(fas_mod_ER))

# predict vs real 
# ggplot(data.amr.test, aes(predict(amr_mod_W), lnAMR, color = tempTest))+
#   geom_point()
# ggplot(data.rmr.test, aes(predict(rmr_mod_W), lnRMR, color = tempTest))+
#   geom_point()
# ggplot(data.as.test, aes(predict(as_mod_W), lnAS, color = tempTest))+
#   geom_point()
# ggplot(data.fas.test, aes(predict(fas_mod_W), lnFAS, color = tempTest))+
#   geom_point()
# plot(resid(fas_mod_W))

    # 
    # # General scaling plots ------
    # AMRmodel_plot1<-ggplot(data=data.amrER, aes(x=lnBWg, y=lnAMR)) +
    #   geom_point(alpha=0.9,  size=1, pch=1, color="grey75")+
    #   geom_line(data=data.plotAMRint_ER[round(data.plotAMRint_ER$tempTest,2)== 5 |
    #                                       round(data.plotAMRint_ER$tempTest,2)== 15|
    #                                       round(data.plotAMRint_ER$tempTest,2)== 25|
    #                                       round(data.plotAMRint_ER$tempTest,2)== 31,],
    #             aes(y = model_predFE, x=lnBWg,  group=tempTest, color = tempTest),
    #             linewidth=0.5, lty=1,alpha=0.8, show.legend=FALSE) +
    #   geom_point(data=data.amr.test, aes(x=lnBWg, y=lnAMR),
    #              alpha=0.9,  size=1, pch=21, show.legend = FALSE, stroke =0.2,
    #              fill = cols.amr[3], color = cols.amr[1])+
    #   geom_line(data=data.plotAMR_warm[round(data.plotAMR_warm$tempTest,2)==25,],
    #             aes(y = model_predFE, x=lnBWg,  group=tempTest), color=cols.amr[1],
    #             linewidth=1, lty=1, show.legend=FALSE) +
    #   annotate("text",  x = -5.2, y = 11.5,
    #            label = bquote(Optimal:~italic(b)[MMR] == change~with~degree*C),
    #            size=4, hjust=0, family="Helvetica", color = "black")+
    #   annotate("text",  x = -5.2, y = 9.8,
    #            label = bquote(Warm:~italic(b)[MMR] == .(AMR_slope_w)),
    #            size=4, hjust=0, family="Helvetica", color = cols.amr[1])+
    #   annotate("text",  x = 2.5, y = -3,
    #            label = bquote(~5*degree*C:~italic(b)[MMR] == .(AMR_slope5)),
    #            size=3.5, hjust=0, family="Helvetica", color = "black")+
    #   annotate("text",  x = 2.5, y = -4,
    #            label = bquote(15*degree*C:~italic(b)[MMR] == .(AMR_slope15)),
    #            size=3.5, hjust=0, family="Helvetica", color = "black")+
    #   annotate("text",  x = 2.5, y = -5,
    #            label = bquote(25*degree*C:~italic(b)[MMR] == .(AMR_slope25)),
    #            size=3.5, hjust=0, family="Helvetica", color = "black")+
    #   annotate("text",  x = 2.5, y = -6,
    #            label = bquote(31*degree*C:~italic(b)[MMR] == .(AMR_slope31)),
    #            size=3.5, hjust=0, family="Helvetica", color = "black")+
    #   # scale_fill_gradient( low = cols.amr[5], high = cols.amr[1])+
    #   scale_color_gradient( low = "grey", high = "black")+
    #   ylim(x = -6.3, 12)+
    #   xlim(x = -6.3, 12)+
    #   annotate("text", label = paste("n = ", nrow(data.amrER), sep=""),
    #            x = -5.2, y = 8.5, size=3, hjust=0, family="Helvetica", color = "black")+
    #   annotate("text", label = paste("n = ", nrow(data.amr.test), sep=""),
    #            x = -5.2, y = 7.6, size=3, hjust=0, family="Helvetica", color = cols.amr[1])
    # ggformat(AMRmodel_plot1, x_title=expression(italic(ln)*Body~weight~(g)), y_title=expression(italic(ln)*MMR~(mg~O[2]~h^-1)), size_text = 12, print = F)
    # 
    # RMRmodel_plot1<-ggplot(data=data.rmrER, aes(x=lnBWg, y=lnRMR)) +
    #   geom_point(alpha=0.9, size=1, pch=1, color="grey75")+
    #   geom_line(data=data.plotRMR_ER[round(data.plotRMR_ER$tempTest,2)==25,],
    #             mapping=aes(y = model_predFE, x=lnBWg,  group=tempTest, color= tempTest),
    #             color="black", linewidth=0.5, lty=1, alpha=0.8, show.legend=FALSE) +
    #   geom_point(data=data.rmr.test, aes(x=lnBWg, y=lnRMR),
    #              alpha=0.9,  size=1, pch=21, show.legend = FALSE, stroke =0.2,
    #              fill = cols.rmr[3], color = cols.rmr[1])+
    #   geom_line(data=data.plotRMR_warm[round(data.plotRMR_warm$tempTest,2)==25,],
    #             aes(y = model_predFE, x=lnBWg,  group=tempTest),
    #             color=cols.rmr[1], linewidth=1, lty=1, show.legend=FALSE) +
    #   annotate("text",  x = -5.2, y = 11.5, label = bquote(Optimal:~italic(b)[RMR] == .(RMR_slope)),
    #            size=4, hjust=0, family="Helvetica", color = "black")+
    #   annotate("text",  x = -5.2, y = 9.8, label = bquote(Warm:~italic(b)[RMR] == .(RMR_slope_w)),
    #            size=4, hjust=0, family="Helvetica", color = cols.rmr[1])+
    #   scale_fill_gradient( low = cols.rmr[5], high = cols.rmr[1])+
    #   scale_color_gradient( low = cols.rmr[5], high = cols.rmr[1])+
    #   ylim(x = -6.3, 12)+
    #   xlim(x = -6.3, 12)+
    #   annotate("text", label = paste("n = ", nrow(data.rmrER), sep=""),
    #            x = -5.2, y = 8.5, size=3, hjust=0, family="Helvetica", color = "black")+
    #   annotate("text", label = paste("n = ", nrow(data.rmr.test), sep=""),
    #            x = -5.2, y = 7.6, size=3, hjust=0, family="Helvetica", color = cols.rmr[1])
    # ggformat(RMRmodel_plot1, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*RMR~(mg~O[2]~h^-1)),size_text = 12, print = F)
    # 
    # # FAS! 
    # FASmodel_plot1<-ggplot(data=data.fasER, aes(x=lnBWg, y=log(FAS))) +
    #   geom_point(alpha=0.9, size=1, pch=1, color="grey75")+
    #   geom_line(data=data.plotFAS_ER[round(data.plotFAS_ER$tempTest,2)==25,],
    #             aes(y = model_predFE, x=lnBWg,  group=tempTest),
    #             color="black", linewidth=0.5, lty=1, alpha=0.8, show.legend=FALSE) +
    #   geom_point(data=data.fas.test, aes(x=lnBWg, y=log(FAS), fill=tempTest, color = tempTest),
    #              alpha=0.9, size=1, pch=21, show.legend = FALSE, stroke =0.2,
    #              fill = cols.fas[3], color = cols.fas[1])+
    #   geom_line(data=data.plotFAS_warm[round(data.plotFAS_warm$tempTest,2)==25,],
    #             aes(y = model_predFE, x=lnBWg,  group=tempTest), color=cols.fas[1], linewidth=1, lty=1, show.legend=FALSE) +
    #   annotate("text",  x = 2, y = 3.9, label = bquote(Optimal:~italic(b)[FAS] == .(FAS_slope)),
    #            size=4, hjust=0, family="Helvetica", color = "black")+
    #   annotate("text",  x = 2, y = 3.5, label = bquote(Warm:~italic(b)[FAS] == .(FAS_slope_w)),
    #            size=4, hjust=0, family="Helvetica", color = cols.fas[1])+
    #   annotate("text", label = paste("n = ", nrow(data.fasER), sep=""),
    #            x = 4, y = 3.1, size=3, hjust=0, family="Helvetica", color = "black")+
    #   annotate("text", label = paste("n = ", nrow(data.fas.test), sep=""),
    #            x = 4, y = 2.85, size=3, hjust=0, family="Helvetica", color = cols.fas[1])+
    #   scale_fill_gradient(high = cols.fas[1], low = cols.fas[5])+
    #   scale_color_gradient(high = cols.fas[1], low = cols.fas[5])+
    #   ylim(0,4)
    # ggformat(FASmodel_plot1, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*FAS),
    #          size_text = 12, print = FALSE)
    # 
    # # AS 
    # ASmodel_plot1<-ggplot(data=data.asER, aes(x=lnBWg, y=lnAS)) +
    #   geom_point(alpha=0.9, size=1, pch=1, color="grey75")+
    #   geom_line(data=data.plotAS_ER[round(data.plotAS_ER$tempTest,2)==25,], 
    #             mapping=aes(y = model_predFE, x=lnBWg,  group=tempTest, color= tempTest),
    #             color="black", linewidth=0.5, lty=1, alpha=0.8, show.legend=FALSE) +
    #   geom_point(data=data.as.test, aes(x=lnBWg, y=lnAS),
    #              alpha=0.9, size=1, pch=21, show.legend = FALSE, stroke =0.2,
    #              fill = cols.as[3], color = cols.as[1])+
    #   geom_line(data=data.plotAS_warm[round(data.plotAS_warm$tempTest,1)==25,],
    #             aes(y = model_predFE, x=lnBWg,  group=tempTest),
    #             color=cols.as[1], linewidth=1, lty=1, show.legend=FALSE) +
    #   annotate("text",  x = -5.2, y = 11.5, label = bquote(Optimal:~italic(b)[AS] == .(AS_slope)),
    #            size=4, hjust=0, family="Helvetica", color = "black")+
    #   annotate("text",  x = -5.2, y = 9.8, label = bquote(Warm:~italic(b)[AS] == .(AS_slope_w)),
    #            size=4, hjust=0, family="Helvetica", color = cols.as[1])+
    #   scale_fill_gradient( low = cols.as[5], high = cols.as[1])+
    #   scale_color_gradient( low = cols.as[5], high = cols.as[1])+
    #   ylim(x = -6.3, 12)+
    #   xlim(x = -6.3, 12)+
    #   annotate("text", label = paste("n = ", nrow(data.asER), sep=""),
    #            x = -5.2, y = 8.5, size=3, hjust=0, family="Helvetica", color = "black")+
    #   annotate("text", label = paste("n = ", nrow(data.as.test), sep=""),
    #            x = -5.2, y = 7.6, size=3, hjust=0, family="Helvetica", color = cols.as[1])
    # ggformat(ASmodel_plot1, x_title=expression(italic(ln)*Body~mass~(g)),
    #          y_title=expression(italic(ln)*AS~(mg~O[2]~h^-1)), size_text = 12, print = F)
    # 
    # 
    # scaling<-cowplot:::plot_grid(AMRmodel_plot1, RMRmodel_plot1,
    #                              ASmodel_plot1, FASmodel_plot1,
    #                               align = "hv",
    #                               axis = "l",
    #                               nrow = 2,
    #                               ncol = 2,
    #                               labels = "AUTO",
    #                              label_x = c(0.18, 0.18),
    #                              label_y = c(0.895, 0.895),
    #                               label_size = 12)
    # # scaling
    # ggsave(filename = paste("./Figures/Figure3.png", sep=""),
    #        plot=scaling, width = 6.8, height = 6.8, units = "in")
    # 
    # 
    # 
    # # Global outcomes by temprature  --------
    # 
    # AMRmodel_plot1_t<-ggplot(data=data.amrER, aes(x=tempTest, y=lnAMR)) +
    #   geom_point(alpha=0.9, size=1, pch=1, color="grey75")+
    #   geom_line(data=data.plotAMRint_ER[
    #                                 round(data.plotAMRint_ER$lnBWg,2)==-2.54 |
    #                                 round(data.plotAMRint_ER$lnBWg,2)== 2.46 | # ~ 10
    #                                 round(data.plotAMRint_ER$lnBWg,2)== 4.46 | # ~ 100; 119.104 grams
    #                                 round(data.plotAMRint_ER$lnBWg,2)== 6.96
    #                                 ,], # ~ 1000
    #             mapping=aes(y = model_predFE, x=tempTest,  group=lnBWg, color= lnBWg),
    #             # color="black",
    #             linewidth=0.5, lty=1, alpha=0.8, show.legend=FALSE) +
    #   geom_point(data=data.amr.test, aes(x=tempTest, y=lnAMR),
    #              alpha=0.9, size=1, pch=21, show.legend = FALSE, stroke =0.2,
    #              fill = cols.amr[3], color = cols.amr[1])+
    #   geom_line(data=data.plotAMR_warm[round(data.plotAMR_warm$lnBWg,2)== 4.88,],
    #             aes(y = model_predFE, x=tempTest,  group=lnBWg),
    #             color=cols.amr[1], linewidth=1, lty=1, show.legend=FALSE) +
    #   scale_fill_gradient( low = cols.amr[5], high = cols.amr[1])+
    #   scale_color_gradient( low = "grey", high = "black")+
    #   ylim(x = -6.3, 12)+
    #   xlim(x = 0,40)
    # ggformat(AMRmodel_plot1_t, x_title="Temperature ºC",
    #          y_title=expression(italic(ln)*MMR~(mg~O[2]~h^-1)), size_text = 12, print = F)
    # 
    # 
    # # RMR 
    # RMRmodel_plot1_t<-ggplot(data=data.rmrER, aes(x=tempTest, y=lnRMR)) +
    #   geom_point(alpha=0.9, size=1, pch=1, color="grey75")+
    #   geom_line(data=data.plotRMR_ER[round(data.plotRMR_ER$lnBWg,2)== 4.78 ,], # ~ 1000
    #             mapping=aes(y = model_predFE, x=tempTest,  group=lnBWg, color= tempTest),
    #             color="black", linewidth=0.5, lty=1, alpha=0.8, show.legend=FALSE) +
    #   geom_point(data=data.rmr.test, aes(x=tempTest, y=lnRMR),
    #              alpha=0.9, size=1, pch=21, show.legend = FALSE, stroke =0.2,
    #              fill = cols.rmr[3], color = cols.rmr[1])+
    #   geom_line(data=data.plotRMR_warm[round(data.plotRMR_warm$lnBWg,2)== 4.88,],
    #             aes(y = model_predFE, x=tempTest,  group=lnBWg),
    #             color=cols.rmr[1], linewidth=1, lty=1, show.legend=FALSE) +
    #   scale_fill_gradient( low = cols.rmr[5], high = cols.rmr[1])+
    #   scale_color_gradient( low = cols.rmr[5], high = cols.rmr[1])+
    #   ylim(x = -6.3, 12)+
    #   xlim(x = 0,40)
    # ggformat(RMRmodel_plot1_t, x_title="Temperature ºC",
    #          y_title=expression(italic(ln)*RMR~(mg~O[2]~h^-1)), size_text = 12, print = F)
    # 
    # 
    # 
    # # AS 
    # ASmodel_plot1_t<-ggplot(data=data.asER, aes(x=tempTest, y=lnAS)) +
    #   geom_point(alpha=0.9, size=1, pch=1, color="grey75")+
    #   geom_line(data=data.plotAS_ER[round(data.plotAS_ER$lnBWg,2)== 4.78 ,], # ~ 1000
    #             mapping=aes(y = model_predFE, x=tempTest,  group=lnBWg, color= tempTest),
    #             color="black", linewidth=0.5, lty=1, alpha=0.8, show.legend=FALSE) +
    #   geom_point(data=data.as.test, aes(x=tempTest, y=lnAS),
    #              alpha=0.9, size=1, pch=21, show.legend = FALSE, stroke =0.2,
    #              fill = cols.as[3], color = cols.as[1])+
    #   geom_line(data=data.plotAS_warm[round(data.plotAS_warm$lnBWg,2)== 4.88,],
    #             aes(y = model_predFE, x=tempTest,  group=lnBWg),
    #             color=cols.as[1], linewidth=1, lty=1, show.legend=FALSE) +
    #   scale_fill_gradient( low = cols.as[5], high = cols.as[1])+
    #   scale_color_gradient( low = cols.as[5], high = cols.as[1])+
    #   ylim(x = -6.3, 12)+
    #   xlim(x = 0,40)+
    #   geom_point(data=data.as.test[grepl(x = data.as.test$Common_name, pattern = "salmon"),],
    #              aes(x=tempTest, y=lnAS),
    #              alpha=1, size=2, pch=21, show.legend = FALSE, stroke =0.2,
    #              fill = "red", color = "black")  
    # ggformat(ASmodel_plot1_t, x_title="Temperature ºC",
    #          y_title=expression(italic(ln)*AS~(mg~O[2]~h^-1)), size_text = 12, print = F)
    # 
    # 
    # # FAS 
    # FASmodel_plot1_t<-ggplot(data=data.fasER, aes(x=tempTest, y=lnFAS)) +
    #   geom_point(alpha=0.9, size=1, pch=1, color="grey75")+
    #   geom_line(data=data.plotFAS_ER[round(data.plotFAS_ER$lnBWg,2)== 4.78 ,], # ~ 1000
    #             mapping=aes(y = model_predFE, x=tempTest,  group=lnBWg, color= tempTest),
    #             color="black", linewidth=0.5, lty=1, alpha=0.8, show.legend=FALSE) +
    #   geom_point(data=data.fas.test, aes(x=tempTest, y=lnFAS),
    #              alpha=0.9, size=1, pch=21, show.legend = FALSE, stroke =0.2,
    #              fill = cols.fas[3], color = cols.fas[1])+
    #   geom_line(data=data.plotFAS_warm[round(data.plotFAS_warm$lnBWg,2)== 4.88,],
    #             aes(y = model_predFE, x=tempTest,  group=lnBWg),
    #             color=cols.fas[1], linewidth=1, lty=1, show.legend=FALSE) +
    #   scale_fill_gradient( low = cols.fas[5], high = cols.fas[1])+
    #   scale_color_gradient( low = cols.fas[5], high = cols.fas[1])+
    #   ylim(x = -1, 6)+
    #   xlim(x = 0,40)+
    #   geom_point(data=data.fas.test[grepl(x = data.fas.test$Common_name, pattern = "salmon"),],
    #              aes(x=tempTest, y=lnFAS),
    #              alpha=1, size=2, pch=21, show.legend = FALSE, stroke =0.2,
    #              fill = "red", color = "black")  
    # ggformat(FASmodel_plot1_t, x_title="Temperature ºC",
    #          y_title=expression(italic(ln)*FAS~(mg~O[2]~h^-1)), size_text = 12, print = F)
    # 
    # 
    # scaling_t<-cowplot:::plot_grid(AMRmodel_plot1_t, RMRmodel_plot1_t,
    #                              ASmodel_plot1_t, FASmodel_plot1_t,
    #                               align = "hv",
    #                               axis = "l",
    #                               nrow = 2,
    #                               ncol = 2,
    #                               labels = "AUTO",
    #                              label_x = c(0.18, 0.18),
    #                              label_y = c(0.895, 0.895),
    #                               label_size = 12)
    # # scaling
    # ggsave(filename = paste("./Figures/Suppl_Figure_t.png", sep=""),
    #        plot=scaling_t, width = 6.8, height = 6.8, units = "in")
    # 
    # 
    # # Final plots with conceptual figure ---------
    # # All fish together:
    # MRmodel_plot_inset<-
    #   ggplot(data = sum_CItable[sum_CItable$var_repeat == "lnBWg" &  sum_CItable$MR == "RMR" &
    #                               sum_CItable$temp_cat == "ER",])+
    #   geom_linerange(aes(ymin = ci5, ymax = ci95, x = MR))+
    #   geom_point(mapping = aes(x = "MMR", y = AMR_slope31), size = 2, pch = 21,
    #              color = "black", stroke = 0.4,fill = "#5d7d7c")+
    #   geom_point(mapping = aes(x = "MMR", y = AMR_slope5), size = 2, pch = 21,
    #              color = "black", stroke = 0.4,fill = "#c4c7cc")+
    #   geom_point(mapping = aes(x = "MMR", y = AMR_slope15), size = 2, pch = 21,
    #              color = "black", stroke = 0.4,fill = "#a2adb5")+
    #   geom_point(mapping = aes(x = "MMR", y = AMR_slope25), size = 2, pch = 21,
    #              color = "black", stroke = 0.4, fill = "#7e959b")+
    #   geom_point(mapping = aes(x = "RMR", y = RMR_slope), size = 2)+
    #   annotate(geom = "text", x = 1.25, y = 0.7, label = "5", size = 2, angle = 45, hjust =0)+
    #   annotate(geom = "text", x = 1.25, y = 0.75, label = "15", size = 2, angle = 45, hjust =0)+
    #   annotate(geom = "text", x = 1.25, y = 0.80, label = "25", size = 2, angle = 45, hjust =0)+
    #   annotate(geom = "text", x = 1.25, y = 0.84, label = "31ºC", size = 2, angle = 45, hjust =0)+
    #   theme_classic()+
    #   coord_flip()+
    #   ylim(0.5, 1)+
    #   xlab("")+
    #   ylab("Slope value")+
    #   geom_hline(yintercept = c(0.75, 1), linetype = "dashed", linewidth = 0.1)+
    #   theme(axis.text = element_text(size = 8, family = "Helvetica"),
    #         axis.title.x = element_text(size = 9, family = "Helvetica"),
    #       panel.background = element_rect(fill = "transparent",
    #                                       colour = NA_character_), # necessary to avoid drawing panel outline
    #       panel.grid.major = element_blank(), # get rid of major grid
    #       panel.grid.minor = element_blank(), # get rid of minor grid
    #       plot.background = element_rect(fill = "transparent",
    #                                      colour = NA_character_), # necessary to avoid drawing plot outline
    #       legend.background = element_rect(fill = "transparent"),
    #       legend.box.background = element_rect(fill = "transparent"),
    #       legend.key = element_rect(fill = "transparent"))
    # 
    # 
    # MRmodel_plot2<-
    #   ggplot(data=data.rmrER, aes(x=lnBWg, y=lnRMR)) +
    #   geom_line(data=data.plotRMR_ER[round(data.plotRMR_ER$tempTest,2)==25,],
    #             aes(y = model_predFE, x=lnBWg,  group=tempTest, color= tempTest),
    #             color="black", linewidth=0.7, lty=1, show.legend=FALSE) +
    #   geom_line(data=data.plotAMRint_ER[round(data.plotAMRint_ER$tempTest,2)==5 | 
    #                                         round(data.plotAMRint_ER$tempTest,2)==15 |
    #                                         round(data.plotAMRint_ER$tempTest,2)==25 |
    #                                         round(data.plotAMRint_ER$tempTest,2)==31,],
    #             aes(y = model_predFE, x=lnBWg,  group=factor(tempTest), color= factor(tempTest)),
    #             size=0.7, lty=1, show.legend=FALSE) +
    #   scale_color_manual(values = c("#c4c7cc", "#a2adb5", "#7e959b","#5d7d7c"))+
    #   ylim(x = -6.3, 12)+
    #   xlim(x = -6.3, 12)+
    #   annotate("text",  x = -5.8, y = 11.8, label = "B", fontface =2,
    #            size=4.5, hjust=1, family="Helvetica", color = "black")+
    #   # annotate("text",  x = 9, y = 8.4, label = expression(paste("MMR")),
    #   #          size=5, hjust=0, family="Helvetica", color = "black", angle = 42)+
    #   # annotate("text",  x = 9, y = 4.1, label = expression(paste("RMR")),
    #   #          size=5, hjust=0, family="Helvetica", color = "black", angle = 42)+
    #   annotate("text",  x = -2.4, y = -5.75, label = expression(paste("OPTIMAL TEMPERATURES")),
    #            size=3.5, hjust=0, family="Helvetica", color = "black", angle = 0)+
    #   xlab(expression(italic(ln)*Body~mass~(g)))+
    #   ylab(expression(italic(ln)*MR~(mg~O[2]~h^-1)))+
    #   theme(axis.text.y=element_text(size=12, colour= 'black'),
    # 		axis.text.x=element_text(size=12, colour= 'black'),
    # 		axis.line.y=element_line(colour = 'black',linewidth=0.5),
    # 		axis.line.x=element_line(colour = 'black',linewidth=0.5),
    # 		axis.ticks.y=element_line(size=0.5),
    # 		panel.background = element_blank(),
    # 		axis.ticks.x.bottom = element_line(linewidth=0.5, colour = "black"),
    # 	  axis.title.y=element_text(size=12),
    # 		axis.title.x=element_text(size=12),
    # 		panel.border = element_rect(linetype = "solid",fill=NA, colour = "black"))+
    #   inset_element(MRmodel_plot_inset, 0.01, 0.6, 0.7, 1)
    # MRmodel_plot2
    #   
    # # warm final plots 
    # # # All fish together:
    # MRmodel_plot_inset_w<-
    #   ggplot(data = sum_CItable[sum_CItable$var_repeat == "lnBWg" &
    #                               c(sum_CItable$MR == "RMR" | sum_CItable$MR == "MMR") &
    #                               sum_CItable$temp_cat == "W",])+
    #   geom_hline(yintercept = c(0.75, 1), linetype = "dashed", linewidth = 0.1)+
    #   geom_linerange(aes(ymin = ci5, ymax = ci95, x = MR, color = MR), show.legend = F)+
    #   geom_point(mapping = aes(x = "MMR", y = AMR_slope_w), size = 2, pch=21, 
    #              color = "black", fill = cols.amr[3], stroke = 0.4)+
    #   geom_point(mapping = aes(x = "RMR", y = RMR_slope_w), size = 2,pch =21,
    #              color = "black", fill = cols.rmr[3], stroke = 0.4)+
    #   theme_classic()+
    #   coord_flip()+
    #   scale_color_manual(values = c(cols.amr[2], cols.rmr[2]))+
    #   ylim(0.5, 1)+
    #   xlab("")+
    #   ylab("Slope value")+
    #   theme(axis.text = element_text(size = 8, family = "Helvetica"),
    #         axis.title.x = element_text(size = 9, family = "Helvetica"),
    #         panel.background = element_rect(fill = "transparent",
    #                                         colour = NA_character_), # necessary to avoid drawing panel outline
    #         panel.grid.major = element_blank(), # get rid of major grid
    #         panel.grid.minor = element_blank(), # get rid of minor grid
    #         plot.background = element_rect(fill = "transparent",
    #                                        colour = NA_character_), # necessary to avoid drawing plot outline
    #         legend.background = element_rect(fill = "transparent"),
    #         legend.box.background = element_rect(fill = "transparent"),
    #         legend.key = element_rect(fill = "transparent"))
    # 
    # MRmodel_plot2_w<-
    #   ggplot(data=data.rmr.test, aes(x=lnBWg, y=lnRMR)) +
    #   geom_line(data=data.plotRMR_warm[round(data.plotRMR_warm$tempTest,2)==25,],
    #             aes(y = model_predFE, x=lnBWg,  group=tempTest, color= tempTest),
    #             linewidth=0.7, lty=1, show.legend=FALSE, color = cols.rmr[2]) +
    #   geom_line(data=data.plotAMR_warm[round(data.plotAMR_warm$tempTest,2)==25,],
    #             aes(y = model_predFE, x=lnBWg,  group=tempTest),
    #             linewidth=0.7, lty=1, show.legend=FALSE, color = cols.amr[2]) +
    #   scale_color_gradient( low = "grey", high = "black")+
    #   ylim(x = -6.3, 12)+
    #   xlim(x = -6.3, 12)+
    #   annotate("text",  x = -5.8, y = 11.8, label = "D", fontface =2,
    #            size=4.5, hjust=1, family="Helvetica", color = "black")+
    #   # annotate("text",  x = 6.2, y = 6, label = expression(paste("MMR")),
    #   #          size=5, hjust=0, family="Helvetica", angle = 30, color = cols.amr[1])+
    #   # annotate("text",  x = 6.7, y = 3.1, label = expression(paste("RMR")),
    #   #          size=5, hjust=0, family="Helvetica", angle = 42, color = cols.rmr[1])+
    #   annotate("text", x = -2.4, y = -5.75, label = expression(paste("WARM TEMPERATURES")),
    #            size=3.5, hjust=0, family="Helvetica", color = "black", angle = 0)+
    #   xlab(expression(italic(ln)*Body~mass~(g)))+
    #   ylab(expression(italic(ln)*MR~(mg~O[2]~h^-1)))+
    #   theme(axis.text.y=element_text(size=12, colour= 'black'),
    # 		axis.text.x=element_text(size=12, colour= 'black'),
    # 		axis.line.y=element_line(colour = 'black',linewidth=0.5),
    # 		axis.line.x=element_line(colour = 'black',linewidth=0.5),
    # 		axis.ticks.y=element_line(size=0.5),
    # 		panel.background = element_blank(),
    # 		axis.ticks.x.bottom = element_line(linewidth=0.5, colour = "black"),
    # 	  axis.title.y=element_text(size=12),
    # 		axis.title.x=element_text(size=12),
    # 		panel.border = element_rect(linetype = "solid",fill=NA, colour = "black"))+
    #   inset_element(MRmodel_plot_inset_w, 0.01, 0.6, 0.7, 1)
    # # MRmodel_plot2_w
    # 
    # # conceptual figures 
    # data_sim<-read.csv(here("./Data/Conceptual_fig_data.csv"))
    # 
    # psim<-ggplot(data_sim[data_sim$simID==3,], aes(x=log(BW_kg), y=log(pred.mmr.mgO2min) ))+
    #   geom_abline(slope = data_sim[data_sim$simID==3,"slope.MMR"][1],
    #               intercept = data_sim[data_sim$simID==3,"int.MMR"][1],
    #               color = "#7e959b", linewidth=0.7)+ # MMR
    #   geom_abline(slope = data_sim[data_sim$simID==3,"slope.SMR"][1],
    #               intercept = data_sim[data_sim$simID==3,"int.SMR"][1],
    #               color = "black", linewidth=0.7)+ # MMR
    #   lims(x = c(-4.5, 4.5), y = c(-4.5, 4.5))+
    #   annotate("text",  x = -4.3, y = 4.4, label = "A", fontface =2,
    #            size=4.5, hjust=1, family="Helvetica", color = "black")+
    #   annotate("text",  x = 0.8, y = 3.5, label = expression(paste("MMR")),
    #        size=5, hjust=0, family="Helvetica", angle = 42, color = "#7e959b")+
    #   annotate("text",  x = 1, y = 0.9, label = expression(paste("RMR")),
    #            size=5, hjust=0, family="Helvetica", angle = 33, color = "black")+
    #   annotate("text",  x = -3.5, y = -4.2, label = expression(paste("OPTIMAL TEMPERATURES")),
    #            size=3.5, hjust=0, family="Helvetica", color = "black", angle = 0)+
    #   annotate("text",  x = -3.8, y = 4.3, label = expression(paste("HYPOTHESIS")),
    #            size=4, hjust=0, family="Helvetica", color = "black", angle = 0)+
    #   annotate("text",  x = -3.8, y = 3.6, label = expression(italic(b)[MMR]~">"~italic(b)[RMR]),
    #            size=3.5, hjust=0, family="Helvetica", color = "black", angle = 0)+
    #   annotate("text",  x = 1.75, y = -0.8,
    #            label = expression(italic(b)[MMR]~"="~1.00),
    #            size=3.5, hjust=0, family="Helvetica", color = "black", angle = 0)+
    #   annotate("text",  x = 1.75, y = -1.5,
    #            label = expression(italic(b)[RMR]~"="~0.89),
    #            size=3.5, hjust=0, family="Helvetica", color = "black", angle = 0)+
    #   annotate("text",  x = 1.75, y = -2.2,
    #            label = expression(italic(b)[AS]~"="~1.05),
    #            size=3,5, hjust=0, family="Helvetica", color = "black", angle = 0)+
    #   annotate("text",  x = 1.75, y = -2.9,
    #            label = expression(italic(b)[FAS]~"="~0.11),
    #            size=3.5, hjust=0, family="Helvetica", color = "black", angle = 0)+
    #   xlab(expression(italic(ln)*Body~mass~(g)))+
    #   ylab(expression(italic(ln)*MR~(mg~O[2]~h^-1)))+
    #   theme(axis.text.y=element_text(size=12, colour= 'black'),
    # 		axis.text.x=element_text(size=12, colour= 'black'),
    # 		axis.line.y=element_line(colour = 'black',size=0.5),
    # 		axis.line.x=element_line(colour = 'black',size=0.5),
    # 		axis.ticks.y=element_line(size=0.5),
    # 		panel.background = element_blank(),
    # 		axis.ticks.x.bottom = element_line(size=0.5, colour = "black"),
    # 	  axis.title.y=element_text(size=12),
    # 		axis.title.x=element_text(size=12),
    # 		panel.border = element_rect(linetype = "solid",fill=NA, colour = "black"))
    # # psim
    #     
    # # conceptual WArm 
    # psim_w<-ggplot(data_sim[data_sim$simID==1,], aes(x=log(BW_kg), y=log(pred.mmr.mgO2min) ))+
    #   geom_abline(slope = data_sim[data_sim$simID==1,"slope.MMR"][1],
    #               intercept = data_sim[data_sim$simID==1,"int.MMR"][1],
    #               color = cols.amr[2], linewidth=0.7)+ # MMR
    #   geom_abline(slope = data_sim[data_sim$simID==1,"slope.SMR"][1],
    #               intercept = data_sim[data_sim$simID==1,"int.SMR"][1],
    #               color = cols.rmr[2], linewidth=0.7)+ # MMR
    #   lims(x = c(-4.5, 4.5), y = c(-4.5, 4.5))+
    #   annotate("text",  x = -4.3, y = 4.4, label = "C", fontface =2,
    #            size=4.5, hjust=1, family="Helvetica", color = "black")+
    #   annotate("text",  x = 0.8, y = 3.3, label = expression(paste("MMR")),
    #        size=5, hjust=0, family="Helvetica", angle = 32, color = cols.amr[2])+
    #   annotate("text",  x = 1, y = 0.9, label = expression(paste("RMR")),
    #            size=5, hjust=0, family="Helvetica", angle = 47, color = cols.rmr[2])+
    #   annotate("text",  x = -3.3, y = -4.2, label = expression(paste("WARM TEMPERATURES")),
    #            size=3.5, hjust=0, family="Helvetica", color = "black", angle = 0)+
    #   annotate("text",  x = -3.8, y = 4.3, label = expression(paste("HYPOTHESIS")),
    #            size=4, hjust=0, family="Helvetica", color = "black", angle = 0)+
    #   annotate("text",  x = -3.8, y = 3.6, label = expression(italic(b)[MMR]~"<"~italic(b)[RMR]),
    #            size=4, hjust=0, family="Helvetica", color = "black", angle = 0)+
    #   annotate("text",  x = 1.75, y = -0.8,
    #            label = expression(italic(b)[MMR]~"="~0.75),
    #            size=3.5, hjust=0, family="Helvetica", color = "black", angle = 0)+
    #   annotate("text",  x = 1.75, y = -1.5,
    #            label = expression(italic(b)[RMR]~"="~0.95),
    #            size=3.5, hjust=0, family="Helvetica", color = "black", angle = 0)+
    #   annotate("text",  x = 1.75, y = -2.2,
    #            label = expression(italic(b)[AS]~"="~0.69),
    #            size=3.5, hjust=0, family="Helvetica", color = "black", angle = 0)+
    #   annotate("text",  x = 1.75, y = -2.9,
    #            label = expression(italic(b)[FAS]~"="~-0.14),
    #            size=3.5, hjust=0, family="Helvetica", color = "black", angle = 0)+
    #   xlab(expression(italic(ln)*Body~mass~(g)))+
    #   ylab(expression(italic(ln)*MR~(mg~O[2]~h^-1)))+
    #   theme(axis.text.y=element_text(size=12, colour= 'black'),
    # 		axis.text.x=element_text(size=12, colour= 'black'),
    # 		axis.line.y=element_line(colour = 'black',size=0.5),
    # 		axis.line.x=element_line(colour = 'black',size=0.5),
    # 		axis.ticks.y=element_line(size=0.5),
    # 		panel.background = element_blank(),
    # 		axis.ticks.x.bottom = element_line(size=0.5, colour = "black"),
    # 	  axis.title.y=element_text(size=12),
    # 		axis.title.x=element_text(size=12),
    # 		panel.border = element_rect(linetype = "solid",fill=NA, colour = "black"))
    # # psim_w
    # # patchwork library 
    # top<-  psim + MRmodel_plot2+ plot_layout(axis_titles = "collect") 
    # # top
    # bottom <- psim_w + MRmodel_plot2_w +  plot_layout(axis_titles = "collect")
    #   
    # # plot_layout(widths = c(2, 1), heights = unit(c(5, 1), c('cm', 'null')))
    # 
    # # save the plots 
    # mainplot<-
    #   cowplot::plot_grid(top, bottom,
    #                    align = "hv", nrow = 2)
    # 
    # ggsave(filename = paste("./Figures/FigureMAIN.png", sep=""),
    #        plot=mainplot, width = 6.8, height = 6.8, units = "in")
    # 
    # # --- misc NOT USED -----
    # # fas_temp<-ggplot(data=data.fas, aes(y=FAS, x = tempTest))+
    # #   geom_point(show.legend = F)+
    # #   ylim(0,20)+
    # #   facet_grid(DemersPelag~test_category3)+
    # #   geom_smooth(method = "lm")
    # # ggformat(fas_temp, y_title = expression(FAS~(MMR/RMR)), x_title = "Temperature ºC", print = F)
    # # fas_temp<-fas_temp+theme(legend.position = "none")
    # # fas_temp
    # 
    # # Not used in manuscript -------------
    # # Overall all fish together, with booted CI:
    # # MRmodel_plot1<-ggplot() +
    # #   geom_ribbon(data=data.plotRMR_ER, mapping = aes(y = model_predFE, x=lnBWg, ymin=CI_2.5, ymax=CI_97.5, group = tempTestK1000), linetype=2, alpha=0.1, fill = "grey30")+
    # #   geom_ribbon(data.plotAMRint_ER, mapping = aes(y = model_predFE, x=lnBWg, ymin=CI_2.5, ymax=CI_97.5, group = tempTestK1000), linetype=2, alpha=0.1, fill = "grey30")+
    # #   geom_line(data=data.plotAMRint_ER[round(data.plotAMRint_ER$tempTestK1000,2)==3.49,], mapping = aes(y = model_predFE, x=lnBWg,  group=tempTestK1000), color="black", linewidth=0.7, lty=1, show.legend=FALSE) +
    # #   geom_line(data=data.plotRMR_ER[round(data.plotRMR_ER$tempTestK1000,2)==3.49,], aes(y = model_predFE, x=lnBWg,  group=tempTestK1000, color= tempTestK1000), color="black", linewidth=0.7, lty=1, show.legend=FALSE) +
    # #   scale_y_continuous(limits = c(-7, 12), breaks = seq(-7,12, 2))+
    # #   scale_x_continuous(limits = c(-7, 12), breaks = seq(-7,12, 2))
    # # ggformat(MRmodel_plot1, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*MR~(mg~O[2]~h^-1)), print = F)
    # # 
    # # MRmodel_plotW<-ggplot() +
    # #   geom_ribbon(data.plotAMR_warm, mapping = aes(y = model_predFE, x=lnBWg, ymin=CI_2.5, ymax=CI_97.5, group = tempTestK1000), linetype=2, alpha=0.1, fill = cols.amr[2])+
    # #   geom_ribbon(data=data.plotRMR_warm, mapping = aes(y = model_predFE, x=lnBWg, ymin=CI_2.5, ymax=CI_97.5, group = tempTestK1000), linetype=2, alpha=0.1, fill = cols.rmr[2])+
    # #   geom_line(data=data.plotAMR_warm[round(data.plotAMR_warm$tempTestK1000,2)==3.49,], mapping = aes(y = model_predFE, x=lnBWg,  group=tempTestK1000), color=cols.amr[1], size=0.7, lty=1, show.legend=FALSE) +
    # #   geom_line(data=data.plotRMR_warm[round(data.plotRMR_warm$tempTestK1000,2)==3.49,], aes(y = model_predFE, x=lnBWg,  group=tempTestK1000, color= tempTestK1000), color=cols.rmr[1], size=0.7, lty=1, show.legend=FALSE) +
    # #   scale_y_continuous(limits = c(-7, 12), breaks = seq(-7,12, 2))+
    # #   scale_x_continuous(limits = c(-7, 12), breaks = seq(-7,12, 2))
    # # ggformat(MRmodel_plotW, x_title=expression(italic(ln)*Body~mass~(g)), y_title=expression(italic(ln)*MR~(mg~O[2]~h^-1)), print = F)
