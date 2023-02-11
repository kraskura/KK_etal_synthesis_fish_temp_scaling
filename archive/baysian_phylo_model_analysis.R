
# library("MCMCglmm") # implements Markov chain Monte Carlo routines for fitting multi-response generalized linear mixed models
# Loading required package: Matrix
# Loading required package: coda
# Loading required package: ape
# https://cran.r-project.org/web/packages/MCMCglmm/vignettes/Overview.pdf
# source the data --------
source("/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Codes/Metadata-temperature/JEB/get_scaling_data_temp.R")
source("/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Codes/MR-fish-metadata/deltaIC_BICdelta_ggformat.R")

data.list<-get_scaling_data_temp(data.amr = "/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Data/MR-fish-metadata-data/Fish_AMR_temp_dataset_jun2021.csv",
                                 data.rmr = "/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Data/MR-fish-metadata-data/Fish_RMR_temp_dataset_jun2021.csv",
                                 ecology.data = "/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Data/MR-fish-metadata-data/Kraskura_JEB_species_ecology.csv", 
                                 onlyTop.above = TRUE)


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





# Get phylo tree  #######################
# 1. rotl & ape -----------

# ************************************
# https://taylorreiter.github.io/2017-07-28-Taxonomy-from-Species-Name-in-R/
# http://ape-package.ird.fr/
# https://cran.r-project.org/web/packages/ape/vignettes/DrawingPhylogenies.pdf
library(rotl)
library(ape)
# ************************************

taxon_search <- tnrs_match_names(names=unique(levels(data.rmrER$species)), context_name="Vertebrates") ## hhave matches
knitr::kable(taxon_search)
tnrs_contexts()

# have a match for all species

ott_in_tree <- ott_id(taxon_search)[is_in_tree(ott_id(taxon_search))]
tr <- tol_induced_subtree(ott_ids = ott_in_tree)
class(tr)

tr$tip.label
tr$tip.label <- strip_ott_ids(tr$tip.label, remove_underscores = TRUE)

labels <- as.data.frame(tr$tip.label)
labels$`tr$tip.label`[which(labels$`tr$tip.label` == "Oncorhynchus mykiss (species in domain Eukaryota)")]<-"Oncorhynchus mykiss"
labels$`tr$tip.label`[which(labels$`tr$tip.label` == "Gadus morhua (species in domain Eukaryota)")]<-"Gadus morhua"
labels$`tr$tip.label`[which(labels$`tr$tip.label` == "Leiocassis longirostris")]<-"Tachysurus dumerili" # Leiocassis longirostris by OTL and FIshbase has Tachysurus dumerili
labels$`tr$tip.label`[which(labels$`tr$tip.label` == "Tachysurus vachellii")]<-"Pseudobagrus vachellii" # Fishbase = "Pseudobagrus vachellii" OTL = "Tachysurus_vachellii"
labels$`tr$tip.label`[which(labels$`tr$tip.label` == "Rhinogobius similis")]<- "Rhinogobius giurinus" # fishbase = "Rhinogobius giurinus" OTL= Rhinogobius similis
labels$`tr$tip.label`[which(labels$`tr$tip.label` == "Chromis atripectoralis")]<- "Chromis atripectoralis" # fishbase = "Rhinogobius giurinus" OTL= Rhinogobius similis
tr$tip.label <- labels$`tr$tip.label`

tr$tip.label 
tr$tip.label  %in% data.rmrER$species # all tip labels are found in dataset 
unique(data.rmrER$species) %in% tr$tip.label # all dataset species are found in the phylo tree

plot(tr)

# >>> look for a tree already exists and that would work for us
# studies_properties()
# studies_find_trees(property = "ot:ottId", value = as.character(ott_id(taxon_search)[1]))
# ott_in_tree <- ott_id(taxon_search)[is_in_tree(ott_id(taxon_search))]
# tr2 <- tol_induced_subtree(ott_ids = ott_in_tree)
# plot(tr2)

# >>> compute branch lengths
tr2<-compute.brlen(tr)
plot(tr2)
tr2$edge.length

A <- ape::vcv.phylo(tr2)
tree <- compute.brlen(tr2)
cor <- vcv(tree, cor = T)
data.rmrER$phylo<-data.rmrER$species


# ************************************
# Figures
plot(tree, cex=0.3, label.offset = 0.05)






# 2. fishphylogeny & fishtree -----------

# ************************************# ************************************
# https://cran.r-project.org/web/packages/fishtree/vignettes/community-analysis.html
# https://cran.r-project.org/web/packages/fishtree/vignettes/comparative-analysis.html
library(fishtree)
loadNamespace("rfishbase")
# ************************************
species.fishtree<-as.character(unique(data.rmrER$species))
fish.tr <- fishtree_phylogeny(species = species.fishtree)
fish.tr$tip.label
fish.tr$tip.label <- strip_ott_ids(fish.tr$tip.label, remove_underscores = TRUE)

labels <- as.data.frame(fish.tr$tip.label)
fish.tr$tip.label  %in% data.rmrER$species
fish.not.found<-which(unique(data.rmrER$species) %in% fish.tr$tip.label  == FALSE)
unique(data.rmrER$species) [fish.not.found]
# [1] Bellapiscis medius     
# [2] Caesio cuning          
# [3] Channa argus           
# [4] Enteromius neumayeri   
# [5] Hemiscyllium ocellatum 
# [6] Neogobius melanostomus 
# [7] Salvelinus alpinus     
# [8] Sinibrama taeniatus    
# [9] Somniosus microcephalus
# [10] Tachysurus dumerili    
# [11] Torpedo marmorata 









# Scaling models ------
# 1. brms --------
# https://paul-buerkner.github.io/brms/articles/brms_phylogenetics.html
# following the section: "A Phylogenetic Model with Repeated Measurements"
library(brms)

# model 1. Simple: NO temperature
model_repeat <- brm(
  lnRMR ~ lnBWg + (1|gr(phylo, cov = A)) + (1|species), 
  data = data.rmrER, 
  family = gaussian(), 
  data2 = list(A = A),
  prior = c(
    prior(normal(0, 10), "b"), # the first number is estimate, the second is SD 
    prior(normal(3, 2.5), "Intercept"),
    prior(student_t(3, 0, 20), "sd"), 
    prior(student_t(3, 0, 20), "sigma")),
  sample_prior = TRUE,
  chains = 2, cores = 2, 
  iter = 10000, warmup = 1000,
  control = list(adapt_delta = 0.99) # << for convergence 
)


# Warning messages:
# 1: There were 1 divergent transitions after warmup. See
# http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
# to find out why this is a problem and how to eliminate them. 
# 2: There were 9439 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
# http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded 
# 3: Examine the pairs() plot to diagnose sampling problems

pairs(model_repeat)
# RESULTS:
# 1. 1 divergent transitions after warmup
summary(model_repeat)
plot(model_repeat)
plot(model_repeat, N = 2, ask = FALSE)
plot(conditional_effects(model_repeat), points = TRUE) 

hyp <- paste(
  "sd_phylo__Intercept^2 /", 
  "(sd_phylo__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0"
)
(hyp <- hypothesis(model_repeat, hyp, class = NULL))

summary(model_repeat)

plot(hyp)
plot(model_repeat, N = 2, ask = FALSE)
plot(conditional_effects(model_repeat), points = TRUE) 


# 2. metfor --------
# https://cran.r-project.org/web/packages/metafor/vignettes/metafor.pdf
# http://environmentalcomputing.net/meta-analysis-3/
# https://www.metafor-project.org/doku.php/diagram
# https://wviechtb.github.io/metafor/reference/rma.mv.html#
library(metafor)
# library(Hmisc)

# !!! Drawback - cannot use all raw data, only the study means? 


# 3. pgls --------
library(pgls)

 

# 4. MCMCglmm -------
















### ----------- BELOW IS CJFAS dataset, play # 1------------

# Import data ---------
source("/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Codes/MR-fish-metadata/deltaIC_BICdelta_ggformat.R")
source("/Users/kristakraskura/Github_repositories/Plots-formatting/ggplot_format.R")
source("/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Codes/MR-fish-metadata/get.scaling.data.R")

data.list<-get.scaling.data(data.amr = "Fish_AMR_temp_dataset_mar2021.csv",
                            data.rmr = "Fish_RMR_temp_dataset_mar2021.csv",
                            wd = "/Users/kristakraskura/Github_repositories/Metabolic-scaling-fish/Data/MR-fish-metadata-data/",
                            # wd.amr.Species.data = "/Users/kristakraskura/Desktop/BOX/UCSB/Research/Metabolic_scaling/ms-AMR-RMR 4-message paper /data analysis/SpeciesIC/AMR/IC Data outputs/",
                            # wd.rmr.Species.data = "/Users/kristakraskura/Desktop/BOX/UCSB/Research/Metabolic_scaling/ms-AMR-RMR 4-message paper /data analysis/SpeciesIC/RMR/IC Data outputs/",
                            colors = c("#0083FF","#f45905","#0083FF","#f45905", "#A9B8F0","#FFC6A5", "#00B498", "#007E66"))

data.amr <- data.frame(data.list[1])
data.rmr <- data.frame(data.list[2])
dataMR <- data.frame(data.list[3])
data.as <- data.frame(data.list[4])

# Best models: from Kraskura_AMR_RMR_combined_y2020.R-----------
#Random species and slopes model by species with free parameter
AMRweight.model17=lmer(lnAMR~lnBWg  + tempTest + (0+lnBWg|species) + (1|species) + (1|species:trial), data=data.as, REML=FALSE)
summary(AMRweight.model17)
logLik(AMRweight.model17)

#Random species and slopes model by species with free parameter
RMRweight.model17=lmer(lnRMR~lnBWg  + tempTest + (0+lnBWg|species) + (1|species) + (1|species:trial), data=data.rmr, REML=FALSE)
summary(RMRweight.model17)
logLik(RMRweight.model17)

# # best AS model that is used for size corrections in AS plots
# ASweight.model17=lmer(lnAS~lnBWg  + tempTest + (0+lnBWg|species) + (1|species) + (1|species:trial), data=data.amr, REML=FALSE)
# summary(ASweight.model17)
# logLik(ASweight.model17)


# http://devillemereuil.legtux.org/wp-content/uploads/2018/08/opm-11.pdf


# phylogeny ----------
#add tree
taxon_search <- tnrs_match_names(names=unique(levels(data.rmr$species)), context_name="Vertebrates") ## hhave matches
knitr::kable(taxon_search)
tnrs_contexts()
# 
# tree <- rtree(n = length(unique(data.rmr$species)))
# tree_l <- compute.brlen(tree, method = 'Grafen',power = 1)
# A <- ape::vcv.phylo(tree_l)


# data$ott_name <- unique_name(taxon_search)
# data$ott_id <- taxon_search$ott_id

ott_in_tree <- ott_id(taxon_search)[is_in_tree(ott_id(taxon_search))]
tr <- tol_induced_subtree(ott_ids = ott_in_tree)
class(tr)

tr$tip.label
tr$tip.label <- strip_ott_ids(tr$tip.label, remove_underscores = TRUE)

labels <- as.data.frame(tr$tip.label)
labels$`tr$tip.label`[32]<-"Oncorhynchus mykiss"
labels$`tr$tip.label`[63]<-"Gadus morhua"
labels$`tr$tip.label`[28]<-"Tachysurus dumerili" # Leiocassis longirostris by OTL and FIshbase has Tachysurus dumerili
labels$`tr$tip.label`[30]<-"Pseudobagrus vachellii" # Fishbase = "Pseudobagrus vachellii" OTL = "Tachysurus_vachellii"
labels$`tr$tip.label`[62]<- "Rhinogobius giurinus" # fishbase = "Rhinogobius giurinus" OTL= Rhinogobius similis
tr$tip.label <- labels$`tr$tip.label`

tr$tip.label %in% data.rmr$species
plot(tr)

# 2. look for a tree that works
studies_properties()
studies_find_trees(property = "ot:ottId", value = as.character(ott_id(taxon_search)[1]))
ott_in_tree <- ott_id(taxon_search)[is_in_tree(ott_id(taxon_search))]
tr2 <- tol_induced_subtree(ott_ids = ott_in_tree)
plot(tr2)


# 1.1. compute branch lengths
tr2<-compute.brlen(tr)
plot(tr2)
tr2$edge.length

A <- ape::vcv.phylo(tr2)
data.rmr$phylo<-data.rmr$species


# 'brms' package
# https://cran.r-project.org/web/packages/brms/vignettes/brms_phylogenetics.html

model_repeat <- brm(
  lnRMR ~ lnBWg + (1|gr(phylo, cov = A)) + (1|species), 
  data = data.rmr, 
  family = gaussian(), 
  data2 = list(A = A),
  prior = c(
    prior(normal(1, 1), "b"),
    prior(normal(3, 2.5), "Intercept"),
    prior(student_t(3, 0, 2.5), "sd"),
    prior(student_t(3, 0, 2.5), "sigma")),
    # prior(student_t(3, 0, 2.5), class = "sd", group = "species"),
    # prior(student_t(3, 0, 2.5), class = "sd", coef="Intercept", group = "species"),
    # prior(student_t(3, 0, 2.5), class = "sd", group = "phylo"),
    # prior(student_t(3, 0, 2.5), class = "sd", coef="Intercept", group = "phylo")),
  sample_prior = TRUE,
  chains = 2, cores = 2, 
  iter = 10000, warmup = 1000,
  control = list(adapt_delta = 0.99) # << for convergence 
)



# get_prior(lnRMR ~ lnBWg + (1|gr(phylo, cov = A)) + (1|species), 
#           data = data.rmr, 
#           family = gaussian(), 
#           data2=list(A=A))


# convergence issues: 
# Warning messages:
#   1: There were XX divergent transitions after warmup. See
# http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
# to find out why this is a problem and how to eliminate them. 

# stan(model_repeat, control = list(adapt_delta = 0.99))



# RESULTS:
# 1. 1 divergent transitions after warmup
summary(model_repeat)

plot(model_repeat)
plot(model_repeat, N = 2, ask = FALSE)
plot(conditional_effects(model_repeat), points = TRUE) 

hyp <- paste(
  "sd_phylo__Intercept^2 /", 
  "(sd_phylo__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0"
)
(hyp <- hypothesis(model_repeat, hyp, class = NULL))

summary(model_repeat)

plot(hyp)
plot(model_repeat, N = 2, ask = FALSE)
plot(conditional_effects(model_repeat), points = TRUE) 




### ------ MMR -----
taxon_search.amr <- tnrs_match_names(names=unique(levels(data.amr$species)), context_name="Vertebrates") ## hhave matches
knitr::kable(taxon_search.amr)
tnrs_contexts()
# 
# tree <- rtree(n = length(unique(data.rar$species)))
# tree_l <- compute.brlen(tree, method = 'Grafen',power = 1)
# A <- ape::vcv.phylo(tree_l)


# data$ott_name <- unique_name(taxon_search)
# data$ott_id <- taxon_search$ott_id

ott_in_tree.amr <- ott_id(taxon_search.amr)[is_in_tree(ott_id(taxon_search.amr))]
tr.amr <- tol_induced_subtree(ott_ids = ott_in_tree.amr)
class(tr.amr)

tr.amr$tip.label
tr.amr$tip.label <- strip_ott_ids(tr.amr$tip.label, remove_underscores = TRUE)

labels <- as.data.frame(tr.amr$tip.label)
labels$`tr.amr$tip.label`[39]<-"Oncorhynchus mykiss"
labels$`tr.amr$tip.label`[34]<-"Tachysurus dumerili" # Leiocassis longirostris by OTL and FIshbase has Tachysurus dumerili
labels$`tr.amr$tip.label`[36]<-"Pseudobagrus vachellii" # Fishbase = "Pseudobagrus vachellii" OTL = "Tachysurus_vachellii"
labels$`tr.amr$tip.label`[61]<- "Rhinogobius giurinus" # fishbase = "Rhinogobius giurinus" OTL= Rhinogobius similis
tr.amr$tip.label <- labels$`tr.amr$tip.label`

tr.amr$tip.label %in% data.amr$species
plot(tr.amr)

# 2. look for a tree that works
studies_properties()
studies_find_trees(property = "ot:ottId", value = as.character(ott_id(taxon_search.amr)[1]))
ott_in_tree.amr <- ott_id(taxon_search.amr)[is_in_tree(ott_id(taxon_search.amr))]
tr2.amr <- tol_induced_subtree(ott_ids = ott_in_tree.amr)
plot(tr2.amr)

library(phylobase)

# dtempamr <- data.amr %>% 
#   group_by(species) %>% 
#   summarise( scaling_a_amr = mean(scaling_coef), scaling_b_amr = mean(scaling_exp), .groups = "drop") %>% 
#   as.data.frame()
# amr_numeric <-dtempamr[,-1]
# rownames(amr_numeric) <- dtempamr[,1]
# tree_data <- phylo4d(tr2.amr)

# labels$`tr$tip.label`<-gsub("_"," ",labels$`tr$tip.label`)
# colnames(labels)<-c("species")

# 1. get tree of life data


# 1.1. compute branch lengths
tr2.amr<-compute.brlen(tr.amr)
plot(tr2.amr)
tr2.amr$edge.length


A.amr <- ape::vcv.phylo(tr2.amr)
data.amr$phylo<-data.amr$species

model_repeat.amr <- brm(
  lnAMR ~ lnBWg + (1|gr(phylo, cov = A)) + (1|species), 
  data = data.amr, 
  family = gaussian(), 
  data2 = list(A = A.amr),
  prior = c(
    prior(normal(1, 1), "b"),
    prior(normal(3, 2.5), "Intercept"),
    prior(student_t(3, 0, 2.5), "sd"),
    prior(student_t(3, 0, 2.5), "sigma")),
  # prior(student_t(3, 0, 2.5), class = "sd", group = "species"),
  # prior(student_t(3, 0, 2.5), class = "sd", coef="Intercept", group = "species"),
  # prior(student_t(3, 0, 2.5), class = "sd", group = "phylo"),
  # prior(student_t(3, 0, 2.5), class = "sd", coef="Intercept", group = "phylo")),
  sample_prior = TRUE,
  chains = 2, cores = 2, 
  iter = 10000, warmup = 1000,
  thin=500,
  control = list(adapt_delta = 0.99) # << for convergence 
)

# Warning message:
#   There were 1 divergent transitions after warmup. Increasing adapt_delta above 0.99 may help. See http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup 

# get_prior(lnRMR ~ lnBWg + (1|gr(phylo, cov = A)) + (1|species), 
#           data = data.rmr, 
#           family = gaussian(), 
#           data2=list(A=A))

summary(model_repeat.amr)

plot(model_repeat.amr)
plot(model_repeat.amr, N = 2, ask = FALSE)
plot(conditional_effects(model_repeat.amr), points = TRUE) 

hyp <- paste(
  "sd_phylo__Intercept^2 /", 
  "(sd_phylo__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0"
)
(hyp <- hypothesis(model_repeat.amr, hyp, class = NULL))

summary(model_repeat.amr)







# MCMCglmm package ---------
# pr <- list(
#   R = list(V = 1, nu = 0.002),
#   G = list(G1 = list(V = 1, nu = 0.002)))
# 
# model <- MCMCglmm(lnRMR ~ lnBWg + tempTest,
#                   random = ~animal + us(1 + lnBWg):species,
#                   pedigree = tr2,
#                   # prior = pr,
#                   data = data.rmr,
#                   family = "gaussian",
#                   verbose = TRUE
# )
# 
# data.aZmr$animal<-data.amr$species
# tr2.I <- inverseA(tr2,nodes="TIPS",scale=TRUE)$Ainv
# 
# prior <- list(G=list(G1=list(V=diag(2), nu=0.002, alpha.mu=c(0,0)),
#                      G2=list(V=diag(2), nu=0.002, alpha.mu=c(0,0))),
#               R=list(V=1,nu=0.02))
# 
# 
# # run! all the way once! 
# model.rmr1 <- MCMCglmm(lnRMR ~ lnBWg + tempTest,
#                       random = ~ us(1 + lnBWg):animal,
#                       ginverse = list(animal = tr2.I),
#                       nodes="TIPS",
#                       rcov = units,
#                       # prior = prior,
#                       family = "gaussian",
#                       data = data.rmr,
#                       verbose = TRUE,
#                       nitt=10000,
#                       burnin=1000,
#                       thin=500
# )
# 
# 
# summary(model.rmr1$Sol)
# summary(model.rmr1)
# 

# PGLS from ape (Killen et al) 
# 
# 
# 
# data.amr$animal<-data.amr$species
# tr2.amrI <- inverseA(tr2.amr,nodes="TIPS",scale=TRUE)$Ainv
# 
# 
# 
# prior2<-list(G=list(G1=list(V=1,nu=0.02),G2=list(V=1,nu=0.02)),
#              R=list(V=1,nu=0.02))
# # rcov : ensures that 
# 
# model.amr <- MCMCglmm(lnAMR ~ lnBWg + tempTest,
#                       random = ~ animal + us(lnBWg):species,
#                       rcov = ~ units,
#                       prior = prior2,
#                       # pedigree = tr2.amr,
#                       ginverse = list(animal = tr2.amrI),
#                       nodes="TIPS",
#                       family = "gaussian",
#                       data = data.amr,
#                       verbose = TRUE,
#                       nitt=500000,
#                       burnin=1000,
#                       thin=500
# )
# 
# summary(model.amr)
# summary(model.amr$Sol)
# par(mfrow=c(8,2), mar=c(2,2,1,0))
# plot(model.amr$Sol)
# plot(model.amr$VCV)
# 
# 
# 
# lambda <- model_simple$VCV[,'phylo']/
#   (model_simple$VCV[,'phylo']+model_simple$VCV[,'units'])

