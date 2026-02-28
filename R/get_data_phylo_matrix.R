# testing: ******
# species.list = unique(levels(data.amrER$species))
# matrix.name = "A.mmr.er"
# dataset.ID = "MMR optimal"

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

