get_phylo_matrix<-function(species.list, matrix.name, tree.name, dataset.ID){

  taxon_search <- tnrs_match_names(names=species.list, context_name="Vertebrates") 
  ott_in_tree <- ott_id(taxon_search)[is_in_tree(ott_id(taxon_search))]
  tr <- tol_induced_subtree(ott_ids = ott_in_tree)
  
  tr$tip.label <- strip_ott_ids(tr$tip.label, remove_underscores = TRUE)
  labels <- as.data.frame(tr$tip.label)

  labels$`tr$tip.label`[which(labels$`tr$tip.label` == "Oncorhynchus mykiss (species in domain Eukaryota)")]<-"Oncorhynchus mykiss"
  labels$`tr$tip.label`[which(labels$`tr$tip.label` == "Gadus morhua (species in domain Eukaryota)")]<-"Gadus morhua"
  labels$`tr$tip.label`[which(labels$`tr$tip.label` == "Leiocassis longirostris")]<-"Tachysurus dumerili" 
  labels$`tr$tip.label`[which(labels$`tr$tip.label` == "Tachysurus vachellii")]<-"Pseudobagrus vachellii" 
  labels$`tr$tip.label`[which(labels$`tr$tip.label` == "Rhinogobius similis")]<- "Rhinogobius giurinus" 
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

  assign(matrix.name, value = A, envir = .GlobalEnv)
  assign(tree.name, value = tree, envir = .GlobalEnv)

}
