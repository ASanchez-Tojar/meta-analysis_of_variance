##############################################################
# Authors: 
# Alfredo Sanchez-Tojar (@ASanchez_Tojar)
# Profile: https://goo.gl/PmpPEB
# Department of Evolutionary Biology, Bielefeld University (GER) 
# Email: alfredo.tojar@gmail.com

# Script first created on the 17th of June 2019

##############################################################
# Description of script and instructions
##############################################################

# This script is to build the phylogeny and estimate the 
# phylogenetic relatedness among the species included in:

# Eyck et al. 2019: Effects of developmental stress on animal
# phenotype and performance: a quantitative review


##############################################################
# Packages needed
##############################################################

pacman::p_load(openxlsx,readxl, ape, fulltext, metafor,rotl,
               treebase,diagram,dplyr)

# Clear memory
rm(list=ls())


##############################################################
# Functions needed
##############################################################

# none

##############################################################
# Importing datasets
##############################################################

# final and clean database 
stress.data <- read.xlsx("data_re-extraction/clean_data/EyckDev_stress_clean_effect_sizes.xlsx",
                         colNames=T,sheet = 1)


# generating list of species
species <- sort(unique(stress.data$scientific.name))


##############################################################
# Formatting species data
##############################################################

# obtaining dataframe listing the Open Tree identifiers potentially 
# matching our list of species.
taxa <- tnrs_match_names(names = species)


# according to the `approximate_match` column, there might be 
# 2 typos in the species list (well actually 1 with the final subset)
# nrow(taxa[taxa$approximate_match==TRUE,])
taxa[taxa$approximate_match==TRUE,]


# fixing those typos
species[species=="Salmo salari"] <- "Salmo salar"
#species[species=="Basal Tarantula"] <- "Brachypelma smithi" #info obtained from the study

stress.data[stress.data$scientific.name=="Salmo salari","scientific.name"] <- "Salmo salar"
#stress.data[stress.data$scientific.name=="Basal Tarantula",] <- "Brachypelma smithi" #info obtained from the study
stress.data[stress.data$scientific.name=="rattus norvegicus","scientific.name"] <- "Rattus norvegicus"

# rerun
taxa.c <- tnrs_match_names(names = species)

# exploring which species return more than one match, and the
# reasons to make sure we retrieve the correct data.
taxa.c[taxa.c$number_matches != 1,]
ott_id_tocheck <- taxa.c[taxa.c$number_matches != 1,"ott_id"]

for(i in 1:length(ott_id_tocheck)){
  print(inspect(taxa.c, ott_id = ott_id_tocheck[i]))
}

# everything seems good so far


##############################################################
# Retrieving phylogenetic relationships
##############################################################

# retrieving phylogenetic relationships among taxa in the form 
# of a trimmed sub-tree
tree <- tol_induced_subtree(ott_ids = taxa.c[["ott_id"]], label_format = "name")
plot(tree, cex=.5, label.offset =.1, no.margin = TRUE)


# Notice that the species names shown in the tree are not exactly 
# the same as the species names that we had in our list. This is 
# because those names had synonyms in the tree of life database, 
# and we are using those names for the plot.


##############################################################
# Dealing with polytomies
##############################################################

# we can check for the existence of polytomies by running the 
# following code. If polytomies exist, the output will be 
# `FALSE`, and vice versa.

is.binary.tree(tree) # there are some polytomies


# to take care of these polytomies, we are goint to use a 
# randomization approach
set.seed(111) #making it replicable, at least for this version of R (i.e. v.3.5.1)
tree_random <- multi2di(tree,random=TRUE)
is.binary.tree(tree_random)


##############################################################
# Final checks
##############################################################

# exploring whether our tree covers all the species we wanted 
# it to include, and makeing sure that the species names in our 
# database match those in the tree. We use the following code.

tree_random$tip.label <- gsub("_"," ", tree_random$tip.label)
intersect(as.character(tree_random$tip.label), as.character(species))
setdiff(species, as.character(tree_random$tip.label)) #listed in our database but not in the tree
setdiff(as.character(tree_random$tip.label),species) # listed in the tree but not in our database


# they are the same species, the "problem" is that synonyms 
# have been used in the tree (see also `taxa.c`). We are going
# to leave all the names as in Open Tree of Life as it seems
# to be the most updated nomenclature

# we'll just fix this one
tree_random$tip.label[49]<-"Oncorhynchus mykiss"

tiff("plots/phylogenetic_tree_pruned.tiff",
     height=20, width=10,
     units='cm', compression="lzw", res=800)

plot(tree_random, cex=.5, label.offset =.1, no.margin = TRUE)

dev.off()

# we can now save the tree
save(tree_random, file = "data_re-extraction/clean_data/tree_random.Rdata")


##############################################################
# Computing branch lengths
##############################################################

# we are computing branch lengths for our tree following 
# Grafen (1989)(https://royalsocietypublishing.org/doi/abs/10.1098/rstb.1989.0106)

# before we now need to make sure that tree labels and database
# use the same nomenclature
setdiff(stress.data$scientific.name, as.character(tree_random$tip.label))
setdiff(as.character(tree_random$tip.label),stress.data$scientific.name)

tree_random.fixed <- tree_random
#tree_random.fixed$tip.label <- gsub("Rhinella marina","Bufo marinus", tree_random.fixed$tip.label)
#tree_random.fixed$tip.label <- gsub("Spea hammondii","Scaphiopus hammondii", tree_random.fixed$tip.label)
tree_random.fixed$tip.label <- gsub("Gulosus aristotelis","Phalacrocorax aristotelis", tree_random.fixed$tip.label)
tree_random.fixed$tip.label <- gsub("Pseudosimochromis pleurospilus","Simochromis pleurospilus", tree_random.fixed$tip.label)
tree_random.fixed$tip.label <- gsub("Coturnix japonica","Coturnix coturnix japonica", tree_random.fixed$tip.label)
tree_random.fixed$tip.label <- gsub("Saproscincus mustelinus","Saproscincus mustelina", tree_random.fixed$tip.label)
#tree_random.fixed$tip.label <- gsub("Oncorhynchus mykiss (species in domain Eukaryota)","Oncorhynchus mykiss", tree_random.fixed$tip.label)
tree_random.fixed$tip.label <- gsub("Litoria ewingii","Litoria ewingi", tree_random.fixed$tip.label)
tree_random.fixed$tip.label <- gsub("Gallus gallus","Gallus gallus domesticus", tree_random.fixed$tip.label)

# nothing done with "Rana sylvatica"/"Lithobates sylvaticus"
# because both nomenclatures are present in our database, 
# however, we need to choose one for our database:

stress.data[stress.data$scientific.name=="Lithobates sylvaticus","scientific.name"] <- "Rana sylvatica"

setdiff(stress.data$scientific.name, as.character(tree_random.fixed$tip.label))
setdiff(as.character(tree_random.fixed$tip.label),stress.data$scientific.name)
# all good!


# compute branch lengths of tree
phylo_branch <- compute.brlen(tree_random.fixed, method = "Grafen", power = 1)


# check tree is ultrametric
is.ultrametric(phylo_branch) # TRUE


##############################################################
# Phylogenetic matrix
##############################################################

# matrix to be included in the models
phylo_cor <- vcv(phylo_branch, cor = T)


# finally, save matrix for future analyses
save(phylo_cor, file = "data_re-extraction/clean_data/phylo_cor.Rdata")


# exporting fixed dataset for analyses
write.xlsx(stress.data,
           "data_re-extraction/clean_data/EyckDev_stress_clean_effect_sizes_sp_corrected.xlsx",
           sheetName="Sheet1",col.names=TRUE, row.names=F,
           append=FALSE, showNA=TRUE, password=NULL)

# saving session information with all packages versions for reproducibility purposes
sink("data_re-extraction/clean_data/phylogeny_reconstruction_R_session.txt")
sessionInfo()
sink()
