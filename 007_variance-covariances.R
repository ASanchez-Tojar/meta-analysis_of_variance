##############################################################
# Authors: 
# Alfredo Sanchez-Tojar (@ASanchez_Tojar)
# Profile: https://goo.gl/PmpPEB
# Department of Evolutionary Biology, Bielefeld University (GER) 
# Email: alfredo.tojar@gmail.com

# Rose O'Dea (@rose_odea)
# Profile: https://www.roseodea.com/
# I-DEEL, UNSW, Sydney (AUSTRALIA) 
# Email: rose.eleanor.o.dea@gmail.com

# Script first created on the 18th of June 2019

##############################################################
# Description of script and instructions
##############################################################

# This script is to build variance-covariance matrices as the
# sampling variances that will be included in the models for 
# the re-analysis of the data collected in:

# Eyck et al. 2019: Effects of developmental stress on animal
# phenotype and performance: a quantitative review

# This is a conservative approach that assumes a correlation
# between the effect size sampling variances with the same 
# study ID. Furthermore, we will also try to account for 
# the correlation between the effect size sampling variances
# with the same control (shared control comparisons).

# Note: a non-conservative approach would assume independent 
# sampling variances.

##############################################################
# Packages needed
##############################################################

pacman::p_load(openxlsx,corpcor,dplyr,metafor,Matrix)

# Clear memory
rm(list=ls())


##############################################################
# Functions needed
##############################################################


# Building a variance-covariance matrix for the sampling variances
# Function adjusted from: https://osf.io/ayxrt/
covMatrix <- function(data, es_var, cor){
  Var1 = sqrt(data[, es_var]) %>% unlist %>% as.numeric
  Var2 = sqrt(data[, es_var]) %>% unlist %>% as.numeric
  tmp <- expand.grid(Var1, Var2)
  tmp$cor <- ifelse(tmp$Var1 == tmp$Var2, 1, cor)
  tmp$cov <- tmp$cor * tmp$Var1 * tmp$Var2
  corMat <- matrix(tmp$cov , nrow = nrow(data), ncol = nrow(data))
  return(corMat)
}

# this function is to make sure that a is matrix positive-definitive
# it uses the function nearPD in the Matrix package to compute the 
# nearest positive definite matrix to an approximate one 
# Function reused from: https://osf.io/ayxrt/
PDfunc <- function(matrix){
  require(Matrix)
  if(corpcor::is.positive.definite(matrix) == T){
    x <- corDiag
  }else{
    x <- Matrix::nearPD(matrix)
  }
  mat <- x$mat
  return(mat)
}


##############################################################
# Importing dataset
##############################################################

# database with the corrected data from our pilot re-extraction
stress.data <- read.xlsx("data_re-extraction/clean_data/EyckDev_stress_clean_effect_sizes_sp_corrected.xlsx",
                         colNames=T,sheet = 1)

#test<-stress.data[stress.data$studyID=="114" | stress.data$studyID=="80",]

# subsetting data
stress.data.ours <- stress.data[!(is.na(stress.data$SMDH.ours)),]
stress.data.HE <- stress.data[!(is.na(stress.data$SMDH.HE)),]

##############################################################
# METAFOR and BRMS: variance-covariance matrices
##############################################################

##############################################################
# lnRR.ours
##############################################################

# subset of data needed
stress.data.ours.lnRR <- stress.data.ours[!(is.na(stress.data.ours$lnRR.sv)),]

##################################
# within-study correlation = 0.5
##################################

# generating a within-studyID variance-covariance matrix 
# assuming a correlation of 0.5 for within-study sampling 
# variances 

# The way the function works is: 1) Split data; 2) lapply to 
# list of each paper; 3) create block diagonal matrix

# test$studyID <- factor(test$studyID, levels=unique(test$studyID))
# splt <- split(test, test$studyID)

stress.data.ours.lnRR$studyID <- factor(stress.data.ours.lnRR$studyID, levels=unique(stress.data.ours.lnRR$studyID))
splt <- split(stress.data.ours.lnRR, stress.data.ours.lnRR$studyID)
covMat <- lapply(splt, function(x) covMatrix(x, "lnRR.sv", cor = 0.5))
varcovar.studyID.lnRR.ours_0.5 <- bldiag(covMat)

# is the matrix positive-definitive? No, run PDfunc to make it
# so. If a matrix is not positive-definite, it won't work
is.positive.definite(varcovar.studyID.lnRR.ours_0.5) # FALSE
varcovar.studyID.lnRR.ours_0.5<-PDfunc(varcovar.studyID.lnRR.ours_0.5)
is.positive.definite(varcovar.studyID.lnRR.ours_0.5)

# all good, save it
save(varcovar.studyID.lnRR.ours_0.5, file = "data_re-extraction/clean_data/varcovar_studyID_lnRR_ours_0.50.Rdata")



##############################################################
# lnVR.ours
##############################################################

# subset of data needed
stress.data.ours.lnVR <- stress.data.ours[!(is.na(stress.data.ours$lnVR.sv)),]

##################################
# within-study correlation = 0.5
##################################

# generating a within-studyID variance-covariance matrix 
# assuming a correlation of 0.5 for within-study sampling 
# variances 

# The way the function works is: 1) Split data; 2) lapply to 
# list of each paper; 3) create block diagonal matrix

# test$studyID <- factor(test$studyID, levels=unique(test$studyID))
# splt <- split(test, test$studyID)

stress.data.ours.lnVR$studyID <- factor(stress.data.ours.lnVR$studyID, levels=unique(stress.data.ours.lnVR$studyID))
splt <- split(stress.data.ours.lnVR, stress.data.ours.lnVR$studyID)
covMat <- lapply(splt, function(x) covMatrix(x, "lnVR.sv", cor = 0.5))
varcovar.studyID.lnVR.ours_0.5 <- bldiag(covMat)

# is the matrix positive-definitive? No, run PDfunc to make it
# so. If a matrix is not positive-definite, it won't work
is.positive.definite(varcovar.studyID.lnVR.ours_0.5) # FALSE
varcovar.studyID.lnVR.ours_0.5<-PDfunc(varcovar.studyID.lnVR.ours_0.5)
is.positive.definite(varcovar.studyID.lnVR.ours_0.5)

# all good, save it
save(varcovar.studyID.lnVR.ours_0.5, file = "data_re-extraction/clean_data/varcovar_studyID_lnVR_ours_0.50.Rdata")



##############################################################
# lnCVR.ours
##############################################################

# subset of data needed
stress.data.ours.lnCVR <- stress.data.ours[!(is.na(stress.data.ours$lnCVR.sv)),]

##################################
# within-study correlation = 0.5
##################################

# generating a within-studyID variance-covariance matrix 
# assuming a correlation of 0.5 for within-study sampling 
# variances 

# The way the function works is: 1) Split data; 2) lapply to 
# list of each paper; 3) create block diagonal matrix

# test$studyID <- factor(test$studyID, levels=unique(test$studyID))
# splt <- split(test, test$studyID)

stress.data.ours.lnCVR$studyID <- factor(stress.data.ours.lnCVR$studyID, levels=unique(stress.data.ours.lnCVR$studyID))
splt <- split(stress.data.ours.lnCVR, stress.data.ours.lnCVR$studyID)
covMat <- lapply(splt, function(x) covMatrix(x, "lnCVR.sv", cor = 0.5))
varcovar.studyID.lnCVR.ours_0.5 <- bldiag(covMat)

# is the matrix positive-definitive? No, run PDfunc to make it
# so. If a matrix is not positive-definite, it won't work
is.positive.definite(varcovar.studyID.lnCVR.ours_0.5) # FALSE
varcovar.studyID.lnCVR.ours_0.5<-PDfunc(varcovar.studyID.lnCVR.ours_0.5)
is.positive.definite(varcovar.studyID.lnCVR.ours_0.5)

# all good, save it
save(varcovar.studyID.lnCVR.ours_0.5, file = "data_re-extraction/clean_data/varcovar_studyID_lnCVR_ours_0.50.Rdata")



##############################################################
# SMDH.ours
##############################################################

##################################
# within-study correlation = 0.5
##################################

# generating a within-studyID variance-covariance matrix 
# assuming a correlation of 0.5 for within-study sampling 
# variances 

# The way the function works is: 1) Split data; 2) lapply to 
# list of each paper; 3) create block diagonal matrix

# test$studyID <- factor(test$studyID, levels=unique(test$studyID))
# splt <- split(test, test$studyID)

stress.data.ours$studyID <- factor(stress.data.ours$studyID, levels=unique(stress.data.ours$studyID))
splt <- split(stress.data.ours, stress.data.ours$studyID)
covMat <- lapply(splt, function(x) covMatrix(x, "SMDH.sv", cor = 0.5))
varcovar.studyID.SMDH.ours_0.5 <- bldiag(covMat)

# is the matrix positive-definitive? No, run PDfunc to make it
# so. If a matrix is not positive-definite, it won't work
is.positive.definite(varcovar.studyID.SMDH.ours_0.5) # FALSE
varcovar.studyID.SMDH.ours_0.5<-PDfunc(varcovar.studyID.SMDH.ours_0.5)
is.positive.definite(varcovar.studyID.SMDH.ours_0.5)

# all good, save it
save(varcovar.studyID.SMDH.ours_0.5, file = "data_re-extraction/clean_data/varcovar_studyID_SMDH_ours_0.50.Rdata")



##############################################################
# MCMCglmm: variance-covariance matrices
##############################################################

##############################################################
# SMDH.ours
##############################################################

# # Losia's Function for creating covarianc matrix for MCMCglmm
# # data is the dataframe
# # GroupID is a column with unique ID joining ES that are not independednt within each study (blocks of correlated ES)
# # es_var is variance of effect sizes
# # obs is the column with unique ID for each data row (effect size)
# 
# VC_Matrix <- function(data, es_var, GroupID, obs){
#   vc_matrix <- matrix(0,nrow = dim(data)[1],ncol = dim(data)[1]) #make empty matrix of the same size as data length
#   rownames(vc_matrix) <- data[ ,obs]
#   colnames(vc_matrix) <- data[ ,obs]
#   # find start and end coordinates for the subsets
#   shared_coord <- which(data[ ,GroupID] %in% data[duplicated(data[ ,GroupID]), GroupID]==TRUE)
#   # matrix of combinations of coordinates for each experiment with shared control
#   combinations <- do.call("rbind", tapply(shared_coord, data[shared_coord,GroupID], function(x) t(combn(x,2))))
#   # calculate covariance values between  values at the positions in shared_list and place them on the matrix
#   for (i in 1:dim(combinations)[1]){
#     p1 <- combinations[i,1]
#     p2 <- combinations[i,2]
#     p1_p2_cov <- 0.5*sqrt(data[p1,es_var]) * sqrt(data[p2,es_var])
#     vc_matrix[p1,p2] <- p1_p2_cov
#     vc_matrix[p2,p1] <- p1_p2_cov
#   }
#   diag(vc_matrix) <- data[ ,es_var]   #add the diagonal
#   return(vc_matrix)
# }


# saving session information with all packages versions for reproducibility purposes
sink("data_re-extraction/clean_data/variance_covariance_R_session.txt")
sessionInfo()
sink()