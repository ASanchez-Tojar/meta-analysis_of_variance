##############################################################
# Authors: 
# Alfredo Sanchez-Tojar (@ASanchez_Tojar)
# Profile: https://goo.gl/PmpPEB
# Department of Evolutionary Biology, Bielefeld University (GER) 
# Email: alfredo.tojar@gmail.com

# Script first created on the 8th of July 2019

##############################################################
# Description of script and instructions
##############################################################

# This script is to first reduce the dataset by reducing the 
# number of lnRR and lnCVR effect sizes available, and then
# to re-analyze the data collected in:

# Eyck et al. 2019: Effects of developmental stress on animal
# phenotype and performance: a quantitative review

# We use the R package 'brms' for the analyses. The idea of 
# these reductions is to show how the borrowing of strength 
# takes place when using bivariate models.

##############################################################
# Packages needed
##############################################################

pacman::p_load(openxlsx,brms,tidybayes,ggplot2,plotly,stringr,dplyr)


# Clear memory
rm(list=ls())


##############################################################
# Functions needed
##############################################################

# none

##############################################################
# Importing dataset
##############################################################

# database with the corrected data
stress.data <- read.xlsx("data_re-extraction/clean_data/EyckDev_stress_clean_effect_sizes_sp_corrected.xlsx",
                         colNames=T,sheet = 1)

# loading phylogenetic matrix "phylo_cor"
load("data_re-extraction/clean_data/phylo_cor.Rdata") #phylo_cor

##############################################################
# Reducing dataset
##############################################################

# # since we are running our analyses in R version 3.5.1, it is
# # still possible to use set.seed() to make the randomization
# # reproducible (though they may have fixed this 'bug' already!)
# 
# set.seed(2+13) #KL+PG
# 
# # we want to make the following data reductions: 25% and 50%
# # We are going to based our reduction in the whole dataset,
# # which contains 711 effect sizes, despite that it changes slightly
# # from 677 to 708 depending on whether it is based on lnRR,
# # lnVR or lnCVR
# 
# twentyfive.per.cent <- round(nrow(stress.data)*0.25,0)
# fifty.per.cent <- round(nrow(stress.data)*0.50,0)
# 
# # let's first reset the row names to have them from 1:nrow()
# # and make choosing them randomly easier, also, it looks neater
# rownames(stress.data) <- NULL
# 
# 
# # randomly choosing 25 and 50% of the effect sizes
# twentyfive.per.cent.list <- base::sample(c(1:nrow(stress.data)),
#                                          twentyfive.per.cent)
# 
# fifty.per.cent.list <- base::sample(c(1:nrow(stress.data)),
#                                     fifty.per.cent)
# 
# 
# # removing those randomly selected effect sizes from the database
# # They way we are setting the randomization makes sure that the
# # same esID is removed for both lnRR and lnCVR. This makes sense
# # because what we are trying to simulate is a scenario where
# # an X% of the data is non ratio scale, and that's why lnRR and
# # lnCVR could not be calculated, something that does not affect
# # the lnVR. Note that the specific % are approximate.
# 
# # duplicating variables to be selectively emptied
# stress.data$lnRR.sc.ours.25 <- stress.data$lnRR.sc.ours
# stress.data$lnRR.sc.ours.50 <- stress.data$lnRR.sc.ours
# stress.data$lnCVR.sc.25 <- stress.data$lnCVR.sc
# stress.data$lnCVR.sc.50 <- stress.data$lnCVR.sc
# 
# # selectively reducing varibles
# stress.data[twentyfive.per.cent.list,"lnRR.sc.ours.25"] <- NA
# stress.data[fifty.per.cent.list,"lnRR.sc.ours.50"] <- NA
# stress.data[twentyfive.per.cent.list,"lnCVR.sc.25"] <- NA
# stress.data[fifty.per.cent.list,"lnCVR.sc.50"] <- NA
# 
# 
# # saving dataset
# write.xlsx(stress.data,
#            "data_re-extraction/clean_data/EyckDev_stress_clean_effect_sizes_sp_corrected_reduction_ours.xlsx",
#            sheetName="Sheet1",col.names=TRUE, row.names=F,
#            append=FALSE, showNA=TRUE, password=NULL)

# database with the reduced data 
stress.data <- read.xlsx("data_re-extraction/clean_data/EyckDev_stress_clean_effect_sizes_sp_corrected_reduction_ours.xlsx",
                         colNames=T,sheet = 1)

##############################################################
# -------------------------- BRMS -------------------------- #
##############################################################

# We run univariate and bivariate multilevel meta-analyses 
# in brms.

# Using brms, we will run univariate models based on the effect 
# sizes we calculated ourselves: lnRR, lnVR and lnCVR.

# Then we will run bivariate models as: c(lnRR,lnVR),
# and c(lnVR,lnCVR)

# The main models will be run based on the effect sizes that
# account for shared control non-independence (see script 004),
# plus, in this case, we based the models on a reduced dataset
# to show the borrowing of strength effect.

# model specifications
adapt_delta_value <- 0.9999
max_treedepth_value <- 20
iterations <- 6000
burnin <- 3000
thinning <- 2

###########################
# UNIVARIATE MODELS: lnRR #
###########################

# sharedcontrol - 25%: 0.999, 20, 6000, 3000, 2:  73 min (corei7)
# sharedcontrol - 50%: 0.999, 20, 6000, 3000, 2:  34 min (corei7)

# subset of data needed just to make it easier for the univariate
# models instead of going through tons of NA's.
stress.data.lnRR.25 <- stress.data[!(is.na(stress.data$lnRR.sc.ours.25)),]
stress.data.lnRR.50 <- stress.data[!(is.na(stress.data$lnRR.sc.ours.50)),]


#################
# 25% reduction #
#################

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_univariate_lnRR_ours_25_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")


brms.univariate.lnRR.ours.25 <- brm(lnRR.sc.ours.25 | se(sqrt(lnRR.sc.sv)) ~ 1 + 
                                      (1|studyID) + (1|esID) +
                                      (1|scientific.name) + (1|speciesID),
                                    data = stress.data.lnRR.25,
                                    family = gaussian(),
                                    cov_ranef = list(scientific.name = phylo_cor),
                                    control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                    chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.univariate.lnRR.ours.25,
     file=filename)


#################
# 50% reduction #
#################

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_univariate_lnRR_ours_50_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")


brms.univariate.lnRR.ours.50 <- brm(lnRR.sc.ours.50 | se(sqrt(lnRR.sc.sv)) ~ 1 + 
                                      (1|studyID) + (1|esID) +
                                      (1|scientific.name) + (1|speciesID),
                                    data = stress.data.lnRR.50,
                                    family = gaussian(),
                                    cov_ranef = list(scientific.name = phylo_cor),
                                    control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                    chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.univariate.lnRR.ours.50,
     file=filename)


############################
# UNIVARIATE MODELS: lnCVR #
############################

# sharedcontrol - 25%: 0.999, 20, 6000, 3000, 2: 4 min (corei7)
# sharedcontrol - 50%: 0.999, 20, 6000, 3000, 2: 3  min (corei7)

# subset of data needed just to make it easier for the univariate
# models instead of going through tons of NA's.
stress.data.lnCVR.25 <- stress.data[!(is.na(stress.data$lnCVR.sc.25)),]
stress.data.lnCVR.50 <- stress.data[!(is.na(stress.data$lnCVR.sc.50)),]


#################
# 25% reduction #
#################

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_univariate_lnCVR_25_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")


brms.univariate.lnCVR.25 <- brm(lnCVR.sc.25 | se(sqrt(lnCVR.sc.sv)) ~ 1 + 
                                       (1|studyID) + (1|esID) +
                                       (1|scientific.name) + (1|speciesID),
                                     data = stress.data.lnCVR.25,
                                     family = gaussian(),
                                     cov_ranef = list(scientific.name = phylo_cor),
                                     control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                     chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.univariate.lnCVR.25,
     file=filename)


#################
# 50% reduction #
#################

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_univariate_lnCVR_50_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")


brms.univariate.lnCVR.50 <- brm(lnCVR.sc.50 | se(sqrt(lnCVR.sc.sv)) ~ 1 + 
                                       (1|studyID) + (1|esID) +
                                       (1|scientific.name) + (1|speciesID),
                                     data = stress.data.lnCVR.50,
                                     family = gaussian(),
                                     cov_ranef = list(scientific.name = phylo_cor),
                                     control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                     chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.univariate.lnCVR.50,
     file=filename)


###################
# BIVARIATE MODELS: 
###################

# specifying models' structure
bf.lnRR.ours.25 <- bf(lnRR.sc.ours.25 | se(sqrt(lnRR.sc.sv)) ~
                        1 + (1|p|studyID) + (1|q|esID) + (1|a|scientific.name) + (1|d|speciesID))

bf.lnRR.ours.50 <- bf(lnRR.sc.ours.50 | se(sqrt(lnRR.sc.sv)) ~
                        1 + (1|p|studyID) + (1|q|esID) + (1|a|scientific.name) + (1|d|speciesID))

bf.lnVR <- bf(lnVR.sc | se(sqrt(lnVR.sc.sv)) ~
                1 + (1|p|studyID) + (1|q|esID) + (1|a|scientific.name) + (1|d|speciesID))

bf.lnCVR.25 <- bf(lnCVR.sc.25 | se(sqrt(lnCVR.sc.sv)) ~
                    1 + (1|p|studyID) + (1|q|esID) + (1|a|scientific.name) + (1|d|speciesID))

bf.lnCVR.50 <- bf(lnCVR.sc.50 | se(sqrt(lnCVR.sc.sv)) ~
                    1 + (1|p|studyID) + (1|q|esID) + (1|a|scientific.name) + (1|d|speciesID))


###################
# cbind(lnRR,lnVR) 
###################

# sharedcontrol - 25%: 0.999, 20, 6000, 3000, 2:  1385 min   (corei7)
# sharedcontrol - 50%: 0.999, 20, 6000, 3000, 2:  896 min (corei7)


#################
# 25% reduction #
#################

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_bivariate_lnRR_ours_lnVR_25_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")


brms.bivariate.lnRR.ours.lnVR.25 <- brm(bf.lnRR.ours.25 + bf.lnVR,
                                        data = stress.data,
                                        cov_ranef = list(scientific.name = phylo_cor),
                                        family = gaussian(),
                                        control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                        chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.bivariate.lnRR.ours.lnVR.25,
     file=filename)


#################
# 50% reduction #
#################

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_bivariate_lnRR_ours_lnVR_50_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")


brms.bivariate.lnRR.ours.lnVR.50 <- brm(bf.lnRR.ours.50 + bf.lnVR,
                                        data = stress.data,
                                        cov_ranef = list(scientific.name = phylo_cor),
                                        family = gaussian(),
                                        control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                        chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.bivariate.lnRR.ours.lnVR.50,
     file=filename)


###################
# cbind(lnVR,lnCVR) 
###################

# sharedcontrol - 25%: 0.999, 20, 6000, 3000, 2:  201 min (corei7)
# sharedcontrol - 50%: 0.999, 20, 6000, 3000, 2:  126 min (corei7)

#################
# 25% reduction #
#################

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_bivariate_lnVR_lnCVR_25_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")


brms.bivariate.lnVR.lnCVR.25 <- brm(bf.lnVR + bf.lnCVR.25,
                                    data = stress.data,
                                    cov_ranef = list(scientific.name = phylo_cor),
                                    family = gaussian(),
                                    control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                    chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.bivariate.lnVR.lnCVR.25,
     file=filename)


#################
# 50% reduction #
#################

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_bivariate_lnVR_lnCVR_50_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")


brms.bivariate.lnVR.lnCVR.50 <- brm(bf.lnVR + bf.lnCVR.50,
                                    data = stress.data,
                                    cov_ranef = list(scientific.name = phylo_cor),
                                    family = gaussian(),
                                    control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                    chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.bivariate.lnVR.lnCVR.50,
     file=filename)


####################################################################################
# saving session information with all packages versions for reproducibility purposes
sink("models/brms/brms_meta-analysis_reductions_R_session.txt")
sessionInfo()
sink()