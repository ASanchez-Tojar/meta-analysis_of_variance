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

# database with the corrected data from our pilot re-extraction
stress.data <- read.xlsx("data_re-extraction/clean_data/EyckDev_stress_clean_effect_sizes_sp_corrected.xlsx",
                         colNames=T,sheet = 1)

# subsetting data
stress.data.ours <- stress.data[!(is.na(stress.data$SMDH.ours)),]

# loading phylogenetic matrix "phylo_cor"
load("data_re-extraction/clean_data/phylo_cor.Rdata") #phylo_cor

##############################################################
# Reducing dataset
##############################################################

# # since we are running our analyses in R version 3.5.1, it is
# # still possible to use set.seed() to make the randomization
# # reproducible (though they may have fixed this bug already!)
# 
# set.seed(2+13) #KL+PG
# 
# # we want to make the following data reductions: 10%, 30% and 50%
# # We are going to based our reduction in the whole dataset,
# # which contains 684 effect sizes, despite that changes slightly
# # from 674 to 681 depending on whether it is based on lnRR,
# # lnVR or lnCVR
# 
# ten.per.cent <- round(length(stress.data.ours$SMDH.sc.ours)*0.10,0)
# thirty.per.cent <- round(length(stress.data.ours$SMDH.sc.ours)*0.30,0)
# fifty.per.cent <- round(length(stress.data.ours$SMDH.sc.ours)*0.50,0)
# 
# # let's first reset the row names to have them from 1:nrow()
# # and make choosing them randomly easier, also, it looks neater
# rownames(stress.data.ours) <- NULL
# 
# 
# # randomly choosing 10, 30 and 50% of the effect sizes 
# ten.per.cent.list <- base::sample(c(1:nrow(stress.data.ours)),
#                                   ten.per.cent)
# 
# thirty.per.cent.list <- base::sample(c(1:nrow(stress.data.ours)),
#                                      thirty.per.cent)
# 
# fifty.per.cent.list <- base::sample(c(1:nrow(stress.data.ours)),
#                                     fifty.per.cent)
# 
# 
# # removing those randomly selected effect sizes from the database
# # They way we are setting the randomization makes sure that the
# # same esID is removed for both lnRR and lnCVR. This makes sense
# # because what we are trying to simulate is a scenario where
# # an X% of the data is non ratio scale, and that's why lnRR and
# # lnCVR could not be calculated, something that does not affect
# # the lnVR
# 
# # duplicating variables to be selectively emptied
# stress.data.ours$lnRR.sc.ours.10 <- stress.data.ours$lnRR.sc.ours
# stress.data.ours$lnRR.sc.ours.30 <- stress.data.ours$lnRR.sc.ours
# stress.data.ours$lnRR.sc.ours.50 <- stress.data.ours$lnRR.sc.ours
# stress.data.ours$lnCVR.sc.ours.10 <- stress.data.ours$lnCVR.sc.ours
# stress.data.ours$lnCVR.sc.ours.30 <- stress.data.ours$lnCVR.sc.ours
# stress.data.ours$lnCVR.sc.ours.50 <- stress.data.ours$lnCVR.sc.ours
# 
# # selectively reducing varibles
# stress.data.ours[ten.per.cent.list,"lnRR.sc.ours.10"] <- NA
# stress.data.ours[thirty.per.cent.list,"lnRR.sc.ours.30"] <- NA
# stress.data.ours[fifty.per.cent.list,"lnRR.sc.ours.50"] <- NA
# stress.data.ours[ten.per.cent.list,"lnCVR.sc.ours.10"] <- NA
# stress.data.ours[thirty.per.cent.list,"lnCVR.sc.ours.30"] <- NA
# stress.data.ours[fifty.per.cent.list,"lnCVR.sc.ours.50"] <- NA
# 
# 
# # saving dataset
# write.xlsx(stress.data.ours,
#            "data_re-extraction/clean_data/EyckDev_stress_clean_effect_sizes_sp_corrected_reduction_ours.xlsx",
#            sheetName="Sheet1",col.names=TRUE, row.names=F,
#            append=FALSE, showNA=TRUE, password=NULL)

# database with the corrected data from our pilot re-extraction
stress.data.ours <- read.xlsx("data_re-extraction/clean_data/EyckDev_stress_clean_effect_sizes_sp_corrected_reduction_ours.xlsx",
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
adapt_delta_value <- 0.999
max_treedepth_value <- 20
iterations <- 6000
burnin <- 3000
thinning <- 2

###########################
# UNIVARIATE MODELS: lnRR #
###########################

# sharedcontrol - 10%: 0.999, 20, 6000, 3000, 2: XX min (corei7)
# sharedcontrol - 30%: 0.999, 20, 6000, 3000, 2: XX min (corei7)
# sharedcontrol - 50%: 0.999, 20, 6000, 3000, 2: XX min (corei7)

# subset of data needed just to make it easier for the univariate
# models instead of going through tons of NA's.
stress.data.ours.lnRR.10 <- stress.data.ours[!(is.na(stress.data.ours$lnRR.sc.ours.10)),]
stress.data.ours.lnRR.30 <- stress.data.ours[!(is.na(stress.data.ours$lnRR.sc.ours.30)),]
stress.data.ours.lnRR.50 <- stress.data.ours[!(is.na(stress.data.ours$lnRR.sc.ours.50)),]


#################
# 10% reduction #
#################

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_univariate_lnRR_ours_10_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")


brms.univariate.lnRR.ours.10 <- brm(lnRR.sc.ours.10 | se(sqrt(lnRR.sc.sv)) ~ 1 + 
                                      (1|studyID) + (1|esID) +
                                      (1|scientific.name) + (1|speciesID),
                                    data = stress.data.ours.lnRR.10,
                                    family = gaussian(),
                                    cov_ranef = list(scientific.name = phylo_cor),
                                    control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                    chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.univariate.lnRR.ours.10,
     file=filename)


#################
# 30% reduction #
#################

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_univariate_lnRR_ours_30_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")


brms.univariate.lnRR.ours.30 <- brm(lnRR.sc.ours.30 | se(sqrt(lnRR.sc.sv)) ~ 1 + 
                                      (1|studyID) + (1|esID) +
                                      (1|scientific.name) + (1|speciesID),
                                    data = stress.data.ours.lnRR.30,
                                    family = gaussian(),
                                    cov_ranef = list(scientific.name = phylo_cor),
                                    control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                    chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.univariate.lnRR.ours.30,
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
                                    data = stress.data.ours.lnRR.50,
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

# sharedcontrol - 10%: 0.999, 20, 6000, 3000, 2: XX min (corei7)
# sharedcontrol - 30%: 0.999, 20, 6000, 3000, 2: XX min (corei7)
# sharedcontrol - 50%: 0.999, 20, 6000, 3000, 2: XX min (corei7)

# subset of data needed just to make it easier for the univariate
# models instead of going through tons of NA's.
stress.data.ours.lnCVR.10 <- stress.data.ours[!(is.na(stress.data.ours$lnCVR.sc.ours.10)),]
stress.data.ours.lnCVR.30 <- stress.data.ours[!(is.na(stress.data.ours$lnCVR.sc.ours.30)),]
stress.data.ours.lnCVR.50 <- stress.data.ours[!(is.na(stress.data.ours$lnCVR.sc.ours.50)),]


#################
# 10% reduction #
#################

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_univariate_lnCVR_ours_10_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")


brms.univariate.lnCVR.ours.10 <- brm(lnCVR.sc.ours.10 | se(sqrt(lnCVR.sc.sv)) ~ 1 + 
                                       (1|studyID) + (1|esID) +
                                       (1|scientific.name) + (1|speciesID),
                                     data = stress.data.ours.lnCVR.10,
                                     family = gaussian(),
                                     cov_ranef = list(scientific.name = phylo_cor),
                                     control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                     chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.univariate.lnCVR.ours.10,
     file=filename)


#################
# 30% reduction #
#################

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_univariate_lnCVR_ours_30_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")


brms.univariate.lnCVR.ours.30 <- brm(lnCVR.sc.ours.30 | se(sqrt(lnCVR.sc.sv)) ~ 1 + 
                                       (1|studyID) + (1|esID) +
                                       (1|scientific.name) + (1|speciesID),
                                     data = stress.data.ours.lnCVR.30,
                                     family = gaussian(),
                                     cov_ranef = list(scientific.name = phylo_cor),
                                     control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                     chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.univariate.lnCVR.ours.30,
     file=filename)


#################
# 50% reduction #
#################

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_univariate_lnCVR_ours_50_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")


brms.univariate.lnCVR.ours.50 <- brm(lnCVR.sc.ours.50 | se(sqrt(lnCVR.sc.sv)) ~ 1 + 
                                       (1|studyID) + (1|esID) +
                                       (1|scientific.name) + (1|speciesID),
                                     data = stress.data.ours.lnCVR.50,
                                     family = gaussian(),
                                     cov_ranef = list(scientific.name = phylo_cor),
                                     control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                     chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.univariate.lnCVR.ours.50,
     file=filename)



###################
# BIVARIATE MODELS: 
###################

# specifying models' structure
bf.lnRR.ours.10 <- bf(lnRR.sc.ours.10 | se(sqrt(lnRR.sc.sv)) ~
                        1 + (1|p|studyID) + (1|q|esID) + (1|a|scientific.name) + (1|d|speciesID))

bf.lnRR.ours.30 <- bf(lnRR.sc.ours.30 | se(sqrt(lnRR.sc.sv)) ~
                        1 + (1|p|studyID) + (1|q|esID) + (1|a|scientific.name) + (1|d|speciesID))

bf.lnRR.ours.50 <- bf(lnRR.sc.ours.50 | se(sqrt(lnRR.sc.sv)) ~
                        1 + (1|p|studyID) + (1|q|esID) + (1|a|scientific.name) + (1|d|speciesID))

bf.lnVR.ours <- bf(lnVR.sc.ours | se(sqrt(lnVR.sc.sv)) ~
                     1 + (1|p|studyID) + (1|q|esID) + (1|a|scientific.name) + (1|d|speciesID))

bf.lnCVR.ours.10 <- bf(lnCVR.sc.ours.10 | se(sqrt(lnCVR.sc.sv)) ~
                         1 + (1|p|studyID) + (1|q|esID) + (1|a|scientific.name) + (1|d|speciesID))

bf.lnCVR.ours.30 <- bf(lnCVR.sc.ours.30 | se(sqrt(lnCVR.sc.sv)) ~
                         1 + (1|p|studyID) + (1|q|esID) + (1|a|scientific.name) + (1|d|speciesID))

bf.lnCVR.ours.50 <- bf(lnCVR.sc.ours.50 | se(sqrt(lnCVR.sc.sv)) ~
                         1 + (1|p|studyID) + (1|q|esID) + (1|a|scientific.name) + (1|d|speciesID))


###################
# cbind(lnRR,lnVR) 
###################

# sharedcontrol - 10%: 0.999, 20, 6000, 3000, 2:  h (corei7)
# sharedcontrol - 30%: 0.999, 20, 6000, 3000, 2:  h (corei7)
# sharedcontrol - 50%: 0.999, 20, 6000, 3000, 2:  h (corei7)

#################
# 10% reduction #
#################

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_bivariate_lnRR_lnVR_ours_10_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")


brms.bivariate.lnRR.lnVR.ours.10 <- brm(bf.lnRR.ours.10 + bf.lnVR.ours,
                                        data = stress.data.ours,
                                        cov_ranef = list(scientific.name = phylo_cor),
                                        family = gaussian(),
                                        control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                        chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.bivariate.lnRR.lnVR.ours.10,
     file=filename)


#################
# 30% reduction #
#################

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_bivariate_lnRR_lnVR_ours_30_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")


brms.bivariate.lnRR.lnVR.ours.30 <- brm(bf.lnRR.ours.30 + bf.lnVR.ours,
                                        data = stress.data.ours,
                                        cov_ranef = list(scientific.name = phylo_cor),
                                        family = gaussian(),
                                        control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                        chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.bivariate.lnRR.lnVR.ours.30,
     file=filename)


#################
# 50% reduction #
#################

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_bivariate_lnRR_lnVR_ours_50_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")


brms.bivariate.lnRR.lnVR.ours.50 <- brm(bf.lnRR.ours.50 + bf.lnVR.ours,
                                        data = stress.data.ours,
                                        cov_ranef = list(scientific.name = phylo_cor),
                                        family = gaussian(),
                                        control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                        chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.bivariate.lnRR.lnVR.ours.50,
     file=filename)


###################
# cbind(lnVR,lnCVR) 
###################

# sharedcontrol - 10%: 0.999, 20, 6000, 3000, 2: 5 h (corei7)
# sharedcontrol - 30%: 0.999, 20, 6000, 3000, 2: 5 h (corei7)
# sharedcontrol - 50%: 0.999, 20, 6000, 3000, 2: 5 h (corei7)

#################
# 10% reduction #
#################

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_bivariate_lnVR_lnCVR_ours_10_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")


brms.bivariate.lnVR.lnCVR.ours.10 <- brm(bf.lnVR.ours + bf.lnCVR.ours.10,
                                         data = stress.data.ours,
                                         cov_ranef = list(scientific.name = phylo_cor),
                                         family = gaussian(),
                                         control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                         chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.bivariate.lnVR.lnCVR.ours.10,
     file=filename)


#################
# 30% reduction #
#################

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_bivariate_lnVR_lnCVR_ours_30_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")


brms.bivariate.lnVR.lnCVR.ours.30 <- brm(bf.lnVR.ours + bf.lnCVR.ours.30,
                                         data = stress.data.ours,
                                         cov_ranef = list(scientific.name = phylo_cor),
                                         family = gaussian(),
                                         control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                         chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.bivariate.lnVR.lnCVR.ours.30,
     file=filename)


#################
# 50% reduction #
#################

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_bivariate_lnVR_lnCVR_ours_50_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")


brms.bivariate.lnVR.lnCVR.ours.50 <- brm(bf.lnVR.ours + bf.lnCVR.ours.50,
                                         data = stress.data.ours,
                                         cov_ranef = list(scientific.name = phylo_cor),
                                         family = gaussian(),
                                         control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                         chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.bivariate.lnVR.lnCVR.ours.50,
     file=filename)


####################################################################################
# saving session information with all packages versions for reproducibility purposes
sink("models/brms/brms_meta-analysis_reductions_R_session.txt")
sessionInfo()
sink()