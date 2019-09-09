##############################################################
# Authors: 
# Alfredo Sanchez-Tojar (@ASanchez_Tojar)
# Profile: https://goo.gl/PmpPEB
# Department of Evolutionary Biology, Bielefeld University (GER) 
# Email: alfredo.tojar@gmail.com

# Script first created on the 20th of May 2019

##############################################################
# Description of script and instructions
##############################################################

# This script is to analyze the data collected in:

# Eyck et al. 2019: Effects of developmental stress on animal
# phenotype and performance: a quantitative review

# We use the R package 'brms' for the analyses.

##############################################################
# Packages needed
##############################################################

pacman::p_load(openxlsx,brms,tidybayes,ggplot2,plotly,stringr,dplyr,tidyr)


# Clear memory
rm(list=ls())


##############################################################
# Functions needed
##############################################################

# none

##############################################################
# Settings
##############################################################

# model specifications
adapt_delta_value <- 0.9999
max_treedepth_value <- 20
iterations <- 6000
burnin <- 3000
thinning <- 2

##############################################################
# Importing dataset
##############################################################

# database with the corrected data
stress.data <- read.xlsx("data_re-extraction/clean_data/EyckDev_stress_clean_effect_sizes_sp_corrected.xlsx",
                         colNames=T,sheet = 1)

# subsets needed
stress.data.lnRR.ours <- stress.data[!(is.na(stress.data$lnRR.sc.ours)),]
stress.data.lnCVR <- stress.data[!(is.na(stress.data$lnCVR.sc.sv)),]
stress.data.SMDH.ours <- stress.data[!(is.na(stress.data$SMDH.sc.ours)),]
stress.data.HE.cohens.biased <- stress.data[!(is.na(stress.data$cohens.biased.HE)),]

# loading phylogenetic matrix "phylo_cor"
load("data_re-extraction/clean_data/phylo_cor.Rdata") #phylo_cor

# adding the corrected trait classes for the meta-regressions
trait.database <- read.xlsx("data_re-extraction/clean_data/EyckDev_stress_clean_effect_sizes_sp_corrected_trait_modification.xlsx",
                            colNames=T,sheet = 1)

trait.database.red <- trait.database[,c("esID","trait.class.2","potential.alternative","agreement")]
stress.data.metareg <- merge(stress.data,trait.database.red,by="esID",all.x=T)

# subsets needed
stress.data.metareg.lnRR.ours <- stress.data.metareg[!(is.na(stress.data.metareg$lnRR.sc.ours)),]
stress.data.metareg.lnCVR <- stress.data.metareg[!(is.na(stress.data.metareg$lnCVR.sc)),]

########################
# Some summary numbers #
########################

# reducing dataset to lnRR and lnCVR only
stress.data.red <- stress.data[!(is.na(stress.data$lnRR.sc.ours)) | !(is.na(stress.data$lnCVR.sc)),]

# counting number of studies, species, etc...
nrow(stress.data.red)
length(unique(stress.data.red$studyID))
length(unique(stress.data.red$speciesID))

stress.data.red %>% 
  group_by(TAXA) %>% 
  summarise(count = n_distinct(speciesID))

# percentage of effect sizes affected by shared control effects
nrow(stress.data.red[stress.data.red$num.shared.control>1,])
length(unique(stress.data.red[stress.data.red$num.shared.control>1,"studyID"]))
round((length(unique(stress.data.red[stress.data.red$num.shared.control>1,"studyID"]))/length(unique(stress.data.red$studyID)))*100,1)
round((nrow(stress.data.red[stress.data.red$num.shared.control>1,])/nrow(stress.data.red))*100,2)

# counts for traits
stress.data.metareg.red <- stress.data.metareg[!(is.na(stress.data.metareg$lnRR.sc.ours)) | !(is.na(stress.data.metareg$lnCVR.sc)),]
stress.data.metareg.red %>% 
  group_by(trait.class.2) %>% 
  summarise(count = n_distinct(esID))

stress.data.metareg.red %>% 
  group_by(trait.class.2) %>% 
  summarise(count = n_distinct(studyID))

##############################################################
# -------------------------- BRMS -------------------------- #
##############################################################

# We ran univariate and bivariate multilevel meta-analyses 
# in brms.

# Using brms, we ran univariate models based on the effect 
# sizes we calculated ourselves: SMDH, lnRR, lnVR and lnCVR.

# Then we ran the following bivariate model: c(lnRR,lnVR)

# The main models were run based on the effect sizes that
# account for shared control non-independence (see script 004).

# some information obtained through reading about brms:
# ### The default prior is an improper flat prior over the reals. https://cran.r-project.org/web/packages/brms/vignettes/brms_overview.pdf
# ### If Rhat is considerably greater than 1 (i.e., > 1.1), the chains have not yet converged and it is necessary to run more iterations and/or set stronger priors. https://cran.r-project.org/web/packages/brms/vignettes/brms_overview.pdf
# ### weakly non-informative priors = student t-distribution:  student_t(3, 0, 10). We also tested that using a stronger prior (e.g. standard normal distribution: normal(0, 1)) has negligible effects on the model results. https://justincally.github.io/SexualSelection/#

# Note: after running the models, we realized that, the high
# rate of effect size sign inversion that we had to do (21%)
# erase/disrupted the mean-var correlation originally observed
# at the raw level. This had two consequences: (1) running a 
# bivariate meta-analysis would not provide any additional information
# compare to the univariate model because the correlation 
# between lnRR and lnVR was essentially zero once the sign 
# inversions were applied, and thus, borrowing of strength did 
# not take place; (2) interpreting the results of lnVR by itself
# was very difficult as one cannot tell whether the effect comes
# from the mean-var relationship or an independent effect on 
# total variance. Due to this, we decided to not present the 
# results based on lnVR to avoid any confussion however, we
# provide the models at the OSF in case any one is interested
# OSF LINK: https://osf.io/yjua8/ 

###########################
# UNIVARIATE MODELS: lnRR #
###########################

# varcovar: 0.99, 15, 6000, 3000, 2: after ca. 20 h, not a single chain had reached 10% (corei7)
# sharedcontrol: 0.99, 15, 4000, 2000, 2: 47 min (corei7)
# sharedcontrol: 0.999, 20, 6000, 3000, 2: 82 min (corei7)
# sharedcontrol: 0.9999, 20, 6000, 3000, 2: 111 min (corei7)

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_univariate_lnRR_ours_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")

brms.univariate.lnRR.ours <- brm(lnRR.sc.ours | se(sqrt(lnRR.sc.sv)) ~ 1 + 
                                   (1|studyID) + (1|esID) +
                                   (1|scientific.name) + (1|speciesID),
                                 data = stress.data.lnRR.ours,
                                 family = gaussian(),
                                 cov_ranef = list(scientific.name = phylo_cor),
                                 #autocor = cor_fixed(varcovar.studyID.lnRR.ours_0.5), # fixed covariance matrix of the response variable for instance to model multivariate effect sizes in meta-analysis (https://rdrr.io/cran/brms/man/cor_fixed.html)
                                 control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                 chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.univariate.lnRR.ours,
     file=filename)


############################
# UNIVARIATE MODELS: lnCVR #
############################

# sharedcontrol: 0.9999, 20, 6000, 3000, 2: 6 min (corei7)

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_univariate_lnCVR_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")

brms.univariate.lnCVR <- brm(lnCVR.sc | se(sqrt(lnCVR.sc.sv)) ~ 1 + 
                               (1|studyID) + (1|esID) +
                               (1|scientific.name) + (1|speciesID),
                             data = stress.data.lnCVR,
                             family = gaussian(),
                             cov_ranef = list(scientific.name = phylo_cor),
                             control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                             chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.univariate.lnCVR,
     file=filename)


###########################
# UNIVARIATE MODELS: SMDH #
###########################

# sharedcontrol: 0.999, 20, 6000, 3000, 2: 13 min (corei7)
# sharedcontrol: 0.9999, 20, 6000, 3000, 2: 14 min (corei7)

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_univariate_SMDH_ours_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")

brms.univariate.SMDH.ours <- brm(SMDH.sc.ours | se(sqrt(SMDH.sc.sv)) ~ 1 + 
                                   (1|studyID) + (1|esID) +
                                   (1|scientific.name) + (1|speciesID),
                                 data = stress.data.SMDH.ours,
                                 family = gaussian(),
                                 cov_ranef = list(scientific.name = phylo_cor),
                                 control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                 chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.univariate.SMDH.ours,
     file=filename)


########################################
# UNIVARIATE MODELS: Cohen's biased HE #
########################################

# This model is run purely to communicate to H.Eyck what changed
# after we corrected the typos found in the dataset

# sharedcontrol: 0.999, 20, 6000, 3000, 2: 9 min (corei7)
# sharedcontrol: 0.9999, 20, 6000, 3000, 2: 9 min (corei7)

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_univariate_cohens_biased_HE_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")

brms.univariate.cohens.biased.HE <- brm(cohens.biased.HE | se(sqrt(cohens.biased.sv)) ~ 1 + 
                                          (1|studyID) + (1|esID) +
                                          (1|speciesID),
                                        data = stress.data.HE.cohens.biased,
                                        family = gaussian(),
                                        control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                        chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.univariate.cohens.biased.HE,
     file=filename)


##############################
# UNIVARIATE META-REGRESSIONS: 
##############################

# #####################
# # TRAIT CLASS: lnRR #
# #####################
# 
# # sharedcontrol: 0.9999, 20, 6000, 3000, 2: 109 min (corei7)
# 
# ptm <- proc.time() # checking the time needed to run the model
# 
# # filename for saving the model, this avoids having to change the
# # text every time
# filename <- paste0("models/brms/brms_univariate_lnRR_ours_trait_",
#                    "sharedcontrol_",
#                    iterations,"iter_",
#                    burnin,"burnin_",
#                    thinning,"thin_",
#                    adapt_delta_value,"delta_",
#                    max_treedepth_value,"treedepth.RData")
# 
# brms.univariate.lnRR.ours.trait <- brm(lnRR.sc.ours | se(sqrt(lnRR.sc.sv)) ~ 
#                                          1 + trait.class.2 +
#                                          (1|studyID) + (1|esID) +
#                                          (1|scientific.name) + (1|speciesID),
#                                        data = stress.data.metareg.lnRR.ours,
#                                        family = gaussian(),
#                                        cov_ranef = list(scientific.name = phylo_cor),
#                                        control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
#                                        chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)
# 
# proc.time() - ptm # checking the time needed to run the model
# 
# save(brms.univariate.lnRR.ours.trait,
#      file=filename)


#####################################
# TRAIT CLASS: lnRR -> NO INTERCEPT #
#####################################

# sharedcontrol: 0.9999, 20, 6000, 3000, 2:  180 min (corei7)

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_univariate_lnRR_ours_trait_no_intercept_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")

brms.univariate.lnRR.ours.trait.no.intercept <- brm(lnRR.sc.ours | se(sqrt(lnRR.sc.sv)) ~ 
                                                      -1 + trait.class.2 +
                                                      (1|studyID) + (1|esID) +
                                                      (1|scientific.name) + (1|speciesID),
                                                    data = stress.data.metareg.lnRR.ours,
                                                    family = gaussian(),
                                                    cov_ranef = list(scientific.name = phylo_cor),
                                                    control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                                    chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.univariate.lnRR.ours.trait.no.intercept,
     file=filename)


# ######################
# # TRAIT CLASS: lnCVR #
# ######################
# 
# # sharedcontrol: 0.9999, 20, 6000, 3000, 2: 8 min (corei7)
# 
# ptm <- proc.time() # checking the time needed to run the model
# 
# # filename for saving the model, this avoids having to change the
# # text every time
# filename <- paste0("models/brms/brms_univariate_lnCVR_trait_",
#                    "sharedcontrol_",
#                    iterations,"iter_",
#                    burnin,"burnin_",
#                    thinning,"thin_",
#                    adapt_delta_value,"delta_",
#                    max_treedepth_value,"treedepth.RData")
# 
# brms.univariate.lnCVR.trait <- brm(lnCVR.sc | se(sqrt(lnCVR.sc.sv)) ~ 
#                                      1 + trait.class.2 +
#                                      (1|studyID) + (1|esID) +
#                                      (1|scientific.name) + (1|speciesID),
#                                    data = stress.data.metareg.lnCVR,
#                                    family = gaussian(),
#                                    cov_ranef = list(scientific.name = phylo_cor),
#                                    control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
#                                    chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)
# 
# proc.time() - ptm # checking the time needed to run the model
# 
# save(brms.univariate.lnCVR.trait,
#      file=filename)


######################################
# TRAIT CLASS: lnCVR -> NO INTERCEPT #
######################################

# sharedcontrol: 0.9999, 20, 6000, 3000, 2: 7 min (corei7)

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_univariate_lnCVR_trait_no_intercept_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")

brms.univariate.lnCVR.trait.no.intercept <- brm(lnCVR.sc | se(sqrt(lnCVR.sc.sv)) ~ 
                                                  -1 + trait.class.2 +
                                                  (1|studyID) + (1|esID) +
                                                  (1|scientific.name) + (1|speciesID),
                                                data = stress.data.metareg.lnCVR,
                                                family = gaussian(),
                                                cov_ranef = list(scientific.name = phylo_cor),
                                                control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                                chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.univariate.lnCVR.trait.no.intercept,
     file=filename)


#####################
# PUBLICATION BIAS: #
#####################

load("models/brms/brms_univariate_lnRR_ours_sharedcontrol_6000iter_3000burnin_2thin_0.9999delta_20treedepth.RData")
load("models/brms/brms_univariate_SMDH_ours_sharedcontrol_6000iter_3000burnin_2thin_0.9999delta_20treedepth.RData")

#######################
# EGGER'S REGRESSIONS #
#######################

# Following Nakagawa and Santos 2012

########################
# lnRR (meta-analysis) #
########################

# model predictions
stress.data.lnRR.ours$prediction<-predict(brms.univariate.lnRR.ours)[,1]

# precision
stress.data.lnRR.ours$precision<-sqrt(1/stress.data.lnRR.ours$lnRR.sc.sv)

# meta-analytic residuals
stress.data.lnRR.ours$MAR<-stress.data.lnRR.ours$lnRR.sc.ours-stress.data.lnRR.ours$prediction # meta-analytic residual!

# mar adjusted by precision
stress.data.lnRR.ours$zMAR<-stress.data.lnRR.ours$MAR*stress.data.lnRR.ours$precision


ptm <- proc.time() # checking the time needed to run the model


# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_Egger_lnRR_ours_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")

brms.Egger.lnRR.ours <- brm(zMAR ~ precision,
                            data = stress.data.lnRR.ours,
                            family = gaussian(),
                            control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                            chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.Egger.lnRR.ours,
     file=filename)


########################
# SMDH (meta-analysis) #
########################

# model predictions
stress.data.SMDH.ours$prediction<-predict(brms.univariate.SMDH.ours)[,1]

# precision
stress.data.SMDH.ours$precision<-sqrt(1/stress.data.SMDH.ours$SMD.sc.sv)

# meta-analytic residuals
stress.data.SMDH.ours$MAR<-stress.data.SMDH.ours$SMDH.sc.ours-stress.data.SMDH.ours$prediction # meta-analytic residual!

# mar adjusted by precision
stress.data.SMDH.ours$zMAR<-stress.data.SMDH.ours$MAR*stress.data.SMDH.ours$precision


ptm <- proc.time() # checking the time needed to run the model


# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_Egger_SMDH_ours_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")

brms.Egger.SMDH.ours <- brm(zMAR ~ precision,
                            data = stress.data.SMDH.ours,
                            family = gaussian(),
                            control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                            chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.Egger.SMDH.ours,
     file=filename)


#################
# TIME-LAG BIAS #
#################

########
# lnRR #
########

# z-transforming year
stress.data.lnRR.ours$year.z <- scale(stress.data.lnRR.ours$year)

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_univariate_lnRR_ours_year_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")

brms.univariate.lnRR.ours.year <- brm(lnRR.sc.ours | se(sqrt(lnRR.sc.sv)) ~ 
                                        1 + year.z +
                                        (1|studyID) + (1|esID) +
                                        (1|scientific.name) + (1|speciesID),
                                      data = stress.data.lnRR.ours,
                                      family = gaussian(),
                                      cov_ranef = list(scientific.name = phylo_cor),
                                      control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                      chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.univariate.lnRR.ours.year,
     file=filename)


########
# SMDH #
########

# sharedcontrol: 0.9999, 20, 6000, 3000, 2:  17 min (corei7)

# z-transforming year
stress.data.SMDH.ours$year.z <- scale(stress.data.SMDH.ours$year)

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_univariate_SMDH_ours_year_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")

brms.univariate.SMDH.ours.year <- brm(SMDH.sc.ours | se(sqrt(SMDH.sc.sv)) ~ 
                                        1 + year.z +
                                        (1|studyID) + (1|esID) +
                                        (1|scientific.name) + (1|speciesID),
                                      data = stress.data.SMDH.ours,
                                      family = gaussian(),
                                      cov_ranef = list(scientific.name = phylo_cor),
                                      control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                      chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.univariate.SMDH.ours.year,
     file=filename)


####################################################################################
# saving session information with all packages versions for reproducibility purposes
sink("models/brms/brms_meta-analysis_R_session.txt")
sessionInfo()
sink()
