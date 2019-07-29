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

# This script is to re-analyze the data collected in:

# Eyck et al. 2019: Effects of developmental stress on animal
# phenotype and performance: a quantitative review

# We use the R package 'brms' for the analyses

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
# Settings
##############################################################

# for the ggplot2 plots
tm <- theme(panel.background = element_blank(),
            axis.line = element_line(size = 0.75),
            axis.text = element_text(size = 10, colour = "black"),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 16, hjust = 0.5))

# model specifications
adapt_delta_value <- 0.9999
max_treedepth_value <- 20
iterations <- 6000
burnin <- 3000
thinning <- 2

##############################################################
# Importing dataset
##############################################################

# database with the corrected data from our pilot re-extraction
stress.data <- read.xlsx("data_re-extraction/clean_data/EyckDev_stress_clean_effect_sizes_sp_corrected.xlsx",
                         colNames=T,sheet = 1)

# counting number of studies, species, etc...
nrow(stress.data)
length(unique(stress.data$studyID))
length(unique(stress.data$speciesID))

stress.data %>% 
  group_by(TAXA) %>% 
  summarise(count = n_distinct(speciesID))

# percentage of effect sizes affected by shared control effects
nrow(stress.data[stress.data$num.shared.control>1,])
length(unique(stress.data[stress.data$num.shared.control>1,"studyID"]))
round((length(unique(stress.data[stress.data$num.shared.control>1,"studyID"]))/length(unique(stress.data$studyID)))*100,1)
round((nrow(stress.data[stress.data$num.shared.control>1,])/nrow(stress.data))*100,2)

# loading phylogenetic matrix "phylo_cor"
load("data_re-extraction/clean_data/phylo_cor.Rdata") #phylo_cor

# # loading variance-covariance matrices
# load("data_re-extraction/clean_data/varcovar_studyID_SMDH_ours_0.50.Rdata")#varcovar.studyID.SMDH.ours_0.5
# load("data_re-extraction/clean_data/varcovar_studyID_lnRR_ours_0.50.Rdata")#varcovar.studyID.lnRR.ours_0.5
# load("data_re-extraction/clean_data/varcovar_studyID_lnVR_ours_0.50.Rdata")#varcovar.studyID.lnVR.ours_0.5
# load("data_re-extraction/clean_data/varcovar_studyID_lnCVR_ours_0.50.Rdata")#varcovar.studyID.lnCVR.ours_0.5


##############################################################
# mean-variance relationship
##############################################################

# Before running the analyses, is there a mean-variance
# relationship in our dataset?

# reducing dataset to ratio scale data
non.negative <- stress.data[stress.data$mean.control>0 &
                              stress.data$mean.treat>0,]

# all means in dataset
all.ln.means <- log(c(non.negative$mean.control,non.negative$mean.treat))
all.ln.SDs <- log(c(non.negative$SD.control,non.negative$SD.treat))

mv.db <- data.frame(means=as.numeric(all.ln.means),
                    sds=as.numeric(all.ln.SDs),
                    origin=as.factor(c(rep("control",length(non.negative$mean.control)),
                                       rep("treatment",length(non.negative$mean.treat)))))

# they are highly correlated: r = 0.89
pearsons <- round(cor(mv.db$means,mv.db$sds),2)

tiff("plots/mean-variance_relationship.tiff",
     height=15, width=15,
     units='cm', compression="lzw", res=800)

ggplot(mv.db, aes(x = means, y = sds, colour = origin)) +
  geom_point() + tm + labs(x = "ln(mean)", y = "ln(SD)") +
  scale_color_manual(values=c(rgb(0,0,0,alpha=0.75),
                              rgb(224/255,224/255,224/255,alpha=0.8))) +
  #geom_text(x=-5, y=10, label=paste0("r = ",pearsons)) +
  annotate(geom="text", x=-3.5, y=10.5, label=paste0("r = ",pearsons),size=5)+
  theme(legend.position = c(0.85, 0.15),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.key = element_rect(fill = alpha("white", 0.0)))

dev.off()

##############################################################
# -------------------------- BRMS -------------------------- #
##############################################################

# We run univariate and bivariate multilevel meta-analyses 
# in brms.

# Using brms, we will run univariate models based on the effect 
# sizes we calculated ourselves: SMDH, lnRR, lnVR and lnCVR.

# Then we will run bivariate models as: c(lnRR,lnVR),
# and c(lnVR,lnCVR)

# The main models will be run based on the effect sizes that
# account for shared control non-independence (see script 004).

# some information obtained through reading about brms:
# ### The default prior is an improper flat prior over the reals. https://cran.r-project.org/web/packages/brms/vignettes/brms_overview.pdf
# ### If Rhat is considerably greater than 1 (i.e., > 1.1), the chains have not yet converged and it is necessary to run more iterations and/or set stronger priors. https://cran.r-project.org/web/packages/brms/vignettes/brms_overview.pdf
# ### weakly non-informative priors = student t-distribution:  student_t(3, 0, 10). We also tested that using a stronger prior (e.g. standard normal distribution: normal(0, 1)) has negligible effects on the model results. https://justincally.github.io/SexualSelection/#


###########################
# UNIVARIATE MODELS: lnRR #
###########################

# varcovar: 0.99, 15, 6000, 3000, 2: after ca. 20 h, not a single chain had reached 10% (corei7)
# sharedcontrol: 0.99, 15, 4000, 2000, 2: 47 min (corei7)
# sharedcontrol: 0.999, 20, 6000, 3000, 2: 82 min (corei7)
# sharedcontrol: 0.9999, 20, 6000, 3000, 2: 111 min (corei7)

# subset of data needed
stress.data.lnRR.ours <- stress.data[!(is.na(stress.data$lnRR.sc.ours)),]

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_univariate_lnRR_ours_",
                   #"0.5varcov_",
                   #"novarcovar_",
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


###########################
# UNIVARIATE MODELS: lnVR #
###########################

# sharedcontrol: 0.999, 20, 6000, 3000, 2: 6 min (corei7)
# sharedcontrol: 0.9999, 20, 6000, 3000, 2: 7 min (corei7)

# subset of data needed
stress.data.lnVR <- stress.data[!(is.na(stress.data$lnVR.sc.sv)),]

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_univariate_lnVR_",
                   #"0.5varcov_",
                   #"novarcovar_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")


brms.univariate.lnVR <- brm(lnVR.sc | se(sqrt(lnVR.sc.sv)) ~ 1 + 
                              (1|studyID) + (1|esID) +
                              (1|scientific.name) + (1|speciesID),
                            data = stress.data.lnVR,
                            family = gaussian(),
                            cov_ranef = list(scientific.name = phylo_cor),
                            #autocor = cor_fixed(varcovar.studyID.lnRR.ours_0.5), # fixed covariance matrix of the response variable for instance to model multivariate effect sizes in meta-analysis (https://rdrr.io/cran/brms/man/cor_fixed.html)
                            control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                            chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.univariate.lnVR,
     file=filename)


###########################
# UNIVARIATE MODELS: lnCVR #
###########################

# sharedcontrol: 0.9999, 20, 6000, 3000, 2: 6 min (corei7)

# subset of data needed
stress.data.lnCVR <- stress.data[!(is.na(stress.data$lnCVR.sc.sv)),]

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_univariate_lnCVR_",
                   #"0.5varcov_",
                   #"novarcovar_",
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
                             #autocor = cor_fixed(varcovar.studyID.lnRR.ours_0.5), # fixed covariance matrix of the response variable for instance to model multivariate effect sizes in meta-analysis (https://rdrr.io/cran/brms/man/cor_fixed.html)
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

# subset of data needed
stress.data.SMDH <- stress.data[!(is.na(stress.data$SMDH.sc.sv)),]

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_univariate_SMDH_ours_",
                   #"0.5varcov_",
                   #"novarcovar_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")


brms.univariate.SMDH.ours <- brm(SMDH.sc.ours | se(sqrt(SMDH.sc.sv)) ~ 1 + 
                                   (1|studyID) + (1|esID) +
                                   (1|scientific.name) + (1|speciesID),
                                 data = stress.data.SMDH,
                                 family = gaussian(),
                                 cov_ranef = list(scientific.name = phylo_cor),
                                 #autocor = cor_fixed(varcovar.studyID.lnRR.ours_0.5), # fixed covariance matrix of the response variable for instance to model multivariate effect sizes in meta-analysis (https://rdrr.io/cran/brms/man/cor_fixed.html)
                                 control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                 chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.univariate.SMDH.ours,
     file=filename)


########################################
# UNIVARIATE MODELS: Cohen's biased HE #
########################################

# This model is run purely to communicate to H.Eyck what changed
# after we corrected his typos

# sharedcontrol: 0.999, 20, 6000, 3000, 2: 9 min (corei7)
# sharedcontrol: 0.9999, 20, 6000, 3000, 2: 9 min (corei7)

# subset of data needed
stress.data.HE.cohens.biased <- stress.data[!(is.na(stress.data$cohens.biased.HE)),]

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_univariate_cohens_biased_HE_",
                   #"0.5varcov_",
                   #"novarcovar_",
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


###################
# BIVARIATE MODELS: 
###################

# when including a varcovariance matrix, syntax has to be adjusted. See the following links
# https://discourse.mc-stan.org/t/brms-multivariate-meta-analysis-syntax/6801
# https://groups.google.com/forum/#!topic/brms-users/KPQaLs-PU4s

# specifying models' structure
bf.lnRR.ours <- bf(lnRR.sc.ours | se(sqrt(lnRR.sc.sv)) ~
                     1 + (1|p|studyID) + (1|q|esID) + (1|a|scientific.name) + (1|d|speciesID))

bf.lnVR <- bf(lnVR.sc | se(sqrt(lnVR.sc.sv)) ~
                1 + (1|p|studyID) + (1|q|esID) + (1|a|scientific.name) + (1|d|speciesID))

bf.lnCVR <- bf(lnCVR.sc | se(sqrt(lnCVR.sc.sv)) ~
                 1 + (1|p|studyID) + (1|q|esID) + (1|a|scientific.name) + (1|d|speciesID))

# By writing |p| and |q| in between we indicate that all varying
# effects of Study and Index should be modeled as correlated.
# This makes sense since we actually have two model parts, one
# for lnRR and one for lnVR

###################
# cbind(lnRR,lnVR) 
###################
#Error: Fixed residual covariance matrices are not implemented when 'rescor' is estimated.

# sharedcontrol: 0.999, 20, 6000, 3000, 2: 33 h (corei7)
# sharedcontrol: 0.9999, 20, 6000, 3000, 2: 32 h (corei7)

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_bivariate_lnRR_ours_lnVR_",
                   #"0.5varcov_",
                   #"novarcovar_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")


brms.bivariate.lnRR.ours.lnVR <- brm(bf.lnRR.ours + bf.lnVR,
                                     data = stress.data,
                                     cov_ranef = list(scientific.name = phylo_cor),
                                     family = gaussian(),
                                     control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                     chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.bivariate.lnRR.ours.lnVR,
     file=filename)


###################
# cbind(lnVR,lnCVR) 
###################
#Error: Fixed residual covariance matrices are not implemented when 'rescor' is estimated.

# sharedcontrol: 0.999, 20, 6000, 3000, 2: 5 h (corei7)
# sharedcontrol: 0.9999, 20, 6000, 3000, 2: 161 min (corei7)

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_bivariate_lnVR_lnCVR_",
                   #"0.5varcov_",
                   #"novarcovar_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")


brms.bivariate.lnVR.lnCVR <- brm(bf.lnVR + bf.lnCVR,
                                 data = stress.data,
                                 cov_ranef = list(scientific.name = phylo_cor),
                                 family = gaussian(),
                                 control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                 chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.bivariate.lnVR.lnCVR,
     file=filename)


###################
# cbind(lnRR,lnCVR) 
###################

# sharedcontrol: 0.9999, 20, 6000, 3000, 2: 31 h (corei7)

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_bivariate_lnRR_ours_lnCVR_",
                   #"0.5varcov_",
                   #"novarcovar_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")


brms.bivariate.lnRR.ours.lnCVR <- brm(bf.lnRR.ours + bf.lnCVR,
                                      data = stress.data,
                                      cov_ranef = list(scientific.name = phylo_cor),
                                      family = gaussian(),
                                      control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                      chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.bivariate.lnRR.ours.lnCVR,
     file=filename)


###################
# TRIVARIATE MODEL: 
###################

# sharedcontrol: 0.9999, 20, 6000, 3000, 2: 61 h (corei7)

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_trivariate_lnRR_ours_lnRR_lnCVR_",
                   #"0.5varcov_",
                   #"novarcovar_",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")


brms.trivariate.lnRR.ours.lnRR.lnCVR <- brm(bf.lnRR.ours + bf.lnVR + bf.lnCVR,
                                            data = stress.data,
                                            cov_ranef = list(scientific.name = phylo_cor),
                                            family = gaussian(),
                                            control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                            chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.trivariate.lnRR.ours.lnRR.lnCVR,
     file=filename)



##############################
# UNIVARIATE META-REGRESSIONS: 
##############################

# adding the corrected trait classes
trait.database <- read.xlsx("data_re-extraction/clean_data/EyckDev_stress_clean_effect_sizes_sp_corrected_trait_modification.xlsx",
                            colNames=T,sheet = 1)

trait.database.red <- trait.database[,c("esID","trait.class.2","potential.alternative","agreement")]

stress.data.metareg <- merge(stress.data,trait.database.red,by="esID",all.x=T)


#####################
# TRAIT CLASS: lnRR #
#####################

# sharedcontrol: 0.9999, 20, 6000, 3000, 2: 109 min (corei7)

# subset of data needed
stress.data.metareg.lnRR.ours <- stress.data.metareg[!(is.na(stress.data.metareg$lnRR.sc.ours)),]

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_univariate_lnRR_ours_trait",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")


brms.univariate.lnRR.ours.trait <- brm(lnRR.sc.ours | se(sqrt(lnRR.sc.sv)) ~ 
                                         1 + trait.class.2 +
                                         (1|studyID) + (1|esID) +
                                         (1|scientific.name) + (1|speciesID),
                                       data = stress.data.metareg.lnRR.ours,
                                       family = gaussian(),
                                       cov_ranef = list(scientific.name = phylo_cor),
                                       control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                       chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.univariate.lnRR.ours.trait,
     file=filename)


#####################
# TRAIT CLASS: lnVR #
#####################

# sharedcontrol: 0.9999, 20, 6000, 3000, 2: 8 min (corei7)

# subset of data needed
stress.data.metareg.lnVR <- stress.data.metareg[!(is.na(stress.data.metareg$lnVR.sc)),]

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_univariate_lnVR_trait",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")


brms.univariate.lnVR.trait <- brm(lnVR.sc | se(sqrt(lnVR.sc.sv)) ~ 
                                    1 + trait.class.2 +
                                    (1|studyID) + (1|esID) +
                                    (1|scientific.name) + (1|speciesID),
                                  data = stress.data.metareg.lnVR,
                                  family = gaussian(),
                                  cov_ranef = list(scientific.name = phylo_cor),
                                  control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                  chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.univariate.lnVR.trait,
     file=filename)


######################
# TRAIT CLASS: lnCVR #
######################

# sharedcontrol: 0.9999, 20, 6000, 3000, 2: 8 min (corei7)

# subset of data needed
stress.data.metareg.lnCVR <- stress.data.metareg[!(is.na(stress.data.metareg$lnCVR.sc)),]

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_univariate_lnCVR_trait",
                   "sharedcontrol_",
                   iterations,"iter_",
                   burnin,"burnin_",
                   thinning,"thin_",
                   adapt_delta_value,"delta_",
                   max_treedepth_value,"treedepth.RData")


brms.univariate.lnCVR.trait <- brm(lnCVR.sc | se(sqrt(lnCVR.sc.sv)) ~ 
                                     1 + trait.class.2 +
                                     (1|studyID) + (1|esID) +
                                     (1|scientific.name) + (1|speciesID),
                                   data = stress.data.metareg.lnCVR,
                                   family = gaussian(),
                                   cov_ranef = list(scientific.name = phylo_cor),
                                   control = list(adapt_delta = adapt_delta_value, max_treedepth = max_treedepth_value),
                                   chains = 4, cores = 4, iter = iterations, warmup = burnin, thin = thinning)

proc.time() - ptm # checking the time needed to run the model

save(brms.univariate.lnCVR.trait,
     file=filename)


################
# TIME-LAG BIAS: 
################

# sharedcontrol: 0.9999, 20, 6000, 3000, 2: 104 min (corei7)

# subset of data needed
stress.data.lnRR.ours <- stress.data[!(is.na(stress.data$lnRR.sc.ours)),]

stress.data.lnRR.ours$year.z <- scale(stress.data.lnRR.ours$year)

ptm <- proc.time() # checking the time needed to run the model

# filename for saving the model, this avoids having to change the
# text every time
filename <- paste0("models/brms/brms_univariate_lnRR_ours_year",
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


####################################################################################
# saving session information with all packages versions for reproducibility purposes
sink("models/brms/brms_meta-analysis_R_session.txt")
sessionInfo()
sink()
