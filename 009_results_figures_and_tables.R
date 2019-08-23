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

pacman::p_load(openxlsx,brms,tidybayes,ggplot2,plotly,stringr,
               dplyr,pander,gt,ggridges,metafor)

#webshot::install_phantomjs()


# Clear memory
rm(list=ls())


##############################################################
# Functions needed
##############################################################

#####################################################################
# Functions obtained from: 
# https://github.com/JustinCally/SexualSelection

# This code is from the Github page of Jared Lander: https://github.com/jaredlander/coefplot/blob/master/R/position.r

# Detect and prevent collisions.
# Powers dodging, stacking and filling.
collidev <- function(data, height = NULL, name, strategy, check.height = TRUE) {
  # Determine height
  if (!is.null(height)) {
    # height set manually
    if (!(all(c("ymin", "ymax") %in% names(data)))) {
      data$ymin <- data$y - height / 2
      data$ymax <- data$y + height / 2
    }
  } else {
    if (!(all(c("ymin", "ymax") %in% names(data)))) {
      data$ymin <- data$y
      data$ymax <- data$y
    }
    
    # height determined from data, must be floating point constant
    heights <- unique(data$ymax - data$ymin)
    heights <- heights[!is.na(heights)]
    height <- heights[1]
  }
  
  # Reorder by x position, relying on stable sort to preserve existing
  # ordering, which may be by group or order.
  data <- data[order(data$ymin), ]
  
  # Check for overlap
  intervals <- as.numeric(t(unique(data[c("ymin", "ymax")])))
  intervals <- intervals[!is.na(intervals)]
  
  if (length(unique(intervals)) > 1 & any(diff(scale(intervals)) < -1e-6)) {
    warning(name, " requires non-overlapping y intervals", call. = FALSE)
    # This is where the algorithm from [L. Wilkinson. Dot plots.
    # The American Statistician, 1999.] should be used
  }
  
  if (!is.null(data$xmax)) {
    plyr::ddply(data, "ymin", strategy, height = height)
  } else if (!is.null(data$x)) {
    data$xmax <- data$x
    data <- plyr::ddply(data, "ymin", strategy, height = height)
    data$x <- data$xmax
    data
  } else {
    stop("Neither x nor xmax defined")
  }
}

# Stack overlapping intervals.
# Assumes that each set has the same horizontal position
pos_stackv <- function(df, height) {
  if (nrow(df) == 1) return(df)
  
  n <- nrow(df) + 1
  x <- ifelse(is.na(df$x), 0, df$x)
  if (all(is.na(df$y))) {
    heights <- rep(NA, n)
  } else {
    heights <- c(0, cumsum(x))
  }
  
  df$xmin <- heights[-n]
  df$xmax <- heights[-1]
  df$x <- df$xmax
  df
}

# Stack overlapping intervals and set height to 1.
# Assumes that each set has the same horizontal position.
pos_fillv <- function(df, height) {
  stacked <- pos_stackv(df, height)
  stacked$xmin <- stacked$xmin / max(stacked$xmax)
  stacked$xmax <- stacked$xmax / max(stacked$xmax)
  stacked$x <- stacked$xmax
  stacked
}

# Dodge overlapping interval.
# Assumes that each set has the same horizontal position.
pos_dodgev <- function(df, height) {
  n <- length(unique(df$group))
  if (n == 1) return(df)
  
  if (!all(c("ymin", "ymax") %in% names(df))) {
    df$ymin <- df$y
    df$ymax <- df$y
  }
  
  d_height <- max(df$ymax - df$ymin)
  
  # df <- data.frame(n = c(2:5, 10, 26), div = c(4, 3, 2.666666,  2.5, 2.2, 2.1))
  # ggplot(df, aes(n, div)) + geom_point()
  
  # Have a new group index from 1 to number of groups.
  # This might be needed if the group numbers in this set don't include all of 1:n
  groupidy <- match(df$group, sort(unique(df$group)))
  
  # Find the center for each group, then use that to calculate xmin and xmax
  df$y <- df$y + height * ((groupidy - 0.5) / n - .5)
  df$ymin <- df$y - d_height / n / 2
  df$ymax <- df$y + d_height / n / 2
  df
}

position_dodgev <- function(height = NULL) {
  ggproto(NULL, PositionDodgeV, height = height)
}

PositionDodgeV <- ggproto(
  "PositionDodgeV", Position,
  required_aes = "y",
  height = NULL,
  setup_params = function(self, data) {
    if (is.null(data$ymin) && is.null(data$ymax) && is.null(self$height)) {
      warning("height not defined. Set with `position_dodgev(height = ?)`",
              call. = FALSE)
    }
    list(height = self$height)
  },
  
  compute_panel = function(data, params, scales) {
    collidev(data, params$height, "position_dodgev", pos_dodgev, check.height = FALSE)
  }
)


##############################################################
# Settings
##############################################################

# none

##############################################################
# Importing dataset
##############################################################

# database with the corrected data from our pilot re-extraction
stress.data <- read.xlsx("data_re-extraction/clean_data/EyckDev_stress_clean_effect_sizes_sp_corrected.xlsx",
                         colNames=T,sheet = 1)

# subsets
stress.data.lnRR.ours <- stress.data[!(is.na(stress.data$lnRR.sc.ours)),]
#stress.data.lnVR <- stress.data[!(is.na(stress.data$lnVR.sc.sv)),]
stress.data.lnCVR <- stress.data[!(is.na(stress.data$lnCVR.sc.sv)),]
stress.data.SMDH.ours <- stress.data[!(is.na(stress.data$SMDH.sc.ours)),]

# adding the corrected trait classes for the meta-regressions
trait.database <- read.xlsx("data_re-extraction/clean_data/EyckDev_stress_clean_effect_sizes_sp_corrected_trait_modification.xlsx",
                            colNames=T,sheet = 1)
trait.database.red <- trait.database[,c("esID","trait.class.2","potential.alternative","agreement")]
stress.data.metareg <- merge(stress.data,trait.database.red,by="esID",all.x=T)

# subsets needed
stress.data.metareg.lnRR.ours <- stress.data.metareg[!(is.na(stress.data.metareg$lnRR.sc.ours)),]
#stress.data.metareg.lnVR <- stress.data.metareg[!(is.na(stress.data.metareg$lnVR.sc)),]
stress.data.metareg.lnCVR <- stress.data.metareg[!(is.na(stress.data.metareg$lnCVR.sc)),]
stress.data.metareg.SMDH.ours <- stress.data.metareg[!(is.na(stress.data.metareg$SMDH.sc.ours)),]

# loading phylogenetic matrix "phylo_cor"
load("data_re-extraction/clean_data/phylo_cor.Rdata") #phylo_cor

##############################################################
# Loading models
##############################################################

# loading all the models at once
temp = list.files(path = "./models/brms/", pattern="*.RData")
for (i in 1:length(temp)) assign(temp[i], load(paste0("models/brms/",temp[i])))

##############################################################
# Did models run well? Spoiler: All seems good!
##############################################################

# ###############
# # Meta-analyses
# 
# # lnRR
# summary(brms.univariate.lnRR.ours) # Rhat = 1.00 at all times, Eff.sample range=(1267,4764)
# plot(brms.univariate.lnRR.ours)
# pairs(brms.univariate.lnRR.ours)
# brms::pp_check(brms.univariate.lnRR.ours, resp = "lnRR.sc.ours") # graphical posterior predictive check
# 
# # lnVR
# summary(brms.univariate.lnVR) # Rhat <= 1.01 at all times (only species ID has an Rhat = 1.01), Eff.sample range=(1255,3701)
# plot(brms.univariate.lnVR)
# pairs(brms.univariate.lnVR)
# brms::pp_check(brms.univariate.lnVR, resp = "lnVR.sc") # graphical posterior predictive check
# 
# # lnCVR
# summary(brms.univariate.lnCVR) # Rhat = 1.00 at all times, Eff.sample range=(1319,4206)
# plot(brms.univariate.lnCVR)
# pairs(brms.univariate.lnCVR)
# brms::pp_check(brms.univariate.lnCVR, resp = "lnCVR.sc") # graphical posterior predictive check
# 
# # SMDH
# summary(brms.univariate.SMDH.ours) # Rhat = 1.00 at all times, Eff.sample range=(719,4782)
# plot(brms.univariate.SMDH.ours)
# pairs(brms.univariate.SMDH.ours)
# brms::pp_check(brms.univariate.SMDH.ours, resp = "SMDH.sc.ours") # graphical posterior predictive check
# 
# # Cohen's biased
# summary(brms.univariate.cohens.biased.HE) # Rhat <= 1.01 at all times, Eff.sample range=(362,2903)
# plot(brms.univariate.cohens.biased.HE)
# pairs(brms.univariate.cohens.biased.HE)
# brms::pp_check(brms.univariate.cohens.biased.HE, resp = "cohens.biased.HE") # graphical posterior predictive check
# 
# # lnRR,lnVR (bivariate)
# summary(brms.bivariate.lnRR.ours.lnVR) # Rhat = 1.00 at all times, Eff.sample range=(606,4232)
# plot(brms.bivariate.lnRR.ours.lnVR)
# brms::pp_check(brms.bivariate.lnRR.ours.lnVR, resp = "lnRRscours") # graphical posterior predictive check
# brms::pp_check(brms.bivariate.lnRR.ours.lnVR, resp = "lnVRsc") # graphical posterior predictive check
# 
# ##################
# # Meta-regressions
# 
# # lnRR
# summary(brms.univariate.lnRR.ours.trait) # Rhat = 1.00 at all times, Eff.sample range=(1175,4823)
# plot(brms.univariate.lnRR.ours.trait)
# #pairs(brms.univariate.lnRR.ours.trait) # too crowdie
# brms::pp_check(brms.univariate.lnRR.ours.trait, resp = "lnRR.sc.ours") # graphical posterior predictive check
# 
# # lnVR
# summary(brms.univariate.lnVR.trait) # Rhat = 1.00 at all times, Eff.sample range=(1398,4256)
# plot(brms.univariate.lnVR.trait)
# brms::pp_check(brms.univariate.lnVR.trait, resp = "lnVR.sc") # graphical posterior predictive check
# 
# # lnCVR
# summary(brms.univariate.lnCVR.trait) # Rhat = 1.00 at all times, Eff.sample range=(1316,4213)
# plot(brms.univariate.lnCVR.trait)
# brms::pp_check(brms.univariate.lnCVR.trait, resp = "lnCVR.sc") # graphical posterior predictive check
# 
# # SMDH
# summary(brms.univariate.SMDH.ours.trait) # Rhat = 1.00 at all times, Eff.sample range=(688,4458)
# plot(brms.univariate.SMDH.ours.trait)
# brms::pp_check(brms.univariate.SMDH.ours.trait, resp = "SMDH.sc.ours") # graphical posterior predictive check
# 
# # Egger: lnRR
# summary(brms.Egger.lnRR.ours) # Rhat = 1.00 at all times, Eff.sample range=(4397,5626)
# plot(brms.Egger.lnRR.ours)
# brms::pp_check(brms.Egger.lnRR.ours, resp = "zMAR") # graphical posterior predictive check
# 
# # Egger: lnRR.trait
# summary(brms.Egger.lnRR.ours.trait) # Rhat = 1.00 at all times, Eff.sample range=(3434,5905)
# plot(brms.Egger.lnRR.ours.trait)
# brms::pp_check(brms.Egger.lnRR.ours.trait, resp = "zMAR") # graphical posterior predictive check
# 
# # Egger: SMDH
# summary(brms.Egger.SMDH.ours) # Rhat = 1.00 at all times, Eff.sample range=(5332,5700)
# plot(brms.Egger.SMDH.ours)
# brms::pp_check(brms.Egger.SMDH.ours, resp = "zMAR") # graphical posterior predictive check
# 
# # Egger: SMDH.trait
# summary(brms.Egger.SMDH.ours.trait) # Rhat = 1.00 at all times, Eff.sample range=(5129,5629)
# plot(brms.Egger.SMDH.ours.trait)
# brms::pp_check(brms.Egger.SMDH.ours.trait, resp = "zMAR") # graphical posterior predictive check
# 
# # Time-lag: lnRR
# summary(brms.univariate.lnRR.ours.year) # Rhat = 1.00 at all times, Eff.sample range=(1193,5390)
# plot(brms.univariate.lnRR.ours.year)
# brms::pp_check(brms.univariate.lnRR.ours.year, resp = "lnRR.sc.ours") # graphical posterior predictive check
# 
# # Time-lag: SMDH
# summary(brms.univariate.SMDH.ours.year) # Rhat = 1.00 at all times, Eff.sample range=(525,3810)
# plot(brms.univariate.SMDH.ours.year)
# brms::pp_check(brms.univariate.SMDH.ours.year, resp = "SMDH.sc.ours") # graphical posterior predictive check


##############################################################
# Heterogeneity
##############################################################

# check the names of the random effects in our models
# get_variables()
# The variables we are interested in from all univariate models
# are: sd_esID__Intercept, sd_scientific.name__Intercept,
# sd_speciesID__Intercept, sd_studyID__Intercept

#####################
# lnRR
#####################

# extracting the posterior distributions from our models
posterior.brms.univariate.lnRR.ours <- posterior_samples(brms.univariate.lnRR.ours)

# WI = weight
WI.lnRR.ours <- na.omit(1/stress.data.lnRR.ours$lnRR.sc.sv)

# s2I = measurement error variance = sigma2m
s2I.lnRR.ours <- sum(WI.lnRR.ours*(length(WI.lnRR.ours)-1))/(sum(WI.lnRR.ours)^2-sum(WI.lnRR.ours^2))

# total variance, including measurement error variance
total_var.lnRR.ours <- posterior.brms.univariate.lnRR.ours$sd_esID__Intercept +
  posterior.brms.univariate.lnRR.ours$sd_scientific.name__Intercept +
  posterior.brms.univariate.lnRR.ours$sd_speciesID__Intercept +
  posterior.brms.univariate.lnRR.ours$sd_studyID__Intercept +
  s2I.lnRR.ours

# total heterogeneity I2
I2_total.lnRR.ours <- (total_var.lnRR.ours-s2I.lnRR.ours)/total_var.lnRR.ours
# round(mean(I2_total.lnRR.ours),3)*100
# round(quantile(I2_total.lnRR.ours, c(0.025, 0.975), na.rm = TRUE),3)*100

# observational level I2
I2_esiD.lnRR.ours <- posterior.brms.univariate.lnRR.ours$sd_esID__Intercept/total_var.lnRR.ours
# round(mean(I2_esiD.lnRR.ours),3)*100
# round(quantile(I2_esiD.lnRR.ours, c(0.025, 0.975), na.rm = TRUE),3)*100

# studyID I2
I2_studyID.lnRR.ours <- posterior.brms.univariate.lnRR.ours$sd_studyID__Intercept/total_var.lnRR.ours
# round(mean(I2_studyID.lnRR.ours),3)*100
# round(quantile(I2_studyID.lnRR.ours, c(0.025, 0.975), na.rm = TRUE),3)*100

# phylogeny I2: notice that s2I is substracted from this calculation as phylogenetic
# relatedness is a "fixed random effect"
I2_phylo.lnRR.ours <- posterior.brms.univariate.lnRR.ours$sd_scientific.name__Intercept/(total_var.lnRR.ours-s2I.lnRR.ours)
# round(mean(I2_speciesID.lnRR.ours),3)*100
# round(quantile(I2_speciesID.lnRR.ours, c(0.025, 0.975), na.rm = TRUE),3)*100

# speciesID I2
I2_speciesID.lnRR.ours <- posterior.brms.univariate.lnRR.ours$sd_speciesID__Intercept/total_var.lnRR.ours
# round(mean(I2_speciesID.lnRR.ours),3)*100
# round(quantile(I2_speciesID.lnRR.ours, c(0.025, 0.975), na.rm = TRUE),3)*100


# #####################
# # lnVR
# #####################
# 
# # extracting the posterior distributions from our models
# posterior.brms.univariate.lnVR <- posterior_samples(brms.univariate.lnVR)
# 
# # WI = weight
# WI.lnVR <- na.omit(1/stress.data.lnVR$lnVR.sc.sv)
# 
# # s2I = measurement error variance = sigma2m
# s2I.lnVR <- sum(WI.lnVR*(length(WI.lnVR)-1))/(sum(WI.lnVR)^2-sum(WI.lnVR^2))
# 
# # total variance, including measurement error variance
# total_var.lnVR <- posterior.brms.univariate.lnVR$sd_esID__Intercept +
#   posterior.brms.univariate.lnVR$sd_scientific.name__Intercept +
#   posterior.brms.univariate.lnVR$sd_speciesID__Intercept +
#   posterior.brms.univariate.lnVR$sd_studyID__Intercept +
#   s2I.lnVR
# 
# # total heterogeneity I2
# I2_total.lnVR <- (total_var.lnVR-s2I.lnVR)/total_var.lnVR
# 
# # observational level I2
# I2_esiD.lnVR <- posterior.brms.univariate.lnVR$sd_esID__Intercept/total_var.lnVR
# 
# # studyID I2
# I2_studyID.lnVR <- posterior.brms.univariate.lnVR$sd_studyID__Intercept/total_var.lnVR
# 
# # phylogeny I2: notice that s2I is substracted from this calculation as phylogenetic
# # relatedness is a "fixed random effect"
# I2_phylo.lnVR <- posterior.brms.univariate.lnVR$sd_scientific.name__Intercept/(total_var.lnVR-s2I.lnVR)
# 
# # speciesID I2
# I2_speciesID.lnVR <- posterior.brms.univariate.lnVR$sd_speciesID__Intercept/total_var.lnVR


#####################
# lnCVR
#####################

# extracting the posterior distributions from our models
posterior.brms.univariate.lnCVR <- posterior_samples(brms.univariate.lnCVR)

# WI = weight
WI.lnCVR <- na.omit(1/stress.data.lnCVR$lnCVR.sc.sv)

# s2I = measurement error variance = sigma2m
s2I.lnCVR <- sum(WI.lnCVR*(length(WI.lnCVR)-1))/(sum(WI.lnCVR)^2-sum(WI.lnCVR^2))

# total variance, including measurement error variance
total_var.lnCVR <- posterior.brms.univariate.lnCVR$sd_esID__Intercept +
  posterior.brms.univariate.lnCVR$sd_scientific.name__Intercept +
  posterior.brms.univariate.lnCVR$sd_speciesID__Intercept +
  posterior.brms.univariate.lnCVR$sd_studyID__Intercept +
  s2I.lnCVR

# total heterogeneity I2
I2_total.lnCVR <- (total_var.lnCVR-s2I.lnCVR)/total_var.lnCVR

# observational level I2
I2_esiD.lnCVR <- posterior.brms.univariate.lnCVR$sd_esID__Intercept/total_var.lnCVR

# studyID I2
I2_studyID.lnCVR <- posterior.brms.univariate.lnCVR$sd_studyID__Intercept/total_var.lnCVR

# phylogeny I2: notice that s2I is substracted from this calculation as phylogenetic
# relatedness is a "fixed random effect"
I2_phylo.lnCVR <- posterior.brms.univariate.lnCVR$sd_scientific.name__Intercept/(total_var.lnCVR-s2I.lnCVR)

# speciesID I2
I2_speciesID.lnCVR <- posterior.brms.univariate.lnCVR$sd_speciesID__Intercept/total_var.lnCVR


#####################
# SMDH
#####################

# extracting the posterior distributions from our models
posterior.brms.univariate.SMDH.ours <- posterior_samples(brms.univariate.SMDH.ours)

# WI = weight
WI.SMDH.ours <- na.omit(1/stress.data.SMDH.ours$SMDH.sc.sv)

# s2I = measurement error variance = sigma2m
s2I.SMDH.ours <- sum(WI.SMDH.ours*(length(WI.SMDH.ours)-1))/(sum(WI.SMDH.ours)^2-sum(WI.SMDH.ours^2))

# total variance, including measurement error variance
total_var.SMDH.ours <- posterior.brms.univariate.SMDH.ours$sd_esID__Intercept +
  posterior.brms.univariate.SMDH.ours$sd_scientific.name__Intercept +
  posterior.brms.univariate.SMDH.ours$sd_speciesID__Intercept +
  posterior.brms.univariate.SMDH.ours$sd_studyID__Intercept +
  s2I.SMDH.ours

# total heterogeneity I2
I2_total.SMDH.ours <- (total_var.SMDH.ours-s2I.SMDH.ours)/total_var.SMDH.ours

# observational level I2
I2_esiD.SMDH.ours <- posterior.brms.univariate.SMDH.ours$sd_esID__Intercept/total_var.SMDH.ours

# studyID I2
I2_studyID.SMDH.ours <- posterior.brms.univariate.SMDH.ours$sd_studyID__Intercept/total_var.SMDH.ours

# phylogeny I2: notice that s2I is substracted from this calculation as phylogenetic
# relatedness is a "fixed random effect"
I2_phylo.SMDH.ours <- posterior.brms.univariate.SMDH.ours$sd_scientific.name__Intercept/(total_var.SMDH.ours-s2I.SMDH.ours)

# speciesID I2
I2_speciesID.SMDH.ours <- posterior.brms.univariate.SMDH.ours$sd_speciesID__Intercept/total_var.SMDH.ours


##############################################################
# R2 marginal
##############################################################

#####################
# lnRR.trait
#####################

# extracting posterior samples from the model
posterior.brms.univariate.lnRR.ours.trait <- posterior_samples(brms.univariate.lnRR.ours.trait)


# building a design matrix for the fixed effects
newdat<-expand.grid(trait.class.2 = stress.data.metareg.lnRR.ours$trait.class.2)
fixeff.design.matrix<-model.matrix(~trait.class.2,data=newdat)


# Estimating R2
vmVarF<-numeric(1000)
for(i in 1:1000){
  Var<-var(as.vector(as.matrix(posterior.brms.univariate.lnRR.ours.trait[i,c(1:6)]) %*% t(fixeff.design.matrix)))
  vmVarF[i]<-Var}

R2m.lnRR.trait<-100*(vmVarF/(vmVarF+posterior.brms.univariate.lnRR.ours.trait$sd_esID__Intercept+
                               posterior.brms.univariate.lnRR.ours.trait$sd_scientific.name__Intercept+
                               posterior.brms.univariate.lnRR.ours.trait$sd_speciesID__Intercept+
                               posterior.brms.univariate.lnRR.ours.trait$sd_studyID__Intercept))

# #####################
# # lnVR.trait
# #####################
# 
# # extracting posterior samples from the model
# posterior.brms.univariate.lnVR.trait <- posterior_samples(brms.univariate.lnVR.trait)
# 
# 
# # building a design matrix for the fixed effects
# newdat<-expand.grid(trait.class.2 = stress.data.metareg.lnVR$trait.class.2)
# fixeff.design.matrix<-model.matrix(~trait.class.2,data=newdat)
# 
# 
# # Estimating R2
# vmVarF<-numeric(1000)
# for(i in 1:1000){
#   Var<-var(as.vector(as.matrix(posterior.brms.univariate.lnVR.trait[i,c(1:6)]) %*% t(fixeff.design.matrix)))
#   vmVarF[i]<-Var}
# 
# R2m.lnVR.trait<-100*(vmVarF/(vmVarF+posterior.brms.univariate.lnVR.trait$sd_esID__Intercept+
#                                posterior.brms.univariate.lnVR.trait$sd_scientific.name__Intercept+
#                                posterior.brms.univariate.lnVR.trait$sd_speciesID__Intercept+
#                                posterior.brms.univariate.lnVR.trait$sd_studyID__Intercept))
# 


#####################
# lnCVR.trait
#####################

# extracting posterior samples from the model
posterior.brms.univariate.lnCVR.trait <- posterior_samples(brms.univariate.lnCVR.trait)


# building a design matrix for the fixed effects
newdat<-expand.grid(trait.class.2 = stress.data.metareg.lnCVR$trait.class.2)
fixeff.design.matrix<-model.matrix(~trait.class.2,data=newdat)


# Estimating R2
vmVarF<-numeric(1000)
for(i in 1:1000){
  Var<-var(as.vector(as.matrix(posterior.brms.univariate.lnCVR.trait[i,c(1:6)]) %*% t(fixeff.design.matrix)))
  vmVarF[i]<-Var}

R2m.lnCVR.trait<-100*(vmVarF/(vmVarF+posterior.brms.univariate.lnCVR.trait$sd_esID__Intercept+
                                posterior.brms.univariate.lnCVR.trait$sd_scientific.name__Intercept+
                                posterior.brms.univariate.lnCVR.trait$sd_speciesID__Intercept+
                                posterior.brms.univariate.lnCVR.trait$sd_studyID__Intercept))


#####################
# SMDH.trait
#####################

# extracting posterior samples from the model
posterior.brms.univariate.SMDH.ours.trait <- posterior_samples(brms.univariate.SMDH.ours.trait)


# building a design matrix for the fixed effects
newdat<-expand.grid(trait.class.2 = stress.data.metareg.SMDH.ours$trait.class.2)
fixeff.design.matrix<-model.matrix(~trait.class.2,data=newdat)


# Estimating R2
vmVarF<-numeric(1000)
for(i in 1:1000){
  Var<-var(as.vector(as.matrix(posterior.brms.univariate.SMDH.ours.trait[i,c(1:6)]) %*% t(fixeff.design.matrix)))
  vmVarF[i]<-Var}

R2m.SMDH.trait<-100*(vmVarF/(vmVarF+posterior.brms.univariate.SMDH.ours.trait$sd_esID__Intercept+
                               posterior.brms.univariate.SMDH.ours.trait$sd_scientific.name__Intercept+
                               posterior.brms.univariate.SMDH.ours.trait$sd_speciesID__Intercept+
                               posterior.brms.univariate.SMDH.ours.trait$sd_studyID__Intercept))


#####################
# lnRR.year
#####################

# extracting posterior samples from the model
posterior.brms.univariate.lnRR.ours.year <- posterior_samples(brms.univariate.lnRR.ours.year)


# building a design matrix for the fixed effects
stress.data.lnRR.ours$year.z <- scale(stress.data.lnRR.ours$year)
newdat<-expand.grid(year.z = stress.data.lnRR.ours$year.z)
fixeff.design.matrix<-model.matrix(~year.z,data=newdat)


# Estimating R2
vmVarF<-numeric(1000)
for(i in 1:1000){
  Var<-var(as.vector(as.matrix(posterior.brms.univariate.lnRR.ours.year[i,c(1:2)]) %*% t(fixeff.design.matrix)))
  vmVarF[i]<-Var}

R2m.lnRR.year<-100*(vmVarF/(vmVarF+posterior.brms.univariate.lnRR.ours.year$sd_esID__Intercept+
                              posterior.brms.univariate.lnRR.ours.year$sd_scientific.name__Intercept+
                              posterior.brms.univariate.lnRR.ours.year$sd_speciesID__Intercept+
                              posterior.brms.univariate.lnRR.ours.year$sd_studyID__Intercept))


#####################
# SMDH.year
#####################

# extracting posterior samples from the model
posterior.brms.univariate.SMDH.ours.year <- posterior_samples(brms.univariate.SMDH.ours.year)


# building a design matrix for the fixed effects
stress.data.SMDH.ours$year.z <- scale(stress.data.SMDH.ours$year)
newdat<-expand.grid(year.z = stress.data.SMDH.ours$year.z)
fixeff.design.matrix<-model.matrix(~year.z,data=newdat)


# Estimating R2
vmVarF<-numeric(1000)
for(i in 1:1000){
  Var<-var(as.vector(as.matrix(posterior.brms.univariate.SMDH.ours.year[i,c(1:2)]) %*% t(fixeff.design.matrix)))
  vmVarF[i]<-Var}

R2m.SMDH.year<-100*(vmVarF/(vmVarF+posterior.brms.univariate.SMDH.ours.year$sd_esID__Intercept+
                              posterior.brms.univariate.SMDH.ours.year$sd_scientific.name__Intercept+
                              posterior.brms.univariate.SMDH.ours.year$sd_speciesID__Intercept+
                              posterior.brms.univariate.SMDH.ours.year$sd_studyID__Intercept))


##############################################################
# Cochrane's Q test
##############################################################

#####################
# lnRR
#####################
lnRR.ours.Q <- rma.mv(yi=lnRR.sc.ours,
                      V = lnRR.sc.sv, 
                      random = list(~ 1 | studyID, ~ 1 | esID, 
                                    ~ 1 | scientific.name, ~ 1 | speciesID), 
                      R = list(scientific.name = phylo_cor),
                      data=stress.data.lnRR.ours)$QE

#####################
# lnCVR
#####################
lnCVR.Q <- rma.mv(yi=lnCVR.sc, 
                  V = lnCVR.sc.sv, 
                  random = list(~ 1 | studyID, ~ 1 | esID, 
                                ~ 1 | scientific.name, ~ 1 | speciesID), 
                  R = list(scientific.name = phylo_cor),
                  data=stress.data.lnCVR)$QE

#####################
# SMDH
#####################
SMDH.ours.Q <- rma.mv(yi=SMDH.sc.ours, 
                      V = SMDH.sc.sv, 
                      random = list(~ 1 | studyID, ~ 1 | esID, 
                                    ~ 1 | scientific.name, ~ 1 | speciesID), 
                      R = list(scientific.name = phylo_cor),
                      data=stress.data.SMDH.ours)$QE


##############################################################
# Creating tables for the results
##############################################################

###########
# TABLE 1
###########

# using the package gt to create fancy tables: https://github.com/rstudio/gt
table.column.names <- c("Effect.size","K","Metaanalytic.mean",
                        "Observational","Study","Species","Phylogeny","Total","Qtest","Eggers")

# effect size names
effect.sizes <- c("lnRR",
                  #"lnVR",
                  "lnCVR",
                  "SMDH")

# sample sizes
ks <- c(nrow(stress.data.lnRR.ours),
        #nrow(stress.data.lnVR),
        nrow(stress.data.lnCVR),
        nrow(stress.data.SMDH.ours))

# point summaries
brms.univariate.lnRR.ours %>%
  spread_draws(b_Intercept) %>%
  mode_hdi() -> point.summaries.univariate.lnRR.ours

# brms.univariate.lnVR %>%
#   spread_draws(b_Intercept) %>%
#   mode_hdi() -> point.summaries.univariate.lnVR

brms.univariate.lnCVR %>%
  spread_draws(b_Intercept) %>%
  mode_hdi() -> point.summaries.univariate.lnCVR

brms.univariate.SMDH.ours %>%
  spread_draws(b_Intercept) %>%
  mode_hdi() -> point.summaries.univariate.SMDH.ours

# point summaries Egger's regressions
brms.Egger.lnRR.ours %>%
  spread_draws(b_Intercept) %>%
  mode_hdi() -> point.summaries.brms.Egger.lnRR.ours

brms.Egger.SMDH.ours %>%
  spread_draws(b_Intercept) %>%
  mode_hdi() -> point.summaries.brms.Egger.SMDH.ours


# meta-analytic modes
meta.modes <- round(c(unlist(point.summaries.univariate.lnRR.ours[["b_Intercept"]]),
                      #unlist(point.summaries.univariate.lnVR[["b_Intercept"]]),
                      unlist(point.summaries.univariate.lnCVR[["b_Intercept"]]),
                      unlist(point.summaries.univariate.SMDH.ours[["b_Intercept"]])),2)

# meta-analytic lower 2.5% CrIs (HDI)
meta.lower <- round(c(unlist(point.summaries.univariate.lnRR.ours[[".lower"]]),
                      #unlist(point.summaries.univariate.lnVR[[".lower"]]),
                      unlist(point.summaries.univariate.lnCVR[[".lower"]]),
                      unlist(point.summaries.univariate.SMDH.ours[[".lower"]])),2)

# meta-analytic upper 97.5% CrIs (HDI)
meta.upper <- round(c(unlist(point.summaries.univariate.lnRR.ours[[".upper"]]),
                      #unlist(point.summaries.univariate.lnVR[[".upper"]]),
                      unlist(point.summaries.univariate.lnCVR[[".upper"]]),
                      unlist(point.summaries.univariate.SMDH.ours[[".upper"]])),2)

# heterogeneities modes
total.I2.modes <- round(c(MCMCglmm::posterior.mode(I2_total.lnRR.ours),
                          #MCMCglmm::posterior.mode(I2_total.lnVR),
                          MCMCglmm::posterior.mode(I2_total.lnCVR),
                          MCMCglmm::posterior.mode(I2_total.SMDH.ours)),3)*100

esID.I2.modes <- round(c(MCMCglmm::posterior.mode(I2_esiD.lnRR.ours),
                         #MCMCglmm::posterior.mode(I2_esiD.lnVR),
                         MCMCglmm::posterior.mode(I2_esiD.lnCVR),
                         MCMCglmm::posterior.mode(I2_esiD.SMDH.ours)),3)*100

studyID.I2.modes <- round(c(MCMCglmm::posterior.mode(I2_studyID.lnRR.ours),
                            #MCMCglmm::posterior.mode(I2_studyID.lnVR),
                            MCMCglmm::posterior.mode(I2_studyID.lnCVR),
                            MCMCglmm::posterior.mode(I2_studyID.SMDH.ours)),3)*100

phylo.I2.modes <- round(c(MCMCglmm::posterior.mode(I2_phylo.lnRR.ours),
                          #MCMCglmm::posterior.mode(I2_phylo.lnVR),
                          MCMCglmm::posterior.mode(I2_phylo.lnCVR),
                          MCMCglmm::posterior.mode(I2_phylo.SMDH.ours)),3)*100

speciesID.I2.modes <- round(c(MCMCglmm::posterior.mode(I2_speciesID.lnRR.ours),
                              #MCMCglmm::posterior.mode(I2_speciesID.lnVR),
                              MCMCglmm::posterior.mode(I2_speciesID.lnCVR),
                              MCMCglmm::posterior.mode(I2_speciesID.SMDH.ours)),3)*100

# heterogeneities lower 2.5% CrIs (HDI)
total.I2.lower <- round(c(bayestestR::hdi(I2_total.lnRR.ours,ci = 0.95)$CI_low,
                          #bayestestR::hdi(I2_total.lnVR,ci = 0.95)$CI_low,
                          bayestestR::hdi(I2_total.lnCVR,ci = 0.95)$CI_low,
                          bayestestR::hdi(I2_total.SMDH.ours,ci = 0.95)$CI_low),3)*100

esID.I2.lower <- round(c(bayestestR::hdi(I2_esiD.lnRR.ours,ci = 0.95)$CI_low,
                         #bayestestR::hdi(I2_esiD.lnVR,ci = 0.95)$CI_low,
                         bayestestR::hdi(I2_esiD.lnCVR,ci = 0.95)$CI_low,
                         bayestestR::hdi(I2_esiD.SMDH.ours,ci = 0.95)$CI_low),3)*100

studyID.I2.lower <- round(c(bayestestR::hdi(I2_studyID.lnRR.ours,ci = 0.95)$CI_low,
                            #bayestestR::hdi(I2_studyID.lnVR,ci = 0.95)$CI_low,
                            bayestestR::hdi(I2_studyID.lnCVR,ci = 0.95)$CI_low,
                            bayestestR::hdi(I2_studyID.SMDH.ours,ci = 0.95)$CI_low),3)*100

phylo.I2.lower <- round(c(bayestestR::hdi(I2_phylo.lnRR.ours,ci = 0.95)$CI_low,
                          #bayestestR::hdi(I2_phylo.lnVR,ci = 0.95)$CI_low,
                          bayestestR::hdi(I2_phylo.lnCVR,ci = 0.95)$CI_low,
                          bayestestR::hdi(I2_phylo.SMDH.ours,ci = 0.95)$CI_low),3)*100

speciesID.I2.lower <- round(c(bayestestR::hdi(I2_speciesID.lnRR.ours,ci = 0.95)$CI_low,
                              #bayestestR::hdi(I2_speciesID.lnVR,ci = 0.95)$CI_low,
                              bayestestR::hdi(I2_speciesID.lnCVR,ci = 0.95)$CI_low,
                              bayestestR::hdi(I2_speciesID.SMDH.ours,ci = 0.95)$CI_low),3)*100

# heterogeneities upper 97.5% CrIs (HDI)
total.I2.upper <- round(c(bayestestR::hdi(I2_total.lnRR.ours,ci = 0.95)$CI_high,
                          #bayestestR::hdi(I2_total.lnVR,ci = 0.95)$CI_high,
                          bayestestR::hdi(I2_total.lnCVR,ci = 0.95)$CI_high,
                          bayestestR::hdi(I2_total.SMDH.ours,ci = 0.95)$CI_high),3)*100

esID.I2.upper <- round(c(bayestestR::hdi(I2_esiD.lnRR.ours,ci = 0.95)$CI_high,
                         #bayestestR::hdi(I2_esiD.lnVR,ci = 0.95)$CI_high,
                         bayestestR::hdi(I2_esiD.lnCVR,ci = 0.95)$CI_high,
                         bayestestR::hdi(I2_esiD.SMDH.ours,ci = 0.95)$CI_high),3)*100

studyID.I2.upper <- round(c(bayestestR::hdi(I2_studyID.lnRR.ours,ci = 0.95)$CI_high,
                            #bayestestR::hdi(I2_studyID.lnVR,ci = 0.95)$CI_high,
                            bayestestR::hdi(I2_studyID.lnCVR,ci = 0.95)$CI_high,
                            bayestestR::hdi(I2_studyID.SMDH.ours,ci = 0.95)$CI_high),3)*100

phylo.I2.upper <- round(c(bayestestR::hdi(I2_phylo.lnRR.ours,ci = 0.95)$CI_high,
                          #bayestestR::hdi(I2_phylo.lnVR,ci = 0.95)$CI_high,
                          bayestestR::hdi(I2_phylo.lnCVR,ci = 0.95)$CI_high,
                          bayestestR::hdi(I2_phylo.SMDH.ours,ci = 0.95)$CI_high),3)*100

speciesID.I2.upper <- round(c(bayestestR::hdi(I2_speciesID.lnRR.ours,ci = 0.95)$CI_high,
                              #bayestestR::hdi(I2_speciesID.lnVR,ci = 0.95)$CI_high,
                              bayestestR::hdi(I2_speciesID.lnCVR,ci = 0.95)$CI_high,
                              bayestestR::hdi(I2_speciesID.SMDH.ours,ci = 0.95)$CI_high),3)*100

# # bayesian R2 = how much variation in the response variable can be explained by our model using a Bayesian generalization of the R2 coefficient
# # bayesian R2 modes 
# bayesian.R2.modes <- round(c(MCMCglmm::posterior.mode(bayes_R2(brms.univariate.lnRR.ours,summary = FALSE)),
#                              MCMCglmm::posterior.mode(bayes_R2(brms.univariate.lnVR,summary = FALSE)),
#                              MCMCglmm::posterior.mode(bayes_R2(brms.univariate.lnCVR,summary = FALSE)),
#                              MCMCglmm::posterior.mode(bayes_R2(brms.univariate.SMDH.ours,summary = FALSE))),3)*100
# 
# # bayesian R2 lower 2.5% CrIs (HDI)
# bayesian.R2.lower <- round(c(bayestestR::hdi(bayes_R2(brms.univariate.lnRR.ours,summary = FALSE),ci = 0.95)$CI_low,
#                              bayestestR::hdi(bayes_R2(brms.univariate.lnVR,summary = FALSE),ci = 0.95)$CI_low,
#                              bayestestR::hdi(bayes_R2(brms.univariate.lnCVR,summary = FALSE),ci = 0.95)$CI_low,
#                              bayestestR::hdi(bayes_R2(brms.univariate.SMDH.ours,summary = FALSE),ci = 0.95)$CI_low),3)*100
# 
# # bayesian R2 upper 97.5% CrIs (HDI)
# bayesian.R2.upper <- round(c(bayestestR::hdi(bayes_R2(brms.univariate.lnRR.ours,summary = FALSE),ci = 0.95)$CI_high,
#                              bayestestR::hdi(bayes_R2(brms.univariate.lnVR,summary = FALSE),ci = 0.95)$CI_high,
#                              bayestestR::hdi(bayes_R2(brms.univariate.lnCVR,summary = FALSE),ci = 0.95)$CI_high,
#                              bayestestR::hdi(bayes_R2(brms.univariate.SMDH.ours,summary = FALSE),ci = 0.95)$CI_high),3)*100

# Q tests
Qtests <- round(c(lnRR.ours.Q,lnCVR.Q,SMDH.ours.Q),0)


# Egger's regression intercept modes
eggers.modes <- round(c(unlist(point.summaries.brms.Egger.lnRR.ours[["b_Intercept"]]),
                        #NA,
                        NA,
                        unlist(point.summaries.brms.Egger.SMDH.ours[["b_Intercept"]])),2)

# Egger's regression intercept lower 2.5% CrIs (HDI)
eggers.lower <- round(c(unlist(point.summaries.brms.Egger.lnRR.ours[[".lower"]]),
                        #NA,
                        NA,
                        unlist(point.summaries.brms.Egger.SMDH.ours[[".lower"]])),2)

# Egger's regression intercept upper 97.5% CrIs (HDI)
eggers.upper <- round(c(unlist(point.summaries.brms.Egger.lnRR.ours[[".upper"]]),
                        #NA,
                        NA,
                        unlist(point.summaries.brms.Egger.SMDH.ours[[".upper"]])),2)



# Building Table 1
table1 <- data.frame(effect.sizes,
                     ks,
                     paste0(sprintf("%.2f",meta.modes)," [", #sprintf allows to keep the number of digits regardless of 0's
                            sprintf("%.2f",meta.lower),",",
                            sprintf("%.2f",meta.upper),"]"),
                     paste0(sprintf("%.1f",esID.I2.modes)," [", #sprintf allows to keep the number of digits regardless of 0's
                            sprintf("%.1f",esID.I2.lower),",",
                            sprintf("%.1f",esID.I2.upper),"]"),
                     paste0(sprintf("%.1f",studyID.I2.modes)," [", #sprintf allows to keep the number of digits regardless of 0's
                            sprintf("%.1f",studyID.I2.lower),",",
                            sprintf("%.1f",studyID.I2.upper),"]"),
                     paste0(sprintf("%.1f",speciesID.I2.modes)," [", #sprintf allows to keep the number of digits regardless of 0's
                            sprintf("%.1f",speciesID.I2.lower),",",
                            sprintf("%.1f",speciesID.I2.upper),"]"),
                     paste0(sprintf("%.1f",phylo.I2.modes)," [", #sprintf allows to keep the number of digits regardless of 0's
                            sprintf("%.1f",phylo.I2.lower),",",
                            sprintf("%.1f",phylo.I2.upper),"]"),
                     paste0(sprintf("%.1f",total.I2.modes)," [", #sprintf allows to keep the number of digits regardless of 0's
                            sprintf("%.1f",total.I2.lower),",",
                            sprintf("%.1f",total.I2.upper),"]"),
                     Qtests,
                     paste0(sprintf("%.2f",eggers.modes)," [", #sprintf allows to keep the number of digits regardless of 0's
                            sprintf("%.2f",eggers.lower),",",
                            sprintf("%.2f",eggers.upper),"]"))#,
# paste0(sprintf("%.1f",bayesian.R2.modes)," [", #sprintf allows to keep the number of digits regardless of 0's
#        sprintf("%.1f",bayesian.R2.lower),",",
#        sprintf("%.1f",bayesian.R2.upper),"]"))


names(table1) <- table.column.names

# some tweaking for the table to look nice
table1[table1$Eggers=="NA [NA,NA]","Eggers"] <- NA

# reducing table for main manuscript
table1.red <- table1[1:2,]

table1.gt <- table1.red %>% 
  gt() %>% 
  cols_label(Effect.size=md("**Effect size**"),
             K=md("**k**"),
             Metaanalytic.mean=md("**Meta-analytic mean**"),
             Observational=md("***I*<sup>2</sup><sub>Obser.</sub> (%)**"),
             Study=md("***I*<sup>2</sup><sub>Study</sub> (%)**"),
             Species=md("***I*<sup>2</sup><sub>Species</sub> (%)**"),
             Phylogeny=md("***I*<sup>2</sup><sub>Phylo</sub> (%)**"),
             Total=md("***I*<sup>2</sup><sub> Total</sub> (%)**"),
             Qtest=md("***Q*<sub>test</sub>**"),
             Eggers=md("**Egger's test**")) %>%#,
  #Bayes.R2=md("**Bayesian R<sup>2</sup> (%)**")) %>%
  cols_align(align = "right") %>%
  tab_source_note(source_note = md("k = number of estimates; *I*<sup>2</sup> = heterogeneity; *Q*<sub>test</sub> = Cochrane's *Q* test; NA = not applicable; Obser. = Observational or residual variance; Phylo = Phylogeny. Egger's test = intercept of an Egger's regression following Nakagawa and Santos (2012). Estimates shown correspond to modes and 95% Highest Posterior Density Intervals")) %>%
  tab_options(table.width=775)


#table1.gt

# saving table
gtsave(table1.gt,filename="table1_meta-analysis.png", path="./tables/")


###########
# TABLE S1
###########

# reducing table for main manuscript
table1.SMDH <- table1[3,]

tableS1.gt <- table1.SMDH %>% 
  gt() %>% 
  cols_label(Effect.size=md("**Effect size**"),
             K=md("**k**"),
             Metaanalytic.mean=md("**Meta-analytic mean**"),
             Observational=md("***I*<sup>2</sup><sub>Obser.</sub> (%)**"),
             Study=md("***I*<sup>2</sup><sub>Study</sub> (%)**"),
             Species=md("***I*<sup>2</sup><sub>Species</sub> (%)**"),
             Phylogeny=md("***I*<sup>2</sup><sub>Phylo</sub> (%)**"),
             Total=md("***I*<sup>2</sup><sub> Total</sub> (%)**"),
             Qtest=md("***Q*<sub>test</sub>**"),
             Eggers=md("**Egger's test**")) %>%#,
  #Bayes.R2=md("**Bayesian R<sup>2</sup> (%)**")) %>%
  cols_align(align = "right") %>%
  tab_source_note(source_note = md("k = number of estimates; *I*<sup>2</sup> = heterogeneity; *Q*<sub>test</sub> = Cochrane's *Q* test; NA = not applicable; Obser. = Observational or residual variance; Phylo = Phylogeny. Egger's test = intercept of an Egger's regression (Nakagawa and Santos 2012). Estimates shown correspond to modes and 95% Highest Posterior Density Intervals")) %>%
  tab_options(table.width=775)


#table1.gt

# saving table
gtsave(tableS1.gt,filename="tableS1_meta-analysis.png", path="./tables/")


###########
# TABLE 2
###########

# using the package gt to create fancy tables: https://github.com/rstudio/gt
table.2.column.names <- c(#"Effect.size",
  "Effect.size.group",
  "Estimates","Posterior.mode")

# effect size names
# effect.sizes.2 <- c("lnRR",paste0("k = ",nrow(stress.data.metareg.lnRR.ours)),
#                     rep("",5))

#es2 <- c(paste0("lnRR (k = ",nrow(stress.data.metareg.lnRR.ours),")"))

effect.sizes.2.group <- c(rep(paste0("lnRR (k = ",nrow(stress.data.metareg.lnRR.ours),")"),7),
                          rep(paste0("lnRR (k = ",nrow(stress.data.metareg.lnRR.ours),"; time-lag bias test)"),3),
                          rep(paste0("lnCVR (k = ",nrow(stress.data.metareg.lnCVR),")"),7))

# names of estimates
estimates <- c("Intercept",
               "Development",
               "Metabolism and Physiology",
               "Morphology",
               "Reproduction",
               "Survival",
               "R2marginal",
               "Intercept",
               "Year of publication",
               "R2marginal",
               "Intercept",
               "Development",
               "Metabolism and Physiology",
               "Morphology",
               "Reproduction",
               "Survival",
               "R2marginal")

# get_variables(brms.univariate.lnRR.ours.trait)
# point summaries
brms.univariate.lnRR.ours.trait %>%
  spread_draws(c(b_Intercept,
                 b_trait.class.2development,
                 b_trait.class.2metabolism_and_physiology,
                 b_trait.class.2morphological,
                 b_trait.class.2reproduction,
                 b_trait.class.2survival)) %>%
  mode_hdi() -> point.summaries.univariate.lnRR.ours.trait

#get_variables(brms.univariate.lnRR.ours.year)
brms.univariate.lnRR.ours.year %>%
  spread_draws(c(b_Intercept,
                 b_year.z)) %>%
  mode_hdi() -> point.summaries.univariate.lnRR.ours.year

brms.univariate.lnCVR.trait %>%
  spread_draws(c(b_Intercept,
                 b_trait.class.2development,
                 b_trait.class.2metabolism_and_physiology,
                 b_trait.class.2morphological,
                 b_trait.class.2reproduction,
                 b_trait.class.2survival)) %>%
  mode_hdi() -> point.summaries.univariate.lnCVR.trait

brms.univariate.SMDH.ours.trait %>%
  spread_draws(c(b_Intercept,
                 b_trait.class.2development,
                 b_trait.class.2metabolism_and_physiology,
                 b_trait.class.2morphological,
                 b_trait.class.2reproduction,
                 b_trait.class.2survival)) %>%
  mode_hdi() -> point.summaries.univariate.SMDH.ours.trait


# meta-analytic modes
meta.modes <- round(c(unlist(point.summaries.univariate.lnRR.ours.trait[["b_Intercept"]]),
                      unlist(point.summaries.univariate.lnRR.ours.trait[["b_trait.class.2development"]]),
                      unlist(point.summaries.univariate.lnRR.ours.trait[["b_trait.class.2metabolism_and_physiology"]]),
                      unlist(point.summaries.univariate.lnRR.ours.trait[["b_trait.class.2morphological"]]),
                      unlist(point.summaries.univariate.lnRR.ours.trait[["b_trait.class.2reproduction"]]),
                      unlist(point.summaries.univariate.lnRR.ours.trait[["b_trait.class.2survival"]]),
                      MCMCglmm::posterior.mode(R2m.lnRR.trait),
                      unlist(point.summaries.univariate.lnRR.ours.year[["b_Intercept"]]),
                      unlist(point.summaries.univariate.lnRR.ours.year[["b_year.z"]]),
                      MCMCglmm::posterior.mode(R2m.lnRR.year),
                      unlist(point.summaries.univariate.lnCVR.trait[["b_Intercept"]]),
                      unlist(point.summaries.univariate.lnCVR.trait[["b_trait.class.2development"]]),
                      unlist(point.summaries.univariate.lnCVR.trait[["b_trait.class.2metabolism_and_physiology"]]),
                      unlist(point.summaries.univariate.lnCVR.trait[["b_trait.class.2morphological"]]),
                      unlist(point.summaries.univariate.lnCVR.trait[["b_trait.class.2reproduction"]]),
                      unlist(point.summaries.univariate.lnCVR.trait[["b_trait.class.2survival"]]),
                      MCMCglmm::posterior.mode(R2m.lnCVR.trait)),2)


# meta-analytic lower 2.5% CrIs (HDI)
meta.lower <- round(c(unlist(point.summaries.univariate.lnRR.ours.trait[["b_Intercept.lower"]]),
                      unlist(point.summaries.univariate.lnRR.ours.trait[["b_trait.class.2development.lower"]]),
                      unlist(point.summaries.univariate.lnRR.ours.trait[["b_trait.class.2metabolism_and_physiology.lower"]]),
                      unlist(point.summaries.univariate.lnRR.ours.trait[["b_trait.class.2morphological.lower"]]),
                      unlist(point.summaries.univariate.lnRR.ours.trait[["b_trait.class.2reproduction.lower"]]),
                      unlist(point.summaries.univariate.lnRR.ours.trait[["b_trait.class.2survival.lower"]]),
                      bayestestR::hdi(R2m.lnRR.trait,ci = 0.95)$CI_low,
                      unlist(point.summaries.univariate.lnRR.ours.year[["b_Intercept.lower"]]),
                      unlist(point.summaries.univariate.lnRR.ours.year[["b_year.z.lower"]]),
                      bayestestR::hdi(R2m.lnRR.year,ci = 0.95)$CI_low,
                      unlist(point.summaries.univariate.lnCVR.trait[["b_Intercept.lower"]]),
                      unlist(point.summaries.univariate.lnCVR.trait[["b_trait.class.2development.lower"]]),
                      unlist(point.summaries.univariate.lnCVR.trait[["b_trait.class.2metabolism_and_physiology.lower"]]),
                      unlist(point.summaries.univariate.lnCVR.trait[["b_trait.class.2morphological.lower"]]),
                      unlist(point.summaries.univariate.lnCVR.trait[["b_trait.class.2reproduction.lower"]]),
                      unlist(point.summaries.univariate.lnCVR.trait[["b_trait.class.2survival.lower"]]),
                      bayestestR::hdi(R2m.lnCVR.trait,ci = 0.95)$CI_low),2)

# meta-analytic upper 97.5% CrIs (HDI)
meta.upper <- round(c(unlist(point.summaries.univariate.lnRR.ours.trait[["b_Intercept.upper"]]),
                      unlist(point.summaries.univariate.lnRR.ours.trait[["b_trait.class.2development.upper"]]),
                      unlist(point.summaries.univariate.lnRR.ours.trait[["b_trait.class.2metabolism_and_physiology.upper"]]),
                      unlist(point.summaries.univariate.lnRR.ours.trait[["b_trait.class.2morphological.upper"]]),
                      unlist(point.summaries.univariate.lnRR.ours.trait[["b_trait.class.2reproduction.upper"]]),
                      unlist(point.summaries.univariate.lnRR.ours.trait[["b_trait.class.2survival.upper"]]),
                      bayestestR::hdi(R2m.lnRR.trait,ci = 0.95)$CI_high,
                      unlist(point.summaries.univariate.lnRR.ours.year[["b_Intercept.upper"]]),
                      unlist(point.summaries.univariate.lnRR.ours.year[["b_year.z.upper"]]),
                      bayestestR::hdi(R2m.lnRR.year,ci = 0.95)$CI_high,
                      unlist(point.summaries.univariate.lnCVR.trait[["b_Intercept.upper"]]),
                      unlist(point.summaries.univariate.lnCVR.trait[["b_trait.class.2development.upper"]]),
                      unlist(point.summaries.univariate.lnCVR.trait[["b_trait.class.2metabolism_and_physiology.upper"]]),
                      unlist(point.summaries.univariate.lnCVR.trait[["b_trait.class.2morphological.upper"]]),
                      unlist(point.summaries.univariate.lnCVR.trait[["b_trait.class.2reproduction.upper"]]),
                      unlist(point.summaries.univariate.lnCVR.trait[["b_trait.class.2survival.upper"]]),
                      bayestestR::hdi(R2m.lnCVR.trait,ci = 0.95)$CI_high),2)



# Building Table 2
table2 <- data.frame(#effect.sizes.2,
  effect.sizes.2.group,
  estimates,
  paste0(sprintf("%.2f",meta.modes)," [", #sprintf allows to keep the number of digits regardless of 0's
         sprintf("%.2f",meta.lower),",",
         sprintf("%.2f",meta.upper),"]"))#,

names(table2) <- table.2.column.names



table2.gt <- table2 %>% 
  gt(groupname_col = "Effect.size.group") %>%
  #tab_stubhead(label = "Effect.size") %>% 
  cols_label(Estimates=md("**Estimates**"),
             Posterior.mode=md("**Mode [95% HPDI]**")) %>%#,
  #Bayes.R2=md("**Bayesian R<sup>2</sup> (%)**")) %>%
  cols_align(align = "right",columns=c("Posterior.mode")) %>%
  cols_align(align = "left",columns=c(#"Effect.size",
    "Estimates")) %>%
  text_transform(locations = cells_data(columns = vars(Estimates)),
                 fn = function(x){ifelse(x=="R2marginal",
                                         md("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<em>R</em><sup>2</sup><sub> marginal</sub> (%) = "),
                                         x)}) %>%
  tab_source_note(source_note = md("k = number of estimates; *R*<sup>2</sup><sub>marginal</sub> = percentage of variance explained by the moderators (Nakagawa and Schielzeth 2013). The reference level (i.e. Intercept) corresponds to \"Behaviour\" (except for the time-lag bias test). Year of publication was z-transformed. Estimates shown correspond to posterior modes and 95% Highest Posterior Density Intervals (HPDI)"))%>%
  tab_options(table.width=397)

#table2.gt

# saving table
gtsave(table2.gt,filename="table2_meta-regressions.png", path="./tables/")



##############################################################
# Creating figures for the results
##############################################################

###########
# FIGURE 1
###########

tiff("plots/Figure1_meta-analysis.tiff",
     height=10, width=10,
     units='cm', compression="lzw", res=800)

#Plot the posterior values from the Bayesian model as density ridges
pd <- position_dodgev(height = 0.25)

post <- data.frame(c(rep(" lnRR",nrow(posterior.brms.univariate.lnRR.ours)),
                     #rep(" lnVR",nrow(posterior.brms.univariate.lnVR)),
                     rep("lnCVR",nrow(posterior.brms.univariate.lnCVR))),
                   #rep("SMDH",nrow(posterior.brms.univariate.SMDH.ours))),
                   c(posterior.brms.univariate.lnRR.ours$b_Intercept,
                     #posterior.brms.univariate.lnVR$b_Intercept,
                     posterior.brms.univariate.lnCVR$b_Intercept))#,
#posterior.brms.univariate.SMDH.ours$b_Intercept))

names(post) <- c("es","estimate")

posterior.plot <- post %>% 
  mutate(es = factor(es, levels = c("lnCVR",
                                    #" lnVR",
                                    " lnRR"))) %>%
  ggplot() + 
  stat_density_ridges(aes(x=estimate, y = es), 
                      alpha = 0.5, scale = 0.6, 
                      position = position_nudge(y = 0.12), 
                      height = 10, 
                      show.legend = F, 
                      quantile_lines = T, 
                      quantiles = 2,
                      fill=rgb(25/255,100/255,205/255, 0.75)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "black") + 
  ylab("")+
  xlab("Effect size")+
  # scale_fill_manual(values = c("Male" = "#ff7f00", "Female" = "#984ea3", "Both" = "#4daf4a"))+
  # scale_color_manual(values = c("Male" = "#ff7f00", "Female" = "#984ea3", "Both" = "#4daf4a"))+
  scale_x_continuous(limits = c(-0.35, 0.20), breaks = c(-0.3,-0.2,-0.1,0,0.1,0.2)) +
  theme_bw() +
  theme(panel.spacing = unit(0.1, "lines"),
        text = element_text(size=16),
        panel.border= element_blank(),
        axis.line=element_line(), 
        axis.ticks.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(), 
        legend.text = element_text(size=16), 
        legend.title=element_text(size=16, face = "bold"),
        axis.title.x = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, hjust = 0.35, margin = margin(r=-5)),
        axis.text.y = element_text(angle = 0, face = "bold",color="black",hjust=-0.2),
        axis.text.x = element_text(color="black"),
        plot.title = element_text(size = 16))

# rgb(1,165/255,0,0.5)


#Add the modes as circles with error bars showing the 95% HDI
p.sum.uni <- data.frame(c(" lnRR",
                          #"lnVR",
                          "lnCVR"),
                        ks[1:2],
                        meta.modes[1:2],
                        meta.lower[1:2],
                        meta.upper[1:2])
names(p.sum.uni) <- c("es","k","mode","lower","upper")


both.plots <- posterior.plot + 
  geom_errorbarh(data = p.sum.uni %>% mutate(es = factor(es, levels = c("lnCVR",
                                                                        #" lnVR",
                                                                        " lnRR"))), 
                 aes(xmin = p.sum.uni$lower,
                     xmax = p.sum.uni$upper, y = es), 
                 height = 0, show.legend = F, position = pd,
                 color=rgb(25/255,100/255,205/255, 0.99))+
  geom_point(data = p.sum.uni %>% mutate(es = factor(es, levels = c("lnCVR",
                                                                    #" lnVR",
                                                                    " lnRR"))),
             aes(x = mode, y = es, size=2), #k
             shape=21, fill = rgb(25/255,100/255,205/255, 0.99), position = pd, show.legend = F) 

both.plots

dev.off()


###########
# FIGURE 2
###########

# First thing is to get the data in the right format. Manually  
# implemented rather than using, for example, marginal_effects()

########################
# lnRR: point summaries
########################

# this way allows to extract mode and HPDI for each level 
newdat.lnRR.ours.trait <- expand.grid(trait.class.2 = levels(as.factor(stress.data.metareg.lnRR.ours$trait.class.2)))

# fixed effects design matrix
xmat.lnRR.ours.trait <- model.matrix(~trait.class.2,data=newdat.lnRR.ours.trait)

fitmatboth.lnRR.ours.trait <- matrix(NA, 
                                     ncol = nrow(posterior.brms.univariate.lnRR.ours.trait), 
                                     nrow = nrow(newdat.lnRR.ours.trait))

posterior.db.lnRR.ours.trait <- matrix(c(posterior.brms.univariate.lnRR.ours.trait$b_Intercept,
                                         posterior.brms.univariate.lnRR.ours.trait$b_trait.class.2development,
                                         posterior.brms.univariate.lnRR.ours.trait$b_trait.class.2metabolism_and_physiology,
                                         posterior.brms.univariate.lnRR.ours.trait$b_trait.class.2morphological,
                                         posterior.brms.univariate.lnRR.ours.trait$b_trait.class.2reproduction,
                                         posterior.brms.univariate.lnRR.ours.trait$b_trait.class.2survival),
                                       ncol = 6)

for(i in 1:nrow(posterior.brms.univariate.lnRR.ours.trait)) {
  fitmatboth.lnRR.ours.trait[,i] <- xmat.lnRR.ours.trait%*%posterior.db.lnRR.ours.trait[i,]
}

# getting point summaries
newdat.lnRR.ours.trait$mode <- round(apply(fitmatboth.lnRR.ours.trait, 1, MCMCglmm::posterior.mode),2)

HPDI.lnRR.ours.trait <- apply(fitmatboth.lnRR.ours.trait, 1, bayestestR::hdi,ci = 0.95)

newdat.lnRR.ours.trait$lower <- round(sapply(HPDI.lnRR.ours.trait, function(x)x[["CI_low"]]),2)
newdat.lnRR.ours.trait$upper <- round(sapply(HPDI.lnRR.ours.trait, function(x)x[["CI_high"]]),2)

#newdat.lnRR.ours.trait


########################
# lnCVR: point summaries
########################

# this way allows to extract mode and HPDI for each level 
newdat.lnCVR.trait <- expand.grid(trait.class.2 = levels(as.factor(stress.data.metareg.lnCVR$trait.class.2)))

# fixed effects design matrix
xmat.lnCVR.trait <- model.matrix(~trait.class.2,data=newdat.lnCVR.trait)

fitmatboth.lnCVR.trait <- matrix(NA, 
                                 ncol = nrow(posterior.brms.univariate.lnCVR.trait), 
                                 nrow = nrow(newdat.lnCVR.trait))

posterior.db.lnCVR.trait <- matrix(c(posterior.brms.univariate.lnCVR.trait$b_Intercept,
                                     posterior.brms.univariate.lnCVR.trait$b_trait.class.2development,
                                     posterior.brms.univariate.lnCVR.trait$b_trait.class.2metabolism_and_physiology,
                                     posterior.brms.univariate.lnCVR.trait$b_trait.class.2morphological,
                                     posterior.brms.univariate.lnCVR.trait$b_trait.class.2reproduction,
                                     posterior.brms.univariate.lnCVR.trait$b_trait.class.2survival),
                                   ncol = 6)

for(i in 1:nrow(posterior.brms.univariate.lnCVR.trait)) {
  fitmatboth.lnCVR.trait[,i] <- xmat.lnCVR.trait%*%posterior.db.lnCVR.trait[i,]
}

# getting point summaries
newdat.lnCVR.trait$mode <- round(apply(fitmatboth.lnCVR.trait, 1, MCMCglmm::posterior.mode),2)

HPDI.lnCVR.trait <- apply(fitmatboth.lnCVR.trait, 1, bayestestR::hdi,ci = 0.95)

newdat.lnCVR.trait$lower <- round(sapply(HPDI.lnCVR.trait, function(x)x[["CI_low"]]),2)
newdat.lnCVR.trait$upper <- round(sapply(HPDI.lnCVR.trait, function(x)x[["CI_high"]]),2)

#newdat.lnCVR.trait

# full dataset with all point summaries
newdat.trait <- rbind(newdat.lnRR.ours.trait,newdat.lnCVR.trait)
newdat.trait$es <- c(rep(" lnRR",6),rep("lnCVR",6))

# adding sample sizes to it
lnRR.ks <- stress.data.metareg.lnRR.ours %>% 
  group_by(trait.class.2) %>% 
  summarise(count = n()) %>%
  select(count)

lnCVR.ks <- stress.data.metareg.lnCVR %>% 
  group_by(trait.class.2) %>% 
  summarise(count = n()) %>%
  select(count)

newdat.trait$k <- c(lnRR.ks$count,lnCVR.ks$count)
newdat.trait$trait.colour <- rep(c("#999999","#E69F00","#56B4E9","#009E73","#D55E00","#CC79A7"),2)

# renaming variables and levels for visualization
newdat.trait <- newdat.trait %>% rename(Trait = trait.class.2)
newdat.trait <- newdat.trait %>%
  mutate(Trait = recode(Trait, 
                        behavioural = "Behaviour",
                        development = "Development",
                        metabolism_and_physiology = "Metabolism and Physiology",
                        morphological = "Morphology",
                        reproduction = "Reproduction",
                        survival = "Survival"))

# table(stress.data.metareg.lnRR.ours$trait.class.2)
# table(stress.data.metareg.lnCVR$trait.class.2)


###########################
# lnRR and lnCVR: posterior
###########################

post.metaregression <- data.frame(c(rep(" lnRR",ncol(fitmatboth.lnRR.ours.trait)*nrow(fitmatboth.lnRR.ours.trait)),
                                    rep("lnCVR",ncol(fitmatboth.lnCVR.trait)*nrow(fitmatboth.lnCVR.trait))),
                                  c(rep("Behaviour",length(fitmatboth.lnRR.ours.trait[1,])),
                                    rep("Development",length(fitmatboth.lnRR.ours.trait[1,])),
                                    rep("Metabolism and Physiology",length(fitmatboth.lnRR.ours.trait[1,])),
                                    rep("Morphology",length(fitmatboth.lnRR.ours.trait[1,])),
                                    rep("Reproduction",length(fitmatboth.lnRR.ours.trait[1,])),
                                    rep("Survival",length(fitmatboth.lnRR.ours.trait[1,])),
                                    rep("Behaviour",length(fitmatboth.lnCVR.trait[1,])),
                                    rep("Development",length(fitmatboth.lnCVR.trait[1,])),
                                    rep("Metabolism and Physiology",length(fitmatboth.lnCVR.trait[1,])),
                                    rep("Morphology",length(fitmatboth.lnCVR.trait[1,])),
                                    rep("Reproduction",length(fitmatboth.lnCVR.trait[1,])),
                                    rep("Survival",length(fitmatboth.lnCVR.trait[1,]))),
                                  c(fitmatboth.lnRR.ours.trait[1,],
                                    fitmatboth.lnRR.ours.trait[2,],
                                    fitmatboth.lnRR.ours.trait[3,],
                                    fitmatboth.lnRR.ours.trait[4,],
                                    fitmatboth.lnRR.ours.trait[5,],
                                    fitmatboth.lnRR.ours.trait[6,],
                                    fitmatboth.lnCVR.trait[1,],
                                    fitmatboth.lnCVR.trait[2,],
                                    fitmatboth.lnCVR.trait[3,],
                                    fitmatboth.lnCVR.trait[4,],
                                    fitmatboth.lnCVR.trait[5,],
                                    fitmatboth.lnCVR.trait[6,]))

names(post.metaregression) <- c("es","Trait","estimate")

#######################
# Actual plotting time

# this should be a colorblind-friendly pallete (from: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/)
#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


#Plot the posterior values from the Bayesian model as density ridges
pd <- position_dodgev(height = 0.4)

tiff("plots/Figure2_meta-regressions.tiff",
     height=10, width=15,
     units='cm', compression="lzw", res=800)

posterior.metaregression.plot <- post.metaregression %>% 
  mutate(es = factor(es, levels = c("lnCVR",
                                    " lnRR"))) %>%
  ggplot() + 
  stat_density_ridges(aes(x = estimate, y = es, fill = Trait), 
                      alpha = 0.65, scale = 0.6, 
                      position = position_nudge(y = 0.20), 
                      height = 10, 
                      show.legend = F, 
                      quantile_lines = T, 
                      quantiles = 2) +
  geom_vline(xintercept = 0, linetype = 2, colour = "black") + 
  ylab("")+
  xlab("Effect size")+
  scale_fill_manual(values = c("Behaviour" = "#999999", "Development" = "#E69F00", "Metabolism and Physiology" = "#56B4E9",
                               "Morphology" = "#009E73", "Reproduction" = "#D55E00", "Survival" = "#CC79A7"))+
  scale_x_continuous(limits = c(-0.45, 0.45), breaks = c(-0.4,-0.2,0,0.2,0.4)) +
  theme_bw() +
  theme(panel.spacing = unit(0.1, "lines"),
        text = element_text(size=16),
        panel.border= element_blank(),
        axis.line=element_line(), 
        axis.ticks.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(), 
        legend.text = element_text(size=7.5), 
        legend.title=element_text(size=7.5, face = "bold"),
        legend.position = c(0.86, 0.87),
        legend.key.size = unit(0.25, "cm"),
        axis.title.x = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, hjust = 0.35, margin = margin(r=-5)),
        axis.text.y = element_text(angle = 0, face = "bold",color="black",hjust=-0.2),
        axis.text.x = element_text(color="black"),
        plot.title = element_text(size = 16))

# rgb(1,165/255,0,0.5)
#posterior.metaregression.plot

both.plots <- posterior.metaregression.plot + 
  geom_errorbarh(data = newdat.trait %>% mutate(es = factor(es, levels = c("lnCVR"," lnRR")),
                                                Trait =  factor(Trait, levels = c("Survival",
                                                                                  "Reproduction",
                                                                                  "Morphology",
                                                                                  "Metabolism and Physiology",
                                                                                  "Development",
                                                                                  "Behaviour"))), 
                 aes(xmin = newdat.trait$lower,
                     xmax = newdat.trait$upper, 
                     y = es, color = Trait), 
                 height = 0, show.legend = F, position = pd) +
  geom_point(data = newdat.trait %>% mutate(es = factor(es, levels = c("lnCVR"," lnRR")),
                                            Trait =  factor(Trait, levels = c("Survival",
                                                                              "Reproduction",
                                                                              "Morphology",
                                                                              "Metabolism and Physiology",
                                                                              "Development",
                                                                              "Behaviour"))),
             aes(x = mode, y = es, size = k, fill = Trait),
             shape=21, color = "grey20", position = pd) +
  scale_color_manual(values = c("Behaviour" = "#999999", "Development" = "#E69F00", "Metabolism and Physiology" = "#56B4E9",
                                "Morphology" = "#009E73", "Reproduction" = "#D55E00", "Survival" = "#CC79A7"))+
  guides(fill = guide_legend(reverse=F, override.aes = list(size = 2)))+
  scale_size(guide = 'none')+
  scale_y_discrete(expand=c(0.075, 0))

both.plots

dev.off()
