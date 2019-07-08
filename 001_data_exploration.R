
##############################################################
# Authors: 
# Alfredo Sanchez-Tojar (@ASanchez_Tojar)
# Profile: https://goo.gl/PmpPEB
# Department of Evolutionary Biology, Bielefeld University (GER) 
# Email: alfredo.tojar@gmail.com

# Nicholas P. Moran
# Profile: https://www.researchgate.net/profile/Nicholas_Moran
# Department of Evolutionary Biology, Bielefeld University (GER) 
# Email: nicholaspatrickmoran@gmail.com

# Rose O'Dea (@rose_odea)
# Profile: https://www.roseodea.com/
# I-DEEL, UNSW, Sydney (AUSTRALIA) 
# Email: rose.eleanor.o.dea@gmail.com

# Script first created on the 25th of March 2019

##############################################################
# Description of script and Instructions
##############################################################

# This script is to explore the meta-analytic data from:

# Eyck et al. 2019: Effects of developmental stress on animal
# phenotype and performance: a quantitative review

# H.Eyck shared with us a dataset that combines the data presented
# in their published table 2:

# https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fbrv.12496&file=brv12496-sup-0002-TableS2.xlsx

# with the corresponding raw data, which is included in their
# published table 4:

# https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fbrv.12496&file=brv12496-sup-0004-TableS4.xlsx

# We will use for the analyses since it contains the raw data of 
# the selection of studies that were used in the meta-analysis 
# (note that there are 865 total effect sizes in this dataset, in 
# contrast to the published data set, which contained 866 effect sizes).

##############################################################
# Packages needed
##############################################################

pacman::p_load(metafor, doBy, metaDigitise, dplyr, tibble,
               ggplot2, plotly, brms, tidybayes, tidyr,
               modelr)

# Clear memory
rm(list=ls())


##############################################################
# Functions needed
##############################################################

# functions needed to estimate SMD (or Cohen's d) from raw
# data

sdpooled <- function(n1,n2,sd1,sd2){
  sp <- sqrt((((n1-1)*sd1^2)+(n2-1)*sd2^2)/(n1+n2-2))
  #sp <- sqrt((((n1)*sd1^2)+(n2)*sd2^2)/(n1+n2+2))
}

means.to.d <- function(x1,x2,sp){
  d <- (x1-x2)/sp
}

v.from.raw<-function(n1,n2,d){
  v<-(1/n1) + (1/n2) + (d^2)/(2*(n1+n2))
}


################################################################################
# Comparing new dataset with the published dataset (table 2)
################################################################################

# importing published dataset
eyck.efs <- read.xlsx("data_re-extraction/brv12496-sup-0002-tables2.xlsx",
                      colNames=T,sheet = 1)

# importing new dataset shared by H.Eyck
eyck.dat <- read.xlsx("data/EyckDev.stress_Data_FULL_TABLE.xlsx",
                      colNames=T,sheet = 1)


# do both datasets include the same number of studies?
length(unique(eyck.efs$sdyID))==length(unique(eyck.dat$studyID))


# do both datasets include the same number of effect sizes?
nrow(eyck.efs)==nrow(eyck.dat)


# how many effect sizes are missing from the new dataset?
nrow(eyck.efs)-nrow(eyck.dat)


# let's find out which effect size is missing
count1<-plyr::count(eyck.efs$sdyID)
count2<-plyr::count(eyck.dat$studyID)

setdiff(count1$freq, count2$freq)
setdiff(count2$freq, count1$freq)

count3<-count1$freq - count2$freq

studyID.missing <- count1[match(c(1),count3),"x"]
#studyID 93 has 1 missing effect size in the the eyck.dat database


# is the missing effect size reported in the published table 4?
# if it is, we could added to the dataset we are using. If not,
# we will have to exclude it from our final analyses

# importing table 4
eyck.table4 <- read.xlsx("data_re-extraction/brv12496-sup-0004-tables4.xlsx",
                      colNames=T,sheet = 1)

eyck.table4.93 <- subset(eyck.table4, studyID==studyID.missing)

eyck.efs93<-subset(eyck.efs, sdyID==studyID.missing)

nrow(eyck.table4.93)==nrow(eyck.efs93)

# it does not seems so as the number of effect sizes differ

# let's find that effect size to make sure we are missing it
eyck.efs93<-subset(eyck.efs, sdyID==studyID.missing)
eyck.dat93<-subset(eyck.dat, studyID==studyID.missing)
names(eyck.efs93)
names(eyck.dat93)

# looking for it
eyck.efs93.cut<-round(eyck.efs93[,c("d","sv")], 2)
eyck.dat93.cut<-round(eyck.dat93[,c("Cohen's.D","sv")],2)
names(eyck.dat93.cut)<-c("d","sv")

setdiff(eyck.dat93.cut$d, eyck.efs93.cut$d)
setdiff(eyck.efs93.cut$d, eyck.dat93.cut$d)

# identifiying all information available for the missind effect size
eyck.efs93[round(eyck.efs93$d,2)==setdiff(eyck.efs93.cut$d, eyck.dat93.cut$d),]


# is the raw data available in table 4?
# in table 2 the missing effect size is the only morphological trait from study 93
eyck.table4.93[eyck.table4.93$trait.class=="morphological"]


# this doesn't find anything, which we visually double-checked.
# Thus, we will ignore this missing effect size for our reanalysis 
# and use the 865 effect sizes available in the dataset that H.Eyck
# shared with us.


################################################################################
# Understanding how effect sizes were calculated
################################################################################

#stress.red <- read.table("data/EyckDev.stress_Data_FULL_TABLE_red.csv",header=T,sep=',')
stress.red <- read.xlsx("data/EyckDev.stress_Data_FULL_TABLE.xlsx",
                        colNames=T,sheet = 1)

# Estimating Cohen's d using equations from W.Viechtbauer course
stress.red$Cohens.2 <- round(means.to.d(stress.red$mean.treat,stress.red$mean.control,sdpooled(stress.red$N.treat,stress.red$N.control,stress.red$SD.1,stress.red$SD)),2)
stress.red$sv.2 <- round(v.from.raw(stress.red$N.treat,stress.red$N.control,means.to.d(stress.red$mean.treat,stress.red$mean.control,sdpooled(stress.red$N.treat,stress.red$N.control,stress.red$SD.1,stress.red$SD))),2)


# Estimating Hedge's g using metafor
x<-escalc(measure="SMD",
                     m1i=mean.treat,n1i=N.treat,sd1i=SD.1,
                     m2i=mean.control,n2i=N.control,sd2i=SD,
                     data=stress.red, append=F,digits=2)

stress.red$yi <- round(x[["yi"]][1:nrow(stress.red)],2)
stress.red$vi <- round(x[["vi"]][1:nrow(stress.red)],2)


# Estimating standardized mean difference with heteroscedastic population 
# variances in the two groups (Bonett, 2008, 2009). In case that's what
# they calculated.
y<-escalc(measure="SMDH",
          m1i=mean.treat,n1i=N.treat,sd1i=SD.1,
          m2i=mean.control,n2i=N.control,sd2i=SD,
          data=stress.red, append=F,digits=2)

stress.red$yi.SMDH <- round(y[["yi"]][1:nrow(stress.red)],2)
stress.red$vi.SMDH <- round(y[["vi"]][1:nrow(stress.red)],2)


# eliminating effect sizes based on test statistics for the time being
stress.red.2 <- stress.red[!(is.na(stress.red$Cohens.2)),]

# proportion of effect sizes derived from test statistics
nrow(stress.red.2)/nrow(stress.red)

# 83% of effect sizes are based on means, sds and ns. These are the 
# effect sizes that we can potentially use to estimate lnRR, lnVR
# and lnCVR.
 

# estimating how many Cohen's d values differ between the two approaches
table(stress.red.2$Cohens.D == stress.red.2$Cohens.2)
table(stress.red.2$Cohens.D == stress.red.2$Cohens.2)[1]/nrow(stress.red.2)
table(abs(stress.red.2$Cohens.D) == abs(stress.red.2$Cohens.2))
table(abs(stress.red.2$Cohens.D) == abs(stress.red.2$Cohens.2))[1]/nrow(stress.red.2)

# There are 459 (64%) effect sizes for which the two approaches differ


# Value are larger for Hedges' g and SMDH, so this is not the answer for 
# why the approaches differ
table(abs(stress.red$Cohens.D) == abs(stress.red$yi))
table(abs(stress.red$Cohens.D) == abs(stress.red$yi.SMDH))


# estimating how many sampling variances differ between the two approaches
table(stress.red.2$sv == stress.red.2$sv.2)
table(stress.red.2$sv == stress.red.2$sv.2)[1]/nrow(stress.red.2)

# Regarding the sampling variance, there are 697 (97%) sampling variances
# that differe between the two approaches.


# Further exploring whether the difference tends to be consistently 
# towards one direction

# subsetting the database for those values that are different between the
# two approaches
stress.red.dif <- stress.red.2[stress.red.2$Cohens.D != stress.red.2$Cohens.2,]


# are their values generally larger?
table(abs(stress.red.dif$Cohens.D) > abs(stress.red.dif$Cohens.2))
table(abs(stress.red.dif$Cohens.D) > abs(stress.red.dif$Cohens.2))[2]/nrow(stress.red.dif)


summary(abs(stress.red.dif[stress.red.dif$Cohens.2<3000,"Cohens.D"]))
summary(abs(stress.red.dif[stress.red.dif$Cohens.2<3000,"Cohens.2"])) #typo in one of the means, it seems

summary(stress.red.dif[stress.red.dif$Cohens.2<3000,"Cohens.D"])
summary(stress.red.dif[stress.red.dif$Cohens.2<3000,"Cohens.2"]) #typo in one of the means, it seems


# absolute difference between effect sizes
summary(abs(stress.red.dif[stress.red.dif$Cohens.2<3000,"Cohens.D"]-stress.red.dif[stress.red.dif$Cohens.2<3000,"Cohens.2"]))

# The effect sizes used in the paper are larger that the ones I have
# calculated in 97% of the cases. Additionally, their effect sizes tend
# to be 0.04 larger at the median level.


# What about the sampling variances?
stress.red.dif.sv <- stress.red.2[stress.red.2$sv != stress.red.2$sv.2,]

table(stress.red.dif.sv$sv < stress.red.dif.sv$sv.2)
summary(stress.red.dif[stress.red.dif$sv.2<100000,"sv"])
summary(stress.red.dif[stress.red.dif$sv.2<100000,"sv.2"]) #typo in one of the means, it seems

# Around 50% of the variances are larger, and vice versa.
# But sampling variances are generally smaller in the database used
# for the paper