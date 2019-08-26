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

# This script is to (re-)calculate effect sizes for:

# Eyck et al. 2019: Effects of developmental stress on animal
# phenotype and performance: a quantitative review

# using a cleaner version of the dataset (more in scripts 001-004).

# Note that some of the effect size statistic calculated in this
# script were never analyzed (e.g. SMD, .HE, etc...)


##############################################################
# Packages needed
##############################################################

pacman::p_load(openxlsx,metafor, doBy, metaDigitise, dplyr, tibble,
               ggplot2, plotly, brms, tidybayes, tidyr,modelr)

# Clear memory
rm(list=ls())


##############################################################
# Functions needed
##############################################################

# function that calculates biased (used in the original meta-analysis)
# and unbiased cohen's d, and their sampling variances
cohensD <- function(mean_control,mean_treatment,sd_control,sd_treatment,n_control,n_treatment,append=T,data){
  sd_pooled <- sqrt(((data[,n_treatment]-1)*data[,sd_treatment]^2+(data[,n_control]-1)*data[,sd_control]^2)/(data[,n_control]+data[,n_treatment]-2))
  correction <- sqrt((data[,n_treatment]+data[,n_control]-2)/(data[,n_treatment]+data[,n_control]))
  cohens.biased <- (data[,mean_treatment]-data[,mean_control])/(sd_pooled*correction)
  cohens.biased.sv <- (((data[,n_control]+data[,n_treatment])/(data[,n_control]*data[,n_treatment]))+((cohens.biased^2)/(2*(data[,n_control]+data[,n_treatment]))))
  cohens.unbiased <- (data[,mean_treatment]-data[,mean_control])/sd_pooled
  cohens.unbiased.sv <- (((data[,n_control]+data[,n_treatment])/(data[,n_control]*data[,n_treatment]))+((cohens.unbiased^2)/(2*(data[,n_control]+data[,n_treatment]))))
  if(append){
    return(cbind(data,cohens.biased,cohens.biased.sv,cohens.unbiased,cohens.unbiased.sv))
  } else{
    return(as.data.frame(cbind(data,cohens.biased,cohens.biased.sv,cohens.unbiased,cohens.unbiased.sv)))
  }
}

se.to.sd <- function(se,n){
  sd <- se*sqrt(n)
  return(sd)
}

# # Code to test that the function works as expected:
# # indeed the original study used biased Cohen's D but 
# # it is unclear how they estimated sampling variance
# # importing dataset shared by H.Eyck
# eyck.dat <- read.xlsx("data/EyckDev.stress_Data_FULL_TABLE.xlsx",
#                       colNames=T,sheet = 1)
# 
# eyck.dat.test <- eyck.dat[1:5,]
# cols <- which(names(eyck.dat.test)=="SD")
# names(eyck.dat.test)[cols] <- paste0("SD", seq.int(1,2))
# 
# cohensD(mean_control="mean.control",mean_treatment="mean.treat",
#         sd_control="SD1",sd_treatment="SD2",
#         n_control="N.control",n_treatment="N.treat",
#         append=T,data=eyck.dat.test)

##############################################################
# Importing datasets
##############################################################

# database with the corrected data
stress.data <- read.xlsx("data_re-extraction/clean_data/EyckDev_stress_clean_raw_data_shared_control.xlsx",
                         colNames=T,sheet = 1)

stress.data$sign.inversion.HE<-as.numeric(as.character(stress.data$sign.inversion.HE))
stress.data$sign.inversion.ours<-as.numeric(as.character(stress.data$sign.inversion.ours))


# imputing values for those two SDs that were reported as 0
# this was found via the outlier exploration we did based
# on funnel plotting (more info below, see also section
# "Double-checking outliers")
stress.data[stress.data$esID==361,"SD.control"] <- 0.14*((sum(se.to.sd(c(0.01,0.02,0.01,0.01),
                                                                       c(41,30,32,32))))/
                                                           (sum(c(0.15,0.14,0.15,0.16,0.17))))

stress.data[stress.data$esID==546,"SD.control"] <- 0.08*((sum(se.to.sd(c(0.01),
                                                                       c(49))))/
                                                           (sum(c(0.08))))


##############################################################
# Calculating Cohen's D biased (used in original publication)
##############################################################

stress.data.for.cohens <- stress.data[!(is.na(stress.data$mean.control)),
                                      c("esID","mean.control","SD.control","N.control","mean.treat","SD.treat","N.treat","sign.inversion.HE","sign.inversion.ours")]

options(scipen = 999)#avoiding scientific notation, it was annoying! :)
stress.data.for.cohens <- cohensD(mean_control="mean.control",mean_treatment="mean.treat",
                                  sd_control="SD.control",sd_treatment="SD.treat",
                                  n_control="N.control",n_treatment="N.treat",
                                  append=T,data=stress.data.for.cohens)


stress.data.for.cohens$cohens.biased.HE <- stress.data.for.cohens$cohens.biased*stress.data.for.cohens$sign.inversion.HE
stress.data.for.cohens$cohens.unbiased.HE <- stress.data.for.cohens$cohens.unbiased*stress.data.for.cohens$sign.inversion.HE
stress.data.for.cohens$cohens.biased.ours <- stress.data.for.cohens$cohens.biased*stress.data.for.cohens$sign.inversion.ours
stress.data.for.cohens$cohens.unbiased.ours <- stress.data.for.cohens$cohens.unbiased*stress.data.for.cohens$sign.inversion.ours

stress.data.cohens.red <- stress.data.for.cohens[,c("esID","cohens.biased.HE","cohens.biased.ours","cohens.biased.sv",
                                                    "cohens.unbiased.HE","cohens.unbiased.ours","cohens.unbiased.sv")]



##############################################################
# Calculating Cohen's D biased: shared control
##############################################################

stress.data.for.cohens.sc <- stress.data[!(is.na(stress.data$mean.control)),
                                         c("esID","mean.control","SD.control","N.control.sc","mean.treat","SD.treat","N.treat","sign.inversion.HE","sign.inversion.ours")]

# renaming the variable so that it is consistent throughtout
colnames(stress.data.for.cohens.sc)[colnames(stress.data.for.cohens.sc)=="N.control.sc"] <- "N.control"


options(scipen = 999)#avoiding scientific notation, it was annoying! :)
stress.data.for.cohens.sc <- cohensD(mean_control="mean.control",mean_treatment="mean.treat",
                                     sd_control="SD.control",sd_treatment="SD.treat",
                                     n_control="N.control",n_treatment="N.treat",
                                     append=T,data=stress.data.for.cohens.sc)


stress.data.for.cohens.sc$cohens.biased.sc.HE <- stress.data.for.cohens.sc$cohens.biased*stress.data.for.cohens.sc$sign.inversion.HE
stress.data.for.cohens.sc$cohens.unbiased.sc.HE <- stress.data.for.cohens.sc$cohens.unbiased*stress.data.for.cohens.sc$sign.inversion.HE
stress.data.for.cohens.sc$cohens.biased.sc.ours <- stress.data.for.cohens.sc$cohens.biased*stress.data.for.cohens.sc$sign.inversion.ours
stress.data.for.cohens.sc$cohens.unbiased.sc.ours <- stress.data.for.cohens.sc$cohens.unbiased*stress.data.for.cohens.sc$sign.inversion.ours

colnames(stress.data.for.cohens.sc)[colnames(stress.data.for.cohens.sc)=="cohens.biased.sv"] <- "cohens.biased.sc.sv"
colnames(stress.data.for.cohens.sc)[colnames(stress.data.for.cohens.sc)=="cohens.unbiased.sv"] <- "cohens.unbiased.sc.sv"


stress.data.cohens.sc.red <- stress.data.for.cohens.sc[,c("esID","cohens.biased.sc.HE","cohens.biased.sc.ours","cohens.biased.sc.sv",
                                                          "cohens.unbiased.sc.HE","cohens.unbiased.sc.ours","cohens.unbiased.sc.sv")]



##############################################################
# Calculating Hedge's g
##############################################################

es.list.SMD <- stress.data[!(is.na(stress.data$mean.control)),"esID"]

# Homocedasticity assumed
stress.data.SMD<-escalc(measure="SMD", 
                        n1i=N.treat, n2i=N.control, 
                        m1i=mean.treat, m2i=mean.control, 
                        sd1i=SD.treat, sd2i=SD.control,
                        data=stress.data,
                        subset=stress.data$esID %in% es.list.SMD,
                        var.names=c("SMD","SMD.sv"), add.measure=FALSE,
                        append=TRUE)

stress.data.SMD <- as.data.frame(stress.data.SMD)
stress.data.SMD$SMD.HE <- stress.data.SMD$SMD*stress.data.SMD$sign.inversion.HE
stress.data.SMD$SMD.ours <- stress.data.SMD$SMD*stress.data.SMD$sign.inversion.ours

stress.data.SMD.red <- stress.data.SMD[,c("esID","SMD.HE","SMD.ours","SMD.sv")]


# Heterocedasticity assumed: heteroscedastic population variances 
# in the two groups (Bonett, 2008, 2009).
stress.data.SMDH<-escalc(measure="SMDH", 
                         n1i=N.treat, n2i=N.control, 
                         m1i=mean.treat, m2i=mean.control, 
                         sd1i=SD.treat, sd2i=SD.control,
                         data=stress.data,
                         subset=stress.data$esID %in% es.list.SMD,
                         var.names=c("SMDH","SMDH.sv"), add.measure=FALSE,
                         append=TRUE)

stress.data.SMDH <- as.data.frame(stress.data.SMDH)
stress.data.SMDH$SMDH.HE <- stress.data.SMDH$SMDH*stress.data.SMDH$sign.inversion.HE
stress.data.SMDH$SMDH.ours <- stress.data.SMDH$SMDH*stress.data.SMDH$sign.inversion.ours

stress.data.SMDH.red <- stress.data.SMDH[,c("esID","SMDH.HE","SMDH.ours","SMDH.sv")]



##############################################################
# Calculating Hedge's g: shared control
##############################################################

es.list.SMD <- stress.data[!(is.na(stress.data$mean.control)),"esID"]

# Homocedasticity assumed
stress.data.SMD.sc<-escalc(measure="SMD", 
                           n1i=N.treat, n2i=N.control.sc, 
                           m1i=mean.treat, m2i=mean.control, 
                           sd1i=SD.treat, sd2i=SD.control,
                           data=stress.data,
                           subset=stress.data$esID %in% es.list.SMD,
                           var.names=c("SMD","SMD.sv"), add.measure=FALSE,
                           append=TRUE)

stress.data.SMD.sc <- as.data.frame(stress.data.SMD.sc)
stress.data.SMD.sc$SMD.sc.HE <- stress.data.SMD.sc$SMD*stress.data.SMD.sc$sign.inversion.HE
stress.data.SMD.sc$SMD.sc.ours <- stress.data.SMD.sc$SMD*stress.data.SMD.sc$sign.inversion.ours

colnames(stress.data.SMD.sc)[colnames(stress.data.SMD.sc)=="SMD.sv"] <- "SMD.sc.sv"

stress.data.SMD.sc.red <- stress.data.SMD.sc[,c("esID","SMD.sc.HE","SMD.sc.ours","SMD.sc.sv")]



# Heterocedasticity assumed: heteroscedastic population variances 
# in the two groups (Bonett, 2008, 2009).
stress.data.SMDH.sc<-escalc(measure="SMDH", 
                            n1i=N.treat, n2i=N.control.sc, 
                            m1i=mean.treat, m2i=mean.control, 
                            sd1i=SD.treat, sd2i=SD.control,
                            data=stress.data,
                            subset=stress.data$esID %in% es.list.SMD,
                            var.names=c("SMDH","SMDH.sv"), add.measure=FALSE,
                            append=TRUE)

stress.data.SMDH.sc <- as.data.frame(stress.data.SMDH.sc)
stress.data.SMDH.sc$SMDH.sc.HE <- stress.data.SMDH.sc$SMDH*stress.data.SMDH.sc$sign.inversion.HE
stress.data.SMDH.sc$SMDH.sc.ours <- stress.data.SMDH.sc$SMDH*stress.data.SMDH.sc$sign.inversion.ours

colnames(stress.data.SMDH.sc)[colnames(stress.data.SMDH.sc)=="SMDH.sv"] <- "SMDH.sc.sv"

stress.data.SMDH.sc.red <- stress.data.SMDH.sc[,c("esID","SMDH.sc.HE","SMDH.sc.ours","SMDH.sc.sv")]



##############################################################
# Calculating lnRR
##############################################################

es.list.lnRR <- stress.data[!(is.na(stress.data$mean.control)) & 
                              stress.data$ratioscale==1,"esID"]

stress.data.lnRR<-escalc(measure="ROM", 
                         n1i=N.treat, n2i=N.control, 
                         m1i=mean.treat, m2i=mean.control, 
                         sd1i=SD.treat, sd2i=SD.control,
                         data=stress.data,
                         subset=stress.data$esID %in% es.list.lnRR,
                         var.names=c("lnRR","lnRR.sv"), add.measure=FALSE,
                         append=TRUE)

stress.data.lnRR <- as.data.frame(stress.data.lnRR)
stress.data.lnRR$lnRR.HE <- stress.data.lnRR$lnRR*stress.data.lnRR$sign.inversion.HE
stress.data.lnRR$lnRR.ours <- stress.data.lnRR$lnRR*stress.data.lnRR$sign.inversion.ours

stress.data.lnRR.red <- stress.data.lnRR[,c("esID","lnRR.HE","lnRR.ours","lnRR.sv")]



##############################################################
# Calculating lnRR: shared control
##############################################################

es.list.lnRR <- stress.data[!(is.na(stress.data$mean.control)) & 
                              stress.data$ratioscale==1,"esID"]

stress.data.lnRR.sc<-escalc(measure="ROM", 
                            n1i=N.treat, n2i=N.control.sc, 
                            m1i=mean.treat, m2i=mean.control, 
                            sd1i=SD.treat, sd2i=SD.control,
                            data=stress.data,
                            subset=stress.data$esID %in% es.list.lnRR,
                            var.names=c("lnRR","lnRR.sv"), add.measure=FALSE,
                            append=TRUE)

stress.data.lnRR.sc <- as.data.frame(stress.data.lnRR.sc)
stress.data.lnRR.sc$lnRR.sc.HE <- stress.data.lnRR.sc$lnRR*stress.data.lnRR.sc$sign.inversion.HE
stress.data.lnRR.sc$lnRR.sc.ours <- stress.data.lnRR.sc$lnRR*stress.data.lnRR.sc$sign.inversion.ours

colnames(stress.data.lnRR.sc)[colnames(stress.data.lnRR.sc)=="lnRR.sv"] <- "lnRR.sc.sv"

stress.data.lnRR.sc.red <- stress.data.lnRR.sc[,c("esID","lnRR.sc.HE","lnRR.sc.ours","lnRR.sc.sv")]



##############################
# Finding outliers
##############################

# running some funnel plots to identify values that are very off

tm <- theme(panel.background = element_blank(),
            axis.line = element_line(size = 0.75),
            axis.text = element_text(size = 10, colour = "black"),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 16, hjust = 0.5))

ggplotly(ggplot(stress.data.cohens.red) + tm +
           geom_point(aes(x = cohens.biased.ours, y = 1/cohens.biased.sv, group = esID))) 

ggplotly(ggplot(stress.data.cohens.red) + tm +
           geom_point(aes(x = cohens.unbiased.ours, y = 1/cohens.biased.sv, group = esID))) 

ggplotly(ggplot(stress.data.SMD.red) + tm +
           geom_point(aes(x = SMD.ours, y = 1/SMD.sv, group = esID))) 

ggplotly(ggplot(stress.data.SMDH.red) + tm +
           geom_point(aes(x = SMDH.ours, y = 1/SMDH.sv, group = esID))) 

ggplotly(ggplot(stress.data.lnRR.red) + tm +
           geom_point(aes(x = lnRR.ours, y = 1/lnRR.sv, group = esID))) 


# when shared control is accounted for

ggplotly(ggplot(stress.data.cohens.sc.red) + tm +
           geom_point(aes(x = cohens.biased.sc.ours, y = 1/cohens.biased.sc.sv, group = esID))) 

ggplotly(ggplot(stress.data.cohens.sc.red) + tm +
           geom_point(aes(x = cohens.unbiased.sc.ours, y = 1/cohens.biased.sc.sv, group = esID))) 

ggplotly(ggplot(stress.data.SMD.sc.red) + tm +
           geom_point(aes(x = SMD.sc.ours, y = 1/SMD.sc.sv, group = esID))) 

ggplotly(ggplot(stress.data.SMDH.sc.red) + tm +
           geom_point(aes(x = SMDH.sc.ours, y = 1/SMDH.sc.sv, group = esID))) 

ggplotly(ggplot(stress.data.lnRR.sc.red) + tm +
           geom_point(aes(x = lnRR.sc.ours, y = 1/lnRR.sc.sv, group = esID))) 


# effect sizes that should be double-checked
esID.to.doublecheck <- c(141,142,173,176,305,823) #141,142,823 were due to typos that were later on fixed, that's why they no longer pops up

# # then
# # C(512,176,173,823)
# 
# stress.data <- stress.data[!(stress.data$esID %in% esID.to.exclude),]
# 
# stress.data.cohens.red.test <- stress.data.cohens.red[!(stress.data.cohens.red$esID %in% esID.to.exclude),]



##############################################################
# Calculating lnVR
##############################################################

es.list.lnVR <- stress.data[!(is.na(stress.data$mean.control)) & 
                              stress.data$prop.data==0,"esID"]

stress.data.lnVR<-escalc(measure="VR", 
                         n1i=N.treat, n2i=N.control, 
                         m1i=mean.treat, m2i=mean.control, 
                         sd1i=SD.treat, sd2i=SD.control,
                         data=stress.data,
                         subset=stress.data$esID %in% es.list.lnVR,
                         var.names=c("lnVR","lnVR.sv"), add.measure=FALSE,
                         append=TRUE)

stress.data.lnVR <- as.data.frame(stress.data.lnVR)
stress.data.lnVR.red <- stress.data.lnVR[,c("esID","lnVR","lnVR.sv")]

# adding the effect sizes for which metafor produces the following error:
# Warning message:
#   In escalc.default(measure = "VR", n1i = N.treat, n2i = N.control,  :
#                       Some 'yi' and/or 'vi' values equal to +-Inf. Recoded to NAs.
esID.to.doublecheck <- c(esID.to.doublecheck,stress.data.lnVR.red[is.na(stress.data.lnVR.red$lnVR),"esID"])



##############################################################
# Calculating lnVR:shared control
##############################################################

es.list.lnVR <- stress.data[!(is.na(stress.data$mean.control)) & 
                              stress.data$prop.data==0,"esID"]

stress.data.lnVR.sc<-escalc(measure="VR", 
                            n1i=N.treat, n2i=N.control.sc, 
                            m1i=mean.treat, m2i=mean.control, 
                            sd1i=SD.treat, sd2i=SD.control,
                            data=stress.data,
                            subset=stress.data$esID %in% es.list.lnVR,
                            var.names=c("lnVR","lnVR.sv"), add.measure=FALSE,
                            append=TRUE)

stress.data.lnVR.sc <- as.data.frame(stress.data.lnVR.sc)

stress.data.lnVR.sc.red <- stress.data.lnVR.sc[,c("esID","lnVR","lnVR.sv")]
names(stress.data.lnVR.sc.red) <- c("esID","lnVR.sc","lnVR.sc.sv")


##############################################################
# Calculating lnCVR
##############################################################

es.list.lnCVR <- stress.data[!(is.na(stress.data$mean.control)) & 
                               stress.data$prop.data==0 & 
                               stress.data$ratioscale==1,"esID"]

stress.data.lnCVR<-escalc(measure="CVR", 
                          n1i=N.treat, n2i=N.control, 
                          m1i=mean.treat, m2i=mean.control, 
                          sd1i=SD.treat, sd2i=SD.control,
                          data=stress.data,
                          subset=stress.data$esID %in% es.list.lnCVR,
                          var.names=c("lnCVR","lnCVR.sv"), add.measure=FALSE,
                          append=TRUE)

stress.data.lnCVR <- as.data.frame(stress.data.lnCVR)
stress.data.lnCVR.red <- stress.data.lnCVR[,c("esID","lnCVR","lnCVR.sv")]

# showing the effect sizes for which metafor produces the following error:
# Warning message:
#   In escalc.default(measure = "VR", n1i = N.treat, n2i = N.control,  :
#                       Some 'yi' and/or 'vi' values equal to +-Inf. Recoded to NAs.
esID.to.doublecheck <- c(esID.to.doublecheck,stress.data.lnCVR.red[is.na(stress.data.lnCVR.red$lnCVR),"esID"])



##############################################################
# Calculating lnCVR
##############################################################

es.list.lnCVR <- stress.data[!(is.na(stress.data$mean.control)) & 
                               stress.data$prop.data==0 & 
                               stress.data$ratioscale==1,"esID"]

stress.data.lnCVR.sc<-escalc(measure="CVR", 
                             n1i=N.treat, n2i=N.control.sc, 
                             m1i=mean.treat, m2i=mean.control, 
                             sd1i=SD.treat, sd2i=SD.control,
                             data=stress.data,
                             subset=stress.data$esID %in% es.list.lnCVR,
                             var.names=c("lnCVR","lnCVR.sv"), add.measure=FALSE,
                             append=TRUE)

stress.data.lnCVR.sc <- as.data.frame(stress.data.lnCVR.sc)

stress.data.lnCVR.sc.red <- stress.data.lnCVR.sc[,c("esID","lnCVR","lnCVR.sv")]
names(stress.data.lnCVR.sc.red) <- c("esID","lnCVR.sc","lnCVR.sc.sv")


##############################################################
# Double-checking outliers
##############################################################

esID.to.doublecheck.final <- sort(unique(esID.to.doublecheck))

# creating and exporting a dataset to make the process easier
# write.xlsx(stress.data[stress.data$esID %in% esID.to.doublecheck.final,],
#            "data_re-extraction/outlier_exploration/list_of_references_to_check_because_of_potential_outliers.xlsx",
#            sheetName="Sheet1",col.names=TRUE, row.names=F,
#            append=FALSE, showNA=TRUE, password=NULL)

# In case of SE/SD reported as 0 in the original publications,
# We are going to use the following SD imputation approach, which 
# is based on the following equation:

# SDj = Meanj * (sum(SDi)/sum(Meani)) 

# This approach assumes that log(SD)/log(Mean) = ca. sum(log(SDi))/sum(log(Meani))

#######################
# dealing with outliers

# esID==141:a typo that was fixed in the original dataset

# esID==142:a typo that was fixed in the original dataset

# stress.data[stress.data$esID==173 & stress.data$esID==176,]
# esID==173 and 176:the response variable is number of vertebrae,
# which seems to be between 58 and 59, and (1) n is large (176-180
# individuals), and the reported SE is very small (0.04-0.05).
# metafor warns about the inclusion of this data points as it 
# identifies them as extreme outliers. Therefore, we have decided
# to exclude them from the final dataset (see code at the bottom)

# esID==305:all values are as shown in original study, plus,
# this effect size is not that much off.

# esID==361:SE.control reported as 0.00 in the original paper.
# To take care of it, we use the means and SDs reported in the 
# paper for that specific response variable (i.e. Mass at 
# hatching (g)), see calculations below. Note: there is a
# value that is very off compare to the rest (mean=0.16,se=0.35)
# we are goint to exclude this one for the calculation.

# stress.data[stress.data$esID==361,"SD.control"] <- 0.14*((sum(se.to.sd(c(0.01,0.02,0.01,0.01),
#                                                                        c(41,30,32,32))))/
#                                                            (sum(c(0.15,0.14,0.15,0.16,0.17))))

# esID==546:SE.control reported as 0.00 in the original paper.
# To take care of it, we use the means and SDs reported in the 
# paper for that specific response variable (i.e. wattle chroma), 
# see calculations below. 

# stress.data[stress.data$esID==546,"SD.control"] <- 0.08*((sum(se.to.sd(c(0.01),
#                                                                        c(49))))/
#                                                            (sum(c(0.08))))

# esID==823:a typo that was fixed in the original dataset


##############################################################
# Adding all effect sizes
##############################################################

stress.data <- merge(stress.data,stress.data.cohens.red,
                     by="esID",all.x=TRUE)

stress.data <- merge(stress.data,stress.data.cohens.sc.red,
                     by="esID",all.x=TRUE)

stress.data <- merge(stress.data,stress.data.SMD.red,
                     by="esID",all.x=TRUE)

stress.data <- merge(stress.data,stress.data.SMDH.red,
                     by="esID",all.x=TRUE)

stress.data <- merge(stress.data,stress.data.SMD.sc.red,
                     by="esID",all.x=TRUE)

stress.data <- merge(stress.data,stress.data.SMDH.sc.red,
                     by="esID",all.x=TRUE)

stress.data <- merge(stress.data,stress.data.lnRR.red,
                     by="esID",all.x=TRUE)

stress.data <- merge(stress.data,stress.data.lnRR.sc.red,
                     by="esID",all.x=TRUE)

stress.data <- merge(stress.data,stress.data.lnVR.red,
                     by="esID",all.x=TRUE)

stress.data <- merge(stress.data,stress.data.lnVR.sc.red,
                     by="esID",all.x=TRUE)

stress.data <- merge(stress.data,stress.data.lnCVR.red,
                     by="esID",all.x=TRUE)

stress.data <- merge(stress.data,stress.data.lnCVR.sc.red,
                     by="esID",all.x=TRUE)


# exporting clean data for effect size estimation but 
# excluding effect sizes = NA (i.e. STATs based effect
# sizes). Also, we are excluding two effect sizes that
# were found to be extreme outliers (see lines:313-319)
write.xlsx(stress.data[!(is.na(stress.data$mean.control)) & stress.data$esID!=173 & stress.data$esID!=176,],
           "data_re-extraction/clean_data/EyckDev_stress_clean_effect_sizes.xlsx",
           sheetName="Sheet1",col.names=TRUE, row.names=F,
           append=FALSE, showNA=TRUE, password=NULL)


# saving session information with all packages versions for reproducibility purposes
sink("data_re-extraction/clean_data/effect_size_calculation_R_session.txt")
sessionInfo()
sink()
