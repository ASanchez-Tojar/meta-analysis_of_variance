
##############################################################
# Authors: 
# Alfredo Sanchez-Tojar (@ASanchez_Tojar)
# Profile: https://goo.gl/PmpPEB
# Department of Evolutionary Biology, Bielefeld University (GER) 
# Email: alfredo.tojar@gmail.com

# Script first created on the 2nd of July 2019

##############################################################
# Description of script and instructions
##############################################################

# This script is to identify studies in which a single control
# group is shared among different tests in the dataset used in:

# Eyck et al. 2019: Effects of developmental stress on animal
# phenotype and performance: a quantitative review

# Sharing a control is a form of non-independence that needs to 
# be dealt with. One way of dealing with this non-independence
# is to adjust sample sizes accordingly. The idea behind this
# is simple, to reduce the sample size of the control groups
# that are used multiple times so that they have less weight
# in the meta-analysis.


##############################################################
# Packages needed
##############################################################

pacman::p_load(openxlsx,dplyr)

# Clear memory
rm(list=ls())


##############################################################
# Functions needed
##############################################################

# none

##############################################################
# Importing dataset
##############################################################

# database coming from 003_data_preparation.R
db <- read.xlsx("data_re-extraction/clean_data/EyckDev_stress_clean_raw_data.xlsx",
                colNames=T,sheet = 1)


# removing rows without raw data as they will not be analyzed
db.red <- db[!(is.na(db$mean.control)),]

# # for testing
# test <- db.red[db.red$studyID==10 | db.red$studyID==100 | db.red$studyID==21,]


##############################################################
# Cleaning variables
##############################################################

# removing the extra space sometimes added by mistake
db.red$specific.trait <- str_remove(db.red$specific.trait, regex(" $"))


# list levels that potentially need some extra cleaning
# All good, they come from different studies
# [29] "beak colour score"                                                                                                                       
# [30] "beak score (colour, hue, brightness, saturation)"  

# From the same study, but their mean.control are different
# so it does not affect our process here
# [53] "cell mediated immunity"                                                                                                                  
# [54] "cell mediated immunity (difference between pre and post immunization)" 

# From the same study and with the same means, so reworded
# below
# [102] "fat content"                                                                                                                             
# [103] "fat content [%]"  

# From the same study, but their mean.control are different
# so it does not affect our process here
# [141] "humoral immunity"                                                                                                                        
# [142] "humoral immunity (difference between pre and post immunization)" 

# From the same study, but their mean.control are different
# so it does not affect our process here
# [371] "wattle chroma"                                                                                                                           
# [372] "wattle chroma (uv)"        

# From the same study, but their mean.control are different
# so it does not affect our process here
# [373] "wattle hue"                                                                                                                              
# [374] "wattle hue (visible)"


# substituting "fat content [%]" by "fat content"
db.red[db.red$specific.trait=="fat content [%]","specific.trait"] <- "fat content"


# # checking the variables related to the species studied
# unique(db.red[order(db.red$common.name),c("common.name","speciesID","scientific.name")])
# unique(db.red[order(db.red$scientific.name),c("common.name","speciesID","scientific.name")])
# unique(db.red[order(db.red$speciesID),c("common.name","speciesID","scientific.name")])

# studyID seems to be properly coded


##############################################################
# Generating new sample sizes
##############################################################

# the following code groups the database by studyID, then by 
# specific.trait and then by mean.control. That is, for each
# specifict.trait of a specific studyID, it checks if a 
# mean.control is shown more than once. If so, we are assuming 
# this is because of a control group is shared across different 
# comparisons. After having gone through the data base (see
# previous scripts), we know this is a fair assumption.
# The code also counts the number of times a mean.control is
# present, and then it calculates a new sample size "N.control.sc"
# that is equal to the original sample size divided by the number
# of times a specific control group was used.


#db.red.dup <- test %>% 
db.red.dup <- db.red %>% 
  group_by(studyID,speciesID,specific.trait,mean.control) %>% 
  mutate(num.shared.control = n()) %>%
  mutate(N.control.sc = N.control/num.shared.control)


# print(tbl_df(db.red.dup[,c("studyID","specific.trait","mean.control","num.shared.control","N.control","N.control.sc")]), 
#       n=nrow(db.red.dup))

# tibble to data.frame
db.red.dup.df <- as.data.frame(db.red.dup)
#db.red.dup.df[c(1:20),]


# additional typo found
db.red.dup.df[db.red.dup.df$TAXA=="Aves ","TAXA"] <- "Aves"


# # what's the number of studies using suffering from shared control?
# length(unique(db.red.dup.df[db.red.dup.df$num.shared.control>1,"studyID"]))
# 
# # what's the % then?
# round((length(unique(db.red.dup.df[db.red.dup.df$num.shared.control>1,"studyID"]))/length(unique(db.red.dup.df$studyID)))*100,1)


# exporting clean data for effect size estimation
write.xlsx(db.red.dup.df,
           "data_re-extraction/clean_data/EyckDev_stress_clean_raw_data_shared_control.xlsx",
           sheetName="Sheet1",col.names=TRUE, row.names=F,
           append=FALSE, showNA=TRUE, password=NULL)


# saving session information with all packages versions for reproducibility purposes
sink("data_re-extraction/clean_data/shared_control_R_session.txt")
sessionInfo()
sink()
