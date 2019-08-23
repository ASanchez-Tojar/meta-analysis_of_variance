
##############################################################
# Authors: 
# Alfredo Sanchez-Tojar (@ASanchez_Tojar)
# Profile: https://goo.gl/PmpPEB
# Department of Evolutionary Biology, Bielefeld University (GER) 
# Email: alfredo.tojar@gmail.com

# Input from:
#
# Nicholas P. Moran
# Profile: https://www.researchgate.net/profile/Nicholas_Moran
# Department of Evolutionary Biology, Bielefeld University (GER) 
# Email: nicholaspatrickmoran@gmail.com

# Script first created on the 18th of April 2019

##############################################################
# Description of script and Instructions
##############################################################

# This script is to double-check and re-extract the data used in:

# Eyck et al. 2019: Effects of developmental stress on animal
# phenotype and performance: a quantitative review

# Alfredo Sanchez-Tojar and Nicholas P. Moran did a pilot 
# re-extraction of some of the studies in Jan 2019. This
# pilot re-extraction revealed some inconsistencies. Thus,
# we decided to double-check and re-extract the all the 
# raw data used in that study

# This script was originally designed with the idea that we
# did not need to double-check all the raw data, which explains
# why we originally used randomization to assign the studies
# and what order. However, at the end, we decided to
# double-check all the studies.

# The output of this script were later used to double-check and 
# re-extract data, adding this way some additional variables,
# and renaming some of the existing variables. Thus, next script
# ("003_data_preparation.R"), uses those versions of the datasets,
# named as "*_EDITED.xlsx".

##############################################################
# Packages needed
##############################################################

pacman::p_load(openxlsx)

# Clear memory
rm(list=ls())


##############################################################
# Functions needed
##############################################################

# none

##############################################################
# Assigning re-extracting efforts
##############################################################

#####################
# Preparing the data
#####################

pilot <- read.xlsx("data_re-extraction/re-extracted/re-extracing_data_Eyck_early-life_stress_pilot.xlsx",
                   colNames=T,sheet = 2)


# study ID of studies re-extracted so far
rextracted <- sort(unique(pilot$studyID))


# updated dataset sent by H.Eyck
stress.red <- read.xlsx("data/EyckDev.stress_Data_FULL_TABLE.xlsx",
                        colNames=T,sheet = 1)

# identifying if we re-extracted any study that is not
# included in the final dataset that we will re-analyse
missing <- setdiff(pilot$studyID,stress.red$studyID)


# excluding this study from the list of studies that
# we re-extracted to keep everything under control,
# clean and sensical
rextracted.adj <- rextracted[-match(c(missing),rextracted)]


################################
# Assigning re-extraction order
################################

# full list of studies
study.ids <- sort(unique(stress.red$studyID))


# list of studies that have not been re-extracted yet
study.ids.adj <- study.ids[!(study.ids %in% rextracted.adj)]
#length(setdiff(rextracted.adj,study.ids.adj)) #as expected

# we have decided to only re-extract data from studies for which
# raw data (i.e. means, var, n) is presented. The reason is that
# our analyses will only include this data because lnRR can only
# be estimated from raw data (i.e. not from test statistics)
test.stastics <- sort(unique(stress.red[(is.na(stress.red$mean.control)) | 
                                          (is.na(stress.red$mean.treat)),c("studyID")]))


# this approach, however, excludes studies were both raw data and statistics
# have been extracted, which means that, to re-extract all raw data
# used in the study we need to re-extract some additional data (see bottom
# of the script).

study.ids.adj.2 <- study.ids.adj[!(study.ids.adj %in% test.stastics)]


# choosing the order of studies randomly
# this was because, originally, we thought we did not need to re-extract
# data from all studies. However, we were not comfortable with the rate
# of typos found and decided to re-extract all raw data used in the 
# original study.
set.seed(111)

list.studies <- sample(study.ids.adj.2,length(study.ids.adj.2))

nick <- list.studies[1:30]
alfredo <- list.studies[31:61]


####################
# Preparing datasets
####################

####################
# 1st round
####################

db.nick <- stress.red[stress.red$studyID %in% nick,]
db.alfredo <- stress.red[stress.red$studyID %in% alfredo,]


# exporting datasets:
# keep in mind that we added some new variables to this database
# during the re-extraction processand and also renamed some of 
# the original variables
write.xlsx(db.nick, "data_re-extraction/EyckDev.stress_Data_FULL_TABLE_Nick_re-extraction_subset.xlsx",
           sheetName="Sheet1",col.names=TRUE, row.names=F,
           append=FALSE, showNA=TRUE, password=NULL)

write.xlsx(db.alfredo, "data_re-extraction/EyckDev.stress_Data_FULL_TABLE_Alfredo_re-extraction_subset.xlsx",
           sheetName="Sheet1",col.names=TRUE, row.names=F,
           append=FALSE, showNA=TRUE, password=NULL)


####################
# 2nd round
####################

# subset of studies that contained both raw data and test statistics

done.so.far <- unique(c(study.ids.adj.2,rextracted.adj))

db.both.data <- stress.red[!(stress.red$studyID%in%done.so.far),]
  
# effect sizes that still need to be re-extractedS         
db.both.data.1 <- db.both.data[!(is.na(db.both.data$mean.control)),]


# splitting work between NPM and AST
nick.2 <- unique(db.both.data.1$studyID)[1:12]
alfredo.2 <- unique(db.both.data.1$studyID)[13:21]


db.nick.2 <- db.both.data.1[db.both.data.1$studyID %in% nick.2,]
db.alfredo.2 <- db.both.data.1[db.both.data.1$studyID %in% alfredo.2,]

write.xlsx(db.nick.2, "data_re-extraction/EyckDev.stress_Data_FULL_TABLE_Nick_re-extraction_subset_2.xlsx",
           sheetName="Sheet1",col.names=TRUE, row.names=F,
           append=FALSE, showNA=TRUE, password=NULL)

write.xlsx(db.alfredo.2, "data_re-extraction/EyckDev.stress_Data_FULL_TABLE_Alfredo_re-extraction_subset_2.xlsx",
           sheetName="Sheet1",col.names=TRUE, row.names=F,
           append=FALSE, showNA=TRUE, password=NULL)


####################
# Kept as original
####################

# effect sizes that do not need to be re-extracted
db.both.data.2 <- db.both.data[(is.na(db.both.data$mean.control)),]

write.xlsx(db.both.data.2, "data_re-extraction/EyckDev.stress_Data_FULL_TABLE_statistics_re-extraction_not_needed.xlsx",
           sheetName="Sheet1",col.names=TRUE, row.names=F,
           append=FALSE, showNA=TRUE, password=NULL)


# final number test
(nrow(pilot[pilot$studyID %in% rextracted.adj,])+nrow(db.nick)+nrow(db.alfredo)+nrow(db.nick.2)+nrow(db.alfredo.2)+nrow(db.both.data.2))==(nrow(stress.red)) # matches


# saving the reduced pilot dataset
write.xlsx(pilot[pilot$studyID %in% rextracted.adj,], 
           "data_re-extraction/re-extracing_data_Eyck_early-life_stress_pilot_reduced.xlsx",
           sheetName="Sheet1",col.names=TRUE, row.names=F,
           append=FALSE, showNA=TRUE, password=NULL)


# saving session information with all packages versions for reproducibility purposes
sink("data_re-extraction/data_re-extraction_R_session.txt")
sessionInfo()
sink()
