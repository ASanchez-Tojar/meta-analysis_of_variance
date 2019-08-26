##############################################################
# Authors: 
#
# Alfredo Sanchez-Tojar (@ASanchez_Tojar)
# Profile: https://goo.gl/PmpPEB
# Bielefeld University
# Email: alfredo.tojar@gmail.com

##############################################################
# Description of script and Instructions
##############################################################

# Script first created the 19th of Aug 2019

# This script is to perform a few extra calculations and create
# some additional plots for our case study.


##############################################################
# Packages needed
##############################################################

# load pacakges
pacman::p_load(dplyr,revtools,ggplot2,tidyr,metafor,openxlsx)

# cleaning up
rm(list=ls())


##############################################################
# Functions needed
##############################################################

# none


##############################################################
# Effect size overview: 
##############################################################

# We use the data from:

# Senior et al. 2016. Heterogeneity in ecological and 
# evolutionary meta-analyses: its magnitude and implications

# to explore what type of effect sizes are used across
# meta-analyses in Ecology and Evolution.


# importing the database published with Senior et al. 2016
db <- read.table("literature_review/Senior_et_al_2016/Data_Package_Part_3_Survey_Data.csv",header=T,sep=",")

# exploring the frequency of each effect size in this database
table(db$Effect.Size)

db$Effect.Size <- as.character(db$Effect.Size)

# recategorizing effect size type following Nakagawa and Santos 2012
db$effect.size.v2 <- ifelse(db$Effect.Size=="beta" |
                              db$Effect.Size=="CV" |
                              db$Effect.Size=="diff" |
                              db$Effect.Size=="exponent" |
                              db$Effect.Size=="mean" |
                              db$Effect.Size=="proportion", 
                            "Others",
                            db$Effect.Size)

table(db$effect.size.v2)

db <- db %>% mutate(effect.size.v2 = factor(effect.size.v2, levels = c("d","lnRR","Zr","r","lnOR","RR", "Others")))

# one can see that the majority of effect sizes are based on means
# (i.e. d + lnRR vs. Zr + r... etc)
table(db$effect.size.v2)


##############################################################
# Meta-analysis of variance in ecology and evolution
##############################################################

# We perform a systematic review to understand the extend of
# the use of meta-analysis of variance in ecology and evolution
# since Nakagawa et al. 2015 (Methods in Ecology and Evolution).
# For that, we performed a search to find all the references
# citing Nakagawa et al. 2015 according to Web of Science 
# (date: 20th Aug, 2019). This script import the results of that
# search, which are then formatted for the title-and-abstract
# screening software Rayyan. Additionally, the results of the 
# title-and-abstract and full-text screening are analyzed to
# provide a temporal overview of the use of meta-analysis of 
# variance in ecology and evolution.

# importing references
citing.refs <- read_bibliography("literature_review/Nakagawa_et_al_2015/Nakagawa_et_al_2015_citations_WoS_20190820.bib")

# reducing fields to the minimum number of fields
# so that all databases have the same columns. Also, these fields
# are the important ones for the screening (see below).
reducing.fields.wos <- c("label","title","author","journal","issn","volume",
                         "number","pages","year","publisher","doi","abstract") #number to issue, 

citing.refs.red <- citing.refs[,reducing.fields.wos]


#################################
# Formatting data for RAYYAN QCRI

# choose only the fields needed for creating a .csv file importable by: https://rayyan.qcri.org

# example of a valid .csv file. The fields are the following:
# key,title,authors,journal,issn,volume,issue,pages,year,publisher,url,abstract
names.rayyan <- c("key","title","authors","journal","issn","volume","issue","pages","year","publisher","url","abstract")
names.rayyan
names(citing.refs.red)

# standardizing fields according to rayyan.example

# what's different between the two?
setdiff(names.rayyan,names(citing.refs.red))
setdiff(names(citing.refs.red),names.rayyan)

citing.refs.red.rayyan <- plyr::rename(citing.refs.red, c("label"="key", "author"="authors", "number"="issue", "doi"="url"))
names(citing.refs.red.rayyan)

setdiff(names.rayyan,names(citing.refs.red.rayyan))
setdiff(names(citing.refs.red.rayyan),names.rayyan)

# reorder
citing.refs.red.rayyan <- citing.refs.red.rayyan[,names.rayyan]


##################
# Creating output
write.csv(citing.refs.red.rayyan[order(citing.refs.red.rayyan$title),],"literature_review/Nakagawa_et_al_2015/Nakagawa_et_al_2015_citing_references_rayyan.csv",row.names=FALSE)
#remember to manually remove the quotes for the column names only in the .csv file before importing into rayyan


#################################
# Exploring the temporal pattern

# importing the extracted data after fulltext screening
data.ref <- read.table("literature_review/Nakagawa_et_al_2015/Nakagawa_et_al_2015_citing_references_fulltext_screening.csv",
                       header=T,sep=",")

summary(data.ref)

# counting numbers
nrow(data.ref)
table(data.ref$t.and.a_decision)
table(data.ref$fulltext_decision)

# subsetting those included in the database
data.ref.included <- data.ref[data.ref$fulltext_decision=="yes" |
                                data.ref$fulltext_decision=="yes_but_no_lnCVR",]

data.ref.included <- data.ref.included[!(is.na(data.ref.included$fulltext_decision)),]

nrow(data.ref.included)

# counting number of studies per year
barplot <- data.ref.included %>% 
  group_by(year) %>% 
  summarise(count = n())

barplot

# data.ref.included$authors # a few papers involving Shinichi, so the real extend should be perhaps consider overestimated

# plotting the number of meta-analyses of variance per year in
# ecology and evolution as per our review
tiff("plots/FigureS1_meta-analysis_variance_numbers.tiff",
     height=10, width=10,
     units='cm', compression="lzw", res=800)

p<-ggplot(data=barplot, aes(x=year, y=count)) +
  geom_bar(stat="identity", fill="steelblue")+
  ylab("Number of meta-analyses of variance")+
  xlab("Year")+
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
        axis.title.x = element_text(hjust = 0.5, size = 13, face = "bold"),
        axis.title.y = element_text(size = 13, hjust = 0.5, margin = margin(r=8),face = "bold"),
        axis.text.y = element_text(angle = 0,color="black",hjust=0),
        axis.text.x = element_text(color="black"),
        plot.title = element_text(size = 16))
p

dev.off()


##############################################################
# Quick exploration of directionless effect sizes
##############################################################

# How many studies were affected by directionless effect sizes?

inversion.unclear.2 <- c("(dihydroxyphenylacetic acid/dopamine) ratio in the nuclues accumbens after observation of a receptive stimulus female", # studyID=107:authors do not have a clear prediciton about this measurement but argue that DOPAC/DA correlates with dopaminergic turnover, which could help Tf males better compete with Tm males                   
                         "(dihydroxyphenylacetic acid/dopamine) ratio in the preoptic area after observation of a receptive stimulus female", # studyID=107:authors do not have a clear prediciton about this measurement but argue that DOPAC/DA correlates with dopaminergic turnover, which could help Tf males better compete with Tm males                       
                         "(dihydroxyphenylacetic acid/dopamine) ratio in the ventral tegumental area after observation of a receptive stimulus female", # studyID=107:authors do not have a clear prediciton about this measurement but argue that DOPAC/DA correlates with dopaminergic turnover, which could help Tf males better compete with Tm males
                         "adrenal weight to body weight ratio", # studyID=37:the authors do not set predictions for any of the ratios. However, we will assume that brain and gonad ratios should be negatively correlated with stress. The remaining ratios will be part of the unclear category
                         "heart weight to body weight ratio", # studyID=37:the authors do not set predictions for any of the ratios. However, we will assume that brain and gonad ratios should be negatively correlated with stress. The remaining ratios will be part of the unclear category
                         "hepatosomatic index (100*liver wet weight/body wet weight)", # studyID=7:see comment for "mean liver glycogen"
                         "hippocampus weight to brain weight ratio", # studyID=37:the authors do not set predictions for any of the ratios. However, we will assume that brain and gonad ratios should be negatively correlated with stress. The remaining ratios will be part of the unclear category
                         "kidney weight to body weight ratio", # studyID=37:the authors do not set predictions for any of the ratios. However, we will assume that brain and gonad ratios should be negatively correlated with stress. The remaining ratios will be part of the unclear category
                         "lung weight to body weight ratio", # studyID=37:the authors do not set predictions for any of the ratios. However, we will assume that brain and gonad ratios should be negatively correlated with stress. The remaining ratios will be part of the unclear category
                         "mean liver glycogen", # studyID=7:authors do not have a prediction regarding this measurement. Also, in the discussion (page 6 (178)) they estated that it is indeed unclear what to expect as many studies find contrasting results.
                         "mean selected temperature", # studyID=101: authors do not have a prediction regarding this measurement, and this variable only seem to make sense in the light of whether incubation (two levels) and "mean selected temperature" match or not. Definition: "Position data were then converted to selected air temperatures by substituting position into the appropriate regression equation, and correcting for minor differences among lanes. From these data, we calculated each individual's mean selected temperature" 
                         "percentage self grooming during aggression test", # studyID=61:the authors do not provide any prediction or information regarding this measurement
                         "pituitary weight to body weight ratio", # studyID=37:the authors do not set predictions for any of the ratios. However, we will assume that brain and gonad ratios should be negatively correlated with stress. The remaining ratios will be part of the unclear category
                         "thermal preference", # studyID=103:same methodology than in studyID=101, and again, authors do not have a prediction regarding this measurement
                         "thermal preference ", # studyID=103:same methodology than in studyID=101, and again, authors do not have a prediction regarding this measurement
                         "thermoregularity precision", # studyID=101: it is a mean of SDs that should not be included? ("standard deviation of their 17 temperature observations (as an index of the precision with which they maintained this mean temperature)")
                         "total granule cell number in brain",  # studyID=59:the authors do not provide any predictions regarding this variable
                         "volume of granular cell layer in brain", # studyID=59:the authors do not provide any predictions regarding this variable
                         "volume of molecular cell layer in brain" # studyID=59:the authors do not provide any predictions regarding this variable
)

studyIDs <- c(107,107,107,37,37,7,37,37,37,7,101,61,37,103,103,101,59,59,59)

length(unique(studyIDs))


# saving versions used for reproducibility purposes
sink("literature_review/Nakagawa_et_al_2015/metaanalysis_variance_review_Nakagawa_et_al_2015_Rpackages_session.txt")
sessionInfo()
sink()
