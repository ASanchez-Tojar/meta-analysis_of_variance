
##############################################################
# Authors: 
# Alfredo Sanchez-Tojar (@ASanchez_Tojar)
# Profile: https://goo.gl/PmpPEB
# Department of Evolutionary Biology, Bielefeld University (GER) 
# Email: alfredo.tojar@gmail.com

# Script first created on the 26th of April 2019

##############################################################
# Description of script and instructions
##############################################################

# This script is to put together the re-extraced/double-checked
# data used in:

# Eyck et al. 2019: Effects of developmental stress on animal
# phenotype and performance: a quantitative review

# Alfredo Sanchez-Tojar and Nicholas P. Moran double-checked 
# all the raw data (i.e. means, sds and n; excluding effect 
# sizes based on F, t or chi-squared) shown in the table that
# H.Eyck shared with us (i.e. EyckDev.stress_Data_FULL_TABLE),
# which is the table that we will use for the analyses since
# it contains the raw data of the selection of studies that 
# were used in the meta-analysis (note that there are 865 total
# effect sizes in this dataset, in contrast to the published
# data set, which contained 866 effect sizes).

# There were different rounds of double-checking and re-extraction:

# 1. Pilot: 10 studies from the published full database (i.e
#           brv12496-sup-0004-tables4.xls). Data re-extracted by 
#           AST and NPM in Jan 2019.

# 2. NPM round 1: 30 studies from the table that H.Eyck shared 
#                 with us (i.e. EyckDev.stress_Data_FULL_TABLE).
.
# 3. AST round 1: 31 studies from the table that H.Eyck shared 
#                 with us (i.e. EyckDev.stress_Data_FULL_TABLE).

# 4. NPM round 2: 12 studies from the table that H.Eyck shared 
#                 with us (i.e. EyckDev.stress_Data_FULL_TABLE).

# 5. AST round 2: 9 studies from the table that H.Eyck shared 
#                 with us (i.e. EyckDev.stress_Data_FULL_TABLE).

# The remaining effect sizes come from F, t or chi-square tests.
# We did not double-checked those as they cannot be used to
# calculate lnRR and relatives, which is the aim of our re-analysis.


##############################################################
# Packages needed
##############################################################

pacman::p_load(openxlsx,plyr,stringr)

# Clear memory
rm(list=ls())


##############################################################
# Functions needed
##############################################################

# none


##############################################################
# Importing datasets
##############################################################

# database with the corrected data from our pilot re-extraction
db.pilot <- read.xlsx("data_re-extraction/re-extracted/re-extracing_data_Eyck_early-life_stress_pilot_reduced_EDITED.xlsx",
                      colNames=T,sheet = 1)


# databases with the corrected data from AST re-extractions
db.AST.1 <- read.xlsx("data_re-extraction/re-extracted/EyckDev.stress_Data_FULL_TABLE_Alfredo_re-extraction_subset_EDITED.xlsx",
                      colNames=T,sheet = 1)

db.AST.2 <- read.xlsx("data_re-extraction/re-extracted/EyckDev.stress_Data_FULL_TABLE_Alfredo_re-extraction_subset_2_EDITED.xlsx",
                      colNames=T,sheet = 1)


# databases with the corrected data from NPM re-extractions
db.NPM.1 <- read.xlsx("data_re-extraction/re-extracted/EyckDev.stress_Data_FULL_TABLE_Nick_re-extraction_subset_EDITED.xlsx",
                      colNames=T,sheet = 1)

db.NPM.2 <- read.xlsx("data_re-extraction/re-extracted/EyckDev.stress_Data_FULL_TABLE_Nick_re-extraction_subset_2_EDITED.xlsx",
                      colNames=T,sheet = 1)


# database with the non-double-checked test statistic

db.statistics <- read.xlsx("data_re-extraction/re-extracted/EyckDev.stress_Data_FULL_TABLE_statistics_re-extraction_not_needed_EDITED.xlsx",
                           colNames=T,sheet = 1)


##############################################################
# Subsetting datasets
##############################################################

# db.statistics contains all the variables we would like to keep.
# We will subset the remaining databases according to the exact
# selection of variables we want. 

names(db.statistics)

# Some checks before we subset:

# AST databases: all good
table(names(db.AST.1)==names(db.AST.2))

# NPM databases: all good
table(names(db.NPM.1)==names(db.NPM.2))

# AST vs. NPM: all good
setdiff(names(db.AST.1),names(db.NPM.1))
setdiff(names(db.NPM.1),names(db.AST.1))

# are the variables from db.statistics in AST and NPM? YES
setdiff(names(db.statistics),names(db.AST.1))
setdiff(names(db.statistics),names(db.NPM.1))

# are the variables from db.statistics in pilot? NO, we will subset it accordingly
setdiff(names(db.statistics),names(db.pilot))
setdiff(names(db.pilot),names(db.statistics))

# Then, deciding on which variables to keep/re-name before
# compiling all the databases in a single database. We only
# need the following list, so we will subset the remaining
# databases according to this list of databases

var.to.exclude <- setdiff(names(db.statistics),names(db.pilot))

var.to.include <- names(db.statistics[,-which(names(db.statistics) %in% var.to.exclude)])


# Actual subsetting
db.statistics.red <- db.statistics[,var.to.include]
db.AST.1.red <- db.AST.1[,var.to.include]
db.AST.2.red <- db.AST.2[,var.to.include]
db.NPM.1.red <- db.NPM.1[,var.to.include]
db.NPM.2.red <- db.NPM.2[,var.to.include]
db.pilot.red <- db.pilot[,var.to.include]


# Checking the order of the variables, as a likely unncessary double-check
table(names(db.statistics.red)==names(db.AST.1.red))
table(names(db.statistics.red)==names(db.AST.2.red))
table(names(db.statistics.red)==names(db.NPM.1.red))
table(names(db.statistics.red)==names(db.NPM.2.red))
table(names(db.statistics.red)==names(db.pilot.red))


##############################################################
# Assembling datasets
##############################################################

# Creating a variable to keep track of data subsets
db.statistics.red$origin <- "STATS"
db.AST.1.red$origin <- "AST"
db.AST.2.red$origin <- "AST"
db.NPM.1.red$origin <- "NPM"
db.NPM.2.red$origin <- "NPM"
db.pilot.red$origin <- "PILOT"


# assembling dataset
db.full <- rbind(db.statistics.red,
                 db.AST.1.red,db.AST.2.red,
                 db.NPM.1.red,db.NPM.2.red,
                 db.pilot.red)


##############################################################
# Subsetting and creating variables
##############################################################

# transforming some variables into factors
var.factor<- c("TAXA","scientific.name","developmental.stressor","method","exposure","sex",
               "trait.class","age","specific.trait","origin")

db.full[,var.factor] <- lapply(db.full[,var.factor] , factor)


# Reducing the number of variables to those that we are interested and those moderators without NA's.
# Many additional variables were not kept in the dataset because it became clear that
# they were modified from this version of the dataset to the final one published (table 2).
# For example, some variables show NA's in this dataset, some levels are miss-labelled, etc

subset.variables <- c("studyID","speciesID","TAXA","common.name","scientific.name","citation",
                      "location","lat","long","av..temp","htmnth","cldmnth",
                      "seasonality","developmental.stressor",
                      "sex","age","trait.class","specific.trait", #note that age was added later, which is why it is not present in list_of_references_to_check_for_group-level_prop_data.xlsx, however, this does not matter for the purposes of that dataset, which is simply double-checking data type
                      "mean.control","SD.control","N.control",
                      "mean.treat","SD.treat","N.treat",
                      "t.value","F.value","origin")


db.full.red <- db.full[,subset.variables]


# some cleaning of moderators
levels(db.full.red$trait.class)
db.full.red$trait.class <- revalue(db.full.red$trait.class, c("behavioural "="behavioural"))
levels(db.full.red$trait.class)


levels(db.full.red$sex)
db.full.red$sex <- revalue(db.full.red$sex, c("Female"="female", "females"="female"))
db.full.red$sex <- revalue(db.full.red$sex, c("Male"="male", "male "="male","males"="male"))
levels(db.full.red$sex)


# creating an effect size ID
db.full.red$esID <- 1:nrow(db.full.red)


##############################################################
# Identifying ratio scale data
##############################################################

# This is because lnRR and lnCVR can only be estimated for ratio
# scale data. In contrast, lnVR could be estimated for non-ratio 
# scale data too.

# Our strategy to find about ratio scale data is to simply assume
# that negative mean values are non-ratio scale data, whereas the
# remaining mean values are ratio scale data.

db.full.red$ratioscale <- ifelse(db.full.red$mean.control<0 | db.full.red$mean.treat<0,
                                 0,
                                 1)

#db.full.red[db.full.red$mean.control<0 | db.full.red$mean.treat<0,c("mean.control","mean.treat","ratioscale")]

# in none of the effect sizes identified as non-ratio scale data
# there were cases where both means were negative
#db.full.red[db.full.red$mean.control<0 & db.full.red$mean.treat<0,]

# Estimating the number and % of effect sizes that are ratio vs. 
# non-ratio scale data

# total numbers
table(db.full.red$ratioscale)

# percentages
round((table(db.full.red$ratioscale)/sum(table(db.full.red$ratioscale)))*100,1)


##############################################################
# Identifying group level proportion data
##############################################################

# This is because lnVR and lnCVR can only be estimated for 
# non-group level proportion data

# Our strategy to find ratio scale data consisted of going 
# through all the levels of the variable "specific.trait" and 
# flagged all those that could potentially refer to group level 
# proportion data. Those identified as potential group level 
# proportion data were then checked in the original papers to
# confirm data type.

# Note that levels for the subset origin==STATS were not checked
# as they cannot be used to calculate lnRR, lnVR or lnCVR.
# First, AST went through all the levels and split them into two
# groups: non.prop.data and unknown.data. Second, NPM went through
# those two groups to also give opinion on which variables should
# be checked in the original papers to confirm data type. NPM
# suggested the following additional levels:
# "dominance interactions won", #it is probably per individual, so fine
# "number of individuals defeated in domince display", #it is probably per individual, so fine
# "thermal preference", #potentially coded as binomial
# "thermal preference ", #potentially coded as binomial
# "visiting frequency", # could be a treatment group level frequency
# "cling to surrogate", #don't know what this means
# "number of assymmetric characters", #don't know what this means
# All those suggestions were included in the unknown.data vector and
# subsequently checked in the original papers.

levels(factor(db.full.red[db.full.red$origin!="STATS",c("specific.trait")]))


# the following vector contains a list of levels of 
# "specific.trait" that can 'safely' be noted as non-proportional
# data. Note that AST first developed the list, and NPM double-checked it
non.prop.data <- c("(dihydroxyphenylacetic acid/dopamine) ratio in the nuclues accumbens after observation of a receptive stimulus female",
                   "(dihydroxyphenylacetic acid/dopamine) ratio in the preoptic area after observation of a receptive stimulus female",
                   "(dihydroxyphenylacetic acid/dopamine) ratio in the ventral tegumental area after observation of a receptive stimulus female",
                   "abdominal length",
                   "abdominal width",
                   "activity level (lines crossed)",
                   "adrenal weight to body weight ratio",
                   "adrenocorticotropic hormone response to anesthesia",
                   "adrenocorticotropic hormone response to stress challenge",
                   "age at death",
                   "age at first mating ",
                   "age at first oviposition",
                   "age at metamorphosis",
                   "age at vaginal opening (puberty) ",
                   "age of first oestrus ",
                   "age post maturation at mating",
                   "aggression towards conspecific",
                   "anal lepidotrichia",
                   "arm length",
                   "arm length ",
                   "aterial pressure",
                   "attempted mounts",
                   "average clutch size",
                   "average flight time under flight stress",
                   "average hatch success per female",
                   "average number of deformed fish per temperature treatment",
                   "average sprint speed",
                   "average total eggs per female",
                   "basal metabolic rate",
                   "baseline cort",
                   "beak colour score",
                   "beak score (colour, hue, brightness, saturation)",
                   "bill coloration length",
                   "bill hue",
                   "bites to mirror",
                   "body depth",
                   "body mass",
                   "body mass ",
                   "body mass at hatching ",
                   "body mass at metamorphosis",
                   "body mass change during first immune challenge",
                   "body shape at hatching (mass/snout-vent length)",
                   "body temperature",
                   "body weight",
                   "body weight ",
                   "body weight at 20 days of gestation (1 year postnatal)",
                   "body weight at 20 days of gestation (150 days postnatal)",
                   "body weight at first oestrus",
                   "body weight at vaginal opening ",
                   "brain weight to body weight ratio",
                   "breeding time in natural conditions",
                   "breeding time in supplemented conditons",
                   "brood fledgling mass",
                   "caudal fin width",
                   "caudal peduncle width",
                   "cephalothorax area",
                   "cephalothorax length",
                   "cephalothorax width",
                   "cerebellum protein carbonyls",
                   "cerebellum superoxide dismutase",
                   "change in length after 1 year",
                   "circulating oestradiol",
                   "climbing",
                   "clutch mass",
                   "clutch size",
                   "clutch size in poor environment",
                   "cognitive performance",
                   "comb and wattle size",
                   "condition factor (100*g cm-3)",
                   "cort 1hr after swimming water maze",
                   "cort 30 mins after handling stress",
                   "corticosterone concentration",
                   "cortisol response to anesthesia",
                   "cortisol response to stress challenge",
                   "culmen length",
                   "cycle length at 140 days ",
                   "cylcle length at 1 year",
                   "days to learn associative task",
                   "days until pupation",
                   "deformed vertebrae",
                   "density of pyramidal neurons in left hipocampus",
                   "development duration",
                   "development rate",
                   "development time",
                   "dihydroxyphenylacetic acid in the nuclues accumbens after observation of a receptive stimulus female",
                   "dihydroxyphenylacetic acid in the preoptic area after observation of a receptive stimulus female",
                   "distance travelled during test session",
                   "dopamine in the nuclues accumbens after observation of a receptive stimulus female",
                   "dopamine in the preoptic area after observation of a receptive stimulus female",
                   "dopamine in the ventral tegumental area after observation of a receptive stimulus female",
                   "dorsal lepidotrichia",
                   "dorsal spines",
                   "dry body mass",
                   "dry weight",
                   "ear tuft",
                   "egg mass",
                   "egg size",
                   "epaxial muscle depth",
                   "exploration",
                   "exploratory behaviour",
                   "eye depth",
                   "fat content",
                   "fat content [%]",
                   "feather mass",
                   "feather tip length",
                   "fecal pellets in response to novel stress",
                   "feeding latency",
                   "female activity, calling and tail quivering in response to male",
                   "femur length",
                   "femur width",
                   "final instar development time ",
                   "forebrain volume",
                   "forebrain weight ",
                   "forehead colour patch area",
                   "forewing length",
                   "fork length",
                   "gape length",
                   "geniohyoideus length",
                   "gonad weight to body weight ratio",
                   "gonopodium length",
                   "granulocyte:lymphocyte ratio",
                   "growth rate",
                   "habenula size",
                   "haematocrit",
                   "haematocrit (ratio of red blood cells to whole blood volume) in natural conditions",
                   "haematocrit (ratio of red blood cells to whole blood volume) in supplemented conditons",
                   "head and bill length",
                   "head length",
                   "head length ",
                   "head width",
                   "head width ",
                   "heart rate",
                   "heart weight to body weight ratio",
                   "hematocrit (% packed cell volume)",
                   "hemibranch area",
                   "hemibranch perimeter",
                   "hepatosomatic index (100*liver wet weight/body wet weight)",
                   "hippocampus weight to brain weight ratio",
                   "hopping distance",
                   "hours until death by starvation following metamorphosis",
                   "hours until metamorphosis after hatching",
                   "humerus length",
                   "humoral immunity",
                   "humoral immunity (difference between pre and post immunization)",
                   "incubation length",
                   "incubation period",
                   "incubation period ",
                   "initial activity",
                   "interferone (ifn)-y (cytokine) concentration in blood", #annoying character in the name
                   "interleukin 10 (cytokine) concentration in blood",
                   "interleukin 4 (cytokine) concentration in blood",
                   "interleukin 6 (cytokine) concentration in blood",
                   "key pecks during song preference tests of familiar songs",
                   "key pecks during song preference tests of unfamiliar songs",
                   "kidney weight to body weight ratio",
                   "latency (s) to leave starting box and enter maze",
                   "latency to escape",
                   "latency to exit begin exploring",
                   "latency to feed from enclosed novel environment",
                   "latency to find food reward after training ",
                   "latency to first appearence in maze",
                   "latency to first egg sac",
                   "latency to immobility in forced swim test",
                   "latency to lay",
                   "latency to lay eggs after pairing",
                   "latency to leave starting box and enter maze",
                   "latency to visit food patch under high threat",
                   "lateral forebrain bundle size",
                   "learning score",
                   "leg 1 femur length",
                   "leg 1 tibia length",
                   "leg length",
                   "length",
                   "length at sexual maturity",
                   "length of the iridescent portion of feather",
                   "life span",
                   "litter size (150 days postnatal)",
                   "litter size at 1 year",
                   "litter size in natural conditions",
                   "litter size in supplemented conditons",
                   "locomotion",
                   "longevity",
                   "lower caudal dermatotrichia",
                   "lower pharyngea jaw keel",
                   "lung weight to body weight ratio",
                   "mass (log change in g/day)",
                   "mass at 20 days ",
                   "mass at hatching ",
                   "mass at maturity",
                   "mass at metamorphosis",
                   "mass growth rate",
                   "maximum sprint speed",
                   "maximum sprinting distance between rests",
                   "mean gill filament length",
                   "mean increase in mass after 100 days",
                   "mean increase in snout-vent length after 100 days",
                   "mean liver glycogen",
                   "mean number of song strophes per choice test",
                   "mean selected temperature",
                   "mean stroke force during swim test",
                   "mean time to metamorphosis",
                   "metabolised energy coefficient",
                   "midbrain glutathione peroxidase",
                   "midbrain non enzymatic antioxidant capacity",
                   "midbrain protein carbonyls",
                   "midbrain superoxide dismutase",
                   "mitogenic response to concanavalin A (cell mediated immunity)",
                   "mitogenic response to polyinosinic-polycytidylic acid (cell mediated immunity)",
                   "nerve cell count",
                   "number of copulations",
                   "number of different copulation partners",
                   "number of eggs in first egg sac",
                   "number of eggs produced",
                   "number of fledglings produced",
                   "number of reproductive attempts",
                   "number of sperm at day 0 of experiment",
                   "number of sperm at day 1 of experiment",
                   "number of stimulations to achieve eliptogenesis",
                   "number of stops over 1m",
                   "number of vertebrae",
                   "oestradiol concentration in diencephalon(pg/g)",
                   "oestradiol concentration in hippocampus (pg/g)",
                   "oestradiol concentration in medial prefrontal cortex (pg/g)",
                   "offspring size",
                   "offspring size  in poor environment",
                   "optical density of 5-ht-immunoreactive fibres and varicosities in anterior hypothalymus",
                   "optical density of 5-ht-immunoreactive fibres and varicosities in anterior hypothalymus (2h after resident intruder exposure)",
                   "optical density of 5-ht-immunoreactive fibres and varicosities in basolateral amygdala",
                   "optical density of 5-ht-immunoreactive fibres and varicosities in basolateral amygdala (2h after resident intruder exposure)",
                   "optical density of 5-ht-immunoreactive fibres and varicosities in dosromedial hypothalamic nucleus",
                   "optical density of 5-ht-immunoreactive fibres and varicosities in dosromedial hypothalamic nucleus (2h after resident intruder exposure)",
                   "optical density of 5-ht-immunoreactive fibres and varicosities in lateral hypothalmic area",
                   "optical density of 5-ht-immunoreactive fibres and varicosities in lateral hypothalmic area (2h after resident intruder exposure)",
                   "optical density of 5-ht-immunoreactive fibres and varicosities in supraoptic nucleus",
                   "optical density of 5-ht-immunoreactive fibres and varicosities in supraoptic nucleus (2h after resident intruder exposure)",
                   "ovarian weight",
                   "ovarian weight ",
                   "ovariole number per ovary",
                   "pectoral lepidotrichia",
                   "pelvis width",
                   "percent time spent thrashing during predator cue",
                   "percent time stock female spent with control and treatment male",
                   "percentage exploration behaviour in aggresion test",
                   "percentage of black colouration",
                   "percentage of time spent in open arms of maze",
                   "percentage self grooming during aggression test",
                   "percentage time social behaviour in aggresion test",
                   "period of first incubation stage",
                   "period of first incubation stage in poor environment",
                   "period of second incubation stage",
                   "period of second incubation stage in poor environment",
                   "pituitary weight to body weight ratio",
                   "plasma lysozyme concentrations",
                   "plasma proteins in natural conditions",
                   "plasma proteins in supplemented conditons",
                   "plasma testosterone levels",
                   "playing behaviour",
                   "preoptic area size",
                   "primary feather length",
                   "pronotum length",
                   "proportion of sperm replenished",
                   "pup weight in natural conditions",
                   "pup weight in supplemented conditons",
                   "pupal mass",
                   "pupal weight",
                   "rate of  raising young (no of successful broods/reproductive lifespan)",
                   "rate of  raising young (no of successful broods/reproductivelifespan) in poor environment",
                   "red blood cell protein carbonyls",
                   "red blood cell superoxide dismutase",
                   "reproductive lifespan",
                   "reproductive lifespan in poor environment",
                   "resting metabolic rate",
                   "scaled mass index",
                   "scapular feather length",
                   "seconds spent begging",
                   "self directed behaviour (e.g. scratching, grooming)",
                   "snout-vent length",
                   "snout-vent length ",
                   "snout-vent length at 20 days ",
                   "snout-vent length at hatching ",
                   "snout-vent length at metamorphosis",
                   "snout-vent length growth (log change in g/day)",
                   "snout-vent length growth rate",
                   "song onset",
                   "song rate",
                   "song repertoire size",
                   "song repoirtoire",
                   "songbout length",
                   "spawning rate (no of spawnings/reproductive lifespan)",
                   "spawning rate (no of spawnings/reproductive lifespan) in poor environment",
                   "speed over 0.25m at 18c ",
                   "speed over 0.25m at 23.5c ",
                   "speed over 0.25m at 29c ",
                   "speed over 1m ",
                   "speed over 25cm ",
                   "sperm velocity ",
                   "spermatophore size",
                   "sprint speed",
                   "sprint speed ",
                   "spur length",
                   "sternohyoideus cross section area",
                   "sternohyoideus length",
                   "stops during running",
                   "straight carapace length  ",
                   "stress recovery (latency to feed after handling)",
                   "sum of asymmetrics",
                   "survival time",
                   "tail length",
                   "tail length ",
                   "tail length at 20 days ",
                   "tail length at hatching",
                   "tail width",
                   "tarsus length",
                   "testosterone concentration of eggs laid",
                   "testosterone in natural conditions",
                   "testosterone in supplemented conditons",
                   "thermoregularity precision",
                   "thorax length",
                   "thorax ratio",
                   "thorax ratio (thorax dry mass/total dry mass)",
                   "throat feather length",
                   "tibiotarsus length",
                   "time moving in novel environment",
                   "time spent climbing in forced swim test",
                   "time spent clinging to surrogate",
                   "time spent immobile in forced swim test",
                   "time spent in open areas of maze test",
                   "time spent solving puzzle on unsuccessful attempt",
                   "time spent swimming in forced swim test",
                   "time swimming during 8 minute trial",
                   "time swimming during 8min trial",
                   "time taken to escape water",
                   "time to adulthood",
                   "time to finish maze",
                   "total biomass produced",
                   "total biomass produced in poor environment",
                   "total circulating androgens",
                   "total eggs produced",
                   "total filament number",
                   "total gill filament length",
                   "total granule cell number in brain",
                   "total immunoglobulin level in natural conditions",
                   "total immunoglobulin level in supplemented conditons",
                   "total length",
                   "total number of pyramidal nuerons in left hippocampus",
                   "total number of young",
                   "total number of young in poor environment",
                   "tumour necrosis factor (tnf)-a (cytokine) concentration in blood", #annoying character in the name
                   "ulna length",
                   "ultraviolet chroma",
                   "upper caudal dermatotrichia",
                   "upper pharyngea jaw depth",
                   "uterine weight",
                   "uterine weight ",
                   "vasopressin-immunoreactive staining in nucleus circularis",
                   "vasopressin-immunoreactive staining in nucleus circularis (2h after resident intruder exposure)",
                   "vasopressin-immunoreactive staining in posterior part of the paraventricular hypothalymus",
                   "vasopressin-immunoreactive staining in posterior part of the paraventricular hypothalymus (2h after resident intruder exposure)",
                   "vasopressin-immunoreactive staining in the lateral hypothalamic area",
                   "vasopressin-immunoreactive staining in the lateral hypothalamic area (2h after resident intruder exposure)",
                   "vasopressin messenger ribonucleic-acid expression in bed nucleus of the stria terminalis ",
                   "ventromedial nucleus of the hypothalymus size",
                   "visual exploration",
                   "volume of central body",
                   "volume of granular cell layer in brain",
                   "volume of molecular cell layer in brain",
                   "volume of protocerebral neurophil",
                   "volume of subdivisions in left hippocampus",
                   "wattle area",
                   "wattle chroma",
                   "wattle chroma (uv)",
                   "wattle hue",
                   "wattle hue (visible)",
                   "weight",
                   "weight at maturity",
                   "weight gain at 20 days gestation (1 year postnatal)",
                   "weight gain at 20 days gestation (150 days postnatal)",
                   "weight of protocerebrum",
                   "weight of protocerebrum ",
                   "wet body mass",
                   "whole brain metabolic capacity (umol/min/g tissue wet weight)", #annoying character in the name
                   "wing chord length",
                   "wing colour patch area",
                   "wing length",
                   "wing length ")


# the following vector contains a list of levels of 
# "specific.trait" that need to be checked before deciding
# whether they are proportional or non-proportional data
unknown.data <- c("after discharge threshold prior to kindling",
                  "cell mediated immune response",
                  "cell mediated immunity",
                  "cell mediated immunity (difference between pre and post immunization)",
                  "cling to surrogate",
                  "desiccation rate (change in mg/h)c",
                  "dominance interactions won",
                  "microbiota similarity (%) from fecal samples",
                  "number of assymmetric characters",
                  "number of individuals defeated in domince display",
                  "offspring matured",
                  "percentage of pop breeding in natural conditions",
                  "percentage of pop breeding in supplemented conditons",
                  "predator response (approached dummy predator)",
                  "pupal encapsulation rate",
                  "tadpole survival",
                  "thermal preference",
                  "thermal preference ",
                  "visiting frequency")


# # testing whether all levels where checked
# length(non.prop.data)+length(unknown.data) == length(levels(factor(db.full.red[db.full.red$origin!="STATS",c("specific.trait")])))
# setdiff(c(non.prop.data,unknown.data),levels(factor(db.full.red[db.full.red$origin!="STATS",c("specific.trait")])))
# setdiff(levels(factor(db.full.red[db.full.red$origin!="STATS",c("specific.trait")])),c(non.prop.data,unknown.data))

# # extra double-check
# # this list was put together by NPM in a very preliminary check in Jan 2019
# potential.prop.data.NPM.Jan2019 <- c("pupal encapsulation rate",
#                                      "predator response (approached dummy predator)",
#                                      #"successful copulations",
#                                      #"number of copulatory mounts",
#                                      "rate of  raising young (no of successful broods/reproductivelifespan) in poor environment",
#                                      "spawning rate (no of spawnings/reproductive lifespan)",
#                                      "spawning rate (no of spawnings/reproductive lifespan) in poor environment",
#                                      "rate of  raising young (no of successful broods/reproductive lifespan)",
#                                      "percentage of pop breeding in natural conditions",
#                                      "percentage of pop breeding in supplemented conditons",
#                                      "thorax ratio",
#                                      "percent time spent thrashing during predator cue",
#                                      "average number of deformed fish per temperature treatment",
#                                      "thorax ratio (thorax dry mass/total dry mass)",
#                                      #"learning trials solved",
#                                      "cell mediated immunity (difference between pre and post immunization)",
#                                      "percentage of black colouration",
#                                      "percent time stock female spent with control and treatment male",
#                                      "percentage exploration behaviour in aggresion test",
#                                      "percentage self grooming during aggression test",
#                                      "percentage time social behaviour in aggresion test",
#                                      #"percentage of entries into open arms",
#                                      "percentage of time spent in open arms of maze",
#                                      #"rate of learning problem solving task",
#                                      #"daily survival rate",
#                                      #"survival rate",
#                                      #"wing aspect ratio",
#                                      #"adult encapsulation rate",
#                                      #"peak flight metabolic rate",
#                                      "tadpole survival")#
#                                      #"negative reaction to novel stimulus")
# 
# 
# setdiff(potential.prop.data.NPM.Jan2019,c(non.prop.data,unknown.data))
# setdiff(unknown.data,potential.prop.data.NPM.Jan2019)
# setdiff(potential.prop.data.NPM.Jan2019,unknown.data)
# 
# # 14th of May 2019
# # NPM and AST have gone through the differences and agreed that the new
# # list is fine, and that those that were included in the JAN2019 but not
# # in the new list do not need to be added to the new list.


##############################################################
# subsetting those levels so that they can be further explored
# in the original papers. No "STATS" subsetting will be applied
# for simplicity.

# exporting the list of references that need to be double-checked 

double.check <- db.full.red[db.full.red$specific.trait%in%unknown.data,]

# write.xlsx(double.check[order(double.check$studyID),],
#            "data_re-extraction/prop_data_exploration/list_of_references_to_check_for_group-level_prop_data.xlsx",
#            sheetName="Sheet1",col.names=TRUE, row.names=F,
#            append=FALSE, showNA=TRUE, password=NULL)


# creating a variable to identify group-level proportional data 
# (and alike)
db.full.red$prop.data <- ifelse(db.full.red$specific.trait %in% non.prop.data,
                                0,
                                NA)

# the following vector contains those levels from non.prop.data
# that were checked in the original papers and identified as 
# non-group-level proportional data
non.prop.data.confirmed <- c("offspring matured", # studyID=12:number of chicks fledged
                             "cling to surrogate", # studyID=36:frequency and duration of contact with an artificial mother-surrogate
                             "tadpole survival", # studyID=43:Survival was expressed as the proportion of larvae per tub that completed development (N=38 tubs/treatment)
                             "pupal encapsulation rate", # studyID=47,111:The degree of encapsulation was analysed as grey values of reflecting light from the implants. N's wrong, F df = 1,112 (same method in both studies)
                             "visiting frequency", # studyID=59:individual level behaviour
                             "cell mediated immune response", # studyID=62,70:thickness of web wing/individual
                             "cell mediated immunity (difference between pre and post immunization)", # studyID=71:thickness of web wing/individual
                             "cell mediated immunity", # studyID=71:thickness of web wing/individual
                             "dominance interactions won", #studyID=68:number measured for each individual level
                             "number of individuals defeated in domince display", #studyID=68:number measured for each individual level
                             "predator response (approached dummy predator)", #studyID=98:number of times that the subject moved within 14cm of the predator
                             "number of assymmetric characters", #studyID=101:number of asymmetrical characters = how many of the 12 characters had left-right differences not equal to zero
                             "desiccation rate (change in mg/h)c", # studyID=103:desiccation rate was assessed as short-term evaporative water loss (change in mass/h) for each individual
                             "thermal preference ", # studyID=103:The selected body temperature of each individual was determined as the mean of the six temperature measurements during the 1-h observation period
                             "thermal preference", # studyID=103:The selected body temperature of each individual was determined as the mean of the six temperature measurements during the 1-h observation period
                             "after discharge threshold prior to kindling" # studyID=118:electric intensity in microA at which an afterdischarge for at least 6 s was observed on the EEG recording
)


# the following vector contains those levels from non.prop.data
# that were checked in the original papers and identified as 
# group-level proportional data
prop.data.confirmed <- c("microbiota similarity (%) from fecal samples", # studyID=20:based on Dice's similarity coefficient, links read: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1415224/ and wikipedia 
                         "percentage of pop breeding in natural conditions", # studyID=92:does not seem to be a mean % across enclosures but a % of individuals from all the individuals from all enclosures
                         "percentage of pop breeding in supplemented conditons" # studyID=92:does not seem to be a mean % across enclosures but a % of individuals from all the individuals from all enclosures
)


##############################################################
# updating the variable "prop.data" based on the data obtained
# via checking the original papers
db.full.red$prop.data <- ifelse(db.full.red$specific.trait %in% non.prop.data.confirmed,
                                0,
                                db.full.red$prop.data)

db.full.red$prop.data <- ifelse(db.full.red$specific.trait %in% prop.data.confirmed,
                                1,
                                db.full.red$prop.data)


# further double-checking the dataset before proceeding with
# the next steps

# focusing on the subset containing raw data, i.e. excluding
# the data for which origin=="STATS"
summary(db.full.red[db.full.red$origin!="STATS",])

db.full.red[db.full.red$origin!="STATS" & is.na(db.full.red$prop.data),]
#db.full.red[db.full.red$origin!="STATS" & is.na(db.full.red$prop.data),"specific.trait"]

# there are prop.data missing for 5 rows because of problems matching the 
# specific.trait names with tricky character (alpha,gamma...). Solution:
# assign their prop.data manually using the following code:

db.full.red[db.full.red$origin!="STATS" & is.na(db.full.red$prop.data),"prop.data"] <- 0
# checking everything worked! It did.
#db.full.red[db.full.red$origin!="STATS" & is.na(db.full.red$prop.data),]


db.full.red[db.full.red$origin!="STATS" & is.na(db.full.red$ratioscale),] 

# these data are ratioscale = NA because they are indeed test stastics data
# Furthemore, they have been assigned a prop.data value of 0. this is because
# this data correspond to the PILOT re-extraction subset. To standardize
# the dataset, this values will be manually set up to NA using the
# following code:

db.full.red[db.full.red$origin!="STATS" & is.na(db.full.red$ratioscale),"prop.data"]<- NA
# checking everything worked! It did.
db.full.red[db.full.red$origin!="STATS" & is.na(db.full.red$ratioscale),] 


summary(db.full.red[db.full.red$origin=="STATS",])
# There are some STATs data with prop.data assigned as 0. To standardize
# the dataset, this values will be manually set up to NA using the
# following code:

db.full.red[db.full.red$origin=="STATS","prop.data"] <- NA
# checking everything worked! It did.
summary(db.full.red[db.full.red$origin=="STATS",])


summary(db.full.red) #all seems good!



##############################################################
# Identifying sign inversions in original dataset
##############################################################

# To be able to meaningfully combine data from multiple traits
# there is generally the need of inverting some of the effect
# size signs. This is because traits differ in the sign of the 
# expected relationship between the trait and the response variable. 
# For example, stress hormones are expected to increase with an
# increase in environmental stress, whereas  body size is
# generally expected to decrease with an increase in environmental
# stress.

# For our re-analysis, we would like to use the same sign decisions 
# used in the original meta-analysis. Since we do not have access to 
# a variable noting which signs where inversed and which were not,
# we will calculate those inversions based on the data and  
# corresponding effect sizes provided by H.Eyck

# importing dataset shared by H.Eyck
eyck.dat <- read.xlsx("data/EyckDev.stress_Data_FULL_TABLE.xlsx",
                      colNames=T,sheet = 1)


# transforming some variables into factors
var.factor.2<- c("TAXA","scientific.name","developmental.stressor","method","exposure","sex",
                 "trait.class","age","specific.trait")

eyck.dat[,var.factor.2] <- lapply(eyck.dat[,var.factor.2] , factor)


# Reducing the number of variables to those selected for our dataset 
# (i.e. those that we are interested in and those moderators without 
# NA's; more info above) plus those additional ones that we need 
# for this section (i.e. effect size)

subset.variables.2 <- c("studyID","speciesID","TAXA","common.name","scientific.name","citation",
                        "location","lat","long","av..temp","htmnth","cldmnth",
                        "seasonality","developmental.stressor",
                        "sex","age","trait.class","specific.trait",
                        "mean.control","SD","N.control",
                        "mean.treat","SD.1","N.treat",
                        "t.value","F.value",
                        "Cohens.D","sv")


eyck.dat.red <- eyck.dat[,subset.variables.2]

# performing the same 'moderator cleaning' performed for our 
# cleaned dataset (i.e. db.full.red; see above)
eyck.dat.red$trait.class <- revalue(eyck.dat.red$trait.class, c("behavioural "="behavioural"))
eyck.dat.red$sex <- revalue(eyck.dat.red$sex, c("Female"="female", "females"="female"))
eyck.dat.red$sex <- revalue(eyck.dat.red$sex, c("Male"="male", "male "="male","males"="male"))

eyck.dat.red[,c("trait.class","sex")] <- lapply(eyck.dat.red[,c("trait.class","sex")] , factor)


# creating a unique identifier that can be used to link data from
# this dataset to our cleaned dataset (i.e. db.full.red)
names(eyck.dat.red)
names(db.full.red)

eyck.dat.red$unique.identifier <- paste(eyck.dat.red$studyID,
                                        eyck.dat.red$speciesID,
                                        eyck.dat.red$specific.trait,
                                        eyck.dat.red$sex,
                                        eyck.dat.red$age,
                                        sep="_")

db.full.red$unique.identifier <- paste(db.full.red$studyID,
                                       db.full.red$speciesID,
                                       db.full.red$specific.trait,
                                       db.full.red$sex,
                                       db.full.red$age,
                                       sep="_")


eyck.dat.red$data.sign.HE <- factor(ifelse(eyck.dat.red$mean.treat-eyck.dat.red$mean.control<0,
                                           "negative",
                                           "positive.or.zero"))


eyck.dat.red$es.sign.HE <- factor(ifelse(eyck.dat.red$Cohens.D<0,
                                         "negative",
                                         "positive.or.zero"))


eyck.dat.red$sign.inversion.HE <- factor(ifelse(is.na(eyck.dat.red$data.sign.HE),
                                                NA,
                                                ifelse(eyck.dat.red$data.sign.HE != eyck.dat.red$es.sign.HE,
                                                       -1,
                                                       1)))



# some double-checking
table(eyck.dat.red$sign.inversion.HE) #there are 44 signs that seemed to have been inverted in the original dataset
levels(factor(eyck.dat.red[eyck.dat.red$sign.inversion.HE==-1 & !(is.na(eyck.dat.red$sign.inversion.HE)),"specific.trait"]))
# Quick double-checking to understand what effect size types were inverted: 
# [1] "aggression towards conspecific": sounds ok                     
# [2] "bites to mirror": sounds ok                                     
# [3] "caudal peduncle width": found to be wrong, and it is corrected below                              
# [4] "cling to surrogate": sounds ok                                 
# [5] "development time": there are two effect sizes called "final instar development time" that have not been inverted, despite that both are in insects, nutritional treatment, adult...                                  
# [6] "exploratory behaviour": suspicious: don't understand why activity in other studies does not seem to be inverted, typo?                         
# [7] "incubation length": sounds ok                                  
# [8] "incubation period": sounds ok                                
# [9] "latency (s) to leave starting box and enter maze": sounds ok   
# [10] "latency to escape": sounds ok                                  
# [11] "latency to lay eggs after pairing": sounds ok                  
# [12] "latency to leave starting box and enter maze": sounds ok       
# [13] "mean time to metamorphosis": sounds ok                         
# [14] "predator response (approached dummy predator)": suspicious       
# [15] "resting metabolic rate": sounds ok                             
# [16] "seconds spent begging": sounds ok                              
# [17] "self directed behaviour (e.g. scratching, grooming)": sounds ok
# [18] "survival time": other studies in this dataset studied longevity, and sign was not reverted, thus, this has to be a typo?                                      
# [19] "time spent clinging to surrogate": sounds ok                  
# [20] "time spent immobile in forced swim test": sounds ok            
# [21] "time spent solving puzzle on unsuccessful attempt": sounds ok  
# [22] "time taken to escape water": sounds ok                         
# [23] "time to finish maze": sounds ok                                
# [24] "volume of protocerebral neurophil": found to be wrong, and it is corrected below

# traits.to.be.inverted <- as.character(levels(factor(eyck.dat.red[eyck.dat.red$sign.inversion==1 & !(is.na(eyck.dat.red$sign.inversion)),"specific.trait"])))
# 
# nrow(db.full.red[db.full.red$origin!="STATS" & db.full.red$specific.trait %in% traits.to.be.inverted,])
# xxx <- db.full.red[db.full.red$origin!="STATS" & db.full.red$specific.trait %in% traits.to.be.inverted,]
# xxx[,c("origin","specific.trait","mean.treat","mean.control")]

# vector containing the unique.identifiers for which effect size sign needs to be inverted
eyck.dat.super.red <- unique(eyck.dat.red[eyck.dat.red$sign.inversion.HE==-1 &
                                            !(is.na(eyck.dat.red$sign.inversion.HE)),
                                          c("unique.identifier","sign.inversion.HE")])

db.full.red.sign <- merge(db.full.red,eyck.dat.super.red,by="unique.identifier",all.x=T)
summary(db.full.red.sign)
# there is one additional effect size. Check why:
# quick.double.check<- merge(eyck.dat.red,eyck.dat.super.red,by="unique.identifier",all.x=T)
# quick.double.check[quick.double.check$sign.inversion.y==1  & !(is.na(quick.double.check$sign.inversion.y)),]
# the reason is that there is a typo in the effect size estimation. 
# The following code fixes that issue:
db.full.red.sign[db.full.red.sign$specific.trait=="caudal peduncle width","sign.inversion.HE"]<-NA
db.full.red.sign[!(is.na(db.full.red.sign$mean.treat)) &
                   is.na(db.full.red.sign$sign.inversion.HE),"sign.inversion.HE"]<-"1"


# # additional, the following sign seems like a typo because,
# # according to other variables in the dataset, it should not 
# # have been reverted, and thus, it is re-reverted using 
# # the following code
# db.full.red.sign[db.full.red.sign$specific.trait=="survival time","sign.inversion.HE"]<-"1"


# checking whether there are duplicated studies with same citation but
# different study ID (such as that found for 109 vs. 119; see below)
unique_studyID_citation.db<-unique(db.full.red.sign[,c("studyID","citation")])
unique_studyID_citation.db.counts <- plyr::count(unique_studyID_citation.db,"citation")
unique_studyID_citation.db.counts[unique_studyID_citation.db.counts$freq>1,]


# by trying to find out why "volume of protocerebral neurophil" was reverted
# we found out that studyID==109 was a duplicate version of studyID==119
# Thus, we here delete studyID==109
db.full.red.sign <- db.full.red.sign[db.full.red.sign$studyID!=109,]
# we also realized that during our duble-check, one of us missed the 
# typo on mean.treat for studyID==109


##############################################################
# Assigning sign inversions by ourselves
##############################################################

# the previous section has shed doubts on whether sign inversions 
# where systematically implemented (and how). Thus, we have decided 
# to go through all the "specific.trait" measured for each study, 
# and decide ourselves which signs should be or not inverted.
db.full.red.sign$specific.trait <- factor(db.full.red.sign$specific.trait)
levels(factor(db.full.red.sign[db.full.red.sign$origin!="STATS",c("specific.trait")]))


##############################################################
# First round
##############################################################

# In this round, we went through all the levels of the variable
# "specific.trait" and, based on the description provided, 
# categorize those levels as: inversion not needed, inversion
# needed, inversion unclear and variable measured unclear.


# the following vector contains a list of levels of 
# "specific.trait" for which a negative correlation
# with stress is expected. Thus, the effect size
# sign for this specific.trait WOULD NOT need to be
# inverted
no.inversion.needed <- c("abdominal length",
                         "abdominal width",
                         "age at death",
                         "arm length",
                         "arm length ",
                         "average clutch size",
                         "average hatch success per female",
                         "average total eggs per female",
                         "beak colour score",
                         "beak score (colour, hue, brightness, saturation)",
                         "bill coloration length",
                         "bill hue",
                         "body depth",
                         "body mass",
                         "body mass ",
                         "body mass at hatching ",
                         "body mass at metamorphosis",
                         "body weight",
                         "body weight ",
                         "body weight at 20 days of gestation (1 year postnatal)",
                         "body weight at 20 days of gestation (150 days postnatal)",
                         "body weight at first oestrus",
                         "body weight at vaginal opening ",
                         "brood fledgling mass",
                         "caudal fin width",
                         "caudal peduncle width",
                         "cell mediated immune response",
                         "cell mediated immunity",
                         "cell mediated immunity (difference between pre and post immunization)",
                         "cephalothorax area",
                         "cephalothorax length",
                         "cephalothorax width",
                         "clutch mass",
                         "clutch size",
                         "clutch size in poor environment",
                         "cognitive performance",
                         "comb and wattle size",
                         "culmen length",
                         "development rate",
                         "dominance interactions won",
                         "dry body mass",
                         "dry weight",
                         "ear tuft",
                         "egg mass",
                         "egg size",
                         "epaxial muscle depth",
                         "fat content",
                         "fat content [%]",
                         "feather mass",
                         "feather tip length",
                         "femur length",
                         "femur width",
                         "forebrain volume",
                         "forebrain weight ",
                         "forehead colour patch area",
                         "forewing length",
                         "gonopodium length",
                         "growth rate",
                         "haematocrit",
                         "haematocrit (ratio of red blood cells to whole blood volume) in natural conditions",
                         "haematocrit (ratio of red blood cells to whole blood volume) in supplemented conditons",
                         "head and bill length",
                         "head length",
                         "head length ",
                         "head width",
                         "head width ",
                         "hematocrit (% packed cell volume)",
                         "hemibranch area",
                         "hemibranch perimeter",
                         "hours until death by starvation following metamorphosis",
                         "humerus length",
                         "latency to immobility in forced swim test", #read about it on the internet and immobility is meant to show hoplesness
                         "learning score",
                         "leg 1 femur length",
                         "leg 1 tibia length",
                         "leg length",
                         "length",
                         "length at sexual maturity",
                         "length of the iridescent portion of feather",
                         "life span",
                         "litter size (150 days postnatal)",
                         "litter size at 1 year",
                         "litter size in natural conditions",
                         "litter size in supplemented conditons",
                         "longevity",
                         "mass at 20 days ",
                         "mass at hatching ",
                         "mass at maturity",
                         "mass at metamorphosis",
                         "mass growth rate",
                         "maximum sprint speed",
                         "maximum sprinting distance between rests",
                         "mean gill filament length",
                         "mean increase in mass after 100 days",
                         "mean increase in snout-vent length after 100 days",
                         "mean number of song strophes per choice test",
                         "mean stroke force during swim test",
                         "number of copulations",
                         "number of different copulation partners",
                         "number of eggs in first egg sac",
                         "number of eggs produced",
                         "number of fledglings produced",
                         "number of individuals defeated in domince display",
                         "number of reproductive attempts",
                         "number of sperm at day 0 of experiment",
                         "number of sperm at day 1 of experiment",
                         "offspring matured",
                         "offspring size",
                         "offspring size  in poor environment",
                         "ovarian weight",
                         "ovarian weight ",
                         "ovariole number per ovary",
                         "pelvis width",
                         "percentage of black colouration",
                         "percentage of pop breeding in natural conditions",
                         "percentage of pop breeding in supplemented conditons",
                         "playing behaviour",
                         "primary feather length",
                         "pronotum length",
                         "proportion of sperm replenished",
                         "pup weight in natural conditions",
                         "pup weight in supplemented conditons",
                         "pupal mass",
                         "pupal weight",
                         "rate of  raising young (no of successful broods/reproductive lifespan)",
                         "rate of  raising young (no of successful broods/reproductivelifespan) in poor environment",
                         "reproductive lifespan",
                         "reproductive lifespan in poor environment",                         
                         "scapular feather length",
                         "snout-vent length",
                         "snout-vent length ",
                         "snout-vent length at 20 days ",
                         "snout-vent length at hatching ",
                         "snout-vent length at metamorphosis",
                         "snout-vent length growth (log change in g/day)",
                         "snout-vent length growth rate",
                         "song rate",
                         "song repertoire size",
                         "song repoirtoire",
                         "songbout length",
                         "spawning rate (no of spawnings/reproductive lifespan)",
                         "spawning rate (no of spawnings/reproductive lifespan) in poor environment",
                         "speed over 0.25m at 18c ",
                         "speed over 0.25m at 23.5c ",
                         "speed over 0.25m at 29c ",
                         "speed over 1m ",
                         "speed over 25cm ",
                         "sperm velocity ",
                         "spermatophore size",
                         "sprint speed",
                         "sprint speed ",
                         "spur length",
                         "straight carapace length  ",
                         "survival time",
                         "tadpole survival",
                         "tail length",
                         "tail length ",
                         "tail length at 20 days ",
                         "tail length at hatching",
                         "tail width",
                         "tarsus length",
                         "thorax length",
                         "throat feather length",
                         "tibiotarsus length",
                         "total biomass produced",
                         "total biomass produced in poor environment",
                         "total eggs produced",
                         "total filament number",
                         "total gill filament length",
                         "total immunoglobulin level in natural conditions",
                         "total immunoglobulin level in supplemented conditons",
                         "total length",
                         "total number of young",
                         "total number of young in poor environment",
                         "ulna length",
                         "uterine weight",
                         "uterine weight ",
                         "wattle area",
                         "wattle chroma",
                         "wattle chroma (uv)",
                         "wattle hue",
                         "wattle hue (visible)",
                         "weight",
                         "weight at maturity",
                         "weight gain at 20 days gestation (1 year postnatal)",
                         "weight gain at 20 days gestation (150 days postnatal)",
                         "weight of protocerebrum ",
                         "wet body mass",
                         "wing chord length",
                         "wing colour patch area",
                         "wing length",
                         "wing length ")


# the following vector contains a list of levels of 
# "specific.trait" for which a positive correlation
# with stress is expected. Thus, the effect size
# sign for this specific.trait WOULD need to be
# inverted
inversion.needed <- c("adrenocorticotropic hormone response to anesthesia",
                      "adrenocorticotropic hormone response to stress challenge",
                      "age at first mating ",
                      "age at first oviposition",
                      "age at metamorphosis",
                      "age at vaginal opening (puberty) ",
                      "age of first oestrus ",
                      "age post maturation at mating",
                      "average number of deformed fish per temperature treatment",
                      "baseline cort",
                      "body mass change during first immune challenge",
                      "cling to surrogate",
                      "cort 1hr after swimming water maze",
                      "cort 30 mins after handling stress",
                      "corticosterone concentration",
                      "cortisol response to anesthesia",
                      "cortisol response to stress challenge",
                      "days to learn associative task",
                      "days until pupation",
                      "deformed vertebrae",
                      "desiccation rate (change in mg/h)c",
                      "development duration",
                      "development time",
                      "fecal pellets in response to novel stress",
                      "final instar development time ",
                      "hours until metamorphosis after hatching",
                      "incubation length",
                      "incubation period",
                      "incubation period ",
                      "latency to find food reward after training ",
                      "latency to lay",
                      "latency to lay eggs after pairing",
                      "mean time to metamorphosis",
                      "number of assymmetric characters",
                      "percentage of time spent in open arms of maze",
                      "song onset",
                      "stops during running",
                      "stress recovery (latency to feed after handling)",
                      "sum of asymmetrics",
                      "time spent clinging to surrogate",
                      "time spent immobile in forced swim test",
                      "time spent solving puzzle on unsuccessful attempt",
                      "time to adulthood",
                      "time to finish maze",
                      "tumour necrosis factor (tnf)-a (cytokine) concentration in blood" #annoying character in the name
                      
)


# the following vector contains a list of levels of 
# "specific.trait" for which is difficult to predict
# whether the association with stress should be 
# positive or negative. Thus, it is UNCLEAR whether
# an effect size sign inversion for this specific.trait 
# would be needed
inversion.unclear <- c("activity level (lines crossed)",
                       "aggression towards conspecific",
                       "aterial pressure",
                       "attempted mounts",
                       "average flight time under flight stress",
                       #"average sprint speed", # sorted when doing the ones with unclear meaning
                       "basal metabolic rate",
                       "bites to mirror",
                       "body shape at hatching (mass/snout-vent length)",
                       "body temperature",
                       "climbing",
                       "distance travelled during test session",
                       "exploration",
                       "exploratory behaviour",
                       "feeding latency",
                       "female activity, calling and tail quivering in response to male",
                       "heart rate",
                       "hopping distance",
                       "initial activity",
                       "latency to escape",
                       #"latency to exit begin exploring", # sorted when doing the ones with unclear meaning
                       "latency to feed from enclosed novel environment",
                       #"latency to leave starting box and enter maze", # sorted when doing the ones with unclear meaning
                       "latency to visit food patch under high threat",
                       "locomotion",
                       "percent time spent thrashing during predator cue",
                       "percentage exploration behaviour in aggresion test",
                       "predator response (approached dummy predator)",
                       "resting metabolic rate",
                       "seconds spent begging",
                       "self directed behaviour (e.g. scratching, grooming)",
                       "testosterone concentration of eggs laid",
                       "testosterone in natural conditions",
                       "testosterone in supplemented conditons",
                       "thorax ratio",
                       "thorax ratio (thorax dry mass/total dry mass)",
                       "time moving in novel environment",
                       #"time spent climbing in forced swim test", # sorted when doing the ones with unclear meaning
                       #"time spent in open areas of maze test", # sorted when doing the ones with unclear meaning
                       #"time spent swimming in forced swim test", # sorted when doing the ones with unclear meaning
                       #"time swimming during 8 minute trial", # sorted when doing the ones with unclear meaning
                       #"time swimming during 8min trial", # sorted when doing the ones with unclear meaning
                       #"time taken to escape water", # sorted when doing the ones with unclear meaning
                       "total circulating androgens",
                       "upper pharyngea jaw depth"
)


# the following vector contains a list of levels of 
# "specific.trait" that we do not understand what
# they really are without going back to the original 
# references and reading about what they are exactly 
# and how to interpret their relationship with stress.

variable.meaning.unclear <- c("(dihydroxyphenylacetic acid/dopamine) ratio in the nuclues accumbens after observation of a receptive stimulus female",                   
                              "(dihydroxyphenylacetic acid/dopamine) ratio in the preoptic area after observation of a receptive stimulus female",                       
                              "(dihydroxyphenylacetic acid/dopamine) ratio in the ventral tegumental area after observation of a receptive stimulus female",             
                              "adrenal weight to body weight ratio",
                              "after discharge threshold prior to kindling",
                              "anal lepidotrichia",
                              "brain weight to body weight ratio",
                              "breeding time in natural conditions",
                              "breeding time in supplemented conditons",
                              "cerebellum protein carbonyls",
                              "cerebellum superoxide dismutase",
                              "change in length after 1 year",
                              "circulating oestradiol",
                              "condition factor (100*g cm-3)",
                              "cycle length at 140 days ",
                              "cylcle length at 1 year",
                              "density of pyramidal neurons in left hipocampus",
                              "dihydroxyphenylacetic acid in the nuclues accumbens after observation of a receptive stimulus female",
                              "dihydroxyphenylacetic acid in the preoptic area after observation of a receptive stimulus female",
                              "dopamine in the nuclues accumbens after observation of a receptive stimulus female",
                              "dopamine in the preoptic area after observation of a receptive stimulus female",
                              "dopamine in the ventral tegumental area after observation of a receptive stimulus female",
                              "dorsal lepidotrichia",
                              "dorsal spines",
                              "eye depth",
                              "geniohyoideus length",
                              "granulocyte:lymphocyte ratio",
                              "habenula size",
                              "heart weight to body weight ratio",
                              "hepatosomatic index (100*liver wet weight/body wet weight)",
                              "hippocampus weight to brain weight ratio",
                              "humoral immunity",
                              "humoral immunity (difference between pre and post immunization)",
                              "interferone (ifn)-y (cytokine) concentration in blood", #annoying character in the name
                              "interleukin 10 (cytokine) concentration in blood",
                              "interleukin 4 (cytokine) concentration in blood",
                              "interleukin 6 (cytokine) concentration in blood",
                              "key pecks during song preference tests of familiar songs",
                              "key pecks during song preference tests of unfamiliar songs",
                              "kidney weight to body weight ratio",
                              "latency to first appearence in maze",
                              "latency to first egg sac",
                              "lateral forebrain bundle size",
                              "lower caudal dermatotrichia",
                              "lower pharyngea jaw keel",
                              "lung weight to body weight ratio",
                              "mass (log change in g/day)",
                              "mean liver glycogen",
                              "mean selected temperature",
                              "metabolised energy coefficient",
                              "microbiota similarity (%) from fecal samples",
                              "midbrain glutathione peroxidase",
                              "midbrain non enzymatic antioxidant capacity",
                              "midbrain protein carbonyls",
                              "midbrain superoxide dismutase",
                              "mitogenic response to concanavalin A (cell mediated immunity)",
                              "mitogenic response to polyinosinic-polycytidylic acid (cell mediated immunity)",
                              "nerve cell count",
                              "number of stimulations to achieve eliptogenesis",
                              "number of stops over 1m",
                              "number of vertebrae",
                              "oestradiol concentration in diencephalon(pg/g)",
                              "oestradiol concentration in hippocampus (pg/g)",
                              "oestradiol concentration in medial prefrontal cortex (pg/g)",
                              "optical density of 5-ht-immunoreactive fibres and varicosities in anterior hypothalymus",
                              "optical density of 5-ht-immunoreactive fibres and varicosities in anterior hypothalymus (2h after resident intruder exposure)",
                              "optical density of 5-ht-immunoreactive fibres and varicosities in basolateral amygdala",
                              "optical density of 5-ht-immunoreactive fibres and varicosities in basolateral amygdala (2h after resident intruder exposure)",
                              "optical density of 5-ht-immunoreactive fibres and varicosities in dosromedial hypothalamic nucleus",
                              "optical density of 5-ht-immunoreactive fibres and varicosities in dosromedial hypothalamic nucleus (2h after resident intruder exposure)",
                              "optical density of 5-ht-immunoreactive fibres and varicosities in lateral hypothalmic area",
                              "optical density of 5-ht-immunoreactive fibres and varicosities in lateral hypothalmic area (2h after resident intruder exposure)",
                              "optical density of 5-ht-immunoreactive fibres and varicosities in supraoptic nucleus",
                              "optical density of 5-ht-immunoreactive fibres and varicosities in supraoptic nucleus (2h after resident intruder exposure)",
                              "pectoral lepidotrichia",
                              "percent time stock female spent with control and treatment male",
                              "percentage self grooming during aggression test",
                              "percentage time social behaviour in aggresion test",
                              "period of first incubation stage",
                              "period of first incubation stage in poor environment",
                              "period of second incubation stage",
                              "period of second incubation stage in poor environment",
                              "pituitary weight to body weight ratio",
                              "plasma lysozyme concentrations",
                              "plasma proteins in natural conditions",
                              "plasma proteins in supplemented conditons",
                              "plasma testosterone levels",
                              "preoptic area size",
                              "pupal encapsulation rate",
                              "red blood cell protein carbonyls",
                              "red blood cell superoxide dismutase",
                              "scaled mass index",
                              "sternohyoideus cross section area",
                              "sternohyoideus length",
                              "thermal preference",
                              "thermal preference ",
                              "thermoregularity precision",
                              "total granule cell number in brain",
                              "total number of pyramidal nuerons in left hippocampus",
                              "ultraviolet chroma",
                              "upper caudal dermatotrichia",
                              "vasopressin-immunoreactive staining in nucleus circularis",
                              "vasopressin-immunoreactive staining in nucleus circularis (2h after resident intruder exposure)",
                              "vasopressin-immunoreactive staining in posterior part of the paraventricular hypothalymus",
                              "vasopressin-immunoreactive staining in posterior part of the paraventricular hypothalymus (2h after resident intruder exposure)",
                              "vasopressin-immunoreactive staining in the lateral hypothalamic area",
                              "vasopressin-immunoreactive staining in the lateral hypothalamic area (2h after resident intruder exposure)",
                              "vasopressin messenger ribonucleic-acid expression in bed nucleus of the stria terminalis ",
                              "ventromedial nucleus of the hypothalymus size",
                              "visiting frequency",
                              "visual exploration",                                                                                                                      
                              "volume of central body",
                              "volume of granular cell layer in brain",
                              "volume of molecular cell layer in brain",
                              "volume of protocerebral neurophil",
                              "volume of subdivisions in left hippocampus",
                              "whole brain metabolic capacity (umol/min/g tissue wet weight)" #annoying character in the name
)

# # quick double-check: everything ok!
# length(levels(factor(db.full.red.sign[db.full.red.sign$origin!="STATS",c("specific.trait")])))==
# length(c(no.inversion.needed,inversion.needed,inversion.unclear,variable.meaning.unclear))
# 
# setdiff(c(no.inversion.needed,inversion.needed,inversion.unclear,variable.meaning.unclear),levels(factor(db.full.red.sign[db.full.red.sign$origin!="STATS",c("specific.trait")])))
# setdiff(levels(factor(db.full.red.sign[db.full.red.sign$origin!="STATS",c("specific.trait")])),c(no.inversion.needed,inversion.needed,inversion.unclear,variable.meaning.unclear))

# intersect(no.inversion.needed,inversion.needed)
# intersect(no.inversion.needed,inversion.unclear)
# intersect(no.inversion.needed,variable.meaning.unclear)
# intersect(inversion.needed,inversion.unclear)
# intersect(inversion.needed,variable.meaning.unclear)
# intersect(inversion.unclear,variable.meaning.unclear)

# # percentage that does not need to be double-checked
# nrow(db.full.red.sign[db.full.red.sign$specific.trait %in% c(no.inversion.needed,inversion.needed),])/nrow(db.full.red.sign[db.full.red.sign$origin!="STATS",])

##############################################################
# Second round
##############################################################

# In this round, we went through the list of levels for which
# either meaning or direction was unclear, and went and double-
# checked each of them in their original publication. Then, we made
# a final decision on whether sign inversion is needed or not, and left
# some variables as direction unclear. We tried to reduce the list
# of unclear ones to the minimum possible, and made some decision
# for some variables, mostly behaviour, for which one could almost
# argue both ways.


# creating two datasets to make the process easier
# write.xlsx(db.full.red.sign[db.full.red.sign$specific.trait %in% variable.meaning.unclear,
#                             c(2:32)],
#            "data_re-extraction/sign_exploration/list_of_references_to_check_because_of_unclear_variable.xlsx",
#            sheetName="Sheet1",col.names=TRUE, row.names=F,
#            append=FALSE, showNA=TRUE, password=NULL)

# write.xlsx(db.full.red.sign[db.full.red.sign$specific.trait %in% inversion.unclear,
#                             c(2:32)],
#            "data_re-extraction/sign_exploration/list_of_references_to_check_because_of_unclear_direction.xlsx",
#            sheetName="Sheet1",col.names=TRUE, row.names=F,
#            append=FALSE, showNA=TRUE, password=NULL)


# probably because of the annoying character, the following variable was not 
# included in the excel, and needed to be double-checked later on:
# "whole brain metabolic capacity (umol/min/g tissue wet weight)"
# first, getting rid off the annoying characters once forever
db.full.red.sign$specific.trait <- as.character(db.full.red.sign$specific.trait)
db.full.red.sign[c(856:858),"specific.trait"] <- "whole brain metabolic capacity (umol/min/g tissue wet weight)"
db.full.red.sign[c(270),"specific.trait"] <- "interferone (ifn)-y (cytokine) concentration in blood"
db.full.red.sign[c(277),"specific.trait"] <- "tumour necrosis factor (tnf)-a (cytokine) concentration in blood"
db.full.red.sign$specific.trait <- as.factor(db.full.red.sign$specific.trait)

# double-checking the whole brain variable
db.full.red.sign[db.full.red.sign$specific.trait=="whole brain metabolic capacity (umol/min/g tissue wet weight)",]


# After checking the original publications, the 
# unclear (both) variables are now assigned as following: 

no.inversion.needed.2 <- c(no.inversion.needed,
                           c("activity level (lines crossed)", # studyID=98:as in other studies (e.g. exploration), we will assume that stressed animals move less
                             "after discharge threshold prior to kindling", # studyID=118:according to the authors low values are indicative of a pro-convulsive state
                             "anal lepidotrichia", # studyID=82:the authors do not provide any predictions on this but we are assuming that stressed individuals should have less - lepidotrichia = fin spines and rays (https://en.wikipedia.org/wiki/Fish_fin)
                             "attempted mounts", # studyID=54:measured as a proxy of sexual behaviour
                             "average flight time under flight stress", # studyID=78:expect stressed animals to fly less, unless predictive adaptive response, as they argue
                             "average sprint speed",
                             "basal metabolic rate", # studyID=73:authors seem to argue that stressed animals should have higher metabolic rate, however, to standardize across studies (see variable: resting metabolic rate, e.g.), we are going to asssume the opposite
                             "bites to mirror", # studyID=98:authors "we interpret ''biting'' as a general measure of social motivation, or the intent to interact with a social partner."
                             "body shape at hatching (mass/snout-vent length)", # studyID=103:authors do not provide any prediction but body shape = (mass^0.3)/length, and therfore, for a specific size, this number goes up when individuals are heavier
                             "body temperature", # studyID=112:in accordance to resting metabolic rate, we will expect a reduction in body temperature due to food shortage
                             "brain weight to body weight ratio", # studyID=37:the authors do not set predictions for any of the ratios. However, we will assume that brain and gonad ratios should be negatively correlated with stress. The remaining ratios will be part of the unclear category
                             "cerebellum superoxide dismutase", # studyID=50:"superoxide dismutase; glutathione peroxidase and non-enzymatic antioxidant capacity are established indicators of antioxidant defences preventing oxidation of cell components"
                             "change in length after 1 year", # studyID=90:proxy of growth
                             "circulating oestradiol",# studyID=108:from the reading of the study, it seems that for females it is good to have higher levels of oestradiol. Note: not sure what should be the control treatment ("'hot' temperatures (32 ?C) resulting in mostly male offspring and relatively 'cold' temperatures (26 ?C) resulting in only female offspring"). EXCLUDE STUDY?
                             "climbing", # studyID=36:authors expect a negative correlation with stress. Also, it is what we have decided to follow in other similar behavioural variables.
                             "condition factor (100*g cm-3)", # studyID=7: individuals condition factor
                             "dihydroxyphenylacetic acid in the nuclues accumbens after observation of a receptive stimulus female", # studyID=107:metabolite of dopamine, which authors argue that relates to sexual behaviour
                             "dihydroxyphenylacetic acid in the preoptic area after observation of a receptive stimulus female", # studyID=107:metabolite of dopamine, which authors argue that relates to sexual behaviour
                             "distance travelled during test session", # studyID=118:distance/time in the open arms is expected to correlate negatively with stress. Since total distance is likely to correlate with distance in open arms, we are going to assume that stressed animals cover less distance
                             "dopamine in the nuclues accumbens after observation of a receptive stimulus female", # studyID=107: authors argue that dopamine relates to sexual behaviour
                             "dopamine in the preoptic area after observation of a receptive stimulus female", # studyID=107: authors argue that dopamine relates to sexual behaviour
                             "dopamine in the ventral tegumental area after observation of a receptive stimulus female", # studyID=107: authors argue that dopamine relates to sexual behaviour
                             "dorsal lepidotrichia", # studyID=82:the authors do not provide any predictions on this but we are assuming that stressed individuals should have less - lepidotrichia = fin spines and rays (https://en.wikipedia.org/wiki/Fish_fin)
                             "dorsal spines", # studyID=82:the authors do not provide any predictions on this but we are assuming that stressed individuals should have less - lepidotrichia = fin spines and rays (https://en.wikipedia.org/wiki/Fish_fin)
                             "exploration", # studyID=36:authors expect a negative correlation with stress. Also, it is what we have decided to follow in other similar behavioural variables.
                             "exploratory behaviour", #studyID=19:defined as the number of entries made into the inner 24 squares of the arena - percentage of entries into the centre of the open field [(central entries/total entries) * 100] as an index of anxiety (Pellow and File 1986; Frye et al. 2000). Authors say that entering the center of the arenta often is a sign of reduced anxiety-type behaviour
                             "eye depth", # studyID=81:structural element surrounding gill: authors interpret increase as negative and consequence of the need of more space for the gills. Authors argue so by saying that body shape can have consequences for hydrodynamics. Nonetheless, to standardize across studies, we are considering that all these type of morphological variables are expecte to be negatively correlated with stress
                             "female activity, calling and tail quivering in response to male", # studyID=89:meant to be a proxy of sexual behaviour
                             "fork length", # studyID=39:body length (see Fig.1)
                             "gape length", # studyID=39:mouth-to-gill distance (see Fig.1)
                             "geniohyoideus length", # studyID=81:structural element surrounding gill: authors interpret increase as negative and consequence of the need of more space for the gills. Authors argue so by saying that body shape can have consequences for hydrodynamics. Nonetheless, to standardize across studies, we are considering that all these type of morphological variables are expecte to be negatively correlated with stress
                             "gonad weight to body weight ratio", # studyID=37:the authors do not set predictions for any of the ratios. However, we will assume that brain and gonad ratios should be negatively correlated with stress. The remaining ratios will be part of the unclear category
                             "habenula size",# studyID=94:"unlikely to contain sex steroid-concentrating neurons". "measured to evaluate whether gonadal sex and incubation temperature influences are specific to putative steroid binding areas". It is a difficult call but we are going to assume that stressed animals should develop smaller brains in general, which means that they should develop smaller brain areas too (difficult call though)
                             "hopping distance", # studyID=78:proxy of fitness where stressed animals should hop less
                             "humoral immunity", # studyID=71:immune response measured as thickness of wing after immunization
                             "humoral immunity (difference between pre and post immunization)", # studyID=71:immune response measured as thickness of wing after immunization
                             "initial activity", # studyID=101:authors do not provide a prediction regarding this trait. Following other studies on mice (maze and swimming test), and to standardize across studies, we are going to assume that stressed animals move less
                             "key pecks during song preference tests of familiar songs", # studyID=64:key pecks is a proxy of motivation to hear male songs, so it could be argue that it is a proxy for sexual behaviour
                             "key pecks during song preference tests of unfamiliar songs", # studyID=64:key pecks is a proxy of motivation to hear male songs, so it could be argue that it is a proxy for sexual behaviour
                             "lateral forebrain bundle size", # studyID=94:"unlikely to contain sex steroid-concentrating neurons". "measured to evaluate whether gonadal sex and incubation temperature influences are specific to putative steroid binding areas". It is a difficult call but we are going to assume that stressed animals should develop smaller brains in general, which means that they should develop smaller brain areas too (difficult call though)
                             "locomotion", # studyID=36:authors expect a negative correlation with stress. Also, it is what we have decided to follow in other similar behavioural variables.
                             "lower caudal dermatotrichia", # studyID=82: the authors do not provide any predictions on this but we are assuming that stressed individuals should have less - lepidotrichia = fin spines and rays (https://en.wikipedia.org/wiki/Fish_fin)
                             "lower pharyngea jaw keel", # studyID=81:structural element surrounding gill: authors interpret increase as negative and consequence of the need of more space for the gills. Authors argue so by saying that body shape can have consequences for hydrodynamics. Nonetheless, to standardize across studies, we are considering that all these type of morphological variables are expecte to be negatively correlated with stress
                             "mass (log change in g/day)", # studyID=103:"Growth in mass was calculated as the difference between an individual's natural log transformed mass at the time of release and hatching divided by the number of days between measurements"
                             "metabolised energy coefficient", # studyID=46:authors do not provide any prediction regarding MEC. However, we have assumed that stressed animals should be less energentically efficient.
                             "microbiota similarity (%) from fecal samples", # studyID=20:authors interpret a reduction in microbiota similarity as altered/disrupted microbiota
                             "midbrain glutathione peroxidase", # studyID=50:"superoxide dismutase; glutathione peroxidase and non-enzymatic antioxidant capacity are established indicators of antioxidant defences preventing oxidation of cell components"
                             "midbrain non enzymatic antioxidant capacity", # studyID=50:"superoxide dismutase; glutathione peroxidase and non-enzymatic antioxidant capacity are established indicators of antioxidant defences preventing oxidation of cell components"
                             "midbrain superoxide dismutase", # studyID=50:"superoxide dismutase; glutathione peroxidase and non-enzymatic antioxidant capacity are established indicators of antioxidant defences preventing oxidation of cell components"
                             "mitogenic response to concanavalin A (cell mediated immunity)", #studyID=84:proxy of immunity
                             "mitogenic response to polyinosinic-polycytidylic acid (cell mediated immunity)", #studyID=84:proxy of immunity
                             "nerve cell count", # studyID=119:authors interpret it as a the more, the better in terms of central nervous system
                             "number of stimulations to achieve eliptogenesis", # studyID=118:according to the authors low values are indicative of susceptibility to kindling epileoptogenesis
                             "number of vertebrae", # studyID=83:triploids have lower counts and are considered to be the low quality group (more deformities, etc). Thus, in this study it seems that it is better to have more vertebrae
                             "oestradiol concentration in diencephalon(pg/g)", # studyID=19:authors interpret low levels as: immune stress during late pregnancy reduced [...] negatively impacted [...] hippocampal function and central neurosteroid formation
                             "oestradiol concentration in hippocampus (pg/g)", # studyID=19:authors interpret low levels as: immune stress during late pregnancy reduced [...] negatively impacted [...] hippocampal function and central neurosteroid formation
                             "oestradiol concentration in medial prefrontal cortex (pg/g)", # studyID=19:authors interpret low levels as: immune stress during late pregnancy reduced [...] negatively impacted [...] hippocampal function and central neurosteroid formation
                             "optical density of 5-ht-immunoreactive fibres and varicosities in anterior hypothalymus",# studyID=61:"Hypothalamic 5-HT seems to diminish aggression, likely via inhibiting local arginine vasopressin actions". The authors consider aggression as a anxiety proxy
                             "optical density of 5-ht-immunoreactive fibres and varicosities in anterior hypothalymus (2h after resident intruder exposure)", # studyID=61:"Hypothalamic 5-HT seems to diminish aggression, likely via inhibiting local arginine vasopressin actions". The authors consider aggression as a anxiety proxy
                             "optical density of 5-ht-immunoreactive fibres and varicosities in basolateral amygdala", # studyID=61:"Hypothalamic 5-HT seems to diminish aggression, likely via inhibiting local arginine vasopressin actions". The authors consider aggression as a anxiety proxy
                             "optical density of 5-ht-immunoreactive fibres and varicosities in basolateral amygdala (2h after resident intruder exposure)", # studyID=61:"Hypothalamic 5-HT seems to diminish aggression, likely via inhibiting local arginine vasopressin actions". The authors consider aggression as a anxiety proxy
                             "optical density of 5-ht-immunoreactive fibres and varicosities in dosromedial hypothalamic nucleus", # studyID=61:"Hypothalamic 5-HT seems to diminish aggression, likely via inhibiting local arginine vasopressin actions". The authors consider aggression as a anxiety proxy
                             "optical density of 5-ht-immunoreactive fibres and varicosities in dosromedial hypothalamic nucleus (2h after resident intruder exposure)", # studyID=61:"Hypothalamic 5-HT seems to diminish aggression, likely via inhibiting local arginine vasopressin actions". The authors consider aggression as a anxiety proxy
                             "optical density of 5-ht-immunoreactive fibres and varicosities in lateral hypothalmic area", # studyID=61:"Hypothalamic 5-HT seems to diminish aggression, likely via inhibiting local arginine vasopressin actions". The authors consider aggression as a anxiety proxy
                             "optical density of 5-ht-immunoreactive fibres and varicosities in lateral hypothalmic area (2h after resident intruder exposure)", # studyID=61:"Hypothalamic 5-HT seems to diminish aggression, likely via inhibiting local arginine vasopressin actions". The authors consider aggression as a anxiety proxy
                             "optical density of 5-ht-immunoreactive fibres and varicosities in supraoptic nucleus", # studyID=61:"Hypothalamic 5-HT seems to diminish aggression, likely via inhibiting local arginine vasopressin actions". The authors consider aggression as a anxiety proxy
                             "optical density of 5-ht-immunoreactive fibres and varicosities in supraoptic nucleus (2h after resident intruder exposure)", # studyID=61:"Hypothalamic 5-HT seems to diminish aggression, likely via inhibiting local arginine vasopressin actions". The authors consider aggression as a anxiety proxy
                             "pectoral lepidotrichia", # studyID=82:the authors do not provide any predictions on this but we are assuming that stressed individuals should have less - lepidotrichia = fin spines and rays (https://en.wikipedia.org/wiki/Fish_fin)
                             "percentage exploration behaviour in aggresion test", # studyID=61:following what we have decided for similar variables (e.g. exploration)
                             "percentage time social behaviour in aggresion test", # studyID=61:mesaure of sociality
                             "percent time spent thrashing during predator cue", # studyID=87:we will assume (as done in similar variables, e.g. exploration), than more thrashin (more activty) is expected from less stressed individuals
                             "percent time stock female spent with control and treatment male", # studyID=66:proxy of male attractiveness
                             "plasma lysozyme concentrations", #studyID=84:proxy of immunity
                             "plasma proteins in natural conditions", # studyID=92:plasma proteins are considered as a measure of individual condition
                             "plasma proteins in supplemented conditons", # studyID=92:plasma proteins are considered as a measure of individual condition
                             "predator response (approached dummy predator)", # studyID=98:number of times that the focal individual got close to the predator dummy
                             "preoptic area size", # studyID=94:"found to be sexually dimorphic" and "are critical as integrative areas for mounting and intromission behaviour in gonadal males and for sexual receptivity in gonadal females". It is a difficult call but we are going to assume that stressed animals should develop smaller brains in general, which means that they should develop smaller brain areas too (difficult call though)
                             "pupal encapsulation rate", # studyID=47,111:"Higher melanization of the capsule corresponds to higher immune activation"
                             "red blood cell superoxide dismutase", #studyID=50:"superoxide dismutase; glutathione peroxidase and non-enzymatic antioxidant capacity are established indicators of antioxidant defences preventing oxidation of cell components"
                             "resting metabolic rate", # studyID=112,47,5: 112=authors predict a reduction in metabolic rate due to food shortage; 47 & 5=authors name resting metabolic rate as expensive, and thus, stressed animals should have lower metabolic rate
                             "scaled mass index", # studyID=97:proxy of body condition
                             "sternohyoideus cross section area", # studyID=81:structural element surrounding gill: authors interpret increase as negative and consequence of the need of more space for the gills. Authors argue so by saying that body shape can have consequences for hydrodynamics. Nonetheless, to standardize across studies, we are considering that all these type of morphological variables are expecte to be negatively correlated with stress
                             "sternohyoideus length", # studyID=81:structural element surrounding gill: authors interpret increase as negative and consequence of the need of more space for the gills. Authors argue so by saying that body shape can have consequences for hydrodynamics. Nonetheless, to standardize across studies, we are considering that all these type of morphological variables are expecte to be negatively correlated with stress
                             "thorax ratio", # studyID=85:defined as thorax mass/total body mass. The larger, the better in terms of flying (dispersion)
                             "time moving in novel environment", # studyID=45:as done for similar variables, we are going to assume that stressed animals move less
                             "time spent swimming in forced swim test", # the more, the more active stress-coping the individual shows (see studyID=61)
                             "time spent in open areas of maze test", # the more visiting to the open arms, the less anxious (more info in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3623971/)
                             "time spent climbing in forced swim test", # the more, the more active stress-coping the individual shows (see studyID=61)
                             "time swimming during 8 minute trial", # the more, the more active stress-coping the individual shows (see studyID=61)
                             "time swimming during 8min trial", # the more, the more active stress-coping the individual shows (see studyID=61)
                             "ultraviolet chroma", # studyID=71:proxy of male attractiveness
                             "upper caudal dermatotrichia", # studyID=82:the authors do not provide any predictions on this but we are assuming that stressed individuals should have less - lepidotrichia = fin spines and rays (https://en.wikipedia.org/wiki/Fish_fin)
                             "upper pharyngea jaw depth", # studyID=81:structural element surrounding gill: authors interpret increase as negative and consequence of the need of more space for the gills. Authors argue so by saying that body shape can have consequences for hydrodynamics. Nonetheless, to standardize across studies, we are considering that all these type of morphological variables are expecte to be negatively correlated with stress
                             "ventromedial nucleus of the hypothalymus size", # studyID=94:"found to be sexually dimorphic" and "are critical as integrative areas for mounting and intromission behaviour in gonadal males and for sexual receptivity in gonadal females". It is a difficult call but we are going to assume that stressed animals should develop smaller brains in general, which means that they should develop smaller brain areas too (difficult call though)
                             "visiting frequency", # studyID=59:the more visiting to the open arms, the less anxious (more info in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3623971/)
                             "visual exploration", # studyID=36:authors argue that stress negatively affects exploration, as found in other studies
                             "volume of central body", # studyID=119:authors interpret it as a the more, the better in terms of central nervous system
                             "volume of protocerebral neurophil", # studyID=119:authors interpret it as a the more, the better in terms of central nervous system
                             "whole brain metabolic capacity (umol/min/g tissue wet weight)" # studyID=99:self-defined
                           )
)


inversion.needed.2 <- c(inversion.needed,
                        c("aggression towards conspecific", # studyID=10:aggression of naive chicks (no previous experience with any other chick) in relation to experimentally increased or not cort levels. Authors predict (based on previous studies) that cort induces aggression, and we follow this expectation
                          "aterial pressure", # studyID=37:authors are not clear about their prediction. We are going to assume that stressed animals should have higher arterial pressure and higher heart rate (as pointed by the authors in the discussion)
                          "breeding time in natural conditions", # studyID=92:the earlier, the better
                          "breeding time in supplemented conditons", # studyID=92:the earlier, the better
                          "cerebellum protein carbonyls", # studyID=50:"protein carbonylation and is considered a reliable proxy of cellular oxidative protein damage"
                          "cycle length at 140 days ", # studyID=117: authors seem to interpret the longer cycle observed in stressed animals as evidence for premature ageing, but it isn't fully clear
                          "cylcle length at 1 year", # studyID=117: authors seem to interpret the longer cycle observed in stressed animals as evidence for premature ageing, but it isn't fully clear
                          "density of pyramidal neurons in left hipocampus", # studyID=115:"It is hypothesised that undernutrition of rats during the major generation and maturation phase of hippocampal pyramidal cells could lead to a diminution in their total numbers."
                          "feeding latency", # studyID=78:proxy of fitness where stressed animals should take longer to feed
                          "granulocyte:lymphocyte ratio", # studyID=13:an increase in the G:L ratio is a general avian stress response
                          "heart rate", # studyID=37:authors are not clear about their prediction. We are going to assume that stressed animals should have higher arterial pressure and higher heart rate (as pointed by the authors in the discussion)
                          "interferone (ifn)-y (cytokine) concentration in blood", # studyID=20:an increase is viewed by the authors as a negative
                          "interleukin 10 (cytokine) concentration in blood", # studyID=20:an increase is viewed by the authors as a negative
                          "interleukin 4 (cytokine) concentration in blood", # studyID=20:an increase is viewed by the authors as a negative
                          "interleukin 6 (cytokine) concentration in blood", # studyID=20:an increase is viewed by the authors as a negative
                          "latency to escape", # studyID=17:"latency to escape onto the [hidden] platform" so that the longer, the less cognitive abilities
                          "latency to exit begin exploring", # following the same logic that for studyID=59
                          "latency to feed from enclosed novel environment", # studyID=45:as done for similar variables, we are going to assume that stressed animals show higher latencies (particularly in a novel environment like this)
                          "latency to first appearence in maze", # studyID=59:the later, the more anxious (more info in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3623971/)
                          "latency to first egg sac", # studyID=120:time to first reproduction
                          "latency to leave starting box and enter maze", # the later, the more anxious (more info in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3623971/)
                          "latency to visit food patch under high threat", # studyID=41:"to determine whether stress during adolescence increased the latency to forage"
                          "midbrain protein carbonyls", #studyID=50:"protein carbonylation and is considered a reliable proxy of cellular oxidative protein damage"
                          "number of stops over 1m", # studyID=103:"was recorded as a measure of escape behaviour" and "greater number of stops accounts for overall slower speed"
                          "period of first incubation stage", # studyID=93:"period between spawning and first food uptake", overall, we can assume that the shorter the better (for the female, at least). However, note that the study is looking a match-mismatches between phenotype (fitness) and environment, which complicates the understanding of what effect is expected under stress
                          "period of first incubation stage in poor environment", # studyID=93:"period between spawning and first food uptake", overall, we can assume that the shorter the better (for the female, at least). However, note that the study is looking a match-mismatches between phenotype (fitness) and environment, which complicates the understanding of what effect is expected under stress
                          "period of second incubation stage", # studyID=93:"# period from end of 1st incubation phase to end of brood care", overall, we can assume that the shorter the better (for the female, at least). However, note that the study is looking a match-mismatches between phenotype (fitness) and environment, which complicates the understanding of what effect is expected under stress
                          "period of second incubation stage in poor environment", # studyID=93:"# period from end of 1st incubation phase to end of brood care", overall, we can assume that the shorter the better (for the female, at least). However, note that the study is looking a match-mismatches between phenotype (fitness) and environment, which complicates the understanding of what effect is expected under stress
                          "plasma testosterone levels", # studyID=69:adult testosterone levels for an experiment where eggs are given or not additional testoterone
                          "red blood cell protein carbonyls", # studyID=50:"protein carbonylation and is considered a reliable proxy of cellular oxidative protein damage"
                          "seconds spent begging", # studyID=10:"defined as frequent vertical movements of the head and pecking on the "beak" of the [adult-like] puppet". Authors expected more begging in stressed individuals, as found in previous studies, and we follow this expectation
                          "self directed behaviour (e.g. scratching, grooming)", # studyID=36:authors expect a positive correlation with stress, and classified this variable in the disturbance category
                          "testosterone concentration of eggs laid", # studyID=8:as for variable "plas testosterone levels", we will assume that stressed animals have more testosterone. This seems to be in accordance with that the authors of this paper say in the introduction
                          "testosterone in natural conditions", # studyID=92:testosterone is meant to correlate negatively with immunity. Thus, as in previous studies, we will assume that stressed animals have more testosterone.
                          "testosterone in supplemented conditons", # studyID=92:testosterone is meant to correlate negatively with immunity. Thus, as in previous studies, we will assume that stressed animals have more testosterone.
                          "thorax ratio (thorax dry mass/total dry mass)", # studyID=78:authors argue that "abdomen size correlates strongly with fecundity in insects", therefore, increased thorax ratio would reflect lower fecundity
                          "time taken to escape water", # the longer, the more passive stress-coping the individual shows (see studyID=61)
                          "total circulating androgens", # studyID=108:from the reading of the study, it seems that for females it is bad to have higher levels of androgens ("In mammals, administration of androgen to neonatal females results in an anovular (sterile) adult"). Note: not sure what should be the control treatment ("'hot' temperatures (32 ?C) resulting in mostly male offspring and relatively 'cold' temperatures (26 ?C) resulting in only female offspring"). EXCLUDE STUDY?
                          "total number of pyramidal nuerons in left hippocampus",# studyID=115:"It is hypothesised that undernutrition of rats during the major generation and maturation phase of hippocampal pyramidal cells could lead to a diminution in their total numbers."
                          "vasopressin-immunoreactive staining in nucleus circularis", # studyID=61:"at the level of the hypothalamus, arginine vasopressin was shown to increase anxiety"
                          "vasopressin-immunoreactive staining in nucleus circularis (2h after resident intruder exposure)", # studyID=61:"at the level of the hypothalamus, arginine vasopressin was shown to increase anxiety"
                          "vasopressin-immunoreactive staining in posterior part of the paraventricular hypothalymus", # studyID=61:"at the level of the hypothalamus, arginine vasopressin was shown to increase anxiety"
                          "vasopressin-immunoreactive staining in posterior part of the paraventricular hypothalymus (2h after resident intruder exposure)", # studyID=61:"at the level of the hypothalamus, arginine vasopressin was shown to increase anxiety"
                          "vasopressin-immunoreactive staining in the lateral hypothalamic area", # studyID=61:"at the level of the hypothalamus, arginine vasopressin was shown to increase anxiety"
                          "vasopressin-immunoreactive staining in the lateral hypothalamic area (2h after resident intruder exposure)", # studyID=61:"at the level of the hypothalamus, arginine vasopressin was shown to increase anxiety"
                          "vasopressin messenger ribonucleic-acid expression in bed nucleus of the stria terminalis ", # studyID=61:"at the level of the hypothalamus, arginine vasopressin was shown to increase anxiety"
                          "volume of subdivisions in left hippocampus"# studyID=115:"It is hypothesised that undernutrition of rats during the major generation and maturation phase of hippocampal pyramidal cells could lead to a diminution in their total numbers."
                        )
)


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


# # double-checking numbers to make sure everything is fine: all is good!
# length(levels(factor(db.full.red.sign[db.full.red.sign$origin!="STATS",c("specific.trait")])))==
#   length(no.inversion.needed.2)+length(inversion.needed.2)+length(inversion.unclear.2)
# setdiff(levels(factor(db.full.red.sign[db.full.red.sign$origin!="STATS",c("specific.trait")])),
#         c(no.inversion.needed.2,inversion.needed.2,inversion.unclear.2))
# setdiff(c(no.inversion.needed.2,inversion.needed.2,inversion.unclear.2),
#         levels(factor(db.full.red.sign[db.full.red.sign$origin!="STATS",c("specific.trait")])))


##############################################################
# Third round
##############################################################

# Assing the signs according to our classification

db.full.red.sign$sign.inversion.ours <- ""

db.full.red.sign[db.full.red.sign$specific.trait %in% no.inversion.needed.2 &
                   !(is.na(db.full.red.sign$mean.treat)),"sign.inversion.ours"]<-1
db.full.red.sign[db.full.red.sign$specific.trait %in% inversion.needed.2 &
                   !(is.na(db.full.red.sign$mean.treat)),"sign.inversion.ours"]<- -1
db.full.red.sign[db.full.red.sign$specific.trait %in% inversion.unclear.2 &
                   !(is.na(db.full.red.sign$mean.treat)),"sign.inversion.ours"]<- NA
db.full.red.sign[(is.na(db.full.red.sign$mean.treat)),"sign.inversion.ours"]<- NA

db.full.red.sign$sign.inversion.ours <- factor(db.full.red.sign$sign.inversion.ours)

summary(db.full.red.sign)


# extracting year of publication from the citation
db.full.red.sign$year <- as.integer(str_extract(db.full.red.sign$citation, regex("(\\d+)")))


# exporting clean data for effect size estimation
write.xlsx(db.full.red.sign,
           "data_re-extraction/clean_data/EyckDev_stress_clean_raw_data.xlsx",
           sheetName="Sheet1",col.names=TRUE, row.names=F,
           append=FALSE, showNA=TRUE, password=NULL)


# saving session information with all packages versions for reproducibility purposes
sink("data_re-extraction/clean_data/data_preparation_R_session.txt")
sessionInfo()
sink()

##############################################################
# Bonus check
##############################################################

# estimating the rate of typos purely affecting means, sds, and
# ns.

db.typos <- read.xlsx("data_re-extraction/re-extracted/001_EyckDev.stress_Data_FULL_TABLE_all_raw_data_re-extracted_VISUAL_manual_comp_typos.xlsx",
                      colNames=T,sheet = 1)

print(paste0(round((table(db.typos$mistake)[2]/nrow(db.typos))*100,0),"% of the effect sizes checked were affected by one or more typos/mistakes"))

