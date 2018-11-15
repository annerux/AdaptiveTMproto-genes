# November 13 2018
# Anne-Ruxandra Carvunis

# The table DATASET_6 was made by Brian using the matlab scripts and the raw data
# This script builds a table containing the resukts of the overexpression screens used to define "adaptive" orfs. 

# these were subject to a very deep inspection from the raw data in November 2018
# the raw data used can be found there:

# EXPID 28, SC-URA+GAL, 34 hours:
# /home/bhsu/yeast_beneficial_data/raw_data/BF1-MET_screen/BF-MET_SC_screen_6144/34hrs/set2/set1
# I checked that the colors on the plate side correspond to this condition(red, brown, purple); that the raw pixel counts in DATASET_6 JPEG_RESULTS table correspond to the raw image analysis; that the FITNESS results averages correspond to the JPEG_RESULTS averages; 

# EXPID 31, SC-URA+GAL+RAF, 34 hours:
#/home/bhsu/yeast_beneficial_data/raw_data/BF1-MET_screen/BF-MET_SC_screen_6144/34hrs/set3/set1
# I checked that the colors on the plate side correspond to this condition(red, brown, purple, blue); that the raw pixel counts in DATASET_6 JPEG_RESULTS table correspond to the raw image analysis; that the FITNESS results averages correspond to the JPEG_RESULTS averages;

# EXPID 33, SC-URA+CASE+ GAL+RAF, 34 hours:
#/home/bhsu/yeast_beneficial_data/raw_data/BF1-MET_screen/BF-MET_SC_screen_6144/34hrs/set3/set3
# I checked that the colors on the plate side correspond to this condition(red, brown, orange, purple, blue);that the raw pixel counts in DATASET_6 JPEG_RESULTS table correspond to the raw image analysis;that the FITNESS results averages correspond to the JPEG_RESULTS averages; NOTE THAT THE EXPS TABLE IN BRIAN_NEW INCORRECTLY ANNOTATES THIS AS SD_D-URA-MET+LYS+GAL+RAF BUT REALLY IT IS SC-URA+CASE+GAL+RAF,

# EXPID 91, SC-URA+CASE+ GAL, 35 hours:
# /home/bhsu/yeast_beneficial_data/raw_data/BF5_CASEGAL_screen/BF_CASEGAL_6144/35hrs/set1/set1
# I checked that the colors on the plate side correspond to this condition(red, brown, orange, purple); that the raw pixel counts in DATASET_6 JPEG_RESULTS table correspond to the raw image analysis; the FITNESS results averages DO NOT correspond to the JPEG_RESULTS averages because of plate manipulation errors that were taken into account by the FLAGS_ALL table.

# EXPID 91, SD-A-URA+GAL+RAF, 58 hours:
# /home/bhsu/yeast_beneficial_data/raw_data/BF5_CASEGAL_screen/BF_CASEGAL_6144/35hrs/set1/set3
# it is very hard to distinguish the colors in the plate sides but i see red brown black blue and purple; I do see an extra hint of red or orange, but since these plates are clearly SD I think it may indicate SD-A or something, it cannot indicate Case. I am confident these plates are SD by looking at the colonies. I checked that the raw pixel counts in DATASET_6 JPEG_RESULTS table correspond to the raw image analysis; the FITNESS results averages DO NOT correspond to the JPEG_RESULTS averages because of plate manipulation errors that were taken into account by the FLAGS_ALL table.


# CONNECT DB
library(RMySQL)
con <- dbConnect(dbDriver("MySQL"),user='carvunis',password='mabiche25',dbname='ANNE',host='paris.csb.pitt.edu')  

general_outdir<-"/Users/annerux/Workspace/cerevisiae/txt/screen_analyses/Protogenes/ProtogeneFitnessFigurePanels/submission/" 
date<-"20181113"



suptablequery<-"select distinct a.orf_name, hours, n, colony_size normalized_cs, q_cs, effect_cs,exp_id from brian_031918.DATASET_6 a, ANNE.SUMMARY_2012 b where a.orf_name = b.orf_name and exp_id in (28,31,91,33, 93) and a.orf_name !='BF_control'"
suptable6conditions<-dbGetQuery(con,suptablequery)
suptable6conditions$exp_environment<-sapply(suptable6conditions$exp_id, function(x)
ifelse(x==33, "N++,C++", ifelse(x==28, "N+, C+", ifelse(x==31,"N+, C++", ifelse(x==91, "N++,C+", ifelse(x==93, "N-,C++"))))) )

suptable6conditions$effect_cs<-sapply(suptable6conditions$effect_cs, function(x)
ifelse(x=="neutral", "unchanged", ifelse(x=="beneficial", "increased", ifelse(x=="deleterious", "decreased", "other"))))

write.csv(suptable6conditions, paste(general_outdir,"Data3.csv",sep=""))

dbWriteTable(con, "BARFLEX_6CONDITIONS_SUPTABLE", suptable6conditions, row.names=FALSE, overwrite=TRUE ) # later renamed 5 conditions

lapply(dbListConnections( dbDriver( drv = "MySQL")), dbDisconnect)

