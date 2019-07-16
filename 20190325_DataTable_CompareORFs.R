# March 25 2019
# Anne-Ruxandra Carvunis

# this script builds a table that contains all the info for redoing all the analyses in the paper, except the 5 condition data results and the data for additional membrane analyses.  


# CONNECT DB
library(RMySQL)
con <- dbConnect(dbDriver("MySQL"),user='carvunis',password='xxx',dbname='ANNE',host='paris.csb.pitt.edu')  

general_outdir<-"/Users/annerux/Workspace/cerevisiae/txt/screen_analyses/Protogenes/ProtogeneFitnessFigurePanels/" 


###########################################################
# list of all ORFs 
allgenes_q<- "select distinct a.orf_name, quality from ANNE.REF_GENES a, ANNE.SUMMARY_2012 b where a.orf_name=b.orf_name"
allgenes<-dbGetQuery(con, allgenes_q)[,1:2]

# Emerging versus Established ORFs
# phylostratigraphy + no string protein-coding signatures
protogenes_q<-"select * from ANNE.PROTOGENES where orf_name like '%y%'  and selected+longer+translated<3"
protogenes<-dbGetQuery(con, protogenes_q)
# identify ORFs with >50% conserved AA in Kud or Bay based on synteny analyses
synteny_table<-read.csv("/Users/annerux/Workspace/cerevisiae/txt/synteny project/20190223_heatmapdata.csv")
too_conserved<-synteny_table[round(synteny_table$Kud.overlap.AA..Identitiy,2)>0.50 & round(synteny_table$Bay.overlap.AA..Identitiy,2)>0.50,]$orf_name
protogenes <- protogenes[! (protogenes[,1] %in% too_conserved),]
protogene<-sapply(allgenes[,1], function(x) ifelse(x %in% protogenes[,1], "Emerging ORFs","Established ORFs"))

# Barflex space
bf_q<-"select distinct orf_name from anne_screen.BARFLEX_SPACE_AGAR"
barflex_space_agar<-dbGetQuery(con,bf_q)[,1]
barflex_space<-sapply(allgenes[,1], function(x) ifelse(x %in% barflex_space_agar, "yes","no"))

## list of genes and proto-genes in the collection
orf_table<-data.frame(allgenes , protogene, barflex_space)
orf_table $protogene <-factor(orf_table $protogene, levels=c("Emerging ORFs","Established ORFs"))

# orf fitness from our SCGAL overexpression experiment
 fitness_q<-" select distinct a.orf_name, a.effect_cs overexpression_relative_fitness from
 brian_031918.DATASET_6 a where exp_id = 28 "
 fitness<-dbGetQuery(con, fitness_q)
 fitness[ fitness[,2]=="neutral",2]<-"unchanged"
 fitness[ fitness[,2]=="deleterious",2]<-"decreased"
 fitness[ fitness[,2]=="beneficial",2]<-"increased"
 orf_table<-merge(orf_table, fitness, by="orf_name",all.x=TRUE) 
 
# fitness from Douglas et al
Douglas_q<-"select distinct orf_name, 20gen overexpression_competitive_fitness from anne_phenotypes.BARFLEX_FITNESS_LIQUID"
Douglas<-dbGetQuery(con, Douglas_q)[,1:2]
Douglas <- subset(Douglas, !duplicated(Douglas[,1]))
orf_table<-merge(orf_table, Douglas, by="orf_name", all.x=TRUE)

# mRNA expression level
mrna_q<-"select distinct orf_name, mrna_rich expr_level from ANNE.ingolia_singlereads"
expr_level<-dbGetQuery(con, mrna_q)
expr_level[,2]<-expr_level[,2]+0.001
expr_level[,2]<-log(expr_level[,2])
orf_table<-merge(orf_table, expr_level, by="orf_name", all.x=TRUE)

# length
length_q<-" select distinct orf_name, (coor2-coor1 +1) Length from ANNE.REF_GENES"
length<-dbGetQuery(con,length_q)
orf_table <-merge(orf_table,length,all.x=TRUE)

# disorder from Carvunis 2012 paper - vsl2 predictions
disorder<-read.csv("/Users/annerux/Workspace/DFCI/PROTOGENES/TABLES/sgrp_anne/percent_disorder.csv")
orf_table <-merge(orf_table, disorder,all.x=TRUE)

# overlap
overlap_q<- "select * from ANNE.OVERLAPPING_ORFS "
overlap<-dbGetQuery(con, overlap_q)[,1]
orf_table $overlap<-sapply(orf_table $orf_name, function(x) ifelse(x %in% overlap,1,0))

# single deletion fitness from costanzo 2016
smf_file<-"/Users/annerux/Workspace/cerevisiae/txt/external_datasets/costanzo/TheCellMapJune52018/strain_ids_and_single_mutant_fitness.csv"
smf<-read.csv(smf_file, header=TRUE)[,1:7]
colnames(smf)<-c("strain_id", "orf_name", "allele_name", "smf_26", "sd_26", "smf_30", "sd_30")
smf_average30<-aggregate(smf$smf_30 ~ smf$orf_name, FUN=mean)
colnames(smf_average30)<-c("orf_name", "loss_fitness")
orf_table <-merge(orf_table, smf_average30, by="orf_name", all.x=TRUE)

# membrane, SP , gravy, aromatocity and GC from Nikos
nikos<-read.csv("/Users/annerux/Workspace/cerevisiae/txt/external_datasets/NikosProteinStructure/20181004_GeneData.csv")
nikos$percent_phobius<-nikos$exp_no_transm_res_PHOBIUS / nikos$length
nikos$percent_tmhmm<-nikos$exp_no_transm_res_TMHMM / nikos$length

nikos<-nikos[,c(1,11, 13, 14, 18, 19)]
orf_table <-merge(orf_table, nikos, by="orf_name", all.x=TRUE)

# isolate data
strain_file<-"/Users/annerux/Workspace/cerevisiae/txt/homology/20181006_proto_strain_data_refgenes.txt"
strains<-read.table(strain_file, header=TRUE)
strains<-strains[strains$missing==0, c("protogene_id", "good_align", "conserved_orf", "pi")]
colnames(strains)<-c("orf_name", "strain_goodalign", "strained_conserved_orf", "strain_pi")
orf_table <-merge(orf_table, strains, by="orf_name", all.x=TRUE)
orf_table[orf_table$length %%3 !=0,c("strain_goodalign", "strained_conserved_orf", "strain_pi")]<-NA


# EXPORT
# text
write.csv(orf_table, paste(general_outdir, "submission/Data1.csv", sep=""))
#DB

dbWriteTable(con, "DATATABLE_ORF_DESC", orf_table, row.names=FALSE, overwrite=TRUE )

lapply(dbListConnections( dbDriver( drv = "MySQL")), dbDisconnect)

