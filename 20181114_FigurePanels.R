# November 14 2018
# Anne-Ruxandra Carvunis

# this script reproduces all the figure panels resulting from computational analyses in Vakirlis, Hsu et al from the supplementary Data tables. 


# LIBRARIES
library(ggplot2)
library(plyr)
library(ggpubr)

  
# FOLDER
general_outdir<-"Workspace/cerevisiae/txt/screen_analyses/Protogenes/ProtogeneFitnessFigurePanels/submission/" # CHANGE THIS LINE TO INDICATE WHERE YOU STORED THE DATA FILES 1-5


# DATA SOURCES

# description of ORFs used in analyses 
orf_table<-read.csv(paste(general_outdir, "Data1.csv",sep=""),header=TRUE)
orf_table $protogene <-factor(orf_table $protogene, levels=c("proto-gene","gene"))
orf_table $overexpression_relative_fitness <-factor(orf_table $overexpression_relative_fitness, levels=c("increased","decreased","unchanged"))

# individual normalized colony sizes for 1 environmental condition (Fig. 3B)
colony_table<-read.csv(paste(general_outdir, "Data2.csv",sep=""),header=TRUE)

# results of overexpression screen in 5 environmental conditions 
fitness_table<-read.csv(paste(general_outdir, "Data3.csv",sep=""),header=TRUE)

# transmembrane analyses (related to Figs 5, S6 and S7)
tm_table<-read.csv(paste(general_outdir, "Data4.csv",sep=""),header=TRUE)
tm_table $type <-factor(tm_table $type, levels=c("genes","proto-genes","intergenes","non-genic ORFs"))
ngo_intergene_table<-read.csv(paste(general_outdir, "Data5.csv",sep=""),header=TRUE)

#######################################
# Figure 2B: fitness cost of ORF loss #
#######################################

# data
df<-orf_table[!is.na(orf_table$loss_fitness) & orf_table$overlap==0,c("protogene","loss_fitness")]

# statistical comparison: Mann Whitney
wilcoxtest<-wilcox.test(df[df $protogene=="gene","loss_fitness"],df[df $protogene=="proto-gene","loss_fitness"])
pval<-wilcoxtest$p.val

# statistical comparison: Fisher's exact with fitness cutoff of 0.9 
n_pg_effect<-nrow(df[df$protogene=="proto-gene" & df$loss_fitness <0.9,]) 
n_pg_noeffect<-nrow(df[df$protogene=="proto-gene" & df$loss_fitness>=0.9,])
n_g_effect<-nrow(df[df$protogene=="gene" & df$loss_fitness <0.9,]) 
n_g_noeffect<-nrow(df[df$protogene=="gene" & df$loss_fitness>=0.9,])
pval<-fisher.test(matrix(c(n_pg_effect, n_pg_noeffect, n_g_effect, n_g_noeffect),2,2))

# cumulative plot
plot_xlabel<-"Fitness upon loss"
plot_ylabel<- "ORFs (cumulative frequency)"
image<-ggplot(df, aes(x= loss_fitness,color=protogene, fill=protogene))+stat_ecdf(geom="step", na.rm=TRUE)+ stat_ecdf(aes(ymin=0,ymax=..y..),
geom="ribbon", na.rm=TRUE)+scale_color_manual(values=c("deepskyblue4","black"))+ scale_fill_manual(values=c("lightblue","NA"))+labs(x=plot_xlabel,y=plot_ylabel)+ geom_vline(xintercept=0.9, colour="red") + theme_bw()+theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8),  	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none") 
ggsave(paste(general_outdir,"Fig2B_CompareGeneProtogenes_CostOfLoss.pdf",sep=""),width=3.5,height=2.5)

#########################################################
# Figure 2C: fixation of ORf sctructure across isolates #
#########################################################

df<-orf_table[!(is.na(orf_table $strain_goodalign)),c("protogene","strain_goodalign", "strained_conserved_orf",  "strain_pi")]

# plot
plot_xlabel<-"% isolates with intact ORF structure"
plot_ylabel<- "ORFs (cumulative frequency)"
image<-ggplot(df, aes(x= strained_conserved_orf/2022*100,color=protogene, fill=protogene))+stat_ecdf(geom="step", na.rm=TRUE)+ stat_ecdf(aes(ymin=0,ymax=..y..),
geom="ribbon", na.rm=TRUE)+scale_color_manual(values=c("deepskyblue4","black"))+ scale_fill_manual(values=c("lightblue","NA"))+labs(x=plot_xlabel,y=plot_ylabel) + geom_vline(xintercept=90, colour="red")+ theme_bw()+theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8),  	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none") + coord_cartesian(xlim=c(0,100), , expand=FALSE)
ggsave(paste(general_outdir,"Fig2C_CompareGeneProtogenes_StrainORFStructure_ECDF.pdf",sep=""),width=3.5,height=2.5)

# intact orf structure in <90% isolates
df$strict<-ifelse(df$strained_conserved_orf/2022*100<90,1,0)
table(df[,c("protogene","strict")])

###################################################
# Figure 2D: nucleotide diversity across isolates #
###################################################

df<-orf_table[!(is.na(orf_table $strain_goodalign)),c("protogene","strain_goodalign", "strained_conserved_orf",  "strain_pi")]

# statistical comaprison: Mann Whitney
wilcoxtest<-wilcox.test(df[df $protogene=="gene","strain_pi"],df[df $protogene=="proto-gene","strain_pi"])
pval<-wilcoxtest$p.val

# plot
mu <- ddply(df, "protogene", summarise, grp.mean=mean(strain_pi,na.rm=TRUE))
plot_xlabel<-"Nucleotide diversity across isolates"
plot_ylabel<- "ORFs (density)"
image<-ggplot(df, aes(x= strain_pi,color=protogene,fill=protogene))+geom_density(na.rm=TRUE)+
   geom_vline(data=mu, aes(xintercept=grp.mean, color=protogene), linetype="dashed")+scale_color_manual(values=c("deepskyblue4","black"))+scale_fill_manual(values=c( "lightblue","NA"))+labs(x=plot_xlabel,y=plot_ylabel)+ theme_bw()+theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8),  	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")
ggsave(paste(general_outdir,"Fig2D_CompareGeneProtogenes_StrainNucleotideDiversity.pdf",sep=""),width=3.5,height=2.5)


######################################################
# Figure 2E: competitive fitness upon overexpression #
######################################################
# data
df<-orf_table[orf_table$barflex_space=="yes",c("protogene","overexpression_competitive_fitness")]

#  statistical comparison: Mann Whitney
wilcoxtest<-wilcox.test(df[df $protogene=="gene","overexpression_competitive_fitness"],df[df $protogene=="proto-gene","overexpression_competitive_fitness"])
pval<-wilcoxtest$p.val

# plot
mu <- ddply(df, "protogene", summarise, grp.mean=mean(overexpression_competitive_fitness,na.rm=TRUE))
plot_xlabel<-"Competitive fitness upon overexpression"
plot_ylabel<- "ORFs (density)"
image<-ggplot(df, aes(x= overexpression_competitive_fitness,color=protogene,fill=protogene))+geom_density(na.rm=TRUE)+   geom_vline(data=mu, aes(xintercept=grp.mean, color=protogene), linetype="dashed")+scale_color_manual(values=c("deepskyblue4", "black"))+scale_fill_manual(values=c("lightblue", "NA"))+labs(x=plot_xlabel,y=plot_ylabel) + theme_bw()+theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8),  	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")   
ggsave(paste(general_outdir,"Fig2E_CompareGeneProtogenes_CompetitiveFitness.pdf",sep=""),width=3.5,height=2.5)


##############################################
# Figure 3B: colony size upon overexpression #
##############################################

df<-colony_table

# category for each colony
df$category<-rep("other", nrow(df))
df[df$orf_name=="BF_control",]$category<-"reference"
df[df$orf_name %in% orf_table[orf_table$protogene=="proto-gene",]$orf_name,]$category<-"proto-gene"
df[df$orf_name %in% orf_table[orf_table$protogene=="gene",]$orf_name,]$category<-"gene"
df $category <-factor(df $category, levels=c("proto-gene","gene", "reference"))

#  statistical comparison: Mann Whitney between genes and proto-genes
wilcoxtest<-wilcox.test(df[df $category=="gene","colony_size"],df[df $category=="proto-gene","colony_size"])
pval<-wilcoxtest$p.val

#  statistical comparison: Mann Whitney between reference strain and proto-genes
wilcoxtest<-wilcox.test(df[df $category=="reference","colony_size"],df[df $category=="proto-gene","colony_size"])
pval<-wilcoxtest$p.val

#Plot
mu <- ddply(df, "category", summarise, grp.mean=mean(colony_size,na.rm=TRUE))
plot_xlabel<-"Normalized colony size"
plot_ylabel<- "Colonies (density)"

image<-ggplot(df, aes(x= colony_size,color=category,fill=category))+geom_density(na.rm=TRUE)+
   geom_vline(data=mu, aes(xintercept=grp.mean, color= category), linetype="dashed")+scale_color_manual(values=c("deepskyblue4","black","grey"))+scale_fill_manual(values=c( "lightblue","NA","NA"))+labs(x=plot_xlabel,y=plot_ylabel)+ theme_bw()+theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8),  	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none") 
ggsave(paste(general_outdir,"Fig3B_Experimental_GeneProtogenesReference_Colonysize.pdf",sep=""),width=3.5,height=2.5)


############################################################################
# Figure 3C: fraction of orfs with relative overexpression fitness effects #
############################################################################

df<-orf_table[orf_table$barflex_space=="yes",c("protogene","overexpression_relative_fitness")]

# calculate proportions and standard error of proportions
a<-data.frame(table(df))
totals<-rep(table(df$protogene),3)
a$groupsize<-totals
sder<-sqrt(a$Freq/a$groupsize*(1-a$Freq/a$groupsize)/a$groupsize)
a$err<-sder
max<-a$Freq/a$groupsize + a$err
a$top<-max
min<-a$Freq/a$groupsize - a$err
 a$bottom<-min
 proportion<-a$Freq/a$groupsize
 a$fraction<-proportion

# Plot 
plot_xlabel<-"Relative fitness"
plot_ylabel<- "Fraction of ORFs"  
image<-ggplot(a,aes(x= overexpression_relative_fitness,y=fraction,ymin=bottom,ymax=top,fill=protogene, color=protogene))+geom_bar(stat="identity",position=position_dodge())+geom_errorbar(position = position_dodge(0.9),width=0.25)+scale_color_manual(values=c("deepskyblue4","black"))+scale_fill_manual(values=c( "lightblue","white"))+labs(x=plot_xlabel,y=plot_ylabel)+ theme_bw()+theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8),  	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none") 
ggsave(paste(general_outdir,"Fig3C_Experimental_GeneProtogenes_RelativeFitness_Fraction.pdf",sep=""),width=3.5,height=2.5)

# range of effect sizes (average of replicates)
increased_colonies<-colony_table[colony_table$orf_name %in% orf_table[orf_table $overexpression_relative_fitness=="increased",]$orf_name,]
mean_colonysizes<-ddply(increased_colonies, "orf_name", summarise, grp.mean=mean(colony_size,na.rm=TRUE))
min(mean_colonysizes $grp.mean)
max(mean_colonysizes $grp.mean)


##############################################################################
# Figure 3D: Odds Ratio of orfs with relative overexpression fitness effects #
##############################################################################

df<-orf_table[orf_table$barflex_space=="yes",c("protogene","overexpression_relative_fitness")]

# calculate odds ratios and confidence intervals and pvalues
a<-data.frame(table(df))
totals<-rep(table(df$protogene),3)
a$groupsize<-totals
a$minus<-a$groupsize-a$Freq
benef<-fisher.test(a[a$overexpression_relative_fitness =="increased",c("Freq","minus")])
neutral<-fisher.test(a[a$overexpression_relative_fitness =="unchanged",c("Freq","minus")])
delet<-fisher.test(a[a$overexpression_relative_fitness =="decreased",c("Freq","minus")])
oddsratio<-data.frame(overexpression_relative_fitness =c("increased","unchanged","decreased"),or=c(benef$estimate,neutral$estimate,delet$estimate))
oddsratio$top<-c(benef$conf.int[2],neutral$conf.int[2],delet$conf.int[2])
oddsratio$bottom<-c(benef$conf.int[1],neutral$conf.int[1],delet$conf.int[1])
oddsratio $overexpression_relative_fitness <-factor(oddsratio $overexpression_relative_fitness, levels=c("increased","decreased","unchanged"))

# plot
plot_xlabel<-"Relative fitness"
plot_ylabel<- "Odds ratio" 
image<-ggplot(oddsratio,aes(x= overexpression_relative_fitness,y=or,ymin=bottom,ymax=top))+geom_point(stat="identity",  shape=21, size=2, stroke=1, fill = "white")+geom_errorbar(width=0.25)+geom_hline(yintercept=1,linetype="dashed",color="red")+ labs(x=plot_xlabel,y=plot_ylabel)+ theme_bw()+theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8),  	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")
ggsave(paste(general_outdir,"Fig3D_Experimental_GeneProtogenes_RelativeFitness_OddsRatio.pdf",sep=""),width=3.5,height=2.5,  useDingbats=FALSE)

###############################################################################################
# Figure 3E: Odds Ratio of orfs with relative overexpression fitness effects - 5 environments #
###############################################################################################

df<-merge(fitness_table, orf_table[,c("orf_name","protogene")], all.x=TRUE)

# calculate odds ratios and confidence intervals and pvalues
a<-data.frame(table(df[,c("exp_environment","protogene", "effect_cs")]))
totals<-table(df[,c("exp_environment","protogene")])
b<-merge(a,totals,by=c("exp_environment","protogene"),all.x=TRUE)
colnames(b)<-c("exp_environment",  "protogene",   "effect_cs", "number", "total")
b$minus<-b$total - b$number
protos<-b[b$protogene=="proto-gene",c("exp_environment" , "protogene",   "effect_cs", "number", "minus")]
genes<-b[b$protogene=="gene",c("exp_environment" , "protogene",   "effect_cs", "number", "minus")]
numbers<-merge(protos,genes,by=c("exp_environment","effect_cs"))
colnames(numbers)<-c("exp_environment" ,   "effect_cs", "proto", "protonumber", "protominus","gene","genenumber","geneminus")
fisher_test_function<-function(x){
	mat<-matrix(as.numeric(c(x[4],x[5],x[7],x[8])),2,2)
	f<-fisher.test(mat)
	return(c(f$estimate, f$conf.int[1],f$conf.int[2], f$p.value))
}
or<-apply(numbers,1,function(x) fisher_test_function(x) )
oddsratios<-data.frame(cbind(numbers[,c("exp_environment","effect_cs")],t(or)))
colnames(oddsratios)<-c("exp_environment",   "effect_cs" , "or", "bottom","top","pval")
oddsratios $effect_cs <-factor(oddsratios $effect_cs, levels=c("increased","decreased","unchanged"))

# plot
forplotting<-oddsratios[oddsratios$effect_cs %in% c("increased","decreased"),]
#forplotting$exp_environment <- factor(forplotting $exp_id, levels = condorder)
plot_xlabel<-"Environmental condition"
plot_ylabel<- "Odds ratio" 
image<-ggplot(forplotting,aes(x= exp_environment,y=or,ymin=bottom,ymax=top))+geom_point(stat="identity", shape=21, size=2, stroke=1, fill = "white") + facet_wrap (~ effect_cs, nrow=2)+geom_errorbar(width=0.25)+geom_hline(yintercept=1,linetype="dashed",color="red")+ labs(x=plot_xlabel,y=plot_ylabel)+ theme_bw()+theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8),  	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")
ggsave(paste(general_outdir,"Fig3E_Experimental_GeneProtogenes_RelativeFitness_OddsRatio_5Environments.pdf",sep=""),width=4,height=3.5,  useDingbats=FALSE)


#############################
# Figure 3F: Adaptive ORFs #
############################

df<-merge(fitness_table, orf_table[,c("orf_name","protogene","loss_fitness", "strained_conserved_orf", "strain_pi")], all.x=TRUE)
adaptive_protogenes<-unique(df[df$protogene=="proto-gene" & df$effect_cs=="increased",]$orf_name)
adaptive_genes<-unique(df[df$protogene=="gene" & df$effect_cs=="increased",]$orf_name)

# counts and error bars
n_adaptive_protogenes<-length(adaptive_protogenes)
n_adaptive_genes<-length(adaptive_genes)
n_protogenes<-length(unique(df[df$protogene=="proto-gene",]$orf_name))
n_genes<-length(unique(df[df$protogene=="gene",]$orf_name))
fraction_protogenes<-n_adaptive_protogenes/n_protogenes
fraction_genes<-n_adaptive_genes/n_genes
err_protogenes<-sqrt(fraction_protogenes * (1-fraction_protogenes)/n_protogenes)
err_genes<-sqrt(fraction_genes * (1-fraction_genes)/n_genes)

# plot
plot_df<-data.frame(groups=c("proto-genes","genes"), adaptive=c(fraction_protogenes, fraction_genes),f_min=c(fraction_protogenes-err_protogenes, fraction_genes-err_genes),f_max=c(fraction_protogenes+err_protogenes, fraction_genes+err_genes))
plot_df $groups <-factor( plot_df $groups, levels=c("proto-genes","genes"))
image<-ggplot(plot_df,aes(x=groups,y= adaptive,ymin= f_min,ymax= f_max,fill=groups, color=groups))+geom_bar(stat="identity",position=position_dodge())+geom_errorbar(position = position_dodge(0.9), width=0.25)+scale_color_manual(values=c("deepskyblue4","black"))+scale_fill_manual(values=c( "lightblue","white"))+labs(x="",y="Adaptive")+ theme_bw()+theme(axis.text.y=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.x=element_text(family="Helvetica", face="bold", colour="black", size=10), 	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none") 
ggsave(paste(general_outdir,"Fig3F_Experimental_GeneProtogene_Adaptive.pdf",sep=""),width=2.5,height=3.5, useDingbats=FALSE)

# statistical comparison: selected effects in adaptive proto-genes
pg_df<-unique(df[df$protogene=="proto-gene", c("orf_name","loss_fitness","strained_conserved_orf", "strain_pi")])
pg_df $adaptive<- ifelse(pg_df$orf_name %in% adaptive_protogenes,1,0)
# fitness upon loss
pg_df $sick<-ifelse(pg_df $loss_fitness<0.9,1,0)
fisher.test(as.matrix(table(pg_df[,c("adaptive", "sick")])))
# orf intactness in isolates
pg_df$intact<-ifelse(pg_df $strained_conserved_orf <0.9 * 2022 ,0,1)
fisher.test(as.matrix(table(pg_df[,c("adaptive", "intact")])))
# nucleotide diversity
wilcox.test(pg_df[pg_df$adaptive==1,"strain_pi"], pg_df[pg_df$adaptive==0,"strain_pi"], na.rm=TRUE)

# how many proto-genes have genetic tradeoffs
adaptive_orfs<-df[df$orf_name %in% c(adaptive_genes, adaptive_protogenes),]
n_adaptive_protogenes_withtradeoff<-length(unique(adaptive_orfs[adaptive_orfs$protogene=="proto-gene" & adaptive_orfs$effect_cs == "decreased",]$orf_name))


#######################################################
# Figure S2: Adaptive proto-genes across environments #
#######################################################

df<-merge(fitness_table, orf_table[,c("orf_name","protogene")], all.x=TRUE)
adaptive_protogenes_df<-df[df$protogene=="proto-gene" & df$effect_cs=="increased",c("orf_name","exp_environment")]

# plot counts per environment (Fig S2A)
plot_x_label<-"Environmental condition"
plot_y_label<-"# adaptive proto-genes"
count_image<-ggplot(adaptive_protogenes_df, aes(exp_environment))+ geom_bar(color='deepskyblue4', fill='lightblue')+theme_bw()+labs(x=plot_x_label,y=plot_y_label) + theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.title.x=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.y=element_text(family="Helvetica", face="bold", colour="black", size=10), legend.position="none") 
ggsave(paste(general_outdir,"FigS2A_AdaptiveProtogenes_PerEnvironment_Counts.pdf",sep=""),width=3.5,height=2.5, useDingbats=FALSE)

# distribution of # environements (Fig S2B)
 plot_x_label<-"# environmental conditions"
plot_y_label<-"# adaptive proto-genes"
distrib<-table(adaptive_protogenes_df$orf_name)
plot_df<-as.data.frame(distrib)
distrib_image<-ggplot(plot_df[plot_df$Freq>0,], aes(x=Freq))+ geom_histogram(binwidth=1, color='deepskyblue4', fill='lightblue')+theme_bw()+labs(x=plot_x_label,y=plot_y_label) + theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.title.x=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.y=element_text(family="Helvetica", face="bold", colour="black", size=10), legend.position="none") 
ggsave(paste(general_outdir,"FigS2B_AdaptiveProtogenes_PerEnvironment_Distribution.pdf",sep=""),width=3.5,height=2.5, useDingbats=FALSE)

#############################################################################################################################
# Figures 4A, S3B-D, S4 A-D: Statistical comparison of features separating neutral from adaptive and deleterious ORFs#
# note that Figure 4A and S3A are the same
############################################################################################################################
#### data
# features
features_column_names<-c("expr_level","Length","Translation","PI","GRAVY","Aromo","CAI","Codon.Bias","GC","INSTABILITY.INDEX..II.","ALIPHATIC.INDEX","negative_charge","positive_charge","percent_disorder","percent_phobius","signal_peptide","tAI_first_10","mfe_first_10_cods", "overlap", "Cys")
pretty_names<-c("Expression level","Length","Translation","Isoelectric point","Hydropathicity","Aromaticity","CAI","Codon bias","GC","Instability", "Aliphaticity", "Negative charge","Positive charge", "Disorder", "Transmembrane","Signal peptide","tAI", "Minimum Free Energy", "Overlap", "Cysteine")
feature_df<-orf_table[orf_table$barflex_space=="yes",c("orf_name","protogene", features_column_names)]

# categories
fitness_df<-merge(fitness_table, orf_table[,c("orf_name","protogene")], all.x=TRUE)
adaptive_orfs<-unique(fitness_df[fitness_df $effect_cs=="increased",]$orf_name)
deleterious_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & fitness_df $effect_cs=="decreased",]$orf_name)
neutral_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & !(fitness_df $orf_name %in% deleterious_orfs),]$orf_name) 
feature_df $fitness_category<-sapply(feature_df $orf_name, function(x) ifelse(x %in% adaptive_orfs,"adaptive",ifelse(x %in% deleterious_orfs, "deleterious", ifelse(x %in% neutral_orfs,"neutral","others"))))

#### function to build a statistics data frame that saves pvalue and effect size per feature 
stats_feature<-function(columnname, orf_categ, fitness_categ){
	# data to analyze
	minidf<-feature_df[feature_df $protogene == orf_categ & feature_df $fitness_category %in% c("neutral", fitness_categ), c("fitness_category",columnname)]
	# pval and effect size
	pval<-wilcox.test(minidf[minidf$fitness_category==fitness_categ, columnname], minidf[minidf$fitness_category=="neutral", columnname])$p.value
	effectsize<-(mean(minidf[minidf$fitness_category==fitness_categ, columnname], na.rm=TRUE) - mean(minidf[minidf$fitness_category=="neutral", columnname], na.rm=TRUE))/sd(minidf[minidf$fitness_category=="neutral", columnname], na.rm=TRUE)
	# return results
	return(c(pval, (-log(pval)),effectsize))
}


#### calculate pvalues and effect sizes
statistical_results<-NULL
molecularfeature<-NULL
protogene_vector<-NULL
fitness_vector<-NULL
pval_vector<-NULL
logpval_vector<-NULL
effectsize_vector<-NULL

for (columnnumber in c(1:length(features_column_names)) ){
	column<-features_column_names[columnnumber]
	prettycolumnname<-pretty_names[columnnumber]
	for( orf in c("gene","proto-gene")){
		for (fitness in c("deleterious","adaptive")){
			molecularfeature<-c(molecularfeature, prettycolumnname)
			protogene_vector <-c(protogene_vector, orf)
			fitness_vector<-c(fitness_vector, fitness)
			statistics<-stats_feature(column, orf, fitness)
			pval_vector<-c(pval_vector, statistics[1])
			logpval_vector <-c(logpval_vector, statistics[2])
			effectsize_vector<-c(effectsize_vector, statistics[3])
		}
	}
}

stat_df<-data.frame(molecularfeature, protogene_vector, fitness_vector,pval_vector, logpval_vector,effectsize_vector, stringsAsFactors=F)

#### Plots

### proto-genes
## adaptive versus neutral
# pvalues - Fig. 4A and S3A (the same!)
d<-stat_df[stat_df$protogene_vector == "proto-gene" & stat_df$fitness_vector =="adaptive",]
d$molecularfeature<-factor(d$molecularfeature, levels=d[order(d$logpval_vector), 'molecularfeature'])
image<-ggplot(d, aes(x= molecularfeature,y= logpval_vector))+geom_col(color='deepskyblue4', fill='lightblue') + geom_hline(yintercept=-log(0.05/length(pretty_names)),color="red")+ geom_hline(yintercept=-log(0.05),linetype="dashed",color="red")+theme_bw()+theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.title.x=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.y=element_blank()	, legend.position="none") +coord_flip()+labs(y= "Adaptive proto-genes (-Log(P))")
ggsave(paste(general_outdir,"Fig4A_Pvalues_features_AdaptiveProtogenes.pdf",sep=""),width=3.5,height=2.5, useDingbats=FALSE)
ggsave(paste(general_outdir,"FigS3A_Pvalues_features_AdaptiveProtogenes.pdf",sep=""),width=3.5,height=2.5, useDingbats=FALSE)
# effect size - Fig. S3B
d<-stat_df[stat_df$protogene_vector == "proto-gene" & stat_df$fitness_vector =="adaptive",]
d$molecularfeature<-factor(d$molecularfeature, levels=d[order(abs(d$effectsize_vector)), 'molecularfeature'])
image<-ggplot(d, aes(x= molecularfeature,y= effectsize_vector))+geom_col(color='deepskyblue4', fill='lightblue')+theme_bw()+theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.title.x=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.y=element_blank()	, legend.position="none") +coord_flip()+labs(y= "Adaptive proto-genes (effect size)")
ggsave(paste(general_outdir,"FigS3B_EffectSize_features_AdaptiveProtogenes.pdf",sep=""),width=3.5,height=2.5, useDingbats=FALSE)
## deleterious versus neutral
# pvalues - Fig. S3C
d<-stat_df[stat_df$protogene_vector == "proto-gene" & stat_df$fitness_vector =="deleterious",]
d$molecularfeature<-factor(d$molecularfeature, levels=d[order(d$logpval_vector), 'molecularfeature'])
image<-ggplot(d, aes(x= molecularfeature,y= logpval_vector))+geom_col(color='deepskyblue4', fill='lightblue') + geom_hline(yintercept=-log(0.05/length(pretty_names)),color="red")+ geom_hline(yintercept=-log(0.05),linetype="dashed",color="red")+theme_bw()+theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.title.x=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.y=element_blank()	, legend.position="none") +coord_flip()+labs(y= "Deleterious proto-genes (-Log(P))")
ggsave(paste(general_outdir,"FigS3C_Pvalues_features_DeleteriousProtogenes.pdf",sep=""),width=3.5,height=2.5, useDingbats=FALSE)
# effect size - Fig. S3D
d<-stat_df[stat_df$protogene_vector == "proto-gene" & stat_df$fitness_vector =="deleterious",]
d$molecularfeature<-factor(d$molecularfeature, levels=d[order(abs(d$effectsize_vector)), 'molecularfeature'])
image<-ggplot(d, aes(x= molecularfeature,y= effectsize_vector))+geom_col(color='deepskyblue4', fill='lightblue')+theme_bw()+theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.title.x=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.y=element_blank()	, legend.position="none") +coord_flip()+labs(y= "Deleterious proto-genes (effect size)")
ggsave(paste(general_outdir,"FigS3D_EffectSize_features_DeleteriousProtogenes.pdf",sep=""),width=3.5,height=2.5, useDingbats=FALSE)
### genes
## adaptive versus neutral
# pvalues - Fig. S4A
d<-stat_df[stat_df$protogene_vector == "gene" & stat_df$fitness_vector =="adaptive",]
d$molecularfeature<-factor(d$molecularfeature, levels=d[order(d$logpval_vector), 'molecularfeature'])
image<-ggplot(d, aes(x= molecularfeature,y= logpval_vector))+geom_col(color='black', fill='white') + geom_hline(yintercept=-log(0.05/length(pretty_names)),color="red")+ geom_hline(yintercept=-log(0.05),linetype="dashed",color="red")+theme_bw()+theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.title.x=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.y=element_blank()	, legend.position="none") +coord_flip()+labs(y= "Adaptive genes (-Log(P))")
ggsave(paste(general_outdir,"FigS4A_Pvalues_features_AdaptiveGenes.pdf",sep=""),width=3.5,height=2.5, useDingbats=FALSE)
# effect size - Fig. S4B
d<-stat_df[stat_df$protogene_vector == "gene" & stat_df$fitness_vector =="adaptive",]
d$molecularfeature<-factor(d$molecularfeature, levels=d[order(abs(d$effectsize_vector)), 'molecularfeature'])
image<-ggplot(d, aes(x= molecularfeature,y= effectsize_vector))+geom_col(color='black', fill='white')+theme_bw()+theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.title.x=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.y=element_blank()	, legend.position="none") +coord_flip()+labs(y= "Adaptive genes (effect size)")
ggsave(paste(general_outdir,"FigS4B_EffectSize_features_AdaptiveGenes.pdf",sep=""),width=3.5,height=2.5, useDingbats=FALSE)
## deleterious versus neutral
# pvalues - Fig. S4C
d<-stat_df[stat_df$protogene_vector == "gene" & stat_df$fitness_vector =="deleterious",]
d$molecularfeature<-factor(d$molecularfeature, levels=d[order(d$logpval_vector), 'molecularfeature'])
image<-ggplot(d, aes(x= molecularfeature,y= logpval_vector))+geom_col(color='black', fill='white') + geom_hline(yintercept=-log(0.05/length(pretty_names)),color="red")+ geom_hline(yintercept=-log(0.05),linetype="dashed",color="red")+theme_bw()+theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.title.x=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.y=element_blank()	, legend.position="none") +coord_flip()+labs(y= "Deleterious genes (-Log(P))")
ggsave(paste(general_outdir,"FigS4C_Pvalues_features_DeleteriousGenes.pdf",sep=""),width=3.5,height=2.5, useDingbats=FALSE)
# effect size - Fig. S4D
d<-stat_df[stat_df$protogene_vector == "gene" & stat_df$fitness_vector =="deleterious",]
d$molecularfeature<-factor(d$molecularfeature, levels=d[order(abs(d$effectsize_vector)), 'molecularfeature'])
image<-ggplot(d, aes(x= molecularfeature,y= effectsize_vector))+geom_col(color='black', fill='white')+theme_bw()+theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.title.x=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.y=element_blank()	, legend.position="none") +coord_flip()+labs(y= "Deleterious genes (effect size)")
ggsave(paste(general_outdir,"FigS4D_EffectSize_features_DeleteriousGenes.pdf",sep=""),width=3.5,height=2.5, useDingbats=FALSE)

## correlation between TM and aromatocity
cor.test(feature_df$percent_phobius, feature_df$Aromo)

##########################################################################################
# Transmembrane (TM) domains in genes and proto-genes as a function of fitness - Fig. 4B #
###########################################################################################
# TM, evolutionary status and fitness category
df<-orf_table[orf_table$barflex_space=="yes",c("orf_name","protogene", "no_helices_PHOBIUS")]
df$is_phobius<-ifelse(df$no_helices_PHOBIUS>0,1,0)
fitness_df<-merge(fitness_table, orf_table[,c("orf_name","protogene")], all.x=TRUE)
adaptive_orfs<-unique(fitness_df[fitness_df $effect_cs=="increased",]$orf_name)
deleterious_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & fitness_df $effect_cs=="decreased",]$orf_name)
neutral_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & !(fitness_df $orf_name %in% deleterious_orfs),]$orf_name) 
df $fitness_category<-sapply(df $orf_name, function(x) ifelse(x %in% adaptive_orfs,"adaptive",ifelse(x %in% deleterious_orfs, "deleterious", ifelse(x %in% neutral_orfs,"neutral","others"))))
df$plottingclass<-paste(df $fitness_category, df $protogene, sep="_")
df$plottingclass<-factor(df$plottingclass, levels=c("adaptive_proto-gene", "neutral_proto-gene", "deleterious_proto-gene", "adaptive_gene", "neutral_gene","deleterious_gene"))

# calculate proportions with at least 1 TM domain
stat_df<-df[,c("plottingclass","is_phobius")]
phobius_helices<-NULL
phobius_helices$total<- as.vector(table(stat_df $plottingclass))
phobius_helices$yes<-as.vector(sapply(levels(stat_df $plottingclass), function(x) sum(stat_df[stat_df $plottingclass ==x,]$is_phobius)))
phobius_helices$fraction<-phobius_helices$yes/phobius_helices$total
phobius_helices$sder<-sqrt(phobius_helices$fraction * (1-phobius_helices$fraction) / phobius_helices$total)
phobius_helices$class<-levels(stat_df $plottingclass)
phobius_helices<- data.frame(phobius_helices)
phobius_helices $protogene<-c(rep("Proto-genes",3),rep("Genes",3))
phobius_helices $fitness<-c("Adaptive","Neutral","Deleterious","Adaptive","Neutral","Deleterious")
phobius_helices $fitness<-factor(phobius_helices $fitness, levels=c("Adaptive","Neutral","Deleterious"))
phobius_helices $protogene <-factor(phobius_helices $protogene, levels=c("Proto-genes","Genes"))

# plot
image<-ggplot(phobius_helices,aes(x= fitness,y=fraction,ymin=(fraction - sder),ymax=(fraction + sder),fill=fitness, color=fitness))+geom_bar(stat="identity",position=position_dodge())+geom_errorbar(position = position_dodge(0.9),width=0.25)+labs(x="",y="Fraction with TM domain")+ scale_color_manual(values=c("darkgoldenrod1", "darkgray", "darkorchid3")) + scale_fill_manual(values=c("gold", "lightgray", "mediumpurple1"))+facet_grid(.~protogene) +theme_bw()+theme(panel.spacing = unit(0, "lines"), axis.text.y=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.x=element_text(family="Helvetica", face="bold", colour="black", size=7), 	axis.title.y=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.x=element_blank(),	legend.position="none") 
ggsave(paste(general_outdir,"Fig4B_PhobiusTM_fitness_barplots.pdf",sep=""),width=3.5,height=2.5,  useDingbats=FALSE)

# statistical comparison : adaptive versus neutral genes Fisher test
fisher.test(matrix(c(phobius_helices[phobius_helices$class=="adaptive_gene",]$yes, phobius_helices[phobius_helices$class=="neutral_gene",]$yes, phobius_helices[phobius_helices$class=="adaptive_gene",]$total - phobius_helices[phobius_helices$class=="adaptive_gene",]$yes, phobius_helices[phobius_helices$class=="neutral_gene",]$total - phobius_helices[phobius_helices$class=="neutral_gene",]$yes),2,2))

# statistical comparison : adaptive versus neutral proto-genes Fisher test
fisher.test(matrix(c(phobius_helices[phobius_helices$class=="adaptive_proto-gene",]$yes, phobius_helices[phobius_helices$class=="neutral_proto-gene",]$yes, phobius_helices[phobius_helices$class=="adaptive_proto-gene",]$total - phobius_helices[phobius_helices$class=="adaptive_proto-gene",]$yes, phobius_helices[phobius_helices$class=="neutral_proto-gene",]$total - phobius_helices[phobius_helices$class=="neutral_proto-gene",]$yes),2,2))

###########################################################################################
# Distribution of %TM residues in genes and proto-genes as a function of fitness - Fig. S5#
###########################################################################################

###data
df<-orf_table[orf_table$barflex_space=="yes",c("orf_name","protogene", "percent_phobius", "percent_tmhmm")]
fitness_df<-merge(fitness_table, orf_table[,c("orf_name","protogene")], all.x=TRUE)
adaptive_orfs<-unique(fitness_df[fitness_df $effect_cs=="increased",]$orf_name)
deleterious_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & fitness_df $effect_cs=="decreased",]$orf_name)
neutral_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & !(fitness_df $orf_name %in% deleterious_orfs),]$orf_name) 
df $class<-sapply(df $orf_name, function(x) ifelse(x %in% adaptive_orfs,"adaptive",ifelse(x %in% deleterious_orfs, "deleterious", ifelse(x %in% neutral_orfs,"neutral","others"))))

### Plots

# proto-genes using phobius
thisdf<-df[df$protogene=="proto-gene" , c("class","percent_phobius")]
thismu<-data.frame(class=c("adaptive","neutral","deleterious"),Mean = c(mean(thisdf[thisdf$class=="adaptive","percent_phobius"],na.rm=TRUE), mean(thisdf[thisdf$class=="neutral","percent_phobius"],na.rm=TRUE),mean(thisdf[thisdf$class=="deleterious","percent_phobius"],na.rm=TRUE) ))
plot_xlabel <-""
plot_ylabel<-"Proto-genes"
image_protoTMphobius<-ggplot(thisdf, aes(x= percent_phobius,color= class,fill= class))+geom_density(na.rm=TRUE)+
   geom_vline(data=thismu, aes(xintercept=Mean, color= class), linetype="dashed")+scale_color_manual(values=c("darkgoldenrod1", "darkorchid3","darkgray"))+xlim(c(0,0.8))+scale_fill_manual(values=c("gold","NA", "NA"))+labs(x=plot_xlabel,y=plot_ylabel) + theme_bw()+theme(axis.text.x=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.y=element_blank(),axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")+  theme(panel.background = element_rect(fill = "lightcyan2"), panel.grid.major = element_line(colour = "white"),   panel.grid.minor = element_line(colour = "white"))

# proto-genes using tmhmm
thisdf<-df[df$protogene=="proto-gene" , c("class","percent_tmhmm")]
thismu<-data.frame(class=c("adaptive","neutral","deleterious"),Mean = c(mean(thisdf[thisdf$class=="adaptive","percent_tmhmm"],na.rm=TRUE), mean(thisdf[thisdf$class=="neutral","percent_tmhmm"],na.rm=TRUE),mean(thisdf[thisdf$class=="deleterious","percent_tmhmm"],na.rm=TRUE) ))
plot_xlabel <-""
plot_ylabel<-""
image_protoTMtmhmm<-ggplot(thisdf, aes(x= percent_tmhmm,color= class,fill= class))+geom_density(na.rm=TRUE)+
   geom_vline(data=thismu, aes(xintercept=Mean, color= class), linetype="dashed")+scale_color_manual(values=c("darkgoldenrod1", "darkorchid3","darkgray"))+xlim(c(0,0.8))+scale_fill_manual(values=c("gold","NA", "NA"))+labs(x=plot_xlabel,y=plot_ylabel) + theme_bw()+theme(axis.text.x=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.y=element_blank(),axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")+  theme(panel.background = element_rect(fill = "lightcyan2"), panel.grid.major = element_line(colour = "white"),   panel.grid.minor = element_line(colour = "white"))

# genes using phobius
thisdf<-df[df$protogene=="gene" , c("class","percent_phobius")]
thismu<-data.frame(class=c("adaptive","neutral","deleterious"),Mean = c(mean(thisdf[thisdf$class=="adaptive","percent_phobius"],na.rm=TRUE), mean(thisdf[thisdf$class=="neutral","percent_phobius"],na.rm=TRUE),mean(thisdf[thisdf$class=="deleterious","percent_phobius"],na.rm=TRUE) ))
plot_xlabel <-"TM residues content (Phobius)"
plot_ylabel<-"Genes"
image_genesTMphobius<-ggplot(thisdf, aes(x= percent_phobius,color= class,fill= class))+geom_density(na.rm=TRUE)+
   geom_vline(data=thismu, aes(xintercept=Mean, color= class), linetype="dashed")+scale_color_manual(values=c("darkgoldenrod1", "darkorchid3","darkgray"))+xlim(c(0,0.8))+scale_fill_manual(values=c("gold","NA", "NA"))+labs(x=plot_xlabel,y=plot_ylabel) + theme_bw()+theme(axis.text.x=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.y=element_blank(),axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")

# genes using tmhmm
thisdf<-df[df$protogene=="gene" , c("class","percent_tmhmm")]
thismu<-data.frame(class=c("adaptive","neutral","deleterious"),Mean = c(mean(thisdf[thisdf$class=="adaptive","percent_tmhmm"],na.rm=TRUE), mean(thisdf[thisdf$class=="neutral","percent_tmhmm"],na.rm=TRUE),mean(thisdf[thisdf$class=="deleterious","percent_tmhmm"],na.rm=TRUE) ))
plot_xlabel <-"TM residues content (TMHMM)"
plot_ylabel<-""
image_genesTMtmhmm<-ggplot(thisdf, aes(x= percent_tmhmm,color= class,fill= class))+geom_density(na.rm=TRUE)+
   geom_vline(data=thismu, aes(xintercept=Mean, color= class), linetype="dashed")+scale_color_manual(values=c("darkgoldenrod1", "darkorchid3","darkgray"))+xlim(c(0,0.8))+scale_fill_manual(values=c("gold","NA", "NA"))+labs(x=plot_xlabel,y=plot_ylabel) + theme_bw()+theme(axis.text.x=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.y=element_blank(),axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")

pdf(paste(general_outdir,"FigS5_TMResiduesContent.pdf",sep=""),width=6,height=5, useDingbats=FALSE)
multiplot(image_protoTMphobius, image_protoTMtmhmm, image_genesTMphobius, image_genesTMtmhmm, cols=2)
dev.off()

######################################################################################################
# Fraction of real and scrambled sequences with putative TM domain - WITH Length Control - Figure 5A #
######################################################################################################
# paramaters
set.seed(35689118)
scramble_color<-"wheat3"

# data
df<-tm_table
df$has_tm<-ifelse(df$no_hel_Phob>0,1,0)

###length controlled data (all proto-genes, and other types are picked to follow proto-gene length distribution)
controlled_df<-df[df$type=="proto-genes",]
# determine the probability distribution for proto-genes
h<-hist(df[df$type=="proto-genes",]$length, breaks= seq(0,650,by=25), plot=FALSE)
protobreaks<-h$breaks
protoprob<-h$counts/nrow(df[df$type=="proto-genes",]) 
# Select 1000 non-genic ORFs following the proto-gene length distribution, with replacement
t<-df[df $type=="non-genic ORFs" & df $mode=="real",]
adj_protobreaks<-c(seq(25,150,by=25), 650)
adj_protoprob<-c(sum(protoprob[1:2]),protoprob[3:6], sum(protoprob[7:length(protoprob)]))
t$lengthbin<-cut(t$length, breaks=adj_protobreaks)
bindistrib<-as.data.frame(table(sample(x=levels(t$lengthbin), 1000, replace=T, prob= adj_protoprob)))
l<-sapply(bindistrib$Var1, function(x) sample(t[t$lengthbin==x,]$names, bindistrib[bindistrib$Var1==x,]$Freq, replace=T))
rownumbers<-NULL
for(x in unlist(l)){rownumbers<-c(rownumbers, which(df$names==x,))}
controlled_df <-rbind(controlled_df, df[rownumbers,])
# Select 1000 genes following the proto-gene length distribution, with replacement
t<-df[df $type=="genes" & df $mode=="real",]
adj_protobreaks<-c(seq(25,300,by=25), 650)
adj_protoprob<-c(sum(protoprob[1:2]), protoprob[3:12], sum(protoprob[13:length(protoprob)]))
t$lengthbin<-cut(t$length, breaks=adj_protobreaks)
t<-t[!(is.na(t$lengthbin)),]
bindistrib<-as.data.frame(table(sample(x=levels(t$lengthbin), 1000, replace=T, prob= adj_protoprob)))
l<-sapply(bindistrib$Var1, function(x) sample(t[t$lengthbin==x,]$names, bindistrib[bindistrib$Var1==x,]$Freq, replace=T))
rownumbers<-NULL
for(x in unlist(l)){rownumbers<-c(rownumbers, which(df$names==x,))}
controlled_df <-rbind(controlled_df, df[rownumbers,])
# Select 1000 intergenes following the proto-gene length distribution, with replacement
t<-df[df $type=="intergenes" & df $mode=="real",]
adj_protobreaks<-c(seq(25,300,by=25 ),650)
adj_protoprob<-c(sum(protoprob[1:2]), protoprob[3:12], sum(protoprob[13:length(protoprob)]))
t$lengthbin<-cut(t$length, breaks=adj_protobreaks)
t<-t[!(is.na(t$lengthbin)),]
bindistrib<-as.data.frame(table(sample(x=levels(t$lengthbin), 1000, replace=T, prob= adj_protoprob)))
l<-sapply(bindistrib$Var1, function(x) sample(t[t$lengthbin==x,]$names, bindistrib[bindistrib$Var1==x,]$Freq, replace=T))
rownumbers<-NULL
for(x in unlist(l)){rownumbers<-c(rownumbers, which(df$names==x,))}
controlled_df <-rbind(controlled_df, df[rownumbers,])

# calculate proportions with at least 1 TM domain in length controlled data
statistics_table<-NULL
statistics_table $total<- as.vector(table(controlled_df $type, controlled_df $mode))
statistics_table $yes<-as.vector(sapply(levels(controlled_df $mode), function(x) sapply(levels(controlled_df $type), function(y) sum(controlled_df[controlled_df $type ==y & controlled_df $mode==x,]$has_tm))))
statistics_table $fraction<-statistics_table $yes/statistics_table $total
statistics_table $sder<-sqrt(statistics_table $fraction * (1-statistics_table $fraction) / statistics_table $total)
statistics_table $evoclass<-rep(levels(controlled_df $type), 2)
statistics_table $seqclass<-c(rep(levels(controlled_df $mode)[1], 4),rep(levels(controlled_df $mode)[2], 4))
statistics_table <- data.frame(statistics_table)
statistics_table $evoclass <-factor(statistics_table $evoclass, levels=c("genes","proto-genes","intergenes","non-genic ORFs"))

# plot
image<-ggplot(statistics_table,aes(x= seqclass,y=fraction,ymin=(fraction - sder),ymax=(fraction + sder),fill=seqclass))+geom_bar(stat="identity",position=position_dodge(), color="black")+geom_errorbar(position = position_dodge(0.9),width=0.25)+labs(x="",y="Fraction with TM domain")+ ylim(0,0.7)+ scale_fill_manual(values=c("pink", scramble_color))  + facet_grid(.~evoclass)+theme_bw()+theme(panel.spacing = unit(0.2, "lines"), axis.text=element_text(family="Helvetica", face="bold", colour="black", size=8),  	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")
ggsave(paste(general_outdir,"Fig5A_FractionWithTM_Lengthcontrol.pdf",sep=""),width=5,height=2.5,  useDingbats=FALSE)

##################################################################################################################
# TM propensity in real and scrambled sequences as a function of Thymine content - WITH Length Control - Figure 5B #
##################################################################################################################
# data
df<-tm_table
df$has_tm<-ifelse(df$no_hel_Phob>0,1,0)

###length controlled data (all proto-genes, and other types are picked to follow proto-gene length distribution)
controlled_df<-df[df$type=="proto-genes",]
# determine the probability distribution for proto-genes
h<-hist(df[df$type=="proto-genes",]$length, breaks= seq(0,650,by=25), plot=FALSE)
protobreaks<-h$breaks
protoprob<-h$counts/nrow(df[df$type=="proto-genes",]) 
# Select 1000 non-genic ORFs following the proto-gene length distribution, with replacement
t<-df[df $type=="non-genic ORFs" & df $mode=="real",]
adj_protobreaks<-c(seq(25,150,by=25), 650)
adj_protoprob<-c(sum(protoprob[1:2]),protoprob[3:6], sum(protoprob[7:length(protoprob)]))
t$lengthbin<-cut(t$length, breaks=adj_protobreaks)
bindistrib<-as.data.frame(table(sample(x=levels(t$lengthbin), 1000, replace=T, prob= adj_protoprob)))
l<-sapply(bindistrib$Var1, function(x) sample(t[t$lengthbin==x,]$names, bindistrib[bindistrib$Var1==x,]$Freq, replace=T))
rownumbers<-NULL
for(x in unlist(l)){rownumbers<-c(rownumbers, which(df$names==x,))}
controlled_df <-rbind(controlled_df, df[rownumbers,])
# Select 1000 genes following the proto-gene length distribution, with replacement
t<-df[df $type=="genes" & df $mode=="real",]
adj_protobreaks<-c(seq(25,300,by=25), 650)
adj_protoprob<-c(sum(protoprob[1:2]), protoprob[3:12], sum(protoprob[13:length(protoprob)]))
t$lengthbin<-cut(t$length, breaks=adj_protobreaks)
t<-t[!(is.na(t$lengthbin)),]
bindistrib<-as.data.frame(table(sample(x=levels(t$lengthbin), 1000, replace=T, prob= adj_protoprob)))
l<-sapply(bindistrib$Var1, function(x) sample(t[t$lengthbin==x,]$names, bindistrib[bindistrib$Var1==x,]$Freq, replace=T))
rownumbers<-NULL
for(x in unlist(l)){rownumbers<-c(rownumbers, which(df$names==x,))}
controlled_df <-rbind(controlled_df, df[rownumbers,])
# Select 1000 intergenes following the proto-gene length distribution, with replacement
t<-df[df $type=="intergenes" & df $mode=="real",]
adj_protobreaks<-c(seq(25,300,by=25 ),650)
adj_protoprob<-c(sum(protoprob[1:2]), protoprob[3:12], sum(protoprob[13:length(protoprob)]))
t$lengthbin<-cut(t$length, breaks=adj_protobreaks)
t<-t[!(is.na(t$lengthbin)),]
bindistrib<-as.data.frame(table(sample(x=levels(t$lengthbin), 1000, replace=T, prob= adj_protoprob)))
l<-sapply(bindistrib$Var1, function(x) sample(t[t$lengthbin==x,]$names, bindistrib[bindistrib$Var1==x,]$Freq, replace=T))
rownumbers<-NULL
for(x in unlist(l)){rownumbers<-c(rownumbers, which(df$names==x,))}
controlled_df <-rbind(controlled_df, df[rownumbers,])
controlled_df $type <-factor(controlled_df $type, levels=c("genes","proto-genes","intergenes","non-genic ORFs"))

# plot fraction of sequences with tm domain and density as a function of T
binsize<-0.1
allbins<-seq(0.05+binsize/2,0.75-binsize/2,binsize)
image_bars<-ggplot(controlled_df,aes(x= T, fill=mode)) +stat_summary_bin(aes(y = has_tm), fun.y="mean", bins=length(allbins),geom="bar", position=position_dodge(), color="darkgray")+stat_summary_bin(aes(y = has_tm), fun.data="mean_se", bins=length(allbins),geom="errorbar", position=position_dodge(),color="gray") + xlim(0.05, 0.7)+scale_fill_manual(values=c("pink", scramble_color))+labs(x="",y="With TM domain")  + facet_wrap(~type,nrow=1)+theme_bw()+theme(panel.spacing = unit(0.2, "lines"), axis.text.y=element_text(family="Helvetica", face="bold", colour="black", size=8), axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none", strip.text=element_blank(), axis.text.x=element_blank()) + geom_density(aes(y=..scaled..), color="darkgrey", fill="lightgrey", alpha = 0.1, size=0.4)

# plot % Tm residues content as a function of T, considering only those sequences that do have a putative tm domain
controlled_df$fractiontm<-controlled_df$no_res_TM_Phob / controlled_df$length
image_points<-ggplot(controlled_df[controlled_df$has_tm==1,],aes(x= T,y= fractiontm,colour=mode, fill=mode))+ geom_point(size=0.2, alpha=0.1) + geom_smooth(method="lm", na.rm=TRUE) + xlim(0.05, 0.7) + scale_color_manual(values=c("pink", scramble_color)) +labs(x="Thymine content",y="TM residues content")+ scale_fill_manual(values=c("pink", scramble_color))  + ylim(0,1)+facet_wrap(~type,nrow=1, strip.position="right")+theme_bw()+theme(panel.spacing = unit(0.2, "lines"), axis.text=element_text(family="Helvetica", face="bold", colour="black", size=8), axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none",strip.text=element_blank())

# composite plot
p_double<-ggarrange(image_bars,image_points,nrow=2)
ggsave(paste(general_outdir,"Fig5B_Relationship_TM_Thymine_Lengthcontrol.pdf",sep=""),width=5,height=3.5,  useDingbats=FALSE)

#####################################################################################################################################
# TM propensity in non-genic ORFs (25 - 75 codons) as a function of their length relative to their containing intergene - Figure 5C #
#####################################################################################################################################
# data
df<-ngo_intergene_table
otherdf<-tm_table

# TM propensity to expect for proto-genes and intergenes of same length:
protovalue<-100*nrow(otherdf[otherdf $type=="proto-genes" & otherdf $mode=="real" & otherdf $no_hel_Phob>0 & otherdf $length>25 & otherdf $length<75,])/nrow(otherdf[otherdf $type=="proto-genes" & otherdf $mode=="real" & otherdf $length>25 & otherdf $length<75,])
intergenevalue<-100*nrow(otherdf[otherdf $type=="intergenes" & otherdf $mode=="real" & otherdf $no_hel_Phob >0 & otherdf $length>25 & otherdf $length<75,])/nrow(otherdf[otherdf $type=="intergenes" & otherdf $mode=="real" & otherdf $length>25 & otherdf $length<75,])

# summarize
t<-NULL
t$bin<-sort(unique(df$bin))
t$total<-as.vector(table(df$bin))
t$yes<-sapply(sort(unique(df$bin)), function(x) sum(ifelse(df[df$bin ==x,]$NGO_no_hel>0,1,0)))
t$fraction<-t$yes/t$total
t$sder<-sqrt(t$fraction * (1-t$fraction) / t$total)
t<-data.frame(t)

# plot
image<-ggplot(t,aes(x= bin, y=fraction *100, ymin=100*(fraction - sder),ymax=100*(fraction + sder)))+ geom_bar(stat="identity", colour="black", fill="pink") +geom_errorbar(width=0.025)+labs(x="% of intergene's length covered by non-genic ORF",y="% of non-genic ORFs \n with TM domain")+geom_hline(yintercept= protovalue,linetype="dashed",color="lightblue", size=1)+ geom_hline(yintercept= intergenevalue,linetype="dashed",color="darkblue", size=1)+theme_bw()+theme( axis.text=element_text(family="Helvetica", face="bold", colour="black", size=8),  	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")
ggsave(paste(general_outdir,"Fig5C_FractionWithTM_NGO_PercentIntergene.pdf",sep=""),width=5,height=2,  useDingbats=FALSE)


###################################################
# Distribution of length and T content - Figure S6#
###################################################
# data
df<-tm_table
# length distribution (log)
image<-ggplot(df[df$mode=="real",], aes(x=log(length), color=type, linetype=type))+geom_density()+scale_color_manual(values=c("black","lightblue","blue","darkblue"))+ scale_linetype_manual(values=c("solid","solid","longdash","dotted"))+theme_bw()+theme( axis.text=element_text(family="Helvetica", face="bold", colour="black", size=8),  	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10))
ggsave(paste(general_outdir,"FigS6A_Density_loglength.pdf",sep=""),width=5,height=2.5,  useDingbats=FALSE)
# T distribution
image<-ggplot(df[df$mode=="real",], aes(x=T, color=type, linetype=type))+geom_density()+scale_color_manual(values=c("black","lightblue","darkblue","darkblue"))+ scale_linetype_manual(values=c("solid","solid","longdash","dotted"))+theme_bw()+theme( axis.text=element_text(family="Helvetica", face="bold", colour="black", size=8),  	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10))
ggsave(paste(general_outdir,"FigS6B_Density_Thymine.pdf",sep=""),width=5,height=2.5,  useDingbats=FALSE)


####################################################################################################
# Fraction of real and scrambled sequences with putative TM domain - NO Length Control - Figure S7A #
####################################################################################################
# data
df<-tm_table
df$has_tm<-ifelse(df$no_hel_Phob>0,1,0)

# calculate proportions with at least 1 TM domain in length controlled data
statistics_table<-NULL
statistics_table $total<- as.vector(table(df $type, df $mode))
statistics_table $yes<-as.vector(sapply(levels(df $mode), function(x) sapply(levels(df $type), function(y) sum(df[df $type ==y & df $mode==x,]$has_tm))))
statistics_table $fraction<-statistics_table $yes/statistics_table $total
statistics_table $sder<-sqrt(statistics_table $fraction * (1-statistics_table $fraction) / statistics_table $total)
statistics_table $evoclass<-rep(levels(df $type), 2)
statistics_table $seqclass<-c(rep(levels(df $mode)[1], 4),rep(levels(df $mode)[2], 4))
statistics_table <- data.frame(statistics_table)
statistics_table $evoclass <-factor(statistics_table $evoclass, levels=c("genes","proto-genes","intergenes","non-genic ORFs"))

# plot
image<-ggplot(statistics_table,aes(x= seqclass,y=fraction,ymin=(fraction - sder),ymax=(fraction + sder),fill=seqclass))+geom_bar(stat="identity",position=position_dodge(), color="black")+geom_errorbar(position = position_dodge(0.9),width=0.25)+labs(x="",y="Fraction with TM domain")+ ylim(0,0.7)+ scale_fill_manual(values=c("pink", scramble_color))  + facet_grid(.~evoclass)+theme_bw()+theme(panel.spacing = unit(0.2, "lines"), axis.text=element_text(family="Helvetica", face="bold", colour="black", size=8),  	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")
ggsave(paste(general_outdir,"FigS7A_FractionWithTM_NoLengthcontrol.pdf",sep=""),width=5,height=2.5,  useDingbats=FALSE)

######################################################################################################################
# TM propensity in real and scrambled sequences as a function of Thymine content - NO Length Control - Figure S7B #
######################################################################################################################
# data
df<-tm_table
df$has_tm<-ifelse(df$no_hel_Phob>0,1,0)
df $type <-factor(df $type, levels=c("genes","proto-genes","intergenes","non-genic ORFs"))

# plot fraction of sequences with tm domain and density as a function of T
binsize<-0.1
allbins<-seq(0.05+binsize/2,0.75-binsize/2,binsize)
image_bars<-ggplot(df,aes(x= T, fill=mode)) +stat_summary_bin(aes(y = has_tm), fun.y="mean", bins=length(allbins),geom="bar", position=position_dodge(), color="darkgray")+stat_summary_bin(aes(y = has_tm), fun.data="mean_se", bins=length(allbins),geom="errorbar", position=position_dodge(),color="gray") + xlim(0.05, 0.8)+scale_fill_manual(values=c("pink", scramble_color))+labs(x="",y="With TM domain")  + facet_wrap(~type,nrow=1)+theme_bw()+theme(panel.spacing = unit(0.2, "lines"), axis.text.y=element_text(family="Helvetica", face="bold", colour="black", size=8), axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none", strip.text=element_blank(), axis.text.x=element_blank()) + geom_density(aes(y=..scaled..), color="darkgrey", fill="lightgrey", alpha = 0.1, size=0.4)

# plot % Tm residues content as a function of T, considering only those sequences that do have a putative tm domain
df$fractiontm<-df$no_res_TM_Phob / df$length
image_points<-ggplot(df[df$has_tm==1,],aes(x= T,y= fractiontm,colour=mode, fill=mode))+ geom_point(size=0.2, alpha=0.1) + geom_smooth(method="lm", na.rm=TRUE) + xlim(0.05, 0.8) + scale_color_manual(values=c("pink", scramble_color)) +labs(x="Thymine content",y="TM residues content")+ scale_fill_manual(values=c("pink", scramble_color))  + ylim(0,1)+facet_wrap(~type,nrow=1, strip.position="right")+theme_bw()+theme(panel.spacing = unit(0.2, "lines"), axis.text=element_text(family="Helvetica", face="bold", colour="black", size=8), axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none",strip.text=element_blank())

# composite plot
p_double<-ggarrange(image_bars,image_points,nrow=2)
ggsave(paste(general_outdir,"FigS7B_Relationship_TM_Thymine_NOLengthcontrol.pdf",sep=""),width=5,height=3.5,  useDingbats=FALSE)

