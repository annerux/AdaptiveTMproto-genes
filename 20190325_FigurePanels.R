# March 25 2019
# Anne-Ruxandra Carvunis

# this script reproduces all the figure panels resulting from computational analyses in Vakirlis, Hsu et al from the supplementary Data tables. 


# LIBRARIES
library(ggplot2)
library(plyr)
library(ggpubr)

# paramaters
set.seed(35689118)
scramble_color<-"wheat3" 
  
# FOLDER
general_outdir<-"/users/annerux/Workspace/cerevisiae/txt/screen_analyses/Protogenes/ProtogeneFitnessFigurePanels/submission/" # CHANGE THIS LINE TO INDICATE WHERE YOU STORED THE DATA FILES 1-5


# DATA SOURCES

# description of ORFs used in analyses 
orf_table<-read.csv(paste(general_outdir, "Data1.csv",sep=""),header=TRUE)
orf_table $protogene <-factor(orf_table $protogene, levels=c("Emerging ORFs","Established ORFs"))
orf_table $overexpression_relative_fitness <-factor(orf_table $overexpression_relative_fitness, levels=c("increased","decreased","unchanged"))

# individual normalized colony sizes for 1 environmental condition (Fig. S3B)
colony_table<-read.csv(paste(general_outdir, "Data2.csv",sep=""),header=TRUE)

# results of overexpression screen in 5 environmental conditions 
fitness_table<-read.csv(paste(general_outdir, "Data3.csv",sep=""),header=TRUE)

# transmembrane analyses (related to Figs 5 and S7)
tm_table<-read.csv(paste(general_outdir, "Data4.csv",sep=""),header=TRUE)
tm_labels<-as.vector(tm_table$type)
tm_labels[tm_table$names %in% orf_table[orf_table$protogene=="Emerging ORFs",]$orf_name ]<-"emerging ORFs"
tm_labels[tm_table$names %in% orf_table[orf_table$protogene=="Established ORFs",]$orf_name ]<-"established ORFs"
tm_table $type<-tm_labels
tm_table $type <-factor(tm_table $type, levels=c("established ORFs","emerging ORFs","artificial ORFs","sORFs"))
ngo_intergene_table<-read.csv(paste(general_outdir, "Data5.csv",sep=""),header=TRUE)


#####################################################################
# Figure 1B: Counting ORFs with deletion and overexpression strains #
######################################################################
# total annotated ORFs:
table(orf_table$protogene)
# with deletion strain
df<-orf_table[!is.na(orf_table$loss_fitness) & orf_table$overlap==0,c("protogene","loss_fitness")]
table(df[,1])
# with overexpression strain
table(orf_table[orf_table$barflex_space=="yes",c("protogene")])

#######################################
# Figure 2A: fitness cost of ORF loss #
#######################################

# data
df<-orf_table[!is.na(orf_table$loss_fitness) & orf_table$overlap==0,c("protogene","loss_fitness")]

# statistical comparison: Mann Whitney
wilcoxtest<-wilcox.test(df[df $protogene=="Established ORFs","loss_fitness"],df[df $protogene=="Emerging ORFs","loss_fitness"])
pval<-wilcoxtest$p.val

# statistical comparison: Fisher's exact with fitness cutoff of 0.9 
n_pg_effect<-nrow(df[df$protogene=="Emerging ORFs" & df$loss_fitness <0.9,]) 
n_pg_noeffect<-nrow(df[df$protogene=="Emerging ORFs" & df$loss_fitness>=0.9,])
n_g_effect<-nrow(df[df$protogene=="Established ORFs" & df$loss_fitness <0.9,]) 
n_g_noeffect<-nrow(df[df$protogene=="Established ORFs" & df$loss_fitness>=0.9,])
pval<-fisher.test(matrix(c(n_pg_effect, n_pg_noeffect, n_g_effect, n_g_noeffect),2,2))

# cumulative plot
plot_xlabel<-"Fitness upon loss"
plot_ylabel<- "ORFs (cumulative frequency)"
image<-ggplot(df, aes(x= loss_fitness,color=protogene, fill=protogene))+stat_ecdf(geom="step", na.rm=TRUE)+ stat_ecdf(aes(ymin=0,ymax=..y..),
geom="ribbon", na.rm=TRUE)+scale_color_manual(values=c("deepskyblue4","black"))+ scale_fill_manual(values=c("lightblue","NA"))+labs(x=plot_xlabel,y=plot_ylabel)+ geom_vline(xintercept=0.9, colour="red") + theme_bw()+theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8),  	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none") 
ggsave(paste(general_outdir,"Fig2A_Compare_FitnessUponLoss_CumulativeFreq.pdf",sep=""),width=2.5,height=2.5)

#########################################################
# Figure 2B: fixation of ORf sctructure across isolates #
#########################################################

df<-orf_table[!(is.na(orf_table $strain_goodalign)),c("protogene","strain_goodalign", "strained_conserved_orf",  "strain_pi")]

# plot
plot_xlabel<-"% isolates with intact ORF"
plot_ylabel<- "ORFs (cumulative frequency)"
image<-ggplot(df, aes(x= strained_conserved_orf/2022*100,color=protogene, fill=protogene))+stat_ecdf(geom="step", na.rm=TRUE)+ stat_ecdf(aes(ymin=0,ymax=..y..),
geom="ribbon", na.rm=TRUE)+scale_color_manual(values=c("deepskyblue4","black"))+ scale_fill_manual(values=c("lightblue","NA"))+labs(x=plot_xlabel,y=plot_ylabel) + geom_vline(xintercept=90, colour="red")+ theme_bw()+theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8),  	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none") + coord_cartesian(xlim=c(0,100), , expand=FALSE)
ggsave(paste(general_outdir,"Fig2B_Compare_StrainORFStructure_CumulativeFreq.pdf",sep=""),width=2.5,height=2.5)

# intact orf structure in <90% isolates
df$strict<-ifelse(df$strained_conserved_orf/2022*100<90,1,0)
table(df[,c("protogene","strict")])

###################################################
# Figure 2C: nucleotide diversity across isolates #
###################################################

df<-orf_table[!(is.na(orf_table $strain_goodalign)),c("protogene","strain_goodalign", "strained_conserved_orf",  "strain_pi")]

# statistical comaprison: Mann Whitney
wilcoxtest<-wilcox.test(df[df $protogene=="Established ORFs","strain_pi"],df[df $protogene=="Emerging ORFs","strain_pi"])
pval<-wilcoxtest$p.val

# plot
mu <- ddply(df, "protogene", summarise, grp.mean=mean(strain_pi,na.rm=TRUE))
plot_xlabel<-"Nucleotide diversity"
plot_ylabel<- "ORFs (density)"

image<-ggplot(df, aes(x= strain_pi,color=protogene,fill=protogene))+geom_density(na.rm=TRUE)+
   geom_vline(data=mu, aes(xintercept=grp.mean, color=protogene), linetype="dashed")+scale_color_manual(values=c("deepskyblue4","black"))+scale_fill_manual(values=c( "lightblue","NA"))+labs(x=plot_xlabel,y=plot_ylabel)+ theme_bw()+theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8),  	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")
ggsave(paste(general_outdir,"Fig2C_Compare_StrainNucleotideDiversity.pdf",sep=""),width=2.5,height=2.5)

################################################################################
# Figure S1A # fitness upon ORF loss - Control for length and expression level #
################################################################################
# data
df<-orf_table[!is.na(orf_table$loss_fitness) & orf_table$overlap==0,c("orf_name","protogene","loss_fitness", "expr_level", "Length")]

### controlled data (all emerging ORFs, all established ORFs and randomly picked estabslished ORFs that follow the length and expr_level distribution of emerging ORFs)
controlled_df<-df
n_toselect<-nrow(df[df $protogene=="Established ORFs" ,]) # number to pick randomly for controlling

#LENGTH
# determine the probability distribution for proto-genes
h<-hist(df[df$protogene=="Emerging ORFs",]$Length, breaks= seq(50,2000,by=50), plot=FALSE)
protobreaks<-h$breaks
protoprob<-h$counts/nrow(df[df$protogene=="Emerging ORFs",]) 
# Select random genes following the proto-gene length distribution, with replacement
t<-df[df $protogene=="Established ORFs" ,]
adj_protobreaks<-c(seq(50,800,by=50), 2000)
adj_protoprob<-c(protoprob[1:15], sum(protoprob[16:length(protoprob)]))
t$lengthbin<-cut(t$Length, breaks=adj_protobreaks)
t<-t[!(is.na(t$lengthbin)),]
bindistrib<-as.data.frame(table(sample(x=levels(t$lengthbin), n_toselect, replace=T, prob= adj_protoprob)))
l<-sapply(bindistrib$Var1, function(x) sample(t[t$lengthbin==x,]$orf_name, bindistrib[bindistrib$Var1==x,]$Freq, replace=T))
rownumbers<-NULL
for(x in unlist(l)){rownumbers<-c(rownumbers, which(df$orf_name==x,))}
lengthcontrolled_df<-df[rownumbers,]
lengthcontrolled_df$protogene <-rep("esta_lengthcontrolled", n_toselect)
controlled_df <-rbind(controlled_df, lengthcontrolled_df)

# EXPR_LEVEL
h<-hist(df[df$protogene=="Emerging ORFs",]$expr_level, breaks= seq(-7,6,by=1), plot=FALSE)
protobreaks<-h$breaks
protoprob<-h$counts/nrow(df[df$protogene=="Emerging ORFs",]) 
# Select random genes following the proto-gene length distribution, with replacement
t<-df[df $protogene=="Established ORFs" ,]
adj_protobreaks<-c((-7), seq(-5,3,by=1), 6)
adj_protoprob<-c(sum(protoprob[1:2]), protoprob[3:10], sum(protoprob[11:length(protoprob)]))
t$expbin<-cut(t$expr_level, breaks=adj_protobreaks)
t<-t[!(is.na(t$expbin)),]
bindistrib<-as.data.frame(table(sample(x=levels(t$expbin), n_toselect, replace=T, prob= adj_protoprob)))
l<-sapply(bindistrib$Var1, function(x) sample(t[t$expbin==x,]$orf_name, bindistrib[bindistrib$Var1==x,]$Freq, replace=T))
rownumbers<-NULL
for(x in unlist(l)){rownumbers<-c(rownumbers, which(df$orf_name==x,))}
expcontrolled_df<-df[rownumbers,]
expcontrolled_df$protogene <-rep("esta_expcontrolled", n_toselect)
controlled_df <-rbind(controlled_df, expcontrolled_df)

# plot fractions with fitness costs
cutoff<-0.9
stat_df<-controlled_df[,c("protogene","loss_fitness")]
stat_df$is_delet<-ifelse(stat_df$loss_fitness<cutoff,1,0)
plotdf<-NULL
plotdf $total<- as.vector(table(stat_df $protogene))
plotdf $yes<-as.vector(sapply(levels(stat_df $protogene), function(x) sum(stat_df[stat_df $protogene ==x,]$is_delet)))
plotdf $fraction<-plotdf $yes/plotdf $total
plotdf $sder<-sqrt(plotdf $fraction * (1-plotdf $fraction) / plotdf $total)
plotdf $class<-levels(stat_df $protogene)
plotdf <- data.frame(plotdf)
plotdf$class<-factor(plotdf$class, levels = c("Emerging ORFs","esta_expcontrolled", "esta_lengthcontrolled", "Established ORFs"))
levels(plotdf$class)<-mapvalues(levels(plotdf$class),from=levels(plotdf$class), to= c("Emerging ORFs", "Established ORFs/matched expression", "Established ORFs/matched length", "All established ORFs"))

image<-ggplot(plotdf,aes(x= class,y=fraction,ymin=(fraction - sder),ymax=(fraction + sder),fill=class, color=class))+geom_bar(stat="identity",position=position_dodge())+geom_errorbar(position = position_dodge(0.9),width=0.25)+coord_flip()+ labs(x="",y="Fraction with Fitness Cost")+ scale_color_manual(values=c("deepskyblue4", "black", "black","black")) + scale_fill_manual(values=c("lightblue", "white", "white","white")) +theme_bw()+theme(panel.spacing = unit(0, "lines"), axis.text.x=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.y=element_text(family="Helvetica", face="bold", colour="black", size=7), 	axis.title.x=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.y=element_blank(),	legend.position="none") 
ggsave(paste(general_outdir,"FigS1A_Compare_FitnessUponLoss_controlled_barplots.pdf",sep=""),width=4.5,height=2.5,  useDingbats=FALSE)

# statistical comparison: Fisher's exact with fitness cutoff ; between emerging ORFs and expr_controled
n_pg_effect<-nrow(controlled_df[controlled_df $protogene=="Emerging ORFs" & controlled_df $loss_fitness <cutoff,]) 
n_pg_noeffect<-nrow(controlled_df[controlled_df $protogene=="Emerging ORFs" & controlled_df $loss_fitness>=cutoff,])
n_g_effect<-nrow(controlled_df[controlled_df $protogene=="esta_expcontrolled" & controlled_df $loss_fitness <cutoff,]) 
n_g_noeffect<-nrow(controlled_df[controlled_df $protogene=="esta_expcontrolled" & controlled_df $loss_fitness>=cutoff,])
pval<-fisher.test(matrix(c(n_pg_effect, n_pg_noeffect, n_g_effect, n_g_noeffect),2,2))

# statistical comparison: Fisher's exact with fitness cutoff ; between emerging ORFs and length_controled
n_pg_effect<-nrow(controlled_df[controlled_df $protogene=="Emerging ORFs" & controlled_df $loss_fitness <cutoff,]) 
n_pg_noeffect<-nrow(controlled_df[controlled_df $protogene=="Emerging ORFs" & controlled_df $loss_fitness>=cutoff,])
n_g_effect<-nrow(controlled_df[controlled_df $protogene=="esta_lengthcontrolled" & controlled_df $loss_fitness <cutoff,]) 
n_g_noeffect<-nrow(controlled_df[controlled_df $protogene=="esta_lengthcontrolled" & controlled_df $loss_fitness>=cutoff,])
pval<-fisher.test(matrix(c(n_pg_effect, n_pg_noeffect, n_g_effect, n_g_noeffect),2,2))


################################################################################
# Figure S1B Fixation of orf structure - Control for length and expression level #
################################################################################
# data
df<-orf_table[!(is.na(orf_table $strain_goodalign)),c("orf_name","protogene","strain_goodalign", "strained_conserved_orf", "expr_level", "Length")]

### controlled data (all emerging ORFs, all established ORFs and randomly picked estabslished ORFs that follow the length and expr_level distribution of emerging ORFs)
controlled_df<-df
n_toselect<-nrow(df[df $protogene=="Established ORFs" ,])

#LENGTH
# determine the probability distribution for proto-genes
h<-hist(df[df$protogene=="Emerging ORFs",]$Length, breaks= seq(50,2000,by=50), plot=FALSE)
protobreaks<-h$breaks
protoprob<-h$counts/nrow(df[df$protogene=="Emerging ORFs",]) 
# Select genes following the proto-gene length distribution, with replacement
t<-df[df $protogene=="Established ORFs" ,]
adj_protobreaks<-c(seq(50,800,by=50), 2000)
adj_protoprob<-c(protoprob[1:15], sum(protoprob[16:length(protoprob)]))
t$lengthbin<-cut(t$Length, breaks=adj_protobreaks)
t<-t[!(is.na(t$lengthbin)),]
bindistrib<-as.data.frame(table(sample(x=levels(t$lengthbin), n_toselect, replace=T, prob= adj_protoprob)))
l<-sapply(bindistrib$Var1, function(x) sample(t[t$lengthbin==x,]$orf_name, bindistrib[bindistrib$Var1==x,]$Freq, replace=T))
rownumbers<-NULL
for(x in unlist(l)){rownumbers<-c(rownumbers, which(df$orf_name==x,))}
lengthcontrolled_df<-df[rownumbers,]
lengthcontrolled_df$protogene <-rep("esta_lengthcontrolled", n_toselect)
controlled_df <-rbind(controlled_df, lengthcontrolled_df)

# EXPR_LEVEL
h<-hist(df[df$protogene=="Emerging ORFs",]$expr_level, breaks= seq(-7,6,by=1), plot=FALSE)
protobreaks<-h$breaks
protoprob<-h$counts/nrow(df[df$protogene=="Emerging ORFs",]) 
# Select genes following the proto-gene length distribution, with replacement
t<-df[df $protogene=="Established ORFs" ,]
adj_protobreaks<-c((-7), seq(-5,3,by=1), 6)
adj_protoprob<-c(sum(protoprob[1:2]), protoprob[3:10], sum(protoprob[11:length(protoprob)]))
t$expbin<-cut(t$expr_level, breaks=adj_protobreaks)
t<-t[!(is.na(t$expbin)),]
bindistrib<-as.data.frame(table(sample(x=levels(t$expbin), n_toselect, replace=T, prob= adj_protoprob)))
l<-sapply(bindistrib$Var1, function(x) sample(t[t$expbin==x,]$orf_name, bindistrib[bindistrib$Var1==x,]$Freq, replace=T))
rownumbers<-NULL
for(x in unlist(l)){rownumbers<-c(rownumbers, which(df$orf_name==x,))}
expcontrolled_df<-df[rownumbers,]
expcontrolled_df$protogene <-rep("esta_expcontrolled", n_toselect)
controlled_df <-rbind(controlled_df, expcontrolled_df)

# plot fractions with fixed structures
cutoff<-0.9
stat_df<-controlled_df[,c("protogene","strain_goodalign", "strained_conserved_orf")]
stat_df$is_delet<-ifelse((stat_df$strained_conserved_orf/2022)>cutoff,1,0)
plotdf<-NULL
plotdf $total<- as.vector(table(stat_df $protogene))
plotdf $yes<-as.vector(sapply(levels(stat_df $protogene), function(x) sum(stat_df[stat_df $protogene ==x,]$is_delet)))
plotdf $fraction<-plotdf $yes/plotdf $total
plotdf $sder<-sqrt(plotdf $fraction * (1-plotdf $fraction) / plotdf $total)
plotdf $class<-levels(stat_df $protogene)
plotdf <- data.frame(plotdf)
plotdf$class<-factor(plotdf$class, levels = c("Emerging ORFs","esta_expcontrolled", "esta_lengthcontrolled", "Established ORFs"))
levels(plotdf$class)<-mapvalues(levels(plotdf$class),from=levels(plotdf$class), to= c("Emerging ORFs", "Established ORFs/matched expression", "Established ORFs/matched length", "All established ORFs"))

# plot the fraction with fitness costs
image<-ggplot(plotdf,aes(x= class,y=fraction,ymin=(fraction - sder),ymax=(fraction + sder),fill=class, color=class))+geom_bar(stat="identity",position=position_dodge())+geom_errorbar(position = position_dodge(0.9),width=0.25)+coord_flip()+ labs(x="",y="Fraction with fixed ORF structure")+ scale_color_manual(values=c("deepskyblue4", "black", "black","black")) + scale_fill_manual(values=c("lightblue", "white", "white","white")) +theme_bw()+theme(panel.spacing = unit(0, "lines"), axis.text.x=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.y=element_text(family="Helvetica", face="bold", colour="black", size=7), 	axis.title.x=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.y=element_blank(),	legend.position="none") 
ggsave(paste(general_outdir,"FigS1B_Compare_StrainORFStructure_controlled_barplots.pdf",sep=""),width=4.5,height=2.5,  useDingbats=FALSE)

# statistical comparison: Fisher's exact with fitness cutoff ; between emerging ORFs and expr_controled
n_pg_effect<-nrow(controlled_df[controlled_df $protogene=="Emerging ORFs" & controlled_df $strained_conserved_orf/2022 <cutoff,]) 
n_pg_noeffect<-nrow(controlled_df[controlled_df $protogene=="Emerging ORFs" & controlled_df $strained_conserved_orf/2022>=cutoff,])
n_g_effect<-nrow(controlled_df[controlled_df $protogene=="esta_expcontrolled" & controlled_df $strained_conserved_orf/2022 <cutoff,]) 
n_g_noeffect<-nrow(controlled_df[controlled_df $protogene=="esta_expcontrolled" & controlled_df $strained_conserved_orf/2022>=cutoff,])
pval<-fisher.test(matrix(c(n_pg_effect, n_pg_noeffect, n_g_effect, n_g_noeffect),2,2))

# statistical comparison: Fisher's exact with fitness cutoff ; between emerging ORFs and length_controled
n_pg_effect<-nrow(controlled_df[controlled_df $protogene=="Emerging ORFs" & controlled_df $strained_conserved_orf/2022 <cutoff,]) 
n_pg_noeffect<-nrow(controlled_df[controlled_df $protogene=="Emerging ORFs" & controlled_df $strained_conserved_orf/2022>=cutoff,])
n_g_effect<-nrow(controlled_df[controlled_df $protogene=="esta_lengthcontrolled" & controlled_df $strained_conserved_orf/2022 <cutoff,]) 
n_g_noeffect<-nrow(controlled_df[controlled_df $protogene=="esta_lengthcontrolled" & controlled_df $strained_conserved_orf/2022>=cutoff,])
pval<-fisher.test(matrix(c(n_pg_effect, n_pg_noeffect, n_g_effect, n_g_noeffect),2,2))


######################################################
# Figure S3A: competitive fitness upon overexpression #
######################################################
# data
df<-orf_table[orf_table$barflex_space=="yes",c("protogene","overexpression_competitive_fitness")]

#  statistical comparison: Mann Whitney
wilcoxtest<-wilcox.test(df[df $protogene=="Established ORFs","overexpression_competitive_fitness"],df[df $protogene=="Emerging ORFs","overexpression_competitive_fitness"])
pval<-wilcoxtest$p.val

# plot
mu <- ddply(df, "protogene", summarise, grp.mean=mean(overexpression_competitive_fitness,na.rm=TRUE))
plot_xlabel<-"Competitive fitness upon overexpression"
plot_ylabel<- "ORFs (density)"
image<-ggplot(df, aes(x= overexpression_competitive_fitness,color=protogene,fill=protogene))+geom_density(na.rm=TRUE)+   geom_vline(data=mu, aes(xintercept=grp.mean, color=protogene), linetype="dashed")+scale_color_manual(values=c("deepskyblue4", "black"))+scale_fill_manual(values=c("lightblue", "NA"))+labs(x=plot_xlabel,y=plot_ylabel) + theme_bw()+theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8),  	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")   
ggsave(paste(general_outdir,"FigS3A_Compare_CompetitiveFitness.pdf",sep=""),width=5,height=2.5)


############################################################################
# Figure S3B: showing the exact colony sizes of increased propo-genes #
############################################################################
df<-colony_table
increased_protogenes <- as.character(orf_table[which(orf_table$protogene=="Emerging ORFs" & orf_table$overexpression_relative_fitness=="increased"), "orf_name"])

# labels for data frame
plot_df <- df[df$orf_name %in% increased_protogenes | df$orf_name == "BF_control",]
plot_df$category<-rep("other", nrow(plot_df))
plot_df[plot_df $orf_name=="BF_control",]$category<-"Reference"
plot_df[plot_df $orf_name!="BF_control",]$category<-"Emerging ORFs"
plot_df $category <-factor(plot_df $category, levels=c("Emerging ORFs","Reference"))

# plot
plot_xlabel<-"Normalized colony size"
medd<-median(plot_df[plot_df $orf_name=="BF_control",]$colony_size)
box_scatter <- ggplot() +
  geom_violin(data = plot_df[which(plot_df $category=="Reference"),],
               aes(x=category,
                   y=colony_size),
              fill="grey") +
  geom_point(data = plot_df[which(plot_df $category=="Emerging ORFs"),],
              aes(x=orf_name,
                  y=colony_size),
             colour="deepskyblue4") +
  geom_hline(data=medd, yintercept=medd, linetype="dashed", color="red") +
  labs(x=NULL,y=plot_xlabel) + 
  theme_bw() +
  theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8),
        axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),
        axis.text.x = element_text(angle=60, hjust=1))
ggsave(paste(general_outdir, "FigS3B_increased_fitness_individualcolonies.pdf",sep=""),width=5,height=2.5)


# ##############################################
# # Figure 3A: number of colonies per category overexpressed#
# ##############################################

df<-colony_table

# category for each colony
df$category<-rep("other", nrow(df))
df[df$orf_name=="BF_control",]$category<-"reference"
df[df$orf_name %in% orf_table[orf_table$protogene=="Emerging ORFs",]$orf_name,]$category<-"Emerging ORFs"
df[df$orf_name %in% orf_table[orf_table$protogene=="Established ORFs",]$orf_name,]$category<-"Established ORFs"
df $category <-factor(df $category, levels=c("Emerging ORFs","Established ORFs", "reference"))
table(df$category)

############################################################################
# Figure 3B: fraction of orfs with relative overexpression fitness effects #
############################################################################

df<-orf_table[orf_table$barflex_space=="yes",c("protogene","overexpression_relative_fitness")]
table(df[,1])

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
ggsave(paste(general_outdir,"Fig3B_ExperimentalResults_RelativeFitness_Fraction.pdf",sep=""),width=3.5,height=2.5)

# range of effect sizes for proto-genes (average of replicates)
increased_colonies<-colony_table[colony_table$orf_name %in% orf_table[orf_table $overexpression_relative_fitness=="increased" & orf_table$protogene=="Emerging ORFs",]$orf_name,]
mean_colonysizes<-ddply(increased_colonies, "orf_name", summarise, grp.mean=mean(colony_size,na.rm=TRUE))
min(mean_colonysizes $grp.mean)
max(mean_colonysizes $grp.mean)


##############################################################################
# Figure 3C: Odds Ratio of orfs with relative overexpression fitness effects #
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
ggsave(paste(general_outdir,"Fig3C_ExperimentalResults_RelativeFitness_OddsRatio.pdf",sep=""),width=3.5,height=2.5,  useDingbats=FALSE)

###############################################################################################
# Figure 3D: Odds Ratio of orfs with relative overexpression fitness effects - 5 environments #
###############################################################################################

df<-merge(fitness_table, orf_table[,c("orf_name","protogene")], all.x=TRUE)

# calculate odds ratios and confidence intervals and pvalues
a<-data.frame(table(df[,c("exp_environment","protogene", "effect_cs")]))
totals<-table(df[,c("exp_environment","protogene")])
b<-merge(a,totals,by=c("exp_environment","protogene"),all.x=TRUE)
colnames(b)<-c("exp_environment",  "protogene",   "effect_cs", "number", "total")
b$minus<-b$total - b$number
protos<-b[b$protogene=="Emerging ORFs",c("exp_environment" , "protogene",   "effect_cs", "number", "minus")]
genes<-b[b$protogene=="Established ORFs",c("exp_environment" , "protogene",   "effect_cs", "number", "minus")]
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
ggsave(paste(general_outdir,"Fig3D_ExperimentalResults_RelativeFitness_OddsRatio_5Environments.pdf",sep=""),width=3.5,height=2.5,  useDingbats=FALSE)


####################################################################
# Figure 3E: Adaptive ORFs - length and expression level controlled#
####################################################################

# make a data frame with adaptive in a column
adaptive_orfs<-unique(fitness_table[fitness_table $effect_cs=="increased",]$orf_name)
df<-orf_table[orf_table$barflex_space=="yes",c("orf_name","protogene","Length","expr_level","loss_fitness", "strained_conserved_orf", "strain_pi")]
df$adaptive<-sapply(df$orf_name, function(x) ifelse(x %in% adaptive_orfs,1,0))
df$sick<- ifelse(df $loss_fitness<0.9,1,0)
df$intact<-ifelse(df $strained_conserved_orf <0.9 * 2022 ,0,1)

# enrichment in adaptive effects
table(df[,c("protogene","adaptive")])s
fisher.test(table(df[,c("protogene","adaptive")]))

# statistical comparison: selected effects in adaptive proto-genes
pg_df<-df[df$protogene=="Emerging ORFs",]
# fitness upon loss
fisher.test(as.matrix(table(pg_df[,c("adaptive", "sick")])))
# orf intactness in isolates
fisher.test(as.matrix(table(pg_df[,c("adaptive", "intact")])))
# nucleotide diversity
wilcox.test(pg_df[pg_df$adaptive==1,"strain_pi"], pg_df[pg_df$adaptive==0,"strain_pi"], na.rm=TRUE)


### controlled data (all proto-genes, all genes and randomly picked genes that follow proto-gene length and expr_level distribution)
controlled_df<-df
n_toselect<-nrow(df[df $protogene=="Established ORFs" ,])

#LENGTH
# determine the probability distribution for proto-genes
h<-hist(df[df$protogene=="Emerging ORFs",]$Length, breaks= seq(50,2000,by=50), plot=FALSE)
protobreaks<-h$breaks
protoprob<-h$counts/nrow(df[df$protogene=="Emerging ORFs",]) 
# Select genes following the proto-gene length distribution, with replacement
t<-df[df $protogene=="Established ORFs" ,]
adj_protobreaks<-c(seq(50,800,by=50), 2000)
adj_protoprob<-c(protoprob[1:15], sum(protoprob[16:length(protoprob)]))
t$lengthbin<-cut(t$Length, breaks=adj_protobreaks)
t<-t[!(is.na(t$lengthbin)),]
bindistrib<-as.data.frame(table(sample(x=levels(t$lengthbin), n_toselect, replace=T, prob= adj_protoprob)))
l<-sapply(bindistrib$Var1, function(x) sample(t[t$lengthbin==x,]$orf_name, bindistrib[bindistrib$Var1==x,]$Freq, replace=T))
rownumbers<-NULL
for(x in unlist(l)){rownumbers<-c(rownumbers, which(df$orf_name==x,))}
lengthcontrolled_df<-df[rownumbers,]
lengthcontrolled_df$protogene <-rep("esta_lengthcontrolled", n_toselect)
controlled_df <-rbind(controlled_df, lengthcontrolled_df)

# EXPR_LEVEL
h<-hist(df[df$protogene=="Emerging ORFs",]$expr_level, breaks= seq(-7,6,by=1), plot=FALSE)
protobreaks<-h$breaks
protoprob<-h$counts/nrow(df[df$protogene=="Emerging ORFs",]) 
# Select genes following the proto-gene length distribution, with replacement
t<-df[df $protogene=="Established ORFs" ,]
adj_protobreaks<-c((-7), seq(-5,3,by=1), 6)
adj_protoprob<-c(sum(protoprob[1:2]), protoprob[3:10], sum(protoprob[11:length(protoprob)]))
t$expbin<-cut(t$expr_level, breaks=adj_protobreaks)
t<-t[!(is.na(t$expbin)),]
bindistrib<-as.data.frame(table(sample(x=levels(t$expbin), n_toselect, replace=T, prob= adj_protoprob)))
l<-sapply(bindistrib$Var1, function(x) sample(t[t$expbin==x,]$orf_name, bindistrib[bindistrib$Var1==x,]$Freq, replace=T))
rownumbers<-NULL
for(x in unlist(l)){rownumbers<-c(rownumbers, which(df$orf_name==x,))}
expcontrolled_df<-df[rownumbers,]
expcontrolled_df$protogene <-rep("esta_expcontrolled", n_toselect)
controlled_df <-rbind(controlled_df, expcontrolled_df)

stat_df<-controlled_df[,c("protogene","adaptive")]
plotdf<-NULL
plotdf $total<- as.vector(table(stat_df $protogene))
plotdf $yes<-as.vector(sapply(levels(stat_df $protogene), function(x) sum(stat_df[stat_df $protogene ==x,]$adaptive)))
plotdf $fraction<-plotdf $yes/plotdf $total
plotdf $sder<-sqrt(plotdf $fraction * (1-plotdf $fraction) / plotdf $total)
plotdf $class<-levels(stat_df $protogene)
plotdf <- data.frame(plotdf)
plotdf$class<-factor(plotdf$class, levels = c("Emerging ORFs","esta_expcontrolled", "esta_lengthcontrolled", "Established ORFs"))
levels(plotdf$class)<-mapvalues(levels(plotdf$class),from=levels(plotdf$class), to= c("Emerging ORFs", "Established ORFs/matched expression", "Established ORFs/matched length", "All established ORFs"))

# plot the fraction adaptive
image<-ggplot(plotdf,aes(x= class,y=fraction,ymin=(fraction - sder),ymax=(fraction + sder),fill=class, color=class))+geom_bar(stat="identity",position=position_dodge())+geom_errorbar(position = position_dodge(0.9),width=0.25)+coord_flip()+ labs(x="",y="Fraction adaptive")+ scale_color_manual(values=c("deepskyblue4", "black", "black","black")) + scale_fill_manual(values=c("lightblue", "white", "white","white")) +theme_bw()+theme(panel.spacing = unit(0, "lines"), axis.text.x=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.y=element_text(family="Helvetica", face="bold", colour="black", size=7), 	axis.title.x=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.y=element_blank(),	legend.position="none") 
ggsave(paste(general_outdir,"Fig3E_AdaptiveORFs_controlled_barplots.pdf",sep=""),width=7,height=2.5,  useDingbats=FALSE)

#######################################################
# Figure S4: Adaptive proto-genes across environments #
#######################################################

df<-merge(fitness_table, orf_table[,c("orf_name","protogene")], all.x=TRUE)
adaptive_protogenes_df<-df[df$protogene=="Emerging ORFs" & df$effect_cs=="increased",c("orf_name","exp_environment")]

# plot counts per environment (Fig S2A)
plot_x_label<-"Environmental condition"
plot_y_label<-"# adaptive emerging ORFs"
count_image<-ggplot(adaptive_protogenes_df, aes(exp_environment))+ geom_bar(color='deepskyblue4', fill='lightblue')+theme_bw()+labs(x=plot_x_label,y=plot_y_label) + theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.title.x=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.y=element_text(family="Helvetica", face="bold", colour="black", size=10), legend.position="none") 
ggsave(paste(general_outdir,"FigS4A_AdaptiveEmergingORFs_PerEnvironment_Counts.pdf",sep=""),width=3.5,height=2.5, useDingbats=FALSE)

# distribution of # environements (Fig S2B)
 plot_x_label<-"# environmental conditions"
plot_y_label<-"# adaptive emerging ORFs"
distrib<-table(adaptive_protogenes_df$orf_name)
plot_df<-as.data.frame(distrib)
distrib_image<-ggplot(plot_df[plot_df$Freq>0,], aes(x=Freq))+ geom_histogram(binwidth=1, color='deepskyblue4', fill='lightblue')+theme_bw()+labs(x=plot_x_label,y=plot_y_label) + theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.title.x=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.y=element_text(family="Helvetica", face="bold", colour="black", size=10), legend.position="none") 
ggsave(paste(general_outdir,"FigS4B_AdaptiveEmerging ORFs_PerEnvironment_Distribution.pdf",sep=""),width=3.5,height=2.5, useDingbats=FALSE)

# show that this cannot happen by chance

# categories
adaptive_orfs<- unique(df[df $effect_cs=="increased",]$orf_name)
deleterious_orfs<-unique(df[!(df $orf_name %in% adaptive_orfs) & df $effect_cs=="decreased",]$orf_name)

#get deleterious proto-genes
delet_protos <- unique(df[df$protogene=="Emerging ORFs" & df$orf_name %in% deleterious_orfs,]$orf_name)
# non-deleterious proto-genes
nondelet_protos <- unique(df[df$protogene=="Emerging ORFs" & !(df$orf_name %in% delet_protos),]$orf_name)

#### calculating after removal of all deleterious

resAll <- c()
freqsAll <- c()
for (i in c(1:10000))
{
  samples <- c()
  for (j in c(18, 14, 14, 9, 12)) # these are the actual number of increased fitness protogenes found per condition
  {
    randSampl <- as.character(sample(nondelet_protos, j))
    samples <- c(samples, randSampl)     
  }
  total_found<-length(unique(samples))    
  picked_repeatedly<-length(unique(samples[duplicated(samples)])) 
  freqsAll <- c(freqsAll,picked_repeatedly/total_found)
  resAll <- c(resAll, picked_repeatedly )
}
summary(resAll)
summary(freqsAll)

### plot

df_plot <- data.frame(freqs=freqsAll)

hist <- ggplot() + 
  geom_histogram(data= df_plot,
           aes(x= freqs),
           na.rm = TRUE,
           bins=20) +
  geom_vline(data=df,
             aes(xintercept=17/28), 
             color = "red",
             linetype="dashed") + 
  labs(x="Expected fraction found more than once",y="count") + 
  theme_bw() + 
  theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8),
        axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),
        legend.title = element_blank(), legend.position="top") 

ggsave(paste(general_outdir, "FigS4C_AdaptiveEmergingORFs_Expectation_underNullModel.pdf",sep=""), width=3.5,height=2.5, useDingbats=FALSE)


##########################################################################################
# Fig. 4A - ISD in Emerging ORFs as a function of fitness #
###########################################################################################
# ISD, evolutionary status and fitness category
df<-orf_table[orf_table$barflex_space=="yes" & orf_table$protogene=="Emerging ORFs",c("orf_name", "percent_disorder")]
fitness_df<-merge(fitness_table, orf_table[,c("orf_name","protogene")], all.x=TRUE)
adaptive_orfs<-unique(fitness_df[fitness_df $effect_cs=="increased",]$orf_name)
deleterious_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & fitness_df $effect_cs=="decreased",]$orf_name)
neutral_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & !(fitness_df $orf_name %in% deleterious_orfs),]$orf_name) 
df $fitness_category<-sapply(df $orf_name, function(x) ifelse(x %in% adaptive_orfs,"adaptive",ifelse(x %in% deleterious_orfs, "deleterious", ifelse(x %in% neutral_orfs,"neutral","others"))))
df$fitness_category <-factor(df$fitness_category, levels=c("adaptive", "neutral", "deleterious"))

# calculate average disorder
df2<-df[,c("fitness_category","percent_disorder")]
diso<-NULL
diso$total<- as.vector(table(df2$fitness_category))
diso$Mean<-as.vector(sapply(levels(df2$fitness_category), function(x) mean(df2[df2$fitness_category ==x,]$percent_disorder, na.rm=TRUE)))
diso$SD<-as.vector(sapply(levels(df2$fitness_category), function(x) sd(df2[df2$fitness_category ==x,]$percent_disorder, na.rm=TRUE)))
diso$sder<-diso$SD/sqrt(diso$total)
diso$class<-levels(df$plottingclass)
diso<-data.frame(diso)
#diso$protogene<-c(rep("Proto-genes",3),rep("Genes",3))
diso$fitness<-c("Adaptive","Neutral","Deleterious")
diso$fitness<-factor(diso$fitness, levels=c("Adaptive","Neutral","Deleterious"))

# plot
image<-ggplot(diso,aes(x= fitness,y=Mean,ymin=(Mean - sder),ymax=(Mean + sder),fill=fitness, color=fitness))+geom_bar(stat="identity",position=position_dodge())+geom_errorbar(position = position_dodge(0.9),width=0.25)+labs(x="",y="Average disorder content")+ scale_color_manual(values=c("darkgoldenrod1", "darkgray", "darkorchid3")) + scale_fill_manual(values=c("gold", "lightgray", "mediumpurple1"))+theme_bw()+theme(panel.spacing = unit(0, "lines"), axis.text.y=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.x=element_text(family="Helvetica", face="bold", colour="black", size=7), 	axis.title.y=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.x=element_blank(),	legend.position="none") 
ggsave(paste(general_outdir,"Fig4A_ISD_fitness_barplots.pdf",sep=""),width=2.5,height=2.5,  useDingbats=FALSE)

# statistical comparison : adaptive versus neutral proto-genes wilxoc test
wilcox.test(df[df$fitness_category =="adaptive",]$percent_disorder, df[df$fitness_category =="neutral",]$percent_disorder)
# statistical comparison : adaptive versus deleterious genes wilxoc test
wilcox.test(df[df$fitness_category =="deleterious",]$percent_disorder, df[df$fitness_category =="adaptive",]$percent_disorder)


##########################################################################################
# Fig. 4B - GCcontent in Emerging ORFs as a function of fitness  #
###########################################################################################

df<-orf_table[orf_table$barflex_space=="yes" & orf_table$protogene=="Emerging ORFs",c("orf_name", "GC")]
fitness_df<-merge(fitness_table, orf_table[,c("orf_name","protogene")], all.x=TRUE)
adaptive_orfs<-unique(fitness_df[fitness_df $effect_cs=="increased",]$orf_name)
deleterious_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & fitness_df $effect_cs=="decreased",]$orf_name)
neutral_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & !(fitness_df $orf_name %in% deleterious_orfs),]$orf_name) 
df $fitness_category<-sapply(df $orf_name, function(x) ifelse(x %in% adaptive_orfs,"adaptive",ifelse(x %in% deleterious_orfs, "deleterious", ifelse(x %in% neutral_orfs,"neutral","others"))))
#df$plottingclass<-paste(df $fitness_category, df $protogene, sep="_")
df$fitness_category <-factor(df$fitness_category, levels=c("adaptive", "neutral", "deleterious"))

# calculate average GC
df2<-df[,c("fitness_category","GC")]
diso<-NULL
diso$total<- as.vector(table(df2$fitness_category))
diso$Mean<-as.vector(sapply(levels(df2$fitness_category), function(x) mean(df2[df2$fitness_category ==x,]$GC, na.rm=TRUE)))
diso$SD<-as.vector(sapply(levels(df2$fitness_category), function(x) sd(df2[df2$fitness_category ==x,]$GC, na.rm=TRUE)))
diso$sder<-diso$SD/sqrt(diso$total)
diso$class<-levels(df$plottingclass)
diso<-data.frame(diso)
diso$fitness<-c("Adaptive","Neutral","Deleterious")
diso$fitness<-factor(diso$fitness, levels=c("Adaptive","Neutral","Deleterious"))

# plot
image<-ggplot(diso,aes(x= fitness,y=Mean,ymin=(Mean - sder),ymax=(Mean + sder),fill=fitness, color=fitness))+geom_bar(stat="identity",position=position_dodge())+geom_errorbar(position = position_dodge(0.9),width=0.25)+labs(x="",y="Average GC content")+ scale_color_manual(values=c("darkgoldenrod1", "darkgray", "darkorchid3")) + scale_fill_manual(values=c("gold", "lightgray", "mediumpurple1"))+theme_bw()+theme(panel.spacing = unit(0, "lines"), axis.text.y=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.x=element_text(family="Helvetica", face="bold", colour="black", size=7), 	axis.title.y=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.x=element_blank(),	legend.position="none") 
ggsave(paste(general_outdir,"Fig4B_GC_fitness_barplots.pdf",sep=""),width=2.5,height=2.5,  useDingbats=FALSE)

# statistical comparison : adaptive versus neutral proto-genes wilxoc test
wilcox.test(df[df$fitness_category =="adaptive",]$GC, df[df$fitness_category =="neutral",]$GC)

##########################################################################################
# Fig. 4C - TM residue content in Emerging ORFs as a function of fitness = TMHMM #
###########################################################################################
# TMHMM, evolutionary status and fitness category
df<-orf_table[orf_table$barflex_space=="yes" & orf_table$protogene=="Emerging ORFs",c("orf_name", "percent_tmhmm")]
fitness_df<-merge(fitness_table, orf_table[,c("orf_name","protogene")], all.x=TRUE)
adaptive_orfs<-unique(fitness_df[fitness_df $effect_cs=="increased",]$orf_name)
deleterious_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & fitness_df $effect_cs=="decreased",]$orf_name)
neutral_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & !(fitness_df $orf_name %in% deleterious_orfs),]$orf_name) 
df $fitness_category<-sapply(df $orf_name, function(x) ifelse(x %in% adaptive_orfs,"adaptive",ifelse(x %in% deleterious_orfs, "deleterious", ifelse(x %in% neutral_orfs,"neutral","others"))))
df$fitness_category <-factor(df$fitness_category, levels=c("adaptive", "neutral", "deleterious"))

# calculate average disorder
df2<-df[,c("fitness_category","percent_tmhmm")]
diso<-NULL
diso$total<- as.vector(table(df2$fitness_category))
diso$Mean<-as.vector(sapply(levels(df2$fitness_category), function(x) mean(df2[df2$fitness_category ==x,]$percent_tmhmm, na.rm=TRUE)))
diso$SD<-as.vector(sapply(levels(df2$fitness_category), function(x) sd(df2[df2$fitness_category ==x,]$percent_tmhmm, na.rm=TRUE)))
diso$sder<-diso$SD/sqrt(diso$total)
diso$class<-levels(df$plottingclass)
diso<-data.frame(diso)
diso$fitness<-c("Adaptive","Neutral","Deleterious")
diso$fitness<-factor(diso$fitness, levels=c("Adaptive","Neutral","Deleterious"))

# plot
image<-ggplot(diso,aes(x= fitness,y=Mean,ymin=(Mean - sder),ymax=(Mean + sder),fill=fitness, color=fitness))+geom_bar(stat="identity",position=position_dodge())+geom_errorbar(position = position_dodge(0.9),width=0.25)+labs(x="",y="Average TM content (TMHMM)")+ scale_color_manual(values=c("darkgoldenrod1", "darkgray", "darkorchid3")) + scale_fill_manual(values=c("gold", "lightgray", "mediumpurple1"))+theme_bw()+theme(panel.spacing = unit(0, "lines"), axis.text.y=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.x=element_text(family="Helvetica", face="bold", colour="black", size=7), 	axis.title.y=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.x=element_blank(),	legend.position="none") 
ggsave(paste(general_outdir,"Fig4C_TMcontent_TMHMM_fitness_barplots.pdf",sep=""),width=2.5,height=2.5,  useDingbats=FALSE)

# statistical comparison : adaptive versus neutral proto-genes wilxoc test
wilcox.test(df[df$fitness_category =="adaptive",]$percent_tmhmm, df[df$fitness_category =="neutral",]$percent_tmhmm)
# statistical comparison : deleterious versus neutral genes wilxoc test
wilcox.test(df[df$fitness_category =="deleterious",]$percent_tmhmm, df[df$fitness_category =="neutral",]$percent_tmhmm)


##########################################################################################
# Fig. 4D - TM residue content in Emerging ORFs as a function of fitness = Phobius #
###########################################################################################
# percent_phobius, evolutionary status and fitness category
df<-orf_table[orf_table$barflex_space=="yes" & orf_table$protogene=="Emerging ORFs",c("orf_name", "percent_phobius")]
fitness_df<-merge(fitness_table, orf_table[,c("orf_name","protogene")], all.x=TRUE)
adaptive_orfs<-unique(fitness_df[fitness_df $effect_cs=="increased",]$orf_name)
deleterious_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & fitness_df $effect_cs=="decreased",]$orf_name)
neutral_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & !(fitness_df $orf_name %in% deleterious_orfs),]$orf_name) 
df $fitness_category<-sapply(df $orf_name, function(x) ifelse(x %in% adaptive_orfs,"adaptive",ifelse(x %in% deleterious_orfs, "deleterious", ifelse(x %in% neutral_orfs,"neutral","others"))))
df$fitness_category <-factor(df$fitness_category, levels=c("adaptive", "neutral", "deleterious"))

# calculate average disorder
df2<-df[,c("fitness_category","percent_phobius")]
diso<-NULL
diso$total<- as.vector(table(df2$fitness_category))
diso$Mean<-as.vector(sapply(levels(df2$fitness_category), function(x) mean(df2[df2$fitness_category ==x,]$percent_phobius, na.rm=TRUE)))
diso$SD<-as.vector(sapply(levels(df2$fitness_category), function(x) sd(df2[df2$fitness_category ==x,]$percent_phobius, na.rm=TRUE)))
diso$sder<-diso$SD/sqrt(diso$total)
diso$class<-levels(df$plottingclass)
diso<-data.frame(diso)
diso$fitness<-c("Adaptive","Neutral","Deleterious")
diso$fitness<-factor(diso$fitness, levels=c("Adaptive","Neutral","Deleterious"))

# plot
image<-ggplot(diso,aes(x= fitness,y=Mean,ymin=(Mean - sder),ymax=(Mean + sder),fill=fitness, color=fitness))+geom_bar(stat="identity",position=position_dodge())+geom_errorbar(position = position_dodge(0.9),width=0.25)+labs(x="",y="Average TM content (Phobius)")+ scale_color_manual(values=c("darkgoldenrod1", "darkgray", "darkorchid3")) + scale_fill_manual(values=c("gold", "lightgray", "mediumpurple1"))+theme_bw()+theme(panel.spacing = unit(0, "lines"), axis.text.y=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.x=element_text(family="Helvetica", face="bold", colour="black", size=7), 	axis.title.y=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.x=element_blank(),	legend.position="none") 
ggsave(paste(general_outdir,"Fig4D_TMcontent_Phobius_fitness_barplots.pdf",sep=""),width=2.5,height=2.5,  useDingbats=FALSE)

# statistical comparison : adaptive versus neutral proto-genes wilxoc test
wilcox.test(df[df$fitness_category =="adaptive",]$percent_phobius, df[df$fitness_category =="neutral",]$percent_phobius)
# statistical comparison : deleterious versus neutral genes wilxoc test
wilcox.test(df[df$fitness_category =="deleterious",]$percent_phobius, df[df$fitness_category =="neutral",]$percent_phobius)

##########################################################################################
# Fig. 4E. Transmembrane (TM) domains in Emerging ORFs as a function of fitness - tmhmm #
###########################################################################################
# TM, evolutionary status and fitness category
df<-orf_table[orf_table$barflex_space=="yes" & orf_table$protogene=="Emerging ORFs",c("orf_name","no_helices_TMHMM")]
df$is_tm<-ifelse(df$no_helices_TMHMM>0,1,0)
fitness_df<-merge(fitness_table, orf_table[,c("orf_name","protogene")], all.x=TRUE)
adaptive_orfs<-unique(fitness_df[fitness_df $effect_cs=="increased",]$orf_name)
deleterious_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & fitness_df $effect_cs=="decreased",]$orf_name)
neutral_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & !(fitness_df $orf_name %in% deleterious_orfs),]$orf_name) 
df $fitness_category<-sapply(df $orf_name, function(x) ifelse(x %in% adaptive_orfs,"adaptive",ifelse(x %in% deleterious_orfs, "deleterious", ifelse(x %in% neutral_orfs,"neutral","others"))))
df$fitness_category <-factor(df$fitness_category, levels=c("adaptive", "neutral", "deleterious"))

# calculate proportions with at least 1 TM domain
stat_df<-df[,c("fitness_category","is_tm")]
domains<-NULL
domains $total<- as.vector(table(stat_df $fitness_category))
domains $yes<-as.vector(sapply(levels(stat_df $fitness_category), function(x) sum(stat_df[stat_df $fitness_category ==x,]$is_tm)))
domains $fraction<-domains $yes/domains $total
domains $sder<-sqrt(domains $fraction * (1-domains $fraction) / domains $total)
domains <- data.frame(domains)
domains $fitness<-c("Adaptive","Neutral","Deleterious")
domains $fitness<-factor(domains $fitness, levels=c("Adaptive","Neutral","Deleterious"))

# plot
image<-ggplot(domains,aes(x= fitness,y=fraction,ymin=(fraction - sder),ymax=(fraction + sder),fill=fitness, color=fitness))+geom_bar(stat="identity",position=position_dodge())+geom_errorbar(position = position_dodge(0.9),width=0.25)+labs(x="",y="With TM domain (TMHMMM)")+ scale_color_manual(values=c("darkgoldenrod1", "darkgray", "darkorchid3")) + scale_fill_manual(values=c("gold", "lightgray", "mediumpurple1"))+theme_bw()+theme(panel.spacing = unit(0, "lines"), axis.text.y=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.x=element_text(family="Helvetica", face="bold", colour="black", size=7), 	axis.title.y=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.x=element_blank(),	legend.position="none") 
ggsave(paste(general_outdir,"Fig4E_TMdomains_TMHMM_fitness_barplots.pdf",sep=""),width=2.5,height=2.5,  useDingbats=FALSE)

# statistical comparison : adaptive versus neutral proto-genes Fisher test
fisher.test(matrix(c(domains[domains $fitness =="Adaptive",]$yes, domains[domains $fitness =="Neutral",]$yes, domains[domains $fitness =="Adaptive",]$total - domains[domains $fitness =="Adaptive",]$yes, domains[domains $fitness =="Neutral",]$total - domains[domains $fitness =="Neutral",]$yes),2,2))

##########################################################################################
# Fig. 4F. Transmembrane (TM) domains in Emerging ORFs as a function of fitness - Phobius #
###########################################################################################
# TM, evolutionary status and fitness category
df<-orf_table[orf_table$barflex_space=="yes" & orf_table$protogene=="Emerging ORFs",c("orf_name","no_helices_PHOBIUS")]
df$is_tm<-ifelse(df$no_helices_PHOBIUS>0,1,0)
fitness_df<-merge(fitness_table, orf_table[,c("orf_name","protogene")], all.x=TRUE)
adaptive_orfs<-unique(fitness_df[fitness_df $effect_cs=="increased",]$orf_name)
deleterious_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & fitness_df $effect_cs=="decreased",]$orf_name)
neutral_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & !(fitness_df $orf_name %in% deleterious_orfs),]$orf_name) 
df $fitness_category<-sapply(df $orf_name, function(x) ifelse(x %in% adaptive_orfs,"adaptive",ifelse(x %in% deleterious_orfs, "deleterious", ifelse(x %in% neutral_orfs,"neutral","others"))))
df$fitness_category <-factor(df$fitness_category, levels=c("adaptive", "neutral", "deleterious"))

# calculate proportions with at least 1 TM domain
stat_df<-df[,c("fitness_category","is_tm")]
domains<-NULL
domains $total<- as.vector(table(stat_df $fitness_category))
domains $yes<-as.vector(sapply(levels(stat_df $fitness_category), function(x) sum(stat_df[stat_df $fitness_category ==x,]$is_tm)))
domains $fraction<-domains $yes/domains $total
domains $sder<-sqrt(domains $fraction * (1-domains $fraction) / domains $total)
domains <- data.frame(domains)
domains $fitness<-c("Adaptive","Neutral","Deleterious")
domains $fitness<-factor(domains $fitness, levels=c("Adaptive","Neutral","Deleterious"))

# plot
image<-ggplot(domains,aes(x= fitness,y=fraction,ymin=(fraction - sder),ymax=(fraction + sder),fill=fitness, color=fitness))+geom_bar(stat="identity",position=position_dodge())+geom_errorbar(position = position_dodge(0.9),width=0.25)+labs(x="",y="With TM domain (Phobius)")+ scale_color_manual(values=c("darkgoldenrod1", "darkgray", "darkorchid3")) + scale_fill_manual(values=c("gold", "lightgray", "mediumpurple1"))+theme_bw()+theme(panel.spacing = unit(0, "lines"), axis.text.y=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.x=element_text(family="Helvetica", face="bold", colour="black", size=7), 	axis.title.y=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.x=element_blank(),	legend.position="none") 
ggsave(paste(general_outdir,"Fig4F_TMdomains_Phobius_fitness_barplots.pdf",sep=""),width=2.5,height=2.5,  useDingbats=FALSE)

# statistical comparison : adaptive versus neutral proto-genes Fisher test
fisher.test(matrix(c(domains[domains $fitness =="Adaptive",]$yes, domains[domains $fitness =="Neutral",]$yes, domains[domains $fitness =="Adaptive",]$total - domains[domains $fitness =="Adaptive",]$yes, domains[domains $fitness =="Neutral",]$total - domains[domains $fitness =="Neutral",]$yes),2,2))


##########################################################################################
# Fig. S5A - ISD in established ORFs as a function of fitness #
###########################################################################################
# ISD, evolutionary status and fitness category
df<-orf_table[orf_table$barflex_space=="yes" & orf_table$protogene=="Established ORFs",c("orf_name", "percent_disorder")]
fitness_df<-merge(fitness_table, orf_table[,c("orf_name","protogene")], all.x=TRUE)
adaptive_orfs<-unique(fitness_df[fitness_df $effect_cs=="increased",]$orf_name)
deleterious_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & fitness_df $effect_cs=="decreased",]$orf_name)
neutral_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & !(fitness_df $orf_name %in% deleterious_orfs),]$orf_name) 
df $fitness_category<-sapply(df $orf_name, function(x) ifelse(x %in% adaptive_orfs,"adaptive",ifelse(x %in% deleterious_orfs, "deleterious", ifelse(x %in% neutral_orfs,"neutral","others"))))
df$fitness_category <-factor(df$fitness_category, levels=c("adaptive", "neutral", "deleterious"))

# calculate average disorder
df2<-df[,c("fitness_category","percent_disorder")]
diso<-NULL
diso$total<- as.vector(table(df2$fitness_category))
diso$Mean<-as.vector(sapply(levels(df2$fitness_category), function(x) mean(df2[df2$fitness_category ==x,]$percent_disorder, na.rm=TRUE)))
diso$SD<-as.vector(sapply(levels(df2$fitness_category), function(x) sd(df2[df2$fitness_category ==x,]$percent_disorder, na.rm=TRUE)))
diso$sder<-diso$SD/sqrt(diso$total)
diso$class<-levels(df$plottingclass)
diso<-data.frame(diso)
diso$fitness<-c("Adaptive","Neutral","Deleterious")
diso$fitness<-factor(diso$fitness, levels=c("Adaptive","Neutral","Deleterious"))

# plot
image<-ggplot(diso,aes(x= fitness,y=Mean,ymin=(Mean - sder),ymax=(Mean + sder),fill=fitness, color=fitness))+geom_bar(stat="identity",position=position_dodge())+geom_errorbar(position = position_dodge(0.9),width=0.25)+labs(x="",y="Average disorder content")+ scale_color_manual(values=c("darkgoldenrod1", "darkgray", "darkorchid3")) + scale_fill_manual(values=c("gold", "lightgray", "mediumpurple1"))+theme_bw()+theme(panel.spacing = unit(0, "lines"), axis.text.y=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.x=element_text(family="Helvetica", face="bold", colour="black", size=7), 	axis.title.y=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.x=element_blank(),	legend.position="none") 
ggsave(paste(general_outdir,"FigS5A_ISD_fitness_barplots.pdf",sep=""),width=2.5,height=2.5,  useDingbats=FALSE)

# statistical comparison : adaptive versus neutral proto-genes wilxoc test
wilcox.test(df[df$fitness_category =="adaptive",]$percent_disorder, df[df$fitness_category =="neutral",]$percent_disorder)
# statistical comparison : adaptive versus deleterious genes wilxoc test
wilcox.test(df[df$fitness_category =="deleterious",]$percent_disorder, df[df$fitness_category =="adaptive",]$percent_disorder)


##########################################################################################
# Fig. S5B - GCcontent in established orfs as a function of fitness  #
###########################################################################################

df<-orf_table[orf_table$barflex_space=="yes" & orf_table$protogene=="Established ORFs",c("orf_name", "GC")]
fitness_df<-merge(fitness_table, orf_table[,c("orf_name","protogene")], all.x=TRUE)
adaptive_orfs<-unique(fitness_df[fitness_df $effect_cs=="increased",]$orf_name)
deleterious_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & fitness_df $effect_cs=="decreased",]$orf_name)
neutral_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & !(fitness_df $orf_name %in% deleterious_orfs),]$orf_name) 
df $fitness_category<-sapply(df $orf_name, function(x) ifelse(x %in% adaptive_orfs,"adaptive",ifelse(x %in% deleterious_orfs, "deleterious", ifelse(x %in% neutral_orfs,"neutral","others"))))
#df$plottingclass<-paste(df $fitness_category, df $protogene, sep="_")
df$fitness_category <-factor(df$fitness_category, levels=c("adaptive", "neutral", "deleterious"))

# calculate average GC
df2<-df[,c("fitness_category","GC")]
diso<-NULL
diso$total<- as.vector(table(df2$fitness_category))
diso$Mean<-as.vector(sapply(levels(df2$fitness_category), function(x) mean(df2[df2$fitness_category ==x,]$GC, na.rm=TRUE)))
diso$SD<-as.vector(sapply(levels(df2$fitness_category), function(x) sd(df2[df2$fitness_category ==x,]$GC, na.rm=TRUE)))
diso$sder<-diso$SD/sqrt(diso$total)
diso$class<-levels(df$plottingclass)
diso<-data.frame(diso)
diso$fitness<-c("Adaptive","Neutral","Deleterious")
diso$fitness<-factor(diso$fitness, levels=c("Adaptive","Neutral","Deleterious"))

# plot
image<-ggplot(diso,aes(x= fitness,y=Mean,ymin=(Mean - sder),ymax=(Mean + sder),fill=fitness, color=fitness))+geom_bar(stat="identity",position=position_dodge())+geom_errorbar(position = position_dodge(0.9),width=0.25)+labs(x="",y="Average GC content")+ scale_color_manual(values=c("darkgoldenrod1", "darkgray", "darkorchid3")) + scale_fill_manual(values=c("gold", "lightgray", "mediumpurple1"))+theme_bw()+theme(panel.spacing = unit(0, "lines"), axis.text.y=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.x=element_text(family="Helvetica", face="bold", colour="black", size=7), 	axis.title.y=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.x=element_blank(),	legend.position="none") 
ggsave(paste(general_outdir,"FigS5B_GC_fitness_barplots.pdf",sep=""),width=2.5,height=2.5,  useDingbats=FALSE)

# statistical comparison : adaptive versus neutral proto-genes wilxoc test
wilcox.test(df[df$fitness_category =="adaptive",]$GC, df[df$fitness_category =="neutral",]$GC)

##########################################################################################
# Fig. S5C - TM residue content in established as a function of fitness = TMHMM #
###########################################################################################
# TMHMM, evolutionary status and fitness category
df<-orf_table[orf_table$barflex_space=="yes" & orf_table$protogene=="Established ORFs",c("orf_name", "percent_tmhmm")]
fitness_df<-merge(fitness_table, orf_table[,c("orf_name","protogene")], all.x=TRUE)
adaptive_orfs<-unique(fitness_df[fitness_df $effect_cs=="increased",]$orf_name)
deleterious_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & fitness_df $effect_cs=="decreased",]$orf_name)
neutral_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & !(fitness_df $orf_name %in% deleterious_orfs),]$orf_name) 
df $fitness_category<-sapply(df $orf_name, function(x) ifelse(x %in% adaptive_orfs,"adaptive",ifelse(x %in% deleterious_orfs, "deleterious", ifelse(x %in% neutral_orfs,"neutral","others"))))
df$fitness_category <-factor(df$fitness_category, levels=c("adaptive", "neutral", "deleterious"))

# calculate average disorder
df2<-df[,c("fitness_category","percent_tmhmm")]
diso<-NULL
diso$total<- as.vector(table(df2$fitness_category))
diso$Mean<-as.vector(sapply(levels(df2$fitness_category), function(x) mean(df2[df2$fitness_category ==x,]$percent_tmhmm, na.rm=TRUE)))
diso$SD<-as.vector(sapply(levels(df2$fitness_category), function(x) sd(df2[df2$fitness_category ==x,]$percent_tmhmm, na.rm=TRUE)))
diso$sder<-diso$SD/sqrt(diso$total)
diso$class<-levels(df$plottingclass)
diso<-data.frame(diso)
diso$fitness<-c("Adaptive","Neutral","Deleterious")
diso$fitness<-factor(diso$fitness, levels=c("Adaptive","Neutral","Deleterious"))

# plot
image<-ggplot(diso,aes(x= fitness,y=Mean,ymin=(Mean - sder),ymax=(Mean + sder),fill=fitness, color=fitness))+geom_bar(stat="identity",position=position_dodge())+geom_errorbar(position = position_dodge(0.9),width=0.25)+labs(x="",y="Average TM content (TMHMM)")+ scale_color_manual(values=c("darkgoldenrod1", "darkgray", "darkorchid3")) + scale_fill_manual(values=c("gold", "lightgray", "mediumpurple1"))+theme_bw()+theme(panel.spacing = unit(0, "lines"), axis.text.y=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.x=element_text(family="Helvetica", face="bold", colour="black", size=7), 	axis.title.y=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.x=element_blank(),	legend.position="none") 
ggsave(paste(general_outdir,"FigS5C_TMcontent_TMHMM_fitness_barplots.pdf",sep=""),width=2.5,height=2.5,  useDingbats=FALSE)

# statistical comparison : adaptive versus neutral proto-genes wilxoc test
wilcox.test(df[df$fitness_category =="adaptive",]$percent_tmhmm, df[df$fitness_category =="neutral",]$percent_tmhmm)
# statistical comparison : deleterious versus neutral genes wilxoc test
wilcox.test(df[df$fitness_category =="deleterious",]$percent_tmhmm, df[df$fitness_category =="neutral",]$percent_tmhmm)


##########################################################################################
# Fig. S5D - TM residue content in established orfs as a function of fitness = Phobius #
###########################################################################################
# percent_phobius, evolutionary status and fitness category
df<-orf_table[orf_table$barflex_space=="yes" & orf_table$protogene=="Established ORFs",c("orf_name", "percent_phobius")]
fitness_df<-merge(fitness_table, orf_table[,c("orf_name","protogene")], all.x=TRUE)
adaptive_orfs<-unique(fitness_df[fitness_df $effect_cs=="increased",]$orf_name)
deleterious_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & fitness_df $effect_cs=="decreased",]$orf_name)
neutral_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & !(fitness_df $orf_name %in% deleterious_orfs),]$orf_name) 
df $fitness_category<-sapply(df $orf_name, function(x) ifelse(x %in% adaptive_orfs,"adaptive",ifelse(x %in% deleterious_orfs, "deleterious", ifelse(x %in% neutral_orfs,"neutral","others"))))
df$fitness_category <-factor(df$fitness_category, levels=c("adaptive", "neutral", "deleterious"))

# calculate average disorder
df2<-df[,c("fitness_category","percent_phobius")]
diso<-NULL
diso$total<- as.vector(table(df2$fitness_category))
diso$Mean<-as.vector(sapply(levels(df2$fitness_category), function(x) mean(df2[df2$fitness_category ==x,]$percent_phobius, na.rm=TRUE)))
diso$SD<-as.vector(sapply(levels(df2$fitness_category), function(x) sd(df2[df2$fitness_category ==x,]$percent_phobius, na.rm=TRUE)))
diso$sder<-diso$SD/sqrt(diso$total)
diso$class<-levels(df$plottingclass)
diso<-data.frame(diso)
diso$fitness<-c("Adaptive","Neutral","Deleterious")
diso$fitness<-factor(diso$fitness, levels=c("Adaptive","Neutral","Deleterious"))

# plot
image<-ggplot(diso,aes(x= fitness,y=Mean,ymin=(Mean - sder),ymax=(Mean + sder),fill=fitness, color=fitness))+geom_bar(stat="identity",position=position_dodge())+geom_errorbar(position = position_dodge(0.9),width=0.25)+labs(x="",y="Average TM content (Phobius)")+ scale_color_manual(values=c("darkgoldenrod1", "darkgray", "darkorchid3")) + scale_fill_manual(values=c("gold", "lightgray", "mediumpurple1"))+theme_bw()+theme(panel.spacing = unit(0, "lines"), axis.text.y=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.x=element_text(family="Helvetica", face="bold", colour="black", size=7), 	axis.title.y=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.x=element_blank(),	legend.position="none") 
ggsave(paste(general_outdir,"FigS5D_TMcontent_Phobius_fitness_barplots.pdf",sep=""),width=2.5,height=2.5,  useDingbats=FALSE)

# statistical comparison : adaptive versus neutral proto-genes wilxoc test
wilcox.test(df[df$fitness_category =="adaptive",]$percent_phobius, df[df$fitness_category =="neutral",]$percent_phobius)
# statistical comparison : deleterious versus neutral genes wilxoc test
wilcox.test(df[df$fitness_category =="deleterious",]$percent_phobius, df[df$fitness_category =="neutral",]$percent_phobius)

##########################################################################################
# Fig. S5E. Transmembrane (TM) domains in established orfs as a function of fitness - tmhmm #
###########################################################################################
# TM, evolutionary status and fitness category
df<-orf_table[orf_table$barflex_space=="yes" & orf_table$protogene=="Established ORFs",c("orf_name","no_helices_TMHMM")]
df$is_tm<-ifelse(df$no_helices_TMHMM>0,1,0)
fitness_df<-merge(fitness_table, orf_table[,c("orf_name","protogene")], all.x=TRUE)
adaptive_orfs<-unique(fitness_df[fitness_df $effect_cs=="increased",]$orf_name)
deleterious_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & fitness_df $effect_cs=="decreased",]$orf_name)
neutral_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & !(fitness_df $orf_name %in% deleterious_orfs),]$orf_name) 
df $fitness_category<-sapply(df $orf_name, function(x) ifelse(x %in% adaptive_orfs,"adaptive",ifelse(x %in% deleterious_orfs, "deleterious", ifelse(x %in% neutral_orfs,"neutral","others"))))
df$fitness_category <-factor(df$fitness_category, levels=c("adaptive", "neutral", "deleterious"))

# calculate proportions with at least 1 TM domain
stat_df<-df[,c("fitness_category","is_tm")]
domains<-NULL
domains $total<- as.vector(table(stat_df $fitness_category))
domains $yes<-as.vector(sapply(levels(stat_df $fitness_category), function(x) sum(stat_df[stat_df $fitness_category ==x,]$is_tm)))
domains $fraction<-domains $yes/domains $total
domains $sder<-sqrt(domains $fraction * (1-domains $fraction) / domains $total)
domains <- data.frame(domains)
domains $fitness<-c("Adaptive","Neutral","Deleterious")
domains $fitness<-factor(domains $fitness, levels=c("Adaptive","Neutral","Deleterious"))

# plot
image<-ggplot(domains,aes(x= fitness,y=fraction,ymin=(fraction - sder),ymax=(fraction + sder),fill=fitness, color=fitness))+geom_bar(stat="identity",position=position_dodge())+geom_errorbar(position = position_dodge(0.9),width=0.25)+labs(x="",y="With TM domain (TMHMMM)")+ scale_color_manual(values=c("darkgoldenrod1", "darkgray", "darkorchid3")) + scale_fill_manual(values=c("gold", "lightgray", "mediumpurple1"))+theme_bw()+theme(panel.spacing = unit(0, "lines"), axis.text.y=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.x=element_text(family="Helvetica", face="bold", colour="black", size=7), 	axis.title.y=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.x=element_blank(),	legend.position="none") 
ggsave(paste(general_outdir,"FigS5E_TMdomains_TMHMM_fitness_barplots.pdf",sep=""),width=2.5,height=2.5,  useDingbats=FALSE)

# statistical comparison : adaptive versus neutral proto-genes Fisher test
fisher.test(matrix(c(domains[domains $fitness =="Adaptive",]$yes, domains[domains $fitness =="Neutral",]$yes, domains[domains $fitness =="Adaptive",]$total - domains[domains $fitness =="Adaptive",]$yes, domains[domains $fitness =="Neutral",]$total - domains[domains $fitness =="Neutral",]$yes),2,2))

##########################################################################################
# Fig. S5F. Transmembrane (TM) domains in established ORFs as a function of fitness - Phobius #
###########################################################################################
# TM, evolutionary status and fitness category
df<-orf_table[orf_table$barflex_space=="yes" & orf_table$protogene=="Established ORFs",c("orf_name","no_helices_PHOBIUS")]
df$is_tm<-ifelse(df$no_helices_PHOBIUS>0,1,0)
fitness_df<-merge(fitness_table, orf_table[,c("orf_name","protogene")], all.x=TRUE)
adaptive_orfs<-unique(fitness_df[fitness_df $effect_cs=="increased",]$orf_name)
deleterious_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & fitness_df $effect_cs=="decreased",]$orf_name)
neutral_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & !(fitness_df $orf_name %in% deleterious_orfs),]$orf_name) 
df $fitness_category<-sapply(df $orf_name, function(x) ifelse(x %in% adaptive_orfs,"adaptive",ifelse(x %in% deleterious_orfs, "deleterious", ifelse(x %in% neutral_orfs,"neutral","others"))))
df$fitness_category <-factor(df$fitness_category, levels=c("adaptive", "neutral", "deleterious"))

# calculate proportions with at least 1 TM domain
stat_df<-df[,c("fitness_category","is_tm")]
domains<-NULL
domains $total<- as.vector(table(stat_df $fitness_category))
domains $yes<-as.vector(sapply(levels(stat_df $fitness_category), function(x) sum(stat_df[stat_df $fitness_category ==x,]$is_tm)))
domains $fraction<-domains $yes/domains $total
domains $sder<-sqrt(domains $fraction * (1-domains $fraction) / domains $total)
domains <- data.frame(domains)
domains $fitness<-c("Adaptive","Neutral","Deleterious")
domains $fitness<-factor(domains $fitness, levels=c("Adaptive","Neutral","Deleterious"))

# plot
image<-ggplot(domains,aes(x= fitness,y=fraction,ymin=(fraction - sder),ymax=(fraction + sder),fill=fitness, color=fitness))+geom_bar(stat="identity",position=position_dodge())+geom_errorbar(position = position_dodge(0.9),width=0.25)+labs(x="",y="With TM domain (Phobius)")+ scale_color_manual(values=c("darkgoldenrod1", "darkgray", "darkorchid3")) + scale_fill_manual(values=c("gold", "lightgray", "mediumpurple1"))+theme_bw()+theme(panel.spacing = unit(0, "lines"), axis.text.y=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.x=element_text(family="Helvetica", face="bold", colour="black", size=7), 	axis.title.y=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.x=element_blank(),	legend.position="none") 
ggsave(paste(general_outdir,"FigS5F_TMdomains_Phobius_fitness_barplots.pdf",sep=""),width=2.5,height=2.5,  useDingbats=FALSE)

# statistical comparison : adaptive versus neutral proto-genes Fisher test
fisher.test(matrix(c(domains[domains $fitness =="Adaptive",]$yes, domains[domains $fitness =="Neutral",]$yes, domains[domains $fitness =="Adaptive",]$total - domains[domains $fitness =="Adaptive",]$yes, domains[domains $fitness =="Neutral",]$total - domains[domains $fitness =="Neutral",]$yes),2,2))



###########################################################################################
# Fig. S6 Distribution of %TM residues in genes and proto-genes as a function of fitness#
###########################################################################################

###data
df<-orf_table[orf_table$barflex_space=="yes",c("orf_name","protogene", "percent_phobius", "percent_tmhmm")]
fitness_df<-merge(fitness_table, orf_table[,c("orf_name","protogene")], all.x=TRUE)
adaptive_orfs<-unique(fitness_df[fitness_df $effect_cs=="increased",]$orf_name)
deleterious_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & fitness_df $effect_cs=="decreased",]$orf_name)
neutral_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & !(fitness_df $orf_name %in% deleterious_orfs),]$orf_name) 
df $class<-sapply(df $orf_name, function(x) ifelse(x %in% adaptive_orfs,"adaptive",ifelse(x %in% deleterious_orfs, "deleterious", ifelse(x %in% neutral_orfs,"neutral","others"))))

### Plots

# emerging using phobius
thisdf<-df[df$protogene=="Emerging ORFs" , c("class","percent_phobius")]
thismu<-data.frame(class=c("adaptive","neutral","deleterious"),Mean = c(mean(thisdf[thisdf$class=="adaptive","percent_phobius"],na.rm=TRUE), mean(thisdf[thisdf$class=="neutral","percent_phobius"],na.rm=TRUE),mean(thisdf[thisdf$class=="deleterious","percent_phobius"],na.rm=TRUE) ))
plot_xlabel <-""
plot_ylabel<-"Emerging ORFs"
image_protoTMphobius<-ggplot(thisdf, aes(x= percent_phobius,color= class,fill= class))+geom_density(na.rm=TRUE)+
   geom_vline(data=thismu, aes(xintercept=Mean, color= class), linetype="dashed")+scale_color_manual(values=c("darkgoldenrod1", "darkorchid3","darkgray"))+xlim(c(0,0.8))+scale_fill_manual(values=c("gold","NA", "NA"))+labs(x=plot_xlabel,y=plot_ylabel) + theme_bw()+theme(axis.text.x=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.y=element_blank(),axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")+  theme(panel.background = element_rect(fill = "lightcyan2"), panel.grid.major = element_line(colour = "white"),   panel.grid.minor = element_line(colour = "white"))

# emerging using tmhmm
thisdf<-df[df$protogene=="Emerging ORFs" , c("class","percent_tmhmm")]
thismu<-data.frame(class=c("adaptive","neutral","deleterious"),Mean = c(mean(thisdf[thisdf$class=="adaptive","percent_tmhmm"],na.rm=TRUE), mean(thisdf[thisdf$class=="neutral","percent_tmhmm"],na.rm=TRUE),mean(thisdf[thisdf$class=="deleterious","percent_tmhmm"],na.rm=TRUE) ))
plot_xlabel <-""
plot_ylabel<-""
image_protoTMtmhmm<-ggplot(thisdf, aes(x= percent_tmhmm,color= class,fill= class))+geom_density(na.rm=TRUE)+
   geom_vline(data=thismu, aes(xintercept=Mean, color= class), linetype="dashed")+scale_color_manual(values=c("darkgoldenrod1", "darkorchid3","darkgray"))+xlim(c(0,0.8))+scale_fill_manual(values=c("gold","NA", "NA"))+labs(x=plot_xlabel,y=plot_ylabel) + theme_bw()+theme(axis.text.x=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.y=element_blank(),axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")+  theme(panel.background = element_rect(fill = "lightcyan2"), panel.grid.major = element_line(colour = "white"),   panel.grid.minor = element_line(colour = "white"))

# genes using phobius
thisdf<-df[df$protogene=="Established ORFs" , c("class","percent_phobius")]
thismu<-data.frame(class=c("adaptive","neutral","deleterious"),Mean = c(mean(thisdf[thisdf$class=="adaptive","percent_phobius"],na.rm=TRUE), mean(thisdf[thisdf$class=="neutral","percent_phobius"],na.rm=TRUE),mean(thisdf[thisdf$class=="deleterious","percent_phobius"],na.rm=TRUE) ))
plot_xlabel <-"TM residues content (Phobius)"
plot_ylabel<-"Established ORFs"
image_genesTMphobius<-ggplot(thisdf, aes(x= percent_phobius,color= class,fill= class))+geom_density(na.rm=TRUE)+
   geom_vline(data=thismu, aes(xintercept=Mean, color= class), linetype="dashed")+scale_color_manual(values=c("darkgoldenrod1", "darkorchid3","darkgray"))+xlim(c(0,0.8))+scale_fill_manual(values=c("gold","NA", "NA"))+labs(x=plot_xlabel,y=plot_ylabel) + theme_bw()+theme(axis.text.x=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.y=element_blank(),axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")

# genes using tmhmm
thisdf<-df[df$protogene=="Established ORFs" , c("class","percent_tmhmm")]
thismu<-data.frame(class=c("adaptive","neutral","deleterious"),Mean = c(mean(thisdf[thisdf$class=="adaptive","percent_tmhmm"],na.rm=TRUE), mean(thisdf[thisdf$class=="neutral","percent_tmhmm"],na.rm=TRUE),mean(thisdf[thisdf$class=="deleterious","percent_tmhmm"],na.rm=TRUE) ))
plot_xlabel <-"TM residues content (TMHMM)"
plot_ylabel<-""
image_genesTMtmhmm<-ggplot(thisdf, aes(x= percent_tmhmm,color= class,fill= class))+geom_density(na.rm=TRUE)+
   geom_vline(data=thismu, aes(xintercept=Mean, color= class), linetype="dashed")+scale_color_manual(values=c("darkgoldenrod1", "darkorchid3","darkgray"))+xlim(c(0,0.8))+scale_fill_manual(values=c("gold","NA", "NA"))+labs(x=plot_xlabel,y=plot_ylabel) + theme_bw()+theme(axis.text.x=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.y=element_blank(),axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")

pdf(paste(general_outdir,"FigS6_TMResiduesContent_distributions.pdf",sep=""),width=4.5,height=5, useDingbats=FALSE)
multiplot(image_protoTMphobius, image_protoTMtmhmm, image_genesTMphobius, image_genesTMtmhmm, cols=2)
dev.off()


####################################################################################################
# Fraction of real and scrambled sequences with putative TM domain - NO Length Control - Figure 5A #
####################################################################################################
# data
df<-tm_table
df$has_tm<-ifelse(df$no_hel_Phob>0,1,0)

# calculate proportions with at least 1 TM domain
statistics_table<-NULL
statistics_table $total<- as.vector(table(df $type, df $mode))
statistics_table $yes<-as.vector(sapply(levels(df $mode), function(x) sapply(levels(df $type), function(y) sum(df[df $type ==y & df $mode==x,]$has_tm))))
statistics_table $fraction<-statistics_table $yes/statistics_table $total
statistics_table $sder<-sqrt(statistics_table $fraction * (1-statistics_table $fraction) / statistics_table $total)
statistics_table $evoclass<-rep(levels(df $type), 2)
statistics_table $seqclass<-c(rep(levels(df $mode)[1], 4),rep(levels(df $mode)[2], 4))
statistics_table <- data.frame(statistics_table)
statistics_table $evoclass <-factor(statistics_table $evoclass, levels=c("established ORFs","emerging ORFs","artificial ORFs","sORFs"))

# plot
image<-ggplot(statistics_table,aes(x= seqclass,y=fraction,ymin=(fraction - sder),ymax=(fraction + sder),fill=seqclass))+geom_bar(stat="identity",position=position_dodge(), color="black")+geom_errorbar(position = position_dodge(0.9),width=0.25)+labs(x="",y="Fraction with TM domain")+ ylim(0,0.7)+ scale_fill_manual(values=c("pink", scramble_color))  + facet_grid(.~evoclass)+theme_bw()+theme(panel.spacing = unit(0.2, "lines"), axis.text=element_text(family="Helvetica", face="bold", colour="black", size=8),  	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")
ggsave(paste(general_outdir,"Fig5A_FractionWithTM_NoLengthcontrol.pdf",sep=""),width=5,height=2.5,  useDingbats=FALSE)

######################################################################################################
# Fraction of real and scrambled sequences with putative TM domain - WITH Length Control - Figure S7A #
######################################################################################################

# data
df<-tm_table
df$has_tm<-ifelse(df$no_hel_Phob>0,1,0)

###length controlled data (all proto-genes, and other types are picked to follow proto-gene length distribution)
controlled_df<-df[df$type=="emerging ORFs",]
# determine the probability distribution for proto-genes
h<-hist(df[df$type=="emerging ORFs",]$length, breaks= seq(0,650,by=25), plot=FALSE)
protobreaks<-h$breaks
protoprob<-h$counts/nrow(df[df$type=="emerging ORFs",]) 
# Select 1000 non-genic ORFs following the proto-gene length distribution, with replacement
t<-df[df $type=="sORFs" & df $mode=="real",]
adj_protobreaks<-c(seq(25,150,by=25), 650)
adj_protoprob<-c(sum(protoprob[1:2]),protoprob[3:6], sum(protoprob[7:length(protoprob)]))
t$lengthbin<-cut(t$length, breaks=adj_protobreaks)
bindistrib<-as.data.frame(table(sample(x=levels(t$lengthbin), 1000, replace=T, prob= adj_protoprob)))
l<-sapply(bindistrib$Var1, function(x) sample(t[t$lengthbin==x,]$names, bindistrib[bindistrib$Var1==x,]$Freq, replace=T))
rownumbers<-NULL
for(x in unlist(l)){rownumbers<-c(rownumbers, which(df$names==x,))}
controlled_df <-rbind(controlled_df, df[rownumbers,])
# Select 1000 genes following the proto-gene length distribution, with replacement
t<-df[df $type=="established ORFs" & df $mode=="real",]
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
t<-df[df $type=="artificial ORFs" & df $mode=="real",]
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
statistics_table $evoclass <-factor(statistics_table $evoclass, levels=c("established ORFs","emerging ORFs","artificial ORFs","sORFs"))

# plot
image<-ggplot(statistics_table,aes(x= seqclass,y=fraction,ymin=(fraction - sder),ymax=(fraction + sder),fill=seqclass))+geom_bar(stat="identity",position=position_dodge(), color="black")+geom_errorbar(position = position_dodge(0.9),width=0.25)+labs(x="",y="Fraction with TM domain")+ ylim(0,0.7)+ scale_fill_manual(values=c("pink", scramble_color))  + facet_grid(.~evoclass)+theme_bw()+theme(panel.spacing = unit(0.2, "lines"), axis.text=element_text(family="Helvetica", face="bold", colour="black", size=8),  	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")
ggsave(paste(general_outdir,"FigS7A_FractionWithTM_Lengthcontrol.pdf",sep=""),width=5,height=2.5,  useDingbats=FALSE)

######################################################################################################################
# TM propensity in real and scrambled sequences as a function of Thymine content - NO Length Control - Figure 5B #
######################################################################################################################
# data
df<-tm_table
df$has_tm<-ifelse(df$no_hel_Phob>0,1,0)
df $type <-factor(df $type, levels=c("established ORFs","emerging ORFs","artificial ORFs","sORFs"))

# plot fraction of sequences with tm domain and density as a function of T
binsize<-0.1
allbins<-seq(0.05+binsize/2,0.75-binsize/2,binsize)
image_bars<-ggplot(df,aes(x= T, fill=mode)) +stat_summary_bin(aes(y = has_tm), fun.y="mean", bins=length(allbins),geom="bar", position=position_dodge(), color="darkgray")+stat_summary_bin(aes(y = has_tm), fun.data="mean_se", bins=length(allbins),geom="errorbar", position=position_dodge(),color="gray") + xlim(0.05, 0.8)+scale_fill_manual(values=c("pink", scramble_color))+labs(x="",y="With TM domain")  + facet_wrap(~type,nrow=1)+theme_bw()+theme(panel.spacing = unit(0.2, "lines"), axis.text.y=element_text(family="Helvetica", face="bold", colour="black", size=8), axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none", strip.text=element_blank(), axis.text.x=element_blank()) + geom_density(aes(y=..scaled..), color="darkgrey", fill="lightgrey", alpha = 0.1, size=0.4)

# plot % Tm residues content as a function of T, considering only those sequences that do have a putative tm domain
df$fractiontm<-df$no_res_TM_Phob / df$length
image_points<-ggplot(df[df$has_tm==1,],aes(x= T,y= fractiontm,colour=mode, fill=mode))+ geom_point(size=0.2, alpha=0.1) + geom_smooth(method="lm", na.rm=TRUE) + xlim(0.05, 0.8) + scale_color_manual(values=c("pink", scramble_color)) +labs(x="Thymine content",y="TM residues content")+ scale_fill_manual(values=c("pink", scramble_color))  + ylim(0,1)+facet_wrap(~type,nrow=1, strip.position="right")+theme_bw()+theme(panel.spacing = unit(0.2, "lines"), axis.text=element_text(family="Helvetica", face="bold", colour="black", size=8), axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none",strip.text=element_blank())

# composite plot
p_double<-ggarrange(image_bars,image_points,nrow=2)
ggsave(paste(general_outdir,"Fig5B_Relationship_TM_Thymine_NOLengthcontrol.pdf",sep=""),width=5,height=3.5,  useDingbats=FALSE)

##################################################################################################################
# TM propensity in real and scrambled sequences as a function of Thymine content - WITH Length Control - Figure S7B #
##################################################################################################################
# data
df<-tm_table
df$has_tm<-ifelse(df$no_hel_Phob>0,1,0)

###length controlled data (all proto-genes, and other types are picked to follow proto-gene length distribution)
controlled_df<-df[df$type=="emerging ORFs",]
# determine the probability distribution for proto-genes
h<-hist(df[df$type=="emerging ORFs",]$length, breaks= seq(0,650,by=25), plot=FALSE)
protobreaks<-h$breaks
protoprob<-h$counts/nrow(df[df$type=="emerging ORFs",]) 
# Select 1000 non-genic ORFs following the proto-gene length distribution, with replacement
t<-df[df $type=="sORFs" & df $mode=="real",]
adj_protobreaks<-c(seq(25,150,by=25), 650)
adj_protoprob<-c(sum(protoprob[1:2]),protoprob[3:6], sum(protoprob[7:length(protoprob)]))
t$lengthbin<-cut(t$length, breaks=adj_protobreaks)
bindistrib<-as.data.frame(table(sample(x=levels(t$lengthbin), 1000, replace=T, prob= adj_protoprob)))
l<-sapply(bindistrib$Var1, function(x) sample(t[t$lengthbin==x,]$names, bindistrib[bindistrib$Var1==x,]$Freq, replace=T))
rownumbers<-NULL
for(x in unlist(l)){rownumbers<-c(rownumbers, which(df$names==x,))}
controlled_df <-rbind(controlled_df, df[rownumbers,])
# Select 1000 genes following the proto-gene length distribution, with replacement
t<-df[df $type=="established ORFs" & df $mode=="real",]
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
t<-df[df $type=="artificial ORFs" & df $mode=="real",]
adj_protobreaks<-c(seq(25,300,by=25 ),650)
adj_protoprob<-c(sum(protoprob[1:2]), protoprob[3:12], sum(protoprob[13:length(protoprob)]))
t$lengthbin<-cut(t$length, breaks=adj_protobreaks)
t<-t[!(is.na(t$lengthbin)),]
bindistrib<-as.data.frame(table(sample(x=levels(t$lengthbin), 1000, replace=T, prob= adj_protoprob)))
l<-sapply(bindistrib$Var1, function(x) sample(t[t$lengthbin==x,]$names, bindistrib[bindistrib$Var1==x,]$Freq, replace=T))
rownumbers<-NULL
for(x in unlist(l)){rownumbers<-c(rownumbers, which(df$names==x,))}
controlled_df <-rbind(controlled_df, df[rownumbers,])
controlled_df $type <-factor(controlled_df $type, levels=c("established ORFs","emerging ORFs","artificial ORFs","sORFs"))

# plot fraction of sequences with tm domain and density as a function of T
binsize<-0.1
allbins<-seq(0.05+binsize/2,0.75-binsize/2,binsize)
image_bars<-ggplot(controlled_df,aes(x= T, fill=mode)) +stat_summary_bin(aes(y = has_tm), fun.y="mean", bins=length(allbins),geom="bar", position=position_dodge(), color="darkgray")+stat_summary_bin(aes(y = has_tm), fun.data="mean_se", bins=length(allbins),geom="errorbar", position=position_dodge(),color="gray") + xlim(0.05, 0.7)+scale_fill_manual(values=c("pink", scramble_color))+labs(x="",y="With TM domain")  + facet_wrap(~type,nrow=1)+theme_bw()+theme(panel.spacing = unit(0.2, "lines"), axis.text.y=element_text(family="Helvetica", face="bold", colour="black", size=8), axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none", strip.text=element_blank(), axis.text.x=element_blank()) + geom_density(aes(y=..scaled..), color="darkgrey", fill="lightgrey", alpha = 0.1, size=0.4)

# plot % Tm residues content as a function of T, considering only those sequences that do have a putative tm domain
controlled_df$fractiontm<-controlled_df$no_res_TM_Phob / controlled_df$length
image_points<-ggplot(controlled_df[controlled_df$has_tm==1,],aes(x= T,y= fractiontm,colour=mode, fill=mode))+ geom_point(size=0.2, alpha=0.1) + geom_smooth(method="lm", na.rm=TRUE) + xlim(0.05, 0.7) + scale_color_manual(values=c("pink", scramble_color)) +labs(x="Thymine content",y="TM residues content")+ scale_fill_manual(values=c("pink", scramble_color))  + ylim(0,1)+facet_wrap(~type,nrow=1, strip.position="right")+theme_bw()+theme(panel.spacing = unit(0.2, "lines"), axis.text=element_text(family="Helvetica", face="bold", colour="black", size=8), axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none",strip.text=element_blank())

# composite plot
p_double<-ggarrange(image_bars,image_points,nrow=2)
ggsave(paste(general_outdir,"FigS7B_Relationship_TM_Thymine_Lengthcontrol.pdf",sep=""),width=5,height=3.5,  useDingbats=FALSE)

#####################################################################################################################################
# TM propensity in sORFs (25 - 75 codons) as a function of their length relative to their containing intergene - Figure 5C #
#####################################################################################################################################
# data
df<-ngo_intergene_table
otherdf<-tm_table

# TM propensity to expect for proto-genes and intergenes of same length:
protovalue<-100*nrow(otherdf[otherdf $type=="emerging ORFs" & otherdf $mode=="real" & otherdf $no_hel_Phob>0 & otherdf $length>25 & otherdf $length<75,])/nrow(otherdf[otherdf $type=="emerging ORFs" & otherdf $mode=="real" & otherdf $length>25 & otherdf $length<75,])
intergenevalue<-100*nrow(otherdf[otherdf $type=="artificial ORFs" & otherdf $mode=="real" & otherdf $no_hel_Phob >0 & otherdf $length>25 & otherdf $length<75,])/nrow(otherdf[otherdf $type=="artificial ORFs" & otherdf $mode=="real" & otherdf $length>25 & otherdf $length<75,])

# summarize
t<-NULL
t$bin<-sort(unique(df$bin))
t$total<-as.vector(table(df$bin))
t$yes<-sapply(sort(unique(df$bin)), function(x) sum(ifelse(df[df$bin ==x,]$NGO_no_hel>0,1,0)))
t$fraction<-t$yes/t$total
t$sder<-sqrt(t$fraction * (1-t$fraction) / t$total)
t<-data.frame(t)

# plot
image<-ggplot(t,aes(x= bin, y=fraction *100, ymin=100*(fraction - sder),ymax=100*(fraction + sder)))+ geom_bar(stat="identity", colour="black", fill="pink") +geom_errorbar(width=0.025)+labs(x="% of intergene's length covered by sORF",y="% of sORFs \n with TM domain")+geom_hline(yintercept= protovalue,linetype="dashed",color="lightblue", size=1)+ geom_hline(yintercept= intergenevalue,linetype="dashed",color="darkblue", size=1)+theme_bw()+theme( axis.text=element_text(family="Helvetica", face="bold", colour="black", size=8),  	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")
ggsave(paste(general_outdir,"Fig5C_FractionWithTM_NGO_PercentIntergene.pdf",sep=""),width=5,height=2,  useDingbats=FALSE)
