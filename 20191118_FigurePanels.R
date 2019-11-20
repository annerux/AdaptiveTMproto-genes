# November 18 2019
# Anne-Ruxandra Carvunis and Nikolaos Vakirlis

# this script reproduces all the figure panels resulting from computational analyses in Vakirlis et al from the Source Data tables. 

# LIBRARIES
library(ggplot2)
library(plyr)
library(ggpubr)
library(effsize)
library(nnet)
library(reshape2)
library(dplyr)
library(tidyr)

# paramaters
set.seed(35689118)
scramble_color<-"wheat3" 
  
# FOLDER
general_outdir<-"/users/annerux/Workspace/cerevisiae/txt/screen_analyses/Protogenes/ProtogeneFitnessFigurePanels/resubmission/" # CHANGE THIS LINE TO INDICATE WHERE YOU STORED THE DATA FILES 1-5 

# DATA SOURCES

# description of ORFs used in analyses 
orf_table<-read.csv(paste(general_outdir, "Data1.csv",sep=""),header=TRUE)
orf_table $protogene <-factor(orf_table $emergence_status, levels=c("Emerging ORFs","Established ORFs"))
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
tm_table $type <-factor(tm_table $type, levels=c("established ORFs","emerging ORFs","iORFs","sORFs"))
ngo_intergene_table<-read.csv(paste(general_outdir, "Data5.csv",sep=""),header=TRUE)


#####################################################################
# Figure 1B: Counting ORFs with deletion and overexpression strains #
######################################################################
### total annotated ORFs:
table(orf_table$emergence_status)

 #  Emerging ORFs Established ORFs 
 #            664             5646 

### with deletion strain
df<-orf_table[!is.na(orf_table$loss_fitness) & orf_table$overlap==0,c("emergence_status","loss_fitness")]
table(df[,1])

 #  Emerging ORFs Established ORFs 
 #            239             4410 
 
### with overexpression strain
table(orf_table[orf_table$barflex_space=="yes",c("emergence_status")])

 #Emerging ORFs Established ORFs 
 #            285             4362
 
 
#### Comparing length and expression level of emerging and established ORFs
 
 wilcox.test(orf_table[orf_table$emergence_status == "Emerging ORFs", ]$Length, orf_table[orf_table$emergence_status == "Established ORFs", ]$Length)
 #W = 239600, p-value < 2.2e-16
 cliff.delta(orf_table[orf_table$emergence_status == "Emerging ORFs", ]$Length, orf_table[orf_table$emergence_status == "Established ORFs", ]$Length)
 #delta estimate: -0.8721752 (large)
 
 wilcox.test(orf_table[orf_table$emergence_status == "Emerging ORFs", ]$expr_level, orf_table[orf_table$emergence_status == "Established ORFs", ]$expr_level)
# W = 527880, p-value < 2.2e-16
 cliff.delta(orf_table[orf_table$emergence_status == "Emerging ORFs", ]$expr_level, orf_table[orf_table$emergence_status == "Established ORFs", ]$expr_level)
#delta estimate: -0.7183866 (large) 
 
### SGD annotation quality for emerging ORFs (also calculated with a more recent version of SGD annotations, see Methods)

table(orf_table[,c("emergence_status", "quality")])[1,4] / sum(table(orf_table[,c("emergence_status", "quality")])[1,])
#0.02710843
             
####################################################################
# Figure 2A: fitness cost of ORF loss (mutants in the laboratory) #
#####################################################################

### data
df<-orf_table[!is.na(orf_table$loss_fitness) & orf_table$overlap==0,c("emergence_status","loss_fitness")]

### statistical comparison: Mann Whitney
wilcoxtest<-wilcox.test(df[df $emergence_status =="Established ORFs","loss_fitness"],df[df $emergence_status =="Emerging ORFs","loss_fitness"])
pval<-wilcoxtest$p.val
#[1] 1.469534e-17
cliff.delta(df[df $emergence_status =="Established ORFs","loss_fitness"],df[df $emergence_status =="Emerging ORFs","loss_fitness"])
# delta estimate: -0.3270942 (small)

### statistical comparison: Fisher's exact with fitness cutoff of 0.9 
n_pg_effect<-nrow(df[df$emergence_status =="Emerging ORFs" & df$loss_fitness <0.9,]) 
# [1] 19
n_pg_noeffect<-nrow(df[df$emergence_status =="Emerging ORFs" & df$loss_fitness>=0.9,])
n_g_effect<-nrow(df[df$emergence_status =="Established ORFs" & df$loss_fitness <0.9,]) 
#[1] 1290
n_g_noeffect<-nrow(df[df$emergence_status =="Established ORFs" & df$loss_fitness>=0.9,])
fisher.test(matrix(c(n_pg_effect, n_pg_noeffect, n_g_effect, n_g_noeffect),2,2))
# data:  
# p-value = 3.553e-15
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
 # 0.1228600 0.3359858
# sample estimates:
# odds ratio 
 # 0.2089081 

### cumulative plot
plot_xlabel<-"Fitness upon loss"
plot_ylabel<- "ORFs (cumulative frequency)"
image<-ggplot(df, aes(x= loss_fitness,color= emergence_status, fill= emergence_status))+stat_ecdf(geom="step", na.rm=TRUE)+ stat_ecdf(aes(ymin=0,ymax=..y..),
geom="ribbon", na.rm=TRUE)+scale_color_manual(values=c("deepskyblue4","black"))+ scale_fill_manual(values=c("lightblue","NA"))+labs(x=plot_xlabel,y=plot_ylabel)+ geom_vline(xintercept=0.9, colour="red") + theme_bw()+theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8),  	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none") 
ggsave(paste(general_outdir,"Fig2A_Compare_FitnessUponLoss_CumulativeFreq.pdf",sep=""),width=2.5,height=2.5)

##############################################################################################################
# Figure S1A # fitness cost of ORF loss (mutants in the laboratory)- Control for length and expression level #
##############################################################################################################

# data
df<-orf_table[!is.na(orf_table$loss_fitness) & orf_table$overlap==0,c("orf_name","emergence_status","loss_fitness", "expr_level", "Length")]

### controlled data (all emerging ORFs, all established ORFs and randomly picked estabslished ORFs that follow the length and expr_level distribution of emerging ORFs)
controlled_df<-df
n_toselect<-nrow(df[df $emergence_status =="Established ORFs" ,]) # number to pick randomly for controlling

# EXPR_LEVEL
h<-hist(df[df$emergence_status =="Emerging ORFs",]$expr_level, breaks= seq(-7,6,by=1), plot=FALSE)
protobreaks<-h$breaks
protoprob<-h$counts/nrow(df[df$emergence_status =="Emerging ORFs",]) 
# Select random genes following the proto-gene length distribution, with replacement
set.seed(35689118)
t<-df[df $emergence_status =="Established ORFs" ,]
adj_protobreaks<-c((-7), seq(-5,3,by=1), 6)
adj_protoprob<-c(sum(protoprob[1:2]), protoprob[3:10], sum(protoprob[11:length(protoprob)]))
t$expbin<-cut(t$expr_level, breaks=adj_protobreaks)
t<-t[!(is.na(t$expbin)),]
bindistrib<-as.data.frame(table(sample(x=levels(t$expbin), n_toselect, replace=T, prob= adj_protoprob)))
l<-sapply(bindistrib$Var1, function(x) sample(t[t$expbin==x,]$orf_name, bindistrib[bindistrib$Var1==x,]$Freq, replace=T))
rownumbers<-NULL
for(x in unlist(l)){rownumbers<-c(rownumbers, which(df$orf_name==x,))}
expcontrolled_df<-df[rownumbers,]
expcontrolled_df$emergence_status <-rep("esta_expcontrolled", n_toselect)
controlled_df <-rbind(controlled_df, expcontrolled_df)

#LENGTH
# determine the probability distribution for proto-genes
h<-hist(df[df$emergence_status =="Emerging ORFs",]$Length, breaks= seq(50,2000,by=50), plot=FALSE)
protobreaks<-h$breaks
protoprob<-h$counts/nrow(df[df$emergence_status =="Emerging ORFs",]) 
# Select random genes following the proto-gene length distribution, with replacement
set.seed(35689118)
t<-df[df $emergence_status =="Established ORFs" ,]
adj_protobreaks<-c(seq(50,800,by=50), 2000)
adj_protoprob<-c(protoprob[1:15], sum(protoprob[16:length(protoprob)]))
t$lengthbin<-cut(t$Length, breaks=adj_protobreaks)
t<-t[!(is.na(t$lengthbin)),]
bindistrib<-as.data.frame(table(sample(x=levels(t$lengthbin), n_toselect, replace=T, prob= adj_protoprob)))
l<-sapply(bindistrib$Var1, function(x) sample(t[t$lengthbin==x,]$orf_name, bindistrib[bindistrib$Var1==x,]$Freq, replace=T))
rownumbers<-NULL
for(x in unlist(l)){rownumbers<-c(rownumbers, which(df$orf_name==x,))}
lengthcontrolled_df<-df[rownumbers,]
lengthcontrolled_df$emergence_status <-rep("esta_lengthcontrolled", n_toselect)
controlled_df <-rbind(controlled_df, lengthcontrolled_df)

# plot fractions with fitness costs
cutoff<-0.9
stat_df<-controlled_df[,c("emergence_status","loss_fitness")]
stat_df$is_delet<-ifelse(stat_df$loss_fitness<cutoff,1,0)
plotdf<-NULL
plotdf $total<- as.vector(table(stat_df $emergence_status))
plotdf $yes<-as.vector(sapply(levels(stat_df $emergence_status), function(x) sum(stat_df[stat_df $emergence_status ==x,]$is_delet)))
plotdf $fraction<-plotdf $yes/plotdf $total
plotdf $sder<-sqrt(plotdf $fraction * (1-plotdf $fraction) / plotdf $total)
plotdf $class<-levels(stat_df $emergence_status)
plotdf <- data.frame(plotdf)
plotdf$class<-factor(plotdf$class, levels = c("Emerging ORFs","esta_expcontrolled", "esta_lengthcontrolled", "Established ORFs"))
levels(plotdf$class)<-mapvalues(levels(plotdf$class),from=levels(plotdf$class), to= c("Emerging ORFs", "Established ORFs/matched expression", "Established ORFs/matched length", "All established ORFs"))

image<-ggplot(plotdf,aes(x= class,y=fraction,ymin=(fraction - sder),ymax=(fraction + sder),fill=class, color=class))+geom_bar(stat="identity",position=position_dodge())+geom_errorbar(position = position_dodge(0.9),width=0.25)+coord_flip()+ labs(x="",y="Fraction with Fitness Cost")+ scale_color_manual(values=c("deepskyblue4", "black", "black","black")) + scale_fill_manual(values=c("lightblue", "white", "white","white")) +theme_bw()+theme(panel.spacing = unit(0, "lines"), axis.text.x=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.y=element_text(family="Helvetica", face="bold", colour="black", size=7), 	axis.title.x=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.y=element_blank(),	legend.position="none") 
ggsave(paste(general_outdir,"FigS1A_Compare_FitnessUponLoss_controlled_barplots.pdf",sep=""),width=4.5,height=2.5,  useDingbats=FALSE)

# statistical comparison: Fisher's exact with fitness cutoff ; between emerging ORFs and expression-controled established ORFs
n_pg_effect<-nrow(controlled_df[controlled_df $emergence_status =="Emerging ORFs" & controlled_df $loss_fitness <cutoff,]) 
n_pg_noeffect<-nrow(controlled_df[controlled_df $emergence_status =="Emerging ORFs" & controlled_df $loss_fitness>=cutoff,])
n_g_effect<-nrow(controlled_df[controlled_df $emergence_status =="esta_expcontrolled" & controlled_df $loss_fitness <cutoff,]) 
n_g_noeffect<-nrow(controlled_df[controlled_df $emergence_status =="esta_expcontrolled" & controlled_df $loss_fitness>=cutoff,])
fisher.test(matrix(c(n_pg_effect, n_pg_noeffect, n_g_effect, n_g_noeffect),2,2))

# data:  
# p-value = 0.05164
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
 # 0.3624975 1.0004706
# sample estimates:
# odds ratio 
 # 0.6189953

# statistical comparison: Fisher's exact with fitness cutoff ; between emerging ORFs and length-controled established ORFs
n_pg_effect<-nrow(controlled_df[controlled_df $emergence_status =="Emerging ORFs" & controlled_df $loss_fitness <cutoff,]) 
n_pg_noeffect<-nrow(controlled_df[controlled_df $emergence_status =="Emerging ORFs" & controlled_df $loss_fitness>=cutoff,])
n_g_effect<-nrow(controlled_df[controlled_df $emergence_status =="esta_lengthcontrolled" & controlled_df $loss_fitness <cutoff,]) 
n_g_noeffect<-nrow(controlled_df[controlled_df $emergence_status =="esta_lengthcontrolled" & controlled_df $loss_fitness>=cutoff,])
fisher.test(matrix(c(n_pg_effect, n_pg_noeffect, n_g_effect, n_g_noeffect),2,2))

# data:  
# p-value = 2.946e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
 # 0.1172831 0.3206336
# sample estimates:
# odds ratio 
 # 0.1993821 

#########################################################
# Figure 2B: fixation of ORf sctructure across isolates #
#########################################################

df<-orf_table[!(is.na(orf_table $strain_goodalign)),c("emergence_status","strain_goodalign", "strained_conserved_orf",  "strain_pi")]

# plot
plot_xlabel<-"% isolates with intact ORF"
plot_ylabel<- "ORFs (cumulative frequency)"
image<-ggplot(df, aes(x= strained_conserved_orf/2022*100,color= emergence_status, fill= emergence_status))+stat_ecdf(geom="step", na.rm=TRUE)+ stat_ecdf(aes(ymin=0,ymax=..y..),
geom="ribbon", na.rm=TRUE)+scale_color_manual(values=c("deepskyblue4","black"))+ scale_fill_manual(values=c("lightblue","NA"))+labs(x=plot_xlabel,y=plot_ylabel) + geom_vline(xintercept=90, colour="red")+ theme_bw()+theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8),  	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none") + coord_cartesian(xlim=c(0,100), , expand=FALSE)
ggsave(paste(general_outdir,"Fig2B_Compare_StrainORFStructure_CumulativeFreq.pdf",sep=""),width=2.5,height=2.5)

# intact orf structure in <90% isolates
df$strict<-ifelse(df$strained_conserved_orf/2022*100<90,1,0)
table(df[,c("emergence_status","strict")])

# emergence_status      0    1
  # Emerging ORFs     274  390
  # Established ORFs 4502 1141
  
table(df[,c("emergence_status","strict")])[2,1]/sum(table(df[,c("emergence_status","strict")])[2,])
#[1] 0.7978026
table(df[,c("emergence_status","strict")])[1,1]/sum(table(df[,c("emergence_status","strict")])[1,])
#[1] 0.4126506

# statistical comparison: Fisher's exact with ORF structure cutoff ; between emerging ORFs and established ORFs
fisher.test(table(df[,c("emergence_status","strict")]))

# data:  table(df[, c("emergence_status", "strict")])
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
 # 0.1500218 0.2112741
# sample estimates:
# odds ratio 
 # 0.1781234 
1/fisher.test(table(df[,c("emergence_status","strict")]))$estimate
#odds ratio 
#  5.614086 

################################################################################
# Figure S1B Fixation of orf structure - Control for length and expression level #
################################################################################
# data
df<-orf_table[!(is.na(orf_table $strain_goodalign)),c("orf_name","emergence_status","strain_goodalign", "strained_conserved_orf", "expr_level", "Length")]

### controlled data (all emerging ORFs, all established ORFs and randomly picked established ORFs that follow the length and expr_level distribution of emerging ORFs)
controlled_df<-df
n_toselect<-nrow(df[df $emergence_status =="Established ORFs" ,])

#LENGTH
# determine the probability distribution for proto-genes
h<-hist(df[df$emergence_status =="Emerging ORFs",]$Length, breaks= seq(50,2000,by=50), plot=FALSE)
protobreaks<-h$breaks
protoprob<-h$counts/nrow(df[df$emergence_status =="Emerging ORFs",]) 
# Select genes following the proto-gene length distribution, with replacement
set.seed(35689118)
t<-df[df $emergence_status =="Established ORFs" ,]
adj_protobreaks<-c(seq(50,800,by=50), 2000)
adj_protoprob<-c(protoprob[1:15], sum(protoprob[16:length(protoprob)]))
t$lengthbin<-cut(t$Length, breaks=adj_protobreaks)
t<-t[!(is.na(t$lengthbin)),]
bindistrib<-as.data.frame(table(sample(x=levels(t$lengthbin), n_toselect, replace=T, prob= adj_protoprob)))
l<-sapply(bindistrib$Var1, function(x) sample(t[t$lengthbin==x,]$orf_name, bindistrib[bindistrib$Var1==x,]$Freq, replace=T))
rownumbers<-NULL
for(x in unlist(l)){rownumbers<-c(rownumbers, which(df$orf_name==x,))}
lengthcontrolled_df<-df[rownumbers,]
lengthcontrolled_df$emergence_status <-rep("esta_lengthcontrolled", n_toselect)
controlled_df <-rbind(controlled_df, lengthcontrolled_df)

# EXPR_LEVEL
h<-hist(df[df$emergence_status =="Emerging ORFs",]$expr_level, breaks= seq(-7,6,by=1), plot=FALSE)
protobreaks<-h$breaks
protoprob<-h$counts/nrow(df[df$emergence_status =="Emerging ORFs",]) 
# Select genes following the proto-gene length distribution, with replacement
set.seed(35689118)
t<-df[df $emergence_status =="Established ORFs" ,]
adj_protobreaks<-c((-7), seq(-5,3,by=1), 6)
adj_protoprob<-c(sum(protoprob[1:2]), protoprob[3:10], sum(protoprob[11:length(protoprob)]))
t$expbin<-cut(t$expr_level, breaks=adj_protobreaks)
t<-t[!(is.na(t$expbin)),]
bindistrib<-as.data.frame(table(sample(x=levels(t$expbin), n_toselect, replace=T, prob= adj_protoprob)))
l<-sapply(bindistrib$Var1, function(x) sample(t[t$expbin==x,]$orf_name, bindistrib[bindistrib$Var1==x,]$Freq, replace=T))
rownumbers<-NULL
for(x in unlist(l)){rownumbers<-c(rownumbers, which(df$orf_name==x,))}
expcontrolled_df<-df[rownumbers,]
expcontrolled_df$emergence_status <-rep("esta_expcontrolled", n_toselect)
controlled_df <-rbind(controlled_df, expcontrolled_df)

# plot fractions with fixed structures
cutoff<-0.9
stat_df<-controlled_df[,c("emergence_status","strain_goodalign", "strained_conserved_orf")]
stat_df$is_delet<-ifelse((stat_df$strained_conserved_orf/2022)>cutoff,1,0)
plotdf<-NULL
plotdf $total<- as.vector(table(stat_df $emergence_status))
plotdf $yes<-as.vector(sapply(levels(stat_df $emergence_status), function(x) sum(stat_df[stat_df $emergence_status ==x,]$is_delet)))
plotdf $fraction<-plotdf $yes/plotdf $total
plotdf $sder<-sqrt(plotdf $fraction * (1-plotdf $fraction) / plotdf $total)
plotdf $class<-levels(stat_df $emergence_status)
plotdf <- data.frame(plotdf)
plotdf$class<-factor(plotdf$class, levels = c("Emerging ORFs","esta_expcontrolled", "esta_lengthcontrolled", "Established ORFs"))
levels(plotdf$class)<-mapvalues(levels(plotdf$class),from=levels(plotdf$class), to= c("Emerging ORFs", "Established ORFs/matched expression", "Established ORFs/matched length", "All established ORFs"))

# plot the fraction with fitness costs
image<-ggplot(plotdf,aes(x= class,y=fraction,ymin=(fraction - sder),ymax=(fraction + sder),fill=class, color=class))+geom_bar(stat="identity",position=position_dodge())+geom_errorbar(position = position_dodge(0.9),width=0.25)+coord_flip()+ labs(x="",y="Fraction with fixed ORF structure")+ scale_color_manual(values=c("deepskyblue4", "black", "black","black")) + scale_fill_manual(values=c("lightblue", "white", "white","white")) +theme_bw()+theme(panel.spacing = unit(0, "lines"), axis.text.x=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.y=element_text(family="Helvetica", face="bold", colour="black", size=7), 	axis.title.x=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.y=element_blank(),	legend.position="none") 
ggsave(paste(general_outdir,"FigS1B_Compare_StrainORFStructure_controlled_barplots.pdf",sep=""),width=4.5,height=2.5,  useDingbats=FALSE)

# statistical comparison: Fisher's exact with fitness cutoff ; between emerging ORFs and expr_controled established ORFs
n_pg_effect<-nrow(controlled_df[controlled_df $emergence_status =="Emerging ORFs" & controlled_df $strained_conserved_orf/2022 <cutoff,]) 
n_pg_noeffect<-nrow(controlled_df[controlled_df $emergence_status =="Emerging ORFs" & controlled_df $strained_conserved_orf/2022>=cutoff,])
n_g_effect<-nrow(controlled_df[controlled_df $emergence_status =="esta_expcontrolled" & controlled_df $strained_conserved_orf/2022 <cutoff,]) 
n_g_noeffect<-nrow(controlled_df[controlled_df $emergence_status =="esta_expcontrolled" & controlled_df $strained_conserved_orf/2022>=cutoff,])
fisher.test(matrix(c(n_pg_effect, n_pg_noeffect, n_g_effect, n_g_noeffect),2,2))

# data:  
# p-value = 3.059e-13
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
 # 1.550602 2.163933
# sample estimates:
# odds ratio 
  # 1.830887 

# statistical comparison: Fisher's exact with fitness cutoff ; between emerging ORFs and length_controled established ORFs
n_pg_effect<-nrow(controlled_df[controlled_df $emergence_status =="Emerging ORFs" & controlled_df $strained_conserved_orf/2022 <cutoff,]) 
n_pg_noeffect<-nrow(controlled_df[controlled_df $emergence_status =="Emerging ORFs" & controlled_df $strained_conserved_orf/2022>=cutoff,])
n_g_effect<-nrow(controlled_df[controlled_df $emergence_status =="esta_lengthcontrolled" & controlled_df $strained_conserved_orf/2022 <cutoff,]) 
n_g_noeffect<-nrow(controlled_df[controlled_df $emergence_status =="esta_lengthcontrolled" & controlled_df $strained_conserved_orf/2022>=cutoff,])
fisher.test(matrix(c(n_pg_effect, n_pg_noeffect, n_g_effect, n_g_noeffect),2,2))

# data:  
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
 # 4.228200 5.944408
# sample estimates:
# odds ratio 
  # 5.010927 

  
###################################################
# Figure 2C: nucleotide diversity across isolates #
###################################################

df<-orf_table[!(is.na(orf_table $strain_goodalign)),c("emergence_status","strain_goodalign", "strained_conserved_orf",  "strain_pi")]

# statistical comaprison: Mann Whitney
wilcoxtest<-wilcox.test(df[df $emergence_status =="Established ORFs","strain_pi"],df[df $emergence_status =="Emerging ORFs","strain_pi"])
wilcoxtest$p.val
#[1] 6.409581e-46
cliff.delta(df[df $emergence_status =="Established ORFs","strain_pi"],df[df $emergence_status =="Emerging ORFs","strain_pi"])
#delta estimate: -0.336976 (medium)

# plot
mu <- ddply(df, "emergence_status", summarise, grp.mean=mean(strain_pi,na.rm=TRUE))
plot_xlabel<-"Nucleotide diversity"
plot_ylabel<- "ORFs (density)"

image<-ggplot(df, aes(x= strain_pi,color= emergence_status,fill= emergence_status))+geom_density(na.rm=TRUE)+
   geom_vline(data=mu, aes(xintercept=grp.mean, color= emergence_status), linetype="dashed")+scale_color_manual(values=c("deepskyblue4","black"))+scale_fill_manual(values=c( "lightblue","NA"))+labs(x=plot_xlabel,y=plot_ylabel)+ theme_bw()+theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8),  	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")
ggsave(paste(general_outdir,"Fig2C_Compare_StrainNucleotideDiversity.pdf",sep=""),width=2.5,height=2.5)

# ###########################################################
# # Figure 3A: number of colonies per category overexpressed#
# ###########################################################
# Each overexpressed ORFs is represented by multiple colonies in the arrays
df<-colony_table

# category for each colony
df$category<-rep("other", nrow(df))
df[df$orf_name=="BF_control",]$category<-"reference"
df[df$orf_name %in% orf_table[orf_table$emergence_status =="Emerging ORFs",]$orf_name,]$category<-"Emerging ORFs"
df[df$orf_name %in% orf_table[orf_table$emergence_status =="Established ORFs",]$orf_name,]$category<-"Established ORFs"
df $category <-factor(df $category, levels=c("Emerging ORFs","Established ORFs", "reference"))
table(df$category)

   # Emerging ORFs Established ORFs        reference 
            # 1157            17078             2301 

############################################################################
# Figure 3B: fraction of orfs with relative overexpression fitness effects #
############################################################################

df<-orf_table[orf_table$barflex_space=="yes",c("emergence_status","overexpression_relative_fitness")]
table(df[,1])

# calculate proportions and standard error of proportions
a<-data.frame(table(df))
totals<-rep(table(df$emergence_status),3)
a$groupsize<-totals
sder<-sqrt(a$Freq/a$groupsize*(1-a$Freq/a$groupsize)/a$groupsize)
a$err<-sder
max<-a$Freq/a$groupsize + a$err
a$top<-max
min<-a$Freq/a$groupsize - a$err
 a$bottom<-min
 proportion<-a$Freq/a$groupsize
 a$fraction<-proportion
a
  # emergence_status overexpression_relative_fitness Freq groupsize         err        top     bottom   fraction
# 1    Emerging ORFs                       increased   14       285 0.012802105 0.06192491 0.03632070 0.04912281
# 2 Established ORFs                       increased   49      4362 0.001595730 0.01282911 0.00963765 0.01123338
# 3    Emerging ORFs                       decreased   55       285 0.023376418 0.21635887 0.16960604 0.19298246
# 4 Established ORFs                       decreased 1869      4362 0.007492682 0.43596586 0.42098050 0.42847318
# 5    Emerging ORFs                       unchanged  216       285 0.025373719 0.78326846 0.73252102 0.75789474
# 6 Established ORFs                       unchanged 2444      4362 0.007515302 0.56780874 0.55277814 0.56029344

# Plot 
plot_xlabel<-"Relative fitness"
plot_ylabel<- "Fraction of ORFs"  
image<-ggplot(a,aes(x= overexpression_relative_fitness,y=fraction,ymin=bottom,ymax=top,fill= emergence_status, color= emergence_status))+geom_bar(stat="identity",position=position_dodge())+geom_errorbar(position = position_dodge(0.9),width=0.25)+scale_color_manual(values=c("deepskyblue4","black"))+scale_fill_manual(values=c( "lightblue","white"))+labs(x=plot_xlabel,y=plot_ylabel)+ theme_bw()+theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8),  	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none") 
ggsave(paste(general_outdir,"Fig3B_ExperimentalResults_RelativeFitness_Fraction.pdf",sep=""),width=3.5,height=2.5)

# range of effect sizes for proto-genes (average of replicates)
increased_colonies<-colony_table[colony_table$orf_name %in% orf_table[orf_table $overexpression_relative_fitness=="increased" & orf_table$protogene=="Emerging ORFs",]$orf_name,]
mean_colonysizes<-ddply(increased_colonies, "orf_name", summarise, grp.mean=mean(colony_size,na.rm=TRUE))
min(mean_colonysizes $grp.mean)
#[1] 1.079289
max(mean_colonysizes $grp.mean)
#[1] 1.190052

##############################################################################
# Figure 3C: Odds Ratio of orfs with relative overexpression fitness effects #
##############################################################################

df<-orf_table[orf_table$barflex_space=="yes",c("emergence_status","overexpression_relative_fitness")]

# calculate odds ratios and confidence intervals and pvalues
a<-data.frame(table(df))
totals<-rep(table(df$emergence_status),3)
a$groupsize<-totals
a$minus<-a$groupsize-a$Freq
benef<-fisher.test(a[a$overexpression_relative_fitness =="increased",c("Freq","minus")])
neutral<-fisher.test(a[a$overexpression_relative_fitness =="unchanged",c("Freq","minus")])
delet<-fisher.test(a[a$overexpression_relative_fitness =="decreased",c("Freq","minus")])
oddsratio<-data.frame(overexpression_relative_fitness =c("increased","unchanged","decreased"),or=c(benef$estimate,neutral$estimate,delet$estimate))
oddsratio$top<-c(benef$conf.int[2],neutral$conf.int[2],delet$conf.int[2])
oddsratio$bottom<-c(benef$conf.int[1],neutral$conf.int[1],delet$conf.int[1])
oddsratio $overexpression_relative_fitness <-factor(oddsratio $overexpression_relative_fitness, levels=c("increased","decreased","unchanged"))

oddsratio
  # overexpression_relative_fitness        or       top    bottom
# 1                       increased 4.5444007 8.4923983 2.2871886
# 2                       unchanged 2.4561852 3.2918093 1.8517137
# 3                       decreased 0.3190501 0.4327424 0.2318243

# inverse to calculate less likely:
1/oddsratio[3,2]
#[1] 3.134304

# plot
plot_xlabel<-"Relative fitness"
plot_ylabel<- "Odds ratio (log)" 

image<-ggplot(oddsratio,aes(x= overexpression_relative_fitness,y=or,ymin=bottom,ymax=top))+geom_point(stat="identity",  shape=21, size=2, stroke=1, fill = "white")+geom_errorbar(width=0.25)+geom_hline(yintercept=1,linetype="dashed",color="red")+ labs(x=plot_xlabel,y=plot_ylabel)+ theme_bw()+theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8),  	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")+scale_y_continuous(trans="log10", breaks=c(0.2,0.5,1,2,5,10))+annotation_logticks(sides="l")
ggsave(paste(general_outdir,"Fig3C_ExperimentalResults_RelativeFitness_OddsRatio_logscale.pdf",sep=""),width=3.5,height=2.5,  useDingbats=FALSE)

######################################################
# Figure S3A: competitive fitness upon overexpression #
######################################################
# data
df<-orf_table[orf_table$barflex_space=="yes",c("emergence_status","overexpression_competitive_fitness")]

#  statistical comparison: Mann Whitney
wilcoxtest<-wilcox.test(df[df $emergence_status =="Established ORFs","overexpression_competitive_fitness"],df[df $emergence_status =="Emerging ORFs","overexpression_competitive_fitness"])
wilcoxtest$p.val
#[1] 5.51811e-32

# plot
mu <- ddply(df, "emergence_status", summarise, grp.mean=mean(overexpression_competitive_fitness,na.rm=TRUE))
plot_xlabel<-"Competitive fitness upon overexpression"
plot_ylabel<- "ORFs (density)"
image<-ggplot(df, aes(x= overexpression_competitive_fitness,color= emergence_status,fill= emergence_status))+geom_density(na.rm=TRUE)+   geom_vline(data=mu, aes(xintercept=grp.mean, color= emergence_status), linetype="dashed")+scale_color_manual(values=c("deepskyblue4", "black"))+scale_fill_manual(values=c("lightblue", "NA"))+labs(x=plot_xlabel,y=plot_ylabel) + theme_bw()+theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8),  	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")   
ggsave(paste(general_outdir,"FigS3A_Compare_CompetitiveFitness.pdf",sep=""),width=5,height=2.5)

############################################################################
# Figure S3B: showing the exact colony sizes of increased propo-genes #
############################################################################
df<-colony_table
increased_protogenes <- as.character(orf_table[which(orf_table$emergence_status =="Emerging ORFs" & orf_table$overexpression_relative_fitness=="increased"), "orf_name"])

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



###############################################################################################
# Figure 3D: Odds Ratio of orfs with relative overexpression fitness effects - 5 environments #
###############################################################################################

df<-merge(fitness_table, orf_table[,c("orf_name","emergence_status")], all.x=TRUE)

# calculate odds ratios and confidence intervals and pvalues
a<-data.frame(table(df[,c("exp_environment","emergence_status", "effect_cs")]))
totals<-table(df[,c("exp_environment","emergence_status")])
b<-merge(a,totals,by=c("exp_environment","emergence_status"),all.x=TRUE)
colnames(b)<-c("exp_environment",  "emergence_status",   "effect_cs", "number", "total")
b$minus<-b$total - b$number
protos<-b[b$emergence_status =="Emerging ORFs",c("exp_environment" , "emergence_status",   "effect_cs", "number", "minus")]
genes<-b[b$emergence_status =="Established ORFs",c("exp_environment" , "emergence_status",   "effect_cs", "number", "minus")]
numbers<-merge(protos,genes,by=c("exp_environment","effect_cs"))
colnames(numbers)<-c("exp_environment" ,   "effect_cs", "proto", "protonumber", "protominus","gene","genenumber","geneminus")
fisher_test_function<-function(x){
	mat<-matrix(as.numeric(c(x[4],x[5],x[7],x[8])),2,2)
	f<-fisher.test(mat)
	return(c(f$estimate, f$conf.int[1],f$conf.int[2], f$p.value))
}
or<-apply(numbers,1,function(x) fisher_test_function(x) )
oddsratios<-data.frame(cbind(numbers[,c("exp_environment","effect_cs")],t(or)))

 oddsratios
   # exp_environment effect_cs        or    bottom        top         pval
# 1           N-,C++ decreased 0.3422379 0.2371567  0.4828383 8.756430e-12
# 2           N-,C++ increased 5.3739590 2.9230355  9.4643426 1.790211e-07
# 3           N-,C++ unchanged 1.9748640 1.4632388  2.7013939 2.374424e-06
# 4           N+, C+ decreased 0.3190501 0.2318243  0.4327424 3.833604e-16
# 5           N+, C+ increased 4.5444007 2.2871886  8.4923983 1.886453e-05
# 6           N+, C+ unchanged 2.4561852 1.8517137  3.2918093 2.273258e-11
# 7          N+, C++ decreased 0.2829955 0.2005648  0.3918857 6.529744e-18
# 8          N+, C++ increased 5.8737662 2.9028051 11.2553585 1.650614e-06
# 9          N+, C++ unchanged 2.6441065 1.9670785  3.5999596 2.376780e-12
# 10          N++,C+ decreased 0.3541566 0.2443548  0.5014300 5.632769e-11
# 11          N++,C+ increased 3.0583496 1.3024172  6.4077607 5.655749e-03
# 12          N++,C+ unchanged 2.3211156 1.6837010  3.2567656 2.221283e-08
# 13         N++,C++ decreased 0.2784128 0.1937974  0.3913094 3.586316e-17
# 14         N++,C++ increased 5.2780961 2.4712347 10.5247646 2.128507e-05
# 15         N++,C++ unchanged 2.7361068 2.0087340  3.7848113 2.828090e-12

colnames(oddsratios)<-c("exp_environment",   "effect_cs" , "or", "bottom","top","pval")
oddsratios $effect_cs <-factor(oddsratios $effect_cs, levels=c("increased","decreased","unchanged"))

# plot
forplotting<-oddsratios[oddsratios$effect_cs %in% c("increased","decreased"),]
plot_xlabel<-"Environmental condition"
plot_ylabel<- "Odds ratio (log)" 
image<-ggplot(forplotting,aes(x= exp_environment,y=or,ymin=bottom,ymax=top))+geom_point(stat="identity", shape=21, size=2, stroke=1, fill = "white") + facet_wrap (~ effect_cs, nrow=2)+geom_errorbar(width=0.25)+geom_hline(yintercept=1,linetype="dashed",color="red")+ labs(x=plot_xlabel,y=plot_ylabel)+ theme_bw()+theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8),  	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")+scale_y_continuous(trans="log10", breaks=c(0.2,0.5,1,2,5,10))+annotation_logticks(sides="l")
ggsave(paste(general_outdir,"Fig3D_ExperimentalResults_RelativeFitness_OddsRatio_5Environments_log.pdf",sep=""),width=3.5,height=2.5,  useDingbats=FALSE)

####################################################################
# Figure 3E: Adaptive ORFs - length and expression level controlled#
####################################################################

# make a data frame with adaptive in a column
adaptive_orfs<-unique(fitness_table[fitness_table $effect_cs=="increased",]$orf_name)
df<-orf_table[orf_table$barflex_space=="yes",c("orf_name","emergence_status","Length","expr_level","loss_fitness", "strained_conserved_orf", "strain_pi")]
df$adaptive<-sapply(df$orf_name, function(x) ifelse(x %in% adaptive_orfs,1,0))
df$sick<- ifelse(df $loss_fitness<0.9,1,0)
df$intact<-ifelse(df $strained_conserved_orf <0.9 * 2022 ,0,1)

# enrichment in adaptive effects
table(df[,c("emergence_status","adaptive")])

# # emergence_status      0    1
  # Emerging ORFs     257   28
  # Established ORFs 4236  126
  
fisher.test(table(df[,c("emergence_status","adaptive")]))

# data:  table(df[, c("emergence_status", "adaptive")])
# p-value = 1.199e-07
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
 # 0.1762793 0.4357862
# sample estimates:
# odds ratio 
 # 0.2731654 


### controlled data (all proto-genes, all established orfs and randomly picked ones that follow emerging orf length and expr_level distribution)
controlled_df<-df
n_toselect<-nrow(df[df $emergence_status =="Established ORFs" ,])

#LENGTH
# determine the probability distribution for proto-genes
h<-hist(df[df$emergence_status =="Emerging ORFs",]$Length, breaks= seq(50,2000,by=50), plot=FALSE)
protobreaks<-h$breaks
protoprob<-h$counts/nrow(df[df$emergence_status =="Emerging ORFs",]) 
# Select genes following the proto-gene length distribution, with replacement
set.seed(35689118)
t<-df[df $emergence_status =="Established ORFs" ,]
adj_protobreaks<-c(seq(50,800,by=50), 2000)
adj_protoprob<-c(protoprob[1:15], sum(protoprob[16:length(protoprob)]))
t$lengthbin<-cut(t$Length, breaks=adj_protobreaks)
t<-t[!(is.na(t$lengthbin)),]
bindistrib<-as.data.frame(table(sample(x=levels(t$lengthbin), n_toselect, replace=T, prob= adj_protoprob)))
l<-sapply(bindistrib$Var1, function(x) sample(t[t$lengthbin==x,]$orf_name, bindistrib[bindistrib$Var1==x,]$Freq, replace=T))
rownumbers<-NULL
for(x in unlist(l)){rownumbers<-c(rownumbers, which(df$orf_name==x,))}
lengthcontrolled_df<-df[rownumbers,]
lengthcontrolled_df$emergence_status <-rep("esta_lengthcontrolled", n_toselect)
controlled_df <-rbind(controlled_df, lengthcontrolled_df)

# EXPR_LEVEL
h<-hist(df[df$emergence_status =="Emerging ORFs",]$expr_level, breaks= seq(-7,6,by=1), plot=FALSE)
protobreaks<-h$breaks
protoprob<-h$counts/nrow(df[df$emergence_status =="Emerging ORFs",]) 
# Select genes following the proto-gene length distribution, with replacement
set.seed(35689118)
t<-df[df $emergence_status =="Established ORFs" ,]
adj_protobreaks<-c((-7), seq(-5,3,by=1), 6)
adj_protoprob<-c(sum(protoprob[1:2]), protoprob[3:10], sum(protoprob[11:length(protoprob)]))
t$expbin<-cut(t$expr_level, breaks=adj_protobreaks)
t<-t[!(is.na(t$expbin)),]
bindistrib<-as.data.frame(table(sample(x=levels(t$expbin), n_toselect, replace=T, prob= adj_protoprob)))
l<-sapply(bindistrib$Var1, function(x) sample(t[t$expbin==x,]$orf_name, bindistrib[bindistrib$Var1==x,]$Freq, replace=T))
rownumbers<-NULL
for(x in unlist(l)){rownumbers<-c(rownumbers, which(df$orf_name==x,))}
expcontrolled_df<-df[rownumbers,]
expcontrolled_df$emergence_status <-rep("esta_expcontrolled", n_toselect)
controlled_df <-rbind(controlled_df, expcontrolled_df)

stat_df<-controlled_df[,c("emergence_status","adaptive")]
plotdf<-NULL
plotdf $total<- as.vector(table(stat_df $emergence_status))
plotdf $yes<-as.vector(sapply(levels(stat_df $emergence_status), function(x) sum(stat_df[stat_df $emergence_status ==x,]$adaptive)))
plotdf $fraction<-plotdf $yes/plotdf $total
plotdf $sder<-sqrt(plotdf $fraction * (1-plotdf $fraction) / plotdf $total)
plotdf $class<-levels(stat_df $emergence_status)
plotdf <- data.frame(plotdf)
plotdf$class<-factor(plotdf$class, levels = c("Emerging ORFs","esta_expcontrolled", "esta_lengthcontrolled", "Established ORFs"))
levels(plotdf$class)<-mapvalues(levels(plotdf$class),from=levels(plotdf$class), to= c("Emerging ORFs", "Established ORFs/matched expression", "Established ORFs/matched length", "All established ORFs"))

# plot the fraction adaptive
image<-ggplot(plotdf,aes(x= class,y=fraction,ymin=(fraction - sder),ymax=(fraction + sder),fill=class, color=class))+geom_bar(stat="identity",position=position_dodge())+geom_errorbar(position = position_dodge(0.9),width=0.25)+coord_flip()+ labs(x="",y="Fraction adaptive")+ scale_color_manual(values=c("deepskyblue4", "black", "black","black")) + scale_fill_manual(values=c("lightblue", "white", "white","white")) +theme_bw()+theme(panel.spacing = unit(0, "lines"), axis.text.x=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.y=element_text(family="Helvetica", face="bold", colour="black", size=7), 	axis.title.x=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.y=element_blank(),	legend.position="none") 
ggsave(paste(general_outdir,"Fig3E_AdaptiveORFs_controlled_barplots.pdf",sep=""),width=7,height=2.5,  useDingbats=FALSE)


# statistical comparison: are there more selected effects in adaptive proto-genes?
pg_df<-df[df$emergence_status =="Emerging ORFs",]
# fitness upon loss
fisher.test(as.matrix(table(pg_df[,c("adaptive", "sick")])))
#p-value = 0.1478
# orf intactness in isolates
fisher.test(as.matrix(table(pg_df[,c("adaptive", "intact")])))
#p-value = 0.5422
# nucleotide diversity
wilcox.test(pg_df[pg_df$adaptive==1,"strain_pi"], pg_df[pg_df$adaptive==0,"strain_pi"], na.rm=TRUE)
#p-value = 0.6249

#######################################################
# Figure S4: Adaptive proto-genes across environments #
#######################################################

df<-merge(fitness_table, orf_table[,c("orf_name","emergence_status")], all.x=TRUE)
adaptive_protogenes_df<-df[df$emergence_status =="Emerging ORFs" & df$effect_cs=="increased",c("orf_name","exp_environment")]

# plot counts per environment (Fig S4A)
plot_x_label<-"Environmental condition"
plot_y_label<-"# adaptive emerging ORFs"
count_image<-ggplot(adaptive_protogenes_df, aes(exp_environment))+ geom_bar(color='deepskyblue4', fill='lightblue')+theme_bw()+labs(x=plot_x_label,y=plot_y_label) + theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.title.x=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.y=element_text(family="Helvetica", face="bold", colour="black", size=10), legend.position="none") 
ggsave(paste(general_outdir,"FigS4A_AdaptiveEmergingORFs_PerEnvironment_Counts.pdf",sep=""),width=3.5,height=2.5, useDingbats=FALSE)

# distribution of # environements (Fig S4B)
 plot_x_label<-"# environmental conditions"
plot_y_label<-"# adaptive emerging ORFs"
distrib<-table(adaptive_protogenes_df$orf_name)
plot_df<-as.data.frame(distrib)
distrib_image<-ggplot(plot_df[plot_df$Freq>0,], aes(x=Freq))+ geom_histogram(binwidth=1, color='deepskyblue4', fill='lightblue')+theme_bw()+labs(x=plot_x_label,y=plot_y_label) + theme(axis.text=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.title.x=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.y=element_text(family="Helvetica", face="bold", colour="black", size=10), legend.position="none") 
ggsave(paste(general_outdir,"FigS4B_AdaptiveEmerging ORFs_PerEnvironment_Distribution.pdf",sep=""),width=3.5,height=2.5, useDingbats=FALSE)

# simulation with a stochastic null model (Fig S4C) (null model: emerging orfs found adaptive are in fact just orfs that are not deleterious, and appear adaptive randomly as technical false positives)

# categories
adaptive_orfs<- unique(df[df $effect_cs=="increased",]$orf_name)
deleterious_orfs<-unique(df[!(df $orf_name %in% adaptive_orfs) & df $effect_cs=="decreased",]$orf_name)

#get deleterious proto-genes
delet_protos <- unique(df[df$emergence_status=="Emerging ORFs" & df$orf_name %in% deleterious_orfs,]$orf_name)
# non-deleterious proto-genes
nondelet_protos <- unique(df[df$emergence_status=="Emerging ORFs" & !(df$orf_name %in% delet_protos),]$orf_name)

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
# ISD, evolutionary status and fitness category (### modified after review to present as violin plots)
df<-orf_table[orf_table$barflex_space=="yes" & orf_table$emergence_status =="Emerging ORFs",c("orf_name", "Length","percent_disorder")]
fitness_df<-merge(fitness_table, orf_table[,c("orf_name","emergence_status")], all.x=TRUE)
adaptive_orfs<-unique(fitness_df[fitness_df $effect_cs=="increased",]$orf_name)
deleterious_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & fitness_df $effect_cs=="decreased",]$orf_name)
neutral_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & !(fitness_df $orf_name %in% deleterious_orfs),]$orf_name) 
df $fitness_category<-sapply(df $orf_name, function(x) ifelse(x %in% adaptive_orfs,"Adaptive",ifelse(x %in% deleterious_orfs, "Deleterious", ifelse(x %in% neutral_orfs,"Neutral","others"))))
df$fitness_category <-factor(df$fitness_category, levels=c("Adaptive", "Neutral", "Deleterious"))

# first check length distributions
wilcox.test(df[df$fitness_category == "Adaptive",]$Length, df[df$fitness_category == "Neutral",]$Length)
#p-value = 0.3434
wilcox.test(df[df$fitness_category == "Adaptive",]$Length, df[df$fitness_category == "Deleterious",]$Length)
#p-value = 0.6737
wilcox.test(df[df$fitness_category == "Deleterious",]$Length, df[df$fitness_category == "Neutral",]$Length)
#p-value = 0.7588

# Plot
df2<-df[,c("fitness_category","percent_disorder")]

image <- ggplot() +
  geom_jitter(data=df2,
              aes(x=fitness_category,
                  y=percent_disorder,
                  color=fitness_category),
              alpha=0.5, size=0.3) +
  geom_violin(data=df2,
              aes(x=fitness_category,
                  y=percent_disorder,
                  fill=fitness_category,
                  color=fitness_category),
              alpha=0.3) +
  stat_summary(data=df2,
               aes(x=fitness_category, y=percent_disorder, color=fitness_category),
               fun.y = "mean",
               geom="point",
               shape=8,
               size=2) +
  labs(x="",y="Disorder content") +
  scale_color_manual(values=c("darkgoldenrod1", "darkgray", "darkorchid3")) +
  scale_fill_manual(values=c("gold", "lightgray", "mediumpurple1")) +
  theme_bw() +
  theme(panel.spacing = unit(0, "lines"), axis.text.y=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.x=element_text(family="Helvetica", face="bold", colour="black", size=8), 	axis.title.y=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.x=element_blank(),	legend.position="none") 

ggsave(paste(general_outdir,"Fig4A_ISD_fitness_violin.pdf",sep=""),width=2.5,height=2.5,  useDingbats=FALSE)

# statistical comparison : adaptive versus neutral proto-genes wilxoc test
wilcox.test(df[df$fitness_category =="Adaptive",]$percent_disorder, df[df$fitness_category =="Neutral",]$percent_disorder)
#p-value = 0.03445
# statistical comparison : adaptive versus deleterious proto-genes wilxoc test
wilcox.test(df[df$fitness_category =="Deleterious",]$percent_disorder, df[df$fitness_category =="Adaptive",]$percent_disorder)
#p-value = 0.01584
# Cliff's delta adaptive vs. neutral protogenes
cliff.delta(df2[df2$fitness_category=="Adaptive", "percent_disorder"],
               df2[df2$fitness_category=="Neutral", "percent_disorder"])
  # delta estimate: -0.2484472 (small)
# Cliff's delta adaptive vs. deleterious protogenes
cliff.delta(df2[df2$fitness_category=="Adaptive", "percent_disorder"],
               df2[df2$fitness_category=="Deleterious", "percent_disorder"])
# delta estimate: -0.3116438 (small)



##########################################################################################
# Fig. S5A - ISD in established ORFs as a function of fitness #
###########################################################################################

# ISD, evolutionary status and fitness category (### modified after review to present as violin plots)
df<-orf_table[orf_table$barflex_space=="yes" & orf_table$emergence_status =="Established ORFs",c("orf_name", "percent_disorder")]
fitness_df<-merge(fitness_table, orf_table[,c("orf_name","protogene")], all.x=TRUE)
adaptive_orfs<-unique(fitness_df[fitness_df $effect_cs=="increased",]$orf_name)
deleterious_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & fitness_df $effect_cs=="decreased",]$orf_name)
neutral_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & !(fitness_df $orf_name %in% deleterious_orfs),]$orf_name) 
df $fitness_category<-sapply(df $orf_name, function(x) ifelse(x %in% adaptive_orfs,"Adaptive",ifelse(x %in% deleterious_orfs, "Deleterious", ifelse(x %in% neutral_orfs,"Neutral","others"))))
df$fitness_category <-factor(df$fitness_category, levels=c("Adaptive", "Neutral", "Deleterious"))

# Plot
df2<-df[,c("fitness_category","percent_disorder")]

image <- ggplot() +
  geom_jitter(data=df2,
              aes(x=fitness_category,
                  y=percent_disorder,
                  color=fitness_category),
              alpha=0.5, size=0.3) +
  geom_violin(data=df2,
              aes(x=fitness_category,
                  y=percent_disorder,
                  fill=fitness_category,
                  color=fitness_category),
              alpha=0.3) +
  stat_summary(data=df2,
               aes(x=fitness_category, y=percent_disorder, color=fitness_category),
               fun.y = "mean",
               geom="point",
               shape=8,
               size=2) +
  labs(x="",y="Disorder content") +
  scale_color_manual(values=c("darkgoldenrod1", "darkgray", "darkorchid3")) +
  scale_fill_manual(values=c("gold", "lightgray", "mediumpurple1")) +
  theme_bw() +
  theme(panel.spacing = unit(0, "lines"), axis.text.y=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.x=element_text(family="Helvetica", face="bold", colour="black", size=8), 	axis.title.y=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.x=element_blank(),	legend.position="none") 

ggsave(paste(general_outdir,"FigS5A_ISD_fitness_violin.pdf",sep=""),width=2.5,height=2.5,  useDingbats=FALSE)

# statistical comparison : adaptive versus neutral established orfs wilxoc test
wilcox.test(df[df$fitness_category =="Adaptive",]$percent_disorder, df[df$fitness_category =="Neutral",]$percent_disorder)
#p-value = 0.3364
# statistical comparison : adaptive versus deleterious established orfs wilxoc test
wilcox.test(df[df$fitness_category =="Deleterious",]$percent_disorder, df[df$fitness_category =="Adaptive",]$percent_disorder)
#p-value = 0.5484

# Cliff's delta adaptive vs. neutral established orfs
cliff.delta(df2[df2$fitness_category=="Adaptive", "percent_disorder"],
               df2[df2$fitness_category=="Neutral", "percent_disorder"])
               #delta estimate: 0.05057812 (negligible)
# Cliff's delta adaptive vs. deleterious established orfs
cliff.delta(df2[df2$fitness_category=="Adaptive", "percent_disorder"],
               df2[df2$fitness_category=="Deleterious", "percent_disorder"])
               #delta estimate: -0.03166468 (negligible)

##########################################################################################
# Fig. 4B - GC content in Emerging ORFs as a function of fitness  #
###########################################################################################

df<-orf_table[orf_table$barflex_space=="yes" & orf_table$emergence_status =="Emerging ORFs",c("orf_name", "GC")]
fitness_df<-merge(fitness_table, orf_table[,c("orf_name","emergence_status")], all.x=TRUE)
adaptive_orfs<-unique(fitness_df[fitness_df $effect_cs=="increased",]$orf_name)
deleterious_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & fitness_df $effect_cs=="decreased",]$orf_name)
neutral_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & !(fitness_df $orf_name %in% deleterious_orfs),]$orf_name) 
df $fitness_category<-sapply(df $orf_name, function(x) ifelse(x %in% adaptive_orfs,"Adaptive",ifelse(x %in% deleterious_orfs, "Deleterious", ifelse(x %in% neutral_orfs,"Neutral","others"))))
df$fitness_category <-factor(df$fitness_category, levels=c("Adaptive", "Neutral", "Deleterious"))

# plot
df2<-df[,c("fitness_category","GC")]
image <- ggplot() +
  geom_jitter(data=df2,
              aes(x=fitness_category,
                  y= GC,
                  color=fitness_category),
              alpha=0.5, size=0.3) +
  geom_violin(data=df2,
              aes(x=fitness_category,
                  y= GC,
                  fill=fitness_category,
                  color=fitness_category),
              alpha=0.3) +
  stat_summary(data=df2,
               aes(x=fitness_category, y= GC, color=fitness_category),
               fun.y = "mean",
               geom="point",
               shape=8,
               size=2) +
  labs(x="",y="GC content") +
  scale_color_manual(values=c("darkgoldenrod1", "darkgray", "darkorchid3")) +
  scale_fill_manual(values=c("gold", "lightgray", "mediumpurple1")) +
  theme_bw() +
  theme(panel.spacing = unit(0, "lines"), axis.text.y=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.x=element_text(family="Helvetica", face="bold", colour="black", size=8), 	axis.title.y=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.x=element_blank(),	legend.position="none") 

ggsave(paste(general_outdir,"Fig4B_GC_fitness_violin.pdf",sep=""),width=2.5,height=2.5,  useDingbats=FALSE)

# statistical comparison : adaptive versus neutral proto-genes wilxoc test
wilcox.test(df[df$fitness_category =="Adaptive",]$GC, df[df$fitness_category =="Neutral",]$GC)
#p-value = 0.04436
# statistical comparison : adaptive versus deleterious proto-genes wilxoc test
wilcox.test(df[df$fitness_category =="Deleterious",]$GC, df[df$fitness_category =="Adaptive",]$GC)
#W = 1152, p-value = 0.3258
# Cliff's delta adaptive vs. neutral protogenes
cliff.delta(df2[df2$fitness_category=="Adaptive", "GC"],
               df2[df2$fitness_category=="Neutral", "GC"])
              #delta estimate: -0.2362189 (small)
# Cliff's delta adaptive vs. deleterious protogenes
cliff.delta(df2[df2$fitness_category=="Adaptive", "GC"],
               df2[df2$fitness_category=="Deleterious", "GC"]) 
               #delta estimate: -0.1272016 (negligible)
  
##########################################################################################
# Fig. S5B - GC content in Established ORFs as a function of fitness  #
###########################################################################################

df<-orf_table[orf_table$barflex_space=="yes" & orf_table$emergence_status =="Established ORFs",c("orf_name", "GC")]
fitness_df<-merge(fitness_table, orf_table[,c("orf_name","protogene")], all.x=TRUE)
adaptive_orfs<-unique(fitness_df[fitness_df $effect_cs=="increased",]$orf_name)
deleterious_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & fitness_df $effect_cs=="decreased",]$orf_name)
neutral_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & !(fitness_df $orf_name %in% deleterious_orfs),]$orf_name) 
df $fitness_category<-sapply(df $orf_name, function(x) ifelse(x %in% adaptive_orfs,"Adaptive",ifelse(x %in% deleterious_orfs, "Deleterious", ifelse(x %in% neutral_orfs,"Neutral","others"))))
df$fitness_category <-factor(df$fitness_category, levels=c("Adaptive", "Neutral", "Deleterious"))

# plot
df2<-df[,c("fitness_category","GC")]
image <- ggplot() +
  geom_jitter(data=df2,
              aes(x=fitness_category,
                  y= GC,
                  color=fitness_category),
              alpha=0.5, size=0.3) +
  geom_violin(data=df2,
              aes(x=fitness_category,
                  y= GC,
                  fill=fitness_category,
                  color=fitness_category),
              alpha=0.3) +
  stat_summary(data=df2,
               aes(x=fitness_category, y= GC, color=fitness_category),
               fun.y = "mean",
               geom="point",
               shape=8,
               size=2) +
  labs(x="",y="GC content") +
  scale_color_manual(values=c("darkgoldenrod1", "darkgray", "darkorchid3")) +
  scale_fill_manual(values=c("gold", "lightgray", "mediumpurple1")) +
  theme_bw() +
  theme(panel.spacing = unit(0, "lines"), axis.text.y=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.x=element_text(family="Helvetica", face="bold", colour="black", size=8), 	axis.title.y=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.x=element_blank(),	legend.position="none") 

ggsave(paste(general_outdir,"FigS5B_GC_fitness_violin.pdf",sep=""),width=2.5,height=2.5,  useDingbats=FALSE)

# statistical comparison : adaptive versus neutral proto-genes wilxoc test
wilcox.test(df[df$fitness_category =="Adaptive",]$GC, df[df$fitness_category =="Neutral",]$GC)
#p-value = 0.1039
# statistical comparison : adaptive versus deleterious proto-genes wilxoc test
wilcox.test(df[df$fitness_category =="Deleterious",]$GC, df[df$fitness_category =="Adaptive",]$GC)
#p-value = 0.687
# Cliff's delta adaptive vs. neutral protogenes
cliff.delta(df2[df2$fitness_category=="Adaptive", "GC"],
               df2[df2$fitness_category=="Neutral", "GC"])
               #delta estimate: -0.08655961 (negligible)
# Cliff's delta adaptive vs. deleterious protogenes
cliff.delta(df2[df2$fitness_category=="Adaptive", "GC"],
               df2[df2$fitness_category=="Deleterious", "GC"]) 
               #delta estimate: 0.02125839 (negligible)               
               
##########################################################################################
# Fig. 4C. Transmembrane (TM) domains in Emerging ORFs as a function of fitness - tmhmm #
###########################################################################################
# TM, evolutionary status and fitness category
df<-orf_table[orf_table$barflex_space=="yes" & orf_table$emergence_status =="Emerging ORFs",c("orf_name","no_helices_TMHMM")]
df$is_tm<-ifelse(df$no_helices_TMHMM>0,1,0)
fitness_df<-merge(fitness_table, orf_table[,c("orf_name","emergence_status")], all.x=TRUE)
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
image<-ggplot(domains,aes(x= fitness,y=fraction,ymin=(fraction - sder),ymax=(fraction + sder),fill=fitness, color=fitness))+geom_bar(stat="identity",position=position_dodge(), alpha=0.5)+geom_errorbar(position = position_dodge(0.9),width=0.25)+labs(x="",y="With TM domain (TMHMMM)")+ scale_color_manual(values=c("darkgoldenrod1", "darkgray", "darkorchid3")) + scale_fill_manual(values=c("gold", "lightgray", "mediumpurple1"))+theme_bw()+theme(panel.spacing = unit(0, "lines"), axis.text.y=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.x=element_text(family="Helvetica", face="bold", colour="black", size=8), 	axis.title.y=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.x=element_blank(),	legend.position="none") 

ggsave(paste(general_outdir,"Fig4C_TMdomains_TMHMM_fitness_barplots.pdf",sep=""),width=2.5,height=2.5,  useDingbats=FALSE)

# statistical comparison : adaptive versus neutral emerging ORFs Fisher test
fisher.test(matrix(c(domains[domains $fitness =="Adaptive",]$yes, domains[domains $fitness =="Neutral",]$yes, domains[domains $fitness =="Adaptive",]$total - domains[domains $fitness =="Adaptive",]$yes, domains[domains $fitness =="Neutral",]$total - domains[domains $fitness =="Neutral",]$yes),2,2))

#data:  
# p-value = 0.002941
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
 # 1.439425 9.084369
# sample estimates:
# odds ratio 
  # 3.519211 
  
# adaptive versus deleterious emerging ORFs Fisher test
fisher.test(matrix(c(domains[domains $fitness =="Adaptive",]$yes, domains[domains $fitness =="Deleterious",]$yes, domains[domains $fitness =="Adaptive",]$total - domains[domains $fitness =="Adaptive",]$yes, domains[domains $fitness =="Deleterious",]$total - domains[domains $fitness =="Deleterious",]$yes),2,2))

# #data:  
# p-value = 0.001429
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
  # 1.613842 12.581587
# sample estimates:
# odds ratio 
  # 4.384247 


##########################################################################################
# Fig. 4D. Transmembrane (TM) domains in Emerging ORFs as a function of fitness - Phobius #
###########################################################################################
# TM, evolutionary status and fitness category
df<-orf_table[orf_table$barflex_space=="yes" & orf_table$emergence_status =="Emerging ORFs",c("orf_name","no_helices_PHOBIUS")]
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
image<-ggplot(domains,aes(x= fitness,y=fraction,ymin=(fraction - sder),ymax=(fraction + sder),fill=fitness, color=fitness))+geom_bar(stat="identity",position=position_dodge(), alpha=0.5)+geom_errorbar(position = position_dodge(0.9),width=0.25)+labs(x="",y="With TM domain (Phobius)")+ scale_color_manual(values=c("darkgoldenrod1", "darkgray", "darkorchid3")) + scale_fill_manual(values=c("gold", "lightgray", "mediumpurple1"))+theme_bw()+theme(panel.spacing = unit(0, "lines"), axis.text.y=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.x=element_text(family="Helvetica", face="bold", colour="black", size=8), 	axis.title.y=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.x=element_blank(),	legend.position="none") 

ggsave(paste(general_outdir,"Fig4D_TMdomains_Phobius_fitness_barplots.pdf",sep=""),width=2.5,height=2.5,  useDingbats=FALSE)

# statistical comparison : adaptive versus neutral proto-genes Fisher test
fisher.test(matrix(c(domains[domains $fitness =="Adaptive",]$yes, domains[domains $fitness =="Neutral",]$yes, domains[domains $fitness =="Adaptive",]$total - domains[domains $fitness =="Adaptive",]$yes, domains[domains $fitness =="Neutral",]$total - domains[domains $fitness =="Neutral",]$yes),2,2))

	# Fisher's Exact Test for Count Data

# data:  
# p-value = 0.02408
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
 # 1.107128 7.242959
# sample estimates:
# odds ratio 
  # 2.731452 

# adaptive versus deleterious emerging ORFs Fisher test
fisher.test(matrix(c(domains[domains $fitness =="Adaptive",]$yes, domains[domains $fitness =="Deleterious",]$yes, domains[domains $fitness =="Adaptive",]$total - domains[domains $fitness =="Adaptive",]$yes, domains[domains $fitness =="Deleterious",]$total - domains[domains $fitness =="Deleterious",]$yes),2,2))

	# Fisher's Exact Test for Count Data

# data:  
# p-value = 0.006659
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
  # 1.38705 10.91592
# sample estimates:
# odds ratio 
  # 3.762787 
  
##########################################################################################
# Fig. S5C. Transmembrane (TM) domains in established orfs as a function of fitness - tmhmm #
###########################################################################################
# TM, evolutionary status and fitness category
df<-orf_table[orf_table$barflex_space=="yes" & orf_table$emergence_status =="Established ORFs",c("orf_name","no_helices_TMHMM")]
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
image<-ggplot(domains,aes(x= fitness,y=fraction,ymin=(fraction - sder),ymax=(fraction + sder),fill=fitness, color=fitness))+geom_bar(stat="identity",position=position_dodge(), alpha=0.5)+geom_errorbar(position = position_dodge(0.9),width=0.25)+labs(x="",y="With TM domain (TMHMMM)")+ scale_color_manual(values=c("darkgoldenrod1", "darkgray", "darkorchid3")) + scale_fill_manual(values=c("gold", "lightgray", "mediumpurple1"))+theme_bw()+theme(panel.spacing = unit(0, "lines"), axis.text.y=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.x=element_text(family="Helvetica", face="bold", colour="black", size=8), 	axis.title.y=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.x=element_blank(),	legend.position="none") 

ggsave(paste(general_outdir,"FigS5C_TMdomains_TMHMM_fitness_barplots.pdf",sep=""),width=2.5,height=2.5,  useDingbats=FALSE)

# statistical comparison : adaptive versus neutral establsihed ORFs Fisher test
fisher.test(matrix(c(domains[domains $fitness =="Adaptive",]$yes, domains[domains $fitness =="Neutral",]$yes, domains[domains $fitness =="Adaptive",]$total - domains[domains $fitness =="Adaptive",]$yes, domains[domains $fitness =="Neutral",]$total - domains[domains $fitness =="Neutral",]$yes),2,2))

# data:  
# p-value = 0.1711
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
 # 0.8563892 2.0653819
# sample estimates:
# odds ratio 
  # 1.344598 

 # adaptive versus deleterious establsihed ORFs Fisher test
fisher.test(matrix(c(domains[domains $fitness =="Adaptive",]$yes, domains[domains $fitness =="Deleterious",]$yes, domains[domains $fitness =="Adaptive",]$total - domains[domains $fitness =="Adaptive",]$yes, domains[domains $fitness =="Deleterious",]$total - domains[domains $fitness =="Deleterious",]$yes),2,2))

# data:  
# p-value = 0.06239
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
 # 0.9511855 2.2771582
# sample estimates:
# odds ratio 
  # 1.488166 

##########################################################################################
# Fig. S5D. Transmembrane (TM) domains in established ORFs as a function of fitness - Phobius #
###########################################################################################
# TM, evolutionary status and fitness category
df<-orf_table[orf_table$barflex_space=="yes" & orf_table$emergence_status =="Established ORFs",c("orf_name","no_helices_PHOBIUS")]
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
image<-ggplot(domains,aes(x= fitness,y=fraction,ymin=(fraction - sder),ymax=(fraction + sder),fill=fitness, color=fitness))+geom_bar(stat="identity",position=position_dodge(), alpha=0.5)+geom_errorbar(position = position_dodge(0.9),width=0.25)+labs(x="",y="With TM domain (Phobius)")+ scale_color_manual(values=c("darkgoldenrod1", "darkgray", "darkorchid3")) + scale_fill_manual(values=c("gold", "lightgray", "mediumpurple1"))+theme_bw()+theme(panel.spacing = unit(0, "lines"), axis.text.y=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.x=element_text(family="Helvetica", face="bold", colour="black", size=8), 	axis.title.y=element_text(family="Helvetica", face="bold", colour="black", size=10), axis.title.x=element_blank(),	legend.position="none") 

ggsave(paste(general_outdir,"FigS5D_TMdomains_Phobius_fitness_barplots.pdf",sep=""),width=2.5,height=2.5,  useDingbats=FALSE)

# statistical comparison : adaptive versus neutral proto-genes Fisher test
fisher.test(matrix(c(domains[domains $fitness =="Adaptive",]$yes, domains[domains $fitness =="Neutral",]$yes, domains[domains $fitness =="Adaptive",]$total - domains[domains $fitness =="Adaptive",]$yes, domains[domains $fitness =="Neutral",]$total - domains[domains $fitness =="Neutral",]$yes),2,2))

# data:  
# p-value = 0.5804
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
 # 0.7228133 1.7552702
# sample estimates:
# odds ratio 
  # 1.139554 

 # adaptive versus deleterious establsihed ORFs Fisher test
fisher.test(matrix(c(domains[domains $fitness =="Adaptive",]$yes, domains[domains $fitness =="Deleterious",]$yes, domains[domains $fitness =="Adaptive",]$total - domains[domains $fitness =="Adaptive",]$yes, domains[domains $fitness =="Deleterious",]$total - domains[domains $fitness =="Deleterious",]$yes),2,2))

# data:  
# p-value = 0.5115
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
 # 0.7290033 1.7554379
# sample estimates:
# odds ratio 
  # 1.144729


###################################################################################################
# Fig. S6 - Distribution of %TM residues in established and emerging ORFs as a function of fitness#
###################################################################################################

###data
df<-orf_table[orf_table$barflex_space=="yes",c("orf_name","emergence_status", "percent_phobius", "percent_tmhmm")]
fitness_df<-merge(fitness_table, orf_table[,c("orf_name","emergence_status")], all.x=TRUE)
adaptive_orfs<-unique(fitness_df[fitness_df $effect_cs=="increased",]$orf_name)
deleterious_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & fitness_df $effect_cs=="decreased",]$orf_name)
neutral_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & !(fitness_df $orf_name %in% deleterious_orfs),]$orf_name) 
df $class<-sapply(df $orf_name, function(x) ifelse(x %in% adaptive_orfs,"adaptive",ifelse(x %in% deleterious_orfs, "deleterious", ifelse(x %in% neutral_orfs,"neutral","others"))))

### stats compare TM content between adaptive and neutral ORFs
cliff.delta(df[df$class =="adaptive" & df$emergence_status =="Emerging ORFs",]$percent_phobius, df[df$class =="neutral" & df$emergence_status =="Emerging ORFs",]$percent_phobius)
#delta estimate: 0.3293866 (small)
cliff.delta(df[df$class =="adaptive" & df$emergence_status =="Emerging ORFs",]$percent_tmhmm, df[df$class =="neutral" & df$emergence_status =="Emerging ORFs",]$percent_tmhmm)
#delta estimate: 0.4340062 (medium)
wilcox.test(df[df$class =="adaptive" & df$emergence_status =="Emerging ORFs",]$percent_phobius, df[df$class =="neutral" & df$emergence_status =="Emerging ORFs",]$percent_phobius)
#p-value = 0.002333
wilcox.test(df[df$class =="adaptive" & df$emergence_status =="Emerging ORFs",]$percent_tmhmm, df[df$class =="neutral" & df$emergence_status =="Emerging ORFs",]$percent_tmhmm)
#p-value = 0.0002186

cliff.delta(df[df$class =="adaptive" & df$emergence_status =="Established ORFs",]$percent_phobius, df[df$class =="neutral" & df$emergence_status =="Established ORFs",]$percent_phobius)
#delta estimate: 0.01539161 (negligible)
cliff.delta(df[df$class =="adaptive" & df$emergence_status =="Established ORFs",]$percent_tmhmm, df[df$class =="neutral" & df$emergence_status =="Established ORFs",]$percent_tmhmm)
#delta estimate: 0.02470608 (negligible)
wilcox.test(df[df$class =="adaptive" & df$emergence_status =="Established ORFs",]$percent_phobius, df[df$class =="neutral" & df$emergence_status =="Established ORFs",]$percent_phobius)
#p-value = 0.6921
wilcox.test(df[df$class =="adaptive" & df$emergence_status =="Established ORFs",]$percent_tmhmm, df[df$class =="neutral" & df$emergence_status =="Established ORFs",]$percent_tmhmm)
#p-value = 0.639

### stats compare TM content between adaptive and deleterious ORFs
cliff.delta(df[df$class =="adaptive" & df$emergence_status =="Emerging ORFs",]$percent_phobius, df[df$class =="deleterious" & df$emergence_status =="Emerging ORFs",]$percent_phobius)
#delta estimate: 0.4094912 (medium)
cliff.delta(df[df$class =="adaptive" & df$emergence_status =="Emerging ORFs",]$percent_tmhmm, df[df$class =="deleterious" & df$emergence_status =="Emerging ORFs",]$percent_tmhmm)
#delta estimate: 0.5342466 (large)
wilcox.test(df[df$class =="adaptive" & df$emergence_status =="Emerging ORFs",]$percent_phobius, df[df$class =="deleterious" & df$emergence_status =="Emerging ORFs",]$percent_phobius)
#p-value = 0.0004982
wilcox.test(df[df$class =="adaptive" & df$emergence_status =="Emerging ORFs",]$percent_tmhmm, df[df$class =="deleterious" & df$emergence_status =="Emerging ORFs",]$percent_tmhmm)
#p-value = 3.455e-05

cliff.delta(df[df$class =="adaptive" & df$emergence_status =="Established ORFs",]$percent_phobius, df[df$class =="deleterious" & df$emergence_status =="Established ORFs",]$percent_phobius)
#delta estimate: 0.02807226 (negligible)
cliff.delta(df[df$class =="adaptive" & df$emergence_status =="Established ORFs",]$percent_tmhmm, df[df$class =="deleterious" & df$emergence_status =="Established ORFs",]$percent_tmhmm)
#delta estimate: 0.06733829 (negligible)
wilcox.test(df[df$class =="adaptive" & df$emergence_status =="Established ORFs",]$percent_phobius, df[df$class =="deleterious" & df$emergence_status =="Established ORFs",]$percent_phobius)
#p-value = 0.4653
wilcox.test(df[df$class =="adaptive" & df$emergence_status =="Established ORFs",]$percent_tmhmm, df[df$class =="deleterious" & df$emergence_status =="Established ORFs",]$percent_tmhmm)
#p-value = 0.1976

### Plots

# emerging using phobius
thisdf<-df[df$emergence_status =="Emerging ORFs" , c("class","percent_phobius")]
thismu<-data.frame(class=c("adaptive","neutral","deleterious"),Mean = c(mean(thisdf[thisdf$class=="adaptive","percent_phobius"],na.rm=TRUE), mean(thisdf[thisdf$class=="neutral","percent_phobius"],na.rm=TRUE),mean(thisdf[thisdf$class=="deleterious","percent_phobius"],na.rm=TRUE) ))
plot_xlabel <-""
plot_ylabel<-"Emerging ORFs"
image_protoTMphobius<-ggplot(thisdf, aes(x= percent_phobius,color= class,fill= class))+geom_density(na.rm=TRUE)+
   geom_vline(data=thismu, aes(xintercept=Mean, color= class), linetype="dashed")+scale_color_manual(values=c("darkgoldenrod1", "darkorchid3","darkgray"))+xlim(c(0,0.8))+scale_fill_manual(values=c("gold","NA", "NA"))+labs(x=plot_xlabel,y=plot_ylabel) + theme_bw()+theme(axis.text.x=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.y=element_blank(),axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")+  theme(panel.background = element_rect(fill = "lightcyan2"), panel.grid.major = element_line(colour = "white"),   panel.grid.minor = element_line(colour = "white"))

# emerging using tmhmm
thisdf<-df[df$emergence_status =="Emerging ORFs" , c("class","percent_tmhmm")]
thismu<-data.frame(class=c("adaptive","neutral","deleterious"),Mean = c(mean(thisdf[thisdf$class=="adaptive","percent_tmhmm"],na.rm=TRUE), mean(thisdf[thisdf$class=="neutral","percent_tmhmm"],na.rm=TRUE),mean(thisdf[thisdf$class=="deleterious","percent_tmhmm"],na.rm=TRUE) ))
plot_xlabel <-""
plot_ylabel<-""
image_protoTMtmhmm<-ggplot(thisdf, aes(x= percent_tmhmm,color= class,fill= class))+geom_density(na.rm=TRUE)+
   geom_vline(data=thismu, aes(xintercept=Mean, color= class), linetype="dashed")+scale_color_manual(values=c("darkgoldenrod1", "darkorchid3","darkgray"))+xlim(c(0,0.8))+scale_fill_manual(values=c("gold","NA", "NA"))+labs(x=plot_xlabel,y=plot_ylabel) + theme_bw()+theme(axis.text.x=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.y=element_blank(),axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")+  theme(panel.background = element_rect(fill = "lightcyan2"), panel.grid.major = element_line(colour = "white"),   panel.grid.minor = element_line(colour = "white"))

# genes using phobius
thisdf<-df[df$emergence_status =="Established ORFs" , c("class","percent_phobius")]
thismu<-data.frame(class=c("adaptive","neutral","deleterious"),Mean = c(mean(thisdf[thisdf$class=="adaptive","percent_phobius"],na.rm=TRUE), mean(thisdf[thisdf$class=="neutral","percent_phobius"],na.rm=TRUE),mean(thisdf[thisdf$class=="deleterious","percent_phobius"],na.rm=TRUE) ))
plot_xlabel <-"TM residues content (Phobius)"
plot_ylabel<-"Established ORFs"
image_genesTMphobius<-ggplot(thisdf, aes(x= percent_phobius,color= class,fill= class))+geom_density(na.rm=TRUE)+
   geom_vline(data=thismu, aes(xintercept=Mean, color= class), linetype="dashed")+scale_color_manual(values=c("darkgoldenrod1", "darkorchid3","darkgray"))+xlim(c(0,0.8))+scale_fill_manual(values=c("gold","NA", "NA"))+labs(x=plot_xlabel,y=plot_ylabel) + theme_bw()+theme(axis.text.x=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.y=element_blank(),axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")

# genes using tmhmm
thisdf<-df[df$emergence_status =="Established ORFs" , c("class","percent_tmhmm")]
thismu<-data.frame(class=c("adaptive","neutral","deleterious"),Mean = c(mean(thisdf[thisdf$class=="adaptive","percent_tmhmm"],na.rm=TRUE), mean(thisdf[thisdf$class=="neutral","percent_tmhmm"],na.rm=TRUE),mean(thisdf[thisdf$class=="deleterious","percent_tmhmm"],na.rm=TRUE) ))
plot_xlabel <-"TM residues content (TMHMM)"
plot_ylabel<-""
image_genesTMtmhmm<-ggplot(thisdf, aes(x= percent_tmhmm,color= class,fill= class))+geom_density(na.rm=TRUE)+
   geom_vline(data=thismu, aes(xintercept=Mean, color= class), linetype="dashed")+scale_color_manual(values=c("darkgoldenrod1", "darkorchid3","darkgray"))+xlim(c(0,0.8))+scale_fill_manual(values=c("gold","NA", "NA"))+labs(x=plot_xlabel,y=plot_ylabel) + theme_bw()+theme(axis.text.x=element_text(family="Helvetica", face="plain", colour="black", size=8), axis.text.y=element_blank(),axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")

pdf(paste(general_outdir,"FigS6_TMResiduesContent_distributions.pdf",sep=""),width=4.5,height=5, useDingbats=FALSE)
multiplot(image_protoTMphobius, image_protoTMtmhmm, image_genesTMphobius, image_genesTMtmhmm, cols=2)
dev.off()



#####################################################
# - Figure 4E # Multinomial logistic regression TMHMM  #
#######################################################

# TM, evolutionary status and fitness category

df<-orf_table[orf_table$barflex_space=="yes",c("orf_name", "emergence_status", "percent_tmhmm")]
fitness_df<-merge(fitness_table, orf_table[,c("orf_name","emergence_status")], all.x=TRUE)
adaptive_orfs<-unique(fitness_df[fitness_df $effect_cs=="increased",]$orf_name)
deleterious_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & fitness_df $effect_cs=="decreased",]$orf_name)
neutral_orfs<-unique(fitness_df[!(fitness_df $orf_name %in% adaptive_orfs) & !(fitness_df $orf_name %in% deleterious_orfs),]$orf_name) 
df $fitness_category<-sapply(df $orf_name, function(x) ifelse(x %in% adaptive_orfs,"Adaptive",ifelse(x %in% deleterious_orfs, "Deleterious", ifelse(x %in% neutral_orfs,"Neutral","others"))))
df$fitness_category <-factor(df$fitness_category, levels=c("Adaptive", "Neutral", "Deleterious"))

# prepare data

df$fitness_category <- relevel(df$fitness_category, ref="Neutral")
df$emergence_status <- relevel(df$emergence_status, ref="Established ORFs")
df$orf_name <- as.character(df$orf_name)


# multinomial logistic regression
# fit
complex <- multinom(data=df, fitness_category~ emergence_status*percent_tmhmm)
simple <- multinom(data=df, fitness_category~ emergence_status +percent_tmhmm)
# coefficients p-values
z_complex <- summary(complex)$coefficients/summary(complex)$standard.errors
p_complex <- (1 - pnorm(abs(z_complex),0,1))*2
z_simple <- summary(simple)$coefficients/summary(simple)$standard.errors
p_simple <- (1 - pnorm(abs(z_simple),0,1))*2
# anova
anova(simple, complex)


# Likelihood ratio tests of Multinomial Models

# Response: fitness_category
                             # Model Resid. df Resid. Dev   Test    Df LR stat.     Pr(Chi)
# 1 emergence_status + percent_tmhmm      9288   7381.947                                  
# 2 emergence_status * percent_tmhmm      9286   7371.279 1 vs 2     2 10.66805 0.004824622


# cross validation to compare the predictive power of the two models
set.seed(35689118)
df_cv <- df
df_cv<-df_cv[sample(nrow(df_cv)),]
df_cv$folds <- cut(seq(1,nrow(df_cv)), breaks=10, labels=FALSE)

cv_df <- data.frame(complex=numeric(), simple=numeric())

for (i in 1:10)
{

  df_train <- df_cv %>%
    filter(folds!=i) %>%
    group_by(fitness_category) %>%
    sample_n(300, replace=TRUE)
             
  df_test <- df_cv %>%
    filter(folds==i)
 
  complex_train <- multinom(data=df_train, fitness_category~ emergence_status*percent_tmhmm)
  simple_train <- multinom(data=df_train, fitness_category~ emergence_status +percent_tmhmm)
 
  complex_test <- predict(complex_train,
                          newdata = df_test,
                          type = "class",
                          se = TRUE)
  simple_test <- predict(simple_train,
                          newdata = df_test,
                          type = "class",
                          se = TRUE)
 
  df_test$pred_comp <- complex_test
  df_test$pred_simp <- simple_test

	# complex Matthews correlation coefficient 
 
  TP <- sum(c(length(which(df_test[df_test$fitness_category=="Neutral","pred_comp"]=="Neutral")),
               length(which(df_test[df_test$fitness_category=="Adaptive","pred_comp"]=="Adaptive")),
               length(which(df_test[df_test$fitness_category=="Deleterious","pred_comp"]=="Deleterious"))))

  FP <- sum(c(length(which(df_test[df_test$fitness_category!="Neutral","pred_comp"]=="Neutral")),
               length(which(df_test[df_test$fitness_category!="Adaptive","pred_comp"]=="Adaptive")),
               length(which(df_test[df_test$fitness_category!="Deleterious","pred_comp"]=="Deleterious"))))
 
  TN <- sum(c(length(which(df_test[df_test$fitness_category!="Neutral","pred_comp"]!="Neutral")),
               length(which(df_test[df_test$fitness_category!="Adaptive","pred_comp"]!="Adaptive")),
               length(which(df_test[df_test$fitness_category!="Deleterious","pred_comp"]!="Deleterious"))))
 
  FN <- sum(c(length(which(df_test[df_test$fitness_category=="Neutral","pred_comp"]!="Neutral")),
               length(which(df_test[df_test$fitness_category=="Adaptive","pred_comp"]!="Adaptive")),
               length(which(df_test[df_test$fitness_category=="Deleterious","pred_comp"]!="Deleterious"))))
 
 	sum1 <- TP+FP; sum2 <-TP+FN ; sum3 <-TN+FP ; sum4 <- TN+FN;
 	denom <- as.double(sum1)*sum2*sum3*sum4
 	if (any(sum1==0, sum2==0, sum3==0, sum4==0)) {denom <- 1}
	mcc_comp <- ((TP*TN)-(FP*FN)) / sqrt(denom)

   # Simple Matthews correlation coefficient 
  TP <- sum(c(length(which(df_test[df_test$fitness_category=="Neutral","pred_simp"]=="Neutral")),
                length(which(df_test[df_test$fitness_category=="Adaptive","pred_simp"]=="Adaptive")),
                length(which(df_test[df_test$fitness_category=="Deleterious","pred_simp"]=="Deleterious"))))
 
  FP <- sum(c(length(which(df_test[df_test$fitness_category!="Neutral","pred_simp"]=="Neutral")),
              length(which(df_test[df_test$fitness_category!="Adaptive","pred_simp"]=="Adaptive")),
              length(which(df_test[df_test$fitness_category!="Deleterious","pred_simp"]=="Deleterious"))))
 
  TN <- sum(c(length(which(df_test[df_test$fitness_category!="Neutral","pred_simp"]!="Neutral")),
              length(which(df_test[df_test$fitness_category!="Adaptive","pred_simp"]!="Adaptive")),
              length(which(df_test[df_test$fitness_category!="Deleterious","pred_simp"]!="Deleterious"))))
 
  FN <- sum(c(length(which(df_test[df_test$fitness_category=="Neutral","pred_simp"]!="Neutral")),
              length(which(df_test[df_test$fitness_category=="Adaptive","pred_simp"]!="Adaptive")),
              length(which(df_test[df_test$fitness_category=="Deleterious","pred_simp"]!="Deleterious"))))

sum1 <- TP+FP; sum2 <-TP+FN ; sum3 <-TN+FP ; sum4 <- TN+FN;
 	denom <- as.double(sum1)*sum2*sum3*sum4
 	if (any(sum1==0, sum2==0, sum3==0, sum4==0)) {denom <- 1}
	mcc_simp <- ((TP*TN)-(FP*FN)) / sqrt(denom)

 	 cv_df <- rbind(cv_df,
                 data.frame(complex=mcc_comp,
                            simple=mcc_simp))
 
}
table(cv_df$complex > cv_df$simple)
#TRUE 
 # 10 
 
 apply(cv_df, 2, mean)
 # complex    simple 
#0.2853441 0.2285268 

wilcox.test(cv_df$complex, cv_df$simple, paired=TRUE)

#data:  cv_df$complex and cv_df$simple
#V = 55, p-value = 0.005889



# Plot

dpercent_tmhmm <- data.frame(emergence_status = rep(c("Established ORFs", "Emerging ORFs"), each = 42),
                             percent_tmhmm = rep(seq(0,1,by=1/41)[1:42], 2))

pp_complex <- cbind(dpercent_tmhmm,
                          predict(complex,
                                  newdata = dpercent_tmhmm,
                                  type = "probs",
                                  se = TRUE))

pp_simple <- cbind(dpercent_tmhmm,
                          predict(simple,
                                  newdata = dpercent_tmhmm,
                                  type = "probs",
                                  se=TRUE))

pp_complex.m <- melt(pp_complex, id.vars = c("emergence_status", "percent_tmhmm"))
pp_complex.m$model <- "Model with interaction"
pp_simple.m <- melt(pp_simple, id.vars = c("emergence_status", "percent_tmhmm"))
pp_simple.m$model <- "Model without interaction"

both_models_df <- bind_rows(pp_complex.m, pp_simple.m)
both_models_df$variable <- factor(both_models_df$variable, levels=c("Adaptive", "Neutral", "Deleterious"), ordered=TRUE)
both_models_df$model <- factor(both_models_df$model, levels=c("Model with interaction", "Model without interaction"), ordered=TRUE)

df$bins <- cut(df$percent_tmhmm, breaks=seq(0,1,by=1/41), include.lowest=TRUE, labels=seq(0,1,by=1/41)[1:41])

real_data <- df %>%
  group_by(emergence_status, bins) %>%
  dplyr::summarize(Neutral = table(fitness_category)[1]/dplyr::n(),
            Adaptive = table(fitness_category)[2]/dplyr::n(),
            Deleterious = table(fitness_category)[3]/dplyr::n())


real_data.m <- melt(real_data, id.vars = c("emergence_status", "bins"))

real_data.m$variable <- factor(real_data.m$variable, levels=c("Adaptive", "Neutral", "Deleterious"), ordered=TRUE)

image <- ggplot() +
  geom_bar(data=real_data.m,
           aes(x = as.numeric(as.character(bins)),
               y = value,
               fill = variable),
           stat="identity",
           alpha=0.5) +
  geom_line(data=both_models_df,
            aes(x=percent_tmhmm,
                y=value,
                colour=variable,
                linetype=model),
            size=1.5) +
  facet_grid(emergence_status ~variable, scales = "free") +
  scale_fill_manual(values=c("gold", "lightgray", "mediumpurple1")) +
  scale_colour_manual(values=c("darkgoldenrod1", "darkgray", "darkorchid3")) +
  ylab("Probability (lines), class frequency (bars)") +
  xlab("TM content") +
  theme_bw() +
  theme(panel.spacing = unit(0, "lines"),
        axis.text.y=element_text(family="Helvetica",
                                 face="plain",
                                 colour="black",
                                 size=7),
        axis.text.x=element_text(family="Helvetica",
                                 face="plain",
                                 colour="black",
                                 size=7),
        axis.title.y=element_text(family="Helvetica",
                                  face="bold",
                                  colour="black",
                                  size=10),
		axis.title.x=element_text(family="Helvetica",
                                  face="bold",
                                  colour="black",
                                  size=10),
        legend.position="bottom",
        legend.title = element_blank()) +
  guides(fill=FALSE, color=FALSE)

ggsave(paste(general_outdir,"Fig4E_MultinomRegressPredProb_2models_TMHMM.pdf",sep=""),width=5,height=4)


# pvalues for the terms of the complex model:
p_complex
#            (Intercept) emergence_statusEmerging ORFs percent_tmhmm emergence_statusEmerging ORFs:percent_tmhmm
#Adaptive              0                  9.064429e-01  9.580590e-01                                 0.001643957
#Deleterious           0                  2.579341e-10  3.469761e-05                                 0.971700577


####################################################################################################
# - Figure 5A - Fraction of real and scrambled sequences with putative TM domain # NO length control #
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
statistics_table $evoclass <-factor(statistics_table $evoclass, levels=c("established ORFs","emerging ORFs","iORFs","sORFs"))
statistics_table
  # total  yes  fraction        sder         evoclass  seqclass
# 1  5645 1296 0.2295837 0.005597593 established ORFs      real
# 2   662  304 0.4592145 0.019368292    emerging ORFs      real
# 3  6177 3499 0.5664562 0.006305374            iORFs      real
# 4 18503 3985 0.2153705 0.003022069            sORFs      real
# 5  5645 2004 0.3550044 0.006368883 established ORFs scrambled
# 6   662  261 0.3942598 0.018993520    emerging ORFs scrambled
# 7  6177 2387 0.3864335 0.006195544            iORFs scrambled
# 8 18503 3256 0.1759715 0.002799440            sORFs scrambled

# plot
image<-ggplot(statistics_table,aes(x= seqclass,y=fraction,ymin=(fraction - sder),ymax=(fraction + sder),fill=seqclass))+geom_bar(stat="identity",position=position_dodge(), color="black")+geom_errorbar(position = position_dodge(0.9),width=0.25)+labs(x="",y="Fraction with TM domain")+ ylim(0,0.7)+ scale_fill_manual(values=c("pink", scramble_color))  + facet_grid(.~evoclass)+theme_bw()+theme(panel.spacing = unit(0.2, "lines"), axis.text=element_text(family="Helvetica", face="bold", colour="black", size=8),  	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")
ggsave(paste(general_outdir,"Fig5A_FractionWithTM_NoLengthcontrol.pdf",sep=""),width=5,height=2.5,  useDingbats=FALSE)

# statistical test for difference in numbers with at least one TM in established, emerging and artificial ORFs
statistics_table_temp <-  statistics_table[1:3,]
statistics_table_temp$no <- statistics_table_temp$total - statistics_table_temp$yes
chisq.test(as.matrix(statistics_table_temp[,c("yes","no")]))
#X-squared = 1392.8, df = 2, p-value < 2.2e-16

######################################################################################################
#  Figure S7A - Fraction of real and scrambled sequences with putative TM domain - WITH Length Control - #
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
set.seed(35689118)
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
set.seed(35689118)
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
set.seed(35689118)
t<-df[df $type=="iORFs" & df $mode=="real",]
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
statistics_table $evoclass <-factor(statistics_table $evoclass, levels=c("established ORFs","emerging ORFs","iORFs","sORFs"))
statistics_table
  # total yes  fraction       sder         evoclass  seqclass
# 1  1000 214 0.2140000 0.01296935 established ORFs      real
# 2   662 304 0.4592145 0.01936829    emerging ORFs      real
# 3  1000 556 0.5560000 0.01571191            iORFs      real
# 4  1000 303 0.3030000 0.01453241            sORFs      real
# 5  1000 127 0.1270000 0.01052953 established ORFs scrambled
# 6   662 261 0.3942598 0.01899352    emerging ORFs scrambled
# 7  1000 381 0.3810000 0.01535705            iORFs scrambled
# 8  1000 244 0.2440000 0.01358175            sORFs scrambled

# plot
image<-ggplot(statistics_table,aes(x= seqclass,y=fraction,ymin=(fraction - sder),ymax=(fraction + sder),fill=seqclass))+geom_bar(stat="identity",position=position_dodge(), color="black")+geom_errorbar(position = position_dodge(0.9),width=0.25)+labs(x="",y="Fraction with TM domain")+ ylim(0,0.7)+ scale_fill_manual(values=c("pink", scramble_color))  + facet_grid(.~evoclass)+theme_bw()+theme(panel.spacing = unit(0.2, "lines"), axis.text=element_text(family="Helvetica", face="bold", colour="black", size=8),  	axis.title=element_text(family="Helvetica", face="bold", colour="black", size=10),	legend.position="none")
ggsave(paste(general_outdir,"FigS7A_FractionWithTM_Lengthcontrol.pdf",sep=""),width=5,height=2.5,  useDingbats=FALSE)

# statistical test for difference in numbers with at least one TM in established, emerging and artificial ORFs
statistics_table_temp <-  statistics_table[1:3,]
statistics_table_temp$no <- statistics_table_temp$total - statistics_table_temp$yes
chisq.test(as.matrix(statistics_table_temp[,c("yes","no")]))

# fischer tests, presence of TM in real vs. scrambled established orfs
fisher.test(table(df[which(df$type=="established ORFs"),c("mode","has_tm")]))
	# Fisher's Exact Test for Count Data

# data:  table(df[which(df$type == "established ORFs"), c("mode", "has_tm")])
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
 # 1.699168 2.007830
# sample estimates:
# odds ratio 
  # 1.846904 
  1/fisher.test(table(df[which(df$type=="established ORFs"),c("mode","has_tm")]))$estimate
#odds ratio 
# 0.5414467 

# established ORFs with at least one TM, comparison of TM content real vs. scrambled
wilcox.test(df[which(df$type=="established ORFs" & df$has_tm==1 & df$mode=="real"),"no_res_TM_Phob"]/df[which(df$type=="established ORFs" & df$has_tm==1 & df$mode=="real"),"length"],
            df[which(df$type=="established ORFs" & df$has_tm==1 & df$mode=="scrambled"),"no_res_TM_Phob"]/df[which(df$type=="established ORFs" & df$has_tm==1 & df$mode=="scrambled"),"length"])
  #W = 1782600, p-value < 2.2e-16
cliff.delta(df[which(df$type=="established ORFs" & df$has_tm==1 & df$mode=="real"),"no_res_TM_Phob"]/df[which(df$type=="established ORFs" & df$has_tm==1 & df$mode=="real"),"length"],
            df[which(df$type=="established ORFs" & df$has_tm==1 & df$mode=="scrambled"),"no_res_TM_Phob"]/df[which(df$type=="established ORFs" & df$has_tm==1 & df$mode=="scrambled"),"length"])
  # delta estimate: 0.3727375 (medium)          
            
# fischer tests, presence of TM in real vs. scrambled
fisher.test(table(df[which(df$type=="iORFs"),c("mode","has_tm")]))
# data:  table(df[which(df$type == "iORFs"), c("mode", "has_tm")])
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
 # 0.4483301 0.5182639
# sample estimates:
# odds ratio 
 # 0.4820637 
 1/fisher.test(table(df[which(df$type=="iORFs"),c("mode","has_tm")]))$estimate
 #odds ratio 
 # 2.074415 
 
fisher.test(table(df[which(df$type=="emerging ORFs"),c("mode","has_tm")]))
	# Fisher's Exact Test for Count Data

# data:  table(df[which(df$type == "emerging ORFs"), c("mode", "has_tm")])
# p-value = 0.01957
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
 # 0.6124143 0.9592527
# sample estimates:
# odds ratio 
 # 0.7666518 
 1/fisher.test(table(df[which(df$type=="emerging ORFs"),c("mode","has_tm")]))$estimate
#odds ratio 
#  1.304373 

######################################################################################################################
# Figure 5B - TM propensity in real and scrambled sequences as a function of Thymine content - NO Length Control -  #
######################################################################################################################
# data
df<-tm_table
df$has_tm<-ifelse(df$no_hel_Phob>0,1,0)
df $type <-factor(df $type, levels=c("established ORFs","emerging ORFs","iORFs","sORFs"))

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

#################################################################################################################
# Figure S7B - TM propensity in real and scrambled sequences as a function of Thymine content - WITH Length Control -  #
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
set.seed(35689118)
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
set.seed(35689118)
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
set.seed(35689118)
t<-df[df $type=="iORFs" & df $mode=="real",]
adj_protobreaks<-c(seq(25,300,by=25 ),650)
adj_protoprob<-c(sum(protoprob[1:2]), protoprob[3:12], sum(protoprob[13:length(protoprob)]))
t$lengthbin<-cut(t$length, breaks=adj_protobreaks)
t<-t[!(is.na(t$lengthbin)),]
bindistrib<-as.data.frame(table(sample(x=levels(t$lengthbin), 1000, replace=T, prob= adj_protoprob)))
l<-sapply(bindistrib$Var1, function(x) sample(t[t$lengthbin==x,]$names, bindistrib[bindistrib$Var1==x,]$Freq, replace=T))
rownumbers<-NULL
for(x in unlist(l)){rownumbers<-c(rownumbers, which(df$names==x,))}
controlled_df <-rbind(controlled_df, df[rownumbers,])
controlled_df $type <-factor(controlled_df $type, levels=c("established ORFs","emerging ORFs","iORFs","sORFs"))

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
#Figure 5C - Figure 5C TM propensity in sORFs (25 - 75 codons) as a function of their length relative to their containing intergene -  #
#####################################################################################################################################
# data
df<-ngo_intergene_table
otherdf<-tm_table

# TM propensity to expect for emerging ORFs and iORFs of same length:
protovalue<-100*nrow(otherdf[otherdf $type=="emerging ORFs" & otherdf $mode=="real" & otherdf $no_hel_Phob>0 & otherdf $length>25 & otherdf $length<75,])/nrow(otherdf[otherdf $type=="emerging ORFs" & otherdf $mode=="real" & otherdf $length>25 & otherdf $length<75,])
intergenevalue<-100*nrow(otherdf[otherdf $type=="iORFs" & otherdf $mode=="real" & otherdf $no_hel_Phob >0 & otherdf $length>25 & otherdf $length<75,])/nrow(otherdf[otherdf $type=="iORFs" & otherdf $mode=="real" & otherdf $length>25 & otherdf $length<75,])

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

# correlation between % of intergene's length covered by ORF and % with TM domain
cor.test(t$bin, t$fraction, method="spearman")

# data:  t$bin and t$fraction
# S = 96.542, p-value = 3.98e-08
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
      # rho 
# 0.9153138
