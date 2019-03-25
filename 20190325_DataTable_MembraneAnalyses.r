# March 1 2019; just reran March 25 2019
# Anne-Ruxandra Carvunis

# this script builds the tables necessary to reproduce membrane analyses in Figure 5 from tables provided by Nikos 
 

general_outdir<-"/Users/annerux/Workspace/cerevisiae/txt/screen_analyses/Protogenes/ProtogeneFitnessFigurePanels/submission/" 

data_tm<-read.csv("/Users/annerux/Workspace/cerevisiae/txt/external_datasets/NikosProteinStructure/20181009_Figure5AB_data.csv")

data_intergenetype<-read.csv("/Users/annerux/Workspace/cerevisiae/txt/external_datasets/NikosProteinStructure/20181009_Figure5C_data.csv")


############
# Data table for Fig 5 A B and sup fig 6 and 7

datatable<-data_tm[data_tm$length>=25,]
levels(datatable$type)<-c(levels(datatable$type),"annotated ORFs", 
"sORFs", "artificial ORFs")
datatable$type[datatable$type=="protogenes"]<-"annotated ORFs"
datatable$type[datatable$type=="genes"]<-"annotated ORFs"
datatable$type[datatable$type=="non-genic ORFs"]<-"sORFs"
datatable$type[datatable$type=="intergenes"]<-"artificial ORFs"
write.csv(datatable, paste(general_outdir, "Data4.csv", sep=""))

###############
# Data table for Fig 5C, examination of TM in the contect of non-genic ORfs as a function of intergene fraction
# only consider ORFs between 25 and 75 aa
datatable<-data_intergenetype[data_intergenetype$NTO_length_aa>25 & data_intergenetype$NTO_length_aa<75,]
# bin
midcut<-function(x,from,to,by){
   ## cut the data into bins...
   x=cut(x,seq(from,to,by),include.lowest=T)
   ## make a named vector of the midpoints, names=binnames
   vec=seq(from+by/2,to-by/2,by)
   names(vec)=levels(x)
   ## use the vector to map the names of the bins to the midpoint values
   unname(vec[x])
}
datatable $fracbin_anne<-midcut(datatable $fract_lengths, 0,1 ,0.05)
datatable<- datatable[,c(1,2,3,4,6,10,c(13:19))]
names(datatable)<-c("nameNGO","nameINTER", "INTER_length_nt", "NGO_length_aa", "NGO_no_hel", "fract_lengths",  "coorStartNGO", "coorEndNGO", "coorStartINTER", "coorEndINTER", "ChromNGO", "ChromINTER","bin")

# EXPORT
# text
write.csv(datatable, paste(general_outdir, "Data5.csv", sep=""))



