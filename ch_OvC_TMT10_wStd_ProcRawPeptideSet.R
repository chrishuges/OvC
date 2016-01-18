# TODO: Add comment
# 
# Author: cshughes
###############################################################################
library(RColorBrewer)
library(lattice)
library(gplots)
library(affy)
library(preprocessCore)
library(limma)
library(cluster)
library(vsn)
library(Biobase)
library(LSD)
###############################################################################
#load in files
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/")
s1apsm<-read.table("./ch_Dec2015_OvC_TMT10_Set1A_PSMs.txt", header=TRUE, sep='\t')
s1apro<-read.table("./ch_Dec2015_OvC_TMT10_Set1A_Proteins.txt", header=TRUE, sep='\t')
s1bpsm<-read.table("./ch_Dec2015_OvC_TMT10_Set1B_PSMs.txt", header=TRUE, sep='\t')
s1bpro<-read.table("./ch_Dec2015_OvC_TMT10_Set1B_Proteins.txt", header=TRUE, sep='\t')
s1cpsm<-read.table("./ch_Dec2015_OvC_TMT10_Set1C_PSMs.txt", header=TRUE, sep='\t')
s1cpro<-read.table("./ch_Dec2015_OvC_TMT10_Set1C_Proteins.txt", header=TRUE, sep='\t')
s2apsm<-read.table("./ch_Dec2015_OvC_TMT10_Set2A_PSMs.txt", header=TRUE, sep='\t')
s2apro<-read.table("./ch_Dec2015_OvC_TMT10_Set2A_Proteins.txt", header=TRUE, sep='\t')
s2bpsm<-read.table("./ch_Dec2015_OvC_TMT10_Set2B_PSMs.txt", header=TRUE, sep='\t')
s2bpro<-read.table("./ch_Dec2015_OvC_TMT10_Set2B_Proteins.txt", header=TRUE, sep='\t')
s2cpsm<-read.table("./ch_Dec2015_OvC_TMT10_Set2C_PSMs.txt", header=TRUE, sep='\t')
s2cpro<-read.table("./ch_Dec2015_OvC_TMT10_Set2C_Proteins.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput/")

##################################################
#initial processing of the data to compile peptides and remove contaminants
##################################################
#dont split replicates yet or do any combination of them
processPSM <- function(psmFile, proteinFile, Rep, ... ){	
	require(limma)
	#initial processing to remove non-informative columns and rows from contaminant proteins or those with many missing values
	pep<-psmFile[,c(5,6,9:11,2,37:47)]
	pro<-proteinFile[,c(3:4)]
	print('Raw data frame dimensions')
	print(dim(pep))
	#remove non-unique peptides or those without quan values
	pep<-subset(pep, Number.of.Proteins==1)
	pep<-subset(pep, !grepl('NoQuanValues',pep$Quan.Info))
	print('Unique Peptides Only')
	print(dim(pep))
	#parse the protein accession
	pep$Accession = sapply(strsplit(as.character(pep$Master.Protein.Accessions), ';'),'[', 1)
	#merge the protein descriptions into the peptide file
	pep.m = merge(pep,pro,by='Accession')
	#get the gene name out
	pep.m$Gene<-sub(".*?GN=(.*?)( .*|$)", "\\1", pep.m$Description)
	#get the organism accession out
	pep.m$Org<-'Human'
	ecoli = grepl('Escherichia', pep.m$Description)
	pep.m[ecoli,21] = 'EColi'
	#parse the peptide column for amino acids
	pep.m$Sequence = toupper(sub('.*?\\.(.*?)(\\..*|$)','\\1',pep.m$Annotated.Sequence))
	#filter information from Description
	pep.m$Descriptions = sub('(.*?)( OS=.*|$)','\\1',pep.m$Description)
	#remove specific proteins
	pep.m = subset(pep.m, !grepl('sp',pep.m$Accession, ignore.case=FALSE))
	#filter the data columns
	print('removing contaminants')
	print(dim(pep.m))
	pep.r = pep.m[,c(1,20,21,23,22,7,9:18)]
	#aggregate the PSMs into peptides
	print('aggregating peptides')
	pep.a = aggregate(cbind(X126,X127N,X127C,X128N,X128C,X129N,X129C,X130N,X130C,X131)~Accession+Gene+Org+Descriptions+Sequence+Confidence,data=pep.r,median,na.action=na.pass,na.rm=TRUE)
	colnames(pep.a) = c('Accession','Gene','Organism','Descriptions','Sequence','Confidence','a1','a2','a3','a4','a5','a6','a7','a8','a9','aStd')
	print(dim(pep.a))
	#replace NaN with NA
	pep.a[,7:16] = round(pep.a[,7:16],2)
	pep.a[,7:16][is.na(pep.a[,7:16])]<-NA
	#filter based on S/N
	pep.f1 = subset(pep.a, rowMeans(pep.a[,7:16], na.rm=TRUE)>=5)
	#filter based on NA
	pep.f2 = subset(pep.f1, !is.na(pep.f1[,16]))
	pep.f3 = subset(pep.f2, rowSums(is.na(pep.f2[,7:15]))<9)
	#add replicate counter
	pep.f3$Rep = Rep
	pep.f3$Acc = paste(pep.f3$Gene,'_',pep.f3$Accession, sep='')
	#output the data
	return(pep.f3)	
}
#run the function
ova1.psm<-processPSM(s1apsm, s1apro, 'a1')
ova2.psm<-processPSM(s1bpsm, s1bpro, 'a2')
ova3.psm<-processPSM(s1cpsm, s1cpro, 'a3')
ovb1.psm<-processPSM(s2apsm, s2apro, 'b1')
ovb2.psm<-processPSM(s2bpsm, s2bpro, 'b2')
ovb3.psm<-processPSM(s2cpsm, s2cpro, 'b3')
#output the data objects
saveRDS(ova1.psm,'ch_OvC_wStd_processedPeptides_a1.rds')
saveRDS(ova2.psm,'ch_OvC_wStd_processedPeptides_a2.rds')
saveRDS(ova3.psm,'ch_OvC_wStd_processedPeptides_a3.rds')
saveRDS(ovb1.psm,'ch_OvC_wStd_processedPeptides_b1.rds')
saveRDS(ovb2.psm,'ch_OvC_wStd_processedPeptides_b2.rds')
saveRDS(ovb3.psm,'ch_OvC_wStd_processedPeptides_b3.rds')


