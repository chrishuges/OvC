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
setwd(dir="/Users/cshughes/Documents/projects/OvC/FFPE/")
fa1psm<-read.table("./ch_08July2015_OvC-TMT10_HpH_repA1_med_PSMs.txt", header=TRUE, sep='\t')
fa1pro<-read.table("./ch_08July2015_OvC-TMT10_HpH_repA1_med_Proteins.txt", header=TRUE, sep='\t')
fa2psm<-read.table("./ch_08July2015_OvC-TMT10_HpH_repA2_med_PSMs.txt", header=TRUE, sep='\t')
fa2pro<-read.table("./ch_08July2015_OvC-TMT10_HpH_repA2_med_Proteins.txt", header=TRUE, sep='\t')
fb2psm<-read.table("./ch_09July2015_OvC-TMT10_HpH_repB2_med_PSMs.txt", header=TRUE, sep='\t')
fb2pro<-read.table("./ch_09July2015_OvC-TMT10_HpH_repB2_med_Proteins.txt", header=TRUE, sep='\t')
fb1psm<-read.table("./ch_09July2015_OvC-TMT10_HpH_repB1_med_PSMs.txt", header=TRUE, sep='\t')
fb1pro<-read.table("./ch_09July2015_OvC-TMT10_HpH_repB1_med_Proteins.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/FFPE/Routput/")

##################################################
#initial processing of the data to compile peptides and remove contaminants
##################################################
#dont split replicates yet or do any combination of them
processPSM <- function(psmFile, proteinFile, replicate, ... ){	
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
	colnames(pep.a) = c('Accession','Gene','Organism','Descriptions','Sequence','Confidence','a1','a2','a3','a4','a5','a6','a7','a8','a9','a10')
	print(dim(pep.a))
	#replace NaN with NA
	pep.a[,7:16] = round(pep.a[,7:16],2)
	pep.a[,7:16][is.na(pep.a[,7:16])]<-NA
	#filter based on S/N
	pep.f1 = subset(pep.a, rowMeans(pep.a[,7:16], na.rm=TRUE)>5)
	#filter based on NA
	pep.f2 = subset(pep.f1, rowSums(is.na(pep.f1[,7:16]))<8)
	#add replicate counter
	pep.f2$Rep = replicate
	#output the data
	return(pep.f2)	
}
#run the function
ovfa1.psm<-processPSM(fa1psm, fa1pro, 'a1')
ovfa2.psm<-processPSM(fa2psm, fa2pro, 'a2')
ovfb1.psm<-processPSM(fb1psm, fb1pro, 'a3')
ovfb2.psm<-processPSM(fb2psm, fb2pro, 'a4')

#output the data objects
saveRDS(ovfa1.psm,'ch_OvC_FFPE_processedPeptides_fa1.rds')
saveRDS(ovfa2.psm,'ch_OvC_FFPE_processedPeptides_fa2.rds')
saveRDS(ovfb1.psm,'ch_OvC_FFPE_processedPeptides_fb1.rds')
saveRDS(ovfb2.psm,'ch_OvC_FFPE_processedPeptides_fb2.rds')




