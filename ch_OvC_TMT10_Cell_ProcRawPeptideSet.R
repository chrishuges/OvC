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
setwd(dir="/Users/cshughes/Documents/projects/OvC/CellLine/")
fa1psm<-read.table("./ch_04Aug2015_OvCcell-TMT6_HpH_repA1_PSMs.txt", header=TRUE, sep='\t')
fa1pro<-read.table("./ch_04Aug2015_OvCcell-TMT6_HpH_repA1_Proteins.txt", header=TRUE, sep='\t')
fa2psm<-read.table("./ch_04Aug2015_OvCcell-TMT6_HpH_repA2_PSMs.txt", header=TRUE, sep='\t')
fa2pro<-read.table("./ch_04Aug2015_OvCcell-TMT6_HpH_repA2_Proteins.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/CellLine/Routput/")

##################################################
#initial processing of the data to compile peptides and remove contaminants
##################################################
#dont split replicates yet or do any combination of them
processPSM <- function(psmFile, proteinFile, replicate, ... ){	
	require(limma)
	#initial processing to remove non-informative columns and rows from contaminant proteins or those with many missing values
	pep<-psmFile[,c(5,6,9:11,2,37:43)]
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
	pep.m[ecoli,17] = 'EColi'
	#parse the peptide column for amino acids
	pep.m$Sequence = toupper(sub('.*?\\.(.*?)(\\..*|$)','\\1',pep.m$Annotated.Sequence))
	#filter information from Description
	pep.m$Descriptions = sub('(.*?)( OS=.*|$)','\\1',pep.m$Description)
	#remove specific proteins
	pep.m = subset(pep.m, !grepl('sp',pep.m$Accession, ignore.case=FALSE))
	#filter the data columns
	print('removing contaminants')
	print(dim(pep.m))
	pep.r = pep.m[,c(1,16,17,19,18,7,9:14)]
	#aggregate the PSMs into peptides
	print('aggregating peptides')
	pep.a = aggregate(cbind(X126,X127,X128,X129,X130,X131)~Accession+Gene+Org+Descriptions+Sequence+Confidence,data=pep.r,median,na.action=na.pass,na.rm=TRUE)
	colnames(pep.a) = c('Accession','Gene','Organism','Descriptions','Sequence','Confidence','a1','a2','a3','a4','a5','a6')
	print(dim(pep.a))
	#replace NaN with NA
	pep.a[,7:12] = round(pep.a[,7:12],2)
	pep.a[,7:12][is.na(pep.a[,7:12])]<-NA
	#filter based on S/N
	pep.f1 = subset(pep.a, rowMeans(pep.a[,7:12], na.rm=TRUE)>5)
	#filter based on NA
	pep.f2 = subset(pep.f1, rowSums(is.na(pep.f1[,7:12]))<8)
	#add replicate counter
	pep.f2$Rep = replicate
	#output the data
	return(pep.f2)	
}
#run the function
ovfa1.psm<-processPSM(fa1psm, fa1pro, 'a1')
ovfa2.psm<-processPSM(fa2psm, fa2pro, 'a2')

#output the data objects
saveRDS(ovfa1.psm,'ch_OvC_CellLine_processedPeptides_fa1.rds')
saveRDS(ovfa2.psm,'ch_OvC_CellLine_processedPeptides_fa2.rds')




