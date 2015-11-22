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

setwd(dir="/Users/cshughes/Documents/projects/OvC/")
spsm<-read.table("./ch_02July2015_OvC-TMT10_HpH_tRep1_PSMs.txt", header=TRUE, sep='\t')
spro<-read.table("./ch_02July2015_OvC-TMT10_HpH_tRep1_Proteins.txt", header=TRUE, sep='\t')
lpsm<-read.table("./ch_02July2015_OvC-TMT10_HpH_tRep1-(01)_PSMs.txt", header=TRUE, sep='\t')
lpro<-read.table("./ch_02July2015_OvC-TMT10_HpH_tRep1-(01)_Proteins.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/Routput/")

#this function takes the proteome discoverer data PSM set and filters it to a usable set
processPSM <- function(psmFile, proteinFile, ... ){	
	require(limma)
	#initial processing to remove non-informative columns and rows from contaminant proteins or those with many missing values
	pep<-psmFile[,c(5,6,10,11,25:27,36:46)]
	pro<-proteinFile[,c(3,4)]
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
	#parse the peptide column for amino acids
	pep.m$Sequence = sub('.*?\\.(.*?)(\\..*|$)','\\1',pep.m$Annotated.Sequence)
	#filter information from Description
	pep.m$Descriptions = sub('(.*?)( OS=.*|$)','\\1',pep.m$Description)
	#output the data
	return(pep.m[,c(1,21,23,22,3,9:19,6:8)])
}

s.pro<-processPSM(spsm,spro)
l.pro<-processPSM(lpsm,lpro)



##########################################
#do the quality control
##########################################
#change first line dependent on data set you are working on
sk = s.pro
sk[sk==0]<-NA
#x[,7:16][x[,7:16]<6]<-NA

#filter for values that are present in all channels
sk = subset(sk, rowSums(is.na(sk[,c(7:16)]))<1)

#MA plot of the raw values
#pre-normalization
M= log2(sk[,7]) - log2(sk[,8])
A= (log2(sk[,7]) + log2(sk[,8]))/2
#plot
pdf('ch_preNormExpr_MAplot_240inj.pdf')
heatscatter(A,
		M,
		cex=1,
		pch=19,
		xlab = 'Average Intensity',
		ylab = 'Intensity Ratio')
lines(lowess(A,M), col=brewer.pal(6,"RdBu")[1],lwd=3,lty=2)
abline(h=0,col=brewer.pal(6,"RdBu")[6],lwd=3)
box(lwd=3)
qqnorm(M)
qqline(M, col = 2)
dev.off()

#plot vs the isolation interference
I= sk$Isolation.Interference.in.Percent
I[is.na(I)]<-0
#
pdf('ch_preNormExpr_IsolationInterference_240inj.pdf')
heatscatter(I,
		M,
		cex=1,
		pch=19,
		xlab = 'Isolation Interference (percent)',
		ylab = 'Intensity Ratio')
abline(h=0,col=brewer.pal(6,"RdBu")[6],lwd=3)
box(lwd=3)
dev.off()

#plot the signal distribution to determine the cut off of noise
pdf('ch_preNormExpr_TotalExpression_240inj.pdf')
for (i in 7:16){
	hist(log10(sk[,i]),
			breaks=200,
			main=paste('total expression',colnames(sk)[i]),
			xlab='s/n value',
			xlim=c(0,3.5))
}
for (i in 7:16){
	l = sk[,i]<20 
	hist(sk[l,i],
			breaks=200,
			main=paste('total expression',colnames(sk)[i]),
			xlab='s/n value',
			xlim=c(0,20))
}
dev.off()

#pre-normalization
M= log2(sk[,7]) - log2(sk[,8])
A= log10(sk$Ion.Inject.Time.in.ms)
#plot
pdf('ch_preNormExpr_MvsInjTime_120inj.pdf')
heatscatter(A,
		M,
		cex=1,
		pch=19,
		xlab = 'log(ion inject time)',
		ylab = 'Intensity Ratio')
abline(h=0,col=brewer.pal(6,"RdBu")[6],lwd=3)
box(lwd=3)
heatscatter(A,
		log2(sk$X126),
		cex=1,
		pch=19,
		xlab = 'ion inject time',
		ylab = 's/n')
abline(h=0,col=brewer.pal(6,"RdBu")[6],lwd=3)
box(lwd=3)
dev.off()
