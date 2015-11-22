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
library(PECA)
###############################################################################
#grab the data
setwd(dir="/Users/cshughes/Documents/projects/OvC/")
ovpsm<-read.table("./ch_09July2015_OvC-TMT10_HpH_repAB_PSMs.txt", header=TRUE, sep='\t')
ovpep<-read.table("./ch_09July2015_OvC-TMT10_HpH_repAB_PeptideGroups.txt", header=TRUE, sep='\t')
ovpro<-read.table("./ch_09July2015_OvC-TMT10_HpH_repAB_Proteins.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/Routput/")

##################################################
#initial processing of the data to compile peptides and remove contaminants
##################################################
#dont split replicates yet or do any combination of them
processPSM <- function(psmFile, proteinFile, ... ){	
	require(limma)
	#initial processing to remove non-informative columns and rows from contaminant proteins or those with many missing values
	pep<-psmFile[,c(5,6,8,10,11,25:27,30,34,37:47)]
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
	#parse the peptide column for amino acids
	pep.m$Sequence = toupper(sub('.*?\\.(.*?)(\\..*|$)','\\1',pep.m$Annotated.Sequence))
	#filter information from Description
	pep.m$Descriptions = sub('(.*?)( OS=.*|$)','\\1',pep.m$Description)
	#remove specific proteins
	#pep.m = subset(pep.m, !grepl('Keratin',pep.m$Descriptions))
	pep.m = subset(pep.m, !grepl('sp',pep.m$Accession, ignore.case=FALSE))
	#filter the data columns
	print('removing contaminants')
	print(dim(pep.m))
	pep.r = pep.m[,c(1,24,26,25,3:4,11:22,7:10)]
	#assign the replicate
	repA = grepl('repA', pep.r$Spectrum.File)
	repB = grepl('repB', pep.r$Spectrum.File)
	pep.r$Rep = 'NA'
	pep.r[repA,23] = 'A'
	pep.r[repB,23] = 'B'
	#assign the organism
	ec = grepl('True', pep.r$Std)
	print(paste('number of e coli peptides = ',table(ec)[2],sep=''))
	hu = grepl('False', pep.r$Std)
	print(paste('number of human peptides = ',table(hu)[2],sep=''))
	pep.r$Source = 'NA'
	pep.r[ec,24] = 'ecoli'
	pep.r[hu,24] = 'human'
	#aggregate the PSMs into peptides
	print('aggregating peptides')
	pep.a = aggregate(cbind(X127N,X127C,X128N,X129C,X130N,X126,X128C,X130C,X131,X129N,Percolator.PEP)~Accession+Gene+Descriptions+Sequence+Rep+Source,data=pep.r,mean,na.action=na.pass,na.rm=TRUE)
	colnames(pep.a) = c('Accession','Gene','Descriptions','Sequence','Rep','Source','t1210a','t1246a','t1215a','t1493a','t1303a','t1811D','t376D','t1842b','t1768a','t1798a','PEP')
	print(dim(pep.a))
	#replace NaN with NA
	pep.a[,7:16] = round(pep.a[,7:16],2)
	pep.a[,7:16][is.na(pep.a[,7:16])]<-NA
	#output the data
	return(pep.a)	
}

ov.psm<-processPSM(ovpsm,ovpro)

##################################################
#do the peptide quant
##################################################
#split by replicates
repA = subset(ov.psm, grepl('A', ov.psm$Rep))
repB = subset(ov.psm, grepl('B', ov.psm$Rep))
#remove ecoli peps
repA = subset(repA, !grepl('ecoli', repA$Source))
repB = subset(repB, !grepl('ecoli', repB$Source))
#do the PECA analysis
exprsA = repA[,c(1:3,4,7:16)]
exprsB = repB[,c(1:3,4,7:16)]
exprs = merge(exprsA,exprsB,by=c('Accession','Gene','Descriptions','Sequence'))


group1<-colnames(exprsA)[10:14]
group2<-colnames(exprsA)[5:9]
xnorm = PECA_df(exprsA,'Gene',group1,group2,normalize='median',test='modt',type='median')
xnorm$Gene = rownames(xnorm)
rownames(xnorm) = NULL
exprsA$pepNum = 1
anno = aggregate(pepNum~Accession+Gene+Descriptions,data=exprsA,sum,na.action=na.pass,na.rm=TRUE)
proA = merge(anno,xnorm,by='Gene')

exprsB = repB[,c(1:3,4,7:16)]
group1<-colnames(exprsB)[10:14]
group2<-colnames(exprsB)[5:9]
xnorm = PECA_df(exprsB,'Gene',group1,group2,normalize='median',test='modt',type='median')
xnorm$Gene = rownames(xnorm)
rownames(xnorm) = NULL
exprsB$pepNum = 1
anno = aggregate(pepNum~Accession+Gene+Descriptions,data=exprsB,sum,na.action=na.pass,na.rm=TRUE)
proB = merge(anno,xnorm,by='Gene')


##########################################################
##PLOTTING
##########################################################

##########################################################
##make a plot for peptide metrics
##########################################################
#use proA as a starter for this
#subset out peptides with greater than 50 peptides assigned
vE = merge(proA,proB,by=c('Gene','Accession','Descriptions'),all=TRUE)
vE$pepSum = ifelse(vE$pepNum.x > vE$pepNum.y, vE$pepNum.x, vE$pepNum.y)
#this text gives you the plot...choose based on the set you want to plot...e.g. sum, repA, or repB
nPeps<-data.frame(prop.table(table(vE[,4]))*100)
up<-subset(nPeps, as.numeric(as.character(nPeps$Var1))>25)
dn<-subset(nPeps, as.numeric(as.character(nPeps$Var1))<=25)
#sum the number of proteins id with more than 50 peptides
x<-sum(up$Freq)
dn[26,1]<-'26'
dn[26,2]<-x
#
mainCol<-brewer.pal(6,"RdBu")
colA<-col2rgb(mainCol[6])
xLabels<-seq(1,26,1)
xValues<-seq(1,26,1)
#
mp<-barplot(dn[,2])
pdf('ch_OvC_PeptideIDMetrics_RepA_ctpTissue.pdf')
barplot(dn[,2],
		col=rgb(colA[1,],colA[2,],colA[3,],100,maxColorValue=255),
		main="Distribution of Peptide Identification Numbers", 
		xlab="Number of Unique Peptides",
		ylim=c(0,max(dn[,2])+2),
		xaxt="n"
)

box(lwd=3)
axis(1,mp,xLabels,las=2,lwd=2)
dev.off()

#show the line plot for both
pdf('ch_OvC_PeptideIDMetrics_RepAB-Lines_ctpTissue.pdf')
plot(dn[,2],
		type='l',
		col='blue'
)
nPeps<-data.frame(prop.table(table(vE[,16]))*100)
up<-subset(nPeps, as.numeric(as.character(nPeps$Var1))>25)
dn<-subset(nPeps, as.numeric(as.character(nPeps$Var1))<=25)
x<-sum(up$Freq)
dn[26,1]<-'26'
dn[26,2]<-x
lines(dn[,2], col='red')
box(lwd=3)
dev.off()

#lines(mp,dn[,2], col='red')



##########################################################
##number of peptides exceeding 0.05 p-value
##########################################################
pdf('ch_OvC_PVal-Distribution_RepAB_ctpTissue.pdf')
vE = merge(proA,proB,by=c('Gene','Accession','Descriptions'),all=TRUE)
vE$pepSum = ifelse(vE$pepNum.x > vE$pepNum.y, vE$pepNum.x, vE$pepNum.y)
vE$mPVal = rowMeans(vE[,c(9,15)], na.rm=TRUE)
vE[,17] = -log10(vE[,17])
vE = subset(vE, vE[,17] >= -log10(0.05))
x = hist(log2(vE[,16]),breaks=20,plot=FALSE)
mp<-barplot(x$density,plot=FALSE)
barplot(x$density,
		space=rep(0.1,17),
		main='distribution of significant hits', 
		xlab='log2(pepNum)',
		xaxt="n"
)
axis(1,x$breaks,seq(0,8.5,0.5),las=2,lwd=2)
dev.off()

##########################################################
##correlation between the two replicates
##########################################################
vE = merge(proA,proB,by=c('Gene','Accession','Descriptions'),all=TRUE)
pdf('ch_OvC_ReplicateCorrelation_ctpTissue.pdf')
heatscatter(vE$slr.x,
		vE$slr.y,
		main = 'replicate correlation',
		xlab = 'replicate A',
		ylab = 'replicate B',
		xlim = c(-4,4),
		ylim = c(-4,4))
abline(h=0,lwd=2,lty=2)
abline(v=0,lwd=2,lty=2)
abline(lm(vE$slr.y~vE$slr.x),lwd=2)
dev.off()











