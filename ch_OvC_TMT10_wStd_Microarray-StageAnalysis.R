# TODO: Add comment
# 
# Author: cshughes
###############################################################################
library(GEOquery)
library(gcrma)
library(genefilter)
library(hgu133plus2probe)
library(hgu133plus2.db)
library(simpleaffy)
library(affyPLM)
library(annotate)
###############################################################################

#build the RNA set to be analyzed
setwd(dir="/Users/cshughes/Documents/projects/OvC/RNAexpression") #change this to whatever directory you have stored the data in
orna = readRDS('ch_OvC_RNA_processedProbes.rds')
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput")
#get out the expression data
eSet = as.data.frame(exprs(orna$eset))
#calculate RLE
eSet$med = apply(eSet,1,function(x) median(x,na.rm=TRUE))
eRLE = as.data.frame(apply(eSet[,1:55],2,function(x) x - eSet$med))
#annotate the data
gene.symbols <- getSYMBOL(row.names(eRLE), "hgu133plus2")
eRLE$Gene = gene.symbols
#remove probes with NA annotation
eSub = subset(eRLE, !is.na(eRLE$Gene))
#aggregate into single gene measurements
eAgg = aggregate(eSub[,1:55],by=list(eSub$Gene), median, na.rm=TRUE)
colnames(eAgg)[1] = 'Gene'

#get out genes you are interested in
gn = c('FOLR1','CRIP1','MSLN','SNCG','CRABP2','LEFTY1','GDF15','QPCT','GPC3','CTH')
#scan through genes to get the expression values out
gnex = data.frame(nrow=55)
for (i in 1:length(gn)){
	eset = t(eAgg[grepl(paste('\\<',gn[i],'\\>',sep=''), eAgg$Gene),2:56])
	colnames(eset) = gn[i]
	gnex = cbind(gnex,eset)
}

#grab the stage data
setwd(dir="/Users/cshughes/Documents/projects/OvC/RNAexpression/")
pData<-read.table("./Uehara_PatientData.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/RNAexpression/Routput")

#bind the data together
gp = cbind(pData,gnex)

#plot the data
#make the colors for the plot
cols1 = brewer.pal(6,'Greens')
cols2 = col2rgb(brewer.pal(6,'Greens'))
#make the plot
pdf('ch_test.pdf')
for (i in 8:17){
	boxplot(gp[grepl('Serous',gp$Histology),i]~gp[grepl('Serous',gp$Histology),4],las=2,main=colnames(gp)[i],border=cols1,outline=FALSE,boxlwd=2,col = rgb(cols2[1,],cols2[2,],cols2[3,],25,maxColorValue=255))
	beeswarm(gp[grepl('Serous',gp$Histology),i]~gp[grepl('Serous',gp$Histology),4],pch=16,col='black',corral="omit",add=TRUE,cex=0.95)
}
dev.off()

pdf('ch_test.pdf')
for (i in 8:17){
	plot(gp[1:36,i],gp[1:36,5])
}
dev.off()

#1:25 is CCC
#40:55 is HGSC

pout = ''
pcor = 0
pall = data.frame()
z=0
gn = as.vector(proX$Gene)
for (i in 1:length(gn)){
	eset = t(eAgg[grepl(paste('\\<',gn[i],'\\>',sep=''), eAgg$Gene),2:56])
	pset = cbind(pData,eset)
	if (gn[i] %in% eAgg$Gene){
		epcor = cor(pset[40:55,7],log2(pset[40:55,3]))
		if (abs(epcor)>0.25){
			z = z +1
			pall[z,1] = gn[i]
			pall[z,2] = cor(pset[40:55,7],log2(pset[40:55,3]))
		}
	}	
}
#reannotate the data frame
colnames(pall) = c('Gene','Corr')
pall = pall[order(-pall$Corr),]


gnex = data.frame(nrow=55)
for (i in 2:length(pall$Gene)){
	eset = as.data.frame(t(eAgg[grepl(paste('\\<',pall[i,1],'\\>',sep=''), eAgg$Gene),2:56]))
	colnames(eset) = pall[i,1]
	gnex = cbind(gnex,eset)
}


pdf('ch_test.pdf')
boxplot(gnex[40:55,130:179],las=2)
beeswarm(gnex[40:55,130:179],pch=16,col='black',corral="omit",add=TRUE,cex=0.5)	
dev.off()


