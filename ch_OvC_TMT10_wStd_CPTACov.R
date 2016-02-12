# TODO: Add comment
# 
# Author: cshughes
###############################################################################
setwd(dir="/Users/cshughes/Documents/projects/cptac/ovarian") #change this to whatever directory you have stored the data in
jhu<-read.table("./TCGA_Ovarian_JHU_Proteome_CDAP.r2.itraq.txt", header=TRUE, sep='\t', na.strings='')
pnnl<-read.table("./TCGA_Ovarian_PNNL_Proteome_CDAP.r2.itraq.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput") #change this to whatever directory you have stored the data in

#get rid of the unshared data
jhu = cbind('Gene'=jhu$Gene, jhu[,which(grepl('Unshared',colnames(jhu)))])
pnnl = cbind('Gene'=pnnl$Gene, pnnl[,which(grepl('Unshared',colnames(pnnl)))])
#get rid of the normal ovarian tissue samples
#jhu = jhu[,which(!grepl('CONTROL',colnames(jhu)))]
#pnnl = pnnl[,which(!grepl('CONTROL',colnames(pnnl)))]
#merge into a single set
cptac = merge(jhu,pnnl,by='Gene',all=TRUE)

#make a plot to compare with our markers for HGS
library(beeswarm)
gn = c('FOLR1','CRIP1','MSLN','SNCG','CRABP2','LEFTY1','GDF15','QPCT','GPC3','CTH')
#map the ovarian study
cpPro = as.data.frame(t(cptac[which(cptac$Gene %in% gn),2:207]))
colnames(cpPro) = cptac[which(cptac$Gene %in% gn),1]
cpPro = cpPro[,c(4,2,6,8,1,9,7,5,3)]
#make the colors for the plot
cols1 = c(rep(brewer.pal(9,'RdBu')[1],5),rep(brewer.pal(9,'RdBu')[9],4))
cols2 = col2rgb(c(rep(brewer.pal(9,'RdBu')[1],5),rep(brewer.pal(9,'RdBu')[9],4)))
#plot the data
pdf('ch_OvC_TMT10_wStd_Human_HGSMarkers_OvC-CPTAC_boxplot.pdf')
boxplot(cpPro,
		las=2,
		border=cols1,
		col = rgb(cols2[1,],cols2[2,],cols2[3,],25,maxColorValue=255),
		outline = FALSE,
		boxlwd = 2,
		ylab = 'log2(CPTAC Expression)',
		main = 'Expression of HGS Markers in Ovarian CPTAC Data'
)
beeswarm(cpPro,pch=16,col='black',corral="omit",add=TRUE,cex=0.35)
dev.off()



#this is maybe cheating a bit above...spectral counts is maybe a better number to use because the itraq value is relative to a pool of HGS samples
###############################################################################
setwd(dir="/Users/cshughes/Documents/projects/cptac/ovarian") #change this to whatever directory you have stored the data in
jhu<-read.table("./TCGA_Ovarian_JHU_Proteome_CDAP.r2.summary.txt", header=TRUE, sep='\t', na.strings='')
pnnl<-read.table("./TCGA_Ovarian_PNNL_Proteome_CDAP.r2.summary.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput") #change this to whatever directory you have stored the data in

#get rid of the unshared data
jhu = cbind('Gene'=jhu$Gene, jhu[,which(grepl('Spectral',colnames(jhu)))])
pnnl = cbind('Gene'=pnnl$Gene, pnnl[,which(grepl('Spectral',colnames(pnnl)))])
#get rid of the normal ovarian tissue samples
jhu = jhu[,which(!grepl('CONTROL',colnames(jhu)))]
pnnl = pnnl[,which(!grepl('CONTROL',colnames(pnnl)))]
#merge into a single set
cptac = merge(jhu,pnnl,by='Gene',all=TRUE)

#make a plot to compare with our markers for HGS
library(beeswarm)
gn = c('FOLR1','CRIP1','MSLN','SNCG','CRABP2','LEFTY1','GDF15','QPCT','GPC3','CTH')
#map the ovarian study
cpPro = as.data.frame(t(cptac[which(cptac$Gene %in% gn),c(2:35,37:64)]))
colnames(cpPro) = cptac[which(cptac$Gene %in% gn),1]
cpPro = cpPro[,c(4,2,6,8,1,9,7,5,3)]
#make the colors for the plot
cols1 = c(rep(brewer.pal(9,'RdBu')[1],5),rep(brewer.pal(9,'RdBu')[9],4))
cols2 = col2rgb(c(rep(brewer.pal(9,'RdBu')[1],5),rep(brewer.pal(9,'RdBu')[9],4)))
#plot the data
pdf('ch_OvC_TMT10_wStd_Human_HGSMarkers_OvC-CPTAC-sc_boxplot.pdf')
boxplot(cpPro,
		las=2,
		border=cols1,
		col = rgb(cols2[1,],cols2[2,],cols2[3,],25,maxColorValue=255),
		outline = FALSE,
		boxlwd = 2,
		ylab = 'log2(CPTAC Expression)',
		main = 'Expression of HGS Markers in Ovarian CPTAC Data'
)
beeswarm(cpPro,pch=16,col='black',corral="omit",add=TRUE,cex=0.35)
dev.off()
#this is more as expected...but highlights how the log2 to pool is maybe not that useful in CPTAC data



###############################################################################
#Make a volcano plot of my data vs cptac high and low expression
###############################################################################
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput") #change this to whatever directory you have stored the data in
phvc = readRDS('ch_OvC_wStd_Proteins_HGSvCCC.rds')

#process the cptac data to extract high and low expressed proteins
setwd(dir="/Users/cshughes/Documents/projects/cptac/ovarian") #change this to whatever directory you have stored the data in
jhu<-read.table("./TCGA_Ovarian_JHU_Proteome_CDAP.r2.summary.txt", header=TRUE, sep='\t', na.strings='')
pnnl<-read.table("./TCGA_Ovarian_PNNL_Proteome_CDAP.r2.summary.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput") #change this to whatever directory you have stored the data in

#get rid of the unshared data
jhu = cbind('Gene'=jhu$Gene, jhu[,which(grepl('Spectral',colnames(jhu)))])
pnnl = cbind('Gene'=pnnl$Gene, pnnl[,which(grepl('Spectral',colnames(pnnl)))])
#get rid of the normal ovarian tissue samples
jhu = jhu[,which(!grepl('CONTROL',colnames(jhu)))]
pnnl = pnnl[,which(!grepl('CONTROL',colnames(pnnl)))]
#merge into a single set
cptac = merge(jhu,pnnl,by='Gene',all=TRUE)
cptac$mSC = apply(cptac[,c(36,65)],1,function(x) median(x,na.rm=TRUE))
#get the top 50 expressed genes
cptac = cptac[order(-cptac$mSC),]
gnUP = as.character(cptac[1:65,1]) #this yields 50 genes with an expression value
#get the bottom 50 expressed genes
cptac = cptac[order(cptac$mSC),]
#cSet = subset(cptac, mSC>0) #dont keep genes with 0 spectral counts as a median
gnDN = as.character(cptac[1:125,1]) #this yields 50 genes with an expression value
#total list of genes
gnALL = c(gnUP,gnDN)

#bring in my protein data
proh = phvc
#only keep markers that are present in the data
markHGSCCC = gnALL[which(gnALL %in% proh$Gene)]
markHGS = gnUP[which(gnUP %in% proh$Gene)]
markCCC = gnDN[which(gnDN %in% proh$Gene)]

#make a plot of the data
proSD<-sd(proh$PROexp, na.rm=TRUE)
#sort out colors
lnCols<-brewer.pal(6,"Blues")
cols<-rev(brewer.pal(6,"RdBu"))

###make an initial plot with all points
pdf('ch_OvC_TMT10_wStd_Human_Proteins_CPTAC-HGSvCCC_Volcano.pdf')
xCol = col2rgb(ifelse(proh$Gene %in% markCCC, cols[1], ifelse(proh$Gene %in% markHGS, cols[6],'gray80')))
xCex = ifelse(proh$Gene %in% markHGSCCC, 2, 1)
plot(proh$PROexp,
		-log10(proh$score),
		col=rgb(xCol[1,],xCol[2,],xCol[3,],95,maxColorValue=255),
		pch=20,
		cex = xCex,
		ylab = '-log10(Adjusted p-value)',
		xlab = 'log2(serous/clear cell)',
		main = 'OvC type comparison',
		xlim = c(-5,5),
		ylim = c(0,5)
)
box(lwd=3)
abline(v = -proSD, col=lnCols[6],lwd=2,lty=2)
abline(v = proSD, col=lnCols[6],lwd=2,lty=2)
abline(h = -log10(0.05), col=lnCols[6],lwd=2,lty=2)
text(4,0.5,paste('n=',nrow(proh),sep=""),cex=1.25)
text(4,2.3,paste('p<0.05'),cex=1.25)

###make a second plot with just the marker points
xCol = ifelse(proh$Gene %in% markCCC, cols[1], ifelse(proh$Gene %in% markHGS, cols[6],'gray60'))
gnRM = !grepl('gray', xCol)
proh.s = proh[gnRM,]
xCol = col2rgb(ifelse(proh.s$Gene %in% markCCC, cols[1], ifelse(proh.s$Gene %in% markHGS, cols[6],'gray80')))
plot(proh.s$PROexp,
		-log10(proh.s$score),
		col=rgb(xCol[1,],xCol[2,],xCol[3,],95,maxColorValue=255),
		pch=20,
		cex = 2,
		ylab = '-log10(Adjusted p-value)',
		xlab = 'log2(serous/clear cell)',
		main = 'OvC type comparison',
		xlim = c(-5,5),
		ylim = c(0,5)
)
box(lwd=3)
abline(v = -proSD, col=lnCols[6],lwd=2,lty=2)
abline(v = proSD, col=lnCols[6],lwd=2,lty=2)
abline(h = -log10(0.05), col=lnCols[6],lwd=2,lty=2)
text(4,0.5,paste('n=',sum(gnRM),sep=""),cex=1.25)
text(4,2.3,paste('p<0.05'),cex=1.25)
text(proh.s$PROexp, -log10(proh.s$score), proh.s$Gene,cex=0.25)
dev.off()





###############################################################################
#Query the high low genes by subtype in Breast Cancer
###############################################################################
#process the cptac data to extract high and low expressed proteins
setwd(dir="/Users/cshughes/Documents/projects/cptac/breast") #change this to whatever directory you have stored the data in
broad<-read.table("./TCGA_Breast_BI_Proteome_CDAP.r2.itraq.txt", header=TRUE, sep='\t')
broadAnno<-read.table("./CPTAC_TCGA_BreastCancer_select_clinical_data_r1.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput") #change this to whatever directory you have stored the data in
#process the CPTAC data
#get the IDs for basal bc only
bas = broadAnno[grepl('Basal',broadAnno$PAM50.mRNA),1]
bas2 = substr(bas,6,13)
#process the expression
#need to parse the column names
bcol = sub('(.*?)(\\.01A.*|$)','\\1',colnames(broad))
bcol = gsub("\\.", "-", bcol)
colnames(broad) = bcol
#get rid of the non-basal samples
broadt = cbind('Gene'=broad$Gene, broad[,which(colnames(broad) %in% bas2)])
#get the median expression
broadt$mSC = apply(broadt[,c(2:27)],1,function(x) median(x,na.rm=TRUE))
#get the proteomic data
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput") #change this to whatever directory you have stored the data in
a1 = readRDS('ch_OvC_wStd_processedPeptides_a1.rds')
a2 = readRDS('ch_OvC_wStd_processedPeptides_a2.rds')
a3 = readRDS('ch_OvC_wStd_processedPeptides_a3.rds')
b1 = readRDS('ch_OvC_wStd_processedPeptides_b1.rds')
b2 = readRDS('ch_OvC_wStd_processedPeptides_b2.rds')
b3 = readRDS('ch_OvC_wStd_processedPeptides_b3.rds')
#use this function for normalization
pepNorm = function(x,...){
	#mark the spike-in control peptides
	spikeins = grepl('EColi', x$Organism)
	#extract the expression values
	xnorm = as.matrix(x[,7:16])
	#make a noise model based on the spike ins
	spfit = vsn2(xnorm[spikeins,],lts.quantile=1)
	meanSdPlot(spfit)
	#apply it to the human peptides
	nkid = predict(spfit, newdata=xnorm)
	#extract only human expression values
	nkid.h = as.data.frame(nkid[!spikeins,])
	#normalize the values to the internal standard sample
	vnorm = apply(nkid.h[,1:9],2, function(x) x - nkid.h$aStd)
	snorm = scale(vnorm, center=TRUE,scale=TRUE)
	#recombine the replicates
	vnormFull = as.data.frame(cbind(x[!spikeins,c(1:6,17)],snorm))
	#output the data
	return(vnormFull)
}
#apply the function
a1n = pepNorm(a1)
a2n = pepNorm(a2)
a3n = pepNorm(a3)
b1n = pepNorm(b1)
b2n = pepNorm(b2)
b3n = pepNorm(b3)
#merge replicates
z1 = rbind(a1n,a2n,a3n)
z2 = rbind(b1n,b2n,b3n)
#aggregate redundant peptides
agg1 = aggregate(cbind(a1,a2,a3,a4,a5,a6,a7,a8,a9)~Gene,data=z1,median,na.action=na.pass,na.rm=TRUE)
agg2 = aggregate(cbind(a1,a2,a3,a4,a5,a6,a7,a8,a9)~Gene,data=z2,median,na.action=na.pass,na.rm=TRUE)
#merge into a gene set
proS = merge(agg1,agg2,by='Gene')
#calculate mean expression for the subtypes
proS$hgs1 = rowMeans(proS[,c(2:4,11:13)],na.rm=TRUE)
proS$ccc1 = rowMeans(proS[,c(5:7,14:16)],na.rm=TRUE)
proS$emc1 = rowMeans(proS[,c(8:10,17:19)],na.rm=TRUE)

#extract the spectral counts in the breast cancer data
gOut = data.frame(Gene=character(),exp=numeric(),anno=character())
qTCGA = function(dataIn,hgsLoc,cccLoc,emcLoc,cptacD,...){
	x = subset(dataIn, !grepl('Homo sapiens', dataIn$Gene))
	set = c(hgsLoc,cccLoc,emcLoc)
	for (l in 1:length(set)){
		#start with DN genes
		xS = x[order(x[,set[l]]),]
		gnDN = xS[1:100,1]
		pSet = cptacD[cptacD$Gene %in% gnDN,c(1,28)]
		pSet$anno = paste(colnames(x)[set[l]],'DN',sep='')
		gOut = rbind(gOut,pSet)
		#do again for UP genes
		xS = x[order(-x[,set[l]]),]
		gnUP = xS[1:100,1]
		pSet = cptacD[cptacD$Gene %in% gnUP,c(1,28)]
		pSet$anno = paste(colnames(x)[set[l]],'UP',sep='')
		gOut = rbind(gOut,pSet)
	}
	#output the data
	return(gOut)	
}
gnSet = qTCGA(proS,20,21,22,broadt)



#plot the data
cols1 = c(brewer.pal(9,'RdBu')[9],brewer.pal(9,'RdBu')[1],brewer.pal(9,'RdBu')[9],brewer.pal(9,'RdBu')[1],brewer.pal(9,'RdBu')[9],brewer.pal(9,'RdBu')[1])
cols2 = col2rgb(c(brewer.pal(9,'RdBu')[9],brewer.pal(9,'RdBu')[1],brewer.pal(9,'RdBu')[9],brewer.pal(9,'RdBu')[1],brewer.pal(9,'RdBu')[9],brewer.pal(9,'RdBu')[1]))
#plot the data
pdf('ch_OvC_TMT10_wStd_Human_HGSupdn_BrC-cptac_boxplot.pdf')
boxplot(gnSet$mSC~gnSet$anno,
		las=2,
		border=cols1,
		col = rgb(cols2[1,],cols2[2,],cols2[3,],25,maxColorValue=255),
		outline = FALSE,
		boxlwd = 2,
		ylab = 'log2(TCGA RNA FPKM)',
		main = 'Expression of HGS Markers in Breast TCGA Data'
)
beeswarm(gnSet$mSC~gnSet$anno,pch=16,col='black',corral="omit",add=TRUE,cex=0.65)
dev.off()

#get the p-values
up = gnSet[grepl('hgs1UP',gnSet$anno),]
dn = gnSet[grepl('hgs1DN',gnSet$anno),]
#sig test
pSig<-wilcox.test(dn$mSC,up$mSC)
#hgs p-value = 0.1105
#ccc p-value = 3.36e-7
#emc p-value = 0.6529








