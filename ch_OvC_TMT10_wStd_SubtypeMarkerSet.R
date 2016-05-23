# TODO: Add comment
# 
# Author: cshughes
###############################################################################
#read in data
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput") #change this to whatever directory you have stored the data in
phvc = readRDS('ch_OvC_wStd_Proteins_HGSvCCC.rds')
phve = readRDS('ch_OvC_wStd_Proteins_HGSvEMC.rds')
pcve = readRDS('ch_OvC_wStd_Proteins_CCCvEMC.rds')



#this is just a basic plot that will allow you to extract subtype markers by comparing with the two others
pSet = merge(phvc,phve,by='Gene')
pdf('ch_OvC_wStd_Proteins_HGSMarkers_Scatter.pdf')
cols<-col2rgb(rev(brewer.pal(6,"RdBu")))
plot(pSet$PROexp.x,pSet$PROexp.y,pch=20,col=rgb(cols[1,1],cols[2,1],cols[3,1],95,maxColorValue=255),cex=0.5)
text(1,-3,paste('r = ',cor(pSet$PROexp.x,pSet$PROexp.y,use='pairwise.complete.obs'),sep=''))
text(pSet$PROexp.x,pSet$PROexp.y,pSet$Gene,cex=0.05)
#abline(lm(pSet$PROexp.y~pSet$PROexp.x),lwd=2,lty=2)
box(lwd=3)
abline(h=0,lwd=2,lty=2)
abline(v=0,lwd=2,lty=2)
dev.off()


###get marker proteins
pSet = merge(phvc,phve,by='Gene')

exp = subset(phvc, PROexp<0.5 & p.fdr<=0.05)
write.table(exp,'ch_OvC_wStd_HGSvCCC_CCCup_proteinSet.txt',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)






#merge the protein and RNA data
hvc = merge(phvc,phve,by='Gene')
#make the plot for all three sets
#bring in marker data
setwd(dir="/Users/cshughes/Documents/projects/OvC/markers/")
genesCC = read.table("./Schwartz_ClearCellTable_geneList.txt", header=FALSE, sep='\t')
genesS = read.table("./Schwartz_SerousTable_geneList.txt", header=FALSE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput/")####CHANGE ME
#combine into a single set
geneCC = as.character(genesCC[,1])
geneS = as.character(genesS[,1])
geneTot = c(geneCC,geneS)
#geneTot = c('CTH','LEFTY1')
#only keep markers that are present in the data
markHGSCCC = geneTot[which(geneTot %in% hvc$Gene)]
markHGS = geneS[which(geneS %in% hvc$Gene)]
markCCC = geneCC[which(geneCC %in% hvc$Gene)]
#make the plot
pdf('ch_OvC_wStd_Proteins_HGSMarkers_Scatter.pdf')
cols<-rev(brewer.pal(6,"RdBu"))
xCol = col2rgb(ifelse(hvc$Gene %in% markCCC, cols[1], ifelse(hvc$Gene %in% markHGS, cols[6],'gray60')))
xCex = ifelse(hvc$Gene %in% markHGSCCC, 2, 1)
plot(hvc$PROexp.x,hvc$PROexp.y,pch=20,col=rgb(xCol[1,],xCol[2,],xCol[3,],95,maxColorValue=255),cex=xCex,xlim = c(-4,4),ylim = c(-4,4))
text(1,-3,paste('r = ',cor(hvc$PROexp.x,hvc$PROexp.y,use='pairwise.complete.obs'),sep=''))
#text(hvc$PROexp,hvc$RNAexp,hvc$Gene,cex=0.25)
abline(lm(hvc$PROexp.y~hvc$PROexp.x),lwd=2,lty=2)
abline(h=0,lwd=2,lty=2)
abline(v=0,lwd=2,lty=2)
box(lwd=3)
xCol = ifelse(hvc$Gene %in% markCCC, cols[1], ifelse(hvc$Gene %in% markHGS, cols[6],'gray60'))
gnRM = !grepl('gray', xCol)
hvc.s = hvc[gnRM,]
xCol = col2rgb(ifelse(hvc.s$Gene %in% markCCC, cols[1], ifelse(hvc.s$Gene %in% markHGS, cols[6],'gray60')))
xCex = ifelse(hvc.s$Gene %in% markHGSCCC, 2, 1)
plot(hvc.s$PROexp.x,hvc.s$PROexp.y,pch=20,col=rgb(xCol[1,],xCol[2,],xCol[3,],95,maxColorValue=255),cex=xCex,xlim = c(-4,4),ylim = c(-4,4))
text(hvc.s$PROexp.x,hvc.s$PROexp.y,hvc.s$Gene,cex=0.5)
text(1,-3,paste('r = ',cor(hvc.s$PROexp.x,hvc.s$PROexp.y,use='pairwise.complete.obs'),sep=''))
dev.off()




###############################################################################
#work for HGS data
###############################################################################
###############################################################################
#this does 1 gene across all cancers in TCGA
###############################################################################

#pull in the TCGA data
library(cgdsr)
mycgds = CGDS("http://www.cbioportal.org/public-portal/")
getCancerStudies(mycgds)[,2]
#cStudies = c(17,77,27,40,97,62,94,51,7,118,88)
#ALWAYS CHECK THESE NUMBERS AS THE LIST SEEMS TO BE FOREVER CHANGING
cStudies = c(78,18,91,28,53,100,80,60,11,56,116)
mRNA = list()
gn = 'MSLN'
for (i in 1:length(cStudies)){
	mycancerstudy = getCancerStudies(mycgds)[cStudies[i],1]
	mycaselist = getCaseLists(mycgds,mycancerstudy)[1,1]
	eDataLoc = which(grepl('rna_seq_v2_mrna$', getGeneticProfiles(mycgds,mycancerstudy)[,1]))[1]
	eDataSel = getGeneticProfiles(mycgds,mycancerstudy)[eDataLoc,1]
	eData = getProfileData(mycgds,gn,eDataSel,mycaselist)
	eData$set = mycancerstudy
	mRNA[[i]] = eData
}
#coolapse into a single frame
tRNA = do.call(rbind, mRNA)
tRNA.s = subset(tRNA, MSLN>0)
cols1 = c(brewer.pal(11,'Set3'))
cols2 = col2rgb(c(brewer.pal(11,'Set3')))
#plot the data
pdf('ch_OvC_TMT10_wStd_Human_HGSMarkers_MSLN_All-TCGA_boxplot.pdf')
boxplot(log2(tRNA.s$MSLN)~tRNA.s$set,
		las=2,
		col = rgb(cols2[1,],cols2[2,],cols2[3,],25,maxColorValue=255),
		border=cols1,
		boxlwd = 2,
		outline=FALSE,
		ylab = 'log2(TCGA RNA FPKM)',
		main = 'Expression of HGS Markers in All TCGA Data'
)
abline(h=median(log2(tRNA.s$MSLN),na.rm=TRUE),lty=2,lwd=3)
beeswarm(log2(tRNA.s$MSLN)~tRNA.s$set,pch=16,col='black',corral="omit",add=TRUE,cex=0.35)
dev.off()


###############################################################################
#this does multiple genes across ovarian cancer in TCGA
###############################################################################
library(beeswarm)
gn = c('FOLR1','CRIP1','MSLN','SNCG','CRABP2','LEFTY1','GDF15','QPCT','GPC3','CTH')
#map the ovarian study
cStudies = 78
mRNA = list()
for (i in 1:length(cStudies)){
	mycancerstudy = getCancerStudies(mycgds)[cStudies[i],1]
	mycaselist = getCaseLists(mycgds,mycancerstudy)[1,1]
	eDataLoc = which(grepl('rna_seq_v2_mrna$', getGeneticProfiles(mycgds,mycancerstudy)[,1]))[1]
	eDataSel = getGeneticProfiles(mycgds,mycancerstudy)[eDataLoc,1]
	eData = getProfileData(mycgds,gn,eDataSel,mycaselist)
	mRNA[[i]] = eData
}
#collapse into a single frame
tRNA = do.call(rbind, mRNA)
tRNA[tRNA<1]<-NA
tRNA.s = tRNA[,c(4,2,8,10,1,7,5,9,6,3)]
cols1 = c(rep(brewer.pal(9,'RdBu')[1],5),rep(brewer.pal(9,'RdBu')[9],5))
cols2 = col2rgb(c(rep(brewer.pal(9,'RdBu')[1],5),rep(brewer.pal(9,'RdBu')[9],5)))
#plot the data
pdf('ch_OvC_TMT10_wStd_Human_HGSMarkers_OvC-TCGA_boxplot.pdf')
boxplot(log2(tRNA.s),
		las=2,
		border=cols1,
		col = rgb(cols2[1,],cols2[2,],cols2[3,],25,maxColorValue=255),
		outline = FALSE,
		boxlwd = 2,
		ylab = 'log2(TCGA RNA FPKM)',
		main = 'Expression of HGS Markers in Ovarian TCGA Data'
)
beeswarm(log2(tRNA.s),pch=16,col='black',corral="omit",add=TRUE,cex=0.35)
dev.off()




###############################################################################
#compare HGS genes with the breast cancer data in TCGA
###############################################################################
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
###############################################################################
#query high and low expressed genes against MSigDB
###############################################################################
hgs = proS[,c(1,20)]
ccc = proS[,c(1,21)]
emc = proS[,c(1,22)]
#function for getting the gene expression values for all subtypes from the TCGA
gOut = data.frame(Gene=character(),exp=numeric(),anno=character())
qTCGA = function(dataIn,hgsLoc,cccLoc,emcLoc,...){
	x = subset(dataIn, !grepl('Homo sapiens', dataIn$Gene))
	set = c(hgsLoc,cccLoc,emcLoc)
	for (l in 1:length(set)){
		#start with DN genes
		xS = x[order(x[,set[l]]),]
		gnDN = xS[1:100,c(1,set[l])]
		colnames(gnDN)[2] = 'exp'
		#query TCGA
		cStudies = 16
		mRNA = list()
		gn = gnDN[,1]
		for (i in 1:length(cStudies)){
			mycancerstudy = getCancerStudies(mycgds)[cStudies[i],1]
			mycaselist = getCaseLists(mycgds,mycancerstudy)[7,1]
			eDataLoc = which(grepl('rna_seq_v2_mrna$', getGeneticProfiles(mycgds,mycancerstudy)[,1]))[1]
			eDataSel = getGeneticProfiles(mycgds,mycancerstudy)[eDataLoc,1]
			eData = getProfileData(mycgds,gn,eDataSel,mycaselist)
			mRNA[[i]] = eData
		}
		#collapse into a single frame
		tRNA = do.call(rbind, mRNA)
		tRNA[tRNA<1]<-NA
		tcDN = as.data.frame(t(tRNA))
		tcDN$exp = apply(tcDN,1,function(x) median(x,na.rm=TRUE))
		tcDN$Gene = row.names(tcDN)
		tcDN$anno = paste(colnames(x)[set[l]],'DN',sep='')
		gOut = rbind(gOut,tcDN[,c(109,108,110)])
		#do again for UP genes
		xS = x[order(-x[,set[l]]),]
		gnUP = xS[1:100,c(1,set[l])]
		colnames(gnUP)[2] = 'exp'
		#query TCGA
		cStudies = 16
		mRNA = list()
		gn = gnUP[,1]
		for (i in 1:length(cStudies)){
			mycancerstudy = getCancerStudies(mycgds)[cStudies[i],1]
			mycaselist = getCaseLists(mycgds,mycancerstudy)[7,1]
			eDataLoc = which(grepl('rna_seq_v2_mrna$', getGeneticProfiles(mycgds,mycancerstudy)[,1]))[1]
			eDataSel = getGeneticProfiles(mycgds,mycancerstudy)[eDataLoc,1]
			eData = getProfileData(mycgds,gn,eDataSel,mycaselist)
			mRNA[[i]] = eData
		}
		#collapse into a single frame
		tRNA = do.call(rbind, mRNA)
		tRNA[tRNA<1]<-NA
		tcUP = as.data.frame(t(tRNA))
		tcUP$exp = apply(tcUP,1,function(x) median(x,na.rm=TRUE))
		tcUP$Gene = row.names(tcUP)
		tcUP$anno = paste(colnames(x)[set[l]],'UP',sep='')
		gOut = rbind(gOut,tcUP[,c(109,108,110)])
	}
	#output the data
	return(gOut)	
}
gnSet = qTCGA(proS,20,21,22)


#plot the data
cols1 = c(brewer.pal(9,'RdBu')[9],brewer.pal(9,'RdBu')[1],brewer.pal(9,'RdBu')[9],brewer.pal(9,'RdBu')[1],brewer.pal(9,'RdBu')[9],brewer.pal(9,'RdBu')[1])
cols2 = col2rgb(c(brewer.pal(9,'RdBu')[9],brewer.pal(9,'RdBu')[1],brewer.pal(9,'RdBu')[9],brewer.pal(9,'RdBu')[1],brewer.pal(9,'RdBu')[9],brewer.pal(9,'RdBu')[1]))
#plot the data
pdf('ch_OvC_TMT10_wStd_Human_HGSupdn_BrC-TCGAcell_boxplot.pdf')
boxplot(log2(gnSet$exp)~gnSet$anno,
		las=2,
		border=cols1,
		col = rgb(cols2[1,],cols2[2,],cols2[3,],25,maxColorValue=255),
		outline = FALSE,
		boxlwd = 2,
		ylab = 'log2(TCGA RNA FPKM)',
		main = 'Expression of HGS Markers in Breast TCGA Data'
)
beeswarm(log2(gnSet$exp)~gnSet$anno,pch=16,col='black',corral="omit",add=TRUE,cex=0.65)
dev.off()

#get the p-values
up = gnSet[grepl('emc1UP',gnSet$anno),]
dn = gnSet[grepl('emc1DN',gnSet$anno),]
#sig test
pSig<-wilcox.test(log2(dn$exp),log2(up$exp))
#hgs p-value = 0.0001614
#ccc p-value = 0.9048
#emc p-value = 0.004959










