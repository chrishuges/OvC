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


