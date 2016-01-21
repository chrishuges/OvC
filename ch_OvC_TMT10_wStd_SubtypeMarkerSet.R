# TODO: Add comment
# 
# Author: cshughes
###############################################################################
#read in data
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput") #change this to whatever directory you have stored the data in
phvc = readRDS('ch_OvC_wStd_Proteins_HGSvCCC.rds')
phve = readRDS('ch_OvC_wStd_Proteins_HGSvEMC.rds')
pcve = readRDS('ch_OvC_wStd_Proteins_CCCvEMC.rds')




pSet = merge(hvc,hve,by='Gene')
pdf('ch_OvC_wStd_Proteins_HGSMarkers_Scatter.pdf')
cols<-col2rgb(rev(brewer.pal(6,"RdBu")))
plot(pSet$PROexp.x,pSet$PROexp.y,pch=20,col=rgb(cols[1,1],cols[2,1],cols[3,1],95,maxColorValue=255),cex=0.5)
text(1,-3,paste('r = ',cor(pSet$PROexp.x,pSet$PROexp.y,use='pairwise.complete.obs'),sep=''))
text(pSet$PROexp.x,pSet$PROexp.y,pSet$Gene,cex=0.25)
abline(lm(pSet$PROexp.y~pSet$PROexp.x),lwd=2,lty=2)
box(lwd=3)
abline(h=0,lwd=2,lty=2)
abline(v=0,lwd=2,lty=2)
dev.off()



###############################################################################
#work for HGS data
###############################################################################
###############################################################################
#this does 1 gene across all cancers
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
#this does multiple genes across ovarian cancer
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








