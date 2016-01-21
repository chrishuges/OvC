# TODO: Add comment
# 
# Author: cshughes
###############################################################################
#get the protein data
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput") #change this to whatever directory you have stored the data in
pro = readRDS('ch_OvC_wStd_RLEset.rds')
#calculate RLE
pro$med = apply(pro[,2:19],1,function(x) median(x,na.rm=TRUE))
pro[,2:19] = apply(pro[,2:19],2,function(x) x - pro$med)
#mean RLE for HGS
pro$hgs = rowMeans(pro[,2:7],na.rm=TRUE)
pro = pro[order(-abs(pro$hgs)),]
proS = subset(pro, !grepl('Homo sapiens', pro$Gene))
proS = subset(proS, hgs > 0.5 | hgs < -0.5)

#get the methylation data
library(cgdsr)
mycgds = CGDS("http://www.cbioportal.org/public-portal/")
getCancerStudies(mycgds)[,2]
#map the ovarian study
cStudies = 79
gn = proS$Gene[1:500]
mRNA = list()
for (i in 1:length(cStudies)){
	mycancerstudy = getCancerStudies(mycgds)[cStudies[i],1]
	mycaselist = getCaseLists(mycgds,mycancerstudy)[1,1]
	eDataLoc = which(grepl('methylation_hm27$', getGeneticProfiles(mycgds,mycancerstudy)[,1]))[1]
	eDataSel = getGeneticProfiles(mycgds,mycancerstudy)[eDataLoc,1]
	eData = getProfileData(mycgds,gn,eDataSel,mycaselist)
	mRNA[[i]] = eData
}
#collapse into a single frame
tRNA = do.call(rbind, mRNA)
xRNA = as.data.frame(t(tRNA))
xRNA$Gene = row.names(xRNA)
row.names(xRNA) = NULL
sRNA = subset(xRNA, rowSums(is.na(xRNA[,1:316]))<50)
sRNA$mMed = apply(sRNA[,1:316],1,function(x) mean(x,na.rm=TRUE))
#recombine to get the RLE values
mSet = merge(proS[,c(1,21)],sRNA[,317:318],by='Gene',sort=FALSE)
mSet = mSet[order(mSet$hgs),]
#choose the previous set of genes
gn = c('FOLR1','CRIP1','MSLN','SNCG','CRABP2','LEFTY1','GDF15','QPCT','GPC3','CTH')
#make some colors
cCol = col2rgb(ifelse(mSet$Gene %in% gn, brewer.pal(9,'Greens')[9],'gray40'))
cCex = ifelse(mSet$Gene %in% gn, 3,1)
#plot the data
pdf('ch_OvC_TMT10_wStd_Human_HGSMarkers_OvC-TCGA_methylation.pdf')
plot(mSet$mMed,
		col = rgb(cCol[1,],cCol[2,],cCol[3,],95,maxColorValue=255),
		ylab = 'methylation level',
		main = 'methylation expression patterns',
		cex = cCex,
		pch = 20
)
text(mSet$mMed,mSet$Gene,cex=0.25)
dev.off()











