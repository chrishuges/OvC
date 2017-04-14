# TODO: Add comment
# 
# Author: cshughes
###############################################################################
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput") #change this to whatever directory you have stored the data in
phvc = readRDS('ch_OvC_wStd_Proteins_HGSvCCC.rds')
phve = readRDS('ch_OvC_wStd_Proteins_HGSvEMC.rds')
pcve = readRDS('ch_OvC_wStd_Proteins_CCCvEMC.rds')

###############################################################################
#Make a volcano plot of the variance between subtypes with overlaid markers
###############################################################################
#bring in the protein data
proh = phve
#make a plot of the data
proSD<-sd(proh$PROexp, na.rm=TRUE)
#sort out colors
lnCols<-brewer.pal(6,"Blues")
cols<-rev(brewer.pal(6,"RdBu"))

###make an initial plot with all points
pdf('ch_OvC_TMT10_wStd_Human_Proteins_HGSvEMC_Volcano.pdf')
xCol = col2rgb(ifelse(proh$PROexp > proSD & -log10(proh$score) > -log10(0.05), cols[6], ifelse(proh$PROexp < -proSD & -log10(proh$score) > -log10(0.05), cols[6],'gray60')))
plot(proh$PROexp,
		-log10(proh$score),
		col=rgb(xCol[1,],xCol[2,],xCol[3,],95,maxColorValue=255),
		pch=20,
		cex = 2,
		ylab = '-log10(PECA score)',
		xlab = 'log2(HGSvEMC)',
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
text(proh$PROexp,-log10(proh$score),proh$Gene,cex=0.25)
dev.off()





#write out the data
write.table(phvc,'ch_oct2016_HGSCvsCCC_wStd_proteinSet.txt',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)


