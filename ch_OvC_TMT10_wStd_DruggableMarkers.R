# TODO: overlay druggable markers onto HGSC vs. CCC comparison
# 
# Author: cshughes
###############################################################################
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput") #change this to whatever directory you have stored the data in

#making the plot
proh = readRDS('ch_OvC_wStd_Proteins_HGSvCCC.rds')
#get the list of candidate genes
setwd(dir="/Users/cshughes/Documents/projects/OvC/markers/")
pog = read.table("./POGMarkers.txt", header=TRUE, sep='\t')
vIHC = read.table("./vpIHC_GeneList.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput/")
#
#combine into a single set
colnames(pog) = c('ensg_id','Gene')
pogM = as.character(pog[,2])
#pogM = as.character(vIHC$Gene)
#only keep markers that are present in the data
markPOG = pogM[which(pogM %in% proh$Gene)]

#make a plot of the data
proSD<-sd(proh$PROexp, na.rm=TRUE)
#sort out colors
lnCols<-brewer.pal(6,"Blues")
cols<-rev(brewer.pal(6,"RdBu"))
pogCols<-brewer.pal(6, "Greens")

###make an initial plot with all points
pdf('ch_OvC_TMT10_wStd_Human_Proteins_HGSvCCC-vIHC_Volcano.pdf')
xCol = col2rgb(ifelse(proh$Gene %in% markPOG, pogCols[5], 'gray80'))
xCex = ifelse(proh$Gene %in% markPOG, 2, 1)
plot(proh$PROexp,
		-log10(proh$score),
		col=rgb(xCol[1,],xCol[2,],xCol[3,],95,maxColorValue=255),
		pch=20,
		cex = xCex,
		ylab = '-log10(Adjusted p-value)',
		xlab = 'log2(HGSC/CCC)',
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
xCol = ifelse(proh$Gene %in% markPOG, pogCols[5], 'gray80')
gnRM = !grepl('gray', xCol)
proh.s = proh[gnRM,]
xCol = col2rgb(ifelse(proh.s$Gene %in% markPOG, pogCols[5], 'gray80'))
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
#text(proh.s$logFC, -log10(proh.s$score), proh.s$Gene)
dev.off()







