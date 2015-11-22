# TODO: find a set of markers that agrees in the data
# 
# Author: cshughes
###############################################################################
#we have 3 data sets, ccle.q, pro, proCell
tis = pro[,c(2,4,7)]
cel = proCell[,c(2,4)]
gen = ccle.q[,c(2:3)]
#merge the data sets
tis.cel = merge(tis,cel,by='Gene')
tcg = merge(tis.cel,gen,by='Gene')
tcg$tisFC = rowMeans(tcg[,2:3], na.rm=TRUE)
tcg = tcg[,c(1,6,4,5)]
colnames(tcg) = c('Gene','tisFC','celFC','genFC')
#remove any values that are not present in all 3 data sets
tcg = subset(tcg, rowSums(is.na(tcg[,2:4]))<1)
#scale the values to be on an equal scale
xnorm = as.data.frame(normalize.quantiles(as.matrix(tcg[,2:4])))
tcg.q = cbind(tcg[,1],xnorm)
colnames(tcg.q) = c('Gene','tisFC','celFC','genFC')
#do a limma test on the values
reps = tcg.q[,c(2:4)]
fit <- lmFit(reps)
fit <- eBayes(fit)
p.value <- fit$p.value
p.adj = p.adjust(p.value, method="BH")
tcg.q$pAdj <- p.adj
tcg.q$meanFC<-rowMeans(tcg.q[,2:4],na.rm=TRUE)


##########################################################
##plot the markers based on the consolidated data
##########################################################
setwd(dir="/Users/cshughes/Documents/projects/OvC/")
genesCC = read.table("./Schwartz_ClearCellTable_geneList.txt", header=FALSE, sep='\t')
genesS = read.table("./Schwartz_SerousTable_geneList.txt", header=FALSE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/cell/Routput/")
#make the data a better list
geneCC = as.character(genesCC[,1])
geneS = as.character(genesS[,1])
geneTot = c(geneCC,geneS)
#make the plotting function 
pdf('ch_OvCall_consolidatedMarkers_vol.pdf')
lnCols<-brewer.pal(6,"Blues")
cols<-rev(brewer.pal(6,"RdBu"))
xCex = ifelse(tcg.q[,1] %in% geneTot, 2, 1)
proSD<-sd(tcg.q[,6], na.rm=TRUE)
xCol = col2rgb(ifelse(-log10(tcg.q[,5]) > -log10(0.05),ifelse(tcg.q[,6] > proSD, cols[6],ifelse(tcg.q[,6] < -proSD, cols[6], 'gray30')),'gray30'))
#plot with all points
plot(tcg.q[,6],
		-log10(tcg.q[,5]),
		col=rgb(xCol[1,],xCol[2,],xCol[3,],95,maxColorValue=255),
		pch=20,
		cex = xCex,
		ylab = '-log10(Adjusted p-value)',
		xlab = 'log2(serous/clear cell)',
		main = 'OvC markers',
		xlim = c(-5,5)
)
box(lwd=3)
abline(v = -proSD, col=lnCols[6],lwd=2,lty=2)
abline(v = proSD, col=lnCols[6],lwd=2,lty=2)
abline(h = -log10(0.05), col=lnCols[6],lwd=2,lty=2)
abline(h = -log10(0.01), col=lnCols[6],lwd=2,lty=2)
text(4.5,0.5,paste('n=',nrow(tcg.q),sep=""),cex=1.25)
text(4.5,1.6,paste('p<0.05'),cex=1.25)
text(4.5,2.3,paste('p<0.01'),cex=1.25)
dev.off()

##########################################################
##get the marker list out
##########################################################
ser = subset(tcg.q, -log10(pAdj) > -log10(0.05) & meanFC > 0)
cc = subset(tcg.q, -log10(pAdj) > -log10(0.05) & meanFC < 0)



