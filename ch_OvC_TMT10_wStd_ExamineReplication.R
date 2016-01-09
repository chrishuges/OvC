# TODO: Add comment
# 
# Author: cshughes
###############################################################################
###read in the data
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput") #change this to whatever directory you have stored the data in
a1e = readRDS('ch_OvC_wStd_normalizedPeptides_EColi_fa1.rds')
a2e = readRDS('ch_OvC_wStd_normalizedPeptides_EColi_fa2.rds')
a3e = readRDS('ch_OvC_wStd_normalizedPeptides_EColi_fa3.rds')
b1e = readRDS('ch_OvC_wStd_normalizedPeptides_EColi_fb1.rds')
b2e = readRDS('ch_OvC_wStd_normalizedPeptides_EColi_fb2.rds')
b3e = readRDS('ch_OvC_wStd_normalizedPeptides_EColi_fa3.rds')
###############################################################################
#look first at how the e-coli samples deviated
###############################################################################
a12e = merge(a1e[,c(1,2,4,5,8:17)],a2e[,c(1,2,4,5,8:17)],by=c('Accession','Gene','Descriptions','Sequence'))
a123e = merge(a12e,a3e[,c(1,2,4,5,8:17)],by=c('Accession','Gene','Descriptions','Sequence'))
b12e = merge(b1e[,c(1,2,4,5,8:17)],b2e[,c(1,2,4,5,8:17)],by=c('Accession','Gene','Descriptions','Sequence'))
b123e = merge(b12e,b3e[,c(1,2,4,5,8:17)],by=c('Accession','Gene','Descriptions','Sequence'))
###make an RLE plot
#merge all replicates
ab123e = merge(a123e,b123e,by=c('Accession','Gene','Descriptions','Sequence'))
#take the deviation from the median
ab123e$RLEmed = apply(ab123e[,5:64],1, function(x) median(x))
vnorm = apply(ab123e[,5:64],2, function(x) x - ab123e$RLEmed)
#make colors for the batches
hCols = brewer.pal(6,'Accent')
cols = c(rep(hCols[1],10),rep(hCols[2],10),rep(hCols[3],10),rep(hCols[4],10),rep(hCols[5],10),rep(hCols[6],10))
#make the plot
pdf('ch_OvC_TMT10_wStd_EColi_Peptides_RLEplot.pdf')
boxplot(vnorm,
		col=cols,
		pch=20,
		cex=0.75,
		las=2,
		xaxt='n',
		ylab = 'Relative Log Expression',
		main = 'RLE plot of EColi Batch Deviation',
		ylim = c(-3.5,3.5))
box(lwd=3)
dev.off()




###############################################################################
#plot the batches in a correlation heatmap with human peptides to look for batch effect
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
	xnorm = as.matrix(x[!spikeins,7:16])
	#make a noise model based on the spike ins
	spfit = justvsn(xnorm)
	nkid.h = as.data.frame(spfit)
	#normalize the values to the internal standard sample
	vnorm = apply(nkid.h[,1:9],2, function(x) x - nkid.h$aStd)
	#recombine the replicates
	vnormFull = as.data.frame(cbind(x[!spikeins,c(1:6,17)],vnorm))
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
#bind all data
a1 = rbind(a1n,a2n,a3n)
b1 = rbind(b1n,b2n,b3n)
#aggregate peptides
aggA = aggregate(cbind(a1,a2,a3,a4,a5,a6,a7,a8,a9)~Accession+Gene+Descriptions+Sequence,data=a1,median,na.action=na.pass,na.rm=TRUE)
aggB = aggregate(cbind(a1,a2,a3,a4,a5,a6,a7,a8,a9)~Accession+Gene+Descriptions+Sequence,data=b1,median,na.action=na.pass,na.rm=TRUE)
#merge the two batches
ab123h = merge(aggA,aggB,by=c('Accession','Gene','Descriptions','Sequence'))
#get the correlation maps
ovCorA = cor(ab123h[,c(5:22)], use='pairwise.complete.obs', method='pearson')
ovCorA[lower.tri(ovCorA)] <- NA
#make the plot
pdf('ch_OvC_TMT10_wStd_Human_Peptides_Batches_HeatMap.pdf')
#make the plot labels and boundaries
xLabels<- names(ovCorA)
mybreaks = seq(0,0.9,by=0.05)
#ColSideColors = c(rep(brewer.pal(6,'PRGn')[1],5),rep(brewer.pal(6,'PRGn')[6],5))
#make the correlation heatmap
heatmap.2(
		ovCorA,
		col= colorRampPalette(brewer.pal(9,"YlOrRd"))(length(mybreaks)-1),
		symkey=FALSE,
		Rowv=FALSE,
		Colv=FALSE,
		dendrogram="none",
		breaks=mybreaks,
		cellnote = round(ovCorA,2),
		#labRow = xLabels,
		#labCol = xLabels,
		notecol = 'black',
		notecex = 0.5,
		#colsep = 1:10,
		#rowsep = 1:10,
		sepwidth = c(0.03, 0.03),
		sepcolor = 'white',
		#ColSideColors=ColSideColors,
		## labels
		main='OvC Subtype Correlation',
		## color key
		key = TRUE,
		keysize = 1,
		density.info = "none",
		scale = "none",
		trace = "none",
		mar=c(8,8),
		cexRow=1.2,
		cexCol=1.2
)
dev.off()



###############################################################################
#plot the replicates in a correlation heatmap with human peptides to look for batch effect
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
	xnorm = as.matrix(x[!spikeins,7:16])
	#make a noise model based on the spike ins
	spfit = justvsn(xnorm)
	nkid.h = as.data.frame(spfit)
	#normalize the values to the internal standard sample
	vnorm = apply(nkid.h[,1:9],2, function(x) x - nkid.h$aStd)
	#recombine the replicates
	vnormFull = as.data.frame(cbind(x[!spikeins,c(1:6,17)],vnorm))
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
#merge the replicates
a12h = merge(a1n[,c(1,2,4,5,8:16)],a2n[,c(1,2,4,5,8:16)],by=c('Accession','Gene','Descriptions','Sequence'))
a123h = merge(a12h,a3n[,c(1,2,4,5,8:16)],by=c('Accession','Gene','Descriptions','Sequence'))
#get the correlation maps
ovCorA = cor(a123h[,c(5:31)], use='pairwise.complete.obs', method='pearson')
ovCorA[lower.tri(ovCorA)] <- NA
#make the plot
pdf('ch_OvC_TMT10_wStd_Human_Peptides_Replicates_HeatMap.pdf')
#make the plot labels and boundaries
xLabels<- names(ovCorA)
mybreaks = seq(0,0.9,by=0.05)
#ColSideColors = c(rep(brewer.pal(6,'PRGn')[1],5),rep(brewer.pal(6,'PRGn')[6],5))
#make the correlation heatmap
heatmap.2(
		ovCorA,
		col= colorRampPalette(brewer.pal(9,"YlOrRd"))(length(mybreaks)-1),
		symkey=FALSE,
		Rowv=FALSE,
		Colv=FALSE,
		dendrogram="none",
		breaks=mybreaks,
		cellnote = round(ovCorA,1),
		#labRow = xLabels,
		#labCol = xLabels,
		notecol = 'black',
		notecex = 0.5,
		#colsep = 1:10,
		#rowsep = 1:10,
		sepwidth = c(0.03, 0.03),
		sepcolor = 'white',
		#ColSideColors=ColSideColors,
		## labels
		main='OvC Subtype Correlation',
		## color key
		key = TRUE,
		keysize = 1,
		density.info = "none",
		scale = "none",
		trace = "none",
		mar=c(8,8),
		cexRow=1.2,
		cexCol=1.2
)
dev.off()



