# TODO: Add comment
# 
# Author: cshughes
###############################################################################
#read in the data objects
setwd(dir="/Users/cshughes/Documents/projects/OvC/FFPE/Routput") #change this to whatever directory you have stored the data in
a1h = readRDS('ch_OvC_FFPE_normalizedPeptides_fa1.rds')
a2h = readRDS('ch_OvC_FFPE_normalizedPeptides_fa2.rds')
b1h = readRDS('ch_OvC_FFPE_normalizedPeptides_fb1.rds')
b2h = readRDS('ch_OvC_FFPE_normalizedPeptides_fb2.rds')

###############################################################################
#compare HGS with CCC using limma ####USE THE PECA ANALYSIS INSTEAD###
###############################################################################
#bind all data
allh = rbind(a1h,a2h,b1h,b2h)
aggh = aggregate(cbind(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)~Accession+Gene+Descriptions+Sequence,data=allh,median,na.action=na.pass,na.rm=TRUE)
#need to reshuffle the patients based on HGS or CCC
aggh.r = aggh[,c(1:4,9,14,10,5,13,8,6,7,12,11)]
#calculate p-values
eset = as.matrix(aggh.r[,c(5:14)])
design <- cbind(SR=c(1,1,1,1,1,0,0,0,0,0),CC=c(0,0,0,0,0,1,1,1,1,1))
fit <- lmFit(eset,design)
cont.matrix <- makeContrasts(SRvsCC=SR-CC, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
aggh.r$logFC = fit2$coef
aggh.r$pVal = fit2$p.value
#add a peptide counter
aggh.r$pepNum = 1
#combine into proteins
proh1 = aggregate(pepNum~Accession+Gene+Descriptions,data=aggh.r,sum,na.action=na.pass,na.rm=TRUE)
proh2 = aggregate(cbind(logFC,pVal)~Accession+Gene+Descriptions,data=aggh.r,median,na.action=na.pass,na.rm=TRUE)
proh = merge(proh1,proh2,by=c('Accession','Gene','Descriptions'))
colnames(proh) = c('Accession','Gene','Descriptions','pepNum','logFC','pVal')
#write out the data into a format usable in excel
write.table(proh,'ch_OvC_FFPE_HGSvCCC_proteinSet.txt',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)

###############################################################################
#compare HGS with CCC using PECA
###############################################################################
library('PECA')
#need to use the non-VSN transformed data here
#read in the data objects
setwd(dir="/Users/cshughes/Documents/projects/OvC/FFPE/Routput") #change this to whatever directory you have stored the data in
a1h = readRDS('ch_OvC_FFPE_processedPeptides_fa1.rds')
a2h = readRDS('ch_OvC_FFPE_processedPeptides_fa2.rds')
b1h = readRDS('ch_OvC_FFPE_processedPeptides_fb1.rds')
b2h = readRDS('ch_OvC_FFPE_processedPeptides_fb2.rds')
#bind all data
allh = rbind(a1h,a2h,b1h,b2h)
aggh = aggregate(cbind(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)~Accession+Gene+Descriptions+Sequence+Organism,data=allh,median,na.action=na.pass,na.rm=TRUE)
#remove any e-coli peptides
spikeins = grepl('EColi', aggh$Organism)
inSet = aggh[!spikeins,c(1:4,6:15)]
#need to reshuffle the patients based on HGS or CCC
pec = inSet[,c(1:4,9,14,10,5,13,8,6,7,12,11)]
#calculate p-values
#assign the groups to compare for stats
group1<-colnames(pec)[c(5:9)]
group2<-colnames(pec)[c(10:14)]
#the data will output as group1/group2 in the fold change
pec.q = PECA_df(pec,'Gene',group1,group2,normalize='median',test='modt',type='median')
#reshape the data frame to reincorporate the meta data
pec.q$Gene = rownames(pec.q)
rownames(pec.q) = NULL
#add a peptide number counter
pec$pepNum = 1
#aggregate into proteins
anno = aggregate(pepNum~Accession+Gene+Descriptions,data=pec,sum,na.action=na.pass,na.rm=TRUE)
prop = merge(pec.q,anno,by='Gene')
#reshape and reannotate dataframe
pecP = prop[,c(1,7,8,9,2:6)]
colnames(pecP) = c('Gene','Accession','Descriptions','pepNum','logFC','t','score','p-value','pVal')
#write out the data into a format usable in excel
saveRDS(pecP,'ch_OvC_FFPE_proteinSet.rds')
write.table(pecP,'ch_OvC_FFPE_HGSvCCC_PECA_proteinSet.txt',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)
#I think the PECA set looks more reliable


###############################################################################
#Make a volcano plot of the variance between subtypes with overlaid markers
###############################################################################
#bring in the protein data
proh = readRDS('ch_OvC_FFPE_proteinSet.rds')
#bring in marker data
setwd(dir="/Users/cshughes/Documents/projects/OvC/markers/")
genesCC = read.table("./Schwartz_ClearCellTable_geneList.txt", header=FALSE, sep='\t')
genesS = read.table("./Schwartz_SerousTable_geneList.txt", header=FALSE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/FFPE/Routput/")
#combine into a single set
geneCC = as.character(genesCC[,1])
geneS = as.character(genesS[,1])
geneTot = c(geneCC,geneS)
#only keep markers that are present in the data
markHGSCCC = geneTot[which(geneTot %in% proh$Gene)]
markHGS = geneS[which(geneS %in% proh$Gene)]
markCCC = geneCC[which(geneCC %in% proh$Gene)]

#make a plot of the data
proSD<-sd(proh$logFC, na.rm=TRUE)
#sort out colors
lnCols<-brewer.pal(6,"Blues")
cols<-rev(brewer.pal(6,"RdBu"))

###make an initial plot with all points
pdf('ch_OvC_TMT10_FFPE_Human_Proteins_HGSvCCC_Volcano.pdf')
xCol = col2rgb(ifelse(proh$Gene %in% markCCC, cols[1], ifelse(proh$Gene %in% markHGS, cols[6],'gray80')))
xCex = ifelse(proh$Gene %in% markHGSCCC, 2, 1)
plot(proh$logFC,
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
plot(proh.s$logFC,
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



###############################################################################
#Make a plot of segregation of subtypes based on differential expression...need to use Limma fold changes for this
###############################################################################
#want to use peptide data for this, the RLE data specifically
#bind replicate data
a12h = rbind(a1h,a2h)
b12h = rbind(b1h,b2h)
#design function for processing the FC values
bindFC = function(x,...){
	pep = aggregate(cbind(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)~Accession+Gene+Descriptions+Sequence,data=x,median,na.action=na.pass,na.rm=TRUE)
	#need to reshuffle the patients based on HGS or CCC
	pep.r = pep[,c(1:4,9,14,10,5,13,8,6,7,12,11)]
	#calculate p-values
	eset = as.matrix(pep.r[,c(5:14)])
	design <- cbind(SR=c(1,1,1,1,1,0,0,0,0,0),CC=c(0,0,0,0,0,1,1,1,1,1))
	fit <- lmFit(eset,design)
	cont.matrix <- makeContrasts(SRvsCC=SR-CC, levels=design)
	fit2 <- contrasts.fit(fit, cont.matrix)
	fit2 <- eBayes(fit2)
	pep.r$logFC = fit2$coef
	pep.r$pVal = fit2$p.value
	#pep.a = aggregate(cbind(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,logFC,pVal)~Accession+Gene+Descriptions,data=pep.r,na.action=na.pass,FUN=median,na.rm=TRUE)
	#colnames(pep.a)[14:15] = c('logFC','pVal')
	#pep.a = subset(pep.a, rowSums(is.na(pep.a[,4:13]))<2)
	#subset anything with a small FC
	#vnorm = subset(pep.r, logFC< -1 | logFC> 1)
	pep.r = pep.r[order(-abs(pep.r$logFC)),]
	vnorm = pep.r[1:1000,]
	#output the data
	return(vnorm)
}
#apply function
ovA = bindFC(a12h)
ovB = bindFC(b12h)
#get the correlation maps
ovCor = cor(ovA[,c(5:14)], use='pairwise.complete.obs', method='pearson')
ovCorB = cor(ovB[,c(5:14)], use='pairwise.complete.obs', method='pearson')
ovCor[lower.tri(ovCor)] <- ovCorB[lower.tri(ovCorB)]
#make the plot
pdf('ch_OvC_TMT10_FFPE_Human_Peptides_HGSvCCC_HeatMap_top500.pdf')
#pdf('ch_test.pdf')
#make the plot labels and boundaries
xLabels<- names(ovCor)
mybreaks = seq(0,0.9,by=0.05)
ColSideColors = c(rep(brewer.pal(6,'PRGn')[1],5),rep(brewer.pal(6,'PRGn')[6],5))
#make the correlation heatmap
heatmap.2(
		ovCor,
		col= colorRampPalette(brewer.pal(9,"YlOrRd"))(length(mybreaks)-1),
		symkey=FALSE,
		Rowv=FALSE,
		Colv=FALSE,
		dendrogram="none",
		breaks=mybreaks,
		cellnote = round(ovCor,2),
		labRow = xLabels,
		labCol = xLabels,
		notecol = 'black',
		notecex = 1.2,
		colsep = 1:10,
		rowsep = 1:10,
		sepwidth = c(0.03, 0.03),
		sepcolor = 'white',
		ColSideColors=ColSideColors,
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


##########################################################
##number of peptides per protein for those exceeding 0.05 p-value
##########################################################
#use proh for this
#transform the p-values
pecP$pValT = -log10(pecP$pVal)
#subset non-significant proteins
pSub = subset(pecP, pecP$pValT >= -log10(0.05))

#make the plot
pdf('ch_OvC_TMT10_FFPE_Human_Proteins_HGSvCCC_pValDistribution.pdf')
x = hist(log2(pSub$pepNum),breaks=20,plot=FALSE)
mp<-barplot(x$density,plot=FALSE)
barplot(x$density,
		space=rep(0.1,17),
		main='distribution of significant hits', 
		xlab='log2(pepNum)',
		xaxt="n"
)
axis(1,x$breaks,seq(0,8.5,0.5),las=2,lwd=2)
dev.off()







