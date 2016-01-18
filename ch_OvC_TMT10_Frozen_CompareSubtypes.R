# TODO: Add comment
# 
# Author: cshughes
###############################################################################
###############################################################################
#compare HGS with CCC using PECA
###############################################################################
library('PECA')
#need to use the non-VSN transformed data here
#read in the data objects
setwd(dir="/Users/cshughes/Documents/projects/OvC/Frozen/Routput") #change this to whatever directory you have stored the data in
a1h = readRDS('ch_OvC_Frozen_processedPeptides_fa1.rds')
a2h = readRDS('ch_OvC_Frozen_processedPeptides_fa2.rds')
b1h = readRDS('ch_OvC_Frozen_processedPeptides_fb1.rds')
b2h = readRDS('ch_OvC_Frozen_processedPeptides_fb2.rds')
#bind all data
allh = rbind(a1h,a2h,b1h,b2h)
aggh = aggregate(cbind(a1,a2,a3,a4,a5,a6,a7,a8)~Accession+Gene+Descriptions+Sequence+Organism,data=allh,median,na.action=na.pass,na.rm=TRUE)
#remove any e-coli peptides
spikeins = grepl('EColi', aggh$Organism)
inSet = aggh[!spikeins,c(1:4,6:13)]
#need to reshuffle the patients based on HGS or CCC
pec = inSet[,c(1:4,9:12,5:8)]
#calculate p-values
#assign the groups to compare for stats
group1<-colnames(pec)[c(5:8)]
group2<-colnames(pec)[c(9:12)]
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
#write out the data into a format usable in excel.
saveRDS(pecP,'ch_OvC_Frozen_proteinSet.rds')
write.table(pecP,'ch_OvC_Frozen_HGSvCCC_PECA_proteinSet.txt',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)
#I think the PECA set looks more reliable


###############################################################################
#Make a volcano plot of the variance between subtypes with overlaid markers
###############################################################################
#bring in the protein data
proh = pecP
#bring in marker data
setwd(dir="/Users/cshughes/Documents/projects/OvC/markers/")
genesCC = read.table("./Schwartz_ClearCellTable_geneList.txt", header=FALSE, sep='\t')
genesS = read.table("./Schwartz_SerousTable_geneList.txt", header=FALSE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/Frozen/Routput/")
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
pdf('ch_OvC_TMT10_Frozen_Human_Proteins_HGSvCCC_Volcano.pdf')
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
text(proh.s$logFC, -log10(proh.s$score), proh.s$Gene)
dev.off()

