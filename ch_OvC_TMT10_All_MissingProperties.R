# TODO: Add comment
# 
# Author: cshughes
###############################################################################
##HPA RNAseq Data
setwd(dir="/Users/cshughes/Documents/projects/OvC/RNAexpression/")
ens<-read.table("./HumanProteomeMap_Science_TableS1.txt", header=TRUE, sep='\t')
hpa<-read.table("./HumanProteinAtlas_FPKM_allsamples.txt", header=TRUE, sep='\t')
#change the gene colname
colnames(ens)[2] = 'Gene'
colnames(hpa)[1] = 'ensg_id'
#merge the two human atlas sets
hpa.m = merge(ens,hpa, by='ensg_id')
#remove isoforms
isoforms = grepl('\\.', hpa.m$Gene)
isoforms2 = grepl('\\-', hpa.m$Gene)
hpa.mi = hpa.m[!isoforms,]
hpa.mi = hpa.mi[!isoforms2,]


##protein data from SP3
setwd(dir="/Users/cshughes/Documents/projects/OvC/FFPE/Routput") #change this to whatever directory you have stored the data in
pro = readRDS('ch_OvC_FFPE_proteinSet.rds')
colnames(pro)[5] = 'PROexp'

setwd(dir="/Users/cshughes/Documents/projects/OvC/Frozen/Routput") #change this to whatever directory you have stored the data in
pro = readRDS('ch_OvC_Frozen_proteinSet.rds')
colnames(pro)[5] = 'PROexp'

setwd(dir="/Users/cshughes/Documents/projects/OvC/CellLine/Routput") #change this to whatever directory you have stored the data in
pro = readRDS('ch_OvC_CellLine_proteinSet.rds')
colnames(pro)[5] = 'PROexp'

#combine with the protein data
prAll = merge(pro[,c(1:5,7)],hpa.mi[,c(1:2,130)],by='Gene',all=TRUE,sort=FALSE)


#####first look at proteins that have an expression value where FPKM is below 1
prSub = subset(prAll, log2(prAll$FPKM.ovary_6a.V233.)>=log2(10) & is.na(prAll$PROexp))
write.table(prSub,'ch_OvC_FFPE_Proteins_highRNAexp.txt',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)








###############################################################################
#Make a volcano plot of the variance between subtypes with overlaid markers
###############################################################################
#bring in the protein data
proh = prSub
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
proSD<-sd(proh$PROexp, na.rm=TRUE)
#sort out colors
lnCols<-brewer.pal(6,"Blues")
cols<-rev(brewer.pal(6,"RdBu"))

###make an initial plot with all points
pdf('ch_OvC_TMT10_FFPE_Human_Proteins_lowRNAexp_Volcano.pdf')
xCol = col2rgb(ifelse(proh$Gene %in% markCCC, cols[1], ifelse(proh$Gene %in% markHGS, cols[6],'gray80')))
xCex = ifelse(proh$Gene %in% markHGSCCC, 2, 1)
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
		xlab = 'log2(HGSC/CCC)',
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
text(proh.s$PROexp, -log10(proh.s$score), proh.s$Gene)
dev.off()



