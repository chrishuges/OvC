# TODO: Add comment
# 
# Author: cshughes
###############################################################################
#prepare the RNA data
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput") #change this to whatever directory you have stored the data in
rna = readRDS('ch_OvC_wStd_RNA_HGSvCCC.rds')

#choose protein set
setwd(dir="/Users/cshughes/Documents/projects/OvC/FFPE/Routput") #change this to whatever directory you have stored the data in
pro = readRDS('ch_OvC_FFPE_proteinSet.rds')
colnames(pro)[5] = 'PROexp'

setwd(dir="/Users/cshughes/Documents/projects/OvC/Frozen/Routput") #change this to whatever directory you have stored the data in
pro = readRDS('ch_OvC_Frozen_proteinSet.rds')
colnames(pro)[5] = 'PROexp'

setwd(dir="/Users/cshughes/Documents/projects/OvC/CellLine/Routput") #change this to whatever directory you have stored the data in
pro = readRDS('ch_OvC_CellLine_proteinSet.rds')
colnames(pro)[5] = 'PROexp'



#merge the protein and RNA data
hvc = merge(pro,rna,by='Gene')
#make the plot for all three sets
#bring in marker data
setwd(dir="/Users/cshughes/Documents/projects/OvC/markers/")
genesCC = read.table("./Schwartz_ClearCellTable_geneList.txt", header=FALSE, sep='\t')
genesS = read.table("./Schwartz_SerousTable_geneList.txt", header=FALSE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/CellLine/Routput/")####CHANGE ME
#combine into a single set
geneCC = as.character(genesCC[,1])
geneS = as.character(genesS[,1])
geneTot = c(geneCC,geneS)
#only keep markers that are present in the data
markHGSCCC = geneTot[which(geneTot %in% phvc$Gene)]
markHGS = geneS[which(geneS %in% phvc$Gene)]
markCCC = geneCC[which(geneCC %in% phvc$Gene)]
#make the plot
pdf('ch_OvC_CellLine_Proteins-RNA_HGSvCCC_CorrScatter.pdf')
cols<-rev(brewer.pal(6,"RdBu"))
xCol = col2rgb(ifelse(hvc$Gene %in% markCCC, cols[1], ifelse(hvc$Gene %in% markHGS, cols[6],'gray60')))
xCex = ifelse(hvc$Gene %in% markHGSCCC, 2, 1)
plot(hvc$PROexp,hvc$RNAexp,pch=20,col=rgb(xCol[1,],xCol[2,],xCol[3,],95,maxColorValue=255),cex=xCex,xlim = c(-4,4),ylim = c(-7,7))
text(1,-3,paste('r = ',cor(hvc$PROexp,hvc$RNAexp,use='pairwise.complete.obs'),sep=''))
#text(hvc$PROexp,hvc$RNAexp,hvc$Gene,cex=0.25)
abline(lm(hvc$RNAexp~hvc$PROexp),lwd=2,lty=2)
abline(h=0,lwd=2,lty=2)
abline(v=0,lwd=2,lty=2)
box(lwd=3)
xCol = ifelse(hvc$Gene %in% markCCC, cols[1], ifelse(hvc$Gene %in% markHGS, cols[6],'gray60'))
gnRM = !grepl('gray', xCol)
hvc.s = hvc[gnRM,]
xCol = col2rgb(ifelse(hvc.s$Gene %in% markCCC, cols[1], ifelse(hvc.s$Gene %in% markHGS, cols[6],'gray60')))
xCex = ifelse(hvc.s$Gene %in% markHGSCCC, 2, 1)
plot(hvc.s$PROexp,hvc.s$RNAexp,pch=20,col=rgb(xCol[1,],xCol[2,],xCol[3,],95,maxColorValue=255),cex=xCex,xlim = c(-4,4),ylim = c(-7,7))
text(hvc.s$PROexp,hvc.s$RNAexp,hvc.s$Gene,cex=0.5)
text(1,-3,paste('r = ',cor(hvc.s$PROexp,hvc.s$RNAexp,use='pairwise.complete.obs'),sep=''))
dev.off()





