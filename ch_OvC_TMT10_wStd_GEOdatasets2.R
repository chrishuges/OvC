# TODO: Add comment
# 
# Author: cshughes
###############################################################################
library(hgu133a.db)
###############################################################################
#get some data
###############################################################################
#Unpack the CEL files
setwd(dir="/Users/cshughes/Documents/projects/OvC/RNAexpression/")
untar("GSE6008_RAW.tar", exdir="gse6008data")
cels = list.files("gse6008data/", pattern = "CEL")
sapply(paste("gse6008data", cels, sep="/"), gunzip)
cels = list.files("gse6008data/", pattern = "CEL")
#go to the directory and create the phenotype annotation file
#$ ls data/*.CEL > data/phenodata.txt
#read the data
celfiles <- read.affy(covdesc="phenodata.txt", path="gse6008data")
celfiles.gcrma <- gcrma(celfiles)
#QC
#cols <- brewer.pal(8, "Set1")
#boxplot(celfiles, col=cols)
#boxplot(celfiles.gcrma, col=cols)
#hist(celfiles, col=cols)
#hist(celfiles.gcrma, col=cols)
#celfiles.qc <- fitPLM(celfiles)
#RLE(celfiles.qc, main="RLE")
#eset <- exprs(celfiles.gcrma)
#distance <- dist(t(eset),method="maximum")
#clusters <- hclust(distance)
#plot(clusters)
#filter the data
celfiles.filtered <- nsFilter(celfiles.gcrma, require.entrez=FALSE, remove.dupEntrez=FALSE)
celfiles.filtered$filter.log
#this is the data set you want to use
#output the data sets
saveRDS(celfiles.filtered,'ch_OvC_RNA2_processedProbes.rds')



###############################################################################
###prepare data to examine expression patterns between the RNA and Protein data
###############################################################################
###work the protein data
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput") #change this to whatever directory you have stored the data in
a1 = readRDS('ch_OvC_wStd_processedPeptides_a1.rds')
a2 = readRDS('ch_OvC_wStd_processedPeptides_a2.rds')
a3 = readRDS('ch_OvC_wStd_processedPeptides_a3.rds')
b1 = readRDS('ch_OvC_wStd_processedPeptides_b1.rds')
b2 = readRDS('ch_OvC_wStd_processedPeptides_b2.rds')
b3 = readRDS('ch_OvC_wStd_processedPeptides_b3.rds')
#bind the data together and preprocess
aIN = rbind(a1,a2,a3)
aSet = PECApre(aIN$Acc,aIN[,c(7:16)],peptides=aIN[,5],iStdtype='divide',iStdloc=aIN[,16],spike='ecoli',spikeloc=aIN[,3],normalize='median')
bIN = rbind(b1,b2,b3)
bSet = PECApre(bIN$Acc,bIN[,c(7:16)],peptides=bIN[,5],iStdtype='divide',iStdloc=bIN[,16],spike='ecoli',spikeloc=bIN[,3],normalize='median')
#grouping of samples
group1 = c('a1.x','a2.x','a3.x','a1.y','a2.y','a3.y')
group2 = c('a4.x','a5.x','a6.x','a4.y','a5.y','a6.y')
#merging batches
pSet = merge(aSet,bSet,by=c('Acc','Sequence'),all=TRUE)
#get stats
oSet = PECAstat(pSet$Acc,pSet[,c(group1,group2)],group1,group2,'modt','median')
colnames(oSet)[1] = 'PROexp'
oSet$Gene = sapply(strsplit(oSet$Acc,'_'),'[', 1)
#output the data sets
saveRDS(oSet,'ch_OvC_wStd_Proteins_HGSvCCC.rds')

###work the RNA data
rna = readRDS('ch_OvC_RNA2_processedProbes.rds')
#get out only the data you want to compare
rna.s = as.data.frame(exprs(rna$eset))[,c(59:99,1:8)]
#do DE analysis
samples <- c(rep('HGS',41),rep('CCC',8))
samples <- as.factor(samples)
design <- model.matrix(~0 + samples)
colnames(design) <- c("CCC", "HGS")
#apply limma
fit <- lmFit(rna.s, design)
contrast.matrix <- makeContrasts(HGS_CCC = HGS - CCC, levels=design)
ovc_fits <- contrasts.fit(fit, contrast.matrix)
ovc_ebFit <- eBayes(ovc_fits)
# coef=1 is CCC_EMC, coef=2 is CCC_HGS, coef=3 is HGS_EMC
topTable(ovc_ebFit, number=20, coef=1)
#annotate the data
#lfc here is a fold change cutoff
probes <- topTable(ovc_ebFit, coef=1, number=50000)
genes <- getSYMBOL(row.names(probes), "hgu133a")
results <- cbind(probes, genes)
rnaSet = aggregate(logFC~genes,data=results,median,na.action=na.pass,na.rm=TRUE)
colnames(rnaSet) = c('Gene','RNAexp')
#output the data sets
saveRDS(rnaSet,'ch_OvC_wStd_RNA2_HGSvCCC.rds')
###these are the processed RNA and Protein data
###############################################################################

###############################################################################
###examine expression patterns between the RNA and Protein data
###############################################################################
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput") #change this to whatever directory you have stored the data in
phvc = readRDS('ch_OvC_wStd_Proteins_HGSvCCC.rds')
phve = readRDS('ch_OvC_wStd_Proteins_HGSvEMC.rds')
pcve = readRDS('ch_OvC_wStd_Proteins_CCCvEMC.rds')
rhvc = readRDS('ch_OvC_wStd_RNA2_HGSvCCC.rds')
#rhve = readRDS('ch_OvC_wStd_RNA2_HGSvEMC.rds')
#rcve = readRDS('ch_OvC_wStd_RNA2_CCCvEMC.rds')

#merge the protein and RNA data
hvc = merge(phvc,rhvc,by='Gene')
#hve = merge(phve,rhve,by='Gene')
#cve = merge(pcve,rcve,by='Gene')
#make the plot for all three sets
#bring in marker data
setwd(dir="/Users/cshughes/Documents/projects/OvC/markers/")
genesCC = read.table("./Schwartz_ClearCellTable_geneList.txt", header=FALSE, sep='\t')
genesS = read.table("./Schwartz_SerousTable_geneList.txt", header=FALSE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput/")
#combine into a single set
geneCC = as.character(genesCC[,1])
geneS = as.character(genesS[,1])
geneTot = c(geneCC,geneS)
#only keep markers that are present in the data
markHGSCCC = geneTot[which(geneTot %in% phvc$Gene)]
markHGS = geneS[which(geneS %in% phvc$Gene)]
markCCC = geneCC[which(geneCC %in% phvc$Gene)]
#make the plot
pdf('ch_OvC_wStd_Proteins-RNA2_HGSvCCC_CorrScatter.pdf')
cols<-rev(brewer.pal(6,"RdBu"))
xCol = col2rgb(ifelse(hvc$Gene %in% markCCC, cols[1], ifelse(hvc$Gene %in% markHGS, cols[6],'gray60')))
xCex = ifelse(hvc$Gene %in% markHGSCCC, 2, 1)
plot(hvc$PROexp,hvc$RNAexp,pch=20,col=rgb(xCol[1,],xCol[2,],xCol[3,],95,maxColorValue=255),cex=xCex,xlim = c(-4,4),ylim = c(-9,9))
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
plot(hvc.s$PROexp,hvc.s$RNAexp,pch=20,col=rgb(xCol[1,],xCol[2,],xCol[3,],95,maxColorValue=255),cex=xCex,xlim = c(-4,4),ylim = c(-9,9))
text(hvc.s$PROexp,hvc.s$RNAexp,hvc.s$Gene,cex=0.5)
text(1,-3,paste('r = ',cor(hvc.s$PROexp,hvc.s$RNAexp,use='pairwise.complete.obs'),sep=''))
dev.off()



