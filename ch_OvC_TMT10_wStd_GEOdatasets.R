# TODO: Add comment
# 
# Author: cshughes
###############################################################################
library(GEOquery)
library(gcrma)
library(genefilter)
library(hgu133plus2probe)
library(hgu133plus2.db)
library(simpleaffy)
library(affyPLM)
library(annotate)
###############################################################################
#get some data
###############################################################################
#Unpack the CEL files
setwd(dir="/Users/cshughes/Documents/projects/OvC/RNAexpression/")
untar("GSE65986_RAW.tar", exdir="gse65986data")
cels = list.files("gse65986data/", pattern = "CEL")
sapply(paste("gse65986data", cels, sep="/"), gunzip)
cels = list.files("gse65986data/", pattern = "CEL")
#go to the directory and create the phenotype annotation file
#$ ls data/*.CEL > data/phenodata.txt
#read the data
celfiles <- read.affy(covdesc="phenodata.txt", path="gse65986data")
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
saveRDS(celfiles.filtered,'ch_OvC_RNA_processedProbes.rds')
#the below analysis compares all subtypes with each other
#do DE analysis
samples <- celfiles.gcrma$Target
samples <- as.factor(samples)
design <- model.matrix(~0 + samples)
colnames(design) <- c("CCC", "EMC", "HGS")
#apply limma
fit <- lmFit(exprs(celfiles.filtered$eset), design)
contrast.matrix <- makeContrasts(CCC_EMC = CCC - EMC, HGS_CCC = HGS - CCC, HGS_EMC = HGS - EMC, levels=design)
ovc_fits <- contrasts.fit(fit, contrast.matrix)
ovc_ebFit <- eBayes(ovc_fits)
# coef=1 is CCC_EMC, coef=2 is CCC_HGS, coef=3 is HGS_EMC
topTable(ovc_ebFit, number=10, coef=1)
#annotate the data
#lfc here is a fold change cutoff
probeset.list <- topTable(ovc_ebFit, coef=1, number=10, lfc=0)
gene.symbols <- getSYMBOL(exprs(celfiles.filtered$eset[row.names(probeset.list),]), "hgu133plus2")
results <- cbind(probeset.list, gene.symbols)
write.table(results, "results.txt", sep="\t", quote=FALSE)

#output the data sets
#saveRDS(results,'ch_OvC_Frozen_processedPeptides_fa1.rds')



###############################################################################
###clustering heat maps by subtype
###############################################################################
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput") #change this to whatever directory you have stored the data in
proG = readRDS('ch_OvC_wStd_Top500Genes_Variance.rds')
proS = subset(proG, hgs1 > 0.5 & ccc1 < -0.5)
proList = proS$Gene
#get the rna expression data
rnaG = exprs(celfiles.filtered$eset)
gene.symbols <- getSYMBOL(row.names(rnaG), "hgu133plus2")
#locate genes from the protein data
geneLoc = which(gene.symbols %in% proList)
#extract them from the RNA data
rnaE = as.data.frame(rnaG[geneLoc,])
rnaE$Gene = getSYMBOL(row.names(rnaE), "hgu133plus2")
rnaA = aggregate(rnaE[,1:55],by=list(rnaE$Gene),median,na.rm=TRUE)
#take the median of each row
rnaA$RLEmed = apply(rnaA[,2:56], 1, function(x) median(x,na.rm=TRUE))
#get the RLE values
vnorm = apply(rnaA[,2:56], 2, function(x) x - rnaA$RLEmed)

#built the heatmap
pdf('ch_OvC_TMT10_wStd_Human_RNA_Subtypes_Top500_HeatMap.pdf')
#make the plot labels and boundaries
xLabels<- c(rep('CCC',25),rep('EMC',14),rep('HGS',16))
mybreaks = seq(-2,2,by=0.05)
ColSideColors = c(rep(brewer.pal(6,'Accent')[1],25),rep(brewer.pal(6,'Accent')[2],14),rep(brewer.pal(6,'Accent')[3],16))
#make the correlation heatmap
heatmap.2(
		as.matrix(vnorm),
		col= colorRampPalette(brewer.pal(6,"RdBu"))(length(mybreaks)-1),
		symkey=TRUE,
		Rowv=TRUE,
		Colv=TRUE,
		dendrogram="both",
		breaks=mybreaks,
		#cellnote = round(ovCor,2),
		labRow = '',
		labCol = xLabels,
		las=2,
		ColSideColors=ColSideColors,
		## labels
		main='OvC Subtype Clustering',
		## color key
		key = TRUE,
		keysize = 1,
		density.info = "none",
		scale = "none",
		trace = "none",
		mar=c(8,8),
		cexRow=1.5,
		cexCol=1.5
)
dev.off()


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
rna = readRDS('ch_OvC_RNA_processedProbes.rds')
#get out only the data you want to compare
rna.s = as.data.frame(exprs(rna$eset))[,c(40:55,1:25)]
#do DE analysis
samples <- c(rep('HGS',16),rep('CCC',25))
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
genes <- getSYMBOL(row.names(probes), "hgu133plus2")
results <- cbind(probes, genes)
rnaSet = aggregate(cbind(logFC,t)~genes,data=results,median,na.action=na.pass,na.rm=TRUE)
colnames(rnaSet) = c('Gene','RNAexp','t')
#output the data sets
saveRDS(rnaSet,'ch_OvC_wStd_RNA_HGSvCCC.rds')
###these are the processed RNA and Protein data
###############################################################################



###############################################################################
###examine expression patterns between the RNA and Protein data
###############################################################################
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput") #change this to whatever directory you have stored the data in
phvc = readRDS('ch_OvC_wStd_Proteins_HGSvCCC.rds')
phve = readRDS('ch_OvC_wStd_Proteins_HGSvEMC.rds')
pcve = readRDS('ch_OvC_wStd_Proteins_CCCvEMC.rds')
rhvc = readRDS('ch_OvC_wStd_RNA_HGSvCCC.rds')
rhve = readRDS('ch_OvC_wStd_RNA_HGSvEMC.rds')
rcve = readRDS('ch_OvC_wStd_RNA_CCCvEMC.rds')

#merge the protein and RNA data
hvc = merge(phvc,rhvc,by='Gene')
hve = merge(phve,rhve,by='Gene')
cve = merge(pcve,rcve,by='Gene')
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
pdf('ch_OvC_wStd_Proteins-RNA_HGSvCCC_CorrScatter.pdf')
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

###plot the other two sets, no markers
pdf('ch_OvC_wStd_Proteins-RNA_HGSvEMC_CorrScatter.pdf')
cols<-col2rgb(rev(brewer.pal(6,"RdBu")))
plot(hve$PROexp,hve$RNAexp,pch=20,col=rgb(cols[1,1],cols[2,1],cols[3,1],95,maxColorValue=255),cex=1.5)
text(1,-3,paste('r = ',cor(hve$PROexp,hve$RNAexp,use='pairwise.complete.obs'),sep=''))
#text(hve$PROexp,hve$RNAexp,hve$Gene,cex=0.25)
abline(lm(hve$RNAexp~hve$PROexp),lwd=2,lty=2)
abline(h=0,lwd=2,lty=2)
abline(v=0,lwd=2,lty=2)
box(lwd=3)
dev.off()

pdf('ch_OvC_wStd_Proteins-RNA_CCCvEMC_CorrScatter.pdf')
cols<-col2rgb(rev(brewer.pal(6,"RdBu")))
plot(cve$PROexp,cve$RNAexp,pch=20,col=rgb(cols[1,1],cols[2,1],cols[3,1],95,maxColorValue=255),cex=1.5)
text(1,-3,paste('r = ',cor(cve$PROexp,cve$RNAexp,use='pairwise.complete.obs'),sep=''))
#text(cve$PROexp,cve$RNAexp,cve$Gene,cex=0.25)
abline(lm(cve$RNAexp~cve$PROexp),lwd=2,lty=2)
box(lwd=3)
abline(h=0,lwd=2,lty=2)
abline(v=0,lwd=2,lty=2)
dev.off()




cor(x$PROexp,x$RNAexp,use='pairwise.complete.obs')


write.table(aSet,'test.txt',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)
l = merge(aNorm,bNorm,by='Acc',all=FALSE)
cor(l$logFC.x,l$logFC.y,use='pairwise.complete.obs')
l$logFC = rowMeans(l[,c(2,7)],na.rm=TRUE)








































###############################################################################
#random code for optimizing the above section
###############################################################################
gset = getGEO('GSE65986', GSEMatrix=TRUE)
#unload the list into an eSet
gunl = gset[[1]]
#access the expression data
gex = exprs(gunl)
#access the annotation
ganno = fData(gunl)
setwd("/Users/cshughes/Documents/projects/OvC/RNAexpression/gse65986data/")
raw.data=ReadAffy(verbose=TRUE, filenames=cels) #From bioconductor
#raw.data=ReadAffy(verbose=TRUE, filenames=cels, cdfname="hgu133plus2cdf") #From bioconductor
#perform GCRMA normalization
data.rma.norm=gcrma(raw.data)
celfiles.filtered <- nsFilter(data.rma.norm, require.entrez=FALSE, remove.dupEntrez=FALSE)
# What got removed and why?
celfiles.filtered$filter.log
samples <- data.rma.norm$Target
#Get the expression estimates for each array
rma=exprs(data.rma.norm)
#Format values to 5 decimal places
rma=as.data.frame(format(rma, digits=5))
#Map probe sets to gene symbols or other annotations
#Extract probe ids, entrez symbols, and entrez ids
probes=row.names(rma)
rma$PROBEID = row.names(rma)
#test call
#select(hgu133plus2.db, "207356_at", c("SYMBOL","ENTREZID"))
GeneAnno = select(hgu133plus2.db, probes, c("SYMBOL","ENTREZID"))
#Combine gene annotations with raw data
rma=merge(GeneAnno,rma,by='PROBEID')
#Write RMA-normalized, mapped data to file
write.table(rma, file = "rma.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
f1 <- kOverA(5, 200)
ffun <- filterfun(f1)
wh1 <- genefilter(exprs(sample.ExpressionSet), ffun)
sum(wh1)