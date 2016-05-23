# TODO: Add comment
# 
# Author: cshughes
###############################################################################

#get the RLE set
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput") #change this to whatever directory you have stored the data in
pro = readRDS('ch_OvC_wStd_RLEset.rds')

##############################
#this is all starter code from the biomaRt users guide
#get some GO information
library(biomaRt)

#choose GO information
listMarts(host="www.ensembl.org")
ensembl=useMart("ENSEMBL_MART_ENSEMBL",host="www.ensembl.org")
ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host = "dec2015.archive.ensembl.org")

#commands
listFilters(ensembl) #attributes you want to obtain
listAttributes(ensembl) #attributes is what you want to retrieve (output)
#values is your input vector
getBM(attributes=c('affy_hg_u133_plus_2', 'entrezgene'), filters = 'affy_hg_u133_plus_2', values = affyids, mart = ensembl)

#for multiple GO and chromosome
go=c("GO:0051330","GO:0000080","GO:0000114")
chrom=c(17,20,"Y")
getBM(attributes= "hgnc_symbol",
		filters=c("go_id","chromosome_name"),
		values=list(go,chrom), mart=ensembl)

##############################
#ids
#GO:0002376 immune system process
#GO:0006955 immune response
#GO:0002684 positive regulation of immune system process
#GO:0050776 regulation of immune response
#GO:0002682 regulation of immune system process
#GO:0006954 inflammatory response
#GO:0050896, 0051869 response to stimulus


#grab for immune response
#go=c("GO:0006955",'GO:0002684','GO:0050776','GO:0002376')
go="GO:0006955"
ires = getBM(attributes= "external_gene_name",filters="go_id",values=go, mart=ensembl)
colnames(ires) = 'Gene'

#merge with the RLE data
p = merge(ires,pro,by='Gene')



#built the heatmap
pdf('ch_OvC_wStd_Proteins-RLE_ImmuneClustering.pdf')
#make the plot labels and boundaries
xLabels<- c(rep('hgsc',6),rep('ccc',6),rep('enoc',6))
mybreaks = seq(-2,2,by=0.05)
ColSideColors = c(rep(brewer.pal(6,'Accent')[1],6),rep(brewer.pal(6,'Accent')[2],6),rep(brewer.pal(6,'Accent')[3],6))
#make the correlation heatmap
heatmap.2(
		as.matrix(x[2:19]),
		col= colorRampPalette(brewer.pal(6,"RdBu"))(length(mybreaks)-1),
		symkey=TRUE,
		Rowv=TRUE,
		Colv=FALSE,
		dendrogram="row",
		breaks=mybreaks,
		#cellnote = round(ovCor,2),
		labRow = '',
		labCol = xLabels,
		las=2,
		ColSideColors=ColSideColors,
		## labels
		main='Immune Clustering',
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


#PCA
p$hgs1 = rowMeans(p[,c(2:7)],na.rm=TRUE)
p$ccc1 = rowMeans(p[,c(8:13)],na.rm=TRUE)
p$emc1 = rowMeans(p[,c(14:19)],na.rm=TRUE)
p$var = apply(p[,20:22], 1, function(x) var(x, na.rm=TRUE))
#order
proV = p[order(-p$var),]
x = proV[rowSums(is.na(proV[,2:19]))<1,c(2:19)]
row.names(x) = proV[rowSums(is.na(proV[,2:19]))<1,1]

pca <- prcomp(t(x))
hCols = brewer.pal(6,'Accent')
pcCols = c(rep(hCols[1],6),rep(hCols[2],6),rep(hCols[3],6))
pcPCH = c(rep(15,3),rep(19,3),rep(8,3),rep(12,3),rep(2,3),rep(11,3))
pdf('ch_OvC_TMT10_wStd_Human_Immune_PCA.pdf')
plot(pca$x[,1:2],
		col=pcCols,
		pch=pcPCH,
		cex = 2
)
box(lwd=3)
dev.off()

