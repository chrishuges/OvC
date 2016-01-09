# TODO: Add comment
# 
# Author: cshughes
###############################################################################


##################################################
#read in the human protein atlas data, for proteins
##################################################
setwd(dir="/Users/cshughes/Documents/projects/PaC/RNAseq/")
ens<-read.table("./HumanProteomeMap_Science_TableS1.txt", header=TRUE, sep='\t')
hpa<-read.table("./HumanProteinAtlas_ProteinExp.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/PaC/Routput/")
#subset out only pancreatic tissue samples
hpa.p = hpa[grepl('pancreatic',hpa$Tumor),]
#get the proteins with high protein expression
hpa.ph = hpa.p[grepl('High',hpa.p$Level),]
#keep only those where greater than 50% of total patients were classified as High
hpa.ph$pHigh = hpa.ph$Count.patients/hpa.ph$Total.patients
hpa.phs = hpa.ph[hpa.ph$pHigh>=0.6,]
#merge with the gene expression data
colnames(hpa.phs)[1] = 'ensg_id'
hpa.rp = merge(hpa.phs,hpa.mi,by='ensg_id')
#plot the FPKM overlaid with the protein data
cols = brewer.pal(6,'Set2')
pdf('ch_HumanProteinAtlas_FPKMwAtlasPro_Pancreas.pdf')
hist(log2(hpa.rp[,62]),
		col = cols[6],
		breaks=75,
		xlim = c(-10,15),
		main = 'Pancreatic Tissue FPKM Distribution',
		xlab = 'log2(FPKM Pancreatic Tissue)',
		ylab = 'Frequency'
)
hist(log2(hpa.mi[,56]),add=TRUE,breaks=200,col=cols[3])
hist(log2(hpa.rp[,62]),add=TRUE,breaks=75,col=cols[6])
dev.off()

