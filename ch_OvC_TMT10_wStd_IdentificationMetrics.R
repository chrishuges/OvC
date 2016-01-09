# TODO: Add comment
# 
# Author: cshughes
###############################################################################
##########################################################
##first need to make the protein set
##########################################################
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput") #change this to whatever directory you have stored the data in
a1h = readRDS('ch_OvC_wStd_processedPeptides_a1.rds')
a2h = readRDS('ch_OvC_wStd_processedPeptides_a2.rds')
a3h = readRDS('ch_OvC_wStd_processedPeptides_a3.rds')
b1h = readRDS('ch_OvC_wStd_processedPeptides_b1.rds')
b2h = readRDS('ch_OvC_wStd_processedPeptides_b2.rds')
b3h = readRDS('ch_OvC_wStd_processedPeptides_b3.rds')
#bind all data
allh = rbind(a1h,a2h,a3h,b1h,b2h,b3h)
aggh = aggregate(cbind(a1,a2,a3,a4,a5,a6,a7,a8,a9,aStd)~Accession+Gene+Descriptions+Sequence+Organism,data=allh,median,na.action=na.pass,na.rm=TRUE)
#remove any e-coli peptides
spikeins = grepl('EColi', aggh$Organism)
inSet = aggh[!spikeins,c(1:4,6:15)]
#no need to reshuffle, just want IDs
pec = inSet[,c(1:4,5:14)]
pec$pepNum = 1
#aggregate into proteins
pecP = aggregate(pepNum~Accession+Gene,data=pec,sum,na.action=na.pass,na.rm=TRUE)

##################################################
#using the human protein atlas data for FPKM correlation
##################################################
setwd(dir="/Users/cshughes/Documents/projects/OvC/RNAexpression/")
ens<-read.table("./HumanProteomeMap_Science_TableS1.txt", header=TRUE, sep='\t')
hpa<-read.table("./HumanProteinAtlas_FPKM_allsamples.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput/")
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
#combine with the protein data
ov.whpa = merge(pecP,hpa.mi,by='Gene',all=FALSE,sort=FALSE)
#map to the protein data
cols = brewer.pal(6,'Set2')
pdf('ch_OvC_TMT10_wStd_Human_Proteins_FPKM_HumanProteinAtlasCoverage.pdf')
hist(log2(hpa.mi[,130]),
		col = cols[3],
		breaks=200,
		xlim = c(-10,15),
		main = 'Ovarian Tissue FPKM Distribution',
		xlab = 'log2(FPKM Pancreatic Tissue)',
		ylab = 'Frequency'
)
hist(log2(ov.whpa[,133]),add=TRUE,breaks=150,col=cols[6])
dev.off()



