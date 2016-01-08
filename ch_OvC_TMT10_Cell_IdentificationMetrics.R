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
setwd(dir="/Users/cshughes/Documents/projects/OvC/CellLine/Routput") #change this to whatever directory you have stored the data in
a1h = readRDS('ch_OvC_CellLine_processedPeptides_fa1.rds')
a2h = readRDS('ch_OvC_CellLine_processedPeptides_fa2.rds')
#bind all data
allh = rbind(a1h,a2h)
aggh = aggregate(cbind(a1,a2,a3,a4,a5,a6)~Accession+Gene+Descriptions+Sequence+Organism,data=allh,median,na.action=na.pass,na.rm=TRUE)
#remove any e-coli peptides
spikeins = grepl('EColi', aggh$Organism)
inSet = aggh[!spikeins,c(1:4,6:11)]
#need to reshuffle the patients based on HGS or CCC
pec = inSet[,c(1:4,5:10)]
#calculate p-values
#assign the groups to compare for stats
group1<-colnames(pec)[c(5:7)]
group2<-colnames(pec)[c(8:10)]
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


##################################################
#using the human protein atlas data for FPKM correlation
##################################################
setwd(dir="/Users/cshughes/Documents/projects/OvC/RNAexpression/")
ens<-read.table("./HumanProteomeMap_Science_TableS1.txt", header=TRUE, sep='\t')
hpa<-read.table("./HumanProteinAtlas_FPKM_allsamples.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/CellLine/Routput/")
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
pdf('ch_OvC_TMT10_CellLine_Human_Proteins_FPKM_HumanProteinAtlasCoverage.pdf')
hist(log2(hpa.mi[,130]),
		col = cols[3],
		breaks=200,
		xlim = c(-10,15),
		main = 'Ovarian Tissue FPKM Distribution',
		xlab = 'log2(FPKM Pancreatic Tissue)',
		ylab = 'Frequency'
)
hist(log2(ov.whpa[,138]),add=TRUE,breaks=100,col=cols[6])
dev.off()







