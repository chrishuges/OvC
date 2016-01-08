# TODO: Add comment
# 
# Author: cshughes
###############################################################################

###############################################################################
#bring in the FFPE data to compare HGS with CCC using PECA
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
pec = inSet[,c(1:4,9,10,5,13,8,6,7,11)]
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
FFPE = prop[,c(1,7,8,9,2:6)]
colnames(FFPE) = c('Gene','Accession','Descriptions','pepNum','logFC','t','score','p-value','pVal')

###############################################################################
#bring in the Frozen data to compare HGS with CCC using PECA
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
FRO = prop[,c(1,7,8,9,2:6)]
colnames(FRO) = c('Gene','Accession','Descriptions','pepNum','logFC','t','score','p-value','pVal')



###############################################################################
#merge the FFPE and Frozen data
###############################################################################
f2 = merge(FFPE, FRO, by=c('Gene'))
#plot the data
#make colors
cols = colorRampPalette(c("white", brewer.pal(9,'YlGnBu')[6], brewer.pal(9,'YlOrRd')[9]))
#make a plot
pdf('ch_OvC_TMT10_FFPEvsFrozen_Human_Proteins_FoldChangeScatter.pdf')
reg = lm(f2$logFC.y~f2$logFC.x)
smoothScatter(f2$logFC.x,
		f2$logFC.y,
		main='FFPE vs Frozen Fold Change',
		xlab='FFPE',
		ylab='Frozen',
		xlim=c(-6,6),
		ylim=c(-6,6),
		colramp = cols,
		#pch = 19,
		cex = 2
)
text(3,-2,paste('n=',nrow(f2),sep=''))
text(3,-3,paste('r2 Set A =',round(cor(f2$logFC.x,f2$logFC.y,use='pairwise.complete.obs'),2),sep=''))
abline(reg,lty=2,lwd=3)
dev.off()









