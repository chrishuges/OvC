# TODO: Add comment
# 
# Author: cshughes
###############################################################################

###############################################################################
#compare the output fold change values for Human proteins
###############################################################################
library('PECA')
#grab the data
setwd(dir="/Users/cshughes/Documents/projects/OvC/Frozen/Routput") #change this to whatever directory you have stored the data in
a1h = readRDS('ch_OvC_Frozen_processedPeptides_fa1.rds')
a2h = readRDS('ch_OvC_Frozen_processedPeptides_fa2.rds')
b1h = readRDS('ch_OvC_Frozen_processedPeptides_fb1.rds')
b2h = readRDS('ch_OvC_Frozen_processedPeptides_fb2.rds')
#bind the two replicates
ah1 = rbind(a1h,a2h)
bh1 = rbind(b1h,b2h)
#make a function to calculate fold change
doQuant = function(x,...){
	#remove any e-coli peptides
	spikeins = grepl('EColi', x$Organism)
	inSet = x[!spikeins,c(1:2,4:5,7:14)]
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
	#output the data
	return(pecP)
}
#apply the function
ph1 = doQuant(ah1)
ph2 = doQuant(bh1)
#merge the two protein sets
pSet = merge(ph1,ph2,by=c('Accession','Gene','Descriptions'))
###plotting
#make colors
cols = colorRampPalette(c("white", brewer.pal(9,'YlGnBu')[6], brewer.pal(9,'YlOrRd')[9]))
#make a plot
pdf('ch_ch_OvC_TMT10_Frozen_Human_ProteinFoldChange_BiologicalReplicateScatter.pdf')
reg = lm(pSet$logFC.y~pSet$logFC.x)
smoothScatter(pSet$logFC.x,pSet$logFC.y,main='Human Peptide Biological Replicates Frozen',
		xlab='repA',
		ylab='repB',
		xlim=c(-6,6),
		ylim=c(-6,6),
		colramp = cols,
		cex = 2
)
text(3,-2,paste('n=',nrow(pSet),sep=''))
text(3,-3,paste('r2 Set A =',round(cor(pSet$logFC.x,pSet$logFC.y,use='pairwise.complete.obs'),2),sep=''))
abline(reg,lty=2,lwd=3)
dev.off()











