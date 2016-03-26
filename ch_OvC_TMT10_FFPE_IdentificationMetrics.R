# TODO: Add comment
# 
# Author: cshughes
###############################################################################
setwd(dir="/Users/cshughes/Documents/projects/OvC/FFPE/Routput") #change this to whatever directory you have stored the data in
a1h = readRDS('ch_OvC_FFPE_normalizedPeptides_fa1.rds')
a2h = readRDS('ch_OvC_FFPE_normalizedPeptides_fa2.rds')
b1h = readRDS('ch_OvC_FFPE_normalizedPeptides_fb1.rds')
b2h = readRDS('ch_OvC_FFPE_normalizedPeptides_fb2.rds')

##########################################################
##first need to make the protein set
##########################################################
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
pec = inSet[,c(1:4,9,14,10,5,13,8,6,7,12,11)]
#calculate p-values
#assign the groups to compare for stats
group1<-colnames(pec)[c(5:9)]
group2<-colnames(pec)[c(10:14)]
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
#write out the data into a format usable in excel
write.table(pecP,'ch_OvC_FFPE_HGSvCCC_PECA_proteinSet.txt',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)

##################################################
#using the human protein atlas data for FPKM correlation
##################################################
setwd(dir="/Users/cshughes/Documents/projects/OvC/RNAexpression/")
ens<-read.table("./HumanProteomeMap_Science_TableS1.txt", header=TRUE, sep='\t')
hpa<-read.table("./HumanProteinAtlas_FPKM_allsamples.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/FFPE/Routput/")
pecP = readRDS('ch_OvC_FFPE_proteinSet.rds')
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
pdf('ch_OvC_TMT10_FFPE_Human_Proteins_FPKM_HumanProteinAtlasCoverage.pdf')
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




##################################################
#make a building Venn to look at ID saturation
##################################################
setwd(dir="/Users/cshughes/Documents/projects/OvC/FFPE/Routput") #change this to whatever directory you have stored the data in
a1h = readRDS('ch_OvC_FFPE_normalizedPeptides_fa1.rds')
a2h = readRDS('ch_OvC_FFPE_normalizedPeptides_fa2.rds')
b1h = readRDS('ch_OvC_FFPE_normalizedPeptides_fb1.rds')
b2h = readRDS('ch_OvC_FFPE_normalizedPeptides_fb2.rds')
#need a protein set with all patients
#bind all data
allh = rbind(a1h,a2h,b1h,b2h)
pepSet = aggregate(cbind(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)~Accession+Gene+Descriptions+Sequence,data=allh,median,na.action=na.pass,na.rm=TRUE)
proSet = aggregate(cbind(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)~Accession+Gene+Descriptions,data=pepSet,median,na.action=na.pass,na.rm=TRUE)
#look for column with most NAs
summary(proSet)
#sample a6 has the most NAs
t1na = is.na(proSet[,9])
t1 = proSet[!t1na,2]
length(t1)#8035
t2na = is.na(proSet[,7])
t2 = c(t1,proSet[!t2na,2])
length(unique(t2))#8109
t3na = is.na(proSet[,12])
t3 = c(t2,proSet[!t3na,2])
length(unique(t3))#8137
t4na = is.na(proSet[,4])
t4 = c(t3,proSet[!t4na,2])
length(unique(t4))#8163
t5na = is.na(proSet[,8])
t5 = c(t4,proSet[!t5na,2])
length(unique(t5))#8164
t6na = is.na(proSet[,10])
t6 = c(t5,proSet[!t6na,2])
length(unique(t6))#8165
t7na = is.na(proSet[,6])
t7 = c(t6,proSet[!t7na,2])
length(unique(t7))#8167
t8na = is.na(proSet[,11])
t8 = c(t7,proSet[!t8na,2])
length(unique(t8))#8167
t9na = is.na(proSet[,5])
t9 = c(t8,proSet[!t9na,2])
length(unique(t9))#8167
t10na = is.na(proSet[,13])
t10 = c(t9,proSet[!t10na,2])
length(unique(t10))#8167


#heat map for IDs
#use proSet
#make a counter for peptide numbers
pepSet$pepNum = 1
proSet2 = aggregate(pepNum~Accession+Gene+Descriptions,data=pepSet,sum,na.action=na.pass,na.rm=TRUE)
lset = cbind(proSet,proSet2)
lset = lset[,c(1:3,17,8,13,9,4,12,7,5,6,11,10)]
lset$missing = ifelse(rowSums(is.na(lset[5:14]))>0,'yes','no')
lset = lset[order(-lset$pepNum),]
xnorm = lset[lset$pepNum==1,5:14]
xnorm = xnorm[order(rowSums(is.na(xnorm))),]
xnorm[xnorm>0] = 1.1
xnorm[is.na(xnorm)] = 0

#make the plot
mybreaks = seq(0,2,by=0.1)
xLabels<- c('hgs1','hgs2','hgs3','hgs4','hgs5','ccc1','ccc2','ccc3','ccc4','ccc5')
ColSideColors = c(rep(brewer.pal(6,'Accent')[1],5),rep(brewer.pal(6,'Accent')[3],5))
pdf('ch_OvC_TMT10_FFPE_Human_Proteins_IDheatmap_Sub.pdf')
heatmap.2(
		as.matrix(xnorm),
		col= colorRampPalette(brewer.pal(11,"RdBu"))(length(mybreaks)-1),
		symkey=FALSE,
		Rowv=FALSE,
		Colv=FALSE,
		dendrogram="none",
		breaks=mybreaks,
		#cellnote = round(ovCor,2),
		labRow = '',
		labCol = xLabels,
		las=2,
		ColSideColors=ColSideColors,
		colsep = 1:10,
		rowsep = 1:10,
		sepwidth = c(0.03, 0.03),
		sepcolor = 'white',
		## labels
		main='test',
		## color key
		key = FALSE,
		keysize = 1,
		density.info = "none",
		scale = "none",
		trace = "none",
		mar=c(8,8),
		cexRow=1.5,
		cexCol=1.5
)
dev.off()


##########################################################
##make a plot for peptide metrics
##########################################################
#remove single peptide hits from the plot
sPeps = subset(pecP, pepNum>0)
#subset out peptides with greater than 50 peptides assigned
nPeps<-data.frame(table(sPeps$pepNum))
up<-subset(nPeps, as.numeric(as.character(nPeps$Var1))>25)
dn<-subset(nPeps, as.numeric(as.character(nPeps$Var1))<=25)
#sum the number of proteins id with more than 50 peptides
x<-sum(up$Freq)
dn[26,1]<-'26'
dn[26,2]<-x
#
mainCol<-brewer.pal(6,"RdBu")
colA<-col2rgb(mainCol[6])
xLabels<-seq(1,26,1)
xValues<-seq(1,26,1)
#make the plot
mp<-barplot(dn[,2])
pdf('ch_OvC_TMT10_FFPE_Human_Proteins_NumberOfPeptides_Barchart.pdf')
barplot(dn[,2],
		col=rgb(colA[1,],colA[2,],colA[3,],100,maxColorValue=255),
		main="Distribution of Peptide Identification Numbers", 
		xlab="Number of Unique Peptides",
		ylim=c(0,max(dn[,2])+50),
		xaxt="n"
)
box(lwd=3)
axis(1,mp,xLabels,las=2,lwd=2)
dev.off()





