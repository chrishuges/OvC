# TODO: Add comment
# 
# Author: cshughes
###############################################################################
#read in the data objects
setwd(dir="/Users/cshughes/Documents/projects/OvC/FFPE/Routput") #change this to whatever directory you have stored the data in
a1h = readRDS('ch_OvC_FFPE_normalizedPeptides_fa1.rds')
a2h = readRDS('ch_OvC_FFPE_normalizedPeptides_fa2.rds')
b1h = readRDS('ch_OvC_FFPE_normalizedPeptides_fb1.rds')
b2h = readRDS('ch_OvC_FFPE_normalizedPeptides_fb2.rds')

setwd(dir="/Users/cshughes/Documents/projects/OvC/FFPE/Routput") #change this to whatever directory you have stored the data in
a1e = readRDS('ch_OvC_FFPE_normalizedPeptides_EColi_fa1.rds')
a2e = readRDS('ch_OvC_FFPE_normalizedPeptides_EColi_fa2.rds')
b1e = readRDS('ch_OvC_FFPE_normalizedPeptides_EColi_fb1.rds')
b2e = readRDS('ch_OvC_FFPE_normalizedPeptides_EColi_fb2.rds')

###############################################################################
#look first at how the e-coli samples deviated
###############################################################################
a12e = merge(a1e[,c(1,2,4,5,8:17)],a2e[,c(1,2,4,5,8:17)],by=c('Accession','Gene','Descriptions','Sequence'))
b12e = merge(b1e[,c(1,2,4,5,8:17)],b2e[,c(1,2,4,5,8:17)],by=c('Accession','Gene','Descriptions','Sequence'))

###make an RLE plot
#merge all replicates
ab12e = merge(a12e,b12e,by=c('Accession','Gene','Descriptions','Sequence'))
#take the deviation from the median
ab12e$RLEmed = apply(ab12e[,5:44],1, function(x) median(x))
vnorm = apply(ab12e[,5:44],2, function(x) x - ab12e$RLEmed)
#make colors for the batches
cols = c(rep(brewer.pal(6,'RdBu')[1],10),rep(brewer.pal(6,'RdBu')[2],10),rep(brewer.pal(6,'RdBu')[5],10),rep(brewer.pal(6,'RdBu')[6],10))
#make the plot
pdf('ch_OvC_TMT10_FFPE_EColi_Peptides_RLEplot.pdf')
boxplot(vnorm,
		col=cols,
		pch=20,
		cex=0.75,
		las=2,
		xaxt='n',
		ylab = 'Relative Log Expression',
		main = 'RLE plot of EColi Batch Deviation',
		ylim = c(-3.5,3.5))
box(lwd=3)
dev.off()


###############################################################################
#compare technical replicates of Human peptides
###############################################################################
a12h = merge(a1h[,c(1,2,4,5,8:17)],a2h[,c(1,2,4,5,8:17)],by=c('Accession','Gene','Descriptions','Sequence'))
b12h = merge(b1h[,c(1,2,4,5,8:17)],b2h[,c(1,2,4,5,8:17)],by=c('Accession','Gene','Descriptions','Sequence'))
###Replicate A
#make a median of all values for one replicate
a12h$med1 = apply(a12h[,5:14],1, function(x) median(x,na.rm=TRUE))
a12h$med2 = apply(a12h[,15:24],1, function(x) median(x,na.rm=TRUE))
vnormA1 = apply(a12h[,5:14],2, function(x) x - a12h$med1)
vnormA2 = apply(a12h[,15:24],2, function(x) x - a12h$med2)
vnormA = cbind(vnormA1,vnormA2)
#take a look at the correlations
trepA = ''
for (i in 1:10){
	x = cor(vnormA[,i],vnormA[,i+10],use='pairwise.complete.obs',method='pearson')
	trepA[i] = as.numeric(x)
}
#for technical replicates A
mean(as.numeric(trepA)) #0.8436

###Replicate B
#make a median of all values for one replicate
b12h$med1 = apply(b12h[,5:14],1, function(x) median(x,na.rm=TRUE))
b12h$med2 = apply(b12h[,15:24],1, function(x) median(x,na.rm=TRUE))
vnormB1 = apply(b12h[,5:14],2, function(x) x - b12h$med1)
vnormB2 = apply(b12h[,15:24],2, function(x) x - b12h$med2)
vnormB = cbind(vnormB1,vnormB2)
#take a look at the correlations
trepB = ''
for (i in 1:10){
	x = cor(vnormB[,i],vnormB[,i+10],use='pairwise.complete.obs',method='pearson')
	trepB[i] = as.numeric(x)
}
#for technical replicates B
mean(as.numeric(trepB)) #0.8553

#make colors
cols = colorRampPalette(c("white", brewer.pal(9,'YlGnBu')[6], brewer.pal(9,'YlOrRd')[9]))
#make a plot
#this is too crazy...make it a smoothscatter
pdf('ch_ch_OvC_TMT10_FFPE_Human_Peptides_TechnicalReplicateScatter.pdf')
reg = lm(vnormA[,14]~vnormA[,4])
smoothScatter(vnormA[,4],vnormA[,14],
		main='Human Peptide Technical Replicates FFPE',
		xlab='repA',
		ylab='repB',
		xlim=c(-6,6),
		ylim=c(-6,6),
		colramp = cols,
		#pch = 19,
		cex = 2
)
text(3,-2,paste('n=',nrow(vnormA),sep=''))
text(3,-3,paste('r2 Set A =',round(mean(as.numeric(trepA)),2),sep=''))
text(3,-4,paste('r2 Set B =',round(mean(as.numeric(trepB)),2),sep=''))
abline(reg,lty=2,lwd=3)
dev.off()




###############################################################################
#compare biological replicates of Human peptides
###############################################################################
abh1 = merge(a1h[,c(1,2,4,5,8:17)],b1h[,c(1,2,4,5,8:17)],by=c('Accession','Gene','Descriptions','Sequence'))
abh2 = merge(a2h[,c(1,2,4,5,8:17)],b2h[,c(1,2,4,5,8:17)],by=c('Accession','Gene','Descriptions','Sequence'))
###Replicate A
#make a median of all values for one replicate
abh1$med1 = apply(abh1[,5:14],1, function(x) median(x,na.rm=TRUE))
abh1$med2 = apply(abh1[,15:24],1, function(x) median(x,na.rm=TRUE))
vnormA1 = apply(abh1[,5:14],2, function(x) x - abh1$med1)
vnormA2 = apply(abh1[,15:24],2, function(x) x - abh1$med2)
vnormA = cbind(vnormA1,vnormA2)
#take a look at the correlations
brepA = ''
for (i in 1:10){
	x = cor(vnormA[,i],vnormA[,i+10],use='pairwise.complete.obs',method='pearson')
	brepA[i] = as.numeric(x)
}
#for biological replicates A
mean(as.numeric(brepA)) #0.7559

###Replicate B
#make a median of all values for one replicate
abh2$med1 = apply(abh2[,5:14],1, function(x) median(x,na.rm=TRUE))
abh2$med2 = apply(abh2[,15:24],1, function(x) median(x,na.rm=TRUE))
vnormB1 = apply(abh2[,5:14],2, function(x) x - abh2$med1)
vnormB2 = apply(abh2[,15:24],2, function(x) x - abh2$med2)
vnormB = cbind(vnormB1,vnormB2)
#take a look at the correlations
brepB = ''
for (i in 1:10){
	x = cor(vnormB[,i],vnormB[,i+10],use='pairwise.complete.obs',method='pearson')
	brepB[i] = as.numeric(x)
}
#for biological replicates B
mean(as.numeric(brepB)) #0.7539

#make colors
cols = colorRampPalette(c("white", brewer.pal(9,'YlGnBu')[6], brewer.pal(9,'YlOrRd')[9]))
#make a plot
pdf('ch_ch_OvC_TMT10_FFPE_Human_Peptides_BiologicalReplicateScatter.pdf')
reg = lm(vnormA[,19]~vnormA[,9])
smoothScatter(vnormA[,9],vnormA[,19],main='Human Peptide Biological Replicates FFPE',
		xlab='repA',
		ylab='repB',
		xlim=c(-6,6),
		ylim=c(-6,6),
		colramp = cols,
		cex = 2
)
text(3,-2,paste('n=',nrow(vnormA),sep=''))
text(3,-3,paste('r2 Set A =',round(mean(as.numeric(brepA)),2),sep=''))
text(3,-4,paste('r2 Set B =',round(mean(as.numeric(brepB)),2),sep=''))
abline(reg,lty=2,lwd=3)
dev.off()


###############################################################################
#compare the output fold change values for Human proteins
###############################################################################
setwd(dir="/Users/cshughes/Documents/projects/OvC/FFPE/Routput") #change this to whatever directory you have stored the data in
a1h = readRDS('ch_OvC_FFPE_processedPeptides_fa1.rds')
a2h = readRDS('ch_OvC_FFPE_processedPeptides_fa2.rds')
b1h = readRDS('ch_OvC_FFPE_processedPeptides_fb1.rds')
b2h = readRDS('ch_OvC_FFPE_processedPeptides_fb2.rds')
#bind the two replicates
ah1 = rbind(a1h,a2h)
bh1 = rbind(b1h,b2h)
#make a function to calculate fold change
doQuant = function(x,...){
	#remove any e-coli peptides
	spikeins = grepl('EColi', x$Organism)
	inSet = x[!spikeins,c(1:2,4:5,7:16)]
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
pdf('ch_ch_OvC_TMT10_FFPE_Human_ProteinFoldChange_BiologicalReplicateScatter.pdf')
reg = lm(pSet$logFC.y~pSet$logFC.x)
smoothScatter(pSet$logFC.x,pSet$logFC.y,main='Human Peptide Biological Replicates FFPE',
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






