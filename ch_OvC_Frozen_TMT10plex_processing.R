# TODO: Add comment
# 
# Author: cshughes
###############################################################################
library(RColorBrewer)
library(lattice)
library(gplots)
library(affy)
library(preprocessCore)
library(limma)
library(cluster)
library(vsn)
library(Biobase)
library(LSD)
###############################################################################
#grab the data
setwd(dir="/Users/cshughes/Documents/projects/OvC/frozen/")
ovpsm<-read.table("./ch_10Sept2015_OvC_TMT8_Frozen_HpH_repAB_PSMs.txt", header=TRUE, sep='\t')
ovpep<-read.table("./ch_10Sept2015_OvC_TMT8_Frozen_HpH_repAB_PeptideGroups.txt", header=TRUE, sep='\t')
ovpro<-read.table("./ch_10Sept2015_OvC_TMT8_Frozen_HpH_repAB_Proteins.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/frozen/Routput/")

##################################################
#initial processing of the data to compile peptides and remove contaminants
##################################################
#dont split replicates yet or do any combination of them
processPSM <- function(psmFile, proteinFile, ... ){	
	require(limma)
	#initial processing to remove non-informative columns and rows from contaminant proteins or those with many missing values
	pep<-psmFile[,c(5,6,9,10,11,24:26,29,36:46)]
	pro<-proteinFile[,c(3:4)]
	print('Raw data frame dimensions')
	print(dim(pep))
	#remove non-unique peptides or those without quan values
	pep<-subset(pep, Number.of.Proteins==1)
	pep<-subset(pep, !grepl('NoQuanValues',pep$Quan.Info))
	print('Unique Peptides Only')
	print(dim(pep))
	#parse the protein accession
	pep$Accession = sapply(strsplit(as.character(pep$Master.Protein.Accessions), ';'),'[', 1)
	#merge the protein descriptions into the peptide file
	pep.m = merge(pep,pro,by='Accession')
	#get the gene name out
	pep.m$Gene<-sub(".*?GN=(.*?)( .*|$)", "\\1", pep.m$Description)
	#parse the peptide column for amino acids
	pep.m$Sequence = toupper(sub('.*?\\.(.*?)(\\..*|$)','\\1',pep.m$Annotated.Sequence))
	#filter information from Description
	pep.m$Descriptions = sub('(.*?)( OS=.*|$)','\\1',pep.m$Description)
	#remove specific proteins
	pep.m = subset(pep.m, !grepl('Keratin',pep.m$Descriptions))
	pep.m = subset(pep.m, !grepl('sp',pep.m$Accession, ignore.case=FALSE))
	#filter the data columns
	print('removing contaminants')
	print(dim(pep.m))
	pep.r = pep.m[,c(1,23,25,24,12:21,10)]
	#assign the replicate
	repA = grepl('10Sept', pep.r$Spectrum.File)
	repB = grepl('18Sept', pep.r$Spectrum.File)
	pep.r$Rep = 'NA'
	pep.r[repA,16] = 'A'
	pep.r[repB,16] = 'B'
	#aggregate the PSMs into peptides
	print('aggregating peptides')
	pep.a = aggregate(cbind(X126,X127N,X127C,X128N,X128C,X129N,X129C,X130N)~Accession+Gene+Descriptions+Sequence+Rep,data=pep.r,mean,na.action=na.pass,na.rm=TRUE)
	colnames(pep.a) = c('Accession','Gene','Descriptions','Sequence','Rep','t1215a','t1246a','t1303a','t1493b','t376C','t1768a','t1798a','t1842a')
	print(dim(pep.a))
	#replace NaN with NA
	pep.a[,6:13] = round(pep.a[,6:13],2)
	pep.a[,6:13][is.na(pep.a[,6:13])]<-NA
	#output the data
	return(pep.a)	
}

ov.psm<-processPSM(ovpsm,ovpro)

##################################################
#do the peptide quant
##################################################
cols = colorRampPalette(c("white", brewer.pal(9,'YlOrRd')))
#split by replicates
repA = subset(ov.psm, grepl('A', ov.psm$Rep))
repB = subset(ov.psm, grepl('B', ov.psm$Rep))
#heat map functions
ovHeat_human <- function(heatIn,...){
	x = as.matrix(cor(heatIn[,6:13], use='pairwise.complete.obs'))
	#make the plot labels and boundaries
	xLabels<- names(heatIn[,6:13])
	mybreaks = seq(0.4,1,by=0.05) 
	#make the correlation heatmap
	heatmap.2(
			x,
			col= colorRampPalette(brewer.pal(9,"YlOrRd"))(length(mybreaks)-1),
			symkey=FALSE,
			Rowv=FALSE,
			Colv=FALSE,
			#hclust=hclustfunc,
			#distfun=distfunc,
			#na.color='black',
			#na.rm=TRUE,
			#symm=FALSE,
			dendrogram="none",
			breaks=mybreaks,
			cellnote = round(x,2),
			labRow = xLabels,
			labCol = xLabels,
			notecol = 'black',
			notecex = 0.5,
			colsep = 1:7,
			rowsep = 1:7,
			sepwidth = c(0.03, 0.03),
			sepcolor = 'white',
			#ColSideColors=ColSideColors,
			## labels
			main='correlation map',
			#xlab='Time Point',
			#ylab='Protein',
			## color key
			key = TRUE,
			keysize = 1,
			density.info = "none",
			scale = "none",
			trace = "none",
			mar=c(8,8),
			cexRow=1.2,
			cexCol=1.2
	)
}
#processing function
processPEP <- function(psmFile,...){
	require(vsn)
	require(LSD)
	require(RColorBrewer)
	#do work
	print('Raw data frame dimensions')
	print(dim(psmFile))
	#take out the ecoli proteins
	hu = subset(psmFile, rowSums(is.na(psmFile[,c(6:13)]))<7)
	print(paste('number of peptides remaining for human = ',nrow(hu),sep=''))
	#remove any PSMs with s/n < 5 on average across channels
	hu = subset(hu, rowMeans(hu[,c(6:13)], na.rm=TRUE) > 2)
	print('after S/N filtering')
	print(paste('number of peptides remaining for human = ',nrow(hu),sep=''))
	#plot the human replication
	pdf(paste('ch_Human_peptideReplicate_',hu$Rep[1],'_corHeat.pdf',sep=''))
	ovHeat_human(hu)
	dev.off()
	#do vsn normalization
	print('doing normalization')
	#hu = (length(spikeins)+1):nrow(huec)
	exprsFile1 = as.matrix(hu[,6:13])
	#spfit = vsn2(exprsFile1[spikeins,],lts.quantile=0.95)
	#exprsFile1[hu,] = normalizeBetweenArrays(exprsFile1[hu,], method = 'quantile')
	#nkid = as.data.frame(predict(spfit, newdata = exprsFile1))
	xnorm = as.data.frame(normalizeVSN(exprsFile1))
	colnames(xnorm)<-names(hu[,6:13])
	pep.q = cbind(hu[,c(1:5)],xnorm)
	#output the model fit
	pdf(paste('ch_PSMnormalization_Rep',hu$Rep[1],'_vsnFit.pdf',sep=''))
	meanSdPlot(exprsFile1, main='preNorm Fit')
	meanSdPlot(normalizeVSN(exprsFile1), main='postNorm Fit')
	M= log2(exprsFile1[,1]) - log2(exprsFile1[,2])
	A= (log2(exprsFile1[,1]) + log2(exprsFile1[,2]))/2
	smoothScatter(A,M,colramp = cols, cex = 2, xlab = 'Average Intensity',ylab = 'Intensity Ratio',main='preNorm expression')
	abline(h=0,col=brewer.pal(6,"RdBu")[6],lwd=3)
	box(lwd=3)
	M= log2(pep.q[,6]) - log2(pep.q[,7])
	A= (log2(pep.q[,6]) + log2(pep.q[,7]))/2
	smoothScatter(A,M,colramp = cols, cex = 2, xlab = 'Average Intensity',ylab = 'Intensity Ratio',main='postNorm expression')
	abline(h=0,col=brewer.pal(6,"RdBu")[6],lwd=3)
	box(lwd=3)
	dev.off()
	#plot the post normalized replication
	pdf(paste('ch_Human_peptideReplicate_',hu$Rep[1],'_NormCorHeat.pdf',sep=''))
	ovHeat_human(pep.q)
	dev.off()
	#output the data
	return(pep.q)	
}

#process the data
ov.pepA = processPEP(repA)
ov.pepB = processPEP(repB)


##################################################
#do the peptide stats
##################################################
#change the first line based on the data set you are using
pepQuan <- function(x,...){
	eset = as.matrix(x[,6:13])
	design <- cbind(SR=c(0,0,0,0,1,1,1,1),CC=c(1,1,1,1,0,0,0,0))
	fit <- lmFit(eset,design)
	cont.matrix <- makeContrasts(SRvsCC=SR-CC, levels=design)
	fit2 <- contrasts.fit(fit, cont.matrix)
	fit2 <- eBayes(fit2)
	x$logFC = fit2$coef
	x$pVal = fit2$p.value
	#output the data
	return(x[,c(1:4,14:15)])
}

ov.qA = pepQuan(ov.pepA)
ov.qB = pepQuan(ov.pepB)

#assign a peptide counter
ov.qA$pepNum = 1
ov.qB$pepNum = 1
#modify the colnames
colnames(ov.qA) = c('Accession','Gene','Descriptions','Sequence','logFC_a','pepPVal_a','pepNum')
colnames(ov.qB) = c('Accession','Gene','Descriptions','Sequence','logFC_b','pepPVal_b','pepNum')
#aggregate into proteins
proA = aggregate(cbind(logFC_a,pepPVal_a)~Accession+Gene+Descriptions,data=ov.qA,median,na.action=na.pass,na.rm=TRUE)
proB = aggregate(cbind(logFC_b,pepPVal_b)~Accession+Gene+Descriptions,data=ov.qB,median,na.action=na.pass,na.rm=TRUE)
#aggregate into proteins
proAp = aggregate(pepNum~Accession+Gene+Descriptions,data=ov.qA,sum,na.action=na.pass,na.rm=TRUE)
proBp = aggregate(pepNum~Accession+Gene+Descriptions,data=ov.qB,sum,na.action=na.pass,na.rm=TRUE)
#need to modify the colnames again for some reason
colnames(proA) = c('Accession','Gene','Descriptions','logFC_a','pepPVal_a')
colnames(proB) = c('Accession','Gene','Descriptions','logFC_b','pepPVal_b')
#merge
proA = merge(proA,proAp, by=c('Accession','Gene','Descriptions'))
proB = merge(proB,proBp, by=c('Accession','Gene','Descriptions'))
#combine into a single protein set
proTisF = merge(proA,proB,by=c('Accession','Gene','Descriptions'),sort=FALSE,all=TRUE)
proTisF$fcStd = rowMeans(proTisF[,c(4,7)], na.rm=TRUE)
proTisF$pvStd = rowMeans(proTisF[,c(5,8)], na.rm=TRUE)
proTisF$pepStd = rowMeans(proTisF[,c(6,9)], na.rm=TRUE)
#output the data table
write.table(proTisF,
		'ch_OvC-Frozen_TMT10plex_proteinSet.txt',
		quote=FALSE,
		sep='\t',
		col.names=TRUE,
		row.names=FALSE)



