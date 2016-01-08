# TODO: processing the TMT cell line data for OvC
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
setwd(dir="/Users/cshughes/Documents/projects/OvC/cell/")
ovpsm<-read.table("./ch_08July2015_OvCcell-TMT6_HpH_All_PSMs.txt", header=TRUE, sep='\t')
ovpep<-read.table("./ch_08July2015_OvCcell-TMT6_HpH_All_PeptideGroups.txt", header=TRUE, sep='\t')
ovpro<-read.table("./ch_08July2015_OvCcell-TMT6_HpH_All_Proteins.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/cell/Routput/")

##################################################
#initial processing of the data to compile peptides and remove contaminants
##################################################
#dont split replicates yet or do any combination of them
processPSM <- function(psmFile, proteinFile, ... ){	
	require(limma)
	#initial processing to remove non-informative columns and rows from contaminant proteins or those with many missing values
	pep<-psmFile[,c(5,6,9,10,24:26,30,34,36:42)]
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
	pep.r = pep.m[,c(1,19,21,20,3:4,11:17,7:10)]
	#aggregate the PSMs into peptides
	print('aggregating peptides')
	pep.a = aggregate(cbind(X126,X127,X128,X129,X130,X131)~Accession+Gene+Descriptions+Sequence,data=pep.r,mean,na.action=na.pass,na.rm=TRUE)
	colnames(pep.a) = c('Accession','Gene','Descriptions','Sequence','OVCAR3','OVCAR5','OVSAHO','JHOC5','OVISE','OVTOKO')
	print(dim(pep.a))
	#replace NaN with NA
	pep.a[,5:10] = round(pep.a[,5:10],2)
	pep.a[,5:10][is.na(pep.a[,5:10])]<-NA
	#output the data
	return(pep.a)	
}

ov.psm<-processPSM(ovpsm,ovpro)

##################################################
#do the peptide quant
##################################################
cols = colorRampPalette(c("white", brewer.pal(9,'YlOrRd')))
#heat map functions
ovHeat_human <- function(heatIn,...){
	x = as.matrix(cor(heatIn[,5:10], use='pairwise.complete.obs'))
	#make the plot labels and boundaries
	xLabels<- names(heatIn[,5:10])
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
			colsep = 1:5,
			rowsep = 1:5,
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
	#remove rows with too many NA
	hu = subset(psmFile, rowSums(is.na(psmFile[,c(5:10)]))<4)
	print(paste('number of peptides remaining for human = ',nrow(hu),sep=''))
	#remove any PSMs with s/n < 5 on average across channels
	hu = subset(hu, rowMeans(hu[,c(5:10)], na.rm=TRUE) > 2)
	print('after S/N filtering')
	print(paste('number of peptides remaining for human = ',nrow(hu),sep=''))
	#plot the human replication
	pdf('ch_Human_peptideReplicate_Cells_corHeat.pdf')
	ovHeat_human(hu)
	dev.off()
	#do vsn normalization
	print('doing normalization')
	#hu = (length(spikeins)+1):nrow(huec)
	exprsFile1 = as.matrix(hu[,5:10])
	#spfit = vsn2(exprsFile1[spikeins,],lts.quantile=0.95)
	#exprsFile1[hu,] = normalizeBetweenArrays(exprsFile1[hu,], method = 'quantile')
	#nkid = as.data.frame(predict(spfit, newdata = exprsFile1))
	xnorm = as.data.frame(normalizeVSN(exprsFile1))
	colnames(xnorm)<-names(hu[,5:10])
	pep.q = cbind(hu[,c(1:4)],xnorm)
	#output the model fit
	pdf('ch_PSMnormalization_vsnFit.pdf')
	meanSdPlot(exprsFile1, main='preNorm Fit')
	meanSdPlot(normalizeVSN(exprsFile1), main='postNorm Fit')
	M= log2(exprsFile1[,1]) - log2(exprsFile1[,2])
	A= (log2(exprsFile1[,1]) + log2(exprsFile1[,2]))/2
	smoothScatter(A,M,colramp = cols, cex = 2, xlab = 'Average Intensity',ylab = 'Intensity Ratio',main='preNorm expression')
	abline(h=0,col=brewer.pal(6,"RdBu")[6],lwd=3)
	box(lwd=3)
	M= log2(pep.q[,7]) - log2(pep.q[,8])
	A= (log2(pep.q[,7]) + log2(pep.q[,8]))/2
	smoothScatter(A,M,colramp = cols, cex = 2, xlab = 'Average Intensity',ylab = 'Intensity Ratio',main='postNorm expression')
	abline(h=0,col=brewer.pal(6,"RdBu")[6],lwd=3)
	box(lwd=3)
	dev.off()
	#plot the post normalized replication
	pdf('ch_Human_peptideReplicate_NormCorHeat.pdf')
	ovHeat_human(pep.q)
	dev.off()
	#output the data
	return(pep.q)	
}

#process the data
ov.pep = processPEP(ov.psm)


##################################################
#do the peptide stats
##################################################
#JHOC5 looks fuuuuuuuucked...really similar to the serous lines, exclude it from the stats
#change the first line based on the data set you are using
pepQuan <- function(x,...){
	eset = as.matrix(x[,c(5:7,9:10)])
	design <- cbind(SR=c(1,1,1,0,0),CC=c(0,0,0,1,1))
	fit <- lmFit(eset,design)
	cont.matrix <- makeContrasts(SRvsCC=SR-CC, levels=design)
	fit2 <- contrasts.fit(fit, cont.matrix)
	fit2 <- eBayes(fit2)
	x$logFC = fit2$coef
	x$pVal = fit2$p.value
	#output the data
	return(x[,c(1:4,11:12)])
}

ov.q = pepQuan(ov.pep)


#assign a peptide counter
ov.q$pepNum = 1
#modify the colnames
colnames(ov.q) = c('Accession','Gene','Descriptions','Sequence','logFC','pepPVal','pepNum')
#aggregate into proteins
proCell = aggregate(cbind(logFC,pepPVal)~Accession+Gene+Descriptions,data=ov.q,median,na.action=na.pass,na.rm=TRUE)
#aggregate into proteins
proCellp = aggregate(pepNum~Accession+Gene+Descriptions,data=ov.q,sum,na.action=na.pass,na.rm=TRUE)
#need to modify the colnames again for some reason
colnames(proCell) = c('Accession','Gene','Descriptions','logFC','pepPVal')
#merge
proCell = merge(proCell,proCellp, by=c('Accession','Gene','Descriptions'))

#output the data table
write.table(proCell,
		'ch_OvCcell_TMT10plex_proteinSet.txt',
		quote=FALSE,
		sep='\t',
		col.names=TRUE,
		row.names=FALSE)




