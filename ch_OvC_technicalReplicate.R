# TODO: Add comment
# 
# Author: cshughes
###############################################################################
#grab the data
setwd(dir="/Users/cshughes/Documents/projects/OvC/")
ovpsm<-read.table("./ch_09July2015_OvC-TMT10_HpH_repAB_PSMs.txt", header=TRUE, sep='\t')
ovpep<-read.table("./ch_09July2015_OvC-TMT10_HpH_repAB_PeptideGroups.txt", header=TRUE, sep='\t')
ovpro<-read.table("./ch_09July2015_OvC-TMT10_HpH_repAB_Proteins.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/Routput/")

##################################################
#determine the average number of PSMs per peptide in each case
##################################################
processPSM <- function(psmFile, proteinFile, ... ){	
	require(limma)
	#initial processing to remove non-informative columns and rows from contaminant proteins or those with many missing values
	pep<-psmFile[,c(5,6,8,10,11,25:27,30,34,37:47)]
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
	pep.m = subset(pep.m, !grepl('Keratin',pep.m$Descriptions) & !grepl('sp',pep.m$Accession, ignore.case=FALSE))
	#filter the data columns
	print('removing contaminants')
	print(dim(pep.m))
	pep.r = pep.m[,c(1,24,26,25,3:4,11:22,7:10)]
	#assign the replicate
	repA = grepl('repA', pep.r$Spectrum.File)
	repB = grepl('repB', pep.r$Spectrum.File)
	pep.r$Rep = 'NA'
	pep.r[repA,23] = 'A'
	pep.r[repB,23] = 'B'
	#assign the organism
	ec = grepl('True', pep.r$Std)
	print(paste('number of e coli peptides = ',table(ec)[2],sep=''))
	hu = grepl('False', pep.r$Std)
	print(paste('number of human peptides = ',table(hu)[2],sep=''))
	pep.r$Source = 'NA'
	pep.r[ec,24] = 'ecoli'
	pep.r[hu,24] = 'human'
	#output the data
	return(pep.r)	
}
#process the data
ov.psm<-processPSM(ovpsm,ovpro)

#make a spectrum counter
ov.psm$specNum = 1
#eliminate e-coli peptides
ecoli = grepl('coli',ov.psm$Source)
ov.psmh = ov.psm[!ecoli,] 
#aggregate to get the total number of psms per peptide
ov.a = aggregate(specNum~Accession+Gene+Descriptions+Sequence+Rep, data=ov.psmh, sum, na.action=na.pass, na.rm=TRUE)
#get the rows for replicate A
repA = grepl('A', ov.a$Rep)
#get the mean of the specNum column for replicate A
mA = mean(ov.a[repA,6], na.rm=TRUE)
#2.74
#aggregate into a single total data set with both replicates
ov.a2 = aggregate(specNum~Accession+Gene+Descriptions+Sequence, data=ov.psmh, sum, na.action=na.pass, na.rm=TRUE)
#get the mean of the specNum column for replicate A
mAB = mean(ov.a2[,5], na.rm=TRUE)
#4.47

#plot the data
pdf('ch_OvC_AverageNumSpectra_perPeptide_ctpTissue.pdf')
boxplot(log2(ov.a$specNum),
		log2(ov.a2$specNum)
)
dev.off()



##################################################
#initial processing of the data to compile peptides and remove contaminants
##################################################
#dont split replicates yet or do any combination of them
processPSM <- function(psmFile, proteinFile, ... ){	
	require(limma)
	#initial processing to remove non-informative columns and rows from contaminant proteins or those with many missing values
	pep<-psmFile[,c(5,6,8,10,11,25:27,30,34,37:47)]
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
	pep.m = subset(pep.m, !grepl('Keratin',pep.m$Descriptions) & !grepl('sp',pep.m$Accession, ignore.case=FALSE))
	#filter the data columns
	print('removing contaminants')
	print(dim(pep.m))
	pep.r = pep.m[,c(1,24,26,25,3:4,11:22,7:10)]
	#assign the replicate
	repA = grepl('repA', pep.r$Spectrum.File)
	repB = grepl('repB', pep.r$Spectrum.File)
	pep.r$Rep = 'NA'
	pep.r[repA,23] = 'A'
	pep.r[repB,23] = 'B'
	#assign the organism
	ec = grepl('True', pep.r$Std)
	print(paste('number of e coli peptides = ',table(ec)[2],sep=''))
	hu = grepl('False', pep.r$Std)
	print(paste('number of human peptides = ',table(hu)[2],sep=''))
	pep.r$Source = 'NA'
	pep.r[ec,24] = 'ecoli'
	pep.r[hu,24] = 'human'
	#aggregate the PSMs into peptides
	print('aggregating peptides')
	pep.a = aggregate(cbind(X127N,X127C,X128N,X129C,X130N,X126,X128C,X130C,X131,X129N,Percolator.PEP)~Accession+Gene+Descriptions+Sequence+Rep+Source,data=pep.r,mean,na.action=na.pass,na.rm=TRUE)
	colnames(pep.a) = c('Accession','Gene','Descriptions','Sequence','Rep','Source','t1210a','t1246a','t1215a','t1493a','t1303a','t1811D','376D','t1842b','t1768a','t1798a','PEP')
	print(dim(pep.a))
	#replace NaN with NA
	pep.a[,7:16] = round(pep.a[,7:16],2)
	pep.a[,7:16][is.na(pep.a[,7:16])]<-NA
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
	x = as.matrix(cor(heatIn[,7:15], use='pairwise.complete.obs'))
	#make the plot labels and boundaries
	xLabels<- names(heatIn[,7:15])
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
			colsep = 1:9,
			rowsep = 1:9,
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
ovHeat_ecoli <- function(heatIn,...){
	x = as.matrix(cor(heatIn[,7:15], use='pairwise.complete.obs'))
	#make the plot labels and boundaries
	xLabels<- names(heatIn[,7:15])
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
			colsep = 1:9,
			rowsep = 1:9,
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
#t1798a looks very weird...exclude for these analyses 
#replicate B seems to follow the cell line difference trend a little better (based on origin)
processPEP <- function(psmFile,...){
	require(vsn)
	require(LSD)
	require(RColorBrewer)
	#do work
	print('Raw data frame dimensions')
	print(dim(psmFile))
	#take out the ecoli proteins
	ec = subset(psmFile, grepl('ecoli',psmFile$Source))
	hu = subset(psmFile, grepl('human', psmFile$Source))
	hu = subset(hu, rowSums(is.na(hu[,c(7:16)]))<5)
	print(paste('number of peptides remaining for human = ',nrow(hu),sep=''))
	#remove any PSMs with s/n < 5 on average across channels
	hu = subset(hu, rowMeans(hu[,c(7:16)], na.rm=TRUE) > 2)
	print('after S/N filtering')
	print(paste('number of peptides remaining for human = ',nrow(hu),sep=''))
	#plot the ecoli replication
	pdf(paste('ch_Ecoli_peptideReplicate_',ec$Rep[1],'_corHeat.pdf',sep=''))
	ovHeat_ecoli(ec)
	dev.off()
	#plot the human replication
	pdf(paste('ch_Human_peptideReplicate_',hu$Rep[1],'_corHeat.pdf',sep=''))
	ovHeat_human(hu)
	dev.off()
	#do vsn normalization
	print('doing normalization')
	#hu = (length(spikeins)+1):nrow(huec)
	exprsFile1 = as.matrix(hu[,7:16])
	#spfit = vsn2(exprsFile1[spikeins,],lts.quantile=0.95)
	#exprsFile1[hu,] = normalizeBetweenArrays(exprsFile1[hu,], method = 'quantile')
	#nkid = as.data.frame(predict(spfit, newdata = exprsFile1))
	xnorm = as.data.frame(normalizeVSN(exprsFile1))
	colnames(xnorm)<-names(hu[,7:16])
	pep.q = cbind(hu[,c(1:5,17)],xnorm)
	#output the model fit
	pdf(paste('ch_PSMnormalization_Rep',hu$Rep[1],'_vsnFit.pdf',sep=''))
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
	eset = as.matrix(x[,7:16])
	design <- cbind(SR=c(0,0,0,0,0,1,1,1,1,1),CC=c(1,1,1,1,1,0,0,0,0,0))
	fit <- lmFit(eset,design)
	cont.matrix <- makeContrasts(SRvsCC=SR-CC, levels=design)
	fit2 <- contrasts.fit(fit, cont.matrix)
	fit2 <- eBayes(fit2)
	x$logFC = fit2$coef
	x$pVal = fit2$p.value
	#output the data
	return(x[,c(1:4,17:18)])
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
proTis = merge(proA,proB,by=c('Accession','Gene','Descriptions'),sort=FALSE,all=TRUE)


##########################################################
##make a plot for peptide metrics
##########################################################
#use proA as a starter for this
#subset out peptides with greater than 50 peptides assigned
vE = proTis
vE$pepSum = ifelse(vE$pepNum.x > vE$pepNum.y, vE$pepNum.x, vE$pepNum.y)
#this text gives you the plot...choose based on the set you want to plot...e.g. sum, repA, or repB
nPeps<-data.frame(table(vE[,9]))
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
#
mp<-barplot(dn[,2])
pdf('ch_OvC_PeptideIDMetrics-RepB_ctpTissue.pdf')
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

#716