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
#want to make a growing venn diagram
##################################################
#process PSMs to peptides as a single pool
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
	pep.a = aggregate(cbind(X127N,X127C,X128N,X129C,X130N,X126,X128C,X130C,X131,X129N,Percolator.PEP)~Accession+Gene+Descriptions+Sequence+Source,data=pep.r,mean,na.action=na.pass,na.rm=TRUE)
	colnames(pep.a) = c('Accession','Gene','Descriptions','Sequence','Source','t1210a','t1246a','t1215a','t1493a','t1303a','t1811D','t376D','t1842b','t1768a','t1798a','PEP')
	print(dim(pep.a))
	#replace NaN with NA
	pep.a[,6:15] = round(pep.a[,6:15],2)
	pep.a[,6:15][is.na(pep.a[,6:15])]<-NA
	#output the data
	return(pep.a)	
}

ov.psm<-processPSM(ovpsm,ovpro)

##################################################
#do the peptide quant
##################################################
#processing function
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
	hu = subset(hu, rowSums(is.na(hu[,c(6:15)]))<8)
	print(paste('number of peptides remaining for human = ',nrow(hu),sep=''))
	#remove any PSMs with s/n < 5 on average across channels
	hu = subset(hu, rowMeans(hu[,c(6:15)], na.rm=TRUE) > 2)
	print('after S/N filtering')
	print(paste('number of peptides remaining for human = ',nrow(hu),sep=''))
	#do vsn normalization
	print('doing normalization')
	exprsFile1 = as.matrix(hu[,6:15])
	xnorm = as.data.frame(normalizeVSN(exprsFile1))
	colnames(xnorm)<-names(hu[,6:15])
	pep.q = cbind(hu[,c(1:5,16)],xnorm)
	#output the data
	return(pep.q)	
}

#process the data
ov.pep = processPEP(ov.psm)

#aggregate into proteins
pro = aggregate(cbind(t1210a,t1246a,t1215a,t1493a,t1303a,t1811D,t376D,t1842b,t1768a,t1798a)~Accession+Gene+Descriptions,data=ov.pep,mean,na.action=na.pass,na.rm=TRUE)



##################################################
#build the Venn numbers
##################################################
t1na = is.na(pro[,13])
t1 = pro[!t1na,2]
length(t1)
t2na = is.na(pro[,6])
t2 = c(t1,pro[!t2na,2])
t3na = is.na(pro[,9])
t3 = c(t2,pro[!t3na,2])
t4na = is.na(pro[,8])
t4 = c(t3,pro[!t4na,2])
t5na = is.na(pro[,5])
t5 = c(t4,pro[!t5na,2])
t6na = is.na(pro[,11])
t6 = c(t5,pro[!t6na,2])
t7na = is.na(pro[,10])
t7 = c(t6,pro[!t7na,2])


##################################################
#make an updated expression plot
##################################################

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

ov.q = pepQuan(ov.pep)

#assign a peptide counter
ov.q$pepNum = 1
#modify the colnames
colnames(ov.q) = c('Accession','Gene','Descriptions','Sequence','logFC','pepPVal','pepNum')
#aggregate into proteins
proA = aggregate(cbind(logFC,pepPVal)~Accession+Gene+Descriptions,data=ov.q,median,na.action=na.pass,na.rm=TRUE)
#aggregate into proteins
proAp = aggregate(pepNum~Accession+Gene+Descriptions,data=ov.q,sum,na.action=na.pass,na.rm=TRUE)
#need to modify the colnames again for some reason
colnames(proA) = c('Accession','Gene','Descriptions','logFC','pepPVal')
#merge
pro = merge(proA,proAp, by=c('Accession','Gene','Descriptions'))

####################
#make the plot
####################
lnCols<-brewer.pal(6,"Blues")
cols<-rev(brewer.pal(6,"RdBu"))
vE = pro
#get the list of candidate genes
setwd(dir="/Users/cshughes/Documents/projects/OvC/")
genesCC = read.table("./Schwartz_ClearCellTable_geneList.txt", header=FALSE, sep='\t')
genesS = read.table("./Schwartz_SerousTable_geneList.txt", header=FALSE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/Routput/")

geneCC = as.character(genesCC[,1])
geneS = as.character(genesS[,1])
geneTot = c(geneCC,geneS)

#make the plotting function 
pdf('ch_OvC_TMT10_expression_volcano_wMarkers.pdf')
proSD<-sd(vE[,4], na.rm=TRUE)
xCol = col2rgb(ifelse(vE[,2] %in% geneCC, cols[1], ifelse(vE[,2] %in% geneS, cols[6],'gray70')))
xCex = ifelse(vE[,2] %in% geneTot, 2, 1)
#plot with all points
plot(vE[,4],
		-log10(vE[,5]),
		col=rgb(xCol[1,],xCol[2,],xCol[3,],150,maxColorValue=255),
		pch=20,
		cex = xCex,
		ylab = '-log10(Adjusted p-value)',
		xlab = 'log2(serous/clear cell)',
		main = 'OvC type comparison',
		xlim = c(-5,5)
)
box(lwd=3)
abline(v = -proSD, col=lnCols[6],lwd=2,lty=2)
abline(v = proSD, col=lnCols[6],lwd=2,lty=2)
abline(h = -log10(0.05), col=lnCols[6],lwd=2,lty=2)
abline(h = -log10(0.01), col=lnCols[6],lwd=2,lty=2)
text(4,0.5,paste('n=',nrow(pro),sep=""),cex=1.25)
text(4,1.6,paste('p<0.05'),cex=1.25)
text(4,2.3,paste('p<0.01'),cex=1.25)
#plot with just subset points
xCol = ifelse(vE[,2] %in% geneCC, cols[1], ifelse(vE[,2] %in% geneS, cols[6],'gray70'))
gnRM = !grepl('gray', xCol)
vB = vE[gnRM,]
xCol = col2rgb(ifelse(vB[,2] %in% geneCC, cols[1], ifelse(vB[,2] %in% geneS, cols[6],'gray70')))
plot(vB[,4],
		-log10(vB[,5]),
		col=rgb(xCol[1,],xCol[2,],xCol[3,],150,maxColorValue=255),
		pch=20,
		cex = 2,
		ylab = '-log10(Adjusted p-value)',
		xlab = 'log2(serous/clear cell)',
		main = 'OvC type comparison',
		xlim = c(-5,5)
)
box(lwd=3)
abline(v = -proSD, col=lnCols[6],lwd=2,lty=2)
abline(v = proSD, col=lnCols[6],lwd=2,lty=2)
abline(h = -log10(0.05), col=lnCols[6],lwd=2,lty=2)
abline(h = -log10(0.01), col=lnCols[6],lwd=2,lty=2)
text(4,0.5,paste('n=',sum(gnRM),sep=""),cex=1.25)
text(4,1.6,paste('p<0.05'),cex=1.25)
text(4,2.3,paste('p<0.01'),cex=1.25)
dev.off()








