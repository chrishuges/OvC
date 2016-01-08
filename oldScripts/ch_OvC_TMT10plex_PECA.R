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
library(PECA)
###############################################################################
#grab the data
setwd(dir="/Users/cshughes/Documents/projects/OvC/FFPESections/")
ovpsm<-read.table("./ch_09July2015_OvC-TMT10_HpH_repAB_PSMs.txt", header=TRUE, sep='\t')
ovpep<-read.table("./ch_09July2015_OvC-TMT10_HpH_repAB_PeptideGroups.txt", header=TRUE, sep='\t')
ovpro<-read.table("./ch_09July2015_OvC-TMT10_HpH_repAB_Proteins.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/FFPESections/Routput/")

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
	pep.m = subset(pep.m, !grepl('sp',pep.m$Accession, ignore.case=FALSE))
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
	colnames(pep.a) = c('Accession','Gene','Descriptions','Sequence','Rep','Source','t1210a','t1246a','t1215a','t1493a','t1303a','t1811D','t376D','t1842b','t1768a','t1798a','PEP')
	print(dim(pep.a))
	#replace NaN with NA
	pep.a[,7:16] = round(pep.a[,7:16],2)
	pep.a[,7:16][is.na(pep.a[,7:16])]<-NA
	#filter based on NA
	pep.f1 = subset(pep.a, rowMeans(pep.a[,7:16], na.rm=TRUE)>2)
	#filter based on S/N
	pep.f2 = subset(pep.f1, rowSums(is.na(pep.f1[,7:16]))<7)
	#output the data
	return(pep.f2)	
}

ov.psm<-processPSM(ovpsm,ovpro)

##################################################
#do the peptide quant
##################################################
ovHeat <- function(heatIn,...){
	x = as.matrix(cor(heatIn, use='pairwise.complete.obs'))
	#make the plot labels and boundaries
	xLabels<- colnames(heatIn)
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

#split by replicates
repA = subset(ov.psm, grepl('A', ov.psm$Rep))
repB = subset(ov.psm, grepl('B', ov.psm$Rep))
#remove ecoli peps
repA = subset(repA, !grepl('ecoli', repA$Source))
repB = subset(repB, !grepl('ecoli', repB$Source))
#combine into a single data frame
pep = merge(repA,repB,by=c('Accession','Gene','Descriptions','Sequence'),all=TRUE)
#comBat processing with imputation
exprs = as.matrix(pep[,c(7:16,20:29)])
xnorm = impute.knn(exprs, k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
samps = colnames(exprs)
batch = rep(1:2, each=10)
pheno = data.frame(samps,batch)
modcombat = model.matrix(~1, data=pheno)
edata = ComBat(dat=xnorm$data, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
vsndata = justvsn(edata)
#do a test plot to see how the data is transformed
pdf('ch_OvC_FFPE_Batches_Heat.pdf')
ovHeat(exprs)
ovHeat(xnorm$data)
ovHeat(edata)
ovHeat(vsndata)
dev.off()
#looks ok I think
#do the PECA analysis
pec = cbind(pep[,1:4],edata)
group1<-colnames(pec)[c(5:9,15:19)]
group2<-colnames(pec)[c(10:14,20:24)]
pec.q = PECA_df(pec,'Gene',group1,group2,normalize='median',test='modt',type='median')
pec.q$Gene = rownames(pec.q)
rownames(pec.q) = NULL
pec$pepNum = 1
anno = aggregate(pepNum~Accession+Gene+Descriptions,data=pec,sum,na.action=na.pass,na.rm=TRUE)
pro = merge(anno,pec.q,by='Gene')
#write out the data
write.table(pro,
		'ch_OvC_wImp-wCombat-wPECA_proteinSet.txt',
		quote=FALSE,
		sep='\t',
		col.names=TRUE,
		row.names=FALSE)

#do a VSN analysis
#split by replicates
repA = subset(ov.psm, grepl('A', ov.psm$Rep))
repB = subset(ov.psm, grepl('B', ov.psm$Rep))
#remove ecoli peps
repA = subset(repA, !grepl('ecoli', repA$Source))
repB = subset(repB, !grepl('ecoli', repB$Source))
#combine into a single data frame
pep = merge(repA,repB,by=c('Accession','Gene','Descriptions','Sequence'),all=TRUE)
exprs = as.matrix(pep[,c(7:16,20:29)])
xnorm = as.matrix(normalizeVSN(exprs))
ximp = impute.knn(xnorm, k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
samps = colnames(exprs)
batch = rep(1:2, each=10)
pheno = data.frame(samps,batch)
modcombat = model.matrix(~1, data=pheno)
edata = ComBat(dat=exprs, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
#do a test plot to see how the data is transformed
pdf('ch_OvC_FFPE_Batches2_Heat.pdf')
ovHeat(exprs)
ovHeat(xnorm)
ovHeat(edata)
dev.off()
#compile to proteins
ov.exprs = cbind(pep.s[,c(1:4)], xnorm)
#quant function
pepQuan <- function(x,...){
	eset = as.matrix(x[,c(5:9,15:19,10:14,20:24)])
	design <- cbind(SR=c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1),CC=c(1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0))
	fit <- lmFit(eset,design)
	cont.matrix <- makeContrasts(SRvsCC=SR-CC, levels=design)
	fit2 <- contrasts.fit(fit, cont.matrix)
	fit2 <- eBayes(fit2)
	x$logFC = fit2$coef
	x$pVal = fit2$p.value
	x.pro = x[,c(1:4,25:26)]
	x.pro$pepNum = 1
	colnames(x.pro)[5:6] = c('logFC','pepPVal')
	proData = aggregate(cbind(logFC,pepPVal)~Accession+Gene+Descriptions,data=x.pro,median,na.action=na.pass,na.rm=TRUE)
	proNum = aggregate(pepNum~Accession+Gene+Descriptions,data=x.pro,sum,na.action=na.pass,na.rm=TRUE)
	pro = merge(proNum,proData, by=c('Accession','Gene','Descriptions'))
	colnames(pro)[5:6] = c('logFC','pepPVal')
	#pro.s = subset(pro, pepNum>1)
	#output the data
	return(pro)
}

ov.q = pepQuan(ov.exprs)


#do the PECA analysis
pec = subset(ov.psm, !grepl('ecoli', ov.psm$Source))
group2<-colnames(pec)[c(7:11)]
group1<-colnames(pec)[c(12:16)]
#the data will output as group1/group2 in the fold change
pec.q = PECA_df(pec,'Gene',group1,group2,normalize='median',test='modt',type='median')
pec.q$Gene = rownames(pec.q)
rownames(pec.q) = NULL
pec$pepNum = 1
anno = aggregate(pepNum~Accession+Gene+Descriptions,data=pec,sum,na.action=na.pass,na.rm=TRUE)
pro = merge(anno,pec.q,by='Gene')

#write out the data
write.table(pro,
		'ch_test.txt',
		quote=FALSE,
		sep='\t',
		col.names=TRUE,
		row.names=FALSE)


x = subset(ov.psm, grepl('P04637', ov.psm$Accession))
ximp = impute.knn(as.matrix(pep[,c(7:16,20:29)]), k = 3, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
pep2 = cbind(pep[,1:6],ximp$data)
y = subset(pep2, grepl('P04637', pep2$Accession))




#split by replicates
repA = subset(ov.psm, grepl('A', ov.psm$Rep))
repB = subset(ov.psm, grepl('B', ov.psm$Rep))
#remove ecoli peps
repA = subset(repA, !grepl('ecoli', repA$Source))
repB = subset(repB, !grepl('ecoli', repB$Source))
#normalize
repA[,c(7:16)] = normalizeVSN(repA[,c(7:16)])
repB[,c(7:16)] = normalizeVSN(repB[,c(7:16)])
#proteins
pA = aggregate(cbind(t1210a,t1246a,t1215a,t1493a,t1303a,t1811D,t376D,t1842b,t1768a,t1798a)~Accession+Gene+Descriptions,data=repA,median,na.action=na.pass,na.rm=TRUE)
pB = aggregate(cbind(t1210a,t1246a,t1215a,t1493a,t1303a,t1811D,t376D,t1842b,t1768a,t1798a)~Accession+Gene+Descriptions,data=repB,median,na.action=na.pass,na.rm=TRUE)
pAB = merge(pA,pB,by=c('Accession','Gene','Descriptions'),all=TRUE)
samps = colnames(pAB[,c(4:23)])
batch = rep(1:2, each=10)
pheno = data.frame(samps,batch)
modcombat = model.matrix(~1, data=pheno)
edata = ComBat(dat=pAB[,c(4:23)], batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

		
		
		
		
ovHeat <- function(heatIn,...){
	x = as.matrix(cor(heatIn[,c(4:23)], use='pairwise.complete.obs'))
	#make the plot labels and boundaries
	xLabels<- colnames(heatIn[,c(4:23)])
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
pdf('ch_test.pdf')
ovHeat(pAB)
dev.off()

##################################################
#do the peptide quant
##################################################
#split by replicates
repA = subset(ov.psm, grepl('A', ov.psm$Rep))
repB = subset(ov.psm, grepl('B', ov.psm$Rep))
#remove ecoli peps
repA = subset(repA, !grepl('ecoli', repA$Source))
repB = subset(repB, !grepl('ecoli', repB$Source))
#do the PECA analysis
exprsA = repA[,c(1:3,4,7:16)]
group1<-colnames(exprsA)[10:14]
group2<-colnames(exprsA)[5:9]
xnorm = PECA_df(exprsA,'Gene',group1,group2,normalize='median',test='modt',type='median')
xnorm$Gene = rownames(xnorm)
rownames(xnorm) = NULL
exprsA$pepNum = 1
anno = aggregate(pepNum~Accession+Gene+Descriptions,data=exprsA,sum,na.action=na.pass,na.rm=TRUE)
proA = merge(anno,xnorm,by='Gene')
#replicate B
exprsB = repB[,c(1:3,4,7:16)]
group1<-colnames(exprsB)[10:14]
group2<-colnames(exprsB)[5:9]
xnorm = PECA_df(exprsB,'Gene',group1,group2,normalize='median',test='modt',type='median')
xnorm$Gene = rownames(xnorm)
rownames(xnorm) = NULL
exprsB$pepNum = 1
anno = aggregate(pepNum~Accession+Gene+Descriptions,data=exprsB,sum,na.action=na.pass,na.rm=TRUE)
proB = merge(anno,xnorm,by='Gene')



write.table(proB,
		'ch_OvC_TMT10plex_proteinSet_B.txt',
		quote=FALSE,
		sep='\t',
		col.names=TRUE,
		row.names=FALSE)


##########################################################
##make volcano plots of expression variance
##########################################################
#make some colors
lnCols<-brewer.pal(6,"Blues")
cols<-rev(brewer.pal(6,"RdBu"))
vE = proA
#get the list of candidate genes
setwd(dir="/Users/cshughes/Documents/projects/OvC/")
genesCC = read.table("./Schwartz_ClearCellTable_geneList.txt", header=FALSE, sep='\t')
genesS = read.table("./Schwartz_SerousTable_geneList.txt", header=FALSE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/Routput/")

geneCC = as.character(genesCC[,1])
geneS = as.character(genesS[,1])
geneTot = c(geneCC,geneS)

#make the plotting function 
pdf('ch_OvC_test.pdf')
proSD<-sd(vE[,5], na.rm=TRUE)
xCol = col2rgb(ifelse(vE[,1] %in% geneCC, cols[1], ifelse(vE[,1] %in% geneS, cols[6],'gray30')))
xCex = ifelse(vE[,1] %in% geneTot, 2, 1)
#plot with all points
plot(vE[,5],
		-log10(vE[,8]),
		col=rgb(xCol[1,],xCol[2,],xCol[3,],95,maxColorValue=255),
		pch=20,
		cex = xCex,
		ylab = '-log10(Adjusted p-value)',
		xlab = 'log2(serous/clear cell)',
		main = 'OvC type comparison',
		xlim = c(-5,5),
		ylim = c(0,10)
)
box(lwd=3)
abline(v = -proSD, col=lnCols[6],lwd=2,lty=2)
abline(v = proSD, col=lnCols[6],lwd=2,lty=2)
abline(h = -log10(0.05), col=lnCols[6],lwd=2,lty=2)
abline(h = -log10(0.01), col=lnCols[6],lwd=2,lty=2)
text(4,0.5,paste('n=',nrow(proA),sep=""),cex=1.25)
text(4,1.6,paste('p<0.05'),cex=1.25)
text(4,2.3,paste('p<0.01'),cex=1.25)
dev.off()
