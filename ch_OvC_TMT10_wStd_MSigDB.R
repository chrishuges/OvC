# TODO: Add comment
# 
# Author: cshughes
###############################################################################
#want to use RLE data, so need to make it
#read in the data
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput") #change this to whatever directory you have stored the data in
a1 = readRDS('ch_OvC_wStd_processedPeptides_a1.rds')
a2 = readRDS('ch_OvC_wStd_processedPeptides_a2.rds')
a3 = readRDS('ch_OvC_wStd_processedPeptides_a3.rds')
b1 = readRDS('ch_OvC_wStd_processedPeptides_b1.rds')
b2 = readRDS('ch_OvC_wStd_processedPeptides_b2.rds')
b3 = readRDS('ch_OvC_wStd_processedPeptides_b3.rds')
#use this function for normalization
pepNorm = function(x,...){
	#mark the spike-in control peptides
	spikeins = grepl('EColi', x$Organism)
	#extract the expression values
	xnorm = as.matrix(x[,7:16])
	#make a noise model based on the spike ins
	spfit = vsn2(xnorm[spikeins,],lts.quantile=1)
	meanSdPlot(spfit)
	#apply it to the human peptides
	nkid = predict(spfit, newdata=xnorm)
	#extract only human expression values
	nkid.h = as.data.frame(nkid[!spikeins,])
	#normalize the values to the internal standard sample
	vnorm = apply(nkid.h[,1:9],2, function(x) x - nkid.h$aStd)
	snorm = scale(vnorm, center=TRUE,scale=TRUE)
	#recombine the replicates
	vnormFull = as.data.frame(cbind(x[!spikeins,c(1:6,17)],snorm))
	#output the data
	return(vnormFull)
}
#apply the function
a1n = pepNorm(a1)
a2n = pepNorm(a2)
a3n = pepNorm(a3)
b1n = pepNorm(b1)
b2n = pepNorm(b2)
b3n = pepNorm(b3)
#merge replicates
z1 = rbind(a1n,a2n,a3n)
z2 = rbind(b1n,b2n,b3n)
#aggregate redundant peptides
agg1 = aggregate(cbind(a1,a2,a3,a4,a5,a6,a7,a8,a9)~Gene,data=z1,median,na.action=na.pass,na.rm=TRUE)
agg2 = aggregate(cbind(a1,a2,a3,a4,a5,a6,a7,a8,a9)~Gene,data=z2,median,na.action=na.pass,na.rm=TRUE)
#merge into a gene set
proS = merge(agg1,agg2,by='Gene')
#calculate mean expression for the subtypes
proS$hgs1 = rowMeans(proS[,c(2:4,11:13)],na.rm=TRUE)
proS$ccc1 = rowMeans(proS[,c(5:7,14:16)],na.rm=TRUE)
proS$emc1 = rowMeans(proS[,c(8:10,17:19)],na.rm=TRUE)
#need t-values for comparisons
reps <- proS[,c(2:4,11:13)]
fit <- lmFit(reps)
eFit <- eBayes(fit)
proS$hgst = eFit$t
reps <- proS[,c(5:7,14:16)]
fit <- lmFit(reps)
eFit <- eBayes(fit)
proS$ccct = eFit$t
reps <- proS[,c(8:10,17:19)]
fit <- lmFit(reps)
eFit <- eBayes(fit)
proS$emct = eFit$t

###############################################################################
#query high and low expressed genes against MSigDB
###############################################################################
hgs = proS[,c(1,20,23)]
ccc = proS[,c(1,21,24)]
emc = proS[,c(1,22,25)]

#get the MSigDB data set of choice
setwd(dir="/Users/cshughes/Documents/projects/OvC/signatures/")
# Read in the data
#x <- scan("msigdb.v5.1.symbols.txt", what="", sep="\n")
x <- scan("compSigs.txt", what="", sep="\n")
# Separate elements by one or more whitepace
sigs <- strsplit(x, "[[:space:]]+")
# Extract the first vector element and set it as the list element name
names(sigs) <- sapply(sigs, `[[`, 1)
#names(y) <- sapply(y, function(x) x[[1]]) # same as above
# Remove the first vector element from each list element
sigs <- lapply(sigs, `[`, -1)
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput")


pro = hgs
#make a depository for the data
gsea = list()
#make a loop over the gene sets
for (i in 1:length(sigs)){
	gs = sigs[[i]]
	eMarker = pro$Gene %in% gs
	if (length(which(eMarker))/length(gs) > 0.5){
		ePopulation = !(pro$Gene %in% gs)
		mMarker = mean(pro[eMarker,2], na.rm=TRUE)
		mPopulation = mean(pro[ePopulation,2], na.rm=TRUE)
		pVal = geneSetTest(eMarker,pro[,3],nsim=10000)
		#pVal = NA
	}
	else{
		mMarker = NA
		mPopluation  = NA
		pVal = NA
	}
	#output the data
	gsea[[i]] = c(mMarker,mPopulation,pVal)
	names(gsea)[i] = names(sigs)[i]
}
#prepare the annotation in a better way
glist = as.data.frame(do.call(rbind,gsea))
glist$Acc = row.names(glist)
row.names(glist) = NULL
colnames(glist) = c('MarkMean','PopMean','pVal','Acc')
glist$meanDiff = glist$MarkMean - glist$PopMean
glist.s = glist[!is.na(glist$meanDiff),]
x = glist.s[order(glist.s$pVal),]


write.table(x,'ch_OvC_wStd_MSigDB_computational_HGS.txt',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)


