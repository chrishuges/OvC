# TODO: Add comment
# 
# Author: cshughes
###############################################################################
setwd(dir="/Users/cshughes/Documents/projects/OvC/signatures/")
# Read in the data
#x <- scan("msigdb.v5.1.symbols.txt", what="", sep="\n")
x <- scan("hallmarkSigs.txt", what="", sep="\n")
# Separate elements by one or more whitepace
sigs <- strsplit(x, "[[:space:]]+")
# Extract the first vector element and set it as the list element name
names(sigs) <- sapply(sigs, `[[`, 1)
#names(y) <- sapply(y, function(x) x[[1]]) # same as above
# Remove the first vector element from each list element
sigs <- lapply(sigs, `[`, -1)
#y <- lapply(y, function(x) x[-1]) # same as above

#get the protein data
setwd(dir="/Users/cshughes/Documents/projects/OvC/FFPE/Routput") #change this to whatever directory you have stored the data in
pro = readRDS('ch_OvC_FFPE_proteinSet.rds')
colnames(pro)[5] = 'PROexp'
#the wStd data
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput") #change this to whatever directory you have stored the data in
phvc = readRDS('ch_OvC_wStd_Proteins_HGSvCCC.rds')
phve = readRDS('ch_OvC_wStd_Proteins_HGSvEMC.rds')
pcve = readRDS('ch_OvC_wStd_Proteins_CCCvEMC.rds')

#pick a dataset
pro = phvc
#make a depository for the data
gsea = list()
#make a loop over the gene sets
for (i in 1:length(sigs)){
	gs = sigs[[i]]
	eMarker = pro[pro$Gene %in% gs,c(1,5)]
	if (nrow(eMarker)/length(gs) > 0.6){
		ePopulation = pro[!(pro$Gene %in% gs),c(1,5)]
		mMarker = mean(eMarker$PROexp, na.rm=TRUE)
		mPopulation = mean(ePopulation$PROexp, na.rm=TRUE)
		pVal = wilcox.test(ePopulation$PROexp,eMarker$PROexp)$p.value
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
glist.s = glist[!is.na(glist$pVal),]
glist.s$fdr = p.adjust(glist.s$pVal, method = "fdr")
glist.f = subset(glist.s, fdr<0.05)
#output the data sets
saveRDS(glist.f,'ch_OvC_wStd_MSigDB-Curated_HGSvCCC.rds')







setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput") #change this to whatever directory you have stored the data in
shvc = readRDS('ch_OvC_wStd_MSigDB-Curated_HGSvCCC.rds')
shve = readRDS('ch_OvC_wStd_MSigDB-Curated_HGSvEMC.rds')

shvc2 = readRDS('ch_OvC_FFPE_MSigDB_HGSvCCC.rds')

pSig = merge(shvc,shve, by='Acc')
#xnorm = apply(pSig[,c(5,10)],2, function(x) x - median(x,na.rm=TRUE))
#row.names(xnorm) = pSig$Acc
pdf('ch_test.pdf')
plot(pSig[,5],pSig[,10])
text(pSig[,5],pSig[,10], pSig$Acc, cex=0.25)
abline(h=0)
abline(v=0)
dev.off()








