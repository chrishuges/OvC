# TODO: Add comment
# 
# Author: cshughes
###############################################################################
setwd(dir="/Users/cshughes/Documents/projects/OvC/signatures/")
# Read in the data
#x <- scan("msigdb.v5.1.symbols.txt", what="", sep="\n")
x <- scan("oncogenicSigs.txt", what="", sep="\n")
# Separate elements by one or more whitepace
sigs <- strsplit(x, "[[:space:]]+")
# Extract the first vector element and set it as the list element name
names(sigs) <- sapply(sigs, `[[`, 1)
#names(y) <- sapply(y, function(x) x[[1]]) # same as above
# Remove the first vector element from each list element
sigs <- lapply(sigs, `[`, -1)
#y <- lapply(y, function(x) x[-1]) # same as above

#the wStd data
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput") #change this to whatever directory you have stored the data in
phvc = readRDS('ch_OvC_wStd_Proteins_HGSvCCC.rds')
phve = readRDS('ch_OvC_wStd_Proteins_HGSvEMC.rds')
pcve = readRDS('ch_OvC_wStd_Proteins_CCCvEMC.rds')
rhvc = readRDS('ch_OvC_wStd_RNA_HGSvCCC.rds')
rhve = readRDS('ch_OvC_wStd_RNA_HGSvEMC.rds')
rcve = readRDS('ch_OvC_wStd_RNA_CCCvEMC.rds')

#pick a dataset
#pro = rhvc[,c(2,1,3)]
pro = phve
#make a depository for the data
gsea = list()
#make a loop over the gene sets

for (i in 1:length(sigs)){
	gs = sigs[[i]]
	eMarker = pro$Gene %in% gs
	if (length(which(eMarker))/length(gs) > 0.5){
		ePopulation = !(pro$Gene %in% gs)
		mMarker = mean(pro[eMarker,1], na.rm=TRUE)
		mPopulation = mean(pro[ePopulation,1], na.rm=TRUE)
		pVal = geneSetTest(eMarker,pro$t,nsim=10000)
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
#glist.s$fdr = p.adjust(glist.s$pVal, method = "fdr")
#glist.f = subset(glist.s, pVal<0.05)
#output the data sets
saveRDS(glist.s,'ch_OvC_Pro_MSigDB-Oncogenic_HGSvEMC.rds')




setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput") #change this to whatever directory you have stored the data in
shvcP = readRDS('ch_OvC_Pro_MSigDB-Oncogenic_HGSvCCC.rds')
shveP = readRDS('ch_OvC_Pro_MSigDB-Oncogenic_HGSvEMC.rds')
scveP = readRDS('ch_OvC_Pro_MSigDB-Oncogenic_CCCvEMC.rds')

#prepare a plot
pSig = merge(shveP,scveP, by='Acc')
cols<-rev(brewer.pal(6,"Set3"))
pCol = col2rgb(ifelse(pSig$pVal.x < 0.05 & pSig$pVal.y < 0.05, cols[3], ifelse(pSig$pVal.x < 0.05 & pSig$pVal.y > 0.05, cols[2], ifelse(pSig$pVal.x > 0.05 & pSig$pVal.y < 0.05, cols[6], NA))))
#make the plot
pdf('ch_OvC_wStd_Proteins_vEMC_GSEA-oncogenic.pdf')
plot(pSig$meanDiff.x,pSig$meanDiff.y,
		col = 'black',
		bg = rgb(pCol[1,],pCol[2,],pCol[3,],250,maxColorValue=255),
		pch = 21,
		cex = 1
)
text(pSig$meanDiff.x,pSig$meanDiff.y, pSig$Acc, cex=0.25)
abline(h=0)
abline(v=0)
box(lwd=3)
dev.off()








