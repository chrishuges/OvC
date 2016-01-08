# TODO: Add comment
# 
# Author: cshughes
###############################################################################
###############################################################################
#load the processed data
###############################################################################
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
	#apply it to the human peptides
	nkid = predict(spfit, newdata=xnorm)
	#extract only human expression values
	nkid.h = as.data.frame(nkid[!spikeins,])
	#normalize the values to the internal standard sample
	vnorm = apply(nkid.h[,1:9],2, function(x) x - nkid.h$aStd)
	#recombine the replicates
	vnormFull = as.data.frame(cbind(x[!spikeins,c(1:6,17)],vnorm))
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
#combine into a single set
aSet = rbind(a1n,a2n,a3n)
bSet = rbind(b1n,b2n,b3n)
#merge all patients
allSet = merge(aSet,bSet,by=c('Accession','Gene','Descriptions','Sequence'),all=TRUE)
allSet.s = allSet[,c(1:4,8:16,20:28)]
colnames(allSet.s) = c('Accession','Gene','Descriptions','Sequence','a1','a2','a3','a4','a5','a6','a7','a8','a9','b1','b2','b3','b4','b5','b6','b7','b8','b9')

eset = as.matrix(allSet.s[,c(8:10,17:19,11:13,20:22)])
design <- cbind(SR=c(1,1,1,1,1,1,0,0,0,0,0,0),CC=c(0,0,0,0,0,0,1,1,1,1,1,1))
fit <- lmFit(eset,design)
cont.matrix <- makeContrasts(SRvsCC=SR-CC, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
allSet.s$logFC = fit2$coef
allSet.s$pVal = fit2$p.value




#aggregate into a single peptide set
#aAgg = aggregate(cbind(a1,a2,a3,a4,a5,a6,a7,a8,a9)~Accession+Gene+Descriptions+Sequence,data=aSet,median,na.action=na.pass,na.rm=TRUE)
#aAgg = aggregate(cbind(a1,a2,a3,a4,a5,a6,a7,a8,a9)~Accession+Gene+Descriptions,data=aSet,median,na.action=na.pass,na.rm=TRUE)
aAgg = aggregate(cbind(a1,a2,a3,a4,a5,a6,a7,a8,a9,b1,b2,b3,b4,b5,b6,b7,b8,b9,logFC,pVal)~Accession+Gene+Descriptions,data=allSet.s,median,na.action=na.pass,na.rm=TRUE)

eset = as.matrix(aAgg[,4:9])
design <- cbind(SR=c(1,1,1,0,0,0),CC=c(0,0,0,1,1,1))
fit <- lmFit(eset,design)
cont.matrix <- makeContrasts(SRvsCC=SR-CC, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
aAgg$logFC = fit2$coef
aAgg$pVal = fit2$p.value

#write out the data into a format usable in excel
write.table(aAgg,'ch_OvC_wStd_CCCvEMC_proteinSet.txt',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)


cor(a1n[,8:16], use='pairwise.complete.obs')
cor(aAgg[,5:13], use='pairwise.complete.obs')
cor(aAgg[,4:21], use='pairwise.complete.obs')

