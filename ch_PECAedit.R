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

#need to adapt PECA because I want to input log values from VSN
###############################################################################
#new PECA function
###############################################################################
PECAnorm = function(genes, accessions = NULL, probeintensities = NULL, iStdtype = FALSE, 
		iStdloc = NULL, spike = FALSE, spikeloc = NULL, samplenames1 = NULL, samplenames2 = NULL, 
		normalize = FALSE, test = NULL, type = NULL, paired = FALSE)
{	
	#divides by the internal standard column
	if (iStdtype == "divide") {
		message("Dividing by internal standard")
		flush.console()
		vnorm = apply(probeintensities,2, function(x) x / iStdloc)
	}
	
	#normalize the ecoli spikes with VSN
	if (spike == 'ecoli') {
		message("Using e-coli rows to normalize")
		flush.console()
		spikeins = grepl('EColi', spikeloc)
		vIN = as.matrix(probeintensities)
		spfit = vsn2(vIN[spikeins,],lts.quantile=1)
		nkid = predict(spfit, newdata=vIN)
		nkid.s = nkid[!spikeins,]
	}
	
	#does median normalization
	if (normalize == "median") {
		message("Performing median normalization")
		flush.console()
		xnorm <- normalizeMedianValues(as.matrix(nkid.s))
	}
	
	#do the statistics
	colnames(xnorm) <- c(samplenames1, samplenames2)
	message("Calculating low-level statistics")
	flush.console()
	if (test == "modt") {
		design <- cbind(G1 = 1, G1vsG2 = c(rep(1, length(samplenames1)),rep(0, length(samplenames2))))
		probeSLR <- as.matrix(cbind(xnorm[, samplenames1],xnorm[, samplenames2]))
		fit <- lmFit(probeSLR, design)
		fit <- eBayes(fit)
		probeSLR <- fit$coefficients[, 2]
		t <- fit$t[, 2]
		df.total <- fit$df.residual[1] + fit$df.prior
		rm(fit)
		gc()
	}
	
	#calculate p-values
	message("Aggregating statistics")
	flush.console()
	pGenes = paste(genes[!spikeins],'_',accessions[!spikeins], sep='')
	#pGenes = genes[!spikeins]
	gene.n <- tapply(t, pGenes, function(x) sum(!is.na(x)))
	print(length(gene.n))
	print(length(probeSLR))
	if (type == "median") {
		geneSLR <- tapply(probeSLR, pGenes, median, na.rm = TRUE)
		t <- tapply(t, pGenes, median, na.rm = TRUE)
	}
	gene.p <- 2 * pt(abs(t), df = df.total, lower.tail = FALSE)
	gene.p2 <- gene.p
	if (type == "median") {
		gene.p2 <- pbeta(gene.p, gene.n/2 + 0.5, gene.n - (gene.n/2 + 0.5) + 1)
	}
	gene.p.fdr <- p.adjust(gene.p2, method = "fdr")
	rm(probeintensities)
	rm(probeSLR)
	gc()
	result <- data.frame(cbind(logFC = geneSLR, t = t, score = gene.p,
					p = gene.p2, p.fdr = gene.p.fdr))
	result$Acc = row.names(result)
	row.names(result) = NULL
	message("Done")
	#output data
	return(result)	
}

PECApre = function(genes, peptides, probeintensities = NULL, iStdtype = FALSE, 
		iStdloc = NULL, spike = FALSE, spikeloc = NULL, normalize = FALSE)
{	
	#divides by the internal standard column
	if (iStdtype == "divide") {
		message("Dividing by internal standard")
		flush.console()
		vnorm = apply(probeintensities,2, function(x) x / iStdloc)
	}
	
	#normalize the ecoli spikes with VSN
	if (spike == 'ecoli') {
		message("Using e-coli rows to normalize")
		flush.console()
		spikeins = grepl('EColi', spikeloc)
		vIN = as.matrix(probeintensities)
		spfit = vsn2(vIN[spikeins,],lts.quantile=1)
		nkid = predict(spfit, newdata=vIN)
		nkid.s = nkid[!spikeins,]
	}
	
	#does median normalization
	if (normalize == "median") {
		message("Performing median normalization")
		flush.console()
		xnorm <- as.data.frame(normalizeMedianValues(as.matrix(nkid.s)))
	}
	
	xnorm$Acc = genes[!spikeins]
	xnorm$Sequence = peptides[!spikeins]
	message("Aggregating replicates")
	flush.console()
	vnorm = aggregate(cbind(a1,a2,a3,a4,a5,a6,a7,a8,a9)~Acc+Sequence,data=xnorm,median,na.action=na.pass,na.rm=TRUE)
	#output data
	return(vnorm)	
}

aIN = rbind(a1,a2,a3)
aSet = PECApre(aIN$Acc,aIN[,c(7:16)],peptides=aIN[,5],iStdtype='divide',iStdloc=aIN[,16],spike='ecoli',spikeloc=aIN[,3],normalize='median')
bIN = rbind(b1,b2,b3)
bSet = PECApre(bIN$Acc,bIN[,c(7:16)],peptides=bIN[,5],iStdtype='divide',iStdloc=bIN[,16],spike='ecoli',spikeloc=bIN[,3],normalize='median')


PECAstat = function(genes, probeintensities = NULL, samplenames1 = NULL, samplenames2 = NULL, test = NULL, type = NULL, paired = FALSE)
{	
	#do the statistics
	colnames(probeintensities) <- c(samplenames1, samplenames2)
	message("Calculating low-level statistics")
	flush.console()
	if (test == "modt") {
		design <- cbind(G1 = 1, G1vsG2 = c(rep(1, length(samplenames1)),rep(0, length(samplenames2))))
		probeSLR <- as.matrix(cbind(probeintensities[, samplenames1],probeintensities[, samplenames2]))
		fit <- lmFit(probeSLR, design)
		fit <- eBayes(fit)
		probeSLR <- fit$coefficients[, 2]
		t <- fit$t[, 2]
		df.total <- fit$df.residual[1] + fit$df.prior
		rm(fit)
		gc()
	}
	
	#calculate p-values
	message("Aggregating statistics")
	flush.console()
	#pGenes = paste(genes[!spikeins],'_',accessions[!spikeins], sep='')
	pGenes = genes
	gene.n <- tapply(t, pGenes, function(x) sum(!is.na(x)))
	print(length(gene.n))
	print(length(probeSLR))
	if (type == "median") {
		geneSLR <- tapply(probeSLR, pGenes, median, na.rm = TRUE)
		t <- tapply(t, pGenes, median, na.rm = TRUE)
	}
	gene.p <- 2 * pt(abs(t), df = df.total, lower.tail = FALSE)
	gene.p2 <- gene.p
	if (type == "median") {
		gene.p2 <- pbeta(gene.p, gene.n/2 + 0.5, gene.n - (gene.n/2 + 0.5) + 1)
	}
	gene.p.fdr <- p.adjust(gene.p2, method = "fdr")
	rm(probeintensities)
	rm(probeSLR)
	gc()
	result <- data.frame(cbind(logFC = geneSLR, t = t, score = gene.p,
					p = gene.p2, p.fdr = gene.p.fdr))
	result$Acc = row.names(result)
	row.names(result) = NULL
	message("Done")
	#output data
	return(result)	
}



#grouping of samples
group1 = c('a1.x','a2.x','a3.x','a1.y','a2.y','a3.y')
group2 = c('a4.x','a5.x','a6.x','a4.y','a5.y','a6.y')
#merging batches
pSet = merge(aSet,bSet,by=c('Acc','Sequence'),all=TRUE)
#get stats
oSet = PECAstat(pSet$Acc,pSet[,c(group1,group2)],group1,group2,'modt','median')








###test run the function
#grouping of samples
group1 = c('a1','a2','a3')
group2 = c('a4','a5','a6')
#apply the function
x = PECAnorm(sn$Acc, probeintensities = sn[,c(group1,group2)], spike = "ecoli", spikeloc = a1[,3], normalize='median', samplenames1 = group1, samplenames2 = group2, test = 'modt', type= 'median')


#check the output data
write.table(oSet,'test.txt',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)
#looks good!









