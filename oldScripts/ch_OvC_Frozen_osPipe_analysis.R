# TODO: processing the OvC Frozen samples using the open source pipeline
# 
# Author: cshughes
###############################################################################

###############################################################################
#Expression analysis
###############################################################################
setwd(dir="/Users/cshughes/Documents/projects/dataPipeline/Routput")
#read the input files
repA = dir("/Users/cshughes/Documents/projects/dataPipeline", pattern="a\\.mzXML", full.names=TRUE)
repB = dir("/Users/cshughes/Documents/projects/dataPipeline", pattern="b\\.mzXML", full.names=TRUE)
#run the quant function
quantA = lapply(repA, do_quant)
quantB = lapply(repB, do_quant)
#get the fraction order
fracA = as.list(sub(".*?\\HpH_(.*?)(\\.mzXML.*|$)", "\\1", repA))
fracB = as.list(sub(".*?\\HpH_(.*?)(\\.mzXML.*|$)", "\\1", repB))


###############################################################################
#Identification analysis
###############################################################################
idA = read.table("/Users/cshughes/Documents/projects/dataPipeline/percRep1f.txt", header=TRUE, sep='\t')
idB = read.table("/Users/cshughes/Documents/projects/dataPipeline/percRep2f.txt", header=TRUE, sep='\t')
#process replicate A...a function would be useful here...make one at some point
idA.s = subset(idA, q.value<=0.01)
idA.s = subset(idA.s, !grepl('sp\\|CONT', idA.s$proteinIds))
idA.s$Sequence = gsub('[^A-Z]','',sapply(strsplit(gsub('UNIMOD','',idA.s$peptide), '\\.'), '[', 2))
idA.s$Frac = sapply(strsplit(as.character(idA.s$PSMId), '\\_'), '[', 7)
idA.s$Scan = sapply(strsplit(as.character(idA.s$PSMId), '\\_'), '[', 9)
qli = lapply(quantA,exprs)
qann <- mapply(cbind, qli, "Frac"=fracA)
qall = as.data.frame(do.call(rbind, lapply(qann, unlist)))
qall$Scan = rownames(qall)
rownames(qall) = NULL
qAnno = qall[,c(11:12)]
qExprs = as.data.frame(do.call(rbind, lapply(qli,unlist)))
colnames(qExprs) = c('X126','X127N','X127C','X128N','X128C','X129N','X129C','X130N','X130C','X131')
rownames(qExprs) = NULL
eData = cbind(qAnno,qExprs)
idA.c = merge(idA.s,eData,by=c('Frac','Scan'))
write.table(idA.c,
		'ch_OvC_Frozen_osPipe_rawSetA.txt',
		quote=FALSE,
		sep='\t',
		col.names=TRUE,
		row.names=FALSE)

#process replicate B
idB.s = subset(idB, q.value<=0.01)
idB.s = subset(idA.s, !grepl('sp\\|CONT', idB.s$proteinIds))
idB.s$Sequence = gsub('[^A-Z]','',sapply(strsplit(gsub('UNIMOD','',idB.s$peptide), '\\.'), '[', 2))
idB.s$Frac = sapply(strsplit(as.character(idB.s$PSMId), '\\_'), '[', 7)
idB.s$Scan = sapply(strsplit(as.character(idB.s$PSMId), '\\_'), '[', 9)
qli = lapply(quantB,exprs)
qann <- mapply(cbind, qli, "Frac"=fracB)
qall = as.data.frame(do.call(rbind, lapply(qann, unlist)))
qall$Scan = rownames(qall)
rownames(qall) = NULL
qAnno = qall[,c(11:12)]
qExprs = as.data.frame(do.call(rbind, lapply(qli,unlist)))
colnames(qExprs) = c('X126','X127N','X127C','X128N','X128C','X129N','X129C','X130N','X130C','X131')
rownames(qExprs) = NULL
eData = cbind(qAnno,qExprs)
idB.c = merge(idB.s,eData,by=c('Frac','Scan'))
write.table(idB.c,
		'ch_OvC_Frozen_osPipe_rawSetB.txt',
		quote=FALSE,
		sep='\t',
		col.names=TRUE,
		row.names=FALSE)

###############################################################################
#combining identifications
###############################################################################
library(PECA)
#read in the processed data
idA.c<-read.table("./ch_OvC_Frozen_osPipe_rawSetA.txt", header=TRUE, sep='\t')
idB.c<-read.table("./ch_OvC_Frozen_osPipe_rawSetB.txt", header=TRUE, sep='\t')
#bind them together
idSet = rbind(idA.c,idB.c)
pec = aggregate(cbind(X126,X127N,X127C,X128N,X128C,X129N,X129C,X130N,X130C,X131)~proteinIds+peptide,data=idSet,mean,na.action=na.pass,na.rm=TRUE)
group2<-colnames(pec)[c(3:7)]
group1<-colnames(pec)[c(8:12)]
#the data will output as group1/group2 in the fold change
pec.q = PECA_df(pec,'proteinIds',group1,group2,normalize='median',test='modt',type='median')
pec.q$proteinIds = rownames(pec.q)
rownames(pec.q) = NULL
pec$pepNum = 1
anno = aggregate(pepNum~proteinIds,data=pec,sum,na.action=na.pass,na.rm=TRUE)
pro = merge(anno,pec.q,by='proteinIds')
write.table(pro,
		'ch_OvC_Frozen_osPipe_proteinSet.txt',
		quote=FALSE,
		sep='\t',
		col.names=TRUE,
		row.names=FALSE)

###############################################################################
#improve the annotations
###############################################################################
library(biomaRt)
#read in data
annoP = read.table("./ch_OvC_Frozen_osPipe_proteinSet.txt", header=TRUE, sep='\t')
#split up the protein accession
annoP$Accession = sapply(strsplit(as.character(annoP$proteinIds), '\\|'), '[', 2)
annoP$uKey = sapply(strsplit(as.character(annoP$proteinIds), '\\|'), '[', 3)
#get the gene names
ensembl=useMart('ENSEMBL_MART_ENSEMBL',host='www.ensembl.org')
listDatasets(ensembl) #pick your organism
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl) #pick your filters
attributes = listAttributes(ensembl) #pick your attributes
attributes[90:100,]#these are the uniprot attributes
attPro = getBM(attributes=c('uniprot_swissprot', 'uniprot_genename', 'wikigene_description'), filters = 'uniprot_swissprot', values = annoP$Accession, mart = ensembl)
#merge back with with the original set
colnames(attPro) = c('Accession','Gene','Description')
proSet = merge(attPro,annoP,by='Accession')
pro = proSet[,c(1:3,5:10)]
colnames(pro) = c('Accession','Gene','Description','NumberofPeptides','log2foldChange','pepTTest','medianPepPVal','proPVal','proQVal')
write.table(pro,
		'ch_OvC_Frozen_osPipe_proteinSet_fullAnno.txt',
		quote=FALSE,
		sep='\t',
		col.names=TRUE,
		row.names=FALSE)




