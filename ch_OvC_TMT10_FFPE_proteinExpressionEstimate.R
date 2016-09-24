# TODO: developing a usable protein expression estimate for individuals
# 
# Author: cshughes
###############################################################################
setwd(dir="/Users/cshughes/Documents/projects/OvC/FFPE/")
fa1psm<-read.table("./ch_08July2015_OvC-TMT10_HpH_repA1_med_PSMs.txt", header=TRUE, sep='\t')
fa1pro<-read.table("./ch_08July2015_OvC-TMT10_HpH_repA1_med_Proteins.txt", header=TRUE, sep='\t')
fa2psm<-read.table("./ch_08July2015_OvC-TMT10_HpH_repA2_med_PSMs.txt", header=TRUE, sep='\t')
fa2pro<-read.table("./ch_08July2015_OvC-TMT10_HpH_repA2_med_Proteins.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/FFPE/Routput/")

##################################################
#initial processing of the data to compile peptides and remove contaminants
##################################################
#dont split replicates yet or do any combination of them
processPSM <- function(psmFile, proteinFile, batch, samplenames, ... ){	
	#reshape the PSM to remove non-informative columns
	pepCols = c('Annotated.Sequence','Modifications','Number.of.Protein.Groups','Master.Protein.Accessions','Quan.Info','X126','X127N','X127C','X128N','X128C','X129N','X129C','X130N','X130C','X131')
	pep<-psmFile[,pepCols]
	#grab more complete annotation information from the Proteins file
	proCols = c('Accession','Description','MW.in.kDa','emPAI')
	pro<-proteinFile[,proCols]
	message('Raw number of PSMs')
	message(nrow(pep))
	#remove non-unique peptides or those without quan values
	#pep<-subset(pep, Number.of.Proteins==1)
	pep<-subset(pep, !grepl('NoQuanValues',pep$Quan.Info))
	message('Unique Peptides Only')
	message(nrow(pep))
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
	#remove contaminant proteins
	message('removing contaminants')
	pep.m = subset(pep.m, !grepl('sp',pep.m$Accession, ignore.case=FALSE))
	#remove the TMT0 PSMs
	pep.m = subset(pep.m, !grepl('TMT\\)',pep.m$Modifications))
	#filter the data columns
	message(nrow(pep.m))
	pepCols = c('Accession','Gene','Descriptions','Sequence','emPAI','X126','X127N','X127C','X128N','X128C','X129N','X129C','X130N','X130C','X131')
	pep.r = pep.m[,pepCols]
	#aggregate the PSMs into peptides
	message('aggregating peptides')
	pep.r$psms = 1
	pep.1 = aggregate(cbind(emPAI,X126,X127N,X127C,X128N,X128C,X129N,X129C,X130N,X130C,X131)~Accession+Gene+Descriptions+Sequence,data=pep.r,median,na.action=na.pass,na.rm=TRUE)
	pep.2 = aggregate(psms~Accession+Gene+Descriptions+Sequence,data=pep.r,sum,na.action=na.pass,na.rm=TRUE)
	pep.a = cbind(pep.1,pep.2[,5])
	colnames(pep.a)[6:15] = samplenames
	colnames(pep.a)[16] = 'psms'
	message(nrow(pep.a))
	#replace NaN with NA
	pep.a[,6:15] = round(pep.a[,6:15],2)
	pep.a[,6:15][is.na(pep.a[,6:15])]<-NA
	#filter based on S/N
	pep.f1 = subset(pep.a, rowMeans(pep.a[,6:15], na.rm=TRUE)>4)
	#filter based on NA
	pep.f2 = subset(pep.f1, rowSums(is.na(pep.f1[,6:15]))<9)
	#add batch label
	pep.f2$batch = batch
	#output the data
	return(pep.f2)	
}

#assign sample names
plex = c('hgsc1','ccc1','ccc2','ccc3','hgsc2','hgsc3','ccc4','ccc5','hgsc4','hgsc5')
#run the function
ovfa1.psm<-processPSM(fa1psm, fa1pro, 'a1', plex)
ovfa2.psm<-processPSM(fa2psm, fa2pro, 'a2', plex)

#output the data objects
saveRDS(ovfa1.psm,'ch_OvC_FFPE_processedPeptides_fa1_forProteinScore.rds')
saveRDS(ovfa2.psm,'ch_OvC_FFPE_processedPeptides_fa2_forProteinScore.rds')


##################################################
####Process to a single data set
#bind the two replicates together
ovBind = rbind(ovfa1.psm,ovfa2.psm)
#aggregate them into a single peptide set
ovPep1 = aggregate(cbind(emPAI,hgsc1,hgsc2,hgsc3,hgsc4,hgsc5,ccc1,ccc2,ccc3,ccc4,ccc5)~Accession+Gene+Descriptions+Sequence,data=ovBind,na.action=na.pass,FUN=mean,na.rm=TRUE)
ovPep2 = aggregate(psms~Accession+Gene+Descriptions+Sequence,data=ovBind,na.action=na.pass,FUN=sum,na.rm=TRUE)
ovPep = cbind(ovPep1,ovPep2[,5])
colnames(ovPep)[16] = 'psms'

##################################################
####calculate peptide scores
pepScore = ovPep
#multiply by psms
pepScore[,6:15] = pepScore[,6:15] * pepScore$psms
#fractional abundance of each peptide
pepScore[,6:15] = apply(pepScore[,6:15], 2, function(x) (x / sum(x, na.rm=TRUE))*100000)

##################################################
####roll into proteins
pepScore$pepNum = 1
ovPro1 = aggregate(cbind(pepNum,psms,hgsc1,hgsc2,hgsc3,hgsc4,hgsc5,ccc1,ccc2,ccc3,ccc4,ccc5)~Accession+Gene+Descriptions,data=pepScore,na.action=na.pass,FUN=sum,na.rm=TRUE)
ovPro2 = aggregate(emPAI~Accession+Gene+Descriptions,data=pepScore,na.action=na.pass,FUN=mean,na.rm=TRUE)
ovPro = merge(ovPro1,ovPro2,by=c('Accession','Gene','Descriptions'))

##################################################
####calculate protein scores
proScore = ovPro
proScore[,6:15] = proScore[,6:15] * proScore$emPAI


##################################################
####sanity checks
proScore[grepl('CRABP2',proScore$Gene),]
proScore[grepl('LEFTY1',proScore$Gene),]
proScore[grepl('MSLN',proScore$Gene),]
proScore[grepl('CTH',proScore$Gene),]
proScore[grepl('CRYAB',proScore$Gene),]

write.table(proScore,'test.txt',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)





##################################################
####ad the HPA RNA data
setwd(dir="/Users/cshughes/Documents/projects/OvC/RNAexpression/")
ens<-read.table("./HumanProteomeMap_Science_TableS1.txt", header=TRUE, sep='\t')
hpa<-read.table("./HumanProteinAtlas_FPKM_allsamples.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/FFPE/Routput/")
#change the gene colname
colnames(ens)[2] = 'Gene'
colnames(hpa)[1] = 'ensg_id'
#merge the two human atlas sets
hpa.m = merge(ens,hpa, by='ensg_id')
#remove isoforms
isoforms = grepl('\\.', hpa.m$Gene)
isoforms2 = grepl('\\-', hpa.m$Gene)
hpa.mi = hpa.m[!isoforms,]
hpa.mi = hpa.mi[!isoforms2,]
#combine with the protein data
prSet = merge(proScore,hpa.mi[,c(2,130:131)],by='Gene',all=FALSE,sort=FALSE)

prSet[grepl('CRABP2',prSet$Gene),]
prSet[grepl('CTH',prSet$Gene),]

#test against GTEx
setwd(dir="/Users/cshughes/Documents/projects/pog/humanReq/RNAseq/")
rna<-read.table("./GTEx_v6_RNAseq_Median_rpkm_byTissue.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/FFPE/Routput/")
rna.ov = rna[,c(1:2,42)]
colnames(rna.ov)[2] = 'Gene'
prgSet = merge(prSet[,c(1,4:15,17:18)],rna.ov[,2:3],by='Gene')
prgSet$meanHGS = rowMeans(prgSet[,4:8],na.rm=TRUE)
prgSet$meanCCC = rowMeans(prgSet[,9:13],na.rm=TRUE)

write.table(prgSet,'test.txt',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)
