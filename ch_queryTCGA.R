# TODO: Add comment
# 
# Author: cshughes
###############################################################################
mycgds = CGDS("http://www.cbioportal.org/public-portal/")
test(mycgds)
getCancerStudies(mycgds)
mycancerstudy = getCancerStudies(mycgds)[2,1]
mycaselist = getCaseLists(mycgds,mycancerstudy)[1,1]
mygeneticprofile = getGeneticProfiles(mycgds,mycancerstudy)[4,1]
getProfileData(mycgds,c('BRCA1','BRCA2'),mygeneticprofile,mycaselist)
myclinicaldata = getClinicalData(mycgds,mycaselist)






df = getProfileData(mycgds, "WT1", c("gbm_tcga_gistic","gbm_tcga_mrna"), "gbm_tcga_all")


##########################################################
##pull out the ovarian cancer data from TCGA
##########################################################
mycgds = CGDS("http://www.cbioportal.org/public-portal/")
getCaseLists(mycgds,'ov_tcga')[,c(1:2)]
getGeneticProfiles(mycgds,'ov_tcga_pub')[,c(1:2)]

mycaselist = getCaseLists(mycgds,'ov_tcga_pub')[17,1]
mygeneticprofile = getGeneticProfiles(mycgds,'ov_tcga_pub')[9,1]


df = getProfileData(mycgds, "WT1",mygeneticprofile,mycaselist)



##########################################################
##match with mutation status
##########################################################
mycgds = CGDS("http://www.cbioportal.org/public-portal/")
mycaselist = getCaseLists(mycgds,'ov_tcga_pub')[17,1]
mygeneticprofile = getGeneticProfiles(mycgds,'ov_tcga_pub')[5,1]
mut = getProfileData(mycgds, c("TP53",'ARID1A','PIK3CA','PTEN','KRAS'),mygeneticprofile,mycaselist)
mut[!is.na(mut$PIK3CA),]
#table(is.na(mut$ARID1A))
#table(is.na(mut$PIK3CA))
#table(is.na(mut$PTEN))
#table(is.na(mut$KRAS))
mygeneticprofile = getGeneticProfiles(mycgds,'ov_tcga_pub')[9,1]
rna = getProfileData(mycgds, c('TP53','ARID1A','PIK3CA','PTEN','KRAS','NID2','LEFTY1','CRABP2','WT1','GDF15','KLHL14','VCAN'),mygeneticprofile,mycaselist)
rna2 = getProfileData(mycgds, c('TP53','NID2','VCAN','SPP1','CRABP2','WT1','KLHL14'),mygeneticprofile,mycaselist)
#get all the patients out who carry ARID1A mutations
ar = rna2[!is.na(mut$ARID1A),][2,]
pt = rna[!is.na(mut$PTEN),]
pk = rna[!is.na(mut$PIK3CA),]
gnAll = rbind(ar,pt,pk)

#very interesting...definitely seem to be able to re-classify the serous cancers and pull out clear cell cases
pdf('ch_OvC_MarkerClustering_TCGA.pdf')
mybreaks = seq(-2,2,by=0.05) 
#make the correlation heatmap
heatmap.2(
		as.matrix(rna2),
		col= colorRampPalette(brewer.pal(9,"RdBu"))(length(mybreaks)-1),
		symkey=TRUE,
		Rowv=TRUE,
		Colv=TRUE,
		#hclust=hclustfunc,
		#distfun=distfunc,
		#na.color='black',
		na.rm=TRUE,
		#symm=FALSE,
		dendrogram="both",
		breaks=mybreaks,
		#cellnote = round(ovCor,2),
		labRow = '',
		#labCol = '',
		#notecol = 'black',
		#notecex = 1.5,
		#colsep = 1:70,
		#rowsep = 1:70,
		sepwidth = c(0.03, 0.03),
		sepcolor = 'white',
		#ColSideColors=ColSideColors,
		## labels
		main='TCGA',
		#xlab='Time Point',
		#ylab='Protein',
		## color key
		key = TRUE,
		keysize = 1,
		density.info = "none",
		scale = "none",
		trace = "none",
		mar=c(8,8),
		cexRow=0.75,
		cexCol=0.5
)
dev.off()

ar = ar[,c(4,3,6,5,1,2,7)]
pdf('ch_OvC_ClearCellinTCGA.pdf')
barplot(as.matrix(ar),las=2)
dev.off()



##########################################################
##cross reference with the tissue data
##########################################################
#get the tissue data
tis = pro[,c(1:5,7:8)]
tis$meanFC = rowMeans(tis[,c(4,6)], na.rm=TRUE)
tis$meanP = rowMeans(tis[,c(5,7)], na.rm=TRUE)
proSD = sd(tis$meanFC, na.rm=TRUE)
#make a gene set for classifying cancers
ser = tis[tis$meanFC> 1.5 & -log10(tis$meanP) > -log10(0.01),]
cc = tis[tis$meanFC< -1.5 & -log10(tis$meanP) > -log10(0.01),]
#make it into a usable list
ser$class = 'serous'
cc$class = 'clear'
tum = rbind(ser[,c(2,10)],cc[,c(2,10)])
#get the data from TCGA
df = getProfileData(mycgds, tum[,1], mygeneticprofile, mycaselist)
#reshape the data
exprs = data.frame(t(df))
exprs = cbind(Row.Names = rownames(exprs), exprs)
rownames(exprs) <- NULL
colnames(exprs)[1] = 'Gene'
#recombined with the class marker
tcga = merge(tum,exprs,by='Gene')
tcga = tcga[order(tcga$class),]
##########################################################
##plot a boxplot of expression for the genes across patients
##########################################################

pdf('ch_test.pdf')
boxplot(df,las=2,cex=0.5)
dev.off()

df = subset(df, WT1>1)
df$NHLRC3 = NULL
#df = data.frame(t(df[,c(1:40,42:70)]))
df = data.frame(t(df))
df = cbind(Row.Names = rownames(df), df)
rownames(df) <- NULL
colnames(df)[1] = 'Gene'
exprs = merge(tum,df,by='Gene',sort=FALSE)
pdf('ch_test.pdf')
#make the plot labels and boundaries
xLabels<- 'WT1'
mybreaks = seq(-1,1,by=0.05) 
#make the correlation heatmap
heatmap.2(
		as.matrix(exprs[,3:ncol(exprs)]),
		col= colorRampPalette(brewer.pal(9,"RdBu"))(length(mybreaks)-1),
		symkey=FALSE,
		Rowv=TRUE,
		Colv=TRUE,
		#hclust=hclustfunc,
		#distfun=distfunc,
		#na.color='black',
		na.rm=TRUE,
		#symm=FALSE,
		dendrogram="row",
		breaks=mybreaks,
		#cellnote = round(ovCor,2),
		labRow = paste(exprs[,1],exprs[,2]),
		labCol = '',
		notecol = 'black',
		notecex = 1.5,
		#colsep = 1:70,
		#rowsep = 1:70,
		sepwidth = c(0.03, 0.03),
		sepcolor = 'white',
		#ColSideColors=ColSideColors,
		## labels
		main='TCGA',
		#xlab='Time Point',
		#ylab='Protein',
		## color key
		key = TRUE,
		keysize = 1,
		density.info = "none",
		scale = "none",
		trace = "none",
		mar=c(8,8),
		cexRow=0.75,
		cexCol=0.5
)
dev.off()

#the TCGA data is all over the map...need to revisit this....



##########################################################
##use the consolidated marker set
##########################################################
mycgds = CGDS("http://www.cbioportal.org/public-portal/")
getCaseLists(mycgds,'ov_tcga')[,c(1:2)]
getGeneticProfiles(mycgds,'ov_tcga_pub')[,c(1:2)]
#extract our case of interest
mycaselist = getCaseLists(mycgds,'ov_tcga_pub')[17,1]
mygeneticprofile = getGeneticProfiles(mycgds,'ov_tcga_pub')[9,1]
#use the consolidated marker set from ch_markerConsolidation.R
##########################################################
##get the marker list out
##########################################################
ser = subset(tcg.q, -log10(pAdj) > -log10(0.05) & meanFC > 0)
ser$class = 'serous'
cc = subset(tcg.q, -log10(pAdj) > -log10(0.05) & meanFC < 0)
cc$class = 'clear'
ser.cc = rbind(ser,cc)
ttcg = c('CXCL11','CXCL10','CXCR3','HMGA2','SOX11','MUC1','MUC16','MCM2','PCNA','SLPI','FAP','ANGPTL2','ANGPTL1')
#query TCGA
#mExprs = getProfileData(mycgds, ser.cc[,1], mygeneticprofile, mycaselist)
mExprs = getProfileData(mycgds, ttcg, mygeneticprofile, mycaselist)
#rework the tcga data to get class back in 
dExprs = data.frame(t(mExprs))
dExprs = cbind(Row.Names = rownames(dExprs), dExprs)
rownames(dExprs) <- NULL
colnames(dExprs)[1] = 'Gene'
#cExprs = merge(ser.cc[,c(1,7)],dExprs,by='Gene',sort=FALSE)
cExprs = dExprs
#make a heatmap of the expression patterns
pdf('ch_OvCall_classMarkers_tcga.pdf')
#make the plot labels and boundaries
mybreaks = seq(-1,1,by=0.05) 
#make the correlation heatmap
heatmap.2(
		as.matrix(cExprs[,2:ncol(cExprs)]),
		col= colorRampPalette(brewer.pal(9,"RdBu"))(length(mybreaks)-1),
		symkey=FALSE,
		Rowv=FALSE,
		Colv=TRUE,
		#hclust=hclustfunc,
		#distfun=distfunc,
		#na.color='black',
		na.rm=TRUE,
		#symm=FALSE,
		dendrogram="column",
		breaks=mybreaks,
		#cellnote = round(ovCor,2),
		#labRow = paste(cExprs[,1],cExprs[,2]),
		labRow = cExprs[,1],
		labCol = '',
		notecol = 'black',
		notecex = 1.5,
		#colsep = 1:70,
		#rowsep = 1:70,
		sepwidth = c(0.03, 0.03),
		sepcolor = 'white',
		#ColSideColors=ColSideColors,
		## labels
		main='TCGA',
		#xlab='Time Point',
		#ylab='Protein',
		## color key
		key = TRUE,
		keysize = 1,
		density.info = "none",
		scale = "none",
		trace = "none",
		mar=c(8,8),
		cexRow=0.75,
		cexCol=0.5
)
dev.off()















