# TODO: Add comment
# 
# Author: cshughes
###############################################################################
setwd(dir="/Users/cshughes/Documents/projects/OvC/RNAExpression") #change this to whatever directory you have stored the data in
prog<-read.table("./ch_prognosismarkers.txt", header=TRUE, sep='\t')
colnames(prog)[2] = 'Gene'

#grab the protein data...want RLE for HGS
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput") #change this to whatever directory you have stored the data in
pro = readRDS('ch_OvC_wStd_RLEset.rds')
#calculate RLE
pro$med = apply(pro[,2:19],1,function(x) median(x,na.rm=TRUE))
pro[,2:19] = apply(pro[,2:19],2,function(x) x - pro$med)
#mean RLE for HGS
pro$hgs = rowMeans(pro[,2:7],na.rm=TRUE)

#merge with the prognosis data
pprog = merge(prog,pro[,c(1,21)],by='Gene')
#make a color set
cols = brewer.pal(9,'RdBu')
#9 is blue
pCol = ifelse(grepl('poor',pprog$Gene.set),cols[1],cols[9])
gn = c('FOLR1','CRIP1','MSLN','SNCG','CRABP2','LEFTY1','GDF15','QPCT','GPC3','CTH')
cCex = ifelse(pprog$Gene %in% gn, 3,1)

pdf('ch_OvC_TMT10_wStd_Human_HGSMarkers_OvC-TCGA_prognosis.pdf')
plot(pprog$hgs,
		-log10(pprog$p.value),
		col = pCol,
		xlim = c(-1.5,1.5),
		xlab = 'RLE in HGS',
		ylab = 'p-value',
		main = 'prognosis markers in HGS',
		pch = 20,
		cex = 2
)
box(lwd=3)
text(pprog$hgs, -log10(pprog$p.value), pprog$Gene, cex = 0.25)
dev.off()









