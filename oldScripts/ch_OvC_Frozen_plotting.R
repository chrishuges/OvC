# TODO: Add comment
# 
# Author: cshughes
###############################################################################


##########################################################
##make a plot for peptide metrics
##########################################################
#subset out peptides with greater than 50 peptides assigned
vE = proTisF
vE$pepSum = ifelse(vE$pepNum.x > vE$pepNum.y, vE$pepNum.x, vE$pepNum.y)
nPeps<-data.frame(table(vE[,13]))
up<-subset(nPeps, as.numeric(as.character(nPeps$Var1))>25)
dn<-subset(nPeps, as.numeric(as.character(nPeps$Var1))<=25)
#sum the number of proteins id with more than 50 peptides
x<-sum(up$Freq)
dn[26,1]<-'26'
dn[26,2]<-x
#
mainCol<-brewer.pal(6,"RdBu")
colA<-col2rgb(mainCol[6])
xLabels<-seq(1,26,1)
xValues<-seq(1,26,1)
#
mp<-barplot(dn[,2])
pdf('ch_OvC-Frozen_PeptideIDMetrics.pdf')
barplot(dn[,2],
		col=rgb(colA[1,],colA[2,],colA[3,],100,maxColorValue=255),
		main="Distribution of Peptide Identification Numbers", 
		xlab="Number of Unique Peptides",
		ylim=c(0,max(dn[,2])+50),
		xaxt="n"
)
box(lwd=3)
axis(1,mp,xLabels,las=2,lwd=2)
dev.off()





