# TODO: Add comment
# 
# Author: cshughes
###############################################################################
library(beeswarm)
library(RColorBrewer)
###############################################################################
#read in the data
setwd(dir="/Users/cshughes/Documents/projects/OvC/IHC/")
hgs<-read.table("./CTH_HGSpanel.txt", header=FALSE, sep='\t')
ccc<-read.table("./CTH_CCCpanel.txt", header=FALSE, sep='\t')
emc<-read.table("./CTH_EMCpanel.txt", header=FALSE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/IHC/Routput")
#add a category
hgs$type = 'hgs'
ccc$type = 'ccc'
emc$type = 'emc'
#bind the data
exp = rbind(hgs,ccc,emc)
colnames(exp)[1] = 'score'


#make some colors
cols1 = c(brewer.pal(6,'Accent')[1],brewer.pal(6,'Accent')[2],brewer.pal(6,'Accent')[3])
cols2 = col2rgb(c(brewer.pal(6,'Accent')[1],brewer.pal(6,'Accent')[2],brewer.pal(6,'Accent')[3]))
#plot the data
pdf('ch_OvC_IHC-TMA_CTH_boxplot.pdf')
boxplot(exp$score~exp$type,
		las=2,
		border=cols1,
		col = rgb(cols2[1,],cols2[2,],cols2[3,],25,maxColorValue=255),
		outline = FALSE,
		boxlwd = 2,
		ylab = 'IHC score',
		main = 'Expression of CTH in OvC TMA'
)
beeswarm(exp$score~exp$type,pch=16,col='black',corral="omit",add=TRUE,cex=0.35)
dev.off()

pdf('ch_OvC_IHC-TMA_CTH_vioplot.pdf')
vioplot(hgs$V1,emc$V1,ccc$V1,
		drawRect=FALSE,
		border=cols1,
		col = rgb(cols2[1,],cols2[2,],cols2[3,],25,maxColorValue=255))
dev.off()


wilcox.test(ccc$V1,emc$V1)



