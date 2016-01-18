# TODO: Add comment
# 
# Author: cshughes
###############################################################################
##################################################
#read in the human protein atlas data, for proteins
##################################################
setwd(dir="/Users/cshughes/Documents/projects/OvC/RNAexpression/")
hpaP<-read.table("./ch_HPAprotein_Data.csv", header=TRUE, sep=',')
setwd(dir="/Users/cshughes/Documents/projects/OvC/wStd/Routput/")
##################################################
#first lets look at just ovarian cancer
##################################################
#subset out only ovarian tissue samples
hpa.p = hpaP[grepl('ovarian',hpaP$Tumor),]
#make the scoring numeric
high = grepl('High',hpa.p$Level)
med = grepl('Medium',hpa.p$Level)
low = grepl('Low',hpa.p$Level)
hpa.p$score = '0'
#assign numbers to levels
hpa.p[high,8] = 9
hpa.p[med,8] = 6
hpa.p[low,8] = 3
#make the expression
hpa.p$exp = as.numeric(hpa.p$score) * hpa.p$Count.patients
#aggregate into a single set per gene...some of the gene IDs are messed up as dates...fix later
hpa.set = aggregate(exp~Gene.name+Total.patients,data=hpa.p,sum,na.rm=TRUE)
#get your gene data out
gn = c('FOLR1','CRIP1','MSLN','SNCG','CRABP2','LEFTY1','GDF15','QPCT','GPC3','CTH')
hpa.sub = hpa.set[hpa.set$Gene %in% gn,]
hpa.t = t(hpa.sub[,2])
colnames(hpa.t) = hpa.sub$Gene.name
hpa.s = hpa.t[,c(2,7,9,1,6,4,8,5,3)]
#plot the protein data
cols1 = c(rep(brewer.pal(9,'RdBu')[1],4),rep(brewer.pal(9,'RdBu')[9],5))
cols2 = col2rgb(c(rep(brewer.pal(9,'RdBu')[1],4),rep(brewer.pal(9,'RdBu')[9],5)))
#plot the data
pdf('ch_OvC_TMT10_wStd_Human_HGSMarkers_OvC-HPA_boxplot.pdf')
barplot(hpa.s,
		las=2,
		col = rgb(cols2[1,],cols2[2,],cols2[3,],25,maxColorValue=255),
		ylab = 'HPA Expression',
		main = 'Expression of HGS Markers in Ovarian HPA Data'
)
dev.off()




##################################################
#now lets do all cancers
##################################################
#subset out only genes of interest
hpa.p = hpaP[grepl('MSLN',hpaP$Gene.name),]
#make the scoring numeric
high = grepl('High',hpa.p$Level)
med = grepl('Medium',hpa.p$Level)
low = grepl('Low',hpa.p$Level)
hpa.p$score = '0'
#assign numbers to levels
hpa.p[high,8] = 9
hpa.p[med,8] = 6
hpa.p[low,8] = 3
#make the expression
hpa.p$exp = as.numeric(hpa.p$score) * hpa.p$Count.patients
#aggregate into a single set per gene...some of the gene IDs are messed up as dates...fix later
hpa.set = aggregate(exp~Tumor,data=hpa.p,sum,na.rm=TRUE)
#plot the protein data
cols1 = c(brewer.pal(10,'Set3'))
cols2 = col2rgb(c(brewer.pal(10,'Set3')))
#plot the data
pdf('ch_OvC_TMT10_wStd_Human_HGSMarkers-MSLN_All-HPA_boxplot.pdf')
barplot(hpa.set$exp,
		las=2,
		names = hpa.set$Tumor,
		col = rgb(cols2[1,],cols2[2,],cols2[3,],25,maxColorValue=255),
		ylab = 'HPA Expression',
		main = 'Expression of HGS Markers in All HPA Data',
		ylim = c(0,80)
)
dev.off()
