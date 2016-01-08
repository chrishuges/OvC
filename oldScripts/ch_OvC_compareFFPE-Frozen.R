# TODO: Add comment
# 
# Author: cshughes
###############################################################################
#proTis is the FFPE set, proTisF is the frozen set
#check the correlation of the expression measurements with both methods
x = merge(proTis,proTisF,by=c('Gene'))
cor(x[,10],x[,21],use='pairwise.complete.obs',method='spearman')
pdf('ch_OvC_FFPEvFroze_ExpressionScatter.pdf')
reg = lm(x[,7]~x[,4])
heatscatter(x[,4],x[,7],main='FFPE Biological replicates',xlab='repA',ylab='repB',xlim=c(-4,4),ylim=c(-4,4))
text(2,-3,paste('n=',nrow(x) - length(which(is.na(x[,4]))),sep=''))
abline(reg,lty=2,lwd=3)
reg = lm(x[,15]~x[,18])
heatscatter(x[,15],x[,18],main='Frozen Biological replicates',xlab='repA',ylab='repB',xlim=c(-4,4),ylim=c(-4,4))
text(2,-3,paste('n=',nrow(x) - length(which(is.na(x[,15]))),sep=''))
abline(reg,lty=2,lwd=3)
reg = lm(x[,10]~x[,21])
heatscatter(x[,10],x[,21],main='FFPE vs Frozen (combined replicates)',xlab='Frozen',ylab='FFPE',xlim=c(-4,4),ylim=c(-4,4))
text(2,-3,paste('n=',nrow(x) - length(which(is.na(x[,10]))),sep=''))
abline(reg,lty=2,lwd=3)
dev.off()





