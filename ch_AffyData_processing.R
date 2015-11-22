# TODO: Add comment
# 
# Author: cshughes
###############################################################################
#grab the data
setwd(dir="/Users/cshughes/Documents/projects/OvC/RNAexpression/")
affy<-read.table("./combined_occc_enoc_u133_gene_expression_data_unique_mapped_gene_combat_processed_mapped_sampleid.txt", header=TRUE, sep='\t')

#grab the annotation
setwd(dir="/Users/cshughes/Documents/projects/OvC/annotation/")
anno<-read.table("./ch_OvC_annotation.txt", header=TRUE, sep='\t')
setwd(dir="/Users/cshughes/Documents/projects/OvC/RNAExpression/Routput/")

#subset the RNA to retain only patients with proteomics data
affyS = affy[,colnames(affy) %in% anno[,2]]

##################################################
#look at some basic data QC
##################################################
#The data has been log2 transformed (using RMA of AFFY R package) and take
#the difference of the data to the mean. We then did combat processing to
#correct the batch effect on the data.

#the correlation between any of the samples looks terrible...
#should it be so poor?
pdf('ch_OvC_RNA-Affy_patientScatter.pdf')
heatscatter(affyS[,3],
		affyS[,2],
		main = paste('Pearson R = ',cor(affyS[,3],affyS[,2],use='pairwise.complete.obs',method='pearson'),sep='')
)
dev.off()
#looks terrible...is this normal?











