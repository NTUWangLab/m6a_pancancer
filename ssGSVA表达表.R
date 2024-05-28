

library(data.table)
library(dplyr)
library(GSVA)

data_all_tumor <- read.csv(paste("~/projects/m6a/data/m6a.csv",sep=''),check.names = F)
riskscore_geneset <- c("ATF6B","ADA","AHNAK","AKR1C4","AHCY","BATF","B2M","CLIC1","CDK1","CDK2","CYP27A1","CALML3","CD19","DUPD1","DBI","FOSB","GPD1","GPN3","HOMER1","CXCL12","HLA-A","IFNL2","IL12B","ITM2A","IL1RN","LTBR","LIG3","MS4A3","MRPL51","NFYB","NGFR","RAN","SLC10A7","ZNF830","DCLRE1B")
tumor_list = c('ACC','BLCA','BRCA','CESC','COAD', 'CHOL', 'ESCA', 'GBM', 'HNSC', 'KIRP', 'KIRC', 'KICH', 'LGG', 'LUAD', 'LUSC', 'LIHC', 'LAML', 'OV', 'PRAD', 'PAAD', 'PCPG', 'READ', 'STAD', 'SKCM', 'THCA', 'TGCT', 'UCEC', 'UCS')

final_data = data.frame()

for (tumor_id in tumor_list)
{
  
  dat_one_tumor <- data_all_tumor[grep(tumor_id,data_all_tumor$CODE),]
  CIBER <- subset(dat_one_tumor,select=-c(Group,CODE))
  
  dat <- fread(paste('~/TCGA_DATA/counts/',tumor_id,'.txt',sep=''), sep = "\t",header = T,stringsAsFactors = F,check.names = F,na.strings="NA",data.table = F)
  dat <- dat[!duplicated(dat[,c(1)]),]
  rownames(dat) <- dat[,c(1)]
  dat <- dat[,-c(1)]
  dat <- na.omit(dat)
  dat <- dat[,which(colnames(dat)%in%CIBER$SampleName)]
  
  
  ssample <- list('risk'=riskscore_geneset)
  ssgsea <- gsva(as.matrix(dat),ssample, method='ssgsea', kcdf='Poisson',abs.ranking=TRUE)
  
  ssgsea <-as.data.frame(t(ssgsea))
  ssgsea$SampleName <- rownames(ssgsea)
  ssgsea$CODE = rep(tumor_id,length(rownames(ssgsea)))
  final_data = rbind(final_data,ssgsea)
}
