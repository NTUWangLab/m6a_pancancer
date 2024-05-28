
setwd("~/projects/Tcell-pancancer/cortri/new")

library(linkET)
library(ggplot2)
library(dplyr)
library(GSVA)
library(data.table)

data_all_tumor <- read.csv(paste("~/projects/Tcell/data/Tcell_tpm_exp.csv",sep=''),check.names = F)
riskscore_geneset <- c("ATF6B","ADA","AHNAK","AKR1C4","AHCY","BATF","B2M","CLIC1","CDK1","CDK2","CYP27A1","CALML3","CD19","DUPD1","DBI","FOSB","GPD1","GPN3","HOMER1","CXCL12","HLA-A","IFNL2","IL12B","ITM2A","IL1RN","LTBR","LIG3","MS4A3","MRPL51","NFYB","NGFR","RAN","SLC10A7","ZNF830","DCLRE1B")
tumor_list = c('ACC','BLCA','BRCA','CESC','COAD', 'CHOL', 'ESCA', 'GBM', 'HNSC', 'KIRP', 'KIRC', 'KICH', 'LGG', 'LUAD', 'LUSC', 'LIHC', 'LAML', 'OV', 'PRAD', 'PAAD', 'PCPG', 'READ', 'STAD', 'SKCM', 'THCA', 'TGCT', 'UCEC', 'UCS')

out_put_list <- c()

for (tumor_id in tumor_list)
{

  dat_one_tumor <- data_all_tumor[grep(tumor_id,data_all_tumor$CODE),]
  CIBER <- subset(dat_one_tumor,select=-c(Group,CODE))
  
  
  
  
  
  ###################ssGSVA
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
  
  out <- merge(CIBER,ssgsea,by='SampleName')
  ssgsea <- data.frame(risk = out$risk)
  out <- subset(out,select=-c(risk,SampleName))
  ###################
  
  
  
  mantel <- mantel_test(ssgsea, out,
                        spec_select = list(ssGSVA = 1:1)) %>% 
    mutate(rd = cut(r, breaks = c(-1, 0.2, 0.4, 1),
                    labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),#对相关系数进行分割，便于映射大小
           pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                    labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))#对P值进行分割，便于映射颜色
  
  
  pdf(paste(tumor_id,'.pdf',sep=''),height=10,width=10)
  p = qcorrplot(correlate(out), type = "lower", diag = FALSE) +#热图绘制
    geom_square() +#热图绘制
    geom_couple(aes(colour = pd, size = rd),data = mantel,curvature = nice_curvature()) +#aes里面是线条格式，data对应的是mantel test 计算结果，curvature控制线条曲率
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"),limits=c(-1,1)) +
    scale_size_manual(values = c(0.5, 1, 2)) +
    scale_colour_manual(values = color_pal(3)) +
    guides(size = guide_legend(title = "r",##guides()函数调整标签样式
                               override.aes = list(colour = "grey35"), 
                               order = 2),
           colour = guide_legend(title = "p", 
                                 override.aes = list(size = 3), 
                                 order = 1),
           fill = guide_colorbar(title = "R", order = 3))
  print(p)
  dev.off()
  eval(parse(text = paste0(tumor_id,'= p')))
  out_put_list <- append(out_put_list,tumor_id)
  
  write.csv(data.frame(r = mantel$r,p=mantel$p,row_1=mantel$spec,row_2=mantel$env),file=paste0(tumor_id,'_ssgsva.csv'))
  write.csv(correlate(out)[[1]],file=paste0(tumor_id,'_tri.csv'))
}

pdf("fig1.pdf",width = 50,height = 70)

eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=5,labels=out_put_list,label_size=50))',sep='')))

dev.off()
