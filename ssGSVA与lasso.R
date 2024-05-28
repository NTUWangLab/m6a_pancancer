rm(list=ls())
file_dir = "~/projects/Tcell-pancancer"
setwd(dir = file_dir)
load('ssGSVA.Rda')
library(ggplot2)
library(ggpubr)
library(stringr)

tumor_list = c('ACC','BLCA','BRCA','CESC','COAD', 'CHOL', 'ESCA', 'GBM', 'HNSC', 'KIRP', 'KIRC', 'KICH', 'LGG', 'LUAD', 'LUSC', 'LIHC', 'LAML', 'OV', 'PRAD', 'PAAD', 'PCPG', 'READ', 'STAD', 'SKCM', 'THCA', 'TGCT', 'UCEC', 'UCS')

other_data = read.csv("~/diff/data/GRADE.csv",,header=T,check.names=F)

dat = data.frame()

for(name in tumor_list)
{
  temp_other_data = other_data[grep(name,other_data$CODE),]
  temp_other_data = subset(temp_other_data,select=-c(CODE))
  temp_final_data = final_data[grep(name,final_data$CODE),]
  out = merge(temp_other_data,temp_final_data,by='SampleName')
  dat = rbind(dat,data.frame(Gene=out$CODE,Group=out[,2],value=out$risk))
}

dat[,3] = as.numeric(dat[,3])
# dat[,2] = gsub(3,'C3',dat[,2])
# dat[,2] = gsub(2,'C2',dat[,2])
# dat[,2] = gsub(1,'C1',dat[,2])
# dat[,2] = gsub(4,'C4',dat[,2])

dat <- na.omit(dat)

pdf(paste("out.pdf",sep=''),width=length(unique(dat[,1])),height = 8)

#lasso
# p <- ggboxplot(dat, x = "Gene", y = "value",
#                fill = "Group", palette  = c('#e60042','#6f0aaa','#1924b1','#fffd00'),
#                shape = "rx",,adjust=2,notch=T,x.text.angle=60)
#cluster
# p <- ggboxplot(dat, x = "Gene", y = "value",
#                color = "Group", palette = c('#e60042','#6f0aaa','#1924b1','#fffd00'),
#                add = c("mean"),x.text.angle=60)
#fenqi
# p <- ggboxplot(dat, x = "Gene", y = "value",
#                color = "Group", palette = c('#e60042','#6f0aaa','#1924b1','#fffd00'),
#                add = c("jitter"),x.text.angle=60)

p <- p + xlab('Tumor')+ylab("Score")
p <- p + theme(axis.text = element_text(size = 30),axis.title=element_text(size=30))
print(p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova"))

dev.off()

