rm(list=ls())
file_dir = "~/projects/m6a/fenqi"
setwd(dir = file_dir)
library(ggplot2)
library(ggpubr)
library(stringr)
library(ggside)
library(dplyr)

tumor_list = c('ACC','BLCA','BRCA','CESC','COAD', 'CHOL', 'ESCA', 'GBM', 'HNSC', 'KIRP', 'KIRC', 'KICH', 'LGG', 'LUAD', 'LUSC', 'LIHC', 'LAML', 'OV', 'PRAD', 'PAAD', 'PCPG', 'READ', 'STAD', 'SKCM', 'THCA', 'TGCT', 'UCEC', 'UCS')

my_data = read.csv('~/projects/m6a/data/m6a.csv',header=T,check.names=F)

other_data = read.csv("~/diff/data/M.csv",,header=T,check.names=F)

temp_my_data = data.frame()
for(i in tumor_list)
{
  if(length(rownames(other_data[grep(i,other_data$CODE),]))<2)
  {next}
  temp_data = my_data[grep(i,my_data$CODE),] %>% subset(select=-c(CODE,Group))

  other_data_temp = other_data[grep(i,other_data$CODE),] %>% subset(select=-c(CODE))
  
  temp_final = merge(temp_data,other_data,by='SampleName')
  temp_my_data = rbind(temp_my_data,temp_final)
}

colnames(temp_my_data)[length(colnames(temp_my_data))]='Group'

gene = 'IGF2BP3'

exp_file <- temp_my_data[,which(colnames(temp_my_data)%in%c(gene,'CODE','Group'))]

exp_file = exp_file[which(str_extract(exp_file$CODE,'[A-Z]+')%in%tumor_list),]

dat <- data.frame(Gene=exp_file$CODE,Group = exp_file$Group,value=exp_file[,1],check.names = F)

dat[,3] = as.numeric(dat[,3])

######log
dat[,3] = log(dat[,3]+1)
######

dat <- na.omit(dat)


pdf(paste(gene,"_M.pdf",sep=''),width=length(unique(dat[,1]))+5,height = 10)
p <- ggboxplot(dat, x = "Gene", y = "value",
                    fill = "Group", palette = c('#253b6e','#d2ecf9','#1f5f8b','#1891ac'),
                    x.text.angle=60)+
     xlab("Gene")+ylab("Gene Expression(log(x+1))")+
     theme(axis.text = element_text(size = 30),axis.title=element_text(size=30))+
     stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
print(p)

# p <- ggplot(dat, aes(x = Gene, y = value, color = Group)) +
#   geom_violin(scale = "width", adjust =1, trim = TRUE) +
#   geom_boxplot(width=0.5,position=position_dodge(0.9)) +
#   labs(x="Tumor",y = "Gene Expression(log(x+1))") +
#   theme_classic()+
#   scale_color_lancet()+
#   theme(legend.key.size = unit(0.2, "inches"),axis.title = element_text(size = 30,color='black'),axis.text.x = element_text(size = 20,angle = 60,vjust = 1,hjust = 1,color='black'),axis.text.y = element_text(size = 30,color='black'))
# print(p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova"))
dev.off()

