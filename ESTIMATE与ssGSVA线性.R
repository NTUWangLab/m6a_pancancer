setwd("~/projects/Tcell-pancancer")
rm(list = ls())
load('ssGSVA.Rda')
library(dplyr)
library(patchwork)
library(ggplotify)
library(psych)
library(stringr)


tumor_list = c('ACC','BLCA','BRCA','CESC','COAD', 'CHOL', 'ESCA', 'GBM', 'HNSC', 'KIRP', 'KIRC', 'KICH', 'LGG', 'LUAD', 'LUSC', 'LIHC', 'LAML', 'OV', 'PRAD', 'PAAD', 'PCPG', 'READ', 'STAD', 'SKCM', 'THCA', 'TGCT', 'UCEC', 'UCS')

other_data = read.csv("~/diff/data/ESTIMATE.csv",,header=T,check.names=F)

temp_data = data.frame()

for (tumor_id in tumor_list)
{
  temp_other_data = other_data[grep(tumor_id,other_data$CODE),] %>% subset(select=-c(CODE,Group))
  temp_final_data = final_data[grep(tumor_id,final_data$CODE),] %>% subset(select=-c(CODE))
  all_name = merge(temp_other_data,temp_final_data,by='SampleName')$SampleName
  temp_other_data = temp_other_data[which(temp_other_data$SampleName%in%all_name),] %>% subset(select=-c(SampleName))
  temp_final_data = temp_final_data[which(temp_final_data$SampleName%in%all_name),] %>% subset(select=-c(SampleName))
  data.corr <- corr.test(temp_final_data, temp_other_data, method="pearson", adjust="fdr")
  tempp_data = data.frame(Row=rep(tumor_id,length(colnames(data.corr$p))),Col=colnames(data.corr$p),rvalue=data.corr$r[1,],pvalue=data.corr$p[1,])
  temp_data = rbind(temp_data,tempp_data)
}




# Row     Col       rvalue       pvalue
# 1    StromalScore   ATF6B -0.091240485 4.300041e-01
# 2    StromalScore     ADA  0.270801059 1.721441e-02
# 3    StromalScore   AHNAK  0.308557224 6.328550e-03
# 4    StromalScore  AKR1C4  0.141196499 2.206219e-01
# 5    StromalScore    AHCY -0.223922389 5.026419e-02
# 6    StromalScore    BATF  0.249722343 2.850313e-02
# 7    StromalScore     B2M  0.174181838 1.297605e-01
# 8    StromalScore   CLIC1 -0.006130566 9.577988e-01
# 9    StromalScore    CDK1 -0.159435151 1.660427e-01
# 10   StromalScore    CDK2 -0.142270447 2.170955e-01
# 11   StromalScore CYP27A1  0.075351859 5.148388e-01



data.r = data.corr$r
data.p = data.corr$p



pdf(paste("out.pdf",sep=''),width = 7,height = 10)
p = ggplot(temp_data, aes(x = Col, y = Row, color = rvalue, size = abs(-log10(pvalue)))) +
  geom_point() +  # 使用点表示数据
  scale_color_gradient2(low = '#d3f6d1', high = '#ff5126') +  # 调整颜色渐变,
  scale_size_continuous(range = c(1,10)) +  # 调整点的大小范围
  labs(title='',x = "", y = "", color = "R", size = "-log10(P)") +  # 添加轴标签
  ggtitle(tumor_id)+
  theme_bw()+
  theme(plot.title = element_text(size=30,color='black',face = "bold"),axis.text.x = element_text(angle = 45,size=20,color='black', hjust = 1),axis.text.y = element_text(size=15,color='black'))
print(p)
dev.off()
