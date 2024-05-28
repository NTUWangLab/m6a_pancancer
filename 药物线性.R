setwd("~/projects/m6a")
rm(list = ls())
library(dplyr)
library(patchwork)
library(ggplotify)
library(psych)
library(ggside)

tumor_list = c('ACC','BLCA','BRCA','CESC','COAD', 'CHOL', 'ESCA', 'GBM', 'HNSC', 'KIRP', 'KIRC', 'KICH', 'LGG', 'LUAD', 'LUSC', 'LIHC', 'LAML', 'OV', 'PRAD', 'PAAD', 'PCPG', 'READ', 'STAD', 'SKCM', 'THCA', 'TGCT', 'UCEC', 'UCS')

other_data = read.csv("~/projects/m6a/Lasso/risk.csv",row.names=1,header=T,check.names=F)

out_put_list <- c()
for (tumor_id in tumor_list)
{
  print(tumor_id)
  dat_input <- read.csv(paste0("drug/data/",tumor_id,".csv"),check.names = T,row.names=1,header=T)
  dat_input <- dat_input[!duplicated(rownames(dat_input)),]
  
  
  status = 1
  #显著性
  dat <- data.frame(check.names = F)
  other_file <- other_data[grep(tumor_id,other_data$CODE),]
  other_file <- subset(other_file,select=-c(CODE))
  allname <- names(which(table(c(rownames(other_file),rownames(dat_input)))==2))
  if(length(allname)!=0)
  {
    for(gene_name in colnames(dat_input))
    {
      for(i in allname)
      {
        dat <- rbind(dat,c(gene_name,other_file[match(i,rownames(other_file)),],dat_input[c(gene_name)][match(i,rownames(dat_input)),]))
      }
      
    }
    colnames(dat) <- c("Gene","Group","value")
    dat[,3] = as.numeric(dat[,3])
    dat <- na.omit(dat)
    xx <- compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova")
  }else{
    status=0
  }
  
  
  #线性
  data_all_tumor <- read.csv(paste("~/projects/m6a/data/m6a.csv",sep=''),check.names = F)
  dat_one_tumor <- data_all_tumor[grep(tumor_id,data_all_tumor$CODE),]
  if (length(colnames(dat_one_tumor))==0)
  {
    next
  }
  rownames(dat_one_tumor) <- make.unique(sub('\n','',dat_one_tumor$SampleName))
  dat_one_tumor <- dat_one_tumor[,!duplicated(colnames(dat_one_tumor))]
  dat_one_tumor <- dat_one_tumor[grep("Tumor",dat_one_tumor$Group),]
  dat_one_tumor <- subset(dat_one_tumor,select=-c(Group,SampleName,CODE))
  one_tumor_sample <- unlist(rownames(dat_one_tumor))
  all_name <- names(which(table(c(rownames(dat_input),one_tumor_sample))==2))
  dat_gene <- dat_one_tumor[match(all_name,rownames(dat_one_tumor)),]
  dat_im <- dat_input[match(all_name,rownames(dat_input)),]
  data.corr <- corr.test(dat_gene, dat_im, method="pearson", adjust="fdr")
  temp_data = data.frame()
  for(i in colnames(as.data.frame(data.corr$r)))
  {
    for(j in rownames(as.data.frame(data.corr$r)))
    {
      eval(parse(text=paste0("temp_row=data.frame(Row='",i,"',Col='",j,"',rvalue=","as.data.frame(data.corr$r)[,which(colnames(as.data.frame(data.corr$r))=='",i,"')][which(rownames(as.data.frame(data.corr$r))=='",j,"')]",",pvalue=","as.data.frame(data.corr$p)[,which(colnames(as.data.frame(data.corr$p))=='",i,"')][which(rownames(as.data.frame(data.corr$p))=='",j,"')])")))
      temp_data = rbind(temp_data,temp_row)
    }
  }
  
  
  
  
  
  #geom_text(aes(label = Label), hjust = 0, vjust = 0, size = 5)
  pdf(paste('drug/corr/',tumor_id,".pdf",sep=''),width = 15,height = 5)
  p = ggplot(temp_data, aes(x = Col, y = Row, color = rvalue, size = abs(-log10(pvalue)))) +
    geom_point() +  # 使用点表示数据
    scale_color_gradient2(low = "green3", high = "red3") +  # 调整颜色渐变
    scale_size_continuous(range = c(1,10)) +  # 调整点的大小范围
    labs(title=tumor_id,x = "", y = "", color = "R", size = "-log10(P)") +  # 添加轴标签
    ggtitle(tumor_id)
    if(status == 1)
    {
      p= p+geom_ysidelabel(aes(label=p,x=x,y=y),data=data.frame(p=xx$p.signif[order(xx$Gene)],x=rep(0,length(xx$p.signif)),y=seq(1:length(xx$p.signif))),inherit.aes=F)+
        theme_bw()+
        theme(plot.title = element_text(size=30,color='black',face = "bold"),axis.text.x = element_text(angle = 45,size=20,color='black', hjust = 1),axis.text.y = element_text(size=15,color='black'))
    }else{
      p=p+theme_bw()+
      theme(plot.title = element_text(size=30,color='black',face = "bold"),axis.text.x = element_text(angle = 45,size=20,color='black', hjust = 1),axis.text.y = element_text(size=15,color='black'))
    }
    print(p)
  dev.off()
  
  #拼图
  eval(parse(text = paste0(tumor_id,'= p')))
  out_put_list <- append(out_put_list,tumor_id)
}
##拼图
pdf("drug/corr/fig1.pdf",width = length(colnames(dat_input))*5.7,height = length(colnames(dat_one_tumor))*length(out_put_list)/20)

eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4))',sep='')))

dev.off()
