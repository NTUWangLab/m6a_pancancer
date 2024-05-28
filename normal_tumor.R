rm(list=ls())
file_dir = "~/projects/m6a"
setwd(dir = file_dir)
library(ggplot2)
library(ggpubr)
library(stringr)
library(ggsci)
need_list = c('data/Group')


plus_symbol='Tumor'
minus_symbol='Normal'

get_symbol = function(data)
{
  temp_symbol = c()
  for(i in unique(data$Gene))
  {
    temp_data = data[which(data$Gene==i),]
    if((sum(temp_data[which(temp_data$Group==plus_symbol),]$value)/length(temp_data[which(temp_data$Group==plus_symbol),]$value))>(sum(temp_data[which(temp_data$Group==minus_symbol),]$value)/length(temp_data[which(temp_data$Group==minus_symbol),]$value)))
    {
      temp_symbol = append(temp_symbol,'+')
    }
    else
    {
      temp_symbol = append(temp_symbol,'-')
    }
  }
  return(temp_symbol)
}


for(file_name in need_list)
{
  tumor_list = c('ACC','BLCA','BRCA','CESC','COAD', 'CHOL', 'ESCA', 'GBM', 'HNSC', 'KIRP', 'KIRC', 'KICH', 'LGG', 'LUAD', 'LUSC', 'LIHC', 'LAML', 'OV', 'PRAD', 'PAAD', 'PCPG', 'READ', 'STAD', 'SKCM', 'THCA', 'TGCT', 'UCEC', 'UCS')
  
  ######拼图列表
  out_put_list <- c()
  ######
  
  final_tumor_list = c()
  
  my_data = read.csv('data/sangerbox-exp.csv',header=T,check.names=F)
  
  other_data = read.csv(paste(file_name,".csv",sep=''),header=T,check.names=F)
  
  p_value_csv <- data.frame()
  
  for(name in tumor_list)
  {
    
    dat <- data.frame(Gene=c(),Group = c(),value=c(),check.names = F)
    
    other_file <- other_data[grep(name,other_data$CODE),]
    
    if (length(rownames(other_file))==0)
    {
      next
    }
    
    other_file <- other_file[!duplicated(other_file$SampleName),]
    
    rownames(other_file) = gsub('\n','',other_file$SampleName)
    
    other_file <- subset(other_file,select=-c(SampleName,CODE))
    
    #分群文件筛选tumor
    # other_file <- other_file[grep("Tumor",other_file$Group),]
    # other_file <- subset(other_file,select=-c(Group))
    
    exp_file <- my_data[grep(name,my_data$CODE),]
    
    colnames(exp_file) = gsub('\n','',colnames(exp_file))
    
    exp_file <- exp_file[!duplicated(exp_file$SampleName),]
    
    rownames(exp_file) = gsub('\n','',exp_file$SampleName)
    
    exp_file <- subset(exp_file,select=-c(Group,CODE,SampleName))
    
    gene_list <- colnames(exp_file)
    
    all_name <- names(which(table(c(rownames(other_file),rownames(exp_file)))==2))
    
    
    if (length(all_name)==0)
    {
      next
    }
    
    for(gene_name in gene_list)
    {
      for(i in all_name)
      {
        dat <- rbind(dat,c(gene_name,other_file[match(i,rownames(other_file)),],exp_file[c(gene_name)][match(i,rownames(exp_file)),]))
      }
      colnames(dat) <- c("Gene","Group","value")
    }
    
    dat[,3] = as.numeric(dat[,3])
    
    ######log
    dat[,3] = log(dat[,3]+1)
    ######
    
    dat <- na.omit(dat)
    
    xx <- compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova")
    
    # ##剔除ns
    # temp_name <- xx[which(xx$p.signif!='ns'),]$Gene
    # dat <- dat[dat[,1]%in%temp_name,]
    ##
    
    p_value <- as.matrix(xx$p)
    
    p_value[is.na(p_value)] <- 1
    
    for(i in 1:length(p_value))
    {
      if(get_symbol(dat)[i]=='-')
      {p_value[i]=-p_value[i]}
    }
    
    
    final_tumor_list <- append(final_tumor_list,name)
    
    print(name)
    
    if(length(final_tumor_list)!=1){p_value_csv <- cbind(p_value_csv,as.data.frame(p_value))}else{p_value_csv <- as.data.frame(p_value)}
    
    pdf(paste('差异/',name,".pdf",sep=''),width=length(unique(dat[,1]))-5,height = 8)
    
    p <- ggboxplot(dat, x = "Gene", y = "value",
                   color = "Group", palette = c('blue3','red3'),
                   x.text.angle=60,add=c('jitter'))
    p <- p + xlab("Gene")+ylab("Gene Expression(log(x+1))")
    p <- p + theme(axis.text = element_text(size = 30),axis.title=element_text(size=30))  +stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")
    
    print(p)
    ##拼图
    eval(parse(text = paste(name,'<- p')))
    out_put_list <- append(out_put_list,name)
    ##
    dev.off()
  }
  rownames(p_value_csv) <- gene_list
  colnames(p_value_csv) <- final_tumor_list
  write.csv(p_value_csv,file=paste("p.csv",sep=''),quote=F)
  
  ######
  pdf(paste("fig1.pdf",sep=''),width = length(unique(dat[,1]))*4+3,height = 8*length(out_put_list)/3.5)
  
  eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4,labels=out_put_list,label_size=40))',sep='')))
  
  dev.off()
  ######
}