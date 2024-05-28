rm(list=ls())
file_dir = "~/projects/m6a"
setwd(dir = file_dir)
library(ggplot2)
library(ggpubr)
library(stringr)

data_list = c('chemokine','Immunoinhibitor','Immunostimulator','MHC','receptor')

for(classify in data_list)
{
  tumor_list = c('ACC','BLCA','BRCA','CESC','COAD', 'CHOL', 'ESCA', 'GBM', 'HNSC', 'KIRP', 'KIRC', 'KICH', 'LGG', 'LUAD', 'LUSC', 'LIHC', 'LAML', 'OV', 'PRAD', 'PAAD', 'PCPG', 'READ', 'STAD', 'SKCM', 'THCA', 'TGCT', 'UCEC', 'UCS')
  out_put_list <- c()
  file_name = paste0('immune_reg/',classify,'/cluster')
  
  if (!dir.exists(file_name)){
    dir.create(paste0('immune_reg/',classify))
    dir.create(file_name)
  } else {
    print("Dir already exists!")
  }
  final_tumor_list = c()
  
  my_data = read.csv(paste0('~/diff/data/immune_reg/',classify,'.csv'),header=T,check.names=F)
  
  other_data = read.csv("~/projects/m6a/cluster/cluster.csv",,header=T,check.names=F)
  
  p_value_csv <- data.frame()
  
  for(name in tumor_list)
  {
    
    dat <- data.frame(check.names = F)
    other_file <- other_data[grep(name,other_data$CODE),]
    if (length(rownames(other_file))==0)
    {
      next
    }
    other_file <- other_file[!duplicated(other_file$SampleName),]
    rownames(other_file) = gsub('\n','',other_file$SampleName)
    other_file <- subset(other_file,select=-c(SampleName,CODE))
    
    exp_file <- my_data[grep(name,my_data$CODE),]
    colnames(exp_file) = gsub('\n','',colnames(exp_file))
    rownames(exp_file) = gsub('\n','',exp_file$SampleName)
    exp_file <- exp_file[!duplicated(exp_file$SampleName),]
    exp_file <- subset(exp_file,select=-c(CODE,SampleName))
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
      
    }
    
    colnames(dat) <- c("Gene","Group","value")
    
    dat[,3] = as.numeric(dat[,3])
    
    dat <- na.omit(dat)
    
    dat[,3] = log(dat[,3]+1)
    xx <-compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova")
    
    # ##剔除ns
    # temp_name <- xx[which(xx$p.signif!='ns'),]$Gene
    # dat <- dat[dat[,1]%in%temp_name,]
    ##
    
    p_value <- as.matrix(xx$p)
    
    p_value[is.na(p_value)] <- 1
    
    
    
    final_tumor_list <- append(final_tumor_list,name)
    print(name)
    if(length(final_tumor_list)!=1){p_value_csv <- cbind(p_value_csv,as.data.frame(p_value))}else{p_value_csv <- as.data.frame(p_value)}
    
    
    pdf(paste(file_name,"/",name,".pdf",sep=''),width=length(unique(dat[,1])),height = 10)
    #填充
    # p <- ggboxplot(dat, x = "Gene", y = "value",
    #                 palette = "jama",
    #                 fill = "Group",x.text.angle=60,x.text.size=20)
    #jitter
    #换样式nejm jco lancet jama  c("#BD2989","#B7E731")
    p <- ggboxplot(dat, x = "Gene", y = "value",
                   color = "Group", palette = c('#f85f73','#1e2022','#f8b500','#34495e'),
                   x.text.angle=60)
     #中间凹陷
    # p <- ggboxplot(dat, x = "Gene", y = "value",
    #                fill = "Group", palette = 'lancet',
    #                shape = "rx",,adjust=2,notch=T,x.text.angle=60)
    
    # p <- ggplot(dat, aes(x = Gene, y = value, fill = Group)) +
    #   # geom_boxplot(outlier.colour="black", outlier.shape=16,
    #   #              outlier.size=2, notch=FALSE)+
    #   theme(axis.text.x = element_text(angle=60))+
    #   geom_jitter(alpha=1,
    #               position=position_jitterdodge(jitter.width = 0.35, 
    #                                             jitter.height = 0, 
    #                                             dodge.width = 0.8))+
    #   geom_boxplot(width=0.5,position=position_dodge(0.9))+
    #   theme_classic()
    
    p <- p + xlab(classify)+ylab("Gene Expression(log(x+1))")
    p <- p + theme(axis.text = element_text(size = 30),axis.title=element_text(size=30))
    #p <- p + theme(axis.text = element_text(size = 15),axis.title=element_text(size=30))
    print(p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova"))
    
    ##拼图
    eval(parse(text = paste(name,'<- p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova")')))
    out_put_list <- append(out_put_list,name)
    ##
    
    dev.off()
  }
  rownames(p_value_csv) <- gene_list
  colnames(p_value_csv) <- final_tumor_list
  write.csv(p_value_csv,file=paste(file_name,"/p.csv",sep=''),quote=F)
  
  ######
  pdf(paste(file_name,"/fig1.pdf",sep=''),width = length(colnames(exp_file))*4,height = 8*length(out_put_list)/3.5)
  
  eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4,labels=out_put_list,label_size=40))',sep='')))
  
  dev.off()
  ######
}

