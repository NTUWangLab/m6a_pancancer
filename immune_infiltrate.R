
rm(list=ls())

setwd("~/projects/m6a")


library(ggplot2)


tumor_list = c('ACC','BLCA','BRCA','CESC','COAD', 'CHOL', 'ESCA', 'GBM', 'HNSC', 'KIRP', 'KIRC', 'KICH', 'LGG', 'LUAD', 'LUSC', 'LIHC', 'LAML', 'OV', 'PRAD', 'PAAD', 'PCPG', 'READ', 'STAD', 'SKCM', 'THCA', 'TGCT', 'UCEC', 'UCS')

data_all_tumor = read.csv(paste("~/projects/m6a/data/m6a.csv",sep=''),check.names = F)

gene_list = colnames(data_all_tumor)[-c(1,2,length(colnames(data_all_tumor)))]



for(tumor_id in tumor_list)
{
  
  
  dat_input = read.csv("~/diff/data/ESTIMATE.csv",check.names = F,header=T)
  dat_input = dat_input[grep(tumor_id,dat_input$CODE),]
  dat_input = dat_input[grep("Tumor",dat_input$Group),]
  dat_input = subset(dat_input,select=-c(Group,CODE,StromalScore,ESTIMATEScore))
  dat_input = dat_input[!duplicated(dat_input$SampleName),]
  
  temp_gene_id = ''
  temp_gene_num = 0
  
  print(tumor_id)
  
  for(gene_id in gene_list)
  {
    
    data_all_tumor = read.csv(paste("~/projects/m6a/data/m6a.csv",sep=''),check.names = F)
    dat_one_tumor = data_all_tumor[grep(tumor_id,data_all_tumor$CODE),]
    dat_one_tumor = dat_one_tumor[,which(colnames(dat_one_tumor)%in%c(gene_id,'SampleName'))]
    
    out= merge(dat_input,dat_one_tumor, by = "SampleName")
    out = subset(out,select=-c(SampleName))
    
    out[,2] = log(out[,2]+1)
    
    fit = lm(out[,2] ~ out[,1])
    
    summary_fit = summary(fit)
    
    r_squared = summary_fit$r.squared
    p_value = summary_fit$coefficients[2, 4]
    
    if(temp_gene_num<r_squared)
    {
      temp_gene_num = r_squared
      temp_gene_id = gene_id
    }
    
    print(paste0(gene_id,' : ',r_squared))
  }
  
  
  
  data_all_tumor = read.csv(paste("~/projects/m6a/data/m6a.csv",sep=''),check.names = F)
  dat_one_tumor = data_all_tumor[grep(tumor_id,data_all_tumor$CODE),]
  dat_one_tumor = dat_one_tumor[,which(colnames(dat_one_tumor)%in%c(temp_gene_id,'SampleName'))]
  
  out= merge(dat_input,dat_one_tumor, by = "SampleName")
  out = rev(subset(out,select=-c(SampleName)))
  out[,1] = log(out[,1]+1)
  fit = lm(out[,2] ~ out[,1])
  
  summary_fit = summary(fit)
  
  r_squared = summary_fit$r.squared
  p_value = summary_fit$coefficients[2, 4]
  
  
  
  p =ggplot(out,aes(x=out[,1],y=out[,2]))+
    geom_point(shape=19) + 
    xlab(paste0(colnames(out)[1],'  Expression')) + 
    ylab(colnames(out)[2])+
    geom_smooth(method = lm)+
    theme_classic()+
    theme(ggside.panel.scale=0.1, ggside.axis.text=element_blank(),axis.text = element_text(color='black',size = 25),axis.title = element_text(color='black',size = 25))+
    annotate('text',size=8, y=max(out[,2])+max(out[,2])/5,x = (min(out[,1])+max(out[,1]))/2,label = paste(tumor_id,'(n=',nrow(out),"),R=",format(sqrt(round(r_squared, 3)),digits = 2), ",P=", format(p_value, scientific=TRUE, digits=2),sep=''),colour="black")+
    geom_xsidedensity(aes(x=out[,1]),fill='#06799f',show.legend = NA,stat = "density")+
    geom_ysidedensity(aes(y=out[,2]),fill='#06799f',show.legend=T,stat = "density")
  
  pdf(paste0(tumor_id,'.pdf'),width=8,height=8)
  print(p)
  dev.off()
}