# Part 1: Generates p-values for comparison of gene expression between tumor and non-tumor samples
rm(list=ls())
file_dir = "~/projects/m6a"
setwd(dir = file_dir)
library(ggplot2)
library(ggpubr)
library(stringr)

tumor_list = c('ACC','BLCA','BRCA','CESC','COAD', 'CHOL', 'ESCA', 'GBM', 'HNSC', 'KIRP', 'KIRC', 'KICH', 'LGG', 'LUAD', 'LUSC', 'LIHC', 'LAML', 'OV', 'PRAD', 'PAAD', 'PCPG', 'READ', 'STAD', 'SKCM', 'THCA', 'TGCT', 'UCEC', 'UCS')

out_put_list <- c()
file_name = 'fenqi/AGE'
final_tumor_list = c()

my_data = read.csv('~/projects/m6a/data/m6a.csv',header=T,check.names=F)
other_data = read.csv("~/diff/data/X.csv",header=T,check.names=F) # X.csv is clinic data

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
  
  exp_file <- exp_file[grep("Tumor",exp_file$Group),]
  
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
    
  }
  
  colnames(dat) <- c("Gene","Group","value")
  
  dat[,3] = as.numeric(dat[,3])
  dat[,3] = log(dat[,3]+1)
  dat <- na.omit(dat)
  xx <-compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova")
  
  p_value <- as.matrix(xx$p)
  
  p_value[is.na(p_value)] <- 1
  
  final_tumor_list <- append(final_tumor_list,name)
  print(name)
  if(length(final_tumor_list)!=1){p_value_csv <- cbind(p_value_csv,as.data.frame(p_value))}else{p_value_csv <- as.data.frame(p_value)}
}
rownames(p_value_csv) <- gene_list
colnames(p_value_csv) <- final_tumor_list
write.csv(p_value_csv,file=paste(file_name,"/p.csv",sep=''),quote=F)

# Part 2: Generates a heatmap of -log10(p-values) with significance annotations
rm(list=ls())
library("pheatmap")
library("jsonlite")

setwd(dir = "~/projects/m6a")
temp = list.files(pattern="*.csv")

df = read.csv(paste('p.csv',sep=''),header=T,row.names=1)
df = replace(df,is.na(df),1)

paletteLength = 1000

getSig <- function(dc) {
  sc <- ' '
  dc = abs(dc)
  if (dc < 0.0001) {sc <- '****'}
  else if (dc < 0.001){sc <- '***'}
  else if (dc < 0.01){sc <- '**'}
  else if (dc < 0.05) {sc <- '*'}
  else{sc <- ''}
  return(sc)
}

sig.mat <- matrix(sapply(as.matrix(df), getSig), nrow=nrow(as.matrix(df)))
str(sig.mat)

tran = function(temp)
{
  temp_p = c()
  for(i in temp)
  {
    if(i > 0)
    {
      if(i < 0.000000001)
      {
        i =  0.000000001
      }
      if(-log10(i)!=0){temp_p = append(temp_p,-log10(i))}else{temp_p = append(temp_p,0.99)}
    }
    else
    {
      if(i > -0.000000001)
      {
        i =  -0.000000001
      }
      temp_p = append(temp_p,log10(abs(i)))
    }
  }
  return(temp_p)
}

df_temp = sapply(df,tran)
rownames(df_temp) = rownames(df)
df = df_temp

myColor <- colorRampPalette(c( 'yellow2', "white", "red3"))(paletteLength)

myBreaks <- c(seq(min(df), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(df)/paletteLength, max(df), length.out=floor(paletteLength/2)))

pdf(paste('pheatmap.pdf',sep=''),length(colnames(df))/2+10,length(rownames(df))/4+3)

xx <- pheatmap(df,
               color=myColor,
               breaks=myBreaks,
               clustering_method="average",main='-log10(P)',number_color='black',border_color = "black",fontsize = 15, cluster_rows=F,cluster_cols=F, cellwidth = 30,cellheight = 15,display_numbers=sig.mat)
print(xx)
dev.off()

# Part 3: Creates boxplots of gene expression for each gene across different tumor types
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
other_data = read.csv("~/diff/data/M.csv",header=T,check.names=F)

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

gene = ''

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
dev.off()
