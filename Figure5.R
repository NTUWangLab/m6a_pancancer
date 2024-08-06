# Figure 5A: Clustering heatmaps and consensus matrices for different samples

setwd("~/projects/m6a/cluster")

library(ggplot2)
library(ggpubr)
library(stringr)
library(ConsensusClusterPlus)
library(pheatmap)
library(ggplotify)

final_class = data.frame()

## Generate clustering heatmaps for various samples
#out_put_list <- c()

for(sample_name in c('ACC','BLCA','BRCA','CESC','COAD', 'CHOL', 'ESCA', 'GBM', 'HNSC', 'KIRP', 'KIRC', 'KICH', 'LGG', 'LUAD', 'LUSC', 'LIHC', 'LAML', 'OV', 'PRAD', 'PAAD', 'PCPG', 'READ', 'STAD', 'SKCM', 'THCA', 'TGCT', 'UCEC', 'UCS'))
{
  
  my_data = read.csv('~/projects/m6a/data/m6a.csv',header=T,check.names=F)
  my_data = my_data[which(my_data$CODE==sample_name),]
  rownames(my_data) = my_data$SampleName
  my_data = subset(my_data,select = -c(SampleName,CODE,Group))
  results = ConsensusClusterPlus(t(as.matrix(my_data)),
                                 maxK=4,
                                 reps=10,
                                 pItem=0.8,
                                 pFeature=1,
                                 title = sample_name,
                                 clusterAlg="pam",
                                 distance="pearson",
                                 seed=123456,
                                 plot="pdf", 
                                 writeTable = TRUE,
                                 tmyPal=c("#fef9ec","#feefe2","#ffe9ce","#ffe5bb","#fbddb2","#fed3a1","#ffbd89","#ffaf7a","#f89f6e","#fa8e61","#f77b4d","#f06c4c","#ea5841","#e04631","#d72f24","#ca1b17","#b4100c","#a70303","#930205")
  )
  n = calcICL(results,title="t",plot = NULL)
  temp_var_k = c()
  n$clusterConsensus = na.omit(n$clusterConsensus)
  
  for(i in 2:4)
  {
    temp_var_k = append(temp_var_k,mean(n$clusterConsensus[n$clusterConsensus[,1]==i,3]))
  }
  
  
  final_matrix = read.csv(paste0(sample_name,'/',sample_name,'.k=',which(temp_var_k==max(temp_var_k))+1,'.consensusMatrix.csv'),row.names = 1)
  
  pdf(paste0('cluster.pdf'))
  heatmap_obj = pheatmap(as.matrix(final_matrix),
                         cluster_rows = TRUE,
                         cluster_cols = TRUE)
  row_order = heatmap_obj$tree_row$order
  col_order = heatmap_obj$tree_col$order
  
  data_matrix_reordered = as.matrix(final_matrix)[row_order, col_order]
  
  class_matrix = read.table(paste0(sample_name,'/',sample_name,'.k=',which(temp_var_k==max(temp_var_k))+1,'.consensusClass.csv'),row.names = 1,sep=',')
  class_matrix$CODE = rep(sample_name,length(rownames(class_matrix)))
  colnames(class_matrix)[1] = 'Cluster'
  
  if(sample_name == 'ACC')
  {
    final_class = class_matrix
  }
  else
  {
    dev.off()
    final_class = rbind(final_class,class_matrix)
  }
  
  results[[which(temp_var_k==max(temp_var_k))+1]]$clrs
  
  hc = hclust(as.dist(1 - results[[which(temp_var_k==max(temp_var_k))+1]]$ml),method='average')
  data_matrix_reordered = results[[which(temp_var_k==max(temp_var_k))+1]]$ml[hc$order, ]
  data_matrix_reordered = rbind(data_matrix_reordered, 0)
  color = colorRampPalette(c("#fef9ec","#feefe2","#ffe9ce","#ffe5bb","#fbddb2","#fed3a1","#ffbd89","#ffaf7a","#f89f6e","#fa8e61","#f77b4d","#f06c4c","#ea5841","#e04631","#d72f24","#ca1b17","#b4100c","#a70303","#930205"))(79)
  
  pdf(paste0(sample_name,'/cluster.pdf'),height=6,width=5.7)
  heatmap(data_matrix_reordered, Colv = as.dendrogram(hc), Rowv = NA, symm = FALSE,
          scale = "none", col = color, na.rm = TRUE,
          labRow = F, labCol = F, mar = c(5, 5), ColSideCol = results[[which(temp_var_k==max(temp_var_k))+1]]$clrs[[1]])
  legend("topright", legend = unique(class_matrix$Cluster), fill = unique(results[[which(temp_var_k==max(temp_var_k))+1]]$clrs[[1]]), horiz = FALSE)
  graphics.off()
}

write.csv(final_class,file='out.csv')

# Figure 5B: Heatmap of significance values with -log10(p-value) transformation

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

myColor <- colorRampPalette(c( 'blue3', "white", "red3"))(paletteLength)

myBreaks <- c(seq(min(df), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(df)/paletteLength, max(df), length.out=floor(paletteLength/2)))

pdf(paste('pheatmap.pdf',sep=''),length(colnames(df))/2+10,length(rownames(df))/4+3)

xx <- pheatmap(df,
               color=myColor,
               breaks=myBreaks,
               clustering_method="average",main='-log10(P)',number_color='black',border_color = "black",fontsize = 15, cluster_rows=F,cluster_cols=F, cellwidth = 30,cellheight = 15,display_numbers=sig.mat)
print(xx)
dev.off()

# Figure 5C: Violin and box plots of gene expression in different tumor types

rm(list=ls())
file_dir = "~/projects/m6a/diff"
setwd(dir = file_dir)
library(ggplot2)
library(ggpubr)
library(stringr)
library(ggside)
library(dplyr)

tumor_list = c('ACC','BLCA','BRCA','CESC','COAD', 'CHOL', 'ESCA', 'GBM', 'HNSC', 'KIRP', 'KIRC', 'KICH', 'LGG', 'LUAD', 'LUSC', 'LIHC', 'LAML', 'OV', 'PRAD', 'PAAD', 'PCPG', 'READ', 'STAD', 'SKCM', 'THCA', 'TGCT', 'UCEC', 'UCS')

my_data = read.csv('~/projects/m6a/data/sangerbox-exp.csv',header=T,check.names=F)
colnames(my_data) = gsub('\n','',colnames(my_data))

gene = ''

exp_file <- my_data[,which(colnames(my_data)%in%c(gene,'CODE','Group'))]

exp_file = exp_file[which(str_extract(exp_file$CODE,'[A-Z]+')%in%tumor_list),]

dat <- data.frame(Gene=exp_file$CODE,Group = exp_file$Group,value=exp_file[,1],check.names = F)

dat[,3] = as.numeric(dat[,3])

######log
dat[,3] = log(dat[,3]+1)
######

dat <- na.omit(dat)

dat$Gene = str_extract(dat$Gene,'[A-Z]+')

pdf(paste(gene,".pdf",sep=''),width=length(unique(dat[,1]))+5,height = 10)

p <- ggplot(dat, aes(x = Gene, y = value, color = Group)) +
  geom_violin(scale = "width", adjust =1, trim = TRUE) +
  geom_boxplot(width=0.5,position=position_dodge(0.9)) +
  labs(x="Tumor",y = "Gene Expression(log(x+1))") +
  theme_classic()+
  scale_color_lancet()+
  theme(legend.key.size = unit(0.2, "inches"),axis.title = element_text(size = 30,color='black'),axis.text.x = element_text(size = 20,angle = 60,vjust = 1,hjust = 1,color='black'),axis.text.y = element_text(size = 30,color='black'))
print(p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova"))
dev.off()
