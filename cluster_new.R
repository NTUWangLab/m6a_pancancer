

setwd("~/projects/m6a/cluster")

library(ggplot2)
library(ggpubr)
library(stringr)
library(ConsensusClusterPlus)
library(pheatmap)
library(ggplotify)

final_class = data.frame()

##拼图
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
  # dev.new()
  # heatmap(data_matrix_reordered, Colv = as.dendrogram(hc), Rowv = NA, symm = FALSE, 
  #         scale = "none", col = color, na.rm = TRUE,
  #         labRow = F, labCol = F, mar = c(5, 5), ColSideCol = results[[which(temp_var_k==max(temp_var_k))+1]]$clrs[[1]])
  # legend("topright", legend = unique(class_matrix$Cluster), fill = unique(results[[which(temp_var_k==max(temp_var_k))+1]]$clrs[[1]]), horiz = FALSE)
  # Sys.sleep(2)
  # eval(parse(text = paste(sample_name,'<- recordPlot()')))
  # out_put_list <- append(out_put_list,sample_name)
}

write.csv(final_class,file='out.csv')


# pdf("fig1.pdf",width = 20,height = 23)
# 
# eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=5))',sep='')))
# 
# dev.off()