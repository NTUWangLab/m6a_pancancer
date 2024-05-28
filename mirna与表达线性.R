setwd("~/projects/m6a/mirna")


tumor_list = c('ACC','BLCA','BRCA','CESC','COAD', 'CHOL', 'ESCA', 'HNSC', 'KIRP', 'KIRC', 'KICH', 'LGG', 'LUAD', 'LUSC', 'LIHC', 'LAML', 'OV', 'PRAD', 'PAAD', 'PCPG', 'READ', 'STAD', 'SKCM', 'THCA', 'TGCT', 'UCEC', 'UCS')

library(dplyr)
library(patchwork)
library(ggplotify)
library(psych)
library(stringr)
out_put_list <- c()
for(tumor_id in tumor_list)
{
  
  
  
  data.r = as.matrix(read.csv(paste0(tumor_id,'_out_r.csv'),header=T,row.names=1))
  data.p = as.matrix(read.csv(paste0(tumor_id,'_out_p.csv'),header=T,row.names=1))
  
  paletteLength = 1000
  
  myColor <- colorRampPalette(c("#7109AA", "white", "#FFFF00"))(paletteLength)
  
  getSig <- function(dc) {
    sc <- ' '
    if (dc < 0.0001) {sc <- '****'}
    else if (dc < 0.001){sc <- '***'}
    else if (dc < 0.01){sc <- '**'}
    else if (dc < 0.05) {sc <- '*'}
    else{sc <- ''
    }
    return(sc)
  }
  
  sig.mat <- matrix(sapply(data.p, getSig), nrow=nrow(data.p))
  str(sig.mat)
  
  paletteLength <- 1000
  
  
  test <- data.r
  myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1),
                seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
  
  
  pdf(paste(tumor_id,".pdf",sep=''),width = length(colnames(data.r))/3+3,height = length(rownames(data.r))/2)
  
  p = pheatmap(data.r, 
               color=myColor,
               breaks=myBreaks,
               clustering_method="average",border_color='white',number_color = "black",fontsize_row = 20,fontsize_col = 20,fontsize=12, cluster_rows=F,cluster_cols=F,cellwidth = 20,cellheight = 20, display_numbers=sig.mat)
  eval(parse(text = paste0(tumor_id,'= as.ggplot(p)')))
  out_put_list <- append(out_put_list,tumor_id)
  dev.off()
}

pdf("fig1.pdf",width = length(colnames(data.r))*1.5,height = length(rownames(data.r))*3+11)

eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4,labels=out_put_list,label_size=40))',sep='')))

dev.off()

