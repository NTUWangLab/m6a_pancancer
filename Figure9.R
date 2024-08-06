# Figure 13A,B: Generate a heatmap of -log10(P) values from a CSV file and save it as a PDF.
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

myColor <- colorRampPalette(c( "orange2", "white", "red3"))(paletteLength)

myBreaks <- c(seq(min(df), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(df)/paletteLength, max(df), length.out=floor(paletteLength/2)))

pdf(paste('pheatmap.pdf',sep=''),length(colnames(df))/2+10,length(rownames(df))/4+3)

xx <- pheatmap(df,
               color=myColor,
               breaks=myBreaks,
               clustering_method="average",main='-log10(P)',number_color='black',border_color = "black",fontsize = 15, cluster_rows=F,cluster_cols=F, cellwidth = 30,cellheight = 15,display_numbers=sig.mat)
print(xx)
dev.off()

# Figure 13C: Generate correlation heatmaps for various tumor types and compile them into a single PDF.
setwd("~/projects/m6a")
rm(list = ls())
library(dplyr)
library(patchwork)
library(ggplotify)
library(psych)

tumor_list = c('ACC','BLCA','BRCA','CESC','COAD', 'CHOL', 'ESCA', 'GBM', 'HNSC', 'KIRP', 'KIRC', 'KICH', 'LGG', 'LUAD', 'LUSC', 'LIHC', 'LAML', 'OV', 'PRAD', 'PAAD', 'PCPG', 'READ', 'STAD', 'SKCM', 'THCA', 'TGCT', 'UCEC', 'UCS')

out_put_list <- c()
for (tumor_id in tumor_list)
{
  print(tumor_id)
  
  dat_input <- read.csv("~/diff/data/gx.csv",check.names = F,header=T)
  dat_input <- dat_input[grep(tumor_id,dat_input$CODE),]
  dat_input <- dat_input[grep("Tumor",dat_input$Group),]
  dat_input <- subset(dat_input,select=-c(Group,CODE))
  dat_input <- dat_input[!duplicated(dat_input$SampleName),]
  rownames(dat_input) <- make.unique(sub('\n','',dat_input$SampleName))
  dat_input <- subset(dat_input,select=-c(SampleName))
  colscluster = length(colnames(dat_input))/3
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
  
  data.r <- data.corr$r
  data.p <- data.corr$p
  paletteLength = 1000
  
  myColor <- colorRampPalette(c("orange2", "white", "red3"))(paletteLength)
  
  getSig <- function(dc) {
    sc <- ' '
    if (dc < 0.0001) {sc <- '****'}
    else if (dc < 0.001){sc <- '***'}
    else if (dc < 0.01){sc <- '**'}
    else if (dc < 0.05) {sc <- '*'}
    else{sc <- ''}
    return(sc)
  }
  
  sig.mat <- matrix(sapply(data.p, getSig), nrow=nrow(data.p))
  str(sig.mat)
  
  paletteLength <- 1000
  
  test <- data.r
  myBreaks <- c(seq(min(test), 0, length.out=ceiling(paletteLength/2) + 1),
                seq(max(test)/paletteLength, max(test), length.out=floor(paletteLength/2)))
  
  pdf(paste('gx/corr/',tumor_id,".pdf",sep=''),width = length(colnames(dat_input)),height = length(colnames(dat_one_tumor))/2.5)
  
  p = pheatmap(data.r, 
               color=myColor,
               breaks=myBreaks,
               clustering_method="average",border_color='white',number_color = "black",fontsize_row = 20,fontsize_col = 20,fontsize=12, cluster_rows=F,cluster_cols=F,cellwidth = 20,cellheight = 20, display_numbers=sig.mat)
  dev.off()
  
  #
  eval(parse(text = paste0(tumor_id,'= as.ggplot(p)')))
  out_put_list <- append(out_put_list,tumor_id)
}
##
pdf("fig1.pdf",width = length(colnames(dat_input))*5,height = length(colnames(dat_one_tumor))*1.5+2)

eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=7,labels=out_put_list,label_size=40))',sep='')))

dev.off()
