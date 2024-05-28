setwd("~/projects/m6a/mirna")


tumor_list = c('ACC','BLCA','BRCA','CESC','COAD', 'CHOL', 'ESCA', 'GBM', 'HNSC', 'KIRP', 'KIRC', 'KICH', 'LGG', 'LUAD', 'LUSC', 'LIHC', 'LAML', 'OV', 'PRAD', 'PAAD', 'PCPG', 'READ', 'STAD', 'SKCM', 'THCA', 'TGCT', 'UCEC', 'UCS')

function_a = function(tumor_id) {
  library(dplyr)
  library(patchwork)
  library(ggplotify)
  library(psych)
  library(stringr)
  library(parallel)
  dat_input <- read.table(paste0("~/MIRNA_DATA/",tumor_id),sep='\t',check.names=F,header=T,row.names=1)
  dat_input = as.data.frame(t(dat_input))   
  na.omit(dat_input)
  dat_input$SampleName = str_sub(rownames(dat_input),end=-2)
  dat_input <- dat_input[!duplicated(dat_input$SampleName),]
  rownames(dat_input) <- make.unique(sub('\n','',dat_input$SampleName))
  dat_input <- subset(dat_input,select=-c(SampleName))
  
  
  data_all_tumor <- read.csv(paste("~/projects/m6a/data/m6a.csv",sep=''),check.names = F)
  dat_one_tumor <- data_all_tumor[grep(tumor_id,data_all_tumor$CODE),]
  rownames(dat_one_tumor) <- make.unique(sub('\n','',dat_one_tumor$SampleName))
  dat_one_tumor <- dat_one_tumor[,!duplicated(colnames(dat_one_tumor))]
  dat_one_tumor <- dat_one_tumor[grep("Tumor",dat_one_tumor$Group),]
  dat_one_tumor <- subset(dat_one_tumor,select=-c(Group,SampleName,CODE))
  one_tumor_sample <- unlist(rownames(dat_one_tumor))
  all_name <- names(which(table(c(rownames(dat_input),one_tumor_sample))==2))
  dat_gene <- dat_one_tumor[match(all_name,rownames(dat_one_tumor)),]
  dat_im <- dat_input[match(all_name,rownames(dat_input)),]
  
  if(length(rownames(dat_gene))==0)
  {
    next
  }
  data.corr <- corr.test(dat_im, dat_gene, method="pearson", adjust="fdr")
  
  
  
  
  
  
  data.r = data.corr$r
  data.p = data.corr$p
  getTotal_p <- function(dc) {
    if(is.na(dc)){return(0)}
    if(dc<0.05) return(1) else return(0)
  }
  
  getTotal_r <- function(dc) {
    if(is.na(dc)){return(0)}
    if(dc<(-0.3)||dc>0.3) return(1) else return(0)
  }
  
  for(num in 10:1000)
  {
    temp.sum = as.data.frame(matrix(sapply(data.p, getTotal_p), nrow=nrow(data.p)))
    rownames(temp.sum) = rownames(data.p)
    name = names(rowSums(temp.sum)[order(-rowSums(temp.sum))])[1:num]
    temp.sum = as.data.frame(matrix(sapply(data.p, getTotal_r), nrow=nrow(data.r)))
    rownames(temp.sum) = rownames(data.r)
    name = name[which(name%in%names(rowSums(temp.sum)[order(-rowSums(temp.sum))])[1:num])]
    if(length(name)>=10)
    {
      break
    }
    print(num)
  }
  name = name[1:10]
  data.r = as.data.frame(data.r)[which(rownames(as.data.frame(data.r))%in%name),]
  data.p = as.data.frame(data.p)[which(rownames(as.data.frame(data.p))%in%name),]
  
  
  write.csv(data.p,file=paste(tumor_id,"_out_p.csv",sep=''),quote=F)
  write.csv(data.r,file=paste(tumor_id,"_out_r.csv",sep=''),quote=F)
  
  
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
  
  
  pdf(paste(tumor_id,".pdf",sep=''),width = length(colnames(dat_input))/3+3,height = length(colnames(dat_one_tumor))/2.5)
  
  p = pheatmap(data.r, 
               color=myColor,
               breaks=myBreaks,
               clustering_method="average",border_color='white',number_color = "black",fontsize_row = 20,fontsize_col = 20,fontsize=12, cluster_rows=F,cluster_cols=F,cellwidth = 20,cellheight = 20, display_numbers=sig.mat)
  eval(parse(text = paste0(tumor_id,'= as.ggplot(p)')))
  eval(parse(text = paste0('save(',tumor_id,",file='",tumor_id,".Rda')")))
  dev.off()
}

clus <- makeCluster(2)
#可以确定当前可用CPU数目
detectCores()

parLapply(clus, tumor_list ,fun = function_a)