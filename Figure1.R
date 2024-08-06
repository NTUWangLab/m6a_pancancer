# Figure 1C
# This figure generates a scatter plot to show the correlation between gene and protein expression levels.

library(dplyr)
library(patchwork)
library(ggplotify)
library(ggplot2)

setwd("~/projects/m6a/protein_exp")

data.corr = read.csv('m6a-pro.csv', header=T, row.names=1)
data.corr[is.na(data.corr)] = 0

temp_data = data.frame()

for(i in colnames(as.data.frame(data.corr))) {
  for(j in rownames(as.data.frame(data.corr))) {
    eval(parse(text=paste0("temp_row=data.frame(Row='", i, "',Col='", j, "',rvalue=", "as.data.frame(data.corr)[,which(colnames(as.data.frame(data.corr))=='", i, "')][which(rownames(as.data.frame(data.corr))=='", j, "')])")))
    temp_data = rbind(temp_data, temp_row)
  }
}

pdf(paste('fig.pdf', sep=''), width=15, height=20)
p = ggplot(temp_data, aes(x = Col, y = Row, color = rvalue, size = rvalue)) +
  geom_point() +
  scale_color_gradient(low = "white", high = "#00A287") +
  scale_size_continuous(range = c(1, 10)) +
  labs(title='', x = "Gene", y = "", color = "Protein_exp", size = "Protein_exp") +
  ggtitle('Protein_exp') +
  theme_bw() +
  theme(plot.title = element_text(size=30, color='black', face = "bold"), axis.text.x = element_text(angle = 45, size=20, color='black', hjust = 1), axis.text.y = element_text(size=15, color='black'))
print(p)
dev.off()

# Figure 1D
# This figure generates box plots to show the differential expression of genes between Tumor and Normal groups.

rm(list=ls())
file_dir = "~/projects/m6a"
setwd(dir = file_dir)
library(ggplot2)
library(ggpubr)
library(stringr)
library(ggsci)
need_list = c('data/Group')

plus_symbol = 'Tumor'
minus_symbol = 'Normal'

get_symbol = function(data) {
  temp_symbol = c()
  for(i in unique(data$Gene)) {
    temp_data = data[which(data$Gene == i),]
    if((sum(temp_data[which(temp_data$Group == plus_symbol),]$value) / length(temp_data[which(temp_data$Group == plus_symbol),]$value)) > (sum(temp_data[which(temp_data$Group == minus_symbol),]$value) / length(temp_data[which(temp_data$Group == minus_symbol),]$value))) {
      temp_symbol = append(temp_symbol, '+')
    } else {
      temp_symbol = append(temp_symbol, '-')
    }
  }
  return(temp_symbol)
}

for(file_name in need_list) {
  tumor_list = c('ACC','BLCA','BRCA','CESC','COAD', 'CHOL', 'ESCA', 'GBM', 'HNSC', 'KIRP', 'KIRC', 'KICH', 'LGG', 'LUAD', 'LUSC', 'LIHC', 'LAML', 'OV', 'PRAD', 'PAAD', 'PCPG', 'READ', 'STAD', 'SKCM', 'THCA', 'TGCT', 'UCEC', 'UCS')
  
  out_put_list <- c()
  final_tumor_list = c()
  
  my_data = read.csv('data/sangerbox-exp.csv', header=T, check.names=F)
  other_data = read.csv(paste(file_name, ".csv", sep=''), header=T, check.names=F)
  p_value_csv <- data.frame()
  
  for(name in tumor_list) {
    dat <- data.frame(Gene=c(), Group=c(), value=c(), check.names = F)
    
    other_file <- other_data[grep(name, other_data$CODE),]
    
    if (length(rownames(other_file)) == 0) {
      next
    }
    
    other_file <- other_file[!duplicated(other_file$SampleName),]
    rownames(other_file) = gsub('\n', '', other_file$SampleName)
    other_file <- subset(other_file, select=-c(SampleName, CODE))
    
    exp_file <- my_data[grep(name, my_data$CODE),]
    
    colnames(exp_file) = gsub('\n', '', colnames(exp_file))
    exp_file <- exp_file[!duplicated(exp_file$SampleName),]
    rownames(exp_file) = gsub('\n', '', exp_file$SampleName)
    exp_file <- subset(exp_file, select=-c(Group, CODE, SampleName))
    
    gene_list <- colnames(exp_file)
    all_name <- names(which(table(c(rownames(other_file), rownames(exp_file))) == 2))
    
    if (length(all_name) == 0) {
      next
    }
    
    for(gene_name in gene_list) {
      for(i in all_name) {
        dat <- rbind(dat, c(gene_name, other_file[match(i, rownames(other_file)),], exp_file[c(gene_name)][match(i, rownames(exp_file)),]))
      }
      colnames(dat) <- c("Gene", "Group", "value")
    }
    
    dat[,3] = as.numeric(dat[,3])
    dat[,3] = log(dat[,3] + 1)
    dat <- na.omit(dat)
    
    xx <- compare_means(value ~ Group, data = dat, group.by = "Gene", method = "anova")
    p_value <- as.matrix(xx$p)
    p_value[is.na(p_value)] <- 1
    
    for(i in 1:length(p_value)) {
      if(get_symbol(dat)[i] == '-') {
        p_value[i] = -p_value[i]
      }
    }
    
    final_tumor_list <- append(final_tumor_list, name)
    
    print(name)
    
    if(length(final_tumor_list) != 1) {
      p_value_csv <- cbind(p_value_csv, as.data.frame(p_value))
    } else {
      p_value_csv <- as.data.frame(p_value)
    }
    
    pdf(paste('diff/', name, ".pdf", sep=''), width=length(unique(dat[,1]))-5, height=8)
    
    p <- ggboxplot(dat, x = "Gene", y = "value", color = "Group", palette = c('blue3', 'red3'), x.text.angle=60, add=c('jitter'))
    p <- p + xlab("Gene") + ylab("Gene Expression(log(x+1))")
    p <- p + theme(axis.text = element_text(size = 30), axis.title=element_text(size=30)) + stat_compare_means(aes(group = Group), label = "p.signif", method = "anova")
    
    print(p)
    eval(parse(text = paste(name, '<- p')))
    out_put_list <- append(out_put_list, name)
    dev.off()
  }
  rownames(p_value_csv) <- gene_list
  colnames(p_value_csv) <- final_tumor_list
  write.csv(p_value_csv, file=paste("p.csv", sep=''), quote=F)
  
  pdf(paste("fig1.pdf", sep=''), width = length(unique(dat[,1]))*4+3, height = 8*length(out_put_list)/3.5)
  eval(parse(text = paste('print(cowplot::plot_grid(', paste(out_put_list, collapse = ","), ', ncol=4, labels=out_put_list, label_size=40))', sep='')))
  dev.off()
}

# Figure 1E,F
# These figures generate violin plots and box plots to compare gene expression levels between Tumor and Normal groups for the specified gene.

rm(list=ls())
file_dir = "~/projects/m6a/diff"
setwd(dir = file_dir)
library(ggplot2)
library(ggpubr)
library(stringr)
library(ggside)
library(dplyr)

tumor_list = c('ACC','BLCA','BRCA','CESC','COAD', 'CHOL', 'ESCA', 'GBM', 'HNSC', 'KIRP', 'KIRC', 'KICH', 'LGG', 'LUAD', 'LUSC', 'LIHC', 'LAML', 'OV', 'PRAD', 'PAAD', 'PCPG', 'READ', 'STAD', 'SKCM', 'THCA', 'TGCT', 'UCEC', 'UCS')

my_data = read.csv('~/projects/m6a/data/sangerbox-exp.csv', header=T, check.names=F)
colnames(my_data) = gsub('\n', '', colnames(my_data))

gene = 'IGF2BP3'

exp_file <- my_data[, which(colnames(my_data) %in% c(gene, 'CODE', 'Group'))]
exp_file = exp_file[which(str_extract(exp_file$CODE, '[A-Z]+') %in% tumor_list),]

dat <- data.frame(Gene=exp_file$CODE, Group=exp_file$Group, value=exp_file[,1], check.names = F)
dat[,3] = as.numeric(dat[,3])
dat[,3] = log(dat[,3] + 1)
dat <- na.omit(dat)
dat$Gene = str_extract(dat$Gene, '[A-Z]+')

pdf(paste(gene, ".pdf", sep=''), width=length(unique(dat[,1]))+5, height=10)
p <- ggplot(dat, aes(x = Gene, y = value, color = Group)) +
  geom_violin(scale = "width", adjust=1, trim=TRUE) +
  geom_boxplot(width=0.5, position=position_dodge(0.9)) +
  labs(x="Tumor", y = "Gene Expression(log(x+1))") +
  theme_classic() +
  scale_color_lancet() +
  theme(legend.key.size = unit(0.2, "inches"), axis.title = element_text(size = 30, color='black'), axis.text.x = element_text(size = 20, angle = 60, vjust = 1, hjust = 1, color='black'), axis.text.y = element_text(size = 30, color='black'))
print(p + stat_compare_means(aes(group = Group), label = "p.signif", method = "anova"))
dev.off()
