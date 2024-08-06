setwd("~/projects/m6a/cortri")

library(ggplotify)

tumor_list = c('ACC', 'BRCA', 'BLCA', 'CESC', 'COAD', 'CHOL', 'ESCA', 'GBM', 'HNSC', 'KIRP', 'KIRC', 'KICH', 'LGG', 'LUAD', 'LUSC', 'LIHC', 'LAML', 'OV', 'PRAD', 'PAAD', 'PCPG', 'READ', 'STAD', 'SKCM', 'THCA', 'TGCT', 'UCEC', 'UCS')

my_data = read.csv('~/projects/m6a/data/m6a.csv', header=T, check.names=F)

out_put_list <- c()

for(i in tumor_list)
{
  res <- my_data[grep(i, my_data$CODE),]
  res <- res[grep("Tumor", res$Group),]
  res <- res[,-1:-2]
  rownames(res) = gsub('\n', '', res[, length(res)])
  res <- res[,1:length(res)-1]
  library(corrplot)
  final <- cor(res)
  
  # Create a PDF file for each tumor type showing a correlation matrix
  pdf(paste(i, ".pdf", sep=''), width=length(colnames(final))/2, height=length(colnames(final))/2)
  
  temp = corrplot(final, type = "upper", order = 'alphabet', 
                  tl.col = "black", tl.srt = 90, tl.cex = 2, method = "square", na.label = i)
  temp = temp$corr
  write.csv(temp, file=paste0(i, '.csv'))
  dev.off()
  
  # Plot the correlation matrix with pie chart symbols and save it as a record plot
  corrplot(final, type = "upper", order = "alphabet", 
           tl.col = "black", tl.cex = 0.5, method = "square", upper = "pie", tl.srt = 90)
  eval(parse(text = paste(i, '<- recordPlot()')))
  out_put_list <- append(out_put_list, c(i, NULL))
}

# Combine all individual plots into one large PDF file
pdf("fig1.pdf", width=20, height=30)

eval(parse(text = paste('print(cowplot::plot_grid(', paste(out_put_list, collapse = ","), ', ncol=5, labels=out_put_list, label_size=40))', sep='')))

dev.off()
