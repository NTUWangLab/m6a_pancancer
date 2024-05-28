rm(list=ls())
library("pheatmap")
library("jsonlite")

setwd(dir = "~/projects/m6a")
temp = list.files(pattern="*.csv")




df = read.csv(paste('差异/p.csv',sep=''),header=T,row.names=1)
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


 
#immune_reg
#myColor <- colorRampPalette(c( 'purple3', "white", '#f9c535'))(paletteLength)
#CIBER
#myColor <- colorRampPalette(c( '#08d9d6', "white", '#ff8260'))(paletteLength)
#gx
#myColor <- colorRampPalette(c( "orange2", "white", "red3"))(paletteLength)
#yzx
#myColor <- colorRampPalette(c( 'yellow2', "white", "red3"))(paletteLength)
#exp
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

