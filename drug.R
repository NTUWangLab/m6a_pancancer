setwd('~/projects/m6a/drug')

rm(list = ls())
library(stringr)
options(stringsAsFactors = F)
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
options(check.names=F)

tumor_list = c('ACC','BLCA','BRCA','CESC','COAD', 'CHOL', 'ESCA', 'GBM', 'HNSC', 'KIRP', 'KIRC', 'KICH', 'LGG', 'LUAD', 'LUSC', 'LIHC', 'LAML', 'OV', 'PRAD', 'PAAD', 'PCPG', 'READ', 'STAD', 'SKCM', 'THCA', 'TGCT', 'UCEC', 'UCS')

GDSC2_Expr = readRDS(file='~/oncopredict/pancancer/GDSC2_Expr (RMA Normalized and Log Transformed).rds')
GDSC2_Res = readRDS(file = "~/oncopredict/pancancer/GDSC2_Res.rds")
GDSC2_Res = exp(GDSC2_Res)

exp = read.csv('~/projects/m6a/data/m6a.csv',header = T,check.names=F)
for (tumor_id in tumor_list)
{
  exp_file <- exp[grep(tumor_id,exp$CODE),]
  
  exp_file <- exp_file[grep("Tumor",exp_file$Group),]
  
  exp_file <- exp_file[,-1:-2]
  
  exp_file <- exp_file[!duplicated(exp_file$SampleName),]
  
  rownames(exp_file) <- gsub('\n','',exp_file[,length(exp_file)])
  
  exp_file <- exp_file[,-length(exp_file)]
  #trace("calcPhenotype",edit=TRUE)
  calcPhenotype(trainingExprData = GDSC2_Expr,
                trainingPtype = GDSC2_Res,
                testExprData = as.matrix(as.data.frame(t(exp_file))),
                batchCorrect = 'eb',  #   "eb" for ComBat  
                powerTransformPhenotype = TRUE,
                removeLowVaryingGenes = 0.2,
                minNumSamples = 1, 
                printOutput = TRUE, 
                removeLowVaringGenesFrom = 'rawData')
  
  temp_text <- read.csv('calcPhenotype_Output/DrugPredictions.csv',header=T,row.names=1,check.names=F)
  temp_text <- na.omit(temp_text)
  temp_text <- rbind(temp_text,Total=colSums(temp_text))
  rownames(temp_text) = sub('\n','',rownames(temp_text))
  temp_text <- data.frame(t(temp_text),check.names=F)
  temp_text <- temp_text[order(temp_text$Total),]
  temp_text <- temp_text[,-length(temp_text)]
  temp_text <- as.data.frame(t(head(temp_text,11)))
  for( i in 1:length(colnames(temp_text)))
  {
    colnames(temp_text)[i] <- sub(str_extract(colnames(temp_text)[i], "_.*"),'',colnames(temp_text)[i])
  }
  write.csv(temp_text,paste(tumor_id,'.csv',sep=''))
}