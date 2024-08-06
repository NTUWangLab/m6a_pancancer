# Figure 6A,B
setwd('~/projects/m6a/Lasso')

Tumor <- c('ACC','BLCA','BRCA','CESC','COAD', 'CHOL', 'ESCA', 'GBM', 'HNSC', 'KIRP', 'KIRC', 'KICH', 'LGG', 'LUAD', 'LUSC', 'LIHC', 'LAML', 'OV', 'PRAD', 'PAAD', 'PCPG', 'READ', 'STAD', 'SKCM', 'THCA', 'TGCT', 'UCEC', 'UCS')

final_list = data.frame()

for(name_fold in Tumor)
{
  all_sur_data <- read.csv('~/history_lasso/lassoPancancer/survival.csv', header=T, check.names=F)
  all_sur_data <- all_sur_data[grep(name_fold, all_sur_data$CODE), ][,-c(1)]
  
  gene_exp <- read.csv('~/projects/m6a/data/m6a.csv')
  
  gene_exp <- gene_exp[grep(name_fold, gene_exp$CODE), ][,-c(1)][grep('Tumor', gene_exp$Group), ][,-c(1)]
  gene_exp$SampleName <- gsub('\n', '', gene_exp$SampleName)
  
  erged_df <- merge(gene_exp, all_sur_data, by = "SampleName")
  erged_df <- erged_df[!duplicated(erged_df$SampleName),]
  erged_df$Time <- erged_df$Time / 365
  erged_df <- na.omit(erged_df)
  
  survival_time <- erged_df$Time
  survival_status <- erged_df$Status
  patient_id <- erged_df$SampleName
  erged_df <- subset(erged_df, select = -Status)
  erged_df <- subset(erged_df, select = -SampleName)
  erged_df <- subset(erged_df, select = -Time)
  erged_df = as.data.frame(lapply(erged_df, as.numeric))
  rownames(erged_df) <- patient_id
  
  if (!dir.exists(name_fold)){
    dir.create(name_fold)
  } else {
    print("Dir already exists!")
  }
  
  standardize <- function(x) {
    rowmean <- apply(x, 1, mean)
    rowsd <- apply(x, 1, sd)  
    rv <- sweep(x, 1, rowmean, "-")
    rv <- sweep(rv, 1, rowsd, "/")
    return(rv)
  }
  
  #########################################
  
  library(survival)
  library(glmnet)
  
  surv_obj <- Surv(event = survival_status, time = survival_time)
  
  fit <- glmnet(as.matrix(erged_df), surv_obj, family = "cox", maxit=1000)
  
  cvfit <- cv.glmnet(as.matrix(erged_df), surv_obj, family = "cox", maxit=1000)
  
  # Save and plot cross-validation fit results
  pdf(paste(name_fold, "/fit.pdf", sep=''), width = 5, height = 4)
  plot(cvfit)
  dev.off()
  
  # Save and plot lambda coefficients
  pdf(paste(name_fold, "/lambda.pdf", sep=''), width = 5, height = 4)
  plot(fit, xvar = "lambda")
  coef_paths <- coef(fit)
  genes <- colnames(as.matrix(erged_df))
  for (i in 1:length(genes)) {
    coef_values <- coef_paths[i, ]  
    lines(log(fit$lambda), coef_values, col = i + 1, lwd = 2)
    gene_label_x <- min(log(fit$lambda))
    gene_label_y <- ifelse(sum(coef_values) <= 0, min(coef_values), max(coef_values))
    text(gene_label_x + (abs(min(coef(fit))) + abs(max(coef(fit)))) / 500, gene_label_y + (abs(min(coef(fit))) + abs(max(coef(fit)))) / 500, labels = genes[i], col = i + 1, cex = 0.7)  # Add gene labels at the maximum coefficient points
  }
  abline(v = log(cvfit$lambda.min), col = "red4", lty = 2)
  legend("topright", legend = genes, col = 2:(length(genes) + 1), lwd = 1.5, cex = 0.25)
  dev.off()
  
  erged_df$Time <- survival_time
  erged_df$Status <- survival_status
  
  if(file.exists(paste(name_fold, '/coef.Rdata', sep='')))
  {
    load(paste(name_fold, '/coef.Rdata', sep=''))
    print('load....')
  } else {
    coef = coef(fit, s = cvfit$lambda.min)
    if(length(coef[which(coef != 0),]) <= 1)
    {
      file.rename(name_fold, paste('z_', name_fold, sep=''))
      print(paste(name_fold, 'coef = 0'))
      next
    }
    save(coef, file = paste(name_fold, '/coef.Rdata', sep=''))
    print('save....')
  }
  
  index = which(coef != 0)
  actCoef = coef[index]
  lassoGene = row.names(coef)[index]
  geneCoef = cbind(Gene = lassoGene, Coef = actCoef)
  write.table(geneCoef, file = paste(name_fold, "/geneCoef.txt", sep=''), sep = "\t", quote = F, row.names = F)
  
  trainFinalGeneExp = erged_df[, lassoGene]
  myFun = function(x) { crossprod(as.numeric(x), actCoef) }
  trainScore = apply(trainFinalGeneExp, 1, myFun)
  outCol = c("Time", "Status", lassoGene)
  risk = as.vector(ifelse(trainScore > median(trainScore), "high", "low"))
  num_out = standardize(as.data.frame(lapply(erged_df[, lassoGene], as.numeric)))
  rownames(num_out) = rownames(erged_df)
  outTab = cbind(erged_df[, c("Time", "Status")], num_out, riskScore = as.vector(trainScore), risk)
  write.table(cbind(id = rownames(outTab), outTab), file = paste(name_fold, "/Risk.txt", sep=''), sep = "\t", quote = F, row.names = F)
  
  ################################################
  
  library(timeROC)
  
  bioROC = function(inputFile = null, rocFile = null) {
    rt = read.table(inputFile, header = T, sep = "\t")
    ROC_rt = timeROC(T = rt$Time, delta = rt$Status,
                      marker = rt$riskScore, cause = 1,
                      weighting = 'aalen',
                      times = c(1, 2, 3), ROC = TRUE)
    pdf(file = rocFile, width = 5, height = 5)
    plot(ROC_rt, time = 1, col = 'green4', title = FALSE, lwd = 2)
    plot(ROC_rt, time = 2, col = 'blue4', add = TRUE, title = FALSE, lwd = 2)
    plot(ROC_rt, time = 3, col = 'red4', add = TRUE, title = FALSE, lwd = 2)
    legend('bottomright',
           c(paste0('AUC at 1 years: ', sprintf("%.03f", ROC_rt$AUC[1])),
             paste0('AUC at 2 years: ', sprintf("%.03f", ROC_rt$AUC[2])),
             paste0('AUC at 3 years: ', sprintf("%.03f", ROC_rt$AUC[3]))),
           col = c("green4", 'blue4', 'red4'), lwd = 2, bty = 'n')
    dev.off()
    return(max(c(ROC_rt$AUC[1], ROC_rt$AUC[2], ROC_rt$AUC[3])))
  }
  roc <- bioROC(inputFile = paste(name_fold, "/Risk.txt", sep=''), rocFile = paste(name_fold, "/ROC.pdf", sep=''))
  
  ######################################
  
  library(survminer)
  bioSurvival = function(inputFile = null, outFile = null) {
    rt = read.table(inputFile, header = T, sep = "\t")                  
    diff = survdiff(Surv(Time, Status) ~ risk, data = rt)
    pValue = 1 - pchisq(diff$chisq, df = 1)
    pValue = signif(pValue, 4)
    pValue = format(pValue, scientific = TRUE)
    fit <- survfit(Surv(Time, Status) ~ risk, data = rt)
    surPlot = ggsurvplot(fit, 
                         conf.int = TRUE,
                         data = rt,
                         pval = paste0("p=", pValue),
                         pval.size = 5,
                         legend.labs = c("High risk", "Low risk"),
                         legend.title = "Risk",
                         xlab = "Time(years)",
                         break.time.by = 1,
                         risk.table.title = "",
                         risk.table = F,
                         ggtheme = theme_classic())
    pdf(file = outFile, width = 5, height = 5)
    print(surPlot)
    dev.off()
    return(pValue)
  }
  pValue <- bioSurvival(inputFile = paste(name_fold, "/Risk.txt", sep=''), outFile = paste(name_fold, "/Survival.pdf", sep=''))
  final_list = rbind(final_list, data.frame(name_fold, ROC = roc, pValue = pValue))
}

write.table(final_list, file = 'final_list.csv', sep = ',', quote = FALSE, row.names = FALSE)
