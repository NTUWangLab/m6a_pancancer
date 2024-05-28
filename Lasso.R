

setwd('~/projects/m6a/Lasso')

Tumor <- c('ACC','BLCA','BRCA','CESC','COAD', 'CHOL', 'ESCA', 'GBM', 'HNSC', 'KIRP', 'KIRC', 'KICH', 'LGG', 'LUAD', 'LUSC', 'LIHC', 'LAML', 'OV', 'PRAD', 'PAAD', 'PCPG', 'READ', 'STAD', 'SKCM', 'THCA', 'TGCT', 'UCEC', 'UCS')


final_list=data.frame()

for(name_fold in Tumor)
{
  all_sur_data <- read.csv('~/history_lasso/lassoPancancer/survival.csv',header=T,check.names=F)
  all_sur_data <- all_sur_data[grep(name_fold,all_sur_data$CODE),][,-c(1)]
  
  gene_exp <- read.csv('~/projects/m6a/data/m6a.csv')
  
  gene_exp <- gene_exp[grep(name_fold,gene_exp$CODE),][,-c(1)][grep('Tumor',gene_exp$Group),][,-c(1)]
  gene_exp$SampleName <- gsub('\n','',gene_exp$SampleName)
  
  erged_df <- merge(gene_exp, all_sur_data, by = "SampleName")
  erged_df <- erged_df[!duplicated(erged_df$SampleName),]
  erged_df$Time <- erged_df$Time/365
  erged_df <- na.omit(erged_df)
  
  survival_time <- erged_df$Time
  survival_status <- erged_df$Status
  patient_id <- erged_df$SampleName
  erged_df <-subset(erged_df, select = -Status)
  erged_df <-subset(erged_df, select = -SampleName)
  erged_df <-subset(erged_df, select = -Time)
  erged_df=as.data.frame(lapply(erged_df,as.numeric))
  rownames(erged_df) <- patient_id
  
  if (!dir.exists(name_fold)){
    dir.create(name_fold)
  } else {
    print("Dir already exists!")
  }
  
  
  standardize <- function(x) {
    rowmean <- apply(x, 1, mean)
    rowsd <- apply(x, 1, sd)  
    rv <- sweep(x, 1, rowmean,"-")  #表达量-均值
    rv <- sweep(rv, 1, rowsd, "/")  #再除以标准差
    return(rv)
  }
  
  #########################################
  
  library(survival)
  library(glmnet)
  
  surv_obj <- Surv(event = survival_status,time = survival_time)
  
  fit <- glmnet(as.matrix(erged_df), surv_obj, family = "cox",maxit=1000)
  
  cvfit <- cv.glmnet(as.matrix(erged_df), surv_obj, family = "cox",maxit=1000)
  
  
  pdf(paste(name_fold,"/fit.pdf",sep=''),width = 5,height = 4)
  plot(cvfit)
  dev.off()
  
  pdf(paste(name_fold,"/lambda.pdf",sep=''),width = 5,height = 4)
  plot(fit, xvar = "lambda")
  coef_paths <- coef(fit)
  genes <- colnames(as.matrix(erged_df))
  for (i in 1:length(genes)) {
    coef_values <- coef_paths[i, ]  # 获取第i个基因的系数路径
    lines(log(fit$lambda), coef_values, col = i + 1, lwd = 2)  # 使用不同的颜色
    gene_label_x <- min(log(fit$lambda))
    gene_label_y <- ifelse(sum(coef_values) <= 0, min(coef_values), max(coef_values))
    text(gene_label_x+(abs(min(coef(fit))) + abs(max(coef(fit))))/500, gene_label_y+(abs(min(coef(fit))) + abs(max(coef(fit))))/500, labels = genes[i], col = i + 1,cex=0.7)  # 在最大系数点上添加基因标签
  }
  abline(v = log(cvfit$lambda.min), col = "red4", lty = 2)
  legend("topright", legend = genes, col = 2:(length(genes) + 1), lwd = 1.5,cex=0.25)
  dev.off()
  
  erged_df$Time <- survival_time
  erged_df$Status <- survival_status
  
  if(file.exists(paste(name_fold,'/coef.Rdata',sep='')))
  {
    load(paste(name_fold,'/coef.Rdata',sep=''))
    print('load....')
  } else {
    coef=coef(fit, s = cvfit$lambda.min)
    if(length(coef[which(coef != 0),])<=1)
    {
      file.rename(name_fold,paste('z_',name_fold,sep=''))
      print(paste(name_fold,'coef = 0'))
      next
    }
    save(coef,file = paste(name_fold,'/coef.Rdata',sep=''))
    print('save....')
  }
  
  index=which(coef != 0)
  actCoef=coef[index]
  lassoGene=row.names(coef)[index]
  geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
  write.table(geneCoef,file=paste(name_fold,"/geneCoef.txt",sep=''),sep="\t",quote=F,row.names=F)
  
  trainFinalGeneExp=erged_df[,lassoGene]
  myFun=function(x){crossprod(as.numeric(x),actCoef)}
  trainScore=apply(trainFinalGeneExp,1,myFun)
  outCol=c("Time","Status",lassoGene)
  risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
  num_out <- standardize(as.data.frame(lapply(erged_df[,lassoGene],as.numeric)))
  rownames(num_out) <- rownames(erged_df)
  outTab=cbind(erged_df[,c("Time","Status")],num_out,riskScore=as.vector(trainScore),risk)
  write.table(cbind(id=rownames(outTab),outTab),file=paste(name_fold,"/Risk.txt",sep=''),sep="\t",quote=F,row.names=F)
  

  
  ################################################
  
  library(timeROC)
  
  bioROC=function(inputFile=null,rocFile=null){
    
    rt=read.table(inputFile,header=T,sep="\t")
    ROC_rt=timeROC(T=rt$Time,delta=rt$Status,
                   marker=rt$riskScore,cause=1,
                   weighting='aalen',
                   times=c(1,2,3),ROC=TRUE)
    pdf(file=rocFile,width=5,height=5)
    plot(ROC_rt,time=1,col='green4',title=FALSE,lwd=2)
    plot(ROC_rt,time=2,col='blue4',add=TRUE,title=FALSE,lwd=2)
    plot(ROC_rt,time=3,col='red4',add=TRUE,title=FALSE,lwd=2)
    legend('bottomright',
           c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
             paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
             paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
           col=c("green4",'blue4','red4'),lwd=2,bty = 'n')
    dev.off()
    return(max(c(ROC_rt$AUC[1],ROC_rt$AUC[2],ROC_rt$AUC[3])))
  }
  roc <- bioROC(inputFile=paste(name_fold,"/Risk.txt",sep=''),rocFile=paste(name_fold,"/ROC.pdf",sep=''))
  
  ######################################
  
  library(survminer)
  bioSurvival=function(inputFile=null,outFile=null){
    rt=read.table(inputFile,header=T,sep="\t")                  
    diff=survdiff(Surv(Time, Status) ~ risk,data = rt)
    pValue=1-pchisq(diff$chisq,df=1)
    pValue=signif(pValue,4)
    pValue=format(pValue, scientific = TRUE)
    fit <- survfit(Surv(Time, Status) ~ risk, data = rt)
    surPlot=ggsurvplot(fit, 
                       conf.int = TRUE,
                       data=rt,
                       pval=paste0("p=",pValue),
                       pval.size=5,
                       legend.labs=c("High risk", "Low risk"),
                       legend.title="Risk",
                       xlab="Time(years)",
                       break.time.by = 1,
                       risk.table.title="",
                       risk.table=F,
                       risk.table.height=.25)
    pdf(file=outFile,onefile = FALSE,width = 5,height =4.5)
    print(surPlot)
    dev.off()
    return(pValue)
  }
  p_value <- bioSurvival(inputFile=paste(name_fold,"/Risk.txt",sep=''),outFile=paste(name_fold,"/survival.pdf",sep=''))
  
  #######################################################################3
  
  library(pheatmap)
  bioRiskPlot=function(inputFile=null,riskScoreFile=null,survStatFile=null){
    rt=read.table(inputFile,sep="\t",header=T,row.names=1,check.names=F)  
    rt=rt[order(rt$riskScore),]    
    riskClass=rt[,"risk"]
    lowLength=length(riskClass[riskClass=="low"])
    highLength=length(riskClass[riskClass=="high"])
    lowMax=max(rt$riskScore[riskClass=="low"])
    line=rt[,"riskScore"]
    line[line>10]=10
    pdf(file=riskScoreFile,width = 8,height = 6)
    plot(line, type="p", pch=20,
         xlab="Patients (increasing risk socre)", ylab="Risk score",
         col=c(rep("blue4",lowLength),rep("red",highLength)) )
    abline(h=lowMax,v=lowLength,lty=2)
    legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("red","blue4"),cex=1.2)
    dev.off()
    color=as.vector(rt$Status)
    color[color==1]="red"
    color[color==0]="blue4"
    pdf(file=survStatFile,width = 8,height = 6)
    plot(rt$Time, pch=19,
         xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
         col=color)
    legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","blue4"),cex=1.2)
    abline(v=lowLength,lty=2)
    dev.off()
  }
  
  bioRiskPlot(inputFile=paste(name_fold,"/Risk.txt",sep=''),riskScoreFile=paste(name_fold,"/riskScore.pdf",sep=''),survStatFile=paste(name_fold,"/survStat.pdf",sep=''))
  
  ##########################################################
  
  risk.out <- read.table(paste(name_fold,'/Risk.txt',sep=''),row.names = 1,header=T)
  risk.out <- risk.out[order(risk.out$risk),]
  risk.out <-subset(risk.out, select = -c(Time,Status,riskScore))
  Group <- risk.out$risk
  risk.out <-subset(risk.out, select = -risk)
  Group <- data.frame(Risk = Group)
  rownames(Group) = rownames(risk.out)
  library(pheatmap)
  pdf(paste(name_fold,'/phep.pdf',sep=''),width=6,height=4)
  mycolors <- colorRampPalette(c("darkgreen","white", "red4"))(15)
  tmp=as.data.frame(lapply(risk.out,as.numeric))
  tmp <- as.matrix(t(tmp))
  colnames(tmp) <- rownames(Group)
  pheatmap(scale(tmp),annotation_col = Group,color = mycolors,cluster_rows=F,show_colnames=F,cluster_cols=F,cellheight = 10,cellwidth = 200/length(colnames(tmp)),border=F)
  dev.off()
  print(paste('p:',p_value,' roc:',roc,sep=''))
  
  if(is.na(roc))
  {
    file.rename(name_fold,paste('z_',name_fold,sep=''))
    print(paste(name_fold,'p o roc not satisify'))
    next
  } else {
    if(as.numeric(p_value) > 0.05||as.numeric(roc)<0.7)
    {
      file.rename(name_fold,paste('z_',name_fold,sep=''))
      print(paste(name_fold,'p o roc not satisify'))
      next
    } else {
      if(ncol(final_list)==0)
      {
        final_list = data.frame(SampleName=rownames(outTab),Risk=outTab$risk,CODE=rep(name_fold,length(outTab$risk)))
      }
      else
      {
        final_list = rbind(final_list,data.frame(SampleName=rownames(outTab),Risk=outTab$risk,CODE=rep(name_fold,length(outTab$risk))))
      }
    }
  }
}
write.csv(final_list,file='risk.csv')
rm(list=ls())
system('ls|grep -E  "z_[A-Z]+"|xargs rm -r')
