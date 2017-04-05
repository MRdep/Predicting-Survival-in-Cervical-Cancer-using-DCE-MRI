# univariate ROC analysis

# function to make Kaplan-Meier plots. Run from batchKM.R.
library(grDevices)
library(survival)
library(rms)
library(colorspace)
library(ggplot2)
library(prodlim)
library(survivalROC)

# 1. Jobs need extract p values and plot as a function of split point.

# ROCanalysis <- function(param){

  computer <- 'lenovo'
  if(computer == 'hils'){computertag <- 'hilarywolfendale'}
  if(computer == 'work'){computertag <- 'Ben Dickie'}
  if(computer == 'lenovo'){computertag <- 'mbcx7bd2'}


  # read in tracer kinetic and clinical parameters
  model.params.ktrans <- read.csv('Ktrans_features.csv')
  model.params.Fp <- read.csv('Fp_features.csv')
  model.params.PS <- read.csv('PS_features.csv')
  model.params.vp <- read.csv('vp_features.csv')
  model.params.ve <- read.csv('ve_features.csv')

  # this
  clinical.params <- read.csv('ClinicalParameters_UpdatedFeb2017.csv')

  model.params <- cbind(model.params.ktrans['Ktrans.median'], model.params.Fp['Fp.median'],model.params.PS['FE.median'],model.params.vp['vp.median'],model.params.ve['ve.median'])
  colnames(model.params) <- c('Median.Ktrans.2CXM', 'Median.Fp.2CXM', 'Median.PS.2CXM', 'Median.vp.2CXM', 'Median.ve.2CXM')
  model.params <- data.frame(model.params)

  classFp2CXM <- vector()
  classFE2CXM <- vector()
  classve2CXM <- vector()
  classvp2CXM <- vector()
  classKtrans2CXM <- vector()

  classKtransTOFTS<- vector()
  classveTOFTS <- vector()

  classKtransETOFTS<- vector()
  classveETOFTS <- vector()

  classAge <- vector()
  classStage <- vector()
  classVol <- vector()
  classStage <- clinical.params$Stage_2level
  classFpandStage <- vector()
  classFpandHistology <- vector()
  classHistology <- vector()
  classNodes <- vector()
  classHistologyandStage <- vector()
  classTreatment <- vector()


  # calculate sample medians
  splitpoint.Fp.2CXM <- sort(model.params$Median.Fp.2CXM)
  splitpoint.PS.2CXM <- sort(model.params$Median.PS.2CXM)
  splitpoint.Ktrans.2CXM <- sort(model.params$Median.Ktrans.2CXM)
  splitpoint.vp.2CXM <- sort(model.params$Median.vp.2CXM)
  splitpoint.ve.2CXM <- sort(model.params$Median.ve.2CXM)

  nsplit <- length(splitpoint.Fp.2CXM)

  splitpoint.Age <- sort(clinical.params$Age)
  splitpoint.Vol <- sort(clinical.params$Vol)

  len <- length(model.params$Median.Fp)

  DFStime <- clinical.params$DFS
  status <- clinical.params$DFSstatus

  pval_Fp <- rep(1, 36)
  pval_PS <- rep(1, 36)
  pval_vp <- rep(1, 36)
  pval_ve <- rep(1, 36)


  # Compute classifications
  # 2cxm parameters
  for(j in 4:(nsplit - 5)){
  for (i in 1:len){
    if(model.params$Median.Fp.2CXM[i] < splitpoint.Fp.2CXM[j]){
      classFp2CXM[i] = 'Fp < splitpoint'
    }
    else{
      classFp2CXM[i] = 'Fp > splitpoint'
    }
  }

  for (i in 1:len){
    if(model.params$Median.PS.2CXM[i] < splitpoint.PS.2CXM[j]){
      classFE2CXM[i] = 'PS < splitpoint'
    }
    else{
      classFE2CXM[i] = 'PS > splitpoint'

    }
  }
  set.seed(123)
  model.params$Median.vp.2CXM <- model.params$Median.vp.2CXM + rnorm(length( model.params$Median.vp.2CXM), 0, 0.0001)
  splitpoint.vp.2CXM <- sort(model.params$Median.vp.2CXM)
  for (i in 1:len){
    if(model.params$Median.vp.2CXM[i] < splitpoint.vp.2CXM[j]){
      classvp2CXM[i] = 'vp < splitpoint'
    }
    else{
      classvp2CXM[i] = 'vp > splitpoint'
    }
  }


  for (i in 1:len){
    if(model.params$Median.ve.2CXM[i] < splitpoint.ve.2CXM[j]){
      classve2CXM[i] = 've < splitpoint'
    }
    else{
      classve2CXM[i] = 've > splitpoint'

    }
  }


  for (i in 1:len){
    if(model.params$Median.Ktrans.2CXM[i] < splitpoint.Ktrans.2CXM[j]){
      classKtrans2CXM[i] = 'Ktrans < splitpoint'
    }
    else{
      classKtrans2CXM[i] = 'Ktrans > splitpoint'
    }
  }
  # pairs(classKtrans2CXM, classFp2CXM, classFE2CXM)

  for (i in 1:len){
    if(clinical.params$Age[i] < splitpoint.Age){
      classAge[i] = 'Age < splitpoint'
    }
    else{
      classAge[i] = 'Age > splitpoint'
    }
  }


  for (i in 1:len){
    if(clinical.params$Vol[i] < splitpoint.Vol){
      classVol[i] = 'Volume < splitpoint'
    }
    else{
      classVol[i] = 'Volume > splitpoint'
    }
  }
  # opposite trend to TK parameters. low is good
  for (i in 1:len){

    if(clinical.params$Stage_2level[i] == 1){

      classStage[i] = "1-2"
    }
    else{
      classStage[i] = "3-4"
    }
  }

  for (i in 1:len){

    if(clinical.params$Histology[i] == 'SCC'){

      classHistology[i] = "SCC"
    }
    else{
      classHistology[i] = "Other"
    }
  }

  for (i in 1:len){
    if(clinical.params$Nodes[i] == '0'){

      classNodes[i] = "Negative"
    }
    else{
      classNodes[i] = "Positive"
    }
  }

  for (i in 1:len){
    if(clinical.params$Treatment[i] == '1'){
      classTreatment[i] = 1
    }
    else{
      classTreatment[i] = 0
    }
  }



    # To univariate cox to compute HR and p value
    cox1 <- coxph(Surv(DFStime,status) ~ classFp2CXM, method = 'breslow')
    sum <- summary(cox1)
    # str(sum$waldtest)
    pval_Fp[j] <- sum$waldtest[3]

    # To univariate cox to compute HR and p value
    cox2 <- coxph(Surv(DFStime,status) ~ classFE2CXM, method = 'breslow')
    sum <- summary(cox2)
    # str(sum$waldtest)
    pval_PS[j] <- sum$waldtest[3]



    # To univariate cox to compute HR and p value
    cox3 <- coxph(Surv(DFStime,status) ~ classvp2CXM, method = 'breslow')
    sum <- summary(cox3)
    # str(sum$waldtest)
    pval_vp[j] <- sum$waldtest[3]



    # To univariate cox to compute HR and p value
    cox4 <- coxph(Surv(DFStime,status) ~ classve2CXM, method = 'breslow')
    sum <- summary(cox4)
    # str(sum$waldtest)
    pval_ve[j] <- sum$waldtest[3]




  }

  par(mfrow = c(2,2), mar = c(4,4,2,1), cex = 1.0)
  plot((1:36/36), runmed(pval_Fp, 2), ylim = c(0, 0.5), xlab = 'Group Split Percentile', ylab = 'P-value', type = 'line', main = expression(italic(F)[p]))
  plot((1:36/36), runmed(pval_PS, 2), ylim = c(0, 0.5), xlab = 'Group Split Percentile', ylab = 'P-value', type = 'line', main = expression(italic(PS)))
  plot((1:36/36), runmed(pval_vp, 2), ylim = c(0, 0.5), xlab = 'Group Split Percentile', ylab = 'P-value', type = 'line', main = expression(italic(v)[p]))
  plot((1:36/36), runmed(pval_ve, 2), ylim = c(0, 0.5), xlab = 'Group Split Percentile', ylab = 'P-value', type = 'line', main = expression(italic(v)[e]))

# }

  # ROC survival
  # FP anf TP assume that high marker values are bad! For TKM parameters, FP = TP and TP = FP. 

  # TK parameters
  png( "SupplementaryFigure2.png", width = 60, height = 22,8, units = 'in', res = 400)
  par(mfrow=c(2,4), cex.lab = 2.5,cex.axis = 2.5, mai = c(2,2,2,2), cex.main = 2.5, cex = 2.5)
  predict.times <- sort(unique(DFStime))
  FP <- array(NA, dim = c((length(DFStime) + 1), (c(length(predict.times) + 1))))
  TP <- array(NA, dim = c((length(predict.times) + 1), (c(length(predict.times) + 1))))
  AUC <- vector()
  Y <- array(NA, dim = c((length(predict.times) + 1), (7)))
  for(j in 1:7){
  if(j == 1){
  marker <- model.params$Median.Fp.2CXM
  param <- expression(italic(F[p]))
  splits <- splitpoint.Fp.2CXM}
    
    if(j == 2){marker <- model.params$Median.PS.2CXM
               param <-expression(italic(PS))
               splits <- splitpoint.PS.2CXM}
    if(j == 3){marker <- model.params$Median.vp.2CXM
               param <- expression(italic(v[p]))
               splits <- splitpoint.vp.2CXM}
    if(j == 4){marker <- model.params$Median.ve.2CXM
               param <- expression(italic(v[e]))
               splits <- splitpoint.ve.2CXM}
    if(j == 5){marker <- model.params$Median.Ktrans.2CXM
               param <- expression(italic(K)^{trans})
               splits <- splitpoint.Ktrans.2CXM}
  if(j == 6){marker <- clinical.params$Age
             param <-'Age'
             splits <- splitpoint.Age}
  if(j == 7){marker <- clinical.params$Vol
             param <-'Tumour volume'
             splits <- splitpoint.Vol}
  
  for(i in 1:length(predict.times)){
  ROC <- survivalROC(DFStime, status, marker, predict.time = predict.times[i], method = "NNE", lambda = NULL, span = 0.0105, window ="symmetric")

  AUC[i] <- 1 - ROC$AUC
  if(j <= 5){
  FP[i,] <- ROC$TP
  TP[i,]<- ROC$FP
  AUC[i] <- 1 - ROC$AUC
  }else{
    FP[i,] <- ROC$FP
    TP[i,]<- ROC$TP 
    AUC[i] <- ROC$AUC
  }
  
  Y[i,j] <- which.max(TP[i,] - FP[i,])
  }
  print('-------------------')
  print(param)
  print('------------------')
  print('AUC')
  print(AUC[18])
  print('Youden threshold (abs)')
  print(splits[Y[i,j]])
  print('Youden threshold (percentile)')
  if(Y[18,j] > 26){Y[18,j] <- 26}
  print((Y[18,j] - 1)/36)
  print('True Positive Rate (Sensitivity)')
  print(TP[18,18])
  print('False Positive Rate (1 - Specificity)')
  print(FP[18,18])
#   print('Number with event')
#   print('18')
#   print('Number without event')
#   print('18')
  print('Number of True Positives')
  print(TP[18,18]*18)
print('Number of False Positives')
print(FP[18,18]*18)
print('Number of True Negatives')
print(18 - 18*FP[18,18])
  print('Number of False Negatives')
  print(18 - 18*TP[18,18])


  
  
  plot(FP[18,1:(length(FP[1,]) -2)], TP[18,1:(length(TP[1,]) - 2)], xlab = '1 - Specificity', ylab = 'Sensitivity', type = 'l', main = param, lwd = 2, xlim = c(0, 1))
  text(0.18, 0.9, paste('AUC = ', format(signif(AUC[18], 2))), cex = 2.5)
  abline(0,1, lwd = 2)
  }
dev.off()
  
