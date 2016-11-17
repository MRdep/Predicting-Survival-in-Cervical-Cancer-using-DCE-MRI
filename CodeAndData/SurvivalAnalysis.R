library(randomForestSRC)
library(survival)
library(rms)
library(prodlim)
library(boot)
library(survcomp)
library(datasets)
library(MASS)

#compile source files needed to run the script 
source('OutcomePartialPlots.R')
source('DicotamiseOnMedian.R')

# set cwd to source file location 
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# # #first make Figures 2 and Supp 2 and write out as pdfs
# # print('----------------------------------------------------------------------------------') 
# # print('Writing Figures 2 and Supp 2 to disk.')
# # print('----------------------------------------------------------------------------------') 

# source('KMplotter.R')
# source('batchKM_2CXM_only.R')

# read in tracer kinetic and clinical parameters
model.params.ktrans <- read.csv('Ktrans_features.csv')
model.params.Fp <- read.csv('Fp_features.csv')
model.params.PS <- read.csv('PS_features.csv')
model.params.vp <- read.csv('vp_features.csv')
model.params.ve <- read.csv('ve_features.csv')

# this
clinical.params <- read.csv('ClinicalParameters_UpdatedJuly2016.csv')

clinical.params <- clinical.params
model.params <- cbind(model.params.ktrans['Ktrans.median'], model.params.Fp['Fp.median'],model.params.PS['FE.median'],model.params.vp['vp.median'],model.params.ve['ve.median'])
colnames(model.params) <- c('Median.Ktrans.2CXM', 'Median.Fp.2CXM', 'Median.PS.2CXM', 'Median.vp.2CXM', 'Median.ve.2CXM')
model.params <- data.frame(model.params)
len <- length(model.params$Median.Fp.2CXM)

# dicotmise variables.
classFp2CXM <- DicotamiseOnMedian(model.params$Median.Fp.2CXM)
classFE2CXM <- DicotamiseOnMedian(model.params$Median.PS.2CXM)
classve2CXM <- DicotamiseOnMedian(model.params$Median.ve.2CXM)
classvp2CXM <- DicotamiseOnMedian(model.params$Median.vp.2CXM)
classKtrans2CXM <- DicotamiseOnMedian(model.params$Median.Ktrans.2CXM)

classAge <- DicotamiseOnMedian(clinical.params$Age)
classVol <- DicotamiseOnMedian(clinical.params$Vol)

classStage <- clinical.params$Stage_2level
classNodes <- clinical.params$Nodes
classHistology<- vector()
classTreatment <- vector()

# changed order
for (i in 1:len){
  if(clinical.params$Histology[i] == 'SCC'){
    classHistology[i] = 0
  } 
  else{
    classHistology[i] = 1
  }
}

# changed order
for (i in 1:len){
  if(clinical.params$Treatment[i] == 1){
    classTreatment[i] = 0
  } 
  if(clinical.params$Treatment[i] == 2){
    classTreatment[i] = 1
  }
}

# extract DFS times. 
DFStime <- clinical.params$DFS
# extract censoring status.
status <- clinical.params$DFSstatus

#make a data frame containing all clinicopathologic and tracer kinetic variables. 
cervixDicot <- data.frame(DFStime, status,classKtrans2CXM, classFp2CXM,classFE2CXM,classve2CXM,classvp2CXM,classAge,classVol,classStage,classNodes,classHistology, classTreatment)
colnames(cervixDicot) <- c('DFStime', 'status','XM.Ktrans', 'XM.Fp', 'XM.FE', 'XM.ve', 'XM.vp', 'Age', 'Vol', 'Stage','Nodes','Histology','Treatment')
cervixDicot_CRTonly <- subset(cervixDicot, classTreatment == 1)

#make a data frame containing only clinicopathologic variables. 
cervixDicot.clinical <- cbind(cervixDicot[,1:2], cervixDicot[,8:13])

# compute Cox model univariate HR
print(coxph(Surv(DFStime, status) ~ XM.Ktrans, data = cervixDicot))
print(coxph(Surv(DFStime, status) ~ XM.Fp, data = cervixDicot))
print(coxph(Surv(DFStime, status) ~ XM.FE, data = cervixDicot))
print(coxph(Surv(DFStime, status) ~ XM.vp, data = cervixDicot))
print(coxph(Surv(DFStime, status) ~ XM.ve, data = cervixDicot))
print(coxph(Surv(DFStime, status) ~ Vol, data = cervixDicot))
print(coxph(Surv(DFStime, status) ~ Stage, data = cervixDicot))
print(coxph(Surv(DFStime, status) ~ Nodes, data = cervixDicot))
print(coxph(Surv(DFStime, status) ~ Histology, data = cervixDicot))
print(coxph(Surv(DFStime, status) ~ Treatment, data = cervixDicot))

#-------------------------------------------------------------------
#Fit RSF model with all clinicopathologic and tracer kinetic variables 100 DFStimes in a bootstrap experiment
#-------------------------------------------------------------------
cervixData <- cervixDicot
cervixData.clinical <- cervixDicot.clinical

#define a variable to hold the number of observations 
n_pats <- nrow(cervixData)
#seed RNG. 
set.seed(123)
# set number of trees.
ntree = 1000
nodesize = 3
nsplit = 2
#fit models once first. 
set.seed(123)
rsf.all.vimp <- rfsrc(Surv(DFStime, status) ~., ntree=ntree, data = cervixData, nodesize = nodesize, bootstrap = c("by.root"),importance = T, nsplit = nsplit)

set.seed(123)
rsf.clinical <- rfsrc(Surv(DFStime, status) ~ ., ntree=ntree, data = cervixData.clinical, nodesize = nodesize, bootstrap = c("by.root"),importance = T, nsplit = nsplit)

oob.predictions <- rsf.clinical$predicted.oob
index.predictions <- order(oob.predictions)
predictor <- vector()
predictor[index.predictions[1:18]] <- 'Short DFS'
predictor[index.predictions[19:36]] <- 'Long DFS'

DFS.actual <- survfit(Surv(DFStime, status)~ predictor, data = cervixData)

oob.survival.functions.pred <- rsf.clinical$survival.oob
oob.survival.functions.longDFS <- oob.survival.functions.pred[index.predictions[1:18],]
oob.survival.functions.shortDFS <- oob.survival.functions.pred[index.predictions[19:36],]

plot(DFS.actual)
lines(rsf.clinical$time.interest, colMeans(oob.survival.functions.shortDFS))
lines(rsf.clinical$time.interest, colMeans(oob.survival.functions.longDFS))

importances.all.vimp <- rsf.all.vimp$importance

set.seed(123)
rsf.top6 <- rfsrc(Surv(DFStime, status) ~ XM.Fp + Age + Treatment + Histology + XM.ve + Nodes, ntree=ntree, data = cervixData, nodesize = nodesize, bootstrap = c("by.root"),importance = T, nsplit = nsplit)

oob.predictions <- rsf.top6$predicted.oob
index.predictions <- order(oob.predictions)
predictor <- vector()
predictor[index.predictions[1:18]] <- 'Short DFS'
predictor[index.predictions[19:36]] <- 'Long DFS'

DFS.actual <- survfit(Surv(DFStime, status) ~ predictor, data = cervixData)

oob.survival.functions.pred <- rsf.top6$survival.oob
oob.survival.functions.longDFS <- oob.survival.functions.pred[index.predictions[1:18],]
oob.survival.functions.shortDFS <- oob.survival.functions.pred[index.predictions[19:36],]

survival.functions.pred <- rsf.top6$survival
survival.functions.longDFS <- survival.functions.pred[index.predictions[1:18],]
survival.functions.shortDFS <- survival.functions.pred[index.predictions[19:36],]

importances.all.vimp <- rsf.all.vimp$importance

n_repeats <- 100

importances.vimp <- matrix(nrow = n_repeats, ncol = (ncol(cervixData) - 2))
importances.clinical <- matrix(nrow = n_repeats, ncol = (ncol(cervixData.clinical) - 2))

colnames(importances.vimp) <- names(importances.all.vimp)
colnames(importances.clinical) <- names(rsf.clinical$importance)

print('Running bootstrap experiment.') 

for(i in 1:n_repeats){
  remainder <- i%%100
  if(remainder == 0){print(paste('MC repeat = ', i, '/1000', sep = ""))}
  # take a bootstrap sample    
  train.sample <- sample(1:n_pats, n_pats, replace = T)  
  data.train <- cervixData[train.sample, ]
  data.train.clinical <- cervixData.clinical[train.sample, ]
  # fit models   
  rsf.all.vimp <- rfsrc(Surv(DFStime, status) ~., ntree=ntree, data = data.train, nodesize = nodesize, bootstrap = c("by.root"),importance = T, nsplit = nsplit)  
  rsf.clinical <- rfsrc(Surv(DFStime, status)~ ., ntree = ntree, data = data.train.clinical, nodesize = nodesize, bootstrap = c("by.root"),importance = T, nsplit = nsplit)
  importances.vimp[i,] <- rsf.all.vimp$importance
  importances.clinical[i,] <- rsf.clinical$importance
}

print('----------------------------------------------------------------------------------') 

bp <- boxplot(importances.vimp, las =2, notch = T, plot = F)
print(bp$stats)
medians <- bp$stats[3,]

bp.clinical <- boxplot(importances.clinical, las =2, notch = T, plot = F)
medians.clinical<- bp.clinical$stats[3,]

print('Computing point estimates and 95% CI on VIMP.') 

bootThetaQuantile <- function(x,i) {
  quantile(x[i], probs=.5)
}

nboot = 100

bootstrapped.medians <- matrix(nrow = nboot, ncol = length(rsf.all.vimp$importance))
bootstrapped.medians.clinical <- matrix(nrow = nboot, ncol = length(rsf.clinical$importance))

# set the signficance criterion for the CI. Bonferroni corrected for the total number of inferences made. Take 2 to subtract out DFStime and status from cervixData. 
criterion = 0.05/(ncol(cervixData) - 2)

lower <- vector()
upper <- vector()

lower.clinical <- vector()
upper.clinical <- vector()

set.seed(123)
for(i in 1:length(importances.all.vimp)){
  boot.out <- boot(data=importances.vimp[,i], statistic=bootThetaQuantile,R=nboot)
  bootstrapped.medians[,i] <- boot.out$t
  ci <- boot.ci(boot.out,conf = 1- criterion, type = 'bca')
  lower[i] <- ci$bca[4] 
  upper[i] <- ci$bca[5] 
}

set.seed(123)
for(i in 1:length(rsf.clinical$importance)){
  boot.out.clinical <- boot(data=importances.clinical[,i], statistic=bootThetaQuantile,R=nboot)
  bootstrapped.medians.clinical[,i] <- boot.out.clinical$t
  ci <- boot.ci(boot.out.clinical,conf = 1- criterion, type = 'bca')
  lower.clinical[i] <- ci$bca[4] 
  upper.clinical[i] <- ci$bca[5] 
}

#output VIMP statistics into a table.
vimp <- rbind(lower, medians, upper)
colnames(vimp) <- names(rsf.all.vimp$importance)

clinical.vimp <- rbind(lower.clinical, medians.clinical, upper.clinical)
colnames(clinical.vimp) <- names(rsf.clinical$importance)

print('----------------------------------------------------------------------------------') 
print('Point estimates and 95% CI on VIMP for model with clinicopathologic variables only.') 
print(clinical.vimp)
print('----------------------------------------------------------------------------------') 
print('Point estimates and 95% CI on VIMP for model with clinicopathologic and tracer kinetic variables.') 
print(vimp)

# write out csv files of VIMP. 
write.csv(vimp, file = 'vimp.csv')
write.csv(clinical.vimp, file = 'clinical_vimp.csv')
print('VIMP files written to disk') 

#make VIMP plot (Figure 3). 
x.labels <- c(expression('2CXM K'^{trans}),
              expression('2CXM F'[p]),
              expression('2CXM PS'),
              expression('2CXM v'[e]),
              expression('2CXM v'[p]),
              expression('Age'),
              expression('Volume'),
              expression('T stage'),
              expression('Nodal Status'),
              expression('Histology'),
              expression('Treatment'))

# set.seed(1234)
# #Make Figure 5. Point estimates and 95% CI on VIMP.
# pdf('Figure5.pdf', width = 7, height = 5)
# par(mai = c(3,1,0.5,0.5), cex.lab = 1.3, cex.axis = 1.3)
# # col = c('red','black','red', 'red', 'red', 'black', 'red', 'red','red', 'black','red', 'red', 'black','red', 'black', 'black')
#  indices <- sort.int(medians, decreasing = T,index.return = T)
# # col <- c('red','chartreuse4','red', 'chartreuse4', 'red', 'chartreuse4', 'red', 'red','chartreuse4', 'chartreuse4','chartreuse4', 'chartreuse4')
# 
# # col <- col[indices$ix]
# plot(medians[indices$ix], xaxt = 'n', yaxt = 'n',xlab = "", ylab = "", pch = 19, ylim = c(-0.001, 0.08), cex = 0.8)
# axis(1, at = 1:length(importances.all.vimp), labels = x.labels[indices$ix], las = 2)
# axis(2, at  = c(0, 0.04, 0.08), labels = c(0, 0.04, 0.08))
# # abline(h  =0, lty =2, lwd = 2)
# abline(h = 0.021, lty = 2, lwd = 2)
# mtext(side = 2, text = "VIMP", line = 2.5, cex = 1.3)
# lower2 <- lower[indices$ix]
# upper2 <- upper[indices$ix]
# for(i in 1:length(importances.all.vimp)){
#   segments(i, lower2[i], x1 = i, y1 = upper2[i], lwd = 2)
# }
# 
# 
# dev.off()
# print('Writing Figure 5 to disk.')
# print('----------------------------------------------------------------------------------')

print('Starting Leave-one-out experiment to assess model accuracy.') 
#compute LOOCV c-indices 
rsf.final <- rfsrc(Surv(DFStime,status) ~ Stage +  Histology + Treatment +  Age + Vol + Nodes + XM.Fp, ntree=ntree,data = cervixData, nodesize = nodesize, bootstrap = c("by.root"),importance = T, nsplit = nsplit)
pred <- vector()
pred.clinical <- vector()
yearProb <- vector()
yearProb.clinical <- vector()

for(i in 1:nrow(cervixData)){
  if(i == 18){print('Half way...')}
  data.train <- cervixData[-i, ]
  data.test <- cervixData[i,]
  data.train.clinical = cervixData.clinical[-i,]
  data.test.clinical = cervixData.clinical[i,]
  rsf.all.vimp <- rfsrc(Surv(DFStime,status) ~ XM.Fp + XM.ve + Nodes + Histology + Treatment + Age, ntree=ntree,data = data.train, nsplit = nsplit)
  rsf.clinical <- rfsrc(Surv(DFStime, status) ~., ntree = ntree, data = data.train.clinical, nsplit = nsplit)
  # pull out the prediction for the text data   
  pred[i] <- predict(rsf.all.vimp, newdata = cervixData)$predicted[i]
  pred.clinical[i] <- predict(rsf.clinical, newdata = cervixData.clinical)$predicted[i]
  yearProb[i] <- predict(rsf.all.vimp, newdata = cervixData)$survival[i,15] 
  yearProb.clinical[i] <- predict(rsf.clinical, newdata = cervixData.clinical)$survival[i,15]
}

c.index.struct <- concordance.index(pred,cervixData$DFStime, cervixData$status, method = 'noether')
c.index.struct.clinical <- concordance.index(pred.clinical,cervixData$DFStime, cervixData$status, method = 'noether')
c.index.final <- c.index.struct$c.index
c.index.clinical <- c.index.struct.clinical$c.index
cindex.comp(c.index.struct,c.index.struct.clinical)

# make 5-year DFS risk prediction partial plots (Figure 5) 
rsf.clinical <- rfsrc(Surv(DFStime,status) ~ Stage +  Histology + Treatment +  Vol + Age + Nodes, ntree=ntree, data = cervixData,nodesize = nodesize, bootstrap = c("by.root"),importance = T, nsplit = nsplit)
rsf.all.vimp <- rfsrc(Surv(DFStime,status) ~ Histology + Treatment + Age + Nodes + XM.Fp + XM.ve, ntree=ntree,data = cervixData,nodesize = nodesize, bootstrap = c("by.root"),importance = T, nsplit = nsplit)

pp.data <- plot.variable(rsf.all.vimp, partial = T,surv.type = "surv",time = DFStime[23], show.plots = F)
pp.data.clinical <- plot.variable(rsf.clinical, partial = T,surv.type = "surv",time = DFStime[23], show.plots = F)

LowFpDFS <- vector()
HighFpDFS <- vector()
LowStageDFS <- vector()
HighStageDFS <- vector()
LowNodesDFS <- vector()
HighNodesDFS <- vector()
LowTreatmentDFS <- vector()
HighTreatmentDFS <- vector()
LowVolumeDFS <- vector()
HighVolumeDFS <- vector()
LowHistologyDFS <- vector()
HighHistologyDFS <- vector()

#contin
print(pp.data)

for(i in 1:6){
  if(pp.data$pData[[i]]$xvar.name == 'Vol'){
    Volnum <- i
  }
  if(pp.data$pData[[i]]$xvar.name == 'Age'){
    Agenum <- i
  }
  if(pp.data$pData[[i]]$xvar.name == 'Treatment'){
    Treatmentnum <- i
  }
  if(pp.data$pData[[i]]$xvar.name == 'XM.Fp'){
    Fpnum <- i
  }
  
  if(pp.data$pData[[i]]$xvar.name == 'XM.ve'){
    venum <- i
  }
  
  if(pp.data$pData[[i]]$xvar.name == 'Nodes'){
    Nodesnum <- i
  }
  if(pp.data$pData[[i]]$xvar.name == 'Histology'){
    Histologynum <- i
  }
  if(pp.data$pData[[i]]$xvar.name == 'Stage'){
    Stagenum <- i
  }
}
print(Agenum)

lowerindices <- seq(1,(36*2),2)
upperindices <- seq(2, (36*2), 2)
print(lowerindices)
#pull out pp info 
# LowTreatmentDFS[1:nrow(cervixData)] <- pp.data$pData[[Treatmentnum]]$yhat[1:nrow(cervixData)]
# HighTreatmentDFS[1:nrow(cervixData)] <- pp.data$pData[[Treatmentnum]]$yhat[(nrow(cervixData) + 1):(nrow(cervixData)*2)]

LowTreatmentDFS[1:nrow(cervixData)] <- pp.data$pData[[Treatmentnum]]$yhat[lowerindices]
HighTreatmentDFS[1:nrow(cervixData)] <- pp.data$pData[[Treatmentnum]]$yhat[upperindices]

LowNodesDFS[1:nrow(cervixData)] <- pp.data$pData[[Nodesnum]]$yhat[lowerindices]
HighNodesDFS[1:nrow(cervixData)] <- pp.data$pData[[Nodesnum]]$yhat[upperindices]

Lower.Age <- pp.data$pData[[Agenum]]$yhat[lowerindices]
Upper.Age <- pp.data$pData[[Agenum]]$yhat[upperindices]
lab.Age <-  expression(paste('Age'))

Lower.Fp <- pp.data$pData[[Fpnum]]$yhat[lowerindices]
Upper.Fp <- pp.data$pData[[Fpnum]]$yhat[upperindices]
lab.Fp <-  expression(paste(italic(F)[p]))

Lower.ve <- pp.data$pData[[venum]]$yhat[lowerindices]
Upper.ve <- pp.data$pData[[venum]]$yhat[upperindices]
lab.ve <-  expression(paste(italic(v)[e]))
lab.Vol <-  expression(paste('MRI tumour volume'))

LowHistologyDFS[1:nrow(cervixData)] <- pp.data$pData[[Histologynum]]$yhat[lowerindices]
HighHistologyDFS[1:nrow(cervixData)] <- pp.data$pData[[Histologynum]]$yhat[upperindices]

# define functions to plot PP 
plotPPdicot <- function(None,With,xlab, ylab, main, labels){
  bp <- data.frame(None, With)
  bp2 <- boxplot(bp, plot = F)
  print(bp)
  plot(c(0.5,1.5), c(mean(None), mean(With)),ylim = labels, pch = 19, xlim = c(0,2), cex = 1, ylab = ylab, xlab = '',main = main, xaxt = 'n', las = 2)
  segments(0.5,bp2$stats[1,1], 0.5, bp2$stats[5,1], lwd = 1.5)
  segments(1.5,bp2$stats[1,2],1.5, bp2$stats[5,2], lwd = 1.5)
  #   axis(side = 2,at = labels, labels = labels)
  axis(side = 1,at = c(0.5,1.5), labels = xlab, cex = 0.8)
}

limits <- c(.25,.75) 
#make 5-year DFS plots (Figure 3) 
print('Writing Figure 5.') 
tiff('Figure3.tiff', width = 8, height = 6, res = 300, units = 'in')
par(mfrow = c(2,3), cex.lab = 1.2, cex.axis = 1.1, oma = c(1,1,1,1))

plotPPdicot(Lower.Fp, Upper.Fp, c('< median', '\u2265 median'), '', lab.Fp, limits)
plotPPdicot(LowTreatmentDFS, HighTreatmentDFS, c('RT', 'CRT'), '3-year DFS probability', expression(paste('Treatment')), limits)
plotPPdicot(LowHistologyDFS, HighHistologyDFS, c('SCC', 'Other'), '', expression(paste('Histology')), limits)  
plotPPdicot(Lower.ve, Upper.ve,  c('< median', '\u2265 median'), '3-year DFS probability', lab.ve, limits)
plotPPdicot(LowNodesDFS, HighNodesDFS, c('-ve', '+ve'), '', expression(paste('Nodal Status')), limits)
plotPPdicot(Lower.Age, Upper.Age,  c('< median', '\u2265 median'), '', lab.Age, limits)

dev.off()

#compute max risk for alterative and null models 
maxrisk.clin <- max(predict(rsf.clinical)$predicted,pred.clinical)
maxrisk.vimp <- max(predict(rsf.all.vimp)$predicted,pred)

#compute c-index for null and alternative models in training data
c.index.struct.train <- concordance.index(predict(rsf.all.vimp)$predicted,cervixData$DFStime, cervixData$status, method = 'noether')
c.index.struct.clinical.train <- concordance.index(predict(rsf.clinical)$predicted,cervixData$DFStime, cervixData$status, method = 'noether')
c.index.final.train <- c.index.struct.train$c.index
c.index.clinical.train <- c.index.struct.clinical.train$c.index

print(paste('c-index of alternative model (train)  = ',c.index.final.train, sep = ""))
print(paste('c-index of null model (train) = ',c.index.clinical.train, sep = ""))

cindex.comp(c.index.struct.train,c.index.struct.clinical.train)
print(paste('P-value  = ',cindex.comp(c.index.struct.train,c.index.struct.clinical.train)$p.value, sep = ""))

#compute c-index for null and alternative models in test (LOO) data
c.index.struct.test <- concordance.index(pred,cervixData$DFStime, cervixData$status, method = 'noether')
c.index.struct.clinical.test <- concordance.index(pred.clinical,cervixData$DFStime, cervixData$status, method = 'noether')
c.index.final.test <- c.index.struct.test$c.index
c.index.clinical.test <- c.index.struct.clinical.test$c.index

print(paste('c-index of null model (test) = ',c.index.clinical.test, sep = ""))
print(paste('c-index of alternative model (test)  = ',c.index.final.test, sep = ""))

cindex.comp(c.index.struct.test,c.index.struct.clinical.test)
print(paste('P-value  = ',cindex.comp(c.index.struct.test,c.index.struct.clinical.test)$p.value, sep = ""))

#print message to alert the user that the warning message printed by R are ok
print('Messages below are nothing to worry about...')