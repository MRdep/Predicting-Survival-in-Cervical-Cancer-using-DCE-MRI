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
source('DicotamiseOnIndex.R')

# set cwd to source file location 
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

##first make Figures 2 and Supp 2 and write out as pdfs
print('----------------------------------------------------------------------------------') 
print('Writing Figures 1 and supplementary Figure 3.')
print('----------------------------------------------------------------------------------') 

source('KMplotter.R')
source('batchKM_2CXM_only.R')

#-------------------------------------------------- 
# read in tracer kinetic and clinical parameters
#-------------------------------------------------- 
model.params.ktrans <- read.csv('Ktrans_features.csv')
model.params.Fp <- read.csv('Fp_features.csv')
model.params.PS <- read.csv('PS_features.csv')
model.params.vp <- read.csv('vp_features.csv')
model.params.ve <- read.csv('ve_features.csv')

clinical.params <- read.csv('ClinicalParameters_UpdatedFeb2017.csv')

clinical.params <- clinical.params
model.params <- cbind(model.params.ktrans['Ktrans.median'], model.params.Fp['Fp.median'],model.params.PS['FE.median'],model.params.vp['vp.median'],model.params.ve['ve.median'])
colnames(model.params) <- c('Median.Ktrans.2CXM', 'Median.Fp.2CXM', 'Median.PS.2CXM', 'Median.vp.2CXM', 'Median.ve.2CXM')
model.params <- data.frame(model.params)
len <- length(model.params$Median.Fp.2CXM)

#-------------------------------------------------- 
# Do some data preprocessing and cleaning
#-------------------------------------------------- 
classFp2CXM <- DicotamiseOnIndex(model.params$Median.Fp.2CXM, 19)
classFE2CXM <- DicotamiseOnIndex(model.params$Median.PS.2CXM, 24)
classvp2CXM <- DicotamiseOnIndex(model.params$Median.vp.2CXM, 22)
classve2CXM <- DicotamiseOnIndex(model.params$Median.ve.2CXM, 18)
classKtrans2CXM <- DicotamiseOnIndex(model.params$Median.Ktrans.2CXM, 26)

classAge <- DicotamiseOnIndex(clinical.params$Age, 18)
classVol <- DicotamiseOnIndex(clinical.params$Vol, 26)

classStage <- clinical.params$Stage_2level
classNodes <- clinical.params$Nodes
classHistology<- vector()
classTreatment <- vector()

# change order of factors
for (i in 1:len){
  if(clinical.params$Histology[i] == 'SCC'){
    classHistology[i] = 0
  } 
  else{
    classHistology[i] = 1
  }
}

# change order of factors
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

# package 0
#-------------------------------------------------------------------
#Compute univariate p-values
#-------------------------------------------------------------------
# compute Cox model univariate HR
print(coxph(Surv(DFStime, status) ~ XM.Ktrans, data = cervixDicot))
print(coxph(Surv(DFStime, status) ~ XM.Fp, data = cervixDicot))
print(coxph(Surv(DFStime, status) ~ XM.FE, data = cervixDicot))
print(coxph(Surv(DFStime, status) ~ XM.vp, data = cervixDicot))
print(coxph(Surv(DFStime, status) ~ XM.ve, data = cervixDicot))
print(coxph(Surv(DFStime, status) ~ Vol, data = cervixDicot))
print(coxph(Surv(DFStime, status) ~ Age, data = cervixDicot))
print(coxph(Surv(DFStime, status) ~ Stage, data = cervixDicot))
print(coxph(Surv(DFStime, status) ~ Nodes, data = cervixDicot))
print(coxph(Surv(DFStime, status) ~ Histology, data = cervixDicot))
print(coxph(Surv(DFStime, status) ~ Treatment, data = cervixDicot))

# package 1
#-------------------------------------------------------------------
#Fit RSF model with all clinicopathologic and tracer kinetic variables. Use VIMP to rank variables!
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

importances.all.vimp <- rsf.all.vimp$importance

set.seed(123)
rsf.top6 <- rfsrc(Surv(DFStime, status) ~ XM.Fp + Age + Treatment + Histology + XM.Ktrans + Nodes, ntree=ntree, data = cervixData, nodesize = nodesize, bootstrap = c("by.root"),importance = T, nsplit = nsplit)

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

print('Running VIMP bootstrap experiment.') 

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

#----------------------------------- 
#output VIMP statistics into a table.
# ----------------------------------- 
vimp <- array(NA, dim = c(length(medians), 2)) 

vimp[, 1] <- names(rsf.all.vimp$importance)
for(i in 1:length(medians)){
vimp[i, 2] <- paste(format(signif(medians[i], 2)), ' (',format(signif(lower[i], 2)),', ',format(signif(upper[i], 2)),')', sep = "")
}
index.ordered <- c(8, 11, 9, 10, 6, 7, 2, 3, 1, 4, 5)
vimp <- vimp[index.ordered,]

clinical.vimp <- rbind(lower.clinical, medians.clinical, upper.clinical)
colnames(clinical.vimp) <- names(rsf.clinical$importance)

# write out csv files of VIMP. 
write.csv(vimp, file = 'vimp.csv')
write.csv(clinical.vimp, file = 'clinical_vimp.csv')

# package 2 
#------------------------------------------------------------
# PLOT KM-curves for predictions from NULL and ALTERNATIVE models
#------------------------------------------------------------

pred.test <- vector()
pred.test.clinical <- vector()
yearProb <- vector()
yearProb.clinical <- vector()

#compute the prediction for the training data 
data.train <- cervixData
rsf.all.vimp <- rfsrc(Surv(DFStime,status) ~ XM.Fp + XM.Ktrans + Nodes + Histology + Treatment + Age, ntree=ntree,data = data.train, nsplit = nsplit)
rsf.clinical <- rfsrc(Surv(DFStime, status) ~., ntree = ntree, data = data.train.clinical, nsplit = nsplit)
pred.train <- predict(rsf.all.vimp, newdata = cervixData)$predicted
pred.train.clinical <- predict(rsf.clinical, newdata = cervixData.clinical)$predicted

#run a LOO crossvalidation, computing the predicted output for the left out individual 
for(i in 1:nrow(cervixData)){
  if(i == 18){print('Half way...')}
  data.train <- cervixData[-i, ]
  data.test <- cervixData[i,]
  data.train.clinical = cervixData.clinical[-i,]
  data.test.clinical = cervixData.clinical[i,]
  rsf.all.vimp <- rfsrc(Surv(DFStime,status) ~ XM.Fp + XM.Ktrans + Nodes + Histology + Treatment + Age, ntree=ntree,data = data.train, nsplit = nsplit)
  rsf.clinical <- rfsrc(Surv(DFStime, status) ~., ntree = ntree, data = data.train.clinical, nsplit = nsplit)
  # pull out the prediction for the test data   
  pred.test[i] <- predict(rsf.all.vimp, newdata = cervixData)$predicted[i]
  pred.test.clinical[i] <- predict(rsf.clinical, newdata = cervixData.clinical)$predicted[i]
}

MortIndexall.train <- vector()
MortIndexclinical.train <- vector()
MortIndexall.test <- vector()
MortIndexclinical.test <- vector()
#find indexes of patients who were predicted to have below median mortality based on :
index_lowPredMort.test <- which(pred.test < median(pred.test))
index_lowPredMort.train <- which(pred.train < median(pred.train))
#find indexes of patients who were predicted to have above median mortality:
index_highPredMort.test <- which(pred.test >= median(pred.test))
index_highPredMort.train <- which(pred.train >= median(pred.train))

MortIndexall.test[index_lowPredMort.test] <- 'L'
MortIndexall.test[index_highPredMort.test] <- 'H'
MortIndexall.train[index_lowPredMort.train] <- 'L'
MortIndexall.train[index_highPredMort.train] <- 'H'
#repeat for clinical only model 
indexclinical_lowPredMort.test <- which(pred.test.clinical < median(pred.test.clinical))
indexclinical_highPredMort.test <- which(pred.test.clinical >= median(pred.test.clinical))
MortIndexclinical.test[indexclinical_lowPredMort.test] <- 'L'
MortIndexclinical.test[indexclinical_highPredMort.test] <- 'H'
indexclinical_lowPredMort.train <- which(pred.train.clinical< median(pred.train.clinical))
indexclinical_highPredMort.train <- which(pred.train.clinical >= median(pred.train.clinical))
MortIndexclinical.train[indexclinical_lowPredMort.train] <- 'L'
MortIndexclinical.train[indexclinical_highPredMort.train] <- 'H'

#output Figure 3
png( "Figure3.png", width = (15), height = (61*2/5), units = 'in', res = 600)
par(mfrow=c(2,1), cex.lab = 1.5,cex.axis = 1.5)
#perform cox regression as multivariate mortality as a predictor
cox1 <- coxph(Surv(DFStime, status) ~ MortIndexall.train, data = cervixDicot)
cox2 <- coxph(Surv(DFStime, status) ~ MortIndexclinical.train, data = cervixDicot)
cox3 <- coxph(Surv(DFStime, status) ~ MortIndexall.test, data = cervixDicot)
cox4 <- coxph(Surv(DFStime, status) ~ MortIndexclinical.test, data = cervixDicot)

KM.modelall.train <- prodlim(Surv(DFStime, status)~ MortIndexall.train, data = cervixDicot)
KM.modelclinical.train <- prodlim(Surv(DFStime, status)~ MortIndexclinical.train, data = cervixDicot)
KM.modelall.test <- prodlim(Surv(DFStime, status)~ MortIndexall.test, data = cervixDicot)
KM.modelclinical.test <- prodlim(Surv(DFStime, status)~ MortIndexclinical.test, data = cervixDicot)

pvalheight = -50
mai = c(2.0,2.0,1.5,1.5)
oma = c(0, 0, 0, 0)

q <- par(mai = mai, oma = oma, xaxs = 'i',yaxs = 'i', ps = 10, las = 1, cex.axis = 4.5, cex.lab = 4.5, cex.main = 3, bty = 'o')

q <- plot(KM.modelclinical.test, 
     data = cervixDicot,
     border = T,
                                       axes = 'F',
                                       atrisk = TRUE,
                                       xlim = c(0,7),
                                       ylim = c(0,1),
                                       lwd = 5,
                                       lty = c(1,2,3,4),
                                       col = c(1,1,1,1),
                                       legend = F,
                                       atRisk.title ='',
                                       atRisk.adj = c(1.2),
                                       atRisk.cex=2.2,
                                       marktime.cex = 2.6,
                                       atRisk.labels=  c('<  group median', '\u2265 group median'),
                                       marktime = TRUE,
                                       atrisk.line = c(5,8),
                                       xlab = NA, 
                                       ylab = NA,
                                       percent = T,
                                       logrank = FALSE,
                                       confint = FALSE)

labelsforLeg1 = '< median predicted mortality'
labelsforLeg2 = '\u2265 median predicted mortality' 


titleForPlot = expression('Null model')

legend(0.0, .175, labelsforLeg1, cex = 3.0, fill = F, border = F, lty = 2, lwd = 3, bty = "n", xpd = T)
legend(0.0, .125, labelsforLeg2, cex = 3.0, fill = F, border = F, lty = 1, lwd = 3, bty = "n", xpd = T)

title(titleForPlot,cex.main = 5.0)
axis(side = 1,lwd = 3, line = -4.3, labels = FALSE, tick = FALSE)
axis(side = 1,lwd = 3, line = -4.3, labels = FALSE)
axis(side = 1,lwd = 3, line = -3.3, tick = FALSE)
axis(side = 2, lwd = 3, line = -3.0, tck = 0, labels = c(0,20,40,60,80,100), at = c(0,0.2,0.4,0.6,0.8,1))
mtext(side = 1, "Years", line = 2.0, cex = 4.5)
mtext(side = 2, "Percentage Recurrence Free", line = 6.4, cex = 4.5, las = 0)

coeffs <- coef(summary(cox4))
pval = as.matrix(coeffs[,5])
mtext(side = 1, paste('P-value',' =',format(signif(pval,2))), at = c(5), line = pvalheight + 15, cex = 4.0)

p <- par(mai = mai,xaxs = 'i',yaxs = 'i', ps = 10, las = 1, cex.axis = 4.5, cex.lab = 4.5, cex.main = 3, bty = 'o')

p <- plot(KM.modelall.test,
     data = cervixDicot,
     border = T,
                                       axes = 'F',
                                       atrisk = TRUE,
                                       xlim = c(0,7),
                                       ylim = c(0,1),
                                       lwd = 5,
                                       lty = c(1,2,3,4),
                                       col = c(1,1,1,1),
                                       legend = F,
                                       atRisk.title ='',
                                       atRisk.adj = c(1.2),
                                       atRisk.cex=2.2,
                                       marktime.cex = 2.6,
                                       atRisk.labels=  c('<  group median', '\u2265 group median'),
                                       marktime = TRUE,
                                       atrisk.line = c(5,8),
                                       xlab = NA, 
                                       ylab = NA,
                                       percent = T,
                                       logrank = FALSE,
                                       confint = FALSE)

labelsforLeg1 = '< median predicted mortality'
labelsforLeg2 = '\u2265 median predicted mortality' 
titleForPlot = expression('Alternative model')

legend(0.0, .175, labelsforLeg1, cex = 3.0, fill = F, border = F, lty = 2, lwd = 3, bty = "n", xpd = T)
legend(0.0, .125, labelsforLeg2, cex = 3.0, fill = F, border = F, lty = 1, lwd = 3, bty = "n", xpd = T)

title(titleForPlot,cex.main = 5.0)

axis(side = 1,lwd = 3, line = -4.3, labels = FALSE, tick = FALSE)
axis(side = 1,lwd = 3, line = -4.3, labels = FALSE)
axis(side = 1,lwd = 3, line = -3.3, tick = FALSE)
axis(side = 2, lwd = 3, line = -3, tck = 0, labels = c(0,20,40,60,80,100), at = c(0,0.2,0.4,0.6,0.8,1))

mtext(side = 1, "Years", line = 2.0, cex = 4.5)
mtext(side = 2, "Percentage Recurrence Free", line = 6.4, cex = 4.5, las = 0)

coeffs <- coef(summary(cox3))
pval = as.matrix(coeffs[,5])
mtext(side = 1, paste('P-value',' =',format(signif(pval,2))), at = c(5), line = pvalheight + 10, cex = 4.0)

dev.off()


# package 3 
#------------------------------------------------------------
# PLOT partial plots with p-values for parameters in the alternative model
#------------------------------------------------------------
# make DFS risk prediction partial plots (Figure 4) 
rsf.clinical <- rfsrc(Surv(DFStime,status) ~ Stage +  Histology + Treatment +  Vol + Age + Nodes, ntree=ntree, data = cervixData,nodesize = nodesize, bootstrap = c("by.root"),importance = T, nsplit = nsplit)
rsf.all.vimp <- rfsrc(Surv(DFStime,status) ~ Histology + Treatment + Age + Nodes + XM.Fp + XM.Ktrans, ntree=ntree,data = cervixData,nodesize = nodesize, bootstrap = c("by.root"),importance = T, nsplit = nsplit)

plotPPdicot <- function(obj, time , variable, xlab, ylab, main, labels){
  
partial.obj <- partial.rfsrc(obj, partial.type = 'mort', partial.xvar = variable, partial.values = sort(unique(rsf.all.vimp$xvar[, variable])),oob = T)

point <- apply(partial.obj$survOutput,2,mean,na.rm = TRUE)
n_comp <- 1
test <- t.test(partial.obj$survOutput[,1],partial.obj$survOutput[,2], alternative = c("two.sided"))
pval <- p.adjust(test$p.value, method = 'bonferroni', n = n_comp)

diff <- sweep(partial.obj$survOutput,2,point,  FUN="-")/sqrt(36)
bp <- sweep(diff, 2, point,  FUN="+")
ndeath = 18
pvallab = paste('P-value = ', format(signif(pval,2)))
bp2 <- boxplot(bp, plot = F)
plot(c(0.5,1.5), c(point[1]/ndeath, point[2]/ndeath),ylim = labels, pch = 19, xlim = c(0,2), cex = 1, ylab = ylab, xlab = '',main = main, xaxt = 'n', las = 2)
segments(0.5, bp2$stats[1,1]/ndeath, 0.5, bp2$stats[5,1]/ndeath, lwd = 1.5)
segments(1.5, bp2$stats[1,2]/ndeath,1.5, bp2$stats[5,2]/ndeath, lwd = 1.5)
axis(side = 1,at = c(0.5,1.5), labels = xlab, cex = 0.9, padj = 0.5)
text(1.3, 0.60, pvallab, cex = 1.1)

}

#make 3-year DFS plots (Supp Figure 2) 
print('Writing Figure 4.') 
png('Figure4.png', width = 8, height = 6, res = 300, units = 'in')
par(mfrow = c(2,3), cex.lab = 1.1, cex.axis = 1.1, oma = c(1,1,1,1))
# limits <- c(.25,.75) 
limits <- c(0.15,0.65) 

plotPPdicot(rsf.all.vimp, rsf.all.vimp$time.interest[16], 'XM.Fp', c('< 50th\npercentile', '\u2265 50th\npercentile'), 'Risk of recurrence/death', expression(paste(italic(F)[p])), limits)
plotPPdicot(rsf.all.vimp, rsf.all.vimp$time.interest[16], 'Treatment', c('RT', 'CRT'), 'Risk of recurrence/death', expression(paste('Treatment')), limits)
plotPPdicot(rsf.all.vimp, rsf.all.vimp$time.interest[16], 'Histology', c('SCC', 'Other'), 'Risk of recurrence/death', expression(paste('Histology')), limits)
plotPPdicot(rsf.all.vimp, rsf.all.vimp$time.interest[16], 'Nodes', c('-ve', '+ve'), 'Risk of recurrence/death', expression(paste('Nodal Status')), limits)
plotPPdicot(rsf.all.vimp, rsf.all.vimp$time.interest[16], 'Age', c('< 47th\npercentile', '\u2265 47th\npercentile'), 'Risk of recurrence/death', expression(paste('Age')), limits)
plotPPdicot(rsf.all.vimp, rsf.all.vimp$time.interest[16], 'XM.Ktrans', c('< 69th\npercentile', '\u2265 69th\npercentile'), 'Risk of recurrence/death', expression(paste(italic(K)^{trans})), limits)

dev.off()

# package 4
#------------------------------------------------------------
# COMPUTE c-indices for null and alternative model...
#------------------------------------------------------------

#compute max risk for alterative and null models 
maxrisk.clin <- max(predict(rsf.clinical)$predicted,pred.train.clinical)
maxrisk.vimp <- max(predict(rsf.all.vimp)$predicted,pred.train)

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
c.index.struct.test <- concordance.index(pred.test,cervixData$DFStime, cervixData$status, method = 'noether')
c.index.struct.clinical.test <- concordance.index(pred.test.clinical,cervixData$DFStime, cervixData$status, method = 'noether')
c.index.final.test <- c.index.struct.test$c.index
c.index.clinical.test <- c.index.struct.clinical.test$c.index

print(paste('c-index of null model (test) = ',c.index.clinical.test, sep = ""))
print(paste('c-index of alternative model (test)  = ',c.index.final.test, sep = ""))

cindex.comp(c.index.struct.test,c.index.struct.clinical.test)
print(paste('P-value  = ',cindex.comp(c.index.struct.test,c.index.struct.clinical.test)$p.value, sep = ""))

#print message to alert the user that the warning message printed by R are ok
print('Messages below are nothing to worry about...')