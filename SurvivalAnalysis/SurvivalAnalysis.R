library(randomForestSRC)
library(survival)
library(rms)
library(prodlim)
library(boot)
library(survcomp)
library(datasets)
library(MASS)

# set cwd to source file location 
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

#compile source files needed to run the script 
source('OutcomePartialPlots.R')
source('DicotamiseOnMedian.R')
source('KMplotter.R')

#first make Figures 3 and 4 and write out as pdfs
print('----------------------------------------------------------------------------------') 
print('Writing Figures 3 and 4 to disk.')
print('----------------------------------------------------------------------------------') 

# Create a folder to hold figures
dir.create('Figures')

source('batchKM.R')

#Next run the Random Survival Forest Analysis and print out figures and VIMP tables. 
#read in data. 
model.params <- read.csv('TracerKineticParameters.csv')
clinical.params <-read.csv('ClinicalParameters.csv')

len <- length(model.params$Median.Fp)

# dicotmise variables.
classFp2CXM <- DicotamiseOnMedian(model.params$Median.Fp.2CXM)
classFE2CXM <- DicotamiseOnMedian(model.params$Median.PS.2CXM)
classve2CXM <- DicotamiseOnMedian(model.params$Median.ve.2CXM)
classvp2CXM <- DicotamiseOnMedian(model.params$Median.vp.2CXM)

classKtransTOFTS <- DicotamiseOnMedian(model.params$Median.Ktrans.Tofts)
classveTOFTS <- DicotamiseOnMedian(model.params$Median.ve.Tofts)

classKtransETOFTS<- DicotamiseOnMedian(model.params$Median.Ktrans.Etofts)
classvpETOFTS<- DicotamiseOnMedian(model.params$Median.vp.Etofts)
classveETOFTS<- DicotamiseOnMedian(model.params$Median.ve.Etofts)

classAge <- DicotamiseOnMedian(clinical.params$Age)
classVol <- DicotamiseOnMedian(clinical.params$Vol)
classStage <- clinical.params$Stage_2level
classNodes <- clinical.params$Nodes

classHistology<- vector()
classTreatment <- vector()

for (i in 1:len){
  if(clinical.params$Histology[i] == 'SCC'){
    classHistology[i] = 1
  } 
  else{
    classHistology[i] = 0
  }
}

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
cervixDicot <- data.frame(DFStime, status,classFp2CXM,classFE2CXM,classve2CXM,classvp2CXM,classKtransTOFTS,classveTOFTS,classKtransETOFTS,classvpETOFTS,classveETOFTS,classAge,classVol,classStage,classNodes,classHistology, classTreatment)
colnames(cervixDicot) <- c('DFStime', 'status','XM.Fp','XM.FE','XM.ve', 'XM.vp','Tofts.Ktrans', 'Tofts.ve','ETofts.Ktrans','ETofts.vp','ETofts.ve', 'Age', 'Vol', 'Stage','Nodes','Histology','Treatment')

#make a data frame containing only clinicopathologic variables. 
cervixDicot.clinical <- cbind(cervixDicot[,1:2], cervixDicot[,12:17])

#-------------------------------------------------------------------
#Fit RSF model with all clinicopathologic and tracer kinetic variables 1000 DFStimes in a bootstrap experiment
#-------------------------------------------------------------------
#define a variable to hold the number of observations 
n_pats <- nrow(cervixDicot)
#seed RNG. 
set.seed(123)
# set number of trees.
ntree = 1000
#fit models once first.  
rsf.all.vimp <- rfsrc(Surv(DFStime, status) ~., ntree=ntree,data = cervixDicot)
rsf.clinical <- rfsrc(Surv(DFStime, status)~ ., ntree = ntree, data = cervixDicot.clinical)
importances.all.vimp <- rsf.all.vimp$importance

n_repeats <- 1000

importances.vimp <- matrix(nrow = n_repeats, ncol = length(importances.all.vimp))
importances.clinical <- matrix(nrow = n_repeats, ncol = ncol(cervixDicot[,12:17]))
colnames(importances.vimp) <- names(importances.all.vimp)
colnames(importances.clinical) <- names(rsf.clinical$importance)

print('Running bootstrap experiment.') 

for(i in 1:n_repeats){
  remainder <- i%%100
  if(remainder == 0){print(paste('MC repeat = ', i, '/1000', sep = ""))}
  # take a bootstrap sample    
  train.sample <- sample(1:n_pats, n_pats, replace = T)  
  data.train <- cervixDicot[train.sample, ]
  data.train.clinical <- cervixDicot.clinical[train.sample, ]
  # fit models   
  rsf.all.vimp <- rfsrc(Surv(DFStime, status) ~., ntree=ntree,data = data.train)  
  rsf.clinical <- rfsrc(Surv(DFStime, status)~ ., ntree = ntree, data = data.train.clinical)
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

nboot = 1000

bootstrapped.medians <- matrix(nrow = nboot, ncol = length(rsf.all.vimp$importance))
bootstrapped.medians.clinical <- matrix(nrow = nboot, ncol = length(rsf.clinical$importance))

# set the signficance criterion for the CI. Bonferroni corrected for the total number of inferences made. Take 2 to subtract out DFStime and status from cervixDicot. 
criterion = 0.05/(ncol(cervixDicot) - 2)

lower <- vector()
upper <- vector()

lower.clinical <- vector()
upper.clinical <- vector()

for(i in 1:length(importances.all.vimp)){
  boot.out <- boot(data=importances.vimp[,i], statistic=bootThetaQuantile,R=nboot)
  bootstrapped.medians[,i] <- boot.out$t
  ci <- boot.ci(boot.out,conf = 1- criterion, type = 'bca')
  lower[i] <- ci$bca[4] 
  upper[i] <- ci$bca[5] 
}

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
x.labels <- c(expression('2CXM F'[p]),
              expression('2CXM F'[E]),
              expression('2CXM v'[p]),
              expression('2CXM v'[e]),
              expression('Tofts K'^{trans}),
              expression('Tofts v'[e]),
              expression('E. Tofts K'^{trans}),
              expression('E. Tofts v'[p]),
              expression('E. Tofts v'[e]),
              expression('Age'),
              expression('Volume'),
              expression('FIGO Stage'),
              expression('Nodal Status'),
              expression('Histology'),
              expression('Treatment'))

#Make Figure 5. Point estimates and 95% CI on VIMP. 
pdf('Figures/Figure5.pdf', width = 7, height = 5)
par(mai = c(3,1,0.5,0.5), cex.lab = 1.3, cex.axis = 1.3)
col = c('black','red', 'red', 'red', 'black', 'red', 'red','red', 'black','red', 'red', 'black','red', 'black', 'black')
indices <- sort.int(medians, decreasing = T,index.return = T)
col <- c('black','red', 'red', 'red', 'black', 'red', 'red','red', 'black','red', 'red', 'black','red', 'black', 'black')

col <- col[indices$ix] 
plot(medians[indices$ix], xaxt = 'n', yaxt = 'n',xlab = "", ylab = "", pch = 19, ylim = c(-0.001, 0.04), col = col, cex = 0.8)
axis(1, at = 1:length(importances.all.vimp), labels = x.labels[indices$ix], las = 2)
axis(2, at  = c(0,0.02, 0.04), labels = c(0, 0.02, 0.04))
abline(h  =0, lty =2, lwd = 2)
abline(h = 0.0155, lty = 2, lwd = 2)
mtext(side = 2, text = "VIMP", line = 2.5, cex = 1.3)
lower2 <- lower[indices$ix]
upper2 <- upper[indices$ix]
for(i in 1:length(importances.all.vimp)){
  segments(i, lower2[i], x1 = i, y1 = upper2[i], lwd = 2, col = col[i])  
}
dev.off()
print('Writing Figure 5 to disk.') 
print('----------------------------------------------------------------------------------') 


print('Starting Leave-one-out experiment to assess model accuracy.') 
#compute LOOCV c-indices 
rsf.final <- rfsrc(Surv(DFStime,status) ~ Stage +  Histology + Treatment +  XM.Fp + Tofts.Ktrans + ETofts.ve, ntree=ntree,data = cervixDicot)
ntree = 10000
pred <- vector()
pred.clinical <- vector()
yearProb <- vector()
yearProb.clinical <- vector()

for(i in 1:36){
  if(i == 18){print('Half way...')}
  data.train <- cervixDicot[-i, ]
  data.test <- cervixDicot[i,]
  data.train.clinical = cervixDicot.clinical[-i,]
  data.test.clinical = cervixDicot.clinical[i,]
  rsf.all.vimp <- rfsrc(Surv(DFStime,status) ~ Stage +  Histology + Treatment +  XM.Fp + Tofts.Ktrans + ETofts.ve, ntree=ntree,data = data.train)
  rsf.clinical <- rfsrc(Surv(DFStime, status) ~., ntree = ntree, data = data.train.clinical,nsplit = 3)
  pred[i] <- predict(rsf.all.vimp, newdata = cervixDicot)$predicted[i]
  pred.clinical[i] <- predict(rsf.clinical, newdata = cervixDicot.clinical)$predicted[i]
  yearProb[i] <- predict(rsf.all.vimp, newdata = cervixDicot)$survival[i,17] 
  yearProb.clinical[i] <- predict(rsf.clinical, newdata = cervixDicot.clinical)$survival[i,17]
}

c.index.struct <- concordance.index(pred,cervixDicot$DFStime, cervixDicot$status, method = 'noether')
c.index.struct.clinical <- concordance.index(pred.clinical,cervixDicot$DFStime, cervixDicot$status, method = 'noether')
c.index.final <- c.index.struct$c.index
c.index.clinical <- c.index.struct.clinical$c.index
cindex.comp(c.index.struct,c.index.struct.clinical)


# make 5-year DFS risk prediction partial plots (Figure 5) 
rsf.clinical <- rfsrc(Surv(DFStime,status) ~ Stage +  Histology + Treatment +  Vol + Age + Nodes, ntree=ntree,data = cervixDicot)
rsf.all.vimp <- rfsrc(Surv(DFStime,status) ~ Stage +  Histology + Treatment +  XM.Fp + Tofts.Ktrans + ETofts.ve, ntree=ntree,data = cervixDicot)

pp.data.het <- plot.variable(rsf.all.vimp, partial = T,surv.type = "surv",time = DFStime[15], show.plots = F)
pp.data.clinical <- plot.variable(rsf.clinical, partial = T,surv.type = "surv",time = DFStime[15], show.plots = F)

LowFpDFS <- vector()
HighFpDFS <- vector()
LowStageDFS <- vector()
HighStageDFS <- vector()
LowEToftsveDFS <- vector()
HighEToftsveDFS <- vector()
LowTreatmentDFS <- vector()
HighTreatmentDFS <- vector()
LowToftsKtransDFS <- vector()
HighToftsKtransDFS <- vector()
LowHistologyDFS <- vector()
HighHistologyDFS <- vector()

LowFpDFS[1:36] <- pp.data.het$pData[[1]]$yhat[1:36]
HighFpDFS[1:36] <- pp.data.het$pData[[1]]$yhat[37:72]
LowStageDFS[1:36] <- pp.data.het$pData[[2]]$yhat[1:36]
HighStageDFS[1:36] <- pp.data.het$pData[[2]]$yhat[37:72]
LowEToftsveDFS[1:36] <- pp.data.het$pData[[3]]$yhat[1:36]
HighEToftsveDFS[1:36] <- pp.data.het$pData[[3]]$yhat[37:72]
LowTreatmentDFS[1:36] <- pp.data.het$pData[[4]]$yhat[1:36]
HighTreatmentDFS[1:36] <- pp.data.het$pData[[4]]$yhat[37:72]
LowToftsKtransDFS[1:36] <- pp.data.het$pData[[5]]$yhat[1:36]
HighToftsKtransDFS[1:36] <- pp.data.het$pData[[5]]$yhat[37:72]
LowHistologyDFS[1:36] <- pp.data.het$pData[[6]]$yhat[1:36]
HighHistologyDFS[1:36] <- pp.data.het$pData[[6]]$yhat[37:72]

#make 5-year DFS plots (Figure 6) 
print('Writing Figure 6.') 
pdf('Figures/Figure6.pdf')
par(mfrow = c(2,3), cex.lab = 1.2, cex.axis = 1.2)
plot2 <- OutcomePartialPlots(LowStageDFS,HighStageDFS, "surv", c('Early','Late'),'FIGO Stage')
plot1 <- OutcomePartialPlots(LowFpDFS,HighFpDFS, "surv", c('Low','High'), expression(paste('2CXM ', italic(F)[p])))
plot4 <- OutcomePartialPlots(LowTreatmentDFS,HighTreatmentDFS, "surv", c('CT','CT + Rx'),'Treatment')
plot6 <- OutcomePartialPlots(LowHistologyDFS,HighHistologyDFS, "surv", c('SCC','Other'),'Histology')
plot3 <- OutcomePartialPlots(LowEToftsveDFS,HighEToftsveDFS, "surv", c('Low','High'),expression(paste('Extended Tofts ', italic(v)[e])))
plot5 <- OutcomePartialPlots(LowToftsKtransDFS,HighToftsKtransDFS, "surv", c('Low','High'),expression(paste('Tofts ', italic(K)^trans)))

dev.off()

# Plot risk in train vs risk in test (Figure 7)

#compute max risk for alterative and null models 
maxrisk.clin <- max(predict(rsf.clinical)$predicted,pred.clinical)
maxrisk.vimp <- max(predict(rsf.all.vimp)$predicted,pred)

#compute c-index for null and alternative models in training data
c.index.struct.train <- concordance.index(predict(rsf.all.vimp)$predicted,cervixDicot$DFStime, cervixDicot$status, method = 'noether')
c.index.struct.clinical.train <- concordance.index(predict(rsf.clinical)$predicted,cervixDicot$DFStime, cervixDicot$status, method = 'noether')
c.index.final.train <- c.index.struct.train$c.index
c.index.clinical.train <- c.index.struct.clinical.train$c.index
cindex.comp(c.index.struct.train,c.index.struct.clinical.train)

#compute c-index for null and alternative models in test (LOO) data
c.index.struct.test <- concordance.index(pred,cervixDicot$DFStime, cervixDicot$status, method = 'noether')
c.index.struct.clinical.test <- concordance.index(pred.clinical,cervixDicot$DFStime, cervixDicot$status, method = 'noether')
c.index.final.test <- c.index.struct.test$c.index
print(paste('c-index of alternative model  = ',c.index.final.test, sep = ""))
c.index.clinical.test <- c.index.struct.clinical.test$c.index
print(paste('c-index of null model  = ',c.index.clinical.test, sep = ""))
cindex.comp(c.index.struct.test,c.index.struct.clinical.test)
print(paste('P-value  = ',cindex.comp(c.index.struct.test,c.index.struct.clinical.test)$p.value, sep = ""))

y.clin <- 100*pred.clinical/maxrisk.clin
x.clin <- 100*predict(rsf.clinical)$predicted/maxrisk.clin

#Fit linear models
y <- 100*pred/maxrisk.vimp
x <- 100*predict(rsf.all.vimp)$predicted/maxrisk.vimp

clin.lm <- lm(y.clin~x.clin)
vimp.lm <- lm(y~x)
r2.clin <- summary(lm(y.clin~x.clin))$r.squared
r2.vimp <- summary(lm(y~x))$r.squared

print('Writing Figure 7.') 
pdf('Figures/Figure7.pdf', width = 7, height = 8)
#define a single column, two rows 
par(mfrow = c(2,1))
#define marker size based on disease-free time 
markersize <- cervixDicot$DFStime*0.4
pchvals <- cervixDicot$status   
index.p <- which(pchvals == 0)
#set the plot marker to symbol 19
pchvals[index.p] <- 19
#plot test vs train Predicted Risk 
plot(x.clin, y.clin,
     col = 1, cex = markersize, pch = pchvals, xlab = 'Predicted Risk (Train)', 
     ylab = 'Predicted Risk (Test)', main = 'Clinical model',xlim = c(0,100),ylim = c(0,100))

#add legend 
legend("bottomright", legend = c("Recurrences", "Censored"), bty = "n",
 lwd = 2, cex = 1.2, col = c("black", "black"), pch = c(1,19), lty = c(NA,NA))

#add line of identify and best fit line 
abline(a = 0, b = 1, lty = 2)
abline(clin.lm)

# add R^2 values 
R2.clin <- format(r2.clin, digits = 2)
labels <- bquote(italic(R)^2 ~ '=' ~ .(R2.clin))
text(10, 90, label=labels, cex = 1.3)

#add annotations 
segments(58, 91.4, x1 = 64, y1 = 91.4)
segments(25.5, 18.18, x1 = 31.5, y1 = 18.18)
text(56, 91.4, 'A')
text(33.5, 18.18, 'B')

#repeat for clinical and tracer kinetic plot 
plot(x,y,
     col = 1, cex = markersize,pch = pchvals, xlab = 'Predicted Risk (Train)', 
     ylab = 'Predicted Risk (Test)', main = 'Clinical and DCE-MRI model',xlim = c(0,100),ylim = c(0,100))

legend("bottomright", legend = c("Recurrences", "Censored"), bty = "n",
 lwd = 2, cex = 1.2, col = c("black", "black"), pch = c(1,19), lty = c(NA,NA))

abline(a = 0, b = 1, lty = 2)
abline(vimp.lm)
R2.vimp <- format(r2.vimp, digits = 2)

labels <- bquote(italic(R)^2 ~ '=' ~ .(R2.vimp))
text(10, 90, label=labels, cex = 1.3)
# 
segments(53, 69, x1 = 59, y1 = 69)
segments(10, 6, x1 = 15, y1 = 6)
text(61, 69, 'C')
text(17, 6, 'D')

#close the pdf device 
dev.off()
#print message to alert the user that the warning message printed by R are ok
print('Messages below are nothing to worry about...')