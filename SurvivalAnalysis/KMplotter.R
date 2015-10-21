# function to make Kaplan-Meier plots. Run from batchKM.R.  
KMplotter <- function(param){
  
  # read in tracer kinetic and clinical parameters
  model.params <- read.csv('TracerKineticParameters.csv')
  clinical.params <-read.csv('ClinicalParameters.csv')
  
  # Define vectors to hold classifications
  classFp2CXM <- vector()
  classFE2CXM <- vector()
  classve2CXM <- vector()
  classvp2CXM <- vector()
  
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
  median.Fp.2CXM <- median(model.params$Median.Fp.2CXM)
  median.PS.2CXM <- median(model.params$Median.PS.2CXM)
  median.Ktrans.Tofts <- median(model.params$Median.Ktrans.Tofts)
  median.Ktrans.Etofts <- median(model.params$Median.Ktrans.Etofts)

  median.Age <- median(clinical.params$Age)
  median.Vol <- median(clinical.params$Vol) 
  
  len <- length(model.params$Median.Fp)
  
  DFStime <- clinical.params$DFS
  status <- clinical.params$DFSstatus

# Compute classifications
# 2cxm parameters 
for (i in 1:len){
  if(model.params$Median.Fp.2CXM[i] < median.Fp.2CXM){
    classFp2CXM[i] = 'Fp < Median'
  } 
  else{
    classFp2CXM[i] = 'Fp > Median'
  }
}

for (i in 1:len){
  if(model.params$Median.PS.2CXM[i] < median.PS.2CXM){
    classFE2CXM[i] = 'FE < Median'
  } 
  else{
    classFE2CXM[i] = 'FE > Median'
  }
}

# Tofts parameters 
for (i in 1:len){
  if(model.params$Median.Ktrans.Tofts[i] < median.Ktrans.Tofts){
    classKtransTOFTS[i] = 'Ktrans < Median'
  } 
  else{
    classKtransTOFTS[i] = 'Ktrans > Median'
  }
}

# ETOFTS params
for (i in 1:len){
  if(model.params$Median.Ktrans.Etofts[i] < median.Ktrans.Etofts){
    classKtransETOFTS[i] = 'Ktrans < Median'
  } 
  else{
    classKtransETOFTS[i] = 'Ktrans > Median'
  }
}

for (i in 1:len){
  if(clinical.params$Age[i] < median.Age){
    classAge[i] = 'Age < Median'
  } 
  else{
    classAge[i] = 'Age > Median'
  }
}


for (i in 1:len){
  if(clinical.params$Vol[i] < median.Vol){
    classVol[i] = 'Volume < Median'
  } 
  else{
    classVol[i] = 'Volume > Median'
  }
}

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
  if(clinical.params$Treatment[i] == 1){
    classTreatment[i] = 1
  } 
  else{
    classTreatment[i] = 0
  }
}

# Define plotting parameters
bottomline <- 7
topline <- 10
mai = c(2.0,3.0,1.5,1.5)

if (param =='Fp')
{
  labelsforLeg = c(expression(paste(italic(F)[p],' ','>',' ','group median')), expression(paste(italic(F)[p],' ','<',' ', 'group median')))
  predictor = as.numeric(factor(classFp2CXM))
  titleForPlot = expression(paste('2CXM ',italic(F)[p]))
  location = 'bottomleft'
  ones <-which(predictor == 2)
  predictor[ones] = 0
  atrisklabs = c('> group median','< group median')
  pvalheight = -37.
  atrisk.line = c(bottomline,topline)
# To univariate cox to compute HR and p value
  cox1 <- coxph(Surv(DFStime,status) ~ classFp2CXM, method = 'breslow')
}

if (param =='FE')
{
  labelsforLeg = c(expression(paste(italic(PS),' ','>',' ','group median')), expression(paste(italic(PS),' ','<',' ', 'group median')))
  predictor = as.numeric(factor(classFE2CXM))
  titleForPlot = expression(paste('2CXM ',italic(PS)))
  location = 'bottomleft'
  ones <-which(predictor == 2)
  predictor[ones] = 0
  atrisklabs = c('> group median','< group median')
  pvalheight = -37.
  atrisk.line = c(bottomline,topline)

  cox1 <- coxph(Surv(DFStime,status) ~ classFE2CXM, method = 'breslow')
}

if (param =='KtransT')
{
  labelsforLeg = c(expression(paste(italic(K)^{trans},' ','>',' ','group median')), expression(paste(italic(K)^{trans},' ','<',' ', 'group median')))
  predictor = as.numeric(factor(classKtransTOFTS))
  titleForPlot = expression(paste('Tofts  ',italic(K)^{trans}))
  location = 'bottomleft'
  ones <-which(predictor == 2)
  predictor[ones] = 0
  atrisklabs = c('> group median','< group median')
  pvalheight = -37.
  atrisk.line = c(bottomline,topline)

  cox1 <- coxph(Surv(DFStime,status) ~ classKtransTOFTS, method = 'breslow')
}

if (param =='KtransET')
{
  labelsforLeg = c(expression(paste(italic(K)^{trans},' ','>',' ','group median')), expression(paste(italic(K)^{trans},' ','<',' ', 'group median')))
  predictor = as.numeric(factor(classKtransETOFTS))
  titleForPlot = expression(paste('Extended Tofts  ',italic(K)^{trans}))
  location = 'bottomleft'
  ones <-which(predictor == 2)
  predictor[ones] = 0
  atrisklabs = c('> group median','< group median')
  pvalheight = -37.
  atrisk.line = c(bottomline,topline)
  cox1 <- coxph(Surv(DFStime,status) ~ classKtransETOFTS, method = 'breslow')
}


if (param =='vol')
{
  labelsforLeg = c('volume < group median', 'volume > group median')
  predictor = as.numeric(factor(classVol))
  titleForPlot = expression('MRI volume')
  location = 'bottomleft'
  atrisklabs = c('< group median','> group median')
  pvalheight = -45.
  atrisk.line = c(bottomline,topline)

  cox1 <- coxph(Surv(DFStime,status) ~ classVol, method = 'breslow')
}

if (param =='age')
{
  labelsforLeg = c('age < group median', 'age > group median')
  predictor = as.numeric(factor(classAge))
  titleForPlot = expression('Age')
  location = 'bottomleft'
  atrisklabs = c('< group median','> group median')
pvalheight = -45.
  atrisk.line = c(bottomline,topline)

  cox1 <- coxph(Surv(DFStime,status) ~ classAge, method = 'breslow')
}


if (param =='stage')
{
  labelsforLeg = c('Stage 1 or 2', 'Stage 3 or 4')
  predictor = as.numeric(factor(classStage))
  titleForPlot = expression('FIGO Stage')
  location = 'bottomright'
  atrisklabs = c('FIGO 1-2','FIGO 3-4')
  pvalheight = -45.
  atrisk.line = c(bottomline,topline)

  cox1 <- coxph(Surv(DFStime,status) ~ classStage, method = 'breslow')
}

if (param =='histology')
{
  labelsforLeg = c('SCC', 'Other')
  predictor = as.numeric(factor(classHistology))
  titleForPlot = expression('Histological subtype')
  location = 'bottomleft'
  ones <-which(predictor == 2)
  predictor[ones] = 0
  atrisklabs = c('    SCC','  Other')
  pvalheight = -45.
  atrisk.line = c(bottomline,topline)

    cox1 <- coxph(Surv(DFStime,status) ~ classHistology, method = 'breslow')
}

if (param =='nodes')
{
  labelsforLeg = c('Negative', 'Positive')
  predictor = as.numeric(factor(classNodes))
  titleForPlot = expression('Nodal Status')
  location = 'bottomleft'
  atrisklabs = c('Negative','Positive')
  pvalheight = -45.
  atrisk.line = c(bottomline,topline)
  cox1 <- coxph(Surv(DFStime,status) ~ classNodes, method = 'breslow')
}

if (param =='treatment')
{
  labelsforLeg = c('ChemoRT', 'RT')
  predictor = as.numeric(factor(classTreatment))
  titleForPlot = expression('Treatment')
  location = 'bottomleft'
  atrisklabs = c('ChemoRT','        RT')
  pvalheight = -45.
  atrisk.line = c(bottomline,topline)
   cox1 <- coxph(Surv(DFStime,status) ~ classTreatment, method = 'breslow')
}

# Make plot
ddist <- datadist(predictor)
options(datadist = 'ddist')

survobj = Surv(DFStime,status)

fits <- prodlim(Surv(DFStime,status) ~ predictor,
        exact = TRUE)
diff <-survdiff(Surv(DFStime,status) ~ predictor)

p <- par(mai = mai, xaxs = 'i',yaxs = 'i', ps = 10, las = 1, cex.axis = 4, cex.lab = 4.0, bty = 'o')

p <- plot(fits,axes = 'F',
     atrisk = TRUE,
     xlim = c(0,7),
     ylim = c(0,1),
     legend.x=c(location),
     lwd = 5,
     lty = c(1,2,3,4),
     col = c(1,1,1,1),
     legend.title="",
     legend.y.intersp=1,
     background.horizontal=seq(0,1,.20),
     legend.cex=4,
     legend.legend = labelsforLeg,
     atRisk.title ='',
     atRisk.adj = c(1.4),
     atRisk.cex=2,
     marktime.cex = 2.3,
    atRisk.labels= atrisklabs,
      marktime = TRUE,
     atrisk.line = atrisk.line,
     xlab = NA, 
     ylab = NA,
     percent = T,
    logrank = FALSE,
    confint = FALSE)
title(titleForPlot,cex.main = 4.0)


 axis(side = 1,lwd = 3, line = -4.3, labels = FALSE, tick = FALSE)
 axis(side = 1,lwd = 3, line = -4.3, labels = FALSE)
 axis(side = 1,lwd = 3, line = -3.3, tick = FALSE)
 axis(side = 2, lwd = 3, line = -2, tck = 0, labels = c(0,20,40,60,80,100), at = c(0,0.2,0.4,0.6,0.8,1))

 mtext(side = 1, "Years", line = 2.0, cex = 3.0)
 mtext(side = 2, "Percentage Recurrence Free", line = 6.4, cex = 3.0, las = 0)

 coeffs <- coef(summary(cox1))
 pval = as.matrix(coeffs[,5])
 mtext(side = 1, paste('P-value',' =',format(signif(pval,2),nsmall = 2)), line = pvalheight, cex = 3.0)

# Return plot object
return(p)
}
