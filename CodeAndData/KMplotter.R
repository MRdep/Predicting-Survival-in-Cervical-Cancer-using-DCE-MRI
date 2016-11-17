# function to make Kaplan-Meier plots. Run from batchKM.R.  
library(grDevices)
library(survival)
library(rms)
library(colorspace)
library(ggplot2)
library(prodlim)

KMplotter <- function(param){
  
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
  clinical.params <- read.csv('ClinicalParameters_UpdatedJuly2016.csv')

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
  median.Fp.2CXM <- median(model.params$Median.Fp.2CXM)
  median.PS.2CXM <- median(model.params$Median.PS.2CXM)
  median.Ktrans.2CXM <- median(model.params$Median.Ktrans.2CXM)
  median.vp.2CXM <- median(model.params$Median.vp.2CXM)
  median.ve.2CXM <- median(model.params$Median.ve.2CXM)

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
    classFE2CXM[i] = 'PS < Median'
  } 
  else{
    classFE2CXM[i] = 'PS > Median'
  
  }
}
  set.seed(123)
  model.params$Median.vp.2CXM <- model.params$Median.vp.2CXM + rnorm(length( model.params$Median.vp.2CXM), 0, 0.0001)
  median.vp.2CXM <- median(model.params$Median.vp.2CXM)
  for (i in 1:len){
    if(model.params$Median.vp.2CXM[i] < median.vp.2CXM){
      classvp2CXM[i] = 'vp < Median'
    } 
    else{
      classvp2CXM[i] = 'vp > Median'
    }
  }
  

  for (i in 1:len){
    if(model.params$Median.ve.2CXM[i] < median.ve.2CXM){
      classve2CXM[i] = 've < Median'
    } 
    else{
      classve2CXM[i] = 've > Median'
      
    }
  }
  
  
  for (i in 1:len){
    if(model.params$Median.PS.2CXM[i] < median.PS.2CXM){
      classFE2CXM[i] = 'PS < Median'
    } 
    else{
      classFE2CXM[i] = 'PS > Median'
      
    }
  }

  for (i in 1:len){
    if(model.params$Median.Ktrans.2CXM[i] < median.Ktrans.2CXM){
      classKtrans2CXM[i] = 'Ktrans < Median'
    } 
    else{
      classKtrans2CXM[i] = 'Ktrans > Median'
    }
  }
  # pairs(classKtrans2CXM, classFp2CXM, classFE2CXM)
  
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

# Define plotting parameters
bottomline <- 7
topline <- 10
mai = c(2.0,3.0,1.5,1.5)

if (param =='Fp')
  
{
  labelsforLeg1 = '\u2265 group median'
  labelsforLeg2 = '< group median'
  predictor = as.numeric(factor(classFp2CXM))
  titleForPlot = expression(paste(italic(F)[p]))
  location = 'bottomleft'
  ones <-which(predictor == 2)
  predictor[ones] = 0
  atrisklabs = c('\u2265 group median', '< group median')
  pvalheight = -50.
  atrisk.line = c(bottomline,topline)
# To univariate cox to compute HR and p value
  cox1 <- coxph(Surv(DFStime,status) ~ classFp2CXM, method = 'breslow')
  print(cox1)
  print(exp(confint(cox1)))
  }

if (param =='FE')
{
  labelsforLeg1 = '\u2265 group median'
  labelsforLeg2 = '< group median'
  predictor = as.numeric(factor(classFE2CXM))
  titleForPlot = expression(paste(italic(PS)))
  location = 'bottomleft'
  ones <-which(predictor == 2)
  predictor[ones] = 0
  atrisklabs = c('\u2265 group median', '< group median')
  pvalheight = -50.
  atrisk.line = c(bottomline,topline)

  cox1 <- coxph(Surv(DFStime,status) ~ classFE2CXM, method = 'breslow')
  print(cox1)
  print(exp(confint(cox1)))
  }

if (param =='Ktrans')
{
  labelsforLeg1 = '\u2265 group median'
  labelsforLeg2 = '< group median'
  predictor = as.numeric(factor(classKtrans2CXM))
  titleForPlot = expression(paste(italic(K)^{trans}))
  location = 'bottomleft'
  ones <-which(predictor == 2)
  predictor[ones] = 0
  atrisklabs =c('\u2265 group median', '< group median')
  pvalheight = -50.
  atrisk.line = c(bottomline,topline)
  cox1 <- coxph(Surv(DFStime,status) ~ classKtrans2CXM, method = 'breslow')
print(cox1)
print(exp(confint(cox1)))
  }

if (param =='2CXMvp'){
  labelsforLeg1 = '\u2265 group median'
  labelsforLeg2 = '< group median'
  predictor = as.numeric(factor(classvp2CXM))
  titleForPlot = expression(paste(italic(v)[p]))
  location = 'bottomleft'
  ones <-which(predictor == 2)
  predictor[ones] = 0
  atrisklabs = c('\u2265 group median', '< group median')
  pvalheight = -50.
  atrisk.line = c(bottomline,topline)
  cox1 <- coxph(Surv(DFStime,status) ~ classvp2CXM, method = 'breslow')
print(cox1)
print(exp(confint(cox1)))


  }

if (param =='2CXMve'){
  labelsforLeg1 = '\u2265 group median'
  labelsforLeg2 = '< group median'
  predictor = as.numeric(factor(classve2CXM))
  titleForPlot = expression(paste(italic(v)[e]))
  location = 'bottomleft'
  ones <-which(predictor == 2)
  predictor[ones] = 0
  atrisklabs = c('\u2265 group median', '< group median')
  pvalheight = -50.
  atrisk.line = c(bottomline,topline)
  cox1 <- coxph(Surv(DFStime,status) ~ classve2CXM, method = 'breslow')
  print(cox1)
  print(exp(confint(cox1)))
  
  }

if (param =='KtransT')
{
  labelsforLeg1 = '\u2265 group median'
  labelsforLeg2 = '< group median'
  predictor = as.numeric(factor(classKtransTOFTS))
  titleForPlot = expression(paste('Tofts  ',italic(K)^{trans}))
  location = 'bottomleft'
  ones <-which(predictor == 2)
  predictor[ones] = 0
  atrisklabs = c('\u2265 group median', '< group median')
  pvalheight = -50.
  atrisk.line = c(bottomline,topline)

  cox1 <- coxph(Surv(DFStime,status) ~ classKtransTOFTS, method = 'breslow')
}

if (param =='KtransET')
{
  labelsforLeg1 = '\u2265 group median'
  labelsforLeg2 = '< group median'
  predictor = as.numeric(factor(classKtransETOFTS))
  titleForPlot = expression(paste('Extended Tofts  ',italic(K)^{trans}))
  location = 'bottomleft'
  ones <-which(predictor == 2)
  predictor[ones] = 0
  atrisklabs = c('\u2265 group median', '<  group median')
  pvalheight = -50.
  atrisk.line = c(bottomline,topline)
  cox1 <- coxph(Surv(DFStime,status) ~ classKtransETOFTS, method = 'breslow')
}


if (param =='vol')
{
  labelsforLeg2 = '\u2265 group median'
  labelsforLeg1 = '< group median'
  predictor = as.numeric(factor(classVol))
  titleForPlot = expression('MRI volume')
  location = 'bottomleft'
  atrisklabs = c( '<  group median', '\u2265 group median')
  pvalheight = -50.
  atrisk.line = c(bottomline,topline)

  cox1 <- coxph(Surv(DFStime,status) ~ classVol, method = 'breslow')
print(cox1)
print(exp(confint(cox1)))
  }
# older age is bad
if (param =='age')
{
  labelsforLeg1 = '< group median'
  labelsforLeg2 = '\u2265 group median' 
  predictor = as.numeric(factor(classAge))
  titleForPlot = expression('Age')
  location = 'bottomleft'
  atrisklabs = c('<  group median', '\u2265 group median')
pvalheight = -50.
  atrisk.line = c(bottomline,topline)
  cox1 <- coxph(Surv(DFStime,status) ~ classAge, method = 'breslow')
print(cox1)
print(exp(confint(cox1)))

  }

# since we have opposite trend to TK params, labels are the other way around to get labels right. 
if (param =='stage')
{
  labelsforLeg1 = c('Stage T1b/T2')
  labelsforLeg2 = c('Stage T3/T4')
  predictor = as.numeric(factor(classStage))
  titleForPlot = expression('T stage')
  location = 'bottomright'
  atrisklabs = c('Stage T1b/T2','Stage T3/T4')
  pvalheight = -50.
  atrisk.line = c(bottomline,topline)
  cox1 <- coxph(Surv(DFStime,status) ~ classStage, method = 'breslow')
print(cox1)
print(exp(confint(cox1)))

}
# SCC histology is good, other is bad
if (param =='histology')
{
  labelsforLeg1 = c('SCC')
  labelsforLeg2 = c('Other')
  predictor = as.numeric(factor(classHistology))
  titleForPlot = expression('Histological subtype')
  location = 'bottomleft'
  ones <-which(predictor == 2)
  predictor[ones] = 0
  atrisklabs = c('    SCC','  Other')
  pvalheight = -50.
  atrisk.line = c(bottomline,topline)

    cox1 <- coxph(Surv(DFStime,status) ~ classHistology, method = 'breslow')
print(cox1)
print(exp(confint(cox1)))

    }

# Positive nodes is bad
if (param =='nodes')
{
  labelsforLeg1 = c('Negative')
  labelsforLeg2 = c('Positive')
  predictor = as.numeric(factor(classNodes))
  titleForPlot = expression('Nodal Status')
  location = 'bottomleft'
  atrisklabs = c('Negative','Positive')
  pvalheight = -50.
  atrisk.line = c(bottomline,topline)
  cox1 <- coxph(Surv(DFStime,status) ~ classNodes, method = 'breslow')
print(cox1)
print(exp(confint(cox1)))

  }

if (param =='treatment')
{
  labelsforLeg1 = c('CRT')
  labelsforLeg2 = c('RT')
  predictor = as.numeric(factor(classTreatment))
  titleForPlot = expression('Treatment')
  location = 'bottomleft'
  atrisklabs = c('  CRT','        RT')
  pvalheight = -50.
  atrisk.line = c(bottomline,topline)
   cox1 <- coxph(Surv(DFStime,status) ~ classTreatment, method = 'breslow')
print(cox1)
print(exp(confint(cox1)))

   }

# Make plot
ddist <- datadist(predictor)
options(datadist = 'ddist')

survobj = Surv(DFStime,status)

fits <- prodlim(Surv(DFStime,status) ~ predictor,
        exact = TRUE)
diff <-survdiff(Surv(DFStime,status) ~ predictor)

p <- par(mai = mai, xaxs = 'i',yaxs = 'i', ps = 10, las = 1, cex.axis = 4.5, cex.lab = 4.5, cex.main = 3, bty = 'o')

p <- plot(fits,
          border = T,
     axes = 'F',
     atrisk = TRUE,
     xlim = c(0,7),
     ylim = c(0,1),
     # legend.x=c(location),
     lwd = 5,
     lty = c(1,2,3,4),
     col = c(1,1,1,1),
     # legend.title="",
     # legend.y.intersp=1,
     # background.horizontal=seq(0,1,.20),
     # legend.cex=4,
     # legend.legend = labelsforLeg,
     legend = F,
     atRisk.title ='',
     atRisk.adj = c(1.4),
     atRisk.cex=2.2,
     marktime.cex = 2.6,
    atRisk.labels= atrisklabs,
      marktime = TRUE,
     atrisk.line = atrisk.line,
     xlab = NA, 
     ylab = NA,
     percent = T,
    logrank = FALSE,
    confint = FALSE)

legend(0.1, 0.25, labelsforLeg1, cex = 4, fill = F, border = F, lty = 1, lwd = 3, bty = "n")
legend(0.1, 0.15, labelsforLeg2, cex = 4, fill = F, border = F, lty = 2, lwd = 3, bty = "n")

title(titleForPlot,cex.main = 8.0)

axis(side = 1,lwd = 3, line = -4.3, labels = FALSE, tick = FALSE)
axis(side = 1,lwd = 3, line = -4.3, labels = FALSE)
axis(side = 1,lwd = 3, line = -3.3, tick = FALSE)
axis(side = 2, lwd = 3, line = -2, tck = 0, labels = c(0,20,40,60,80,100), at = c(0,0.2,0.4,0.6,0.8,1))

mtext(side = 1, "Years", line = 2.0, cex = 4.5)
mtext(side = 2, "Percentage Recurrence Free", line = 6.4, cex = 4.5, las = 0)

coeffs <- coef(summary(cox1))
pval = as.matrix(coeffs[,5])
mtext(side = 1, paste('P-value',' =',format(signif(pval,2),nsmall = 2)), at = c(5), line = pvalheight - 3, cex = 5.0)

# Return plot object
return(p)

}
