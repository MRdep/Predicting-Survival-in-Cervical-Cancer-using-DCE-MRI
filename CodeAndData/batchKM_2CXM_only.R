# batchKM.R 
source('KMplotter.R')
# runs Kaplan-Meier analysis for all clinicopathologic and tracer kinetic variables and makes Figures 3 and 4 

# Make Figure 3 
tiff( "Figure2.tiff", width = 15, height = 56, units = 'in', res = 400)
 
par(mfrow=c(5,1), cex.lab = 1.5,cex.axis = 1.5)
KMplotter('Fp')
KMplotter('FE')
KMplotter('Ktrans')
KMplotter('2CXMvp')
KMplotter('2CXMve')
dev.off()

#Make Figure 4 
tiff( "SupplementaryFigure2.tiff", width = 24*1.5, height = 24*1.5, units = 'in', res = 400)
par(mfrow=c(3,2))
KMplotter('stage')
KMplotter('treatment')
KMplotter('histology')
KMplotter('nodes')
KMplotter('age')
KMplotter('vol')
dev.off()

