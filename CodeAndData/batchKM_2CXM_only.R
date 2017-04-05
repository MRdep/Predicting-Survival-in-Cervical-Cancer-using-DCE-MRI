# batchKM.R 
source('KMplotter.R')
# runs Kaplan-Meier analysis for all clinicopathologic and tracer kinetic variables and makes Figures 3 and 4 

# Make Figure 1
png( "Figure1.png", width = (30), height = (32.6), units = 'in', res = 400)
par(mfrow=c(3,2), cex.lab = 1.5,cex.axis = 1.5)
KMplotter('treatment')
KMplotter('nodes')
KMplotter('age')
KMplotter('vol')
KMplotter('Fp')
KMplotter('Ktrans')
dev.off()

#Make Figure Supp Figure 2
png( "SupplementaryFigure3.png",width = 15, height = 54.8, units = 'in', res = 400)
par(mfrow=c(5,1))
KMplotter('stage')
KMplotter('histology')
KMplotter('FE')
KMplotter('2CXMvp')
KMplotter('2CXMve')

dev.off()

