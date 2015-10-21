# batchKM.R 
# runs Kaplan-Meier analysis for all clinicopathologic and tracer kinetic variables and makes Figures 3 and 4 

# Make Figure 3 
 pdf( "Figures/Figure3.pdf", width = 30, height = 20)
 par(mfrow=c(2,2), cex.lab = 1.5,cex.axis = 1.5)
KMplotter('KtransT')
KMplotter('KtransET')
KMplotter('Fp')
KMplotter('FE')
 dev.off()

#Make Figure 4 
 pdf( "Figures/Figure4.pdf", width = 30, height = 30)
 par(mfrow=c(3,2))
KMplotter('stage')
KMplotter('treatment')
KMplotter('histology')
KMplotter('nodes')
KMplotter('age')
KMplotter('vol')
 dev.off()


