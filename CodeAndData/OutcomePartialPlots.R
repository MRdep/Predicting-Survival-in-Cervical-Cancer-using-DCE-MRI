#Plots predicted DFS plots

OutcomePartialPlots <- function(Low,High, outcome,names, xlab){
  
ylim = c(.20,.70)
ylab = '5-year DFS probability'

boxplot(data.frame(Low, High), ylim = ylim, ylab = ylab, xlab = xlab, names =  names)

}