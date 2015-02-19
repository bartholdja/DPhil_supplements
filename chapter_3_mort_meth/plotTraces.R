# Date: 07 Feb 2015
# Author: Julia Barthold
# Comment: 
# - this script gets sourced by Code S2.R, see Code S2.R for further details
# - plots the traces for the simulation

tracesPlot <- do.call('rbind', lapply(1:length(outParallel), function(x){
  outParallel[[x]]$theta
}))
tracesLamPlot <- do.call('rbind', lapply(1:length(outParallel), function(x){
  outParallel[[x]]$gamma
}))

indKeep <- rep(0, nrow(outParallel[[1]]$theta))
indKeep[keep] <- 1
indKeepNonPar <- rep(indKeep, nsim)

perc <- 0.5

par(mfrow = c(2, 10), mar = c(2, 2, 1, 1))
for (i in 1:ncol(outParallel[[1]]$the)) {
  for (j in 1:length(outParallel)) {
    if (j == 1) {
        ylimRange <- NULL
        if(range(tracesPlot[indKeepNonPar == 1, i])[1] < 0 & 
             range(tracesPlot[indKeepNonPar == 1, i])[2] < 0) {
          ylimRange[1] <- range(tracesPlot[indKeepNonPar == 1, i])[1] * (1 + perc)
          ylimRange[2] <- range(tracesPlot[indKeepNonPar == 1, i])[2] * (1 - perc)
        }
        if(range(tracesPlot[indKeepNonPar == 1, i])[1] < 0 & 
             range(tracesPlot[indKeepNonPar == 1, i])[2] > 0) {
          ylimRange <- range(tracesPlot[indKeepNonPar == 1, i]) * (1 + perc)
        }
        if(range(tracesPlot[indKeepNonPar == 1, i])[1] > 0 & 
             range(tracesPlot[indKeepNonPar == 1, i])[2] > 0) {
          ylimRange[1] <- range(tracesPlot[indKeepNonPar == 1, i])[1] * (1 - perc)
          ylimRange[2] <- range(tracesPlot[indKeepNonPar == 1, i])[2] * (1 + perc)
        }
                
      plot(outParallel[[j]]$theta[keep, i], type = 'l', 
           main = outParallel[[1]]$names[i], 
           ylim = ylimRange,
           col = brewer.pal(ncol(thetaStart), "BrBG")[-3][j], ylab = "")
    } else {
      lines(outParallel[[j]]$theta[keep, i], type = 'l', 
            main = outParallel[[1]]$names[i], 
            col = brewer.pal(ncol(thetaStart), "BrBG")[-3][j])
    } 
  }
}
for(i in 1:2) {
  for (j in 1:length(outParallel)){
    if (j == 1) {
      ylimRangeLam <- NULL
      if(range(tracesLamPlot[indKeepNonPar == 1, i])[1] < 0 & 
           range(tracesLamPlot[indKeepNonPar == 1, i])[2] < 0) {
        ylimRangeLam[1] <- range(tracesLamPlot[indKeepNonPar == 1, i])[1] * (1 + perc)
        ylimRangeLam[2] <- range(tracesLamPlot[indKeepNonPar == 1, i])[2] * (1 - perc)
      }
      if(range(tracesLamPlot[indKeepNonPar == 1, i])[1] < 0 & 
           range(tracesLamPlot[indKeepNonPar == 1, i])[2] > 0) {
        ylimRangeLam <- range(tracesLamPlot[indKeepNonPar == 1, i]) * (1 + perc)
      }
      if(range(tracesLamPlot[indKeepNonPar == 1, i])[1] > 0 & 
           range(tracesLamPlot[indKeepNonPar == 1, i])[2] > 0) {
        ylimRangeLam[1] <- range(tracesLamPlot[indKeepNonPar == 1, i])[1] * (1 - perc)
        ylimRangeLam[2] <- range(tracesLamPlot[indKeepNonPar == 1, i])[2] * (1 + perc)
      }
      
      plot(outParallel[[j]]$gamma[keep, i], type = 'l', main = 
             colnames(outParallel[[1]]$gamma)[i], 
           ylim = ylimRangeLam,
           ylab = "")
    } else {
      lines(outParallel[[j]]$gamma[keep, i], main = 
              colnames(outParallel[[1]]$gamma)[i],
            col = brewer.pal(ncol(thetaStart), "BrBG")[-3][j])
    }
  }
}

