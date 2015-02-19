# Date: 07 Feb 2015
# Author: Julia Barthold
# Comment: 
# - this script gets sourced by Code S2.R, see Code S2.R for further details
# - plots the dispersal pdf from the simulation

op <- par(oma = c(0,0,0,0), tcl = -0.25, cex.axis = 0.75, mar = c(3.6,3.6,1.1,0.1),
          mfrow = c(1,1))

# plot
temp = spOut

# cut objects to only the relevant ages
temp$x <- temp$x[temp$x <= 10]
temp$disp[[1]] <- temp$disp[[1]][ , 1:length(temp$x)]
upperXLim <- 10
upperYLim <- 0.6
intY <- 0.1
intX <- 2

plot(temp$x, temp$dips[[1]][1, ], type = "n",
     ylab = "", 
     ylim = c(0, seq(0, upperYLim, intY)[length(seq(0, upperYLim, intY))]),
     bty = "n", axes = FALSE, 
     xlim = c(0,seq(intX-startAge, upperXLim - startAge, 
                    intX)[length(seq(intX-startAge, 
                                     upperXLim - startAge, intX))]),
     xlab = "")

#Axes
axis(1, at = c(0, seq(intX-startAge, upperXLim - startAge, intX)), 
     labels = c(0, seq(intX-startAge, upperXLim - startAge, intX)))
axis(2, las = 1, at=seq(0, upperYLim, intY), 
     labels = format(seq(0, upperYLim, intY)))
abline(h= seq(0, upperYLim, intY),col="#e5e5e560")

# pdf dispersal
polygon(c(temp$x, rev(temp$x)), c(temp$disp[[1]][2, ], rev(temp$disp[[1]][3, ])),
        col = "white", border = NA)
polygon(c(temp$x, rev(temp$x)), c(temp$disp[[1]][2, ], rev(temp$disp[[1]][3, ])),
        col = adjustcolor(mcol, alpha.f= 0.20), border = NA)
lines(temp$x, temp$disp[[1]][2, ], col = mcol, lwd = 0.4)
lines(temp$x, temp$disp[[1]][3, ], col = mcol, lwd = 0.4)

# Real dispersal rates
if(propDeath %in% ls()) {
lines(temp$x, dgamma(temp$x, lambda[1], lambda[2]), col = mcol)
}
mtext("pdf of ages at dipsersal",
      side = 2, outer = FALSE, cex = 0.8, line = 2.5,col="black")
mtext("Age (years)", side = 1, outer = FALSE, cex = 0.8, line = 2.5, col = "black")
if("propDeath" %in% ls()) {
  mtext(paste("Proportion known death:", propDeath, ""), side = 3, outer = FALSE,
        cex = 1, line = 0, col = "black")
}

par <- op
