# Date: 07 Feb 2015
# Author: Julia Barthold
# Comment: 
# - this script gets sourced by Code S2.R, see Code S2.R for further details
# - plots the mortality rate

source("simFcts.R")
# Mortality
# setting up some objects
startAge <- 0 
cutOff <- 0.05  # plot objects until (1- cutOff) of synthetic cohort are dead

# color
fcol = "#d7191c"
mcol = "#2b83ba"

op <- par(oma = c(0,0,0,0), tcl = -0.25, cex.axis = 0.75, mar = c(3.6,3.6,1.1,0.1),
          mfrow = c(1,1))

# plot
temp = spOut

# Subset to the required data (e.g. birth to 95% mortality etc.)
maleStart <- 1
femaleStart <- 1

if(cutOff == 0.05) maleEnd = temp$cuts95$m; femaleEnd = temp$cuts95$f
if(cutOff == 0.01) maleEnd = temp$cuts99$m; femaleEnd = temp$cuts99$f

# cut objects to only the relevant ages
temp$xm <- temp$x[maleStart:maleEnd]
temp$xf <- temp$x[femaleStart:femaleEnd]
temp$xm <- temp$xm-min(temp$xm)  # This starts age from 0 if not running from birth.
temp$xf <- temp$xf-min(temp$xf)  # This starts age from 0 if not running from birth.
temp$mort$m <- temp$mort$m[ , maleStart:maleEnd]
temp$mort$f <- temp$mort$f[ , femaleStart:femaleEnd]
temp$surv$m <- temp$surv$m[ , maleStart:maleEnd]
temp$surv$f <- temp$surv$f[ , femaleStart:femaleEnd]
temp$pdf$m <- temp$pdf$m[grep("ave", rownames(temp$pdf$m)), maleStart:maleEnd]
temp$pdf$f <- temp$pdf$f[grep("ave", rownames(temp$pdf$f)), femaleStart:femaleEnd]

upperXLim <- 20
upperYLim <- 1
intY <- 0.2
intX <- 2

plot(temp$xf, temp$mort$f[1, ], type = "n",
     ylab = "", 
     ylim = c(0, seq(0, upperYLim, intY)[length(seq(0, upperYLim, intY))]),
     bty = "n", axes = FALSE, 
     xlim = c(0,seq(intX-startAge, upperXLim - startAge, 
                    intX)[length(seq(intX-startAge, 
                                     upperXLim - startAge, intX))]),
     xlab = "")

# Axes
axis(1, at = c(0, seq(intX-startAge, upperXLim - startAge, intX)), 
     labels = c(0, seq(intX-startAge, upperXLim - startAge, intX)))
axis(2, las = 1, at=seq(0, upperYLim, intY), 
     labels = format(seq(0, upperYLim, intY)))
abline(h= seq(0, upperYLim, intY),col="#e5e5e560")

# Mortality trajectories
# females
polygon(c(temp$xf, rev(temp$xf)), c(temp$mort$f[2, ], rev(temp$mort$f[3, ])),
        col = "white", border = NA)
polygon(c(temp$xf, rev(temp$xf)), c(temp$mort$f[2, ], rev(temp$mort$f[3, ])),
        col = adjustcolor(fcol, alpha.f= 0.20), border = NA)
lines(temp$xf, temp$mort$f[2, ], col = fcol, lwd = 0.4)
lines(temp$xf, temp$mort$f[3, ], col = fcol, lwd = 0.4)

# males
polygon(c(temp$xm, rev(temp$xm)), c(temp$mort$m[2, ], rev(temp$mort$m[3, ])),
        col = "white", border = NA)
polygon(c(temp$xm, rev(temp$xm)), c(temp$mort$m[2, ], rev(temp$mort$m[3, ])),
        col = adjustcolor(mcol, alpha.f= 0.50), border = NA)

lines(temp$xm, temp$mort$m [2, ], col = mcol, lwd = 0.4)
lines(temp$xm, temp$mort$m[3, ], col = mcol, lwd = 0.4)

if(propDeath %in% ls()) {
# Real mortality rates
lines(seq(0, 20, 0.1)[1:femaleEnd], CalcMort(thetaReal[1, ], 
                                             seq(0, 20, 0.1)[1:femaleEnd]), col = fcol)
lines(seq(0, 20, 0.1)[1:maleEnd], CalcMort(thetaReal[2, ],
                                           seq(0, 20, 0.1)[1:maleEnd]), col = mcol)
}
mtext(expression(paste("Mortality rate (", mu[x], ")")), 
      side = 2, outer = FALSE, cex = 0.8, line = 2.5,col="black")
mtext("Age (years)", side = 1, outer = FALSE, cex = 0.8, line = 2.5, col = "black")
if("propDeath" %in% ls()) {
mtext(paste("Proportion known death:", propDeath, ""), side = 3, outer = FALSE,
      cex = 1, line = 0, col = "black")
}
par <- op
