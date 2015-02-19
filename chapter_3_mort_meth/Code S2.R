# Code to fit a Siler (Gompertz bathtub) mortality model to re-sighting data 
# of a species with male natal dispersal

# 07 February 2015
# Comments:
# - to start set working directory to the folder that contains simFcts. R, 
#   plotTraces.R, plotMortRate.R, and plotDispPdf.R
# - runs the simulation analysis as described in the manuscript
# - in a loop from 1:12, this scripts sources Code S1.R to simulate data of varying
#   initial birth cohort sample sizes, proportions of individuals that died 
#   before they could be sexed, and proportions of observed deaths (see
#   manuscript for details)
# - Code S1.R and this file (Code S2.R) both sources a function file (simFcts.R)
# - for illustration of results plotTraces.R plots the traces,
#   plotMortRate.R plots the mortality rates for males (blue) and females (pink),
#   and plotDispPdf.R plots the pdf of age at dispersal for males

rm(list = ls())
setwd("...")

library(msm)
library(RColorBrewer)

for (i in 1:2) {
source("Code S1.R")
source("simFcts.R")  # source functions


# Define priors and starting values
thetaPriorMean <- matrix(c(-3, 0.2, 0, -4, 0.01), 2, 5, byrow = TRUE, 
                         dimnames = list(c("f", "m"), 
                                          c("a0", "a1", "c", "b0", "b1")))
thetaPriorSd <- matrix(rep(c(0.5, 0.25, 0.25, 0.5, 0.25), each = 2), 2, 5,  
                       dimnames = dimnames(thetaPriorMean))
thetaStart <-  matrix(rep(c(-2, 0.2, 0, -5, 0.01), 2), 2, 5, byrow = TRUE)
thetaLow <- thetaPriorSd * 0
thetaLow[, c(1, 4)] <- -Inf
nthe <- 10

gammaPrior1 <- c(8, 2)
gammaPrior2 <- c(2, 1)
gamStart <- c(4, 0.1)

# ages
xl <- dat$xli
xl[indObsDi == 1] <- dat$xi[indObsDi == 1] 
xl[which(xl == 0)] <- 0.000001
xt <- rep(0, ni)
xt[indIm == 1] <- dat$xdi[indIm == 1]

# dispersal indicators and indices
idPotEm <- which(covs[, 'm'] == 1 & indIm == 0 & unsexed == 0 & 
                   xl > minDx & indObsDi == 0)
idIm <- which("indIm" == 1)
dispStateNow <- rep(0, nrow(dat))
dispStateNow[idPotEm] <- rbinom(length(idPotEm), 1, 0.5)
idEmNow <- which(dispStateNow == 1)

# sex indicators and indices 
probFem <- 0.5
if (i %in% seq(1, 11, 2)) {
  unsexed <- rep(0, length(unsexed)) 
  sexInd <- "ks"} else {  # known sexes
  sexInd <- "us"}  # unknown sexes
idNoSex <- which(unsexed == 1)
sex <- rbinom(length(idNoSex), 1, probFem)  # success = female
covsNow <- covs
covsNow[idNoSex, ][covsNow[idNoSex, 1] != sex, ] <- 
  abs(covsNow[idNoSex, ][covsNow[idNoSex, 1] != sex, ] - 1)

# run model
#niter <-15000
#burn <- 5000
#thin <- 20
#keep <- seq(burn, niter, thin)

niter <- 1000
burn <- 200
thin <- 20
keep <- seq(burn, niter, thin)

nsim <- 4
ncpus <- 4
require(snowfall)
sfInit(parallel = TRUE, cpus = ncpus)
sfExport(list = c(ls(), ".Random.seed"))
sfLibrary(msm, warn.conflicts = FALSE)
outParallel <- sfClusterApplyLB(1:nsim, RunMCMC)
sfStop() 

out <- ExtractParalOut(outParallel)
quantList <- CalcDemoQuant(out)

spOut <- list(theta = out$theta, coeffs = out$coeffs, mort = quantList$mort,
              surv = quantList$surv, pdf = quantList$pdf,
              cuts99 = quantList$cuts99, cuts95 = quantList$cuts95,
              x = quantList$x, gam = out$gam, disp = quantList$disp)

# plot results
# traces
source("plotTraces.R")
# mortality rate
source("plotMortRate.R")
# pdf ages at dispersal
source("/Users/Viktualia/Dropbox/Projects/008_LionSexDiffMort/code/plotBuildingBlocks/dispPdfFig001.R")
}
system("say Ich bin fertig")   
