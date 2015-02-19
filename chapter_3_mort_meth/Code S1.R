# 07 Feb 2015
# Comments: 
# - this script gets sourced by Code S2.R, see Code S2.R for further details
# - code to simulate demographic census data similar to the ones 
#   from wild lion populations

source("simFcts.R")  # source functions

# initial birth cohort sample size
if (i <= 6) {n <- 500} else {n <- 2000}
# minimum age at dispersal
minDx <- 1.75
# proportion observed deaths
propDeathVec <- rep(c(0.01, 0.01, 0.05, 0.05, 0.1, 0.1), 2)
propDeath <- propDeathVec[i]
# proportion of individuals that remained unsexed at age of death
propUnsex <- 0.3
# female Siler mortality parameters
thf <- matrix(c(-1.4, 0.65, 0.07, -3.8, 0.2), 1, 5, 
              dimnames = list("f", c("a0", "a1", "c", "b0", "b1")))
# male Siler mortality parameters
thm <- matrix(c(-1.2, 0.7, 0.16, -3.5, 0.23), 1, 5, 
              dimnames = list("m", c("a0", "a1", "c", "b0", "b1")))
thetaReal <- matrix(c(thf, thm), 2, 5, byrow = TRUE, 
                    dimnames = list(c("f", "m"), 
                                    c("a0", "a1", "c", "b0", "b1")))
# vector of ages
xv <- seq(0, 25, 0.1)
# cdf of age at death (f: females, m: males)
Fxf <- 1 - CalcSurv(thf, xv)
Fxm <- 1 - CalcSurv(thm, xv)
# sex ratio
sr <- 0.5
# draw sex indicator
sexf <- rev(sort(rbinom(n, 1, sr)))
# number of females
nf <- sum(sexf)
# number of males
nm <- n - nf
# draw ages at death from cdf of age at death
xi <- n * 0
xi[sexf == 1] <- xv[findInterval(runif(nf), Fxf)]
xi[sexf == 0] <- xv[findInterval(runif(nm), Fxm)]

# gamma parameters of age at disperal for males
gammaPars <- c(10, 4)

# residents dispersal age
xdi <- minDx + rgamma(n, gammaPars[1], gammaPars[2])
# disperser index
idDips <- which(xdi < xi & sexf == 0)

# potential immigrant male age at death
xi2 <- xv[findInterval(runif(nm), Fxm)]
# potential immigrant male dispersal age
xdi2 <- minDx + rgamma(nm, gammaPars[1], gammaPars[2])
# indicator of actual immigrant
indImTemp <- which(xdi2 < xi2)
# immigrant male age at death
xi2 <- xi2[indImTemp]
# immigrant male age at immigration
xdi2 <- xdi2[indImTemp]
# number of immigrants
nm2 <- length(xi2)

# build data frame
Xi <- c(xi, xi2)
Xdi <- c(xdi, xdi2)
datRaw <- data.frame(Xi,  Xdi)

# indeces
# immigrant
indIm <- c(rep(0, n), rep(1, nm2))
# female
sex <- c(sexf, rep(0, nm2))
# disperser
indDisp <- rep(0, n + nm2)
indDisp[idDips] <- 1
# died before sexed
unsexed <- rep(0, n + nm2)
unsexed[which(Xi <= 1)] <- rbinom(length(which(Xi<= 1)), 1, propUnsex)
# observed death
indObsDi <- rbinom(length(Xi), 1, propDeath)
# indicator unobserved death 
idNoDeath <- which(indObsDi == 0)
# covariate object sex
covs <- cbind(sex, 1 - sex)
colnames(covs) <- c("f", "m")

# final data frame
# last seen age
Xli <- Xi
Xli[indDisp == 1] <- Xdi[indDisp == 1]
dat <- data.frame(xli = Xli, xdi = Xdi, xi = Xi)
ni <- nrow(dat)
rm(list = setdiff(ls(), c("dat", "ni", "covs", "minDx", "indObsDi",
                          "indIm", "unsexed", "i", "propDeath", "thetaReal")))
