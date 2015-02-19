# Comments: contains function to fit the Siler mortality model
# x: last seen age, xt: truncation age, xd: age at dispersal

# Functions:
# mortality rate
CalcMort <- function(th, ...) UseMethod("CalcMort")

CalcMort.matrix <- function(th, x) {
  exp(th[, 1] - th[, 2] * x) + th[, 3] + exp(th[, 4] + th[, 5] * x)
}

CalcMort.numeric <- function(th, x) {
  exp(th[1] - th[2] * x) + th[3] + exp(th[4] + th[5] * x)
}

# survivorship
CalcSurv <- function(th, ...) UseMethod("CalcSurv")

CalcSurv.matrix <- function(th, x) {
  exp(exp(th[, 1])/th[, 2] * (exp(-th[, 2] * x) - 1) - th[, 3] * x + 
        exp(th[, 4])/th[, 5] * (1 - exp(th[, 5] * x)))
}

CalcSurv.numeric <- function(th, x) {
  exp(exp(th[1])/th[2] * (exp(-th[2] * x) - 1) - th[3] * x + 
        exp(th[4])/th[5] * (1 - exp(th[5] * x)))
}

# fill covariate  matrices
CalcCovPars <- function(pars, covs) {
  covs %*% pars
}

# mortality parameter likelihood
CalcMortPdf <- function(thCov, x) {
  log(CalcMort(thCov, x) * CalcSurv(thCov, x))
}

# mortality parameter posterior
CalcMortPost <- function(th, thCov, mortPdf, idEm) {
  mortPost <- mortPdf
  mortPost[idEm] <- log(CalcSurv(thCov[idEm, ], xl[idEm]))
  sum(mortPost - log(CalcSurv(thCov, xt))) + 
    sum(dtnorm(c(th), c(thetaPriorMean), c(thetaPriorSd),
               low = c(thetaLow), log = TRUE))
}

# age at death posterior
CalcAgePost <- function(thCov, mortPdf, idEm) {
  agePost <- mortPdf + CalcMortPdf(thetaPriorMean[1, ], xl)
  agePost[idEm] <- log(CalcSurv(thCov[idEm, ], xl[idEm])) + 
    log(CalcSurv(thetaPriorMean[1, ], xl[idEm]))
  return(agePost)
}

# dispersal parameter likelihood
CalcDispLike <- function(gam, idEm, thCov) { 
  dispLike <- xl * 0 
  # Immigrants
  dispLike[idIm] <- dgamma(xt[idIm] - minDx, gam[1], gam[2], log = TRUE) - 
    log(CalcSurv(thCov[idIm, ], xt[idIm]))
  # Non-dispersers
  idNotEm <- idPotEm[!(idPotEm %in% idEm)]
  dispLike[idNotEm] <- pgamma(xl[idNotEm] - minDx, gam[1], gam[2], 
                              lower.tail = FALSE, log = TRUE)
  # Dispersers
  dispLike[idEm] <- dgamma(xl[idEm] - minDx, gam[1], gam[2], log = TRUE)
  return(dispLike)
}

# dispersal parameter posterior
CalcGamPost <- function(gam, dispLike) {
  sum(dispLike) + 
    sum(dtnorm(gam, gammaPrior1, gammaPrior2, log = TRUE, 
               low = 0))  
}  

# dispersal state prior
CalcDispPrior <- function(idEm) {
  dispPr <- xl * 0
  dispPr[idEm] <- dgamma(xl[idEm] - minDx, gammaPrior1[1], gammaPrior1[2], 
           log = TRUE)
  idNotEm <- idPotEm[!(idPotEm %in% idEm)]
  dispPr[idNotEm] <- pgamma(xl[idNotEm] - minDx, gammaPrior1[1], gammaPrior1[2], 
           log = TRUE, lower.tail = FALSE)
  return(dispPr)
}

# dispersal state posterior
CalcDispPost <- function(dispLike, idEm) {
  dispPost <- dispLike +  CalcDispPrior(idEm)
  return(dispPost)
}

# sex prior 
CalcSexPrior <- function(covs) {
  (covs[, 1]) * probFem + (1 - covs[, 1]) * (1 - probFem)
}
# sex posterior
CalcSexPost <- function(covs, mortPdf) {
  mortPdf + log(CalcSexPrior(covs))
}

# function to update jumps
UpdateJumps <- function(jumps, updMat, iter, iterUpd, updTarg) {
  updRate <- apply(updMat[iter - ((iterUpd - 1):0), ], 2, sum) / iterUpd  
  updRate[updRate == 0] <- 1e-2
  jumps <- jumps * 
    matrix(updRate, nrow(jumps), ncol(jumps)) / updTarg
  return(jumps)
}

# the MCMC
RunMCMC <- function(sim) {
  theNow <- thetaStart
  gamNow <- gamStart
  if (sim > 1) {
    rm(".Random.seed", envir = .GlobalEnv); runif(1)
    theNow <- matrix(rtnorm(nthe, thetaStart[1:nthe], 0.1, 
                             low = thetaLow[1:nthe]), 2, 5) # adapt the low
  }
  if (sim > 1) {
    rm(".Random.seed", envir = .GlobalEnv); runif(1)
    gamNow[1:2] <- rtnorm(2, mean = gamStart[1:2], sd = 0.1, 
                             low = 0)
  }
  theCovNow <- CalcCovPars(theNow, covsNow)
  mortPdfNow <- CalcMortPdf(theCovNow, xl)
  mortPostNow <- CalcMortPost(theNow, theCovNow, mortPdfNow, idEmNow)
  sexPostNow <- CalcSexPost(covsNow, mortPdfNow)
  dispLikeNow <- CalcDispLike(gamNow, idEmNow, theCovNow)
  gamPostNow <- CalcGamPost(gamNow, dispLikeNow)
  dispPostNow <- CalcDispPost(dispLikeNow, idEmNow)
  agePostNow <- CalcAgePost(theCovNow, mortPdfNow, idEmNow)
  postNow <- dispPostNow + agePostNow
  
  # storing output
  namesThe <- paste(rep(colnames(thetaPriorMean), each = 2), 
                    rep(rownames(thetaPriorMean), 5), 
                    sep = '.')
  theMat <- matrix(0, niter, nthe, dimnames = list(NULL, namesThe))
  theUpdMat <- theMat
  theMat[1, ] <- c(theNow)
  gamMat <- matrix(0, niter, 2, dimnames = list(NULL, c("gam1", "gam2")))
  gamUpdMat <- gamMat
  gamMat[1, ] <- gamNow
  theJumps <- matrix(rep(0.1, nthe), 2, 5, dimnames = dimnames(theNow))
  gamJumps <- matrix(c(0.01, 0.01), 1, 2)
  iterUpd <- 100
  updTarg <- 0.25
  UpdJumps <- TRUE
  dispMat <- dispStateNow[idPotEm]
  sexMat <- matrix(sex, 1, length(idNoSex))
  
  # MCMC:
  Start <- Sys.time()
  for (iter in 2:niter) {
    # 1. Update mortality parameters:
    for (pp in 1:nthe) {
      theNew <- theNow
      theNew[pp] <- rtnorm(1, theNow[pp], theJumps[pp], low = thetaLow[pp])
      theCovNew <- CalcCovPars(theNew, covsNow)
      mortPdfNew <- CalcMortPdf(theCovNew, xl)
      mortPostNew <- CalcMortPost(theNew, theCovNew, mortPdfNew, idEmNow)
      postRatio <- exp(mortPostNew - mortPostNow)
      if (!is.na(postRatio) & postRatio > runif(1)) {
        theNow <- theNew
        theCovNow <- theCovNew
        mortPdfNow <- mortPdfNew
        mortPostNow <- mortPostNew
        if (UpdJumps) theUpdMat[iter, pp] <- 1
      }
    }
    if (any(theUpdMat[iter, ] == 1)) {
      dispLikeNow <- CalcDispLike(gamNow, idEmNow, theCovNow)
      gamPostNow <- CalcGamPost(gamNow, dispLikeNow)
      agePostNow <- CalcAgePost(theCovNow, mortPdfNow, idEmNow)
    }
    
    # 2. Propose sex for unsexed individuals
    covsNew <- covsNow
    sexNew <- rbinom(length(idNoSex), 1, probFem)  # success = female
    covsNew[idNoSex, ][covsNew[idNoSex, 1] != sexNew, ] <- 
      abs(covsNew[idNoSex, ][covsNew[idNoSex, 1] != sexNew, ] - 1)
    theCovNew <- CalcCovPars(theNow, covsNew)
    mortPdfNew <- CalcMortPdf(theCovNew, xl)
    mortPostNew <- CalcMortPost(theNow, theCovNew, mortPdfNew, idEmNow)
    sexPostNew <- CalcSexPost(theCovNew, mortPdfNew)
    
    r <- exp(sexPostNew - sexPostNow)[idNoSex]
    
    idUpd2 <- idNoSex[r > runif(length(idNoSex))]
    idUpd2 <-idUpd2[!(is.na(idUpd2))]
    if (length(idUpd2) > 0 ) {
      mortPdfNow[idUpd2] <- mortPdfNew[idUpd2]
      covsNow[idUpd2, ] <- covsNew[idUpd2, ]
      theCovNow[idUpd2, ] <- theCovNew[idUpd2, ]
      sexPostNow[idUpd2] <- sexPostNew[idUpd2]
      sex <- covsNow[idNoSex, 1]
    }
    mortPostNow <- CalcMortPost(theNow, theCovNow, mortPdfNow, idEmNow)
    dispLikeNow <- CalcDispLike(gamNow,idEmNow, theCovNow)
    gamPostNow <- CalcGamPost(gamNow, dispLikeNow)
    dispPostNow <- CalcDispPost(dispLikeNow, idEmNow)
    agePostNow <- CalcAgePost(theCovNow, mortPdfNow, idEmNow)
    postNow <- dispPostNow + agePostNow
    
    # 3. update dispersal parameters:
    for (pp in 1:2) {
      gamNew <- gamNow
      gamNew[pp] <- rtnorm(1, gamNow[pp], gamJumps[pp], 
                           lower = c(0, 0)[pp])
      dispLikeNew <- CalcDispLike(gamNew,idEmNow, theCovNow)
      gamPostNew <- CalcGamPost(gamNew, dispLikeNew)
      r <- exp(gamPostNew - gamPostNow)
      if (!is.na(r)) {
        if (r > runif(1)) {
          gamNow <- gamNew
          dispLikeNow <- dispLikeNew
          gamPostNow <- gamPostNew
          if (UpdJumps) gamUpdMat[iter, pp] <- 1
        }
      }
    }  
    if (any(gamUpdMat[iter, ] == 1)) {
      dispPostNow <- CalcDispPost(dispLikeNow, idEmNow)
      postNow <- dispPostNow + agePostNow
    }
    
    # 4. Update dispersal state, ages at death and ages at dispersal:
    dispStateNew <- dispStateNow
    dispStateNew[idPotEm] <- rbinom(length(idPotEm), 1, 0.5)
    idEmNew <- which(dispStateNew == 1)
    dispLikeNew <- CalcDispLike(gamNow,idEmNew, theCovNow)
    dispPostNew <- CalcDispPost(dispLikeNew, idEmNew)
    agePostNew <- CalcAgePost(theCovNow, mortPdfNow, idEmNew)
    postNew <- dispPostNew + agePostNew
    r <- exp(postNew - postNow)
    z <- runif(ni)
    idUpd <- which(r > z)
    idUpd <-idUpd[!(is.na(idUpd))] 
    if (length(idUpd) > 0) {
      dispStateNow[idUpd] <- dispStateNew[idUpd]
      dispLikeNow[idUpd] <- dispLikeNew[idUpd]
      dispPostNow[idUpd] <- dispPostNew[idUpd]
      agePostNow[idUpd] <- agePostNew[idUpd]
      postNow[idUpd] <- postNew[idUpd]
    }
    idEmNow <- which(dispStateNow == 1)
    gamPostNow <- CalcGamPost(gamNow, dispLikeNow)
    mortPostNow <- CalcMortPost(theNow, theCovNow, mortPdfNow, idEmNow)
    
    # 5. Dynamic Metropolis to update jumps:
    if (UpdJumps) {
      if (is.element(iter/iterUpd,c(1:10))) {
        theJumps <- UpdateJumps(theJumps, theUpdMat, iter, iterUpd, updTarg)
        gamJumps <- UpdateJumps(gamJumps, gamUpdMat, iter, iterUpd, updTarg)
      }
    }
    
    # 6. Fill in the output matrices:
    theMat[iter, ] <- c(theNow)
    gamMat[iter, ] <- gamNow
    if (iter %in% keep) {
      dispMat <- rbind(dispMat, dispStateNow[idPotEm])
      sexMat <- rbind(sexMat, sex)
    }
  }
  End <- Sys.time()
  parMat <- cbind(theMat, gamMat)
  coeffs <- cbind(apply(parMat[keep, ], 2, mean), apply(parMat[keep, ], 2, sd), 
                  t(apply(parMat[keep, ], 2, quantile, c(0.025, 0.975))))
  colnames(coeffs) <- c("Mean", "SE", "2.5%", "97.5%")
  out <- list(theta = theMat, gamma = gamMat, disp = dispMat, 
              keep = keep, names = namesThe, potEm = idPotEm, coeffs = coeffs)
  return(out)
}

# function to extract thinned sequences from multiple 
# runs and calculate coefficients
ExtractParalOut <- function(out) {
  # joint output objects
  thetaOut <- do.call('rbind', lapply(1:length(out), function(x){
    out[[x]]$theta[keep, ]
  }))
  gamOut <- do.call('rbind', lapply(1:length(out), function(x){
    out[[x]]$gamma[keep, ]
  }))
  coeffs <- rbind(cbind(apply(thetaOut, 2, mean), apply(thetaOut, 2, sd), 
                        t(apply(thetaOut, 2, quantile, c(0.025, 0.975)))),
                  cbind(apply(gamOut, 2, mean), apply(gamOut, 2, sd), 
                        t(apply(gamOut, 2, quantile, c(0.025, 0.975)))))
  colnames(coeffs) <- c("Mean", "SE", "2.5%", "97.5%")
  parList <- list(theta = thetaOut, coeffs = coeffs, gam = gamOut)
  dimnames(parList$theta) <- list(NULL, paste(rep(colnames(thetaPriorMean), each = 2),
                                              rep(rownames(thetaPriorMean), 5), 
                                              sep = '.'))
  
  
  rownames(parList$coeffs)[1:10] <- colnames(parList$theta)
  return(parList)
}

CalcDemoQuant <- function(out, ...) {
  # Construct survival and mortality curves:
  dx <- 0.1
  xv <- seq(0, 100, dx)      
  
  mortList <- lapply(c("f", "m"), function(ss) {
    ids <- grep(ss, colnames(out$theta))
    mort <- apply(out$theta[, ids], 1, function(th) CalcMort(th, xv))
    mortave <- apply(mort, 1, mean)
    mortci <- apply(mort, 1, quantile, c(0.025, 0.975))
    mortfin <- rbind(mortave, mortci)
    return(mortfin)
  })
  
  survList <- lapply(c("f", "m"), function(ss) {
    ids <- ids <- grep(ss, colnames(out$theta))
    surv <- apply(out$theta[, ids], 1, function(th) CalcSurv(th, xv))
    survave <- apply(surv, 1, mean)
    survci <- apply(surv, 1, quantile, c(0.025, 0.975))
    survfin <- rbind(survave, survci)
    return(survfin)
  })
  
  pdfList <- lapply(c("f", "m"), function(ss) {
    ids <- ids <- grep(ss, colnames(out$theta))
    pdfm <- apply(out$theta[, ids], 1, 
                  function(th) CalcSurv(th, xv) * CalcMort(th, xv))
    pdfave <- apply(pdfm, 1, mean)
    pdfci <- apply(pdfm, 1, quantile, c(0.025, 0.975))
    pdffin <- rbind(pdfave, pdfci)
    return(pdffin)
  })
  
  dispList <- list()
  pdf <- apply(out$gam, 1, function(gam) {
    dgamma(xv, gam[1], gam[2])})
  pdfave <- apply(pdf, 1, mean)
  pdfci <- apply(pdf, 1, quantile, c(0.025, 0.975))
  dispList[[1]] <- rbind(pdfave, pdfci)
  
  plotCut99 <- list()
  for (i in 1:length(survList)) {
    plotCut99[[i]] <- which(survList[[i]][1, ] < 0.01)[1]
    if (is.na(plotCut99[[i]])) plotCut99[[i]] <- length(xv)
  }
  plotCut95 <- list()
  for (i in 1:length(survList)) {
    plotCut95[[i]] <- which(survList[[i]][1, ] < 0.05)[1]
    if (is.na(plotCut95[[i]])) plotCut95[[i]] <- length(xv)
  }
  
  names(survList) <- names(mortList) <- names(pdfList) <- 
    names(plotCut95) <- names(plotCut99) <- c("f", "m")
  return(list(mort = mortList, surv = survList, pdf = pdfList, x = xv, 
              cuts99 = plotCut99, cuts95 = plotCut95, disp = dispList))
}


