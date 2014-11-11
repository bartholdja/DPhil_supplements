# Author: Julia A. Barthold
# Contact email: bartholdja@gmail.com
# Date: 7 Jan 2014
# Description: Code to run two-sex IPM
# Comment: The code is based on a two-sex IPM written by S. Schindler in MatLab. 
#   The ipm functions are based on code provided by the authors of 
#   JAE 2010, 79, 1226-1240. Thanks to them and to the people whose code contributions
#   they acknowledged. Throught the script "Fem" or "fem" denotes objects pertaining
#   to females and "Mal" or "mal" those pertaining to males.

rm(list = ls())
setwd("/Users/Viktualia/Dropbox/Projects/003_Lions")
load("Paper/JAE/forDataRepository/parameters.Rdata")

##########
# The functions for the IPM

# Survival function S(z,t)
sFunFem2sex = function(z, intercept, zSlope) {
  u = exp(intercept + zSlope * z) # predict untransformed values from binomial
  return(u / (1+u)) # transform
}

sFunMal2sex = function(z, intercept, zSlope) {
  u = exp(intercept + zSlope * z)
  return(u / (1 + u))
}

# Development function G(z'|z)  zz = z' for R code
gFunFem2sex = function(z, zz, interceptMu, zSlopeMu, interceptVar, zSlopeVar) {
  muZ = interceptMu + zSlopeMu * z
  sigmaZsq = interceptVar
  sigmaZ = sqrt(sigmaZsq)
  temp1 = sqrt(2 * pi) * sigmaZ
  temp2 = ((zz-muZ) ^ 2) / (2 * sigmaZsq)
  return(exp(-temp2) / temp1)
}

gFunMal2sex = function(z, zz, interceptMu, zSlopeMu, interceptVar) {
  muZ = interceptMu + zSlopeMu * z
  sigmaZsq = rep(as.vector(interceptVar), length(z))
  sigmaZ = sqrt(sigmaZsq)
  temp1 = sqrt(2 * pi) * sigmaZ
  temp2 = ((zz-muZ) ^ 2) / (2 * sigmaZsq)
  return(exp(-temp2) / temp1)
}

# Recruitment function R(z,t): in this case number of offspring
rFunFem2sex = function(z, intercept, zSlope) { 
  v = exp(intercept + zSlope * z)
  return(v)
}

# Reproduced if mated for females
repIfMatedFem2sex = function(z, interceptFem, zSlopeFem) {
  u = exp(interceptFem + zSlopeFem * z)
  return(u / (1 + u))
}

# Mating probability
maleSizeAdvantage = function(x, rho) {0.5 * exp(x/50 * rho)}

# Inheritance function D(z'|z,x)  zz = z' for R code.  Note same structure as G
# in this case -- could have different probability distribution if needed
# 
dFunFem2sex = function(z, zz, interceptMu, zSlopeMu, interceptVa) {
  muZ = interceptMu + zSlopeMu * z
  sigmaZsq = as.vector(interceptVa)
  sigmaZ = sqrt(sigmaZsq)
  temp1 = sqrt(2 * pi) * sigmaZ
  temp2 = ((zz-muZ) ^ 2) / (2 * sigmaZsq)
  return(exp(-temp2) / temp1)
}

# D(x'|z,x)
dFunMal2sex = function(z, zz, interceptMu, zSlopeMu, interceptVa) {
  muZ = interceptMu + zSlopeMu * z
  sigmaZsq = as.vector(interceptVa)
  sigmaZ = sqrt(sigmaZsq)
  temp1 = sqrt(2 * pi) * sigmaZ
  temp2 = ((zz-muZ) ^ 2) / (2 * sigmaZsq)
  return(exp(-temp2) / temp1)
}

# The 'big matrix' M of size n x n
bigmatrixFem2sex = function(n, sParams, rParams, gParams, dParams, maxsize, minsize) {
  
  # boundary points b and mesh points y
  b = minsize + c(0 : n) * (maxsize - minsize) / n  # minsize + "200 bins" * "width of interval"
  y = b[1:n] + (maxsize - minsize) / (2 * n)  # left boundary points + "half a width of interval"
  
  # create S, R, G and D matrices
  S <- diag(sFunFem2sex(y, sParams[1], sParams[2]))
  RoffNumb = matrix(rFunFem2sex(y, rParams[3], rParams[4]), n, n) # number of offspring
  RprobRepIfMated = matrix(repIfMatedFem2sex(y, rParams[1], rParams[2]), n, n) # probability of reproduced if mated rows = female character size, columns = male charactersize
  G <- t(outer(y ,y, gFunFem2sex, gParams[1], gParams[2], gParams[3])) # calculates the p of G(z'|z,t) for every z and z' combination
  
  # scale G so columns sum to 1
  G <- G / matrix(as.vector(apply(G, 2, sum)), nrow = n, ncol = n, byrow = TRUE) # because the probability that if you are z size at t that you grow to the same or another size at t+1 has to be 1
  G[apply(G, 2, is.na)] = 0
  
  # kernel to store offspring sizes depending on the mean couple size
  # i = x (female size), j = y (male size), z = x' (female offspring size)
  D <- array(NA, dim = c(150, 150, 150))
  # function to return mean parent size
  meanParentSize <- function(a,b) {
    # function to return mean parent size
    c <- (a+b)/2; return(c)}
  
  # fill the kernel (temp == matrix with mean couple size, dFunMal2sex performs
  # pointwise on each matrix element and returns the probability that a mean couple size
  # "z" produces an offspring of size "zz")
  for( l in 1 : dim(D)[3]) {
    # i is row
    # j is column
    # l is layer
    temp <- outer(y,y, meanParentSize)
    D[ , , l] <- dFunFem2sex(z = temp, zz = y[l],
                             interceptMu = dParams[1], zSlopeMu = dParams[2],
                             interceptVa = dParams[3])
  }
  rm(l, temp)
  
  # Normalise so that offspring size distribution out of one
  # mean couple size sums to one (sum up in z direction and divide each i, j element 
  # in z direction by the sum of all i,j elements in z direction)
  forNorm <- apply(D, FUN = sum, MARGIN = c(1,2))
  
  # loop to divide each matrix in z direction by the matrix that stores the sums
  for (k in 1: dim(D)[3]) {
    D[ , , k] <- D[ , , k]/forNorm
  }
  
  return(list(S = S, RoffNumb = RoffNumb, RprobRepIfMated = RprobRepIfMated, G = G, D = D, meshpts = y))
}

bigmatrixMal2sex<-function(n, sParams, rParams, gParams, dParams, maxsize, minsize) {
  
  # boundary points b and mesh points y
  b <- minsize + c(0 : n) * (maxsize - minsize) / n
  y <- b[1:n] + (maxsize - minsize) / (2 * n)
  
  # create S, R, G and D matrices
  S <- diag(sFunMal2sex(y, sParams[1], sParams[2]))
  G <- t(outer(y ,y, gFunMal2sex, gParams[1], gParams[2], gParams[3]))
  
  # scale G so columns sum to 1
  G <- G / matrix(as.vector(apply(G, 2, sum)), nrow = n, ncol = n, byrow = TRUE)
  G[apply(G, 2, is.na)] = 0
  
  # kernel to store offspring sizes depending on the mean couple size
  # i = x (female size), j = y (male size), z = y' (male offspring size)
  D <- array(NA, dim = c(150, 150, 150))
  
  meanParentSize <- function(a,b) {
    # function to return mean parent size
    c <- (a+b)/2; return(c)}
  
  # fill the kernel (temp == matrix with mean couple size, dFunMal2sex performs
  # pointwise on each matrix element and returns the probability that a mean couple size
  # "z" produces an offspring of size "zz")
  for( l in 1 : dim(D)[3]) {
    # i is row
    # j is column
    # l is layer
    temp <- outer(y,y, meanParentSize)
    D[ , , l] <- dFunMal2sex(z = temp, zz = y[l],
                             interceptMu = dParams[1], zSlopeMu = dParams[2],
                             interceptVa = dParams[3])
  }
  rm(l, temp)
  
  # Normalise so that offspring size distribution out of one
  # mean couple size sums to one (sum up in z direction and divide each i, j element 
  # in z direction by the sum of all i,j elements in z direction)
  forNorm <- apply(D, FUN = sum, MARGIN = c(1,2))  # sums all the i,j elements in z direction
  
  # loop to divide each matrix in z direction by the matrix that stores the sums
  for (k in 1: dim(D)[3]) {
    D[ , , k] <- D[ , , k]/forNorm
  }
  
  return(list(S = S, G = G, D = D, meshpts = y))
}

##########
# Defining initial parameters
minsize = 0  
maxsize = 160
nn = 150    # nn is the size of the matrix
sexRatioAtRecruitment = 0.45  # number of females to males who ever made it to their 1st birthday in Serengeti study 1368/(1368+1680)
epsilon = 0.000001 # difference in population size when to break off iterations
delta = 0.001 # difference in population size of survivors in R0 calculation
numberIterations = 150 # number of iterations when to break off if epsilon is not reached
sizeThresholdRepFem = 87  # size of smallest 2.5 year old female
sizeThresholdRepMal = 100  # size of smallest 2.5 year old male
rho = 0.25  # determines the slope of male-size advantage in the mating function
meshpts = # meshpts of the matrix
  (minsize + c(0:nn) * (maxsize - minsize) / nn)[1:nn] + (maxsize - minsize) / (2 * nn)

parameters = c(sParamsFem, rParamsFem, gParamsFem, dParamsFem[1:3], sexRatioAtRecruitment, sParamsMal, gParamsMal, 
               dParamsMal[4:6], rho)

# create objects to store model outputs
lambda =    vector("list", length(parameters)+1)
popDisFem = vector("list", length(parameters)+1)
popDisMal = vector("list", length(parameters)+1)
popFem = vector("list", length(parameters)+1)
mxFem = vector("list", length(parameters)+1)
mxDaughtFem = vector("list", length(parameters)+1)
mxSonsFem = vector("list", length(parameters)+1)
lxFem = vector("list", length(parameters)+1)
R0Fem = NULL
R0DaughtFem = NULL
R0SonsFem = NULL
generationTimeFem = NULL


names(parameters) = c("survInterceptFem", "survSlopeFem", "probRepInterceptFem", "probRepSlopeFem",
                      "numbOffInterceptFem", "numbOffSlopeFem", "growthInterceptFem", "growthSlopeFem", 
                      "growthVarInterceptFem","inherInterceptFem", "inherSlopeFem" , "inherVarFem", "sexRatioBirth", "survInterceptMal", 
                      "survSlopeMal", "growthInterceptMal", "growthSlopeMal",
                      "growthVarIntercept", "inherSlopeMal", "inherInterceptMal", "inherVarMal", "rho")

##########
# Run the two-sex IPM and calculate output objects 
for(k in 1:(length(parameters)+1)) {
  # for the whole IPM female size (x) increases with increasing i,
  # male size (y) increases with increasing j, and offspring sizes (x' or y')
  # increase with increasing z in the inheritance kernel D
  
  lambdaTemp = NULL
  params = parameters
  
  # upwards perturbations, one parameter at a time, after the first iteration
  # first iteration = unperturbed model
  if (k > 1) {
    if(params[k-1] > 0) params[k-1] = params[k-1] * 1.01
    if(params[k-1] < 0) params[k-1] = params[k-1] * 0.99
  }
  
  sParamsFem = params[1:2]
  rParamsFem = params[3:6]
  gParamsFem = params[7:9]
  dParamsFem = params[10:12]
  sexRatioAtRecruitment = params[13]
  
  sParamsMal = params[14:15]
  gParamsMal = params[16:18]
  dParamsMal = params[19:21]
  rho = params[22] # varies the strength of the male size advantage
  
  # Female kernels
  MFem2sex = bigmatrixFem2sex(n=nn,sParams=sParamsFem,rParamsFem,gParamsFem, dParamsFem, 
                              maxsize, minsize)
  
  # Male kernels
  MMal2sex = bigmatrixMal2sex(nn, sParamsMal, rParamsMal, gParamsMal, dParamsMal,
                              maxsize, minsize)
  
  SandGMal = MMal2sex$G %*% MMal2sex$S
  
  # substitute survival of small females with survival of small males
  # the substitute threshold is chosen dynamically as to guarantee a 
  # smooth survival function for females
  if(k == 1) {
    # determine the unperturbed survival substitute only in the first loop
    # find the threshold
    substitute.threshold <- max(which(diag(MMal2sex$S) >= diag(MFem2sex$S)))
    # substitute 
    femSurvivalSubstitute <- diag(MMal2sex$S)[1:substitute.threshold]
  }
  
  # replace survival values of small females with survival values of small males
  diag(MFem2sex$S)[1:substitute.threshold] <- femSurvivalSubstitute
  
  # Female survival and growth kernel
  SandGFem = MFem2sex$G %*% MFem2sex$S
  
  # empty matrix for mating probabilities
  M = matrix(0, length(MFem2sex$meshpts), length(MFem2sex$meshpts))
  # find meshpts for minimum sizes for reproduction
  
  indexThresholdFem = min(which(MFem2sex$meshpts >= sizeThresholdRepFem))
  indexThresholdMal = min(which(MMal2sex$meshpts >= sizeThresholdRepMal))
  
  # mating function (mating system = male size advantage)
  for (j in indexThresholdMal:ncol(M)) {
    # fills the column M [indexThreshold, ] with probability of mating as a 
    # function of male character sizes
    M[indexThresholdFem, j] = maleSizeAdvantage(
      MMal2sex$meshpts[j], rho = rho)
  } 
  
  # copy the row M[indexThresholdFem, ] for all rows bigger than indexThresholdFem
  for(i in indexThresholdFem:nrow(M)) {M[i, ] = M[indexThresholdFem, ]}
  
  # normalisation so that one column sums up to 1 (all females have a mating probability of 1)
  M[indexThresholdFem:nrow(M), ] = M[indexThresholdFem:nrow(M), ]/ rowSums(M[indexThresholdFem:nrow(M), ])
  
  # Number of offspring depends only on female sizes
  offspringNumb = MFem2sex$RoffNumb
  
  # Initial population
  nf = matrix(1, nrow = numberIterations, ncol = nn) # all sizes = 1
  nm = matrix(1, nrow = numberIterations, ncol = nn) # all sizes = 1
  
  i = 1
  # get densities
  wholePop = sum(nf[i, ]) + sum(nm[i, ])
  nf[i, ] = nf[i, ]/wholePop
  nm[i, ] = nm[i, ]/wholePop
  
  
  # first offspring generation
  i = 2 # "second" iteration
  
  # number of matings
  numbMatings = M * matrix(nf[i-1, ], length(nf[i-1, ]), 1) %*% nm[i-1, ] ## Nominater normalisation constant (without integration)
  
  # normalization such that sum of all couples amounts to number of
  # female density that are in their reproductive age
  if (sum(numbMatings) > 0) {
    numbMatingsNorm = numbMatings*sum(nf[i-1,which(MFem2sex$meshpts>= sizeThresholdRepFem)])/sum(numbMatings)}
  
  # number of offspring per realised mating combination = number of matings *
  # probability of female reproduction if mating occurred * number of offspring if reproduced
  numbOffPerCouple = numbMatingsNorm * MFem2sex$RprobRepIfMated * MFem2sex$RoffNumb
  
  # obtain female offspring size distribution at recruitment
  offspringFem = rep(NA, length(MFem2sex$meshpts))
  
  # multiply the matrix that stores the probablity that a couple of size x and y
  # produces an offspring of size x' (MFem2sex$D[ , , l]) with the number of offspring
  # produced by a couple of size x and y. Sum the probabilities to obtain number of
  # offspring that recruit to the population with size x'.
  for(l in 1:dim(MFem2sex$D)[3]){
    offspringFem[l] = sum(MFem2sex$D[ , , l] * numbOffPerCouple)
  }
  
  # same for male offspring
  offspringMal = rep(NA, length(MMal2sex$meshpts))
  
  for(l in 1:dim(MMal2sex$D)[3]){
    offspringMal[l] = sum(MMal2sex$D[ , , l] * numbOffPerCouple)
  }
  
  # put all together
  nf[i, ] = SandGFem %*% as.matrix(nf[i-1, ], nn, 1) + sexRatioAtRecruitment * offspringFem
  nm[i, ] = SandGMal %*% as.matrix(nm[i-1, ], nn, 1) + (1-sexRatioAtRecruitment) * offspringMal
  
  wholePop = sum(nf[i, ] + nm[i, ])
  
  lambdaTemp[i-1] = wholePop # because N(t-1) = 1 therefore N(t)/N(t-1) = N(t)
  
  # transform into densitites again
  nf[i, ] = nf[i, ]/wholePop
  nm[i, ] = nm[i, ]/wholePop
  
  # iterate IPM further
  for (i in 3:numberIterations) {
    print(paste(i,k))
    
    if(sum(abs(nf[i-1, ] - nf[i-2, ])) > epsilon | sum(abs(nm[i-1, ] - nm [i-2, ])) > epsilon) {
      
      # number of matings
      numbMatings = M * matrix(nf[i-1, ], length(nf[i-1, ]), 1) %*% nm[i-1, ] ## Nominater normalisation constant (without integration)
      
      # normalization such that sum of all couples amounts to number of
      # female density that are in their reproductive age
      if (sum(numbMatings) > 0) {
        numbMatingsNorm = numbMatings*sum(nf[i-1,which(MFem2sex$meshpts>= sizeThresholdRepFem)])/sum(numbMatings)}
      
      # number of offspring per realised mating combination = number of matings *
      # probability of female reproduction if mating occurre * number of offspring if reproduced
      numbOffPerCouple = numbMatingsNorm * MFem2sex$RprobRepIfMated * MFem2sex$RoffNumb
      
      #obtain female offspring size distribution at recruitment
      offspringFem = rep(NA, length(MFem2sex$meshpts))
      
      # multiply the matrix that stores the probablity that a couple of size x and y
      # produces an offspring of size x' (MFem2sex$D[ , , l]) with the number of offspring
      # produced by a couple of size x and y. Sum the probabilities to obtain number of
      # offspring that recruit to the population with size x'.
      for(l in 1:dim(MFem2sex$D)[3]){
        offspringFem[l] = sum(MFem2sex$D[ , , l] * numbOffPerCouple)
      }
      
      # same for male offspring
      offspringMal = rep(NA, length(MMal2sex$meshpts))
      
      for(l in 1:dim(MMal2sex$D)[3]){
        offspringMal[l] = sum(MMal2sex$D[ , , l] * numbOffPerCouple)
      }    
      
      # put all together
      nf[i, ] = SandGFem %*% as.matrix(nf[i-1, ], nn, 1) + sexRatioAtRecruitment * offspringFem
      nm[i, ] = SandGMal %*% as.matrix(nm[i-1, ], nn, 1) + (1-sexRatioAtRecruitment) * offspringMal
      
      wholePop = sum(nf[i, ] + nm[i, ])
      
      lambdaTemp[i-1] = wholePop # because N(t-1) = 1 therefore N(t)/N(t-1) = N(t)
      
      # scale back to density
      nf[i, ] = nf[i, ]/wholePop
      nm[i, ] = nm[i, ]/wholePop
      
    } else {break()}
  }
  # store the outputs
  lambda[[k]] = lambdaTemp
  popDisFem[[k]] = nf[1:(i-1), ]  # it's i-1 because the criterium for breaking is only
  # evaluated after i jumped to the next higher number
  popDisMal[[k]] = nm[1:(i-1), ]
  
  # calculate LRS
  # female cohort, scaled to 1
  startCohFem = sexRatioAtRecruitment * offspringFem/ sum(sexRatioAtRecruitment * offspringFem) # normal vector
  
  # size-specific number of offspring per time step 
  repFem = MFem2sex$RprobRepIfMated * MFem2sex$RoffNumb 
  
  # set reproduction for sizes smaller than threshold to 0
  repFem[which(MFem2sex$meshpts < sizeThresholdRepFem), ] = 0
  
  # male cohort at stable size distribution
  malCoh <- popDisMal[[k]][nrow(popDisMal[[k]]), ]
  
  # abort either if no survivors or when 150 years passed
  
  # let female cohort of size 1 die, grow, and reproduce and count offspring until less than 0.001 still alive
  
  for (j in 1: numberIterations){
    print(j)
    if(j == 1) {popTempFem = matrix(NA, nrow = numberIterations, ncol = nn)
                popTempFem[1, ] = startCohFem
                mxTempDaughtersFem <- NULL
                mxTempSonsFem = NULL
    }
    if(j != 1) { # j = 2 corresponds to t+1
      popTempFem[j, ] = matrix(popTempFem[j-1, ], 1, length(popTempFem[j-1, ])) %*% t(SandGFem)
      
      # number of matings
      numbMatings = (t(t(popTempFem[j-1, ])) %*% malCoh) * M 
      
      # normalization such that sum of all couples amounts to number of
      # female density that are in their reproductive age
      if (sum(numbMatings) > 0) {
        numbMatingsNorm = numbMatings*sum(popTempFem[j-1, which(MFem2sex$meshpts>= sizeThresholdRepFem)])/sum(numbMatings)}
      
      # number of offspring per realised mating combination = number of matings *
      # probability of female reproduction if mating occurre * number of offspring if reproduced
      numbOffPerCouple = numbMatingsNorm * MFem2sex$RprobRepIfMated * MFem2sex$RoffNumb
      
      
      mxTempDaughtersFem[j-1] <- sum(numbOffPerCouple) * sexRatioAtRecruitment
      mxTempSonsFem[j-1] <- sum(numbOffPerCouple) * (1 - sexRatioAtRecruitment)
    }
    if(sum(popTempFem[j, ]) <= delta) break() 
  }
  
  popFem[[k]] = popTempFem[1:j, ]
  
  mxDaughtFem[[k]] = mxTempDaughtersFem[1:(j-1)]
  
  mxSonsFem[[k]] = mxTempSonsFem[1:(j-1)]
  
  lxFem[[k]] = rowSums(popFem[[k]])
  
  R0DaughtFem[k] = sum(mxDaughtFem[[k]])
  
  R0SonsFem[k] = sum(mxSonsFem[[k]])
  
  generationTimeFem[k] = log(R0DaughtFem[k])/log(lambda[[k]] [length(lambda[[k]])])
  
}

# mean character sizes
meanFem = sum(popDisFem[[1]][nrow(popDisFem[[1]]), ] * MFem2sex$meshpts)/sum(popDisFem[[1]][nrow(popDisFem[[1]]), ])
meanMal = sum(popDisMal[[1]][nrow(popDisMal[[1]]), ] * MMal2sex$meshpts)/sum(popDisMal[[1]][nrow(popDisMal[[1]]), ])

# data frame storing the unperturbed and perturbed lambdas
returnStableLambda = function(x) {u = x[length(x)]; return(u)}
lambdas.temp = lapply(lambda, returnStableLambda)
lambdas = data.frame(do.call(rbind, lambdas.temp))
names(lambdas) = "Lambda.pert.anal"

# population densities for sex ratios
popDensFem = NULL
for (i in 1:(length(parameters)+1)){
  popDensFem[i] = rowSums(popDisFem[[i]])[nrow(popDisFem[[i]])]
}

popDensMal = NULL
for (i in 1:(length(parameters)+1)){
  popDensMal[i] = rowSums(popDisMal[[i]])[nrow(popDisMal[[i]])]
}

# matrix storing the unperturbed and perturbed overall sex ratio's
sexRatioDens = t(as.matrix(cbind(popDensFem, popDensMal)))