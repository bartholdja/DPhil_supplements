# Code to analytically calculate the elasticities of biometric heritability
# with respect to changes in the parameter of kernel component functions
# of an age-stage based population model.

# 19 Feb 2015
# Comments: 
# 1) We distinguish between age and age classes. Age is in discrete time
# steps of 1 year (1,2, ..., maximum age). 
# The first age, a = 0, contains the individuals aged 0 to 0.9999 years old.
# The other age classes contain the individuals of older ages analogously.
# Age classes can contain multiple ages, e.g. # "newborns" (age 0), 
# "yearlings" (age 1), "subadults" (ages 2 and 3), etc.
# Functions in separate R script: heritabiliyFctsXXX.R
# 2) to run the analysis for Soay sheep or roe deer set "modelsystem"
# to "soay" or "roe", respectively (line 32)
# 3) this script sources functions from "heritabilityFcts.R"

rm(list=ls())

if ("Viktualia" %in% list.files("/Users/")) {
  setwd("/Users/Viktualia/Dropbox/Projects/006_Heritability/")
} else {
  setwd("...")
}

vers <- "005"
# source functions file
source("Rscripts/heritabilityFcts.R")

# Part (I) ##############################################################
# Set IPM starting values by model system (soay, roe):
modelsystem <- "soay"  # set to "roe" or "soay"
if (modelsystem == "soay") {
  n.age <- 9  # Number of ages in bigmatrix (note four age-classes, just generating expanded matrix)
  nn <-   100  # Number of stage classes
  minsize <- 2.9  # smallest observed size, 2.9 kg
  maxsize <- 34.2  # largest observed size, 34.2 kgb
  
  # parameter list (values from supplementary material Coulson et al. 2010 JAE)
  # four age classes: lambs, yearlings, prime-aged adults, senescent individuals
  # lambs: 0 to 0.999
  # yearlings 1 to 1.999
  # prime-aged adults: 2 to 6.999
  # senescent individuals: 7 to 8.999
  # Parameter lists are always in the the same order:
  # intercepts (lambs, yearlings, adults, senescents)
  # then slopes (lambs, yearlings, adults, senescents)
  surv.params <- c(-3.5652572, -2.298034, 0.1900181, -5.6054482,
                   0.2801004, 0.2552713, 0.1030863, 0.2841259)
  fec.params <- c(-1000,-1000,-6.333452,-10.871048, 0,0,0,0.2860888)
  rep.params <- c(-4.2702348,-5.6249537,-0.473731707,-3.48380923,
                  0.1571317,0.2404116,0.004208452,0.09215088 )
  growth.params.mean <- c(7.6561033,8.2248039,6.1479333,4.4693660, 
                          0.7236385,0.7108617,0.7504923,0.7929918)
  growth.params.var <- c(2.39297,1.1556, -1.1010, 10.0938, 
                         0.04811, 0.1063, 0.1617, -0.2845)
  inheritance.params.mean <- c(1.0891092,6.4214297,9.0065115,13.53349697,
                               0.5015036,0.2941641,0.1651943,-0.04821571)
  inheritance.params.var <- c(5.36041, -1.1985, -1.3743, 2.28464, 
                              -0.06596,0.2729, 0.2863, 0.08917)
  age <- c(1, 2, 3, 3, 3, 3, 3, 3, 4)
  n.ageClass <- length(unique(age))
  sexratio=1
  # upper and lower integration limits
  L <- 0; U <- 1.1*maxsize
  # boundary points b and mesh points y
  b <- L+c(0:nn)*(U-L)/nn
  y <- 0.5*(b[1:nn]+b[2:(nn+1)])
}

if (modelsystem == "roe") {
  minsize <- 1 # minimum possible size of a roe deer
  maxsize <- 40 # maximum possible size of a roe deer
  nn <- 100 # number of stage classes within the model
  n.age <- 12 # number of ages -- senescents pooled into final age
  max.age <- 50 # maximum possible age for a roe deer
  
  # parameter list
  # three age classes: yearlings, adults, senescents
  # Parameter lists are always in the the same order:
  # intercepts (yearlings, adults (8-11), > 11), slopes (yearlings, adults (8-11), > 11)
  
  surv.params <-c(-1.84906+0.69922,-1.84906,-1.84906-1.65167,0.14415,0.14415,0.14415)
  fec.params <- c(-1000,-3.326,-3.326-0.617,0.147,0.147,0.147) # yearlings don't produce
  rep.params <- c(-1000,-1.798,-1.798-0.111,0.120,0.120,0.120)
  growth.params.mean <- c( 7.93611+2.58202, 7.93611, 7.93611-0.53762, 0.68245,0.68245,0.68245)
  growth.params.var <- c(2.2645+1.2598,2.2645,2.2645+0.6738,0,0,0)
  inheritance.params.mean<- c(0, 4.726-0.904, 4.726-0.904+ 0.979, 0.501, 0.501, 0.501)
  inheritance.params.var <- c(0,3.6339+0.3521,3.6339+0.3521,0,0,0)
 # Define age classes                
  age=c(1,2,2,2,2,2,2,3,3,3,3,3)
  n.ageClass <- length(unique(age))
  sexratio=2
  # upper and lower integration limits
  L<-0.9*minsize; U<-1.1*maxsize
  # boundary points b and mesh points y
  b<-L+c(0:nn)*(U-L)/nn
  y<-0.5*(b[1:nn]+b[2:(nn+1)])
}

# Part (II) ##############################################################
# Compute the elasticities

# Populate the IPM matrices for all ages
sub.mats <- NULL
for (i in 1:n.age){
  # params = (SI,SS,GMI,GMS,GVI,GVS,FI,FS,RI,RS,DMI,DMS,DVI,DVS)
  params <- c(surv.params[age[i]],surv.params[age[i]+max(age)],growth.params.mean[age[i]],growth.params.mean[age[i]+max(age)],
              growth.params.var[age[i]],growth.params.var[age[i]+max(age)],fec.params[age[i]],fec.params[age[i]			+max(age)],rep.params[age[i]],rep.params[age[i]+max(age)],inheritance.params.mean[age[i]],inheritance.params.mean[age[i]+max(age)],inheritance.params.var[age[i]],inheritance.params.var[age[i]+max(age)])
  sub.mats[[i]] <- bigmatrix(y,nn,params,sexratio)
}

# create Zh, e, and I
Zh=diag(y)  # diagonal matrix of mid-point trait values of each size class
e=rep(1,nn)  # all-ones vector of length equal to number of size classes (nn)
I=diag(nn)  # Identity matrix of dimensions nn

# survivorship matrices by age
La <- life.history(sub.mats,n.age)

# fertility and survival matrices by age
Fa<-Pa <- array(NA,c(nn,nn,n.age))
for ( i in c(1:n.age)){                
  Fa[,,i]=sub.mats[[i]]$D%*%sub.mats[[i]]$R 
  Pa[,,i]=sub.mats[[i]]$G%*%sub.mats[[i]]$S 
}

# find r by solving the renewal matrix A(r)
u=uniroot(findr,c(0,2),Fa,La,Pa,n.age, tol=0.0000001)
r=u$root

# compute Ar and Ao
z = exp(-r);

Ao=Fa[,,1]%*%La[,,1]
Ar = z*Fa[,,1]%*%La[,,1]
mala=sub.mats[[1]]$R %*%La[,,1]
for( ii in 2:(n.age-1)){
  Ao = Ao + Fa[,,ii]%*%La[,,ii]
  Ar = Ar + z^ii*Fa[,,ii]%*%La[,,ii]
  mala=mala+sub.mats[[ii]]$R %*%La[,,ii]}
Ar = Ar + z^n.age*Fa[,,n.age]%*%solve(I - z*Pa[,,n.age])%*%La[,,n.age]
Ao = Ao + Fa[,,n.age]%*%solve(I - Pa[,,n.age])%*%La[,,n.age]
mala=mala+sub.mats[[n.age]]$R %*%solve(I - Pa[,,n.age])%*%La[,,n.age]

# Find the eigenvector u (stable birth cohort size distribution)
goat = get.eigen.stuff(Ar)  # CHECK -lambda=1 goat[[1]] is 1
u=goat[[2]]
v=goat[[3]]
v = abs(v);
v= v/(v%*%u)[1]
 
# life time reproduction
K=e%*%Ao%*%u

# H, used for brevity
H=Ao/K[1,1]

# parents at own birth phenotype variance (Var_Xo)
Xop=(e%*%mala%*%Zh%*%u)/K[1,1]
SXoap2=e%*%mala%*%Zh%*%Zh%*%u
Var_Xo=SXoap2[1,1]/K[1,1]-(Xop[1,1])^2
  
# covariance between parent phenotype at own birth (Xo) and 
# offspring phenotype at birth (Y)
cov_YXosimp=e%*%Zh%*%H%*%(I-u%*%e%*%mala/K[1,1])%*%Zh%*%u

# compute biometric heritability
b=cov_YXosimp/Var_Xo
b*2

# Delta a
dAo_rec=dAo_Rec_sensitivity(fec.params, rep.params,age,sub.mats,Pa,Fa,La,y, n.age,nn,r, sexratio)
dAo_sur=dAo_Sur_sensitivity(surv.params,age,sub.mats,Pa,Fa,La,y, n.age,nn,r)
dAo_gro=dAo_Gro_sensitivity(age,sub.mats,Pa,Fa,La,y, n.age,nn,r)
dAo_inh=dAo_Inh_sensitivity(age,sub.mats,Pa,Fa,La,y, n.age,nn,r)

# delta u and Delta H
fer_int_u_H=deltau_DeltaH(dAo_rec[[1]],dAo_rec[[1+4]],Fa,La,u,v,r,Ao,Ar,K, H,e,n.age,nn,max.age)
fer_slp_u_H=deltau_DeltaH(dAo_rec[[3]],dAo_rec[[3+4]],Fa,La,u,v,r,Ao,Ar,K, H,e,n.age,nn,max.age)
twi_int_u_H=deltau_DeltaH(dAo_rec[[2]],dAo_rec[[2+4]],Fa,La,u,v,r,Ao,Ar,K, H,e,n.age,nn,max.age)
twi_slp_u_H=deltau_DeltaH(dAo_rec[[4]],dAo_rec[[4+4]],Fa,La,u,v,r,Ao,Ar,K, H,e,n.age,nn,max.age)
sur_int_u_H=deltau_DeltaH(dAo_sur[[1]],dAo_sur[[1+2]],Fa,La,u,v,r,Ao,Ar,K, H,e,n.age,nn,max.age)
sur_slp_u_H=deltau_DeltaH(dAo_sur[[2]],dAo_sur[[2+2]],Fa,La,u,v,r,Ao,Ar,K, H,e,n.age,nn,max.age)
grom_int_u_H=deltau_DeltaH(dAo_gro[[1]],dAo_gro[[1+4]],Fa,La,u,v,r,Ao,Ar,K, H,e,n.age,nn,max.age)
grom_slp_u_H=deltau_DeltaH(dAo_gro[[3]],dAo_gro[[3+4]],Fa,La,u,v,r,Ao,Ar,K, H,e,n.age,nn,max.age)
grov_int_u_H=deltau_DeltaH(dAo_gro[[2]],dAo_gro[[2+4]],Fa,La,u,v,r,Ao,Ar,K, H,e,n.age,nn,max.age)
grov_slp_u_H=deltau_DeltaH(dAo_gro[[4]],dAo_gro[[4+4]],Fa,La,u,v,r,Ao,Ar,K, H,e,n.age,nn,max.age)
inhm_int_u_H=deltau_DeltaH(dAo_inh[[1]],dAo_inh[[1+4]],Fa,La,u,v,r,Ao,Ar,K, H,e,n.age,nn,max.age)
inhm_slp_u_H=deltau_DeltaH(dAo_inh[[3]],dAo_inh[[3+4]],Fa,La,u,v,r,Ao,Ar,K, H,e,n.age,nn,max.age)
inhv_int_u_H=deltau_DeltaH(dAo_inh[[2]],dAo_inh[[2+4]],Fa,La,u,v,r,Ao,Ar,K, H,e,n.age,nn,max.age)
inhv_slp_u_H=deltau_DeltaH(dAo_inh[[4]],dAo_inh[[4+4]],Fa,La,u,v,r,Ao,Ar,K, H,e,n.age,nn,max.age)

# Elasticities of heritability to changes in the intercepts (int)
# and slopes (slp) of kernel functions (fer: probability of reproduction,
# twi: twinning rate, sur: survival, grom: growth mean, grov: growth variance,
# inhm: inheritance mean, inhv: inheritance variance)
Eb_fer_int=Elas_heri(dAo_rec[[1+8]],fer_int_u_H,e,Zh,u,H,nn,cov_YXosimp,Var_Xo,mala,K,Xop)
Eb_fer_slp=Elas_heri(dAo_rec[[3+8]],fer_slp_u_H,e,Zh,u,H,nn,cov_YXosimp,Var_Xo,mala,K,Xop)
Eb_twi_int=Elas_heri(dAo_rec[[2+8]],twi_int_u_H,e,Zh,u,H,nn,cov_YXosimp,Var_Xo,mala,K,Xop)
Eb_twi_slp=Elas_heri(dAo_rec[[4+8]],twi_slp_u_H,e,Zh,u,H,nn,cov_YXosimp,Var_Xo,mala,K,Xop)
Eb_sur_int=Elas_heri(dAo_sur[[1+4]],sur_int_u_H,e,Zh,u,H,nn,cov_YXosimp,Var_Xo,mala,K,Xop)
Eb_sur_slp=Elas_heri(dAo_sur[[2+4]],sur_slp_u_H,e,Zh,u,H,nn,cov_YXosimp,Var_Xo,mala,K,Xop)
Eb_grom_int=Elas_heri(dAo_gro[[1+8]],grom_int_u_H,e,Zh,u,H,nn,cov_YXosimp,Var_Xo,mala,K,Xop)
Eb_grom_slp=Elas_heri(dAo_gro[[3+8]],grom_slp_u_H,e,Zh,u,H,nn,cov_YXosimp,Var_Xo,mala,K,Xop)
Eb_grov_int=Elas_heri(dAo_gro[[2+8]],grov_int_u_H,e,Zh,u,H,nn,cov_YXosimp,Var_Xo,mala,K,Xop)
Eb_grov_slp=Elas_heri(dAo_gro[[4+8]],grov_slp_u_H,e,Zh,u,H,nn,cov_YXosimp,Var_Xo,mala,K,Xop)
Eb_inhm_int=Elas_heri(dAo_inh[[1+8]],inhm_int_u_H,e,Zh,u,H,nn,cov_YXosimp,Var_Xo,mala,K,Xop)
Eb_inhm_slp=Elas_heri(dAo_inh[[3+8]],inhm_slp_u_H,e,Zh,u,H,nn,cov_YXosimp,Var_Xo,mala,K,Xop)
Eb_inhv_int=Elas_heri(dAo_inh[[2+8]],inhv_int_u_H,e,Zh,u,H,nn,cov_YXosimp,Var_Xo,mala,K,Xop)
Eb_inhv_slp=Elas_heri(dAo_inh[[4+8]],inhv_slp_u_H,e,Zh,u,H,nn,cov_YXosimp,Var_Xo,mala,K,Xop)

# Eb_fer_intagec=c(Eb_fer_int)

Eb0=c(Eb_fer_int,Eb_fer_slp,Eb_twi_int,Eb_twi_slp,Eb_sur_int,Eb_sur_slp,Eb_grom_int,Eb_grom_slp,Eb_grov_int,Eb_grov_slp,Eb_inhm_int,Eb_inhm_slp,Eb_inhv_int,Eb_inhv_slp)  

# sum the elasticities across ages within age classes
CalcElasAgeClass <- function(age, elasAge) {
  temp <- NULL
  for (i in 1:max(age)) {
    temp[i] <- sum(elasAge[age == i])
  }
  return(temp)
}

eb_AC_fer_int <- CalcElasAgeClass(age, Eb_fer_int)
eb_AC_fer_slp <- CalcElasAgeClass(age, Eb_fer_slp)
eb_AC_twi_int <- CalcElasAgeClass(age, Eb_twi_int)
eb_AC_twi_slp <- CalcElasAgeClass(age, Eb_twi_slp)
eb_AC_sur_int <- CalcElasAgeClass(age, Eb_sur_int)
eb_AC_sur_slp <- CalcElasAgeClass(age, Eb_sur_slp)
eb_AC_grom_int <- CalcElasAgeClass(age, Eb_grom_int)
eb_AC_grom_slp <- CalcElasAgeClass(age, Eb_grom_slp)
eb_AC_grov_int <- CalcElasAgeClass(age, Eb_grov_int)
eb_AC_grov_slp <- CalcElasAgeClass(age, Eb_grov_slp)
eb_AC_inhm_int <- CalcElasAgeClass(age, Eb_inhm_int)
eb_AC_inhm_slp <- CalcElasAgeClass(age, Eb_inhm_slp)
eb_AC_inhv_int <- CalcElasAgeClass(age, Eb_inhv_int)
eb_AC_inhv_slp <- CalcElasAgeClass(age, Eb_inhv_slp)  
