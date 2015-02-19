# This script gets sourced by heritability.R

# Part I: Kernel component functions and 
# Part II: Functions to compute the elasticities of heritabilities to changes
# in the parameter of kernel component functions.

# Part I:
# Survival function
# Intercept (a) and slope (b) to determine probability of survival 
# from t to t+1 at mid-point trait values (x).
sx<-function(x,a,b) {
  u<-exp(a+b*x)
  return((u/(1+u)))
}

# Growth kernel
# Intercepts and slopes for the mean (a,b) and for the variance (c,d)
# components of the kernel.
# The kernel contains the probabilities to grow from mid-point trait value (x) at t 
# to mid-point trait value (y) at t+1.
# Probabilities are normally distributed around
# the mean (mux) with variance (sigmax2).
gxy<-function(x,y,a,b,c,d) {
  sigmax2 <- c + d*x
  sigmax <- sqrt(sigmax2)
  mux <- a+b*x
  fac1 <- sqrt(2*pi)*sigmax
  fac2 <- ((y-mux)^2)/(2*sigmax2)
  return(exp(-fac2)/fac1)
}

# Reproduction function
# Intercept and slope (a, b) to determine the probability that a reproducing
# female gives birth to a singleton or to twins.
# Intercept and slope (c, d) to determine the
# probability that an female reproduces at mid-point trait value (x).
# Both functions are multiplied to obtain the expected number of offspring
# produced by an individual of trait value x.
rx<-function(x,a,b,c,d){ 
  u1 <- exp(a+b*x) # u1 twinning rate/fecundity-> fec.params
  u1 <- u1/(1+u1)
  u2 <- exp(c+d*x) # u2 probability to reproduce -> rep.params
  u2 <- u2/(1+u2)
  nkids <- u2*(1+u1)
  return(nkids)
}

# inheritance kernel
# intercepts and slopes for the mean (a,b) and for the variance (c,d)
# components of the kernel
# The kernel contains the probability that a female at mid-point trait value (x) at t
# produces an offspring of mid-point trait value y at t+1. 
# Probabilities are normally distributed around
# the mean (mux) with variance (sigmax2).
dxy<-function(x,y,a,b,c,d) {
  mux <- a+b*x
  sigmax2 <- c+d*x
  sigmax <- sqrt(sigmax2)
  fac1<-sqrt(2*pi)*sigmax
  fac2<-((y-mux)^2)/(2*sigmax2)
  f<-exp(-fac2)/fac1
  return(f)
}

############## The 'big matrix' M of size n x n
# y : vector of phenotype mid-point trait values
# n : number of phenotype classes
# p : vector storing the parameters estimated from data
bigmatrix<-function(y,n,p,sexratio) {
  
  # create S, R, G and D matrices
  S <- diag(sx(y,p[1],p[2]))
  R <- diag(rx(y,p[7],p[8],p[9],p[10]))
  G <- t(outer(y,y,gxy,p[3],p[4],p[5],p[6]))
  D <- t(outer(y,y,dxy,p[11],p[12],p[13],p[14]))
  # scale D and G so columns sum to 1
  G <- G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  D <- D/matrix(as.vector(apply(D,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  G <- ifelse(is.nan(G),0,G)
  D <- ifelse(is.nan(D),0,D)
  if(sexratio==2){R=R/2 ### if only one sex in the IPM and sex ratio=0.5
  }
  return(list(S=S,R=R,G=G,D=D,meshpts=y))
}

# For large matrices it is faster to calculate the dominant eigenvalue 
# and eigenvectors associated with it via iteration rather than using eigen.
get.eigen.stuff <- function(mat){ 
  sz <- dim(mat)[1]
  t.now <- runif(sz)
  t.now <- t.now/sum(t.now)
  t.next <- mat%*%t.now
  t.next <- t.next/sum(t.next)
  i <- 0
  while (sum(abs(t.next-t.now))>0.0000001){
    i <- i+1
    # print(i)
    t.now <- t.next
    t.next <- mat%*%t.now
    lambda <- sum(t.next)/sum(t.now)
    t.next <- t.next/sum(t.next)
  }
  r.now <- runif(sz)
  r.now <- r.now/sum(r.now)
  r.next <- r.now%*%mat
  r.next <- r.next/sum(r.next)
  while (sum(abs(r.next-r.now))>0.0000001){
    r.now <- r.next
    r.next <- r.now%*%mat
    r.next <- r.next/sum(r.next)
  }
  return(list(lambda,t.next,r.next))
}

# Tracks a cohort from birth to death recording number of individuals in
# each phenotype class as a function of the phenotype class they were born into (lx tensor). 
life.history <- function(sub.mats,steps){
  sz <- dim(sub.mats[[1]]$S)[1]
  lx <- array(NA,c(sz,sz,(steps+1)))
  lx[,,1] <- diag(1,nrow=sz)
  for (i in 2:(steps+1)){
    lx.mat <- sub.mats[[i-1]]$G%*%sub.mats[[i-1]]$S
    lx[,,i] <- lx.mat%*%lx[,,i-1]
  }
  return(lx)}

# Part II:
# Function to solve the renewal matrix A(r) to find r.
# x : variable to solve for, here r
# Fa, La, Pa : Fertility, survivorship, and survival matrices for ages a, 
# see manuscript for further details.
# n.age : number of ages -- senescents pooled into final class
# y1 : renewal matrix A(r)
findr<- function(x,Fa,La,Pa,n.age) {
  z = exp(-x) 
  y1 = z*Fa[,,1]%*%La[,,1]
  I=diag(dim(Fa[,,1])[1])  
  # sum over all ages younger than last age class to get A(r)
  for (ii in 2:(n.age-1)){
    y1 = y1 + z^ii*Fa[,,ii]%*%La[,,ii]
  }
  # sum over all ages in last age class
  # for closing term see S. K. Kumar et al. (in prep., section A.6)
  y1 = y1 + z^n.age*Fa[,,n.age]%*%solve(I - z*Pa[,,n.age])%*%La[,,n.age]
  ei = get.eigen.stuff(y1)
  r = log(ei[[1]])
  return(r)
}

# Functions to compute delta a
# Equ. 29, 14 November
dAo_Rec_sensitivity<- function(fec.params, rep.params,age,sub.mats,Pa,Fa,La,y, n.age,nn,r,sr){
  dAo_fert_slp= dAo_twin_slp=dAo_fert_int= dAo_twin_int=D1_fert_slp= D1_twin_slp=D1_fert_int= D1_twin_int=dmala_fert_int=dmala_twin_int=dmala_fert_slp=dmala_twin_slp<- array(NA,c(nn,nn,n.age))
  
  for (i in c(1:n.age)){
    #look at each age specific contribution
    u1 <- exp(fec.params[age[i]]+fec.params[age[i]+max(age)]*y)
    u1 <- u1/(1+u1)
    u2 <- exp(rep.params[age[i]]+rep.params[age[i]+max(age)]*y)
    u2 <- u2/(1+u2)
    
    dR_fert_int=u2*(1-u2)*(1+u1) # u2*(1+u1); u2 prob to reproduce, u1 twinning rate
    dR_twin_int=u2*(u1)*(1-u1)
    dR_fert_slp=y*u2*(1-u2)*(1+u1)
    dR_twin_slp=u2*y*u1*(1-u1)
    
    if (sr==2){dR_fert_int=u2*(1-u2)*(1+u1)/2 # u2*(1+u1); u2 prob to reproduce, u1 twinning rate
               dR_twin_int=u2*(u1)*(1-u1)/2
               dR_fert_slp=y*u2*(1-u2)*(1+u1)/2
               dR_twin_slp=u2*y*u1*(1-u1)/2}
    
    dFa_fert_int=sub.mats[[i]]$D%*%diag(dR_fert_int)
    dFa_twin_int=sub.mats[[i]]$D%*%diag(dR_twin_int)
    dFa_fert_slp=sub.mats[[i]]$D%*%diag(dR_fert_slp)
    dFa_twin_slp=sub.mats[[i]]$D%*%diag(dR_twin_slp)
    
    if(i<n.age){
      dAo_fert_int[,,i]=dFa_fert_int%*%La[,,i]
      dAo_twin_int[,,i]=dFa_twin_int%*%La[,,i]
      dAo_fert_slp[,,i]=dFa_fert_slp%*%La[,,i]
      dAo_twin_slp[,,i]=dFa_twin_slp%*%La[,,i]

	dmala_fert_int[,,i]=diag(dR_fert_int)%*%La[,,i]
      dmala_twin_int[,,i]=diag(dR_twin_int)%*%La[,,i]
      dmala_fert_slp[,,i]=diag(dR_fert_slp)%*%La[,,i]
      dmala_twin_slp[,,i]=diag(dR_twin_slp)%*%La[,,i]

      
      D1_fert_int[,,i]=exp(-r*i)*dFa_fert_int%*%La[,,i]
      D1_twin_int[,,i]=exp(-r*i)*dFa_twin_int%*%La[,,i]
      D1_fert_slp[,,i]=exp(-r*i)*dFa_fert_slp%*%La[,,i]
      D1_twin_slp[,,i]=exp(-r*i)*dFa_twin_slp%*%La[,,i]
    }else{
      Hs=solve(I-Pa[,,n.age]) 
      dAo_fert_int[,,n.age]=dFa_fert_int%*%Hs%*%La[,,n.age]
      dAo_twin_int[,,i]=dFa_twin_int%*%Hs%*%La[,,n.age]
      dAo_fert_slp[,,i]=dFa_fert_slp%*%Hs%*%La[,,n.age]
      dAo_twin_slp[,,i]=dFa_twin_slp%*%Hs%*%La[,,n.age]
	
	dmala_fert_int[,,i]=diag(dR_fert_int)%*%Hs%*%La[,,n.age]
      dmala_twin_int[,,i]=diag(dR_twin_int)%*%Hs%*%La[,,n.age]
      dmala_fert_slp[,,i]=diag(dR_fert_slp)%*%Hs%*%La[,,n.age]
      dmala_twin_slp[,,i]=diag(dR_twin_slp)%*%Hs%*%La[,,n.age]
      
      D1_fert_int[,,n.age]=exp(-r*i)*dFa_fert_int%*%Hs%*%La[,,n.age]
      D1_twin_int[,,i]=exp(-r*i)*dFa_twin_int%*%Hs%*%La[,,n.age]
      D1_fert_slp[,,i]=exp(-r*i)*dFa_fert_slp%*%Hs%*%La[,,n.age]
      D1_twin_slp[,,i]=exp(-r*i)*dFa_twin_slp%*%Hs%*%La[,,n.age]
    }
  }
  return(list(dAo_fert_int,dAo_twin_int,dAo_fert_slp,dAo_twin_slp,D1_fert_int,D1_twin_int,D1_fert_slp,D1_twin_slp,dmala_fert_int,dmala_twin_int,dmala_fert_slp,dmala_twin_slp))
}

dAo_Sur_sensitivity<- function(surv.params,age,sub.mats,Pa,Fa,La,y, n.age, nn,r){
  dAo_sur_slp= dAo_sur_int=D1_sur_slp= D1_sur_int=dmala_sur_int=dmala_sur_slp<- array(NA,c(nn,nn,n.age))
  Hw=solve(I-exp(-r)*Pa[,,n.age])
  Hs=solve(I-Pa[,,n.age]) 
  for (i in c(1:n.age)){
    # look at each age specific contribution
    # survival function
    s <- exp(surv.params[age[i]]+surv.params[age[i]+max(age)]*y)
    s <- s/(1+s)
    
    # derivatives of survival function with respect to intercept (int) and slope (slp)
    dS_sur_int=s*(1-s)
    dS_sur_slp=y*s*(1-s)
    
    # calculate the derivatives of Pa (Ga remains unchanged)
    dPa_sur_int=sub.mats[[i]]$G%*%diag(dS_sur_int)
    dPa_sur_slp=sub.mats[[i]]$G%*%diag(dS_sur_slp)
    
    # estimate the successive variations in La
    # see S.K. Kumar et al (in prep, section A5.5)
    if(i<n.age){
      Psucc=diag(nn)
      dLa_sur=D11_sur=dmala_sur=0
      for (b in (i+1):(n.age-1)){
        
        del=sub.mats[[b]]$D%*%sub.mats[[b]]$R%*%Psucc
       
        
        D11_sur=D11_sur+exp(-r*b)*del
        dLa_sur=dLa_sur+del
	dmala_sur=dmala_sur+sub.mats[[b]]$R%*%Psucc

 	Psucc=sub.mats[[b]]$G%*%sub.mats[[b]]$S%*%Psucc
      }
      
      
      dAo_sur_int[,,i]=dLa_sur%*%dPa_sur_int%*%La[,,i]+Fa[,,n.age]%*%Hs%*%Psucc%*%dPa_sur_int%*%La[,,i]
      dAo_sur_slp[,,i]=dLa_sur%*%dPa_sur_slp%*%La[,,i]+Fa[,,n.age]%*%Hs%*%Psucc%*%dPa_sur_slp%*%La[,,i]

	  dmala_sur_int[,,i]=dmala_sur%*%dPa_sur_int%*%La[,,i]+sub.mats[[n.age]]$R%*%Hs%*%Psucc%*%dPa_sur_int%*%La[,,i]
      dmala_sur_slp[,,i]=dmala_sur%*%dPa_sur_slp%*%La[,,i]+sub.mats[[n.age]]$R%*%Hs%*%Psucc%*%dPa_sur_slp%*%La[,,i]	
      
      D1_sur_int[,,i]=D11_sur%*%dPa_sur_int%*%La[,,i]+exp(-r*n.age)*(Fa[,,n.age]%*%Hw%*%Psucc%*%dPa_sur_int%*%La[,,i])
      D1_sur_slp[,,i]=D11_sur%*%dPa_sur_slp%*%La[,,i]+exp(-r*n.age)*(Fa[,,n.age]%*%Hw%*%Psucc%*%dPa_sur_slp%*%La[,,i])
    }else{
      dAo_sur_int[,,i]=Fa[,,n.age]%*%Hs%*%dPa_sur_int%*%Hs%*%La[,,n.age]
      dAo_sur_slp[,,i]=Fa[,,n.age]%*%Hs%*%dPa_sur_slp%*%Hs%*%La[,,n.age]

      dmala_sur_int[,,i]=sub.mats[[n.age]]$R%*%Hs%*%dPa_sur_int%*%Hs%*%La[,,n.age]
      dmala_sur_slp[,,i]=sub.mats[[n.age]]$R%*%Hs%*%dPa_sur_slp%*%Hs%*%La[,,n.age]
      
      D1_sur_int[,,i]=exp(-r*(n.age+1))*Fa[,,n.age]%*%Hw%*%dPa_sur_int%*%Hw%*%La[,,n.age]
      D1_sur_slp[,,i]=exp(-r*(n.age+1))*Fa[,,n.age]%*%Hw%*%dPa_sur_slp%*%Hw%*%La[,,n.age]
      
    }
  }
  return(list(dAo_sur_int,dAo_sur_slp,D1_sur_int,D1_sur_slp,dmala_sur_int,dmala_sur_slp))
}


dAo_Inh_sensitivity<- function(age,sub.mats,Pa,Fa,La,y, n.age,nn,r){
  dAo_inhm_slp= dAo_inhv_slp=dAo_inhm_int= dAo_inhv_int=D1_inhm_slp= D1_inhv_slp=D1_inhm_int=D1_inhv_int<- array(NA,c(nn,nn,n.age))
  
  for (i in c(1:(n.age))){
    #look at each age specific contribution
    dD_inhm_int=dD_inhm_slp=dD_inhv_int=dD_inhv_slp=array(NA,c(nn,nn))
    for (j in c(1:nn)){
      mu=sum(y*(sub.mats[[i]]$D[,j]))  # appendix section A.4, mu1?, equ. 34, but check!!!
      sigma=sqrt(sum(((y-mu)^2)*(sub.mats[[i]]$D[,j])))  # appendix section A.4, sigma1?, equ. 34, but check!!!
      
      for (k in c(1:nn)){
        dD_inhm_int[k,j]=((y[k]-mu)/(sigma^2))*(sub.mats[[i]]$D[k,j])  # derivative with respect to the intercept of the mean of the inheritance function
        dD_inhm_slp[k,j]=y[j]*dD_inhm_int[k,j]
        dD_inhv_int[k,j]=(1/(2*sigma))*((sub.mats[[i]]$D[k,j])/sigma) * ( (((y[k]-mu)^2)/(sigma^2)) - 1 )
        dD_inhv_slp[k,j]=y[j]*dD_inhv_int[k,j]
      }
    }
    dD_inhm_int<- ifelse(is.nan(dD_inhm_int),0,dD_inhm_int)
    dD_inhm_slp<- ifelse(is.nan(dD_inhm_slp),0,dD_inhm_slp)
    dD_inhv_int<- ifelse(is.nan(dD_inhv_int),0,dD_inhv_int)
    dD_inhv_slp<- ifelse(is.nan(dD_inhv_slp),0,dD_inhv_slp)
    
    dFa_inhm_int=dD_inhm_int%*%sub.mats[[i]]$R
    dFa_inhv_int=dD_inhv_int%*%sub.mats[[i]]$R
    dFa_inhm_slp=dD_inhm_slp%*%sub.mats[[i]]$R
    dFa_inhv_slp=dD_inhv_slp%*%sub.mats[[i]]$R
    
    if(i<n.age){
      
      dAo_inhm_int[,,i]=dFa_inhm_int%*%La[,,i]
      dAo_inhv_int[,,i]=dFa_inhv_int%*%La[,,i]
      dAo_inhm_slp[,,i]=dFa_inhm_slp%*%La[,,i]
      dAo_inhv_slp[,,i]=dFa_inhv_slp%*%La[,,i]
      
      D1_inhm_int[,,i]=exp(-r*i)*dFa_inhm_int%*%La[,,i]
      D1_inhv_int[,,i]=exp(-r*i)*dFa_inhv_int%*%La[,,i]
      D1_inhm_slp[,,i]=exp(-r*i)*dFa_inhm_slp%*%La[,,i]
      D1_inhv_slp[,,i]=exp(-r*i)*dFa_inhv_slp%*%La[,,i]
      
    }else {
      Hs=solve(I-Pa[,,n.age]) 
      dAo_inhm_int[,,n.age]=dFa_inhm_int%*%Hs%*%La[,,n.age]
      dAo_inhv_int[,,n.age]=dFa_inhv_int%*%Hs%*%La[,,n.age]
      dAo_inhm_slp[,,n.age]=dFa_inhm_slp%*%Hs%*%La[,,n.age]
      dAo_inhv_slp[,,n.age]=dFa_inhv_slp%*%Hs%*%La[,,n.age]
      
      D1_inhm_int[,,n.age]=exp(-r*i)*dFa_inhm_int%*%Hs%*%La[,,n.age]
      D1_inhv_int[,,n.age]=exp(-r*i)*dFa_inhv_int%*%Hs%*%La[,,n.age]
      D1_inhm_slp[,,n.age]=exp(-r*i)*dFa_inhm_slp%*%Hs%*%La[,,n.age]
      D1_inhv_slp[,,n.age]=exp(-r*i)*dFa_inhv_slp%*%Hs%*%La[,,n.age]
      
    }
  }
  
  return(list(dAo_inhm_int,dAo_inhv_int,dAo_inhm_slp,dAo_inhv_slp,D1_inhm_int,D1_inhv_int,D1_inhm_slp,D1_inhv_slp,array(0,c(nn,nn,n.age)),array(0,c(nn,nn,n.age)),array(0,c(nn,nn,n.age)),array(0,c(nn,nn,n.age))))
}

dAo_Gro_sensitivity<- function(age,sub.mats,Pa,Fa,La,y, n.age,nn,r){
  dAo_grom_slp= dAo_grov_slp=dAo_grom_int= dAo_grov_int= D1_grom_slp=D1_grov_slp=D1_grom_int= D1_grov_int=dmala_grom_int=dmala_grov_int=dmala_grom_slp=dmala_grov_slp<- array(NA,c(nn,nn,n.age))
  
  for (i in c(1:n.age)){
    # look at each age specific contribution
    dG_grom_int=dG_grom_slp=dG_grov_int=dG_grov_slp=array(NA,c(nn,nn))
    for (j in c(1:nn)){
      mu=sum(y*(sub.mats[[i]]$G[,j]))
      sigma=sqrt(sum(((y-mu)^2)*(sub.mats[[i]]$G[,j])))
      
      for (k in c(1:nn)){
        dG_grom_int[k,j]=((y[k]-mu)/(sigma^2))*(sub.mats[[i]]$G[k,j])
        dG_grom_slp[k,j]=y[j]*dG_grom_int[k,j]
        dG_grov_int[k,j]=(1/(2*sigma))*((sub.mats[[i]]$G[k,j])/sigma) * ( (((y[k]-mu)^2)/(sigma^2)) - 1 )
        dG_grov_slp[k,j]=y[j]*dG_grov_int[k,j]
      }
    }
    dG_grom_int<- ifelse(is.nan(dG_grom_int),0,dG_grom_int)
    dG_grom_slp<- ifelse(is.nan(dG_grom_slp),0,dG_grom_slp)
    dG_grov_int<- ifelse(is.nan(dG_grov_int),0,dG_grov_int)
    dG_grov_slp<- ifelse(is.nan(dG_grov_slp),0,dG_grov_slp)
    
    dPa_grom_int=dG_grom_int%*%sub.mats[[i]]$S
    dPa_grov_int=dG_grov_int%*%sub.mats[[i]]$S
    dPa_grom_slp=dG_grom_slp%*%sub.mats[[i]]$S
    dPa_grov_slp=dG_grov_slp%*%sub.mats[[i]]$S
    
    # estimates the successive variations in La   
    Hw=solve(I-exp(-r)*Pa[,,n.age])
    Hs=solve(I-Pa[,,n.age]) 
    
    if(i<n.age){
      Psucc=diag(nn)
      dLa_gro=D11_gro=dmala_gro=0
      for (b in c((i+1):(n.age-1))){
        del=sub.mats[[b]]$D%*%sub.mats[[b]]$R%*%Psucc
        
        
        D11_gro=D11_gro+exp(-r*b)*del
        dLa_gro=dLa_gro+del
	dmala_gro=dmala_gro+sub.mats[[b]]$R%*%Psucc

	Psucc=sub.mats[[b]]$G%*%sub.mats[[b]]$S%*%Psucc
      }
      
      dAo_grom_int[,,i]=dLa_gro%*%dPa_grom_int%*%La[,,i]+Fa[,,n.age]%*%Hs%*%Psucc%*%dPa_grom_int%*%La[,,i]
      dAo_grov_int[,,i]=dLa_gro%*%dPa_grov_int%*%La[,,i]+Fa[,,n.age]%*%Hs%*%Psucc%*%dPa_grov_int%*%La[,,i]
      dAo_grom_slp[,,i]=dLa_gro%*%dPa_grom_slp%*%La[,,i]+Fa[,,n.age]%*%Hs%*%Psucc%*%dPa_grom_slp%*%La[,,i]
      dAo_grov_slp[,,i]=dLa_gro%*%dPa_grov_slp%*%La[,,i]+Fa[,,n.age]%*%Hs%*%Psucc%*%dPa_grov_slp%*%La[,,i]

      dmala_grom_int[,,i]=dmala_gro%*%dPa_grom_int%*%La[,,i]+sub.mats[[n.age]]$R%*%Hs%*%Psucc%*%dPa_grom_int%*%La[,,i]
      dmala_grov_int[,,i]=dmala_gro%*%dPa_grov_int%*%La[,,i]+sub.mats[[n.age]]$R%*%Hs%*%Psucc%*%dPa_grov_int%*%La[,,i]
      dmala_grom_slp[,,i]=dmala_gro%*%dPa_grom_slp%*%La[,,i]+sub.mats[[n.age]]$R%*%Hs%*%Psucc%*%dPa_grom_slp%*%La[,,i]
      dmala_grov_slp[,,i]=dmala_gro%*%dPa_grov_slp%*%La[,,i]+sub.mats[[n.age]]$R%*%Hs%*%Psucc%*%dPa_grov_slp%*%La[,,i]
      
      D1_grom_int[,,i]=D11_gro%*%dPa_grom_int%*%La[,,i]+exp(-r*n.age)*(Fa[,,n.age]%*%Hw%*%Psucc%*%dPa_grom_int%*%La[,,i])
      D1_grov_int[,,i]=D11_gro%*%dPa_grov_int%*%La[,,i]+exp(-r*n.age)*(Fa[,,n.age]%*%Hw%*%Psucc%*%dPa_grov_int%*%La[,,i])
      D1_grom_slp[,,i]=D11_gro%*%dPa_grom_slp%*%La[,,i]+exp(-r*n.age)*(Fa[,,n.age]%*%Hw%*%Psucc%*%dPa_grom_slp%*%La[,,i])
      D1_grov_slp[,,i]=D11_gro%*%dPa_grov_slp%*%La[,,i]+exp(-r*n.age)*(Fa[,,n.age]%*%Hw%*%Psucc%*%dPa_grov_slp%*%La[,,i])
      
    }else{  # only closing term because nothing changes before n.age
      dAo_grom_int[,,i]=Fa[,,n.age]%*%Hs%*%dPa_grom_int%*%Hs%*%La[,,n.age]
      dAo_grov_int[,,i]=Fa[,,n.age]%*%Hs%*%dPa_grov_int%*%Hs%*%La[,,n.age]
      dAo_grom_slp[,,i]=Fa[,,n.age]%*%Hs%*%dPa_grom_slp%*%Hs%*%La[,,n.age]
      dAo_grov_slp[,,i]=Fa[,,n.age]%*%Hs%*%dPa_grov_slp%*%Hs%*%La[,,n.age]

      dmala_grom_int[,,i]=sub.mats[[n.age]]$R%*%Hs%*%dPa_grom_int%*%Hs%*%La[,,n.age]
      dmala_grov_int[,,i]=sub.mats[[n.age]]$R%*%Hs%*%dPa_grov_int%*%Hs%*%La[,,n.age]
      dmala_grom_slp[,,i]=sub.mats[[n.age]]$R%*%Hs%*%dPa_grom_slp%*%Hs%*%La[,,n.age]
      dmala_grov_slp[,,i]=sub.mats[[n.age]]$R%*%Hs%*%dPa_grov_slp%*%Hs%*%La[,,n.age]
      
      D1_grom_int[,,i]=exp(-r*(n.age+1))*Fa[,,n.age]%*%Hw%*%dPa_grom_int%*%Hw%*%La[,,n.age]
      D1_grov_int[,,i]=exp(-r*(n.age+1))*Fa[,,n.age]%*%Hw%*%dPa_grov_int%*%Hw%*%La[,,n.age]
      D1_grom_slp[,,i]=exp(-r*(n.age+1))*Fa[,,n.age]%*%Hw%*%dPa_grom_slp%*%Hw%*%La[,,n.age]
      D1_grov_slp[,,i]=exp(-r*(n.age+1))*Fa[,,n.age]%*%Hw%*%dPa_grov_slp%*%Hw%*%La[,,n.age]
      
    }
    
  }
  return(list(dAo_grom_int,dAo_grov_int,dAo_grom_slp,dAo_grov_slp,D1_grom_int,D1_grov_int,D1_grom_slp,D1_grov_slp,dmala_grom_int,dmala_grov_int,dmala_grom_slp,dmala_grov_slp))
}

deltau_DeltaH<- function(Sdela,D1,Fa,La,u,v,r,Ao,Ar,K, H,e,n.age,nn,max.age){
  delu=array(NA,c(nn,n.age))
  DelK=rep(0,n.age)
  DelH=array(NA,c(nn,nn,n.age))
  
  D2= exp(-r)*Fa[,,1]%*%La[,,1]
  for( ii in 2:(n.age-1)){
    D2= D2+exp(-r*ii)*ii*Fa[,,ii]%*%La[,,ii]
  }
  # closing terms
  Hw=solve(I-exp(-r)*Pa[,,n.age])
  D2= D2+n.age*Ar+exp(-r*(n.age+1))*Fa[,,n.age]%*%Hw%*%Pa[,,n.age]%*%Hw%*%La[,,n.age]
  
  Z=u%*%v
  I=diag(nn)
  for (i in (1:n.age)){
    r1=(v%*%D1[,,i]%*%u)/(v%*%D2%*%u)
    DelA=D1[,,i]-r1[1,1]*D2
    
    delu[,i]=solve(I-(Ar-Z))%*%(I-Z)%*%DelA%*%u
    
    DelK[i]=e%*%Sdela[,,i]%*%u + e%*%Ao%*%delu[,i]
    DelH[,,i]=-(DelK[i]/K[1,1])*H+Sdela[,,i]/K[1,1]
  }
  return(list(delu, DelH,DelK))
}  

# Function to calculate the elasticity of heritability from 
# delu: delta u
# DelH: delta H
# Zh: Zhat (diagonal matrix containing the mid-point values of each phenotype class)
# u: stable birth cohort phenotype distribution
Elas_heri<-function(dmala,deluDelH,e,Zh,u,H,nn, Covyxo,Varxo,mala,K,Xop){
  I=diag(nn)
  delu=deluDelH[[1]]
  DelH=deluDelH[[2]] 
  DelK=deluDelH[[3]] 
  Eb=rep(0,n.age)
  for (i in (1:n.age)){
    co1=e%*%Zh%*%DelH[,,i]%*%(I-u%*%e%*%mala/K[1,1])%*%Zh%*%u
    co2=e%*%Zh%*%H%*%(I-u%*%e%*%mala/K[1,1])%*%Zh%*%delu[,i]
    co3=e%*%Zh%*%H%*%(u%*%e%*%dmala[,,i]/K[1,1]+(delu[,i]*K[1,1]-DelK[i]*u)%*%e%*%mala/K[1,1]^2)%*%Zh%*%u
    DelCov=co1+co2-co3  # Equ. 40, 17 mars
    
    DEXoap2=e%*%((dmala[,,i]*K[1,1]-DelK[i]*mala)/K[1,1]^2)%*%Zh%*%Zh%*%u+e%*%(mala/K[1,1])%*%Zh%*%Zh%*%delu[,i]
    DEXop=e%*%((dmala[,,i]*K[1,1]-DelK[i]*mala)/K[1,1]^2)%*%Zh%*%u+e%*%(mala/K[1,1])%*%Zh%*%delu[,i]
    DelVar=DEXoap2[1,1]-2*DEXop[1,1]*Xop[1,1] # Equ. 41, 17 mars
    
    Eb[i]=DelCov[1,1]/ Covyxo-DelVar/Varxo  # Equ. 30, 14 November
  }
  return(Eb)
}

# create megamatrix for phenotype-/age-structured IPM
megamatrix <- function(subMats, nn, n.age) {
  ## function to populate the megamatrix for a phenotype age structured IPM
  MM <- matrix(0, nn * n.age, nn * n.age)
  for (j in 1:n.age) {
    MM[1:nn, (j*nn-(nn-1)):(j*nn)] <- subMats[[j]]$D %*% subMats[[j]]$R
    if (j <n.age) {
      MM[(j * nn + 1):((j+1)*nn), ((j*nn) - (nn-1)):(j*nn)] <- 
        subMats[[j]]$G %*% subMats[[j]]$S
    } else {
      MM[((j-1) * nn + 1):(j*nn), ((j*nn) - (nn-1)):(j*nn)] <- 
        subMats[[j]]$G %*% subMats[[j]]$S
    }
  }
  return(MM)
}

