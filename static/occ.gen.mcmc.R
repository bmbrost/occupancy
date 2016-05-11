occ.gen.mcmc <- function(y,J,W,X,n.mcmc,alpha.tune=NULL,beta.tune=NULL){

#
#  Mevin Hooten (20111031), Last Updated: 20131029
#
#  W: covariates for detection probability (p)
#  X: covariates for occupancy probability (psi)
#
#

####
####  Libraries and Subroutines
####

expit <- function(logit){
  exp(logit)/(1+exp(logit)) 
}
logit <- function(expit){
  log(expit/(1-expit))
}

####
####  Create Variables 
####

n=length(y)
qX=dim(X)[2]
qW=dim(W)[2]
beta.save=matrix(0,qX,n.mcmc)
alpha.save=matrix(0,qW,n.mcmc)
z.mean=rep(0,n)
N.save=rep(0,n.mcmc)

####
####  Priors and Starting Values 
####

beta.mn=rep(0,qX)
alpha.mn=rep(0,qW)
beta.var=rep(1.5,qX)
alpha.var=rep(1.5,qW)
z=rep(0,n)
z[y>0]=1

beta=as.vector(glm(z ~ 0+X,family=binomial())$coefficients)
if(is.null(beta.tune)){beta.tune=.5*abs(beta)}

alpha=as.vector(glm(cbind(y[y>1],J[y>1]-y[y>1]) ~ 0+W[y>1,],family=binomial())$coefficients)
if(is.null(alpha.tune)){alpha.tune=5*abs(alpha)}

####
####  Begin MCMC Loop 
####

for(k in 1:n.mcmc){
  if(k%%1000==0) cat(k," "); flush.console()

  ####
  ####  Sample beta 
  ####
  
  beta.star=rnorm(qX,beta,beta.tune)
  mh1=sum(dbinom(z,1,expit(X%*%beta.star),log=TRUE))+sum(dnorm(beta.star,beta.mn,sqrt(beta.var),log=TRUE))
  mh2=sum(dbinom(z,1,expit(X%*%beta),log=TRUE))+sum(dnorm(beta,beta.mn,sqrt(beta.var),log=TRUE))
  mh=exp(mh1-mh2)
  if(mh > runif(1)){
    beta=beta.star
  }
  psi=expit(X%*%beta)

  ####
  ####  Sample p 
  ####

  alpha.star=rnorm(qW,alpha,alpha.tune)
  mh1=sum(dbinom(y[z==1],J[z==1],expit(W[z==1,]%*%alpha.star),log=TRUE))+sum(dnorm(alpha.star,alpha.mn,sqrt(alpha.var),log=TRUE))
  mh2=sum(dbinom(y[z==1],J[z==1],expit(W[z==1,]%*%alpha),log=TRUE))+sum(dnorm(alpha,alpha.mn,sqrt(alpha.var),log=TRUE))
  mh=exp(mh1-mh2)
  if(mh > runif(1)){
    alpha=alpha.star
  }
  p=expit(W%*%alpha)
 
  ####
  ####  Sample z 
  ####

  num.tmp=psi*(1-p)^J 
  psi.tmp=num.tmp/(num.tmp+(1-psi))
  z[y==0]=rbinom(sum(y==0),1,psi.tmp[y==0])

  ####
  ####  Save Samples 
  ####

  beta.save[,k]=beta
  alpha.save[,k]=alpha
  z.mean=z.mean+z/n.mcmc
  N.save[k]=sum(z)

}
cat("\n")

####
####  Write Output 
####

list(beta.save=beta.save,alpha.save=alpha.save,N.save=N.save,z.mean=z.mean,n.mcmc=n.mcmc)

}
