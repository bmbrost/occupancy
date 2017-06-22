###
### Code to simulate community occupancy data where the mean
### probabilities of detection and occurrence are correlated.
### See Royle and Dorazio (2008) and Dorazio et al. (2011)
### for further details.
###

setwd("~/git/Occupancy/")

rm(list=ls())

# options(error=recover)
# options(error=stop)

library(mvtnorm)

expit <- function(logit){
		exp(logit)/(1+exp(logit)) 
	}
	
logit <- function(expit){
	log(expit/(1-expit))
}


#########################################################################
###  Simulate community occupancy data
#########################################################################

n <- 15  # total number of species in community
R <- 45  # total number of sites
J <- rpois(R,40)  # number of surveys conducted  at each site
hist(J)


###
### Simulate correlated mean probabilities of detection and occurrence (intercepts)
###

mu.0 <- c(-0.5,0.5)  # mean of alpha_0 and beta_0
rho <- 0.5  # correlation between intercepts
Sigma <- diag(2)*c(0.75,0.5)  # variance-covariance of intercepts
Sigma[1,2] <- Sigma[2,1] <- Sigma[1,1]*Sigma[2,2]*rho

theta.0 <- rmvnorm(n,mu.0,Sigma)  # alpha_0, beta_0
alpha.0 <- theta.0[,1]
beta.0 <- theta.0[,2]


###
### Occupancy design matrix and coefficients
###

# Occupancy varies by site and species
X <- matrix(cbind(rnorm(R),rnorm(R)),R,2)  # site-level covariates
qX <- ncol(X)

mu.beta <- c(-0.5,1)  # mean of species-level occupancy coefficients
sigma.beta <- 0.5  # standard deviation of speciesl level occupancy coefficients
beta <- matrix(c(rnorm(n,mu.beta[1],sigma.beta),  # species-level coefficients 
	rnorm(n,mu.beta[2],sigma.beta)),2,n,byrow=TRUE)  
X <- matrix(X[,1],R,1)
beta <- matrix(beta[1,],1,n)

psi <- sapply(1:R,function(x) pnorm(X[x,]%*%beta+beta.0))  # occupancy probability by site and species
hist(psi[,])


###
### Detection design matrix and coefficients
###

# Detection varies by species only
p <- expit(alpha.0)  # detection probability by species
hist(p)


###
### Latent occupancy state
###

z <- apply(psi,2,function(x) rbinom(n,1,x))
# boxplot(c(psi),c(z))


###
### Observations
###

Y <- sapply(1:R,function(x) rbinom(n,J[x],p*z[,x]))
# idx <- which(z==0)
# Y[idx] <- NA
# plot(p,colMeans(sapply(1:n,function(x) Y[x,]/J),na.rm=TRUE));abline(a=0,b=1,lty=2)

idx <- 10
round(cbind(psi=psi[,idx],z=z[,idx],p=c(p),fo=Y[,idx]/J[idx]),2)


#########################################################################
###  Fit community occupancy model assuming correlated intercepts
############################################################################

###
### Fit model
###

source("community/occ.community.correlated.mcmc.R")
priors <- list(S0=diag(2),nu=3,sigma.0=1.5,sigma.mu.beta=1.5,r=2,q=1)
# hist(sqrt(1/rgamma(1000,1,,2)))
start <- list(alpha.0=alpha.0,beta.0=beta.0,beta=beta,z=ifelse(Y>0,1,0),
	mu.0=mu.0,mu.beta=rep(0,qX),Sigma=diag(2),sigma.beta=1)
tune <- list(alpha.0=0.15,beta.0=0.25)
out1 <- occ.community.correlated.mcmc(Y,J,W=NULL,X,priors,start,tune,n.mcmc=10000,adapt=TRUE)
out1$tune
out1$keep


###
### Examine output
###

idx <- 11
matplot(cbind(out1$alpha.0[,idx],out1$beta.0[,idx]),type="l",lty=1)
abline(h=c(alpha.0[idx],beta.0[idx]),col=1:2,lty=2)
matplot(t(out1$beta[idx,,]),type="l");abline(h=beta[,idx],col=1:qX,lty=2)

matplot(out1$mu.0,type="l");abline(h=mu.0,col=1:2,lty=2)
abline(h=c(mean(alpha.0),mean(beta.0)),col=3,lty=2)

matplot(out1$mu.beta,type="l");abline(h=mu.beta,col=1:qX,lty=2)
abline(h=rowMeans(beta),col=3,lty=2)

hist(out1$Sigma[1,1,],breaks=50);abline(v=Sigma[1,1],lty=2,col=2)
abline(v=var(alpha.0),col=3,lty=2)
hist(out1$Sigma[2,2,],breaks=50);abline(v=Sigma[2,2],lty=2,col=2)
abline(v=var(beta.0),col=3,lty=2)
hist(out1$Sigma[1,2,],breaks=50);abline(v=Sigma[1,2],lty=2,col=2)
abline(v=cov(beta.0,alpha.0),col=3,lty=2)

hist(out1$sigma.beta);abline(v=sigma.beta,lty=2,col=2)
abline(v=apply(beta,1,sd),col=3,lty=2)

boxplot(out1$z.mean~z)
out1$z.mean[z==0]



#########################################################################
###  Fit community occupancy model assuming independent intercepts
#########################################################################

###
### Fit model
###

X.tmp <- cbind(1,X)
beta.tmp <- rbind(beta.0,beta)

source("community/occ.community.mcmc.R")

priors <- list(sigma.mu.beta=1.5,sigma.mu.alpha=1.5,r=2,q=1)
# hist(sqrt(1/rgamma(1000,1,,2)))
start <- list(beta=beta.tmp,alpha=alpha.0,mu.beta=rep(0,ncol(X.tmp)),mu.alpha=0,z=ifelse(Y>0,1,0),
	sigma.alpha=1,sigma.beta=1)

tune <- list(alpha=0.15)
out2 <- occ.community.mcmc(Y,J,W=NULL,X.tmp,priors,start,tune,n.mcmc=10000,adapt=TRUE)
out2$tune
out2$keep


###
### Examine output
###

idx <- 1
matplot(out2$alpha[,idx],type="l");abline(h=alpha.0[idx],col=1:2,lty=2)
matplot(t(out2$beta[idx,,]),type="l");abline(h=beta.tmp[,idx],col=1:3,lty=2)


matplot(out2$mu.beta,type="l");abline(h=c(mu.0[2],mu.beta),col=1:3,lty=2)
abline(h=rowMeans(beta.tmp),col=4,lty=2)
plot(out2$mu.alpha,type="l");abline(h=mu.0[1],col=1,lty=2)
abline(h=mean(alpha.0),col=3,lty=2)

hist(out2$sigma.beta);abline(v=sigma.beta,lty=2,col=2)
hist(out2$sigma.alpha);abline(v=sqrt(Sigma[1,1]),lty=2,col=2)

boxplot(out2$z.mean~z)
out2$z.mean[z==0]

