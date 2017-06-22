###
### Code to simulate community occupancy data with optional
### correlated detection and occurrence probabilities, and to fit  
### the uncorrelated and correlated models. See Royle and Dorazio (2008)
### and Dorazio et al. (2011) for further details.
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

###
### General characteristics of simulated dataset
###

n <- 100  # total number of species in community
R <- 100  # total number of sites
J <- rpois(R,40)  # number of surveys conducted  at each site
hist(J)


###
### Simulate detection and occupancy coefficients
###

# Mean, standard deviation, and correlation of detection and occupancy coefficients
mu.alpha <- -0.5  # mean of detection intercept
mu.beta <- c(0.5,1)  # mean of occupancy coefficients
sigma.alpha <- 0.75  # standard deviation of species-level detection intercept
sigma.beta <- 0.5  # standard deviation of species-level occupancy coefficients
rho <- 0  # uncorrelated intercepts
rho <- 0.5  # correlation between intercepts

# Variance-covariance matrix
Sigma <- diag(c(sigma.alpha,rep(sigma.beta,2)))  # variance of detection and occupancy coefficients
Sigma[1,2] <- Sigma[2,1] <- Sigma[1,1]*Sigma[2,2]*rho  # covariance of intercepts

# Simulated coefficients
theta <- rmvnorm(n,c(mu.alpha,mu.beta),Sigma)  # alpha_0, beta_0, and beta_1
alpha <- matrix(theta[,1],1,n)  # alpha_0
beta <- matrix(theta[,2:3],2,n)  # beta_0 and beta_1


###
### Detection design matrix and coefficients
###

# Detection varies by species only
W <- matrix(1,n,1)  # intercept-only
qW <- ncol(W)
p <- expit(alpha)  # detection probability by species
hist(p)


###
### Occupancy design matrix and coefficients
###

# Occupancy varies by site and species
X <- matrix(cbind(1,rnorm(R)),R,2)  # site-level covariates
qX <- ncol(X)
psi <- sapply(1:R,function(x) pnorm(X[x,]%*%beta))  # occupancy probability by site and species
hist(psi[,])


###
### Latent occupancy state
###

z <- apply(psi,2,function(x) rbinom(n,1,x))
boxplot(c(psi)~c(z))


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
###  Fit community occpuancy model assuming independent detection and occupancy intercepts
#########################################################################

###
###  Fit model
###

source("community/occ.community.mcmc.R")
priors <- list(sigma.mu.alpha=1.5,sigma.mu.beta=1.5,r=2,q=1)
# hist(sqrt(1/rgamma(1000,1,,2)))
start <- list(alpha=alpha,beta=beta,mu.alpha=rep(0,qW),mu.beta=rep(0,qX),z=z,
	sigma.alpha=1.5,sigma.beta=1.5)
tune <- list(alpha=0.15)
out1 <- occ.community.mcmc(Y,J,W,X,priors,start,tune,n.mcmc=10000,adapt=TRUE)
out1$tune
out1$keep


###
### Examine output
###

idx <- 100
dev.new()
matplot(out1$alpha[,,idx],type="l",lty=1);abline(h=alpha[,idx],col=1:qW,lty=2)
matplot(out1$beta[,,idx],type="l",lty=1);abline(h=beta[,idx],col=1:qX,lty=2)

matplot(out1$mu.beta,type="l");abline(h=mu.beta,col=1:2,lty=2)
abline(h=rowMeans(beta),col=3:4,lty=2)
plot(out1$mu.alpha,type="l");abline(h=mu.alpha,col=1,lty=2)
abline(h=mean(alpha),col=3,lty=2)

hist(out1$sigma.beta);abline(v=sigma.beta,lty=2,col=2);abline(v=apply(beta,1,var),col=3,lty=2)
hist(out1$sigma.alpha);abline(v=sigma.alpha,lty=2,col=2);abline(v=apply(alpha,1,var),col=3,lty=2)

boxplot(out1$z.mean~z)
out1$z.mean[z==0]


#########################################################################
###  Fit community occupancy model assuming correlated detection and occupancy intercepts
############################################################################

###
### Fit model
###

source("community/occ.community.correlated.mcmc.R")
priors <- list(S0=diag(2),nu=3,sigma.mu.0=1.5,sigma.mu.beta=1.5,r=2,q=1)
# hist(sqrt(1/rgamma(1000,1,,2)))
start <- list(alpha=alpha,beta=beta,z=z,mu.alpha=rep(0,qW),mu.beta=rep(0,qX),
	Sigma=diag(2)*1.5,sigma.beta=1)
tune <- list(alpha=0.15,beta=0.3)
out2 <- occ.community.correlated.mcmc(Y,J,W,X,priors,start,tune,n.mcmc=10000,adapt=TRUE)
out2$tune
out2$keep


###
### Examine output
###

idx <- 38
matplot(out2$alpha[,,idx],type="l",lty=1);abline(h=alpha[,idx],col=1:qW,lty=2)
matplot(out2$beta[,,idx],type="l",lty=1);abline(h=beta[,idx],col=1:qX,lty=2)

matplot(out2$mu.alpha,type="l");abline(h=mu.alpha,col=1:qW,lty=2);abline(h=rowMeans(alpha),col=3,lty=2)
matplot(out2$mu.beta,type="l");abline(h=mu.beta,col=1:qX,lty=2);abline(h=rowMeans(beta),col=3,lty=2)

hist(out2$Sigma[1,1,],breaks=50);abline(v=Sigma[1,1],lty=2,col=2)
abline(v=var(alpha[1,]),col=3,lty=2)
hist(out2$Sigma[2,2,],breaks=50);abline(v=Sigma[2,2],lty=2,col=2)
abline(v=var(beta[1,]),col=3,lty=2)
hist(out2$Sigma[1,2,],breaks=50);abline(v=Sigma[1,2],lty=2,col=2)
abline(v=cov(alpha[1,],beta[1,]),col=3,lty=2)

hist(out2$sigma.beta);abline(v=sigma.beta,lty=2,col=2)
abline(v=sd(beta[-1,]),col=3,lty=2)

boxplot(out2$z.mean~z)
out2$z.mean[z==0]


