###
### Code to simulate community occupancy data with (optionally)
### correlated detection and occurrence probabilities. 
###

setwd("~/git/Occupancy/")

rm(list=ls())

# options(error=recover)
# options(error=stop)

library(mvtnorm)
library(lattice)

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

n <- 15  # total number of species in community
R <- 45  # total number of sites
# T <- 1  # each unit surveyed once
T <- 10  # units surveyed repeatedly through time, where T is the number of sampling periods
	# i.e., the temporal covariate model of Royle and Dorazio (2008; Page 396)
J <- rpois(R*T,50)  # number of surveys conducted in each unit and sampling period
hist(J)


###
### Simulate detection and occupancy coefficients
###

# Mean value of species-specific coefficients
mu.alpha <- -1  # mean of detection intercept
mu.beta <- c(0.5,1)  # mean of occupancy coefficients

# Standard deviation among species-specific coefficients
sigma.alpha <- 1  # standard deviation of detection intercepts
sigma.beta <- 0.5  # standard deviation of occupancy coefficients

# Correlation between intercepts
# rho <- 0  # uncorrelated intercepts
rho <- 0.5  # correlation between intercepts

# Variance-covariance matrix
Sigma <- diag(c(sigma.alpha,rep(sigma.beta,2)))  # variance of detection and occupancy coefficients
Sigma[1,2] <- Sigma[2,1] <- Sigma[1,1]*Sigma[2,2]*rho  # covariance of intercepts

# Simulated coefficients by species
tmp <- rmvnorm(n,c(mu.alpha,mu.beta),Sigma)  # alpha_0, beta_0, and beta_1
alpha <- matrix(tmp[,1],1,n,byrow=TRUE)  # detection coefficients
beta <- matrix(tmp[,2:3],2,n,byrow=TRUE)  # occupancy coefficients


###
### Occupancy design matrix, probabilities, and latent state
###

# Design matrix - covariates by unit and primary sampling period
X <- matrix(1,nrow=R*T,ncol=2)
X[,2] <- scale(rnorm(R*T,rep(seq(0,1,length.out=T),each=R)))
qX <- ncol(X)  # number of covariates
# plot(1:nrow(X),X[,2],col=rep(1:R,each=T))

# Occupancy probability by unit, sampling period, and species (columns)
psi <- t(apply(X,1,function(x) pnorm(x%*%beta)))
# hist(psi)

# Occupancy latent state
z <- t(apply(psi,1,function(x) rbinom(n,1,x)))
# boxplot(c(psi)~c(z))


###
### Detection: Design matrix and probabilities
###

# Detection varies by species only
W <- matrix(1,R*T,1)  # intercept-only
qW <- ncol(W)
p <- expit(W%*%alpha)  # detection probability by species
hist(p)


###
### Observations
###

Y <- apply(z*p,2,function(x) rbinom(R*T,J,x))
# boxplot(c(Y)~c(z))
# table(Y[z==0])
# table(Y[z==1])

# Check frequency of occurrence
fo <- Y/J*ifelse(z==1,1,NA)
plot(p,fo);abline(a=0,b=1,lty=2)


#########################################################################
###  Fit community occpuancy model where detection and occupancy intercepts are assumed independent
#########################################################################

###
###  Fit model
###

source("community/occ.community.mcmc.R")
priors <- list(sigma.mu.alpha=1.5,sigma.mu.beta=1.5,r=2,q=1)  # prior distribution parameters
# hist(sqrt(1/rgamma(1000,1,,2)))
start <- list(alpha=alpha,beta=beta,mu.alpha=rep(0,qW),mu.beta=rep(0,qX),z=z,  # starting values
	sigma.alpha=1.5,sigma.beta=1.5)
tune <- list(alpha=0.05)
out1 <- occ.community.mcmc(Y,J,W,X,priors,start,tune,n.mcmc=5000,adapt=TRUE)
out1$tune
out1$keep


###
### Examine output
###

idx <- 1
matplot(out1$alpha[,,idx],type="l",lty=1);abline(h=alpha[,idx],col=1:qW,lty=2)
matplot(out1$beta[,,idx],type="l",lty=1);abline(h=beta[,idx],col=1:qX,lty=2)

plot(out1$mu.alpha,type="l");abline(h=mu.alpha,col=1,lty=2);abline(h=mean(alpha),col=3,lty=2)
matplot(out1$mu.beta,type="l");abline(h=mu.beta,col=1:2,lty=2);abline(h=rowMeans(beta),col=3,lty=2)

hist(out1$sigma.beta);abline(v=sigma.beta,lty=2,col=2);abline(v=sd(beta),col=3,lty=2)
hist(out1$sigma.alpha);abline(v=sigma.alpha,lty=2,col=2);abline(v=sd(alpha),col=3,lty=2)

boxplot(out1$z.mean~z)
out1$z.mean[z==0]

# Diversity by sampling period
diversity <- data.frame(richness=c(out1$richness),hill0=c(out1$hill0),hill1=c(out1$hill1),
	hill2=c(out1$hill2),time=rep(1:T,each=R*out1$n.mcmc))
bwplot(time~richness,data=diversity)  # richness based on latent occurence state (z)
bwplot(time~hill0,data=diversity)  # richness based on Hill numbers (q=0)
bwplot(time~hill1,data=diversity)  # Shannon diversity based on Hill numbers (q=1)
bwplot(time~hill2,data=diversity)  # Simpson diversity based on Hill numbers (q=2)


#########################################################################
###  Fit community occpuancy model where the covariance between detection and occupancy is estimated
#########################################################################

###
###  Fit model
###

source("community/occ.community.correlated.mcmc.R")
priors <- list(S0=diag(2),nu=3,sigma.mu.0=1.5,sigma.mu.beta=1.5,r=2,q=1)  # prior distribution parameters
# hist(sqrt(1/rgamma(1000,1,,2)))
start <- list(alpha=alpha,beta=beta,z=z,mu.alpha=rep(0,qW),mu.beta=rep(0,qX),  # starting values
	Sigma=diag(2)*1.5,sigma.beta=1)
tune <- list(alpha=0.05,beta=0.1)
out2 <- occ.community.correlated.mcmc(Y,J,W,X,priors,start,tune,n.mcmc=5000,adapt=TRUE)
out2$tune
out2$keep


###
### Examine output
###

idx <- 5
matplot(out2$alpha[,,idx],type="l",lty=1);abline(h=alpha[,idx],col=1:qW,lty=2)
matplot(out2$beta[,,idx],type="l",lty=1);abline(h=beta[,idx],col=1:qX,lty=2)

matplot(out2$mu.alpha,type="l");abline(h=mu.alpha,col=1:qW,lty=2);abline(h=rowMeans(alpha),col=3,lty=2)
matplot(out2$mu.beta,type="l");abline(h=mu.beta,col=1:qX,lty=2);abline(h=rowMeans(beta),col=3,lty=2)

hist(out2$Sigma[1,1,],breaks=50);abline(v=Sigma[1,1],lty=2,col=2);abline(v=var(alpha[1,]),col=3,lty=2)
hist(out2$Sigma[2,2,],breaks=50);abline(v=Sigma[2,2],lty=2,col=2);abline(v=var(beta[1,]),col=3,lty=2)
hist(out2$Sigma[1,2,],breaks=50);abline(v=Sigma[1,2],lty=2,col=2);abline(v=cov(alpha[1,],beta[1,]),col=3,lty=2)

hist(out2$sigma.beta);abline(v=sigma.beta,lty=2,col=2);abline(v=sd(beta[-1,]),col=3,lty=2)

boxplot(out2$z.mean~z)
hist(out2$z.mean[z==0])
hist(out2$z.mean[z==1])

# Diversity by sampling period
diversity <- data.frame(richness=c(out2$richness),hill0=c(out2$hill0),hill1=c(out2$hill1),
	hill2=c(out2$hill2),time=rep(1:T,each=R*out2$n.mcmc))
bwplot(time~richness,data=diversity)  # richness based on latent occurence state (z)
bwplot(time~hill0,data=diversity)  # richness based on Hill numbers (q=0)
bwplot(time~hill1,data=diversity)  # Shannon diversity based on Hill numbers (q=1)
bwplot(time~hill2,data=diversity)  # Simpson diversity based on Hill numbers (q=2)




