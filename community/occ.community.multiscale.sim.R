###
### Code to simulate multiscale community occupancy data where detection and 
### occurrence probabilities are assumed independent.
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
###  Simulate multiscale community occupancy data
#########################################################################

###
### General characteristics of simulated dataset
###

n <- 15  # total number of species in community
R <- 5  # total number of units
M <- rpois(R,10)  # number subunits per unit
M
# T <- 1  # each unit surveyed once
T <- 20  # units surveyed repeatedly through time, where T is the number of sampling periods
	# i.e., the temporal covariate model of Royle and Dorazio (2008; Page 396)
J <- rpois(sum(M)*T,50)  # number of surveys conducted at each subunit and sampling period
hist(J);table(J)


###
### Simulate detection, use, and occupancy coefficients
###

# Mean value of species-specific coefficients
mu.alpha <- -1  # mean of detection intercept
mu.gamma <- c(1,1.5)  # mean of use coefficients
mu.beta <- c(0.5,1)  # mean of occupancy coefficients

# Standard deviation among species-specific coefficients
sigma.alpha <- 1  # standard deviation of detection intercepts
sigma.gamma <- 0.5  # standard deviation of use coefficients
sigma.beta <- 0.75  # standard deviation of occupancy coefficients

# Correlation between intercepts
rho <- 0  # uncorrelated intercepts
# rho <- 0.5  # correlation between intercepts

# Simulate coefficients by species
alpha <- matrix(rnorm(n,mu.alpha,sigma.alpha),1,n)  # detection coefficients
gamma <- t(rmvnorm(n,mu.gamma,sigma.gamma^2*diag(2)))  # use coefficients
beta <- t(rmvnorm(n,mu.beta,sigma.beta^2*diag(2)))  # occupancy coefficients

# # Variance-covariance matrix
# Sigma <- diag(c(sigma.alpha,rep(sigma.gamma,2)))  # variance of detection and occupancy coefficients
# Sigma[1,2] <- Sigma[2,1] <- Sigma[1,1]*Sigma[2,2]*rho  # covariance of intercepts

# # Simulated coefficients
# tmp <- rmvnorm(n,c(mu.alpha,mu.gamma),Sigma)  # alpha_0, gamma_0, and gamma_1
# alpha <- matrix(tmp[,1],1,n,byrow=TRUE)  # alpha_0
# gamma <- matrix(tmp[,2:3],2,n,byrow=TRUE)  # gamma_0 and gamma_1


###
### Occupancy: Design matrix, probabilities, and latent state
###

# Design matrix - covariates by unit and primary sampling period
X <- matrix(1,nrow=R*T,ncol=2)
X[,2] <- scale(rnorm(R*T,rep(seq(-5,5,length.out=T),R)))
qX <- ncol(X)  # number of covariates
# plot(1:nrow(X),X[,2],col=rep(1:R,each=T))  # check temporal trend

# Occupancy probability by unit, sampling period, and species (columns)
psi <- t(apply(X,1,function(x) pnorm(x%*%beta)))
# plot(psi[2,],pnorm(X[2,]%*%beta))  # check occupancy probabilities
# hist(psi)

# Assignment of rows in X to units and primary sampling periods
groups <- list(X=data.frame(unit=rep(1:R,each=T),time=rep(1:T,R)))

# Occupancy latent state
z <- t(apply(psi,1,function(x) rbinom(n,1,x)))
# boxplot(c(psi)~c(z))


###
### Use: Design matrix, probabilities, and latent state
###

# Design matrix - covariates by unit, subunit, and primary sampling period
U <- matrix(1,nrow=sum(M)*T,ncol=2)
# U[,2] <- sample(c(0,1),nrow(U),replace=TRUE)
U[,2] <- scale(rnorm(nrow(U)))
qU <- ncol(U)  # number of covariates

# Use probability by unit, subunit, primary sampling period, and species (columns)
theta <- t(apply(U,1,function(x) pnorm(x%*%gamma)))
# plot(theta[2,],pnorm(U[2,]%*%gamma))  # check probabilities of use
# hist(theta)

# Assignment of rows in U to units, subunits, and primary sampling periods
groups$U <- data.frame(unit=rep(unlist(sapply(1:R,function(x) rep(x,M[x]))),each=T),
	subunit=rep(unlist(sapply(1:R,function(x) 1:M[x])),each=T),time=rep(1:T,sum(M)))

# Create indicator variable that maps latent occupancy state (z) to use state (a)
z.map <- match(paste(groups$U$unit,groups$U$time),paste(groups$X$unit,groups$X$time))

# Use state
a <- t(apply(z[z.map,]*theta,1,function(x) rbinom(n,1,x)))
# boxplot(c(z[z.map,]*theta)~c(a))
# a[z[z.map,]==0]

###
### Detection: Design matrix and probabilities
###

# Detection varies by species only
W <- matrix(1,sum(M)*T,1)  # intercept-only
qW <- ncol(W)
p <- expit(W%*%alpha)  # detection probability by species
hist(p)


###
### Observations
###

Y <- apply(a*p,2,function(x) rbinom(sum(M)*T,J,x))
# boxplot(c(Y)~c(a))
# table(Y[a==0])
# table(Y[z[z.map,]==0])

# Check frequency of occurrence
fo <- Y/J*ifelse(a==1,1,NA)
plot(p,fo);abline(a=0,b=1,lty=2)


#########################################################################
###  Fit community occpuancy model assuming independent detection and occupancy intercepts
#########################################################################

###
###  Fit model
###

source("community/occ.community.multiscale.mcmc.R")
start <- list(alpha=alpha,beta=beta,gamma=gamma,z=z,a=a,  # starting values
	mu.beta=mu.beta,mu.gamma=mu.gamma,mu.alpha=mu.alpha,
	sigma.beta=sigma.beta,sigma.gamma=sigma.gamma,sigma.alpha=sigma.alpha)  
priors <- list(sigma.mu.alpha=1.5,sigma.mu.gamma=1.5,  # prior distribution parameters
	sigma.mu.beta=1.5,r=2,q=1)
# hist(sqrt(1/rgamma(1000,1,,2)))
tune <- list(alpha=0.05)
out1 <- occ.community.multiscale.mcmc(Y,J,groups,W,U,X,priors,start,tune,n.mcmc=10000,adapt=TRUE)  # fit model
out1$tune
out1$keep


###
### Examine output
###

idx <- 10
p[1,idx]  # generally poor mixing of gamma for low p
# cbind(z[z.map,idx],a[,idx],Y[,idx])
matplot(out1$alpha[,,idx],type="l",lty=1);abline(h=alpha[,idx],col=1:qW,lty=2)
matplot(out1$gamma[,,idx],type="l",lty=1);abline(h=gamma[,idx],col=1:qX,lty=2)
matplot(out1$beta[,,idx],type="l",lty=1);abline(h=beta[,idx],col=1:qX,lty=2)

matplot(out1$mu.alpha,type="l");abline(h=mu.alpha,col=1:qW,lty=2);abline(h=rowMeans(alpha),col=3,lty=2)
matplot(out1$mu.gamma,type="l");abline(h=mu.gamma,col=1:qU,lty=2);abline(h=rowMeans(gamma),col=3,lty=2)
matplot(out1$mu.beta,type="l");abline(h=mu.beta,col=1:qX,lty=2);abline(h=rowMeans(beta),col=3,lty=2)

hist(out1$sigma.alpha);abline(v=sigma.alpha,lty=2,col=2)
abline(v=sd(alpha),col=3,lty=2)

hist(out1$sigma.gamma);abline(v=sigma.gamma,lty=2,col=2)
abline(v=sd(gamma),col=3,lty=2)

hist(out1$sigma.beta);abline(v=sigma.beta,lty=2,col=2)
abline(v=sd(beta),col=3,lty=2)

boxplot(out1$a.mean~a)
out1$a.mean[a==0]
out1$a.mean[a==1]

boxplot(out1$z.mean~z)
out1$z.mean[z==0]
out1$z.mean[z==1]

summary(out1$richness)
hist(out1$richness[,10]);abline(v=sum(z[10,]),lty=2,col=2)

# hist(out3$Sigma[1,1,],breaks=50);abline(v=Sigma[1,1],lty=2,col=2)
# abline(v=var(alpha[1,]),col=3,lty=2)
# hist(out3$Sigma[2,2,],breaks=50);abline(v=Sigma[2,2],lty=2,col=2)
# abline(v=var(beta[1,]),col=3,lty=2)
# hist(out3$Sigma[1,2,],breaks=50);abline(v=Sigma[1,2],lty=2,col=2)
# abline(v=cov(alpha[1,],beta[1,]),col=3,lty=2)
