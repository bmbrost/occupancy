setwd("~/Documents/git/Occupancy/")
rm(list=ls())

expit <- function(y){
	exp(y)/(1+exp(y)) 
}


###
### Simulate 'single-season' occupancy data
###

n <- 100  # number of individuals
J <- 8  # number of samples per individual

# Heterogeneity in occupancy
X <- matrix(cbind(1,rnorm(n)),n,2)  # design matrix for occupancy
qX <- ncol(X)
beta <- matrix(c(0,1.5),2,1)  # coefficients for occupancy
psi <- expit(X%*%beta)  # occupancy probability
hist(psi)

# Heterogeneity in detection
W <- array(1,dim=c(n,2,J))  # design matrix for detection
qW <- dim(W)[2]
for(i in 1:J){
	W[,2,i] <- rnorm(n)
}
alpha <- matrix(c(0.5,0.5),2,1)  # coefficients for detection
p <- apply(W,3,function(x) expit(x%*%alpha))  # detection probability
summary(p)

# State process and observations
z <- rbinom(n,1,psi)  # simulated occupancy state
Y <- sapply(1:J,function(x) rbinom(n,1,z*p[,x]))  # simulated observations


###
### Fit standard occupancy model
###

source("static/occ.mcmc.R")
start <- list(beta=beta,alpha=alpha,z=z)  # starting values
priors <- list(mu.beta=rep(0,qX),mu.alpha=rep(0,qW),  # prior distribution parameters
	sigma.beta=10,sigma.alpha=10)
out1 <- occ.mcmc(Y,W,X,priors,start,1000,alpha.tune=0.1,beta.tune=0.1)  # fit model

# Examine output
matplot(out1$beta,type="l");abline(h=beta,col=1:2,lty=2)  # posterior for beta
matplot(out1$alpha,type="l");abline(h=alpha,col=1:2,lty=2)  # posterior for alpha
apply(out1$beta,2,mean)  # posterior means for beta
apply(out1$alpha,2,mean)  # posterior means for alpha
boxplot(out1$z.mean~z)  # true occupancy versus estimated occupancy
barplot(table(out1$N));sum(z)  # posterior of number in 'occupied' state

