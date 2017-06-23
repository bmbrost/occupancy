###
### Code to simulate community occupancy data where sites are surveyed
### repeatedly over time, and detection and occurrence probabilities are
### optionally correlated. Also includes code to fit model with temporally-
### varying covariates, i.e., the temporal covariate model of Royle 
###	and Dorazio (2008). 
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

n <- 15  # total number of species in community
R <- 45  # total number of sites
T <- 10  # number of primary sampling periods (i.e., years over which sites are surveyed)
J <- matrix(rpois(R*T,10),R,T)  # number of surveys conducted  at each site and sampling period
hist(J)


###
### Simulate detection and occupancy coefficients
###

# Note: detection varies by species, occupancy varies by site and species

# Mean, standard deviation, and correlation of detection and occupancy coefficients
mu.alpha <- -1  # mean of detection intercept
mu.beta <- c(0.5,1)  # mean of occupancy coefficients
sigma.alpha <- 1  # standard deviation of species-level detection intercept
sigma.beta <- 0.5  # standard deviation of species-level occupancy coefficients
# rho <- 0  # uncorrelated intercepts
rho <- 0.5  # correlation between intercepts

# Variance-covariance matrix
Sigma <- diag(c(sigma.alpha,rep(sigma.beta,2)))  # variance of detection and occupancy coefficients
Sigma[1,2] <- Sigma[2,1] <- Sigma[1,1]*Sigma[2,2]*rho  # covariance of intercepts

# Simulated coefficients
theta <- rmvnorm(n,c(mu.alpha,mu.beta),Sigma)  # alpha_0, beta_0, and beta_1
alpha <- matrix(theta[,1],1,n,byrow=TRUE)  # alpha_0
beta <- matrix(theta[,2:3],2,n,byrow=TRUE)  # beta_0 and beta_1
apply(alpha,1,mean)
apply(beta,1,mean)


###
### Detection design matrix and coefficients
###

# Detection varies by species only
W <- matrix(1,n,1)  # intercept-only
qW <- ncol(W)
p <- expit(alpha)  # detection probability by species
hist(p)


###
### Missing data for site x sampling period combinations
###

idx <- sample(1:(T*R),50)
idx <- ifelse(1:(T*R)%in%idx,NA,1)


###
### Occupancy design matrix and coefficients
###

# Occupancy varies by site, species, and primary sampling period
X <- array(1,c(R,2,T))  # Occupancy covariates by site (rows) and sampling periods (third dimension)
X[,1,] <- X[,1,]*idx
X[,2,] <- scale(rnorm(R*T,rep(seq(0,1,length.out=T),each=R)))*idx
plot(1:length(X[,2,]),X[,2,])
qX <- ncol(X)  # number of covariates
psi <- apply(X,c(1,3),function(x) pnorm(x%*%beta))  # occupancy prob. by species, site, and sampling period

plot(psi[,4,2],pnorm(X[4,,2]%*%beta))  # check occupancy probabilities
hist(psi[,,])


###
### Latent occupancy state
###

z <- apply(psi,c(2,3),function(x) rbinom(n,1,x))
boxplot(c(psi)~c(z))


###
### Observations
###

Y <- array(0,c(n,R,T))
for (i in 1:T){
	# i <- 1
	Y[,,i] <- sapply(1:R,function(x) rbinom(n,J[x,i],p*z[,x,i]))
}

# Check simulation
fo <- Y  # frequency of occurence
idx <- which(z==0)
fo[idx] <- NA
for(i in 1:T){  # calculate frequency of occurrence
	J.tmp <- matrix(J[,i],n,R,byrow=TRUE)
	fo[,,i] <- fo[,,i]/J.tmp
}
fo <- apply(fo,3,rowMeans,na.rm=TRUE)
plot(rep(p,T),fo);abline(a=0,b=1,lty=2)

site.idx <- 2  # site indicator
sample.idx <- 3  # sampling period indicator
round(cbind(psi=psi[,site.idx,sample.idx],z=z[,site.idx,sample.idx],p=c(p),
	fo=Y[,site.idx,sample.idx]/J[site.idx,sample.idx]),2)


#########################################################################
###  Fit community occpuancy model assuming independent detection and occupancy intercepts
#########################################################################

###
###  Fit model
###

source("community/occ.community.temporal.mcmc.R")
priors <- list(S0=diag(2),nu=3,sigma.mu.0=1.5,sigma.mu.beta=1.5,r=2,q=1)
# hist(sqrt(1/rgamma(1000,1,,2)))
start <- list(alpha=alpha,beta=beta,z=z,mu.alpha=rep(0,qW),mu.beta=rep(0,qX),
	Sigma=diag(2)*1.5,sigma.beta=1)
tune <- list(alpha=0.05,beta=0.1)
out3 <- occ.community.temporal.mcmc(Y,J,W,X,priors,start,tune,n.mcmc=10000,adapt=TRUE)
out3$tune
out3$keep


###
### Examine output
###

idx <- 6
matplot(out3$alpha[,,idx],type="l",lty=1);abline(h=alpha[,idx],col=1:qW,lty=2)
matplot(out3$beta[,,idx],type="l",lty=1);abline(h=beta[,idx],col=1:qX,lty=2)

matplot(out3$mu.alpha,type="l");abline(h=mu.alpha,col=1:qW,lty=2);abline(h=rowMeans(alpha),col=3,lty=2)
matplot(out3$mu.beta,type="l");abline(h=mu.beta,col=1:qX,lty=2);abline(h=rowMeans(beta),col=3,lty=2)

hist(out3$Sigma[1,1,],breaks=50);abline(v=Sigma[1,1],lty=2,col=2)
abline(v=var(alpha[1,]),col=3,lty=2)
hist(out3$Sigma[2,2,],breaks=50);abline(v=Sigma[2,2],lty=2,col=2)
abline(v=var(beta[1,]),col=3,lty=2)
hist(out3$Sigma[1,2,],breaks=50);abline(v=Sigma[1,2],lty=2,col=2)
abline(v=cov(alpha[1,],beta[1,]),col=3,lty=2)

hist(out3$sigma.beta);abline(v=sigma.beta,lty=2,col=2)
abline(v=sd(beta[-1,]),col=3,lty=2)

boxplot(out3$z.mean~z)
out3$z.mean[z==0]
out3$z.mean[z==1]

plot(p,apply(out3$z.mean,1,mean,na.rm=TRUE))

plot(1:T,apply(out3$z.mean,3,mean,na.rm=TRUE))


