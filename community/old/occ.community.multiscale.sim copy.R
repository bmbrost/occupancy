###
### Code to simulate multiscale community occupancy data where subunits are surveyed
### repeatedly over time, subunits are nested within units, and detection and 
### occurrence probabilities are optionally correlated. Also includes code to fit 
### model with temporally-varying covariates (see Royle and Dorazio 2008). 
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
T <- 20  # number of primary sampling periods (i.e., years over which units are surveyed)
J <- rpois(sum(M)*T,50)  # number of surveys conducted at each subunit and sampling period
hist(J);table(J)


###
### Simulate detection, use, and occupancy coefficients
###

# Note: detection varies by species, whereas use and occupancy varies by site and species

# Mean value of species-specific coefficients
mu.alpha <- -1  # mean of detection intercept
mu.gamma <- c(1,1.5)  # mean of use coefficients
mu.beta <- c(0.5,1)  # mean of occupancy coefficients

# Standard deviation among species-specific coefficients
sigma.alpha <- 1  # standard deviation of detection intercept
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
U[,2] <- sample(c(0,1),nrow(U),replace=TRUE)
# U[,2] <- scale(rnorm(nrow(U)))
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
W <- matrix(1,n,1)  # intercept-only
qW <- ncol(W)
p <- expit(alpha)  # detection probability by species
hist(p)


###
### Observations
###

Y <- t(sapply(1:nrow(a),function(x) rbinom(n,J[x],a[x,]*p)))  # observations
# boxplot(c(Y)~c(a))
# table(Y[a==0])
# table(Y[z[z.map,]==0])

# Check frequency of occurrence
fo <- Y/J*ifelse(a==1,1,NA)
plot(p,apply(fo,2,mean,na.rm=TRUE));abline(a=0,b=1,lty=2)


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
p[idx]  # generally poor mixing of gamma for low p
# cbind(z[z.map,idx],a[,idx],Y[,idx])
matplot(out1$alpha[,idx],type="l",lty=1);abline(h=alpha[,idx],col=1:qW,lty=2)
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

plot(p,apply(out1$z.mean,2,mean,na.rm=TRUE))

# hist(out3$Sigma[1,1,],breaks=50);abline(v=Sigma[1,1],lty=2,col=2)
# abline(v=var(alpha[1,]),col=3,lty=2)
# hist(out3$Sigma[2,2,],breaks=50);abline(v=Sigma[2,2],lty=2,col=2)
# abline(v=var(beta[1,]),col=3,lty=2)
# hist(out3$Sigma[1,2,],breaks=50);abline(v=Sigma[1,2],lty=2,col=2)
# abline(v=cov(alpha[1,],beta[1,]),col=3,lty=2)









# #########################################################################
# ###  Simulate multiscale community occupancy data
# #########################################################################

# ###
# ### General characteristics of simulated dataset
# ###

# n <- 15  # total number of species in community
# R <- 5  # total number of regions
# M <- rpois(R,10)  # number sites per region
# T <- 10  # number of primary sampling periods (i.e., years over which sites are surveyed)
# J <- matrix(rpois(sum(M)*T,10),sum(M),T)  # number of surveys conducted  at each site and sampling period
# hist(J)


# ###
# ### Simulate detection, use, and occupancy coefficients
# ###

# # Note: detection varies by species, whereas use and occupancy varies by site and species

# # Mean, standard deviation, and correlation of detection, use, and occupancy coefficients
# mu.alpha <- -1  # mean of detection intercept
# mu.gamma <- c(1,1.5)  # mean of use coefficients
# mu.beta <- c(0.5,1)  # mean of occupancy coefficients
# sigma.alpha <- 1  # standard deviation of species-level detection intercept
# sigma.gamma <- 0.5  # standard deviation of species-level use coefficients
# sigma.beta <- 0.75  # standard deviation of species-level occupancy coefficients
# rho <- 0  # uncorrelated intercepts
# # rho <- 0.5  # correlation between intercepts

# # Variance-covariance matrix
# Sigma <- diag(c(sigma.alpha,rep(sigma.gamma,2)))  # variance of detection and occupancy coefficients
# Sigma[1,2] <- Sigma[2,1] <- Sigma[1,1]*Sigma[2,2]*rho  # covariance of intercepts

# # Simulated coefficients
# tmp <- rmvnorm(n,c(mu.alpha,mu.gamma),Sigma)  # alpha_0, gamma_0, and gamma_1
# alpha <- matrix(tmp[,1],1,n,byrow=TRUE)  # alpha_0

# gamma <- array(1,dim=c(2,n,R))
# gamma[1,,] <- rnorm(n*R,mu.gamma[1],sigma.gamma)  # gamma_0
# gamma[2,,] <- rnorm(n*R,mu.gamma[2],sigma.gamma)  # gamma_1
# # gamma <- matrix(tmp[,2:3],2,n,byrow=TRUE)  # gamma_0 and gamma_1

# beta <- matrix(c(rnorm(n,mu.beta[1],sigma.beta),rnorm(n,mu.beta[2],sigma.beta)),2,n,byrow=TRUE)  # beta_0 and beta_1
# apply(alpha,1,mean)
# apply(gamma,1,mean)
# apply(beta,1,mean)

# ###
# ### Detection design matrix
# ###

# # Detection varies by species only
# W <- matrix(1,n,1)  # intercept-only
# qW <- ncol(W)
# p <- expit(alpha)  # detection probability by species
# hist(p)


# ###
# ### Missing data for site x sampling period combinations
# ###

# idx <- 1
# # idx <- sample(1:(T*R),50)
# # idx <- ifelse(1:(T*R)%in%idx,NA,1)


# ###
# ### Occupancy --- design matrix, probability, and latent state
# ###

# # Design matrix - occupancy varies by species, region, and primary sampling period
# X <- array(1,c(R,2,T))  # covariates for regions (rows) and primary sampling periods (third dimension) 
# X[,1,] <- X[,1,]*idx
# X[,2,] <- scale(rnorm(T*R,rep(seq(-5,5,length.out=T),each=R)))*idx
# # plot(1:length(X[,2,]),X[,2,])  # check trend over sampling period
# # plot(1:T,apply(X,c(2,3),mean)[2,],col=2)
# qX <- ncol(X)  # number of covariates

# # Occupancy probability by species (rows), region (columns), and sampling period (third dimension)
# psi <- apply(X,c(1,3),function(x) pnorm(x%*%beta))  
# plot(psi[,1,4],pnorm(X[1,,4]%*%beta))  # check occupancy probabilities
# hist(psi[,,])

# # Latent occupancy state
# z <- apply(psi,c(2,3),function(x) rbinom(n,1,x))
# boxplot(c(psi)~c(z))

# psi[,,1]
# z[,,1]

# ###
# ### Use --- design matrix, probability, and latent state
# ###

# # U <- sapply(1:R,function(x) array(1,c(M[x],2,T))*idx)  # list over regions
# # U <- array(1,c(T,2,sum(M)))

 # # Design matrix - use varies by species, region, site, and primary sampling period
# groups <- list(U=data.frame(region=unlist(sapply(1:R,function(x) rep(x,M[x])))))  # region indicator variable
# U <- array(1,c(sum(M),2,T))  # covariates by site (rows) and sampling periods (third dimension)
# U[,1,] <- U[,1,]*idx
# U[,2,] <- scale(rnorm(sum(M)*T,rep(seq(-5,5,length.out=T),each=sum(M))))*idx
# # plot(1:length(U[,2,]),U[,2,])  # check trend over sampling period
# # plot(1:T,apply(U,c(2,3),mean)[2,],col=2)
# qU <- ncol(U)  # number of covariates

# # Use probability by region (list), species (row), site (column), and sampling period (third dimension)
# theta <- array(0,dim=c(n,sum(M),T))
# for(i in 1:T){
	# for(j in 1:sum(M)){
		# theta[,j,i] <- expit(U[j,,i]%*%gamma[,,groups$U$region[j]])  # probability
	# }
# }

# # Laten use state
# z.full <- array(apply(z,3,function(x) x[,groups$U$region]),dim=c(n,sum(M),T))
# a <- apply(z.full*theta,c(2,3),function(x) rbinom(n,1,x))

# # a <- theta*0
# # for(i in 1:T){
	# # for(j in 1:sum(M)){
		# # z.tmp <- z[,groups$U$region[j],i]
		# # a[,j,i] <- rbinom(n,1,z.tmp*theta[,j,i])  # latent state
	# # }
# # }

# a[,,2]
# z[,,2]
	
# boxplot(c(theta)~c(a))
# boxplot(c(theta[,,1])~c(a[,,1]))
# plot(theta[,,10],a[,,10])

# # head(theta)
# # groups
# # sapply(1:sum(M),function(x) apply(U[x,,],2,function(xx) expit(t(xx)%*%gamma[,,groups$U$region[x]])))
# # sapply(1:T,function(x) apply(U[,,x],1,function(xx) expit()))
# # U
# # theta <- sapply(1:R,function(x)  
	# # apply(U[groups$U$region==x,,],c(1,3),function(xx) expit(xx%*%gamma[,,x])))

# # # Latent use state
# # a <- sapply(1:R,function(x) apply(theta[[x]],c(2,3),function(xx) rbinom(n,1,xx)))  
# # boxplot(c(theta[[2]])~c(a[[2]]))


# ###
# ### Observations
# ###

# Y <- array(0,c(n,sum(M),T))
# for (i in 1:T){
	# # i <- 1
	# for(j in 1:sum(M)){
		# Y[,,i] <- rbinom(n,J[j,i],p*a[,j,i])
	# }
# }

# # Check simulation
# fo <- Y  # frequency of occurence
# fo <- fo*ifelse(z.full==1,1,NA)
# J.full <- array(apply(J,2,function(x) matrix(x,n,sum(M),byrow=TRUE)),dim=c(n,sum(M),T))
# fo <- fo/J.full
# apply(fo[,,1],1,mean,na.rm=TRUE)
# apply(fo,3,rowMeans,na.rm=TRUE)
# rowMeans(fo,na.rm=TRUE)
# p
# z[,,1]
# J.full[,,1]

# for(i in 1:T){  # calculate frequency of occurrence
	# z.tmp <- z[,,i]
	# z.tmp <- z.tmp[,groups$U$region]
	# z.tmp <- ifelse(z.tmp==1,1,NA)
	# J.tmp <- matrix(J[,i],n,sum(M),byrow=TRUE)
	# fo[,,i] <- fo[,,i]/J.tmp*z.tmp
# }

# plot(p,rowMeans(fo,na.rm=TRUE));abline(a=0,b=1,lty=2)
# fo <- apply(fo,3,rowMeans,na.rm=TRUE)
# plot(rep(p,T),fo);abline(a=0,b=1,lty=2)

# site.idx <- 2  # site indicator
# sample.idx <- 3  # sampling period indicator
# round(cbind(psi=psi[,site.idx,sample.idx],z=z[,site.idx,sample.idx],p=c(p),
	# fo=Y[,site.idx,sample.idx]/J[site.idx,sample.idx]),2)


# #########################################################################
# ###  Fit community occpuancy model assuming independent detection and occupancy intercepts
# #########################################################################

# ###
# ###  Fit model
# ###

# source("community/occ.community.temporal.mcmc.R")
# priors <- list(S0=diag(2),nu=3,sigma.mu.0=1.5,sigma.mu.beta=1.5,r=2,q=1)
# # hist(sqrt(1/rgamma(1000,1,,2)))
# start <- list(alpha=alpha,beta=beta,z=z,mu.alpha=rep(0,qW),mu.beta=rep(0,qX),
	# Sigma=diag(2)*1.5,sigma.beta=1)
# tune <- list(alpha=0.05,beta=0.1)
# out3 <- occ.community.temporal.mcmc(Y,J,W,X,priors,start,tune,n.mcmc=10000,adapt=TRUE)
# out3$tune
# out3$keep


# ###
# ### Examine output
# ###

# idx <- 6
# matplot(out3$alpha[,,idx],type="l",lty=1);abline(h=alpha[,idx],col=1:qW,lty=2)
# matplot(out3$beta[,,idx],type="l",lty=1);abline(h=beta[,idx],col=1:qX,lty=2)

# matplot(out3$mu.alpha,type="l");abline(h=mu.alpha,col=1:qW,lty=2);abline(h=rowMeans(alpha),col=3,lty=2)
# matplot(out3$mu.beta,type="l");abline(h=mu.beta,col=1:qX,lty=2);abline(h=rowMeans(beta),col=3,lty=2)

# hist(out3$Sigma[1,1,],breaks=50);abline(v=Sigma[1,1],lty=2,col=2)
# abline(v=var(alpha[1,]),col=3,lty=2)
# hist(out3$Sigma[2,2,],breaks=50);abline(v=Sigma[2,2],lty=2,col=2)
# abline(v=var(beta[1,]),col=3,lty=2)
# hist(out3$Sigma[1,2,],breaks=50);abline(v=Sigma[1,2],lty=2,col=2)
# abline(v=cov(alpha[1,],beta[1,]),col=3,lty=2)

# hist(out3$sigma.beta);abline(v=sigma.beta,lty=2,col=2)
# abline(v=sd(beta[-1,]),col=3,lty=2)

# boxplot(out3$z.mean~z)
# out3$z.mean[z==0]
# out3$z.mean[z==1]

# plot(p,apply(out3$z.mean,1,mean,na.rm=TRUE))

# plot(1:T,apply(out3$z.mean,3,mean,na.rm=TRUE))


