setwd("/Users/brost/Documents/research/belugas")
rm(list=ls())

# options(error=recover)
# options(error=stop)

#########################################################################
#########################################################################
##### Simulate static occupancy
#########################################################################
#########################################################################

#####
#####  Simulate two data sources, fully heterogeneous model (covariates on occupancy and detection)
#####

n <- 1000
J <- 3

# Covariates for occupancy (site)
X <- matrix(cbind(1,rnorm(n)),n,2)
beta <- matrix(c(0,3),2,1)
# psi <- exp(X%*%beta)/(1+exp(X%*%beta))
psi <- pnorm(X%*%beta)
summary(psi)

# Covariates for detectability (site and visit)
# W <- array(cbind(1,rnorm(n)),dim=c(n,2,1))
W <- array(NA,dim=c(n,2,J))
W[,1,] <- 1
for(i in 1:J){
	W[,2,i] <- rnorm(nrow(W[,,i]))
}

alpha1 <- matrix(c(0,-5),2,1)
p1 <- apply(W,3,function(x) pnorm(x%*%alpha1))
summary(p1)

alpha2 <- matrix(c(0,-1),2,1)
p2 <- apply(W,3,function(x) pnorm(x%*%alpha2))
summary(p2)

# z <- rbinom(n,1,psi)
# y <- rep(0,n)
# y[z==1] <- rbinom(sum(z==1),J,p[z==1])

z <- rbinom(n,1,psi)
Y1 <- sapply(1:J,function(x) rbinom(n,1,z*p1[,x]))
Y2 <- sapply(1:J,function(x) rbinom(n,1,z*p2[,x]))

#  Fit occupancy model to two data sources
source("static_occupancy/occ.probit.2.mcmc.R")
out1 <- occ.probit.2.mcmc(Y1,Y2,rep(J,n),rep(J,n),W,W,X,10000)
matplot(out1$beta,type="l");abline(h=beta,col=1:2,lty=2)
matplot(out1$alpha1,type="l");abline(h=alpha1,col=1:2,lty=2)
matplot(out1$alpha2,type="l");abline(h=alpha2,col=1:2,lty=2)
apply(out1$beta,2,mean)
apply(out1$alpha1,2,mean)
apply(out1$alpha2,2,mean)



#####
#####  Simulate two data sources at varying temporal scales,
#####  homogenous model (intercept only for occupany and detection)
#####

n <- 1000
J1 <- 4 # Visits per "site" for data source 1
J2 <- 2 # Visits per "site" for data source 2

# Covariates for occupancy (site)
X <- matrix(1,n,1)
beta <- matrix(0,1,1)
psi <- pnorm(X%*%beta)
summary(psi)

# Covariates for detectability (site and visit)
W1 <- array(1,dim=c(n,1,J1))
alpha1 <- matrix(0,1,1)
p1 <- apply(W1,3,function(x) pnorm(x%*%alpha1))
summary(p1)

W2 <- array(1,dim=c(n,1,J2))
alpha2 <- matrix(0,1,1)
p2 <- apply(W2,3,function(x) pnorm(x%*%alpha2))
summary(p2)

z <- rbinom(n,1,psi)
Y1 <- sapply(1:J1,function(x) rbinom(n,1,z*p1[,x]))
Y2 <- sapply(1:J2,function(x) rbinom(n,1,z*p2[,x]))

# Add NA values to simulate missing data
Y1[1000,3:4] <- NA
W1[1000,1,3:4] <- NA
Y2[1000,2] <- NA
W2[1000,1,2] <- NA

#  Fit occupancy model to two data sources
source("static_occupancy/occ.probit.2.mcmc.R")
out1 <- occ.probit.2.mcmc(Y1,Y2,W1,W2,X,500)
matplot(out1$beta,type="l");abline(h=beta,col=1:2,lty=2)
matplot(out1$alpha1,type="l");abline(h=alpha1,col=1:2,lty=2)
matplot(out1$alpha2,type="l");abline(h=alpha2,col=1:2,lty=2)
apply(out1$beta,2,mean)
apply(out1$alpha1,2,mean)
apply(out1$alpha2,2,mean)
apply(out1$beta,2,quantile,c(0.025,0.975))
apply(out1$alpha1,2,quantile,c(0.025,0.975))
apply(out1$alpha2,2,quantile,c(0.025,0.975))

plot(z,out1$z.mean)


#####
#####  Simulate two data sources, homogenous model (intercept only for occupany and detection)
#####

n <- 1000
J <- 5

# Covariates for occupancy (site)
X <- matrix(1,n,1)
beta <- matrix(0,1,1)
psi <- pnorm(X%*%beta)
summary(psi)

# Covariates for detectability (site and visit)
W <- array(1,dim=c(n,1,J))

alpha1 <- matrix(0,1,1)
p1 <- apply(W,3,function(x) pnorm(x%*%alpha1))
summary(p1)

alpha2 <- matrix(0,1,1)
p2 <- apply(W,3,function(x) pnorm(x%*%alpha2))
summary(p2)

z <- rbinom(n,1,psi)
Y1 <- sapply(1:J,function(x) rbinom(n,1,z*p1[,x]))
Y2 <- sapply(1:J,function(x) rbinom(n,1,z*p2[,x]))

#  Fit occupancy model to two data sources
source("static_occupancy/occ.probit.2.mcmc.R")
out1 <- occ.probit.2.mcmc(Y1,Y2,rep(J,n),rep(J,n),W,W,X,5000)
matplot(out1$beta,type="l");abline(h=beta,col=1:2,lty=2)
matplot(out1$alpha1,type="l");abline(h=alpha1,col=1:2,lty=2)
matplot(out1$alpha2,type="l");abline(h=alpha2,col=1:2,lty=2)
apply(out1$beta,2,mean)
apply(out1$alpha1,2,mean)
apply(out1$alpha2,2,mean)
apply(out1$beta,2,quantile,c(0.025,0.975))
apply(out1$alpha1,2,quantile,c(0.025,0.975))
apply(out1$alpha2,2,quantile,c(0.025,0.975))



#####
#####  Simulate one data source, fully heterogeneous model (covariates on occupancy and detection)
#####

n <- 1000 #Number of sites
J <- 3 #Number of visits

# Covariates for occupancy (site)
X <- matrix(cbind(1,rnorm(n)),n,2)
beta <- matrix(c(0,3),2,1)
# psi <- exp(X%*%beta)/(1+exp(X%*%beta))
psi <- pnorm(X%*%beta)
summary(psi)

# Covariates for detectability (site and visit)
# W <- array(cbind(1,rnorm(n)),dim=c(n,2,1))
W <- array(NA,dim=c(n,2,J))
W[,1,] <- 1
for(i in 1:J){
	W[,2,i] <- rnorm(nrow(W[,,i]))
}
alpha <- matrix(c(0,-3),2,1)
# p <- exp(W[,,1]%*%alpha)/(1+exp(W[,,1]%*%alpha))
# p <- pnorm(W[,,1]%*%alpha)
p <- apply(W,3,function(x) pnorm(x%*%alpha))
summary(p)

# z <- rbinom(n,1,psi)
# y <- rep(0,n)
# y[z==1] <- rbinom(sum(z==1),J,p[z==1])

z <- rbinom(n,1,psi)
Y <- sapply(1:J,function(x) rbinom(n,1,z*p[,x]))
# y <- as.numeric(apply(y,1,sum))

# Fit occupancy model to one data source
source("static_occupancy/occ.probit.1.mcmc.R")
out1 <- occ.probit.1.mcmc(Y,rep(J,n),W,X,1000)
matplot(out1$beta.save,type="l");abline(h=beta,col=1:2,lty=2)
matplot(out1$alpha.save,type="l");abline(h=alpha,col=1:2,lty=2)
apply(out1$beta.save,2,mean)
apply(out1$alpha.save,2,mean)



#####
#####  Simulate one data source, homogenous model (intercept only for occupany and detection)
#####

n <- 1000 #Number of sites
J <- 5 #Number of visits

# Covariates for occupancy (site)
X <- matrix(1,n,1)
beta <- matrix(0,1,1)
psi <- pnorm(X%*%beta)
summary(psi)

# Covariates for detectability (site and visit)
W <- array(1,dim=c(n,1,J))
alpha <- matrix(0,1,1)
p <- apply(W,3,function(x) pnorm(x%*%alpha))
summary(p)

z <- rbinom(n,1,psi)
Y <- sapply(1:J,function(x) rbinom(n,1,z*p[,x]))

# Fit occupancy model to one data source
source("static_occupancy/occ.probit.1.mcmc.R")
out1 <- occ.probit.1.mcmc(Y,rep(J,n),W,X,1000)
matplot(out1$beta.save,type="l");abline(h=beta,col=1:2,lty=2)
matplot(out1$alpha.save,type="l");abline(h=alpha,col=1:2,lty=2)
apply(out1$beta.save,2,mean)
apply(out1$alpha.save,2,mean)