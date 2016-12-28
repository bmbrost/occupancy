setwd("~/Documents/git/Occupancy/")
rm(list=ls())

library(lattice)

expit <- function(y){
	exp(y)/(1+exp(y)) 
}


###
### Simulate multi-scale occupancy data
###

N <- 100  # number of sample units
J <- 8  # number of subunits per sample unit
K <- 5 # number of replicates per subunit

# Heterogeneity in occupancy (unit level)
X <- matrix(cbind(1,rnorm(N)),N,2)  # design matrix for occupancy
qX <- ncol(X)
beta <- matrix(c(0.95,1.5),2,1)  # coefficients for occupancy
psi <- expit(X%*%beta)  # occupancy probability
hist(psi)

# Heterogeneity in use (subunit level)
U <- cbind(1,rnorm(N*J))  # ordered by unit then subunit
qU <- ncol(U)
gamma <- matrix(c(-0.5,1),2,1)  # coefficients for use
theta <- expit(U%*%gamma)
hist(theta)

# Heterogeneity in detection
W <- cbind(1,rnorm(K*J*N))  # ordered by unit then subunit then replicate
qW <- ncol(W)
alpha <- matrix(c(0.25,1),2,1)  # coefficients for detection
p <- expit(W%*%alpha)  # detection probability
hist(p)

# Assignment of rows in X, U, and W to units, subunits, and replicates
groups <- list(X=data.frame(unit=1:N),U=data.frame(unit=rep(1:N,each=J),subunit=rep(1:J,N)),
	W=data.frame(unit=rep(1:N,each=J*K),subunit=rep(rep(1:J,each=K),N),replicate=rep(1:K,J*N)))
# groups$replicate <- ave(groups[,1],groups[,1],groups[,2],FUN=seq_along)

# Create indicator variable that maps latent 'occupancy' state (z) to 'use' state (a)
z.map <- match(groups$U$unit,groups$X$unit)

# Create indicator variable that maps latent 'use' state (a) to observations (y)
a.map <- match(paste(groups$W$unit,groups$W$subunit),paste(groups$U$unit,groups$U$subunit))

# Simulate state process, use, and observations
z <- rbinom(N,1,psi)  # latent occupancy state
a <- rbinom(N*J,1,z[z.map]*theta)  # use state
y <- rbinom(N*J*K,1,a[a.map]*p)  # observations

# Examine simulated data
table(tapply(a,z.map,sum)[z==1])  # number of "used" subunits across occupied units
table(tapply(a,z.map,sum)[z==0])  # number of "used" subunits across unoccupied units

table(tapply(y,a.map,sum)[a==1])  # number of detections across "used" subunits
table(tapply(y,a.map,sum)[a==0])  # number of detections across "unused" subunits


###
### Fit multiscale occupancy model
###

source("multiscale/occ.multiscale.mcmc.R")
start <- list(z=z,a=a,beta=beta,gamma=gamma,alpha=alpha)  # starting values
priors <- list(mu.beta=rep(0,qX),mu.gamma=rep(0,qU),  # prior distribution parameters
	mu.alpha=rep(0,qW),sigma.beta=2,sigma.gamma=2,sigma.alpha=2)  
tune <- list(beta=0.7,gamma=0.35,alpha=0.2)
out1 <- occ.multiscale.mcmc(y,groups,W,U,X,priors,start,tune,5000,adapt=TRUE)  # fit model

# Examine output
matplot(out1$beta,type="l");abline(h=beta,col=1:2,lty=2)  # posterior for beta
matplot(out1$gamma,type="l");abline(h=gamma,col=1:2,lty=2)  # posterior for gamma
matplot(out1$alpha,type="l");abline(h=alpha,col=1:2,lty=2)  # posterior for alpha
apply(out1$beta,2,mean)  # posterior means for beta
apply(out1$gamma,2,mean)  # posterior means for gamma
apply(out1$alpha,2,mean)  # posterior means for alpha
boxplot(out1$z.mean~z)  # true occupancy versus estimated occupancy
boxplot(out1$a.mean~a)  # true occupancy versus estimated occupancy
