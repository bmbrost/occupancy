setwd("~/Documents/git/Occupancy/")
rm(list=ls())


###
### Simulate 'single-season' occupancy data
###

n <- 100  # number of individuals
J <- 20  # number of samples per individual

# Heterogeneity in occupancy
X <- matrix(cbind(1,rnorm(n)),n,2)  # design matrix for occupancy
qX <- ncol(X)
beta <- matrix(c(0,1.5),2,1)  # coefficients for occupancy
psi <- pnorm(X%*%beta)  # occupancy probability
hist(psi)

# Heterogeneity in detection
W <- array(1,dim=c(n,2,J))  # design matrix for detection
qW <- dim(W)[2]
for(i in 1:J){
	W[,2,i] <- rnorm(n)
}
alpha <- matrix(c(0.5,0.5),2,1)  # coefficients for detection
p <- apply(W,3,function(x) pnorm(x%*%alpha))  # detection probability
summary(p)

# State process and observations
z <- rbinom(n,1,psi)  # simulated occupancy state
Y <- sapply(1:J,function(x) rbinom(n,1,z*p[,x]))  # simulated observations

# Fit standard occupancy model
source("static/occ.probit.1.mcmc.R")
start <- list(beta=rep(0,qX),alpha=rep(0,qW))  # starting values
priors <- list(mu.beta=rep(0,qX),mu.alpha=rep(0,qW),  # prior distribution parameters
	Sigma.beta=diag(qX)*10,Sigma.alpha=diag(qW)*10)
out1 <- occ.probit.1.mcmc(Y,W,X,priors,start,10000)  # fit model

# Examine output
matplot(out1$beta.save,type="l");abline(h=beta,col=1:2,lty=2)  # posterior for beta
matplot(out1$alpha.save,type="l");abline(h=alpha,col=1:2,lty=2)  # posterior for alpha
apply(out1$beta.save,2,mean)  # posterior means for beta
apply(out1$alpha.save,2,mean)  # posterior means for alpha
boxplot(out1$z.mean~z)  # true occupancy versus estimated occupancy
barplot(table(out1$N.save));sum(z)  # posterior of number in 'occupied' state


###
### Simulate 'single-season' occupancy with false positives
###

pi <- 0.09  # probability of false positive
controls <- rbinom(50,1,pi)  # negative control data set
controls <- list(positive=sum(controls),N=length(controls))  # summarize negative controls

# Add false positives to dataset
Y.tilde <- Y
z0 <- which(z==0)
Y.tilde[z0,] <- rbinom(J*length(z0),1,pi)

# Fit false-positive occupancy model
source("fp/occ.probit.fp.mcmc.R")
start <- list(beta=rep(0,qX),alpha=rep(0,qW))  # starting values
priors <- list(mu.beta=rep(0,qX),mu.alpha=rep(0,qW),  # prior distribution parameters
	Sigma.beta=diag(qX)*10,Sigma.alpha=diag(qW)*10,a=1,b=1)
out2 <- occ.probit.fp.mcmc(Y.tilde,W,X,controls,priors,start,10000)  # fit model

# Examine output
matplot(out2$beta.save,type="l");abline(h=beta,col=1:2,lty=2)  # posterior for beta
matplot(out2$alpha.save,type="l");abline(h=alpha,col=1:2,lty=2)  # posterior for alpha
apply(out2$beta.save,2,mean)  # posterior means for beta
apply(out2$alpha.save,2,mean)  # posterior means for alpha
boxplot(out2$z.mean~z)  # true occupancy versus estimated occupancy
barplot(table(out2$N.save));sum(z)  # posterior of number in 'occupied' state
hist(out2$pi,breaks=100);abline(v=pi,lty=2,col=2)  # posterior for pi

# Fit model ignoring false positives
source("static/occ.probit.1.mcmc.R")
start <- list(beta=rep(0,qX),alpha=rep(0,qW))  # starting values
priors <- list(mu.beta=rep(0,qX),mu.alpha=rep(0,qW),  # prior distribution parameters
	Sigma.beta=diag(qX)*10,Sigma.alpha=diag(qW)*10)
out3 <- occ.probit.1.mcmc(Y.tilde,W,X,priors,start,10000)  # fit model


###
### Compare models
###

library(lattice)

n.mcmc <- 10000
est <- data.frame(post=c(
	out1$beta[,1],out2$beta[,1],out3$beta[,1],out1$beta[,2],out2$beta[,2],out3$beta[,2],
	out1$alpha[,1],out2$alpha[,1],out3$alpha[,1],out1$alpha[,2],out2$alpha[,2],out3$alpha[,2]),
	param=c(rep("beta0",n.mcmc*3),rep("beta1",n.mcmc*3),
	rep("alpha0",n.mcmc*3),rep("alpha1",n.mcmc*3)),
	model=rep(c(rep("no fp",n.mcmc),rep("fp",n.mcmc),rep("ignore fp",n.mcmc)),4))
est$model <- ordered(est$model,levels=c("no fp","fp","ignore fp"))

bwplot(post~model|param,data=est,panel=function(x,y,...){
	panel.violin(x,y,col="lightgray",...)		
	panel.abline(h=c(alpha,beta)[panel.number()],lty=2)
})
	
