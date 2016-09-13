setwd("~/Documents/git/Occupancy/")
rm(list=ls())

library(lattice)

expit <- function(y){
	exp(y)/(1+exp(y)) 
}


###
### Simulate 'single-season' occupancy data
###

n <- 100  # number of individuals
J <- 20  # number of samples per individual

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
	# W[,3,i] <- sample(c(0,1),n,prob=c(0.7,0.3),replace=TRUE)
}
alpha <- matrix(c(0,0.5),2,1)  # coefficients for detection
# alpha <- matrix(c(-1,1),2,1)  # coefficients for detection
# alpha <- matrix(c(0,0.5,5),3,1)  # coefficients for detection
p <- apply(W,3,function(x) expit(x%*%alpha))  # detection probability
summary(p)
# summary(p[W[,3,1]==1,1])

# State process and observations
z <- rbinom(n,1,psi)  # simulated occupancy state
Y <- sapply(1:J,function(x) rbinom(n,1,z*p[,x]))  # simulated observations

# Add false positives to dataset
phi <- 0.09  # probability of false positive
Q <- matrix(rbinom(n*J,1,phi),n,J)  # false positive indicator variables
Y.tilde <- Y+Q  # add false positives to data set
Y.tilde[Y.tilde==2] <- 1


###
### Fit standard occupancy model to dataset without false positives
###

source("static/occ.mcmc.R")
start <- list(beta=beta,alpha=alpha)  # starting values
priors <- list(mu.beta=rep(0,qX),mu.alpha=rep(0,qW),  # prior distribution parameters
	sigma.beta=10,sigma.alpha=10)
tune <- list(beta=0.35,alpha=0.1)
out1 <- occ.mcmc(Y,W,X,priors,start,tune,10000,adapt=TRUE)  # fit model

# Examine output
matplot(out1$beta,type="l");abline(h=beta,col=1:2,lty=2)  # posterior for beta
matplot(out1$alpha,type="l");abline(h=alpha,col=1:2,lty=2)  # posterior for alpha
apply(out1$beta,2,mean)  # posterior means for beta
apply(out1$alpha,2,mean)  # posterior means for alpha
boxplot(out1$z.mean~z)  # true occupancy versus estimated occupancy
barplot(table(out1$N));sum(z)  # posterior of number in 'occupied' state


###
### Fit false positive occupancy model to dataset with false positives
###

source("fp/occ.fp.marginal.lik.mcmc.R")
start <- list(beta=beta,alpha=alpha,z=z,phi=phi)  # starting values
priors <- list(mu.beta=rep(0,qX),mu.alpha=rep(0,qW),  # prior distribution parameters
	sigma.beta=10,sigma.alpha=10)
tune <- list(beta=0.35,alpha=0.1)
out2 <- occ.fp.marginal.lik.mcmc(Y.tilde,W,X,priors,start,tune,10000,adapt=TRUE)  # fit model

# Examine output
matplot(out2$beta,type="l");abline(h=beta,col=1:2,lty=2)  # posterior for beta
matplot(out2$alpha,type="l");abline(h=alpha,col=1:2,lty=2)  # posterior for alpha
apply(out2$beta,2,mean)  # posterior means for beta
apply(out2$alpha,2,mean)  # posterior means for alpha
boxplot(out2$z.mean~z)  # true occupancy versus estimated occupancy
barplot(table(out2$N));sum(z)  # posterior of number in 'occupied' state
hist(out2$pi,breaks=100);abline(v=pi,lty=2,col=2)  # posterior for pi


###
### Fit standard occupancy model to dataset with false positives
###

source("static/occ.mcmc.R")
start <- list(beta=beta,alpha=alpha)  # starting values
priors <- list(mu.beta=rep(0,qX),mu.alpha=rep(0,qW),  # prior distribution parameters
	sigma.beta=10,sigma.alpha=10)
tune <- list(beta=0.35,alpha=0.1)
out3 <- occ.mcmc(Y.tilde,W,X,priors,start,tune,10000,adapt=TRUE)  # fit model

# Examine output
matplot(out3$beta,type="l");abline(h=beta,col=1:2,lty=2)  # posterior for beta
matplot(out3$alpha,type="l");abline(h=alpha,col=1:2,lty=2)  # posterior for alpha
apply(out3$beta,2,mean)  # posterior means for beta
apply(out3$alpha,2,mean)  # posterior means for alpha
boxplot(out3$z.mean~z)  # true occupancy versus estimated occupancy
barplot(table(out3$N));sum(z)  # posterior of number in 'occupied' state
hist(out3$pi,breaks=100);abline(v=pi,lty=2,col=2)  # posterior for pi


###
### Compare models
###

n.mcmc <- 10000

# Stack posteriors from all models
post <- unlist(lapply(list(out1,out2,out3),function(x) c(c(x$beta),c(x$alpha))))
est <- data.frame(post,param=rep(rep(c("beta0","beta1","alpha0","alpha1"),each=n.mcmc),3),
	model=rep(c("no fp","fp","ignore fp"), each=n.mcmc*4))
est$model <- ordered(est$model,levels=c("no fp","fp","ignore fp"))

bwplot(post~model|param,data=est,scales=list(relation="free",y=list(rot=0)),ylab="Posterior",
	panel=function(x,y,...){
		panel.violin(x,y,col="lightgray",...)		
		panel.abline(h=c(alpha,beta)[panel.number()],lty=2,col=1)
})


###
### Simple simulation study
###

sim <- data.frame("iteration"=NA,"model"=NA,"parameter"=NA,"q025"=NA,"q25"=NA,
	"q75"=NA,"q975"=NA,"mean"=NA)
sim <- sim[-1,]

sim.sum <- function(out,mod,i){
	beta.sum <- t(apply(out$beta,2,function(x) 
		c(quantile(x,c(0.025,0.25,0.75,0.975)),mean(x))))
	alpha.sum <- t(apply(out$alpha,2,function(x) 
		c(quantile(x,c(0.025,0.25,0.75,0.975)),mean(x)))) 
	tmp <- data.frame(i,mod,c("beta0","beta1","alpha0","alpha1"),rbind(beta.sum,alpha.sum))
	names(tmp) <- c("iteration","model","parameter","q025","q25","q75","q975","mean")
	tmp
}

n.mcmc <- 10000  # number of MCMC iterations to perform for each model
n.sim <- 1000
for(i in 1:n.sim){

	cat("\n")
	print(i)
	cat("\n")

	###
	### Simulate 'single-season' occupancy data
	###
	
	n <- 100  # number of individuals
	J <- 20  # number of samples per individual
	
	# Heterogeneity in occupancy
	X <- matrix(cbind(1,rnorm(n)),n,2)  # design matrix for occupancy
	qX <- ncol(X)
	beta <- matrix(c(0,1.5),2,1)  # coefficients for occupancy
	psi <- expit(X%*%beta)  # occupancy probability
	# hist(psi)
	
	# Heterogeneity in detection
	W <- array(1,dim=c(n,2,J))  # design matrix for detection
	qW <- dim(W)[2]
	for(j in 1:J){
		W[,2,j] <- rnorm(n)
	}

	alpha <- matrix(c(-1,1),2,1)  # coefficients for detection; unbiased beta p\in{0.02,0.9}
	# alpha <- matrix(c(0.25,0.25),2,1)  # coefficients for detection; biased beta p\in{0.3,0.7}
	# alpha <- matrix(c(0,-0.5),2,1)  # coefficients for detection; biased beta p\in{0.2,0.8}
	# alpha <- matrix(c(2,0.5),2,1)  # coefficients for detection; unbiased beta p\in{0.65,0.97}
	# alpha <- matrix(c(1,1),2,1)  # coefficients for detection; unbiased beta p\in{0.2,0.95}
	p <- apply(W,3,function(x) expit(x%*%alpha))  # detection probability
	# summary(p)
	
	# State process and observations
	z <- rbinom(n,1,psi)  # simulated occupancy state
	Y <- sapply(1:J,function(x) rbinom(n,1,z*p[,x]))  # simulated observations
	

	###
	### Add false positives to dataset
	###

	phi <- 0.09  # probability of false positive
	Q <- matrix(rbinom(n*J,1,phi),n,J)  # false positive indicator variables
	Y.tilde <- Y+Q  # add false positives to data set
	Y.tilde[Y.tilde==2] <- 1

	
	###
	### Fit standard occupancy model to dataset without false positives
	###
	
	source("static/occ.mcmc.R")
	start <- list(beta=beta,alpha=alpha)  # starting values
	priors <- list(mu.beta=rep(0,qX),mu.alpha=rep(0,qW),  # prior distribution parameters
		sigma.beta=2,sigma.alpha=2)
	tune <- list(beta=0.45,alpha=0.1)
	out1 <- occ.mcmc(Y,W,X,priors,start,tune,n.mcmc,adapt=TRUE)  # fit model
	

	###
	### Fit false-positive occupancy model to dataset with false positives
	###
	
	source("fp/occ.fp.marginal.lik.mcmc.R")
	start <- list(beta=beta,alpha=alpha,z=z,phi=phi)  # starting values
	priors <- list(mu.beta=rep(0,qX),mu.alpha=rep(0,qW),  # prior distribution parameters
		sigma.beta=2,sigma.alpha=2)
	tune <- list(beta=0.45,alpha=0.1)
	out2 <- occ.fp.marginal.lik.mcmc(Y.tilde,W,X,priors,start,tune,n.mcmc,adapt=TRUE)  # fit model
		
	
	###
	### Fit standard occupancy model to dataset with false positives (ignore false positives)
	###

	source("static/occ.mcmc.R")
	start <- list(beta=beta,alpha=alpha)  # starting values
	priors <- list(mu.beta=rep(0,qX),mu.alpha=rep(0,qW),  # prior distribution parameters
		sigma.beta=2,sigma.alpha=2)
	tune <- list(beta=0.45,alpha=0.1)
	out3 <- occ.mcmc(Y.tilde,W,X,priors,start,tune,n.mcmc,adapt=TRUE)  # fit model

	
	###
	### Other modeling options to pursue
	###
	
	# Threshold values of y_i to determine unoccupied individuals
	# Threshold values of y_i to create presence/absence data and model with logistic regression
	# Others...

	###
	### Summarize results
	###

	sim <- rbind(sim,sim.sum(out1,mod="no fp",i))	
	sim <- rbind(sim,sim.sum(out2,mod="fp",i))	
	sim <- rbind(sim,sim.sum(out3,mod="ignore fp",i))	
}
sim$model <- ordered(sim$model,levels=c("no fp","fp","ignore fp"))

bwplot(mean~model|parameter,data=sim,scales=list(relation="free",y=list(rot=0)),ylab="Posterior",
	panel=function(x,y,...){
		panel.violin(x,y,col="lightgray",...)		
		panel.abline(h=c(alpha,beta)[panel.number()],lty=2,col=1)
})

# Calculate coverage of 95% credible intervals
cov <- data.frame(model=sim$model,parameter=sim$parameter,cov=NA)
param.tmp <- rep(c(beta,alpha),n.sim*3)
cov$cov <- ifelse(sim$q025<param.tmp&sim$q975>param.tmp,1,0)

tapply(cov$cov,list(cov$model,cov$parameter),sum)/n.sim  # coverage
