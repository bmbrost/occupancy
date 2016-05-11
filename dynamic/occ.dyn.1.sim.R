setwd("/Users/brost/Documents/research/belugas")
rm(list=ls())

# options(error=recover)
# options(error=stop)

#########################################################################
#########################################################################
##### Simulate dynamic occupancy
#########################################################################
#########################################################################

#####
#####  Simulate one data source with dynamic occupancy
#####  Fully heterogeneous model (covariates on colonization, persistence, and detection)
#####  Multiple sites permitted
#####


T <- c(500,750,250,250,1000) #Number of primary periods per site
T <- c(2000) #Number of primary periods per site
# T <- c(5,10) #Number of primary periods
N <- length(T) #Number of sites
site <- unlist(sapply(1:N,function(x) rep(x,T[x])))
y_1.idx <- c(t(sapply(unique(site),function(x) min(which(site==x)))))
y_T.idx <- c(t(sapply(unique(site),function(x) max(which(site==x)))))
z_0.idx <- y_1.idx + 0:(N-1) 
z_1.idx <- z_0.idx + 1
z_T.idx <- y_T.idx + 1:N
J <- 5 #Number of visits per primary period

# Covariates for colonization
X <- matrix(cbind(1,rnorm(sum(T))),sum(T),2)
beta <- matrix(c(-0.75,0.05),2,1)
gamma <- pnorm(X%*%beta)
summary(gamma)

# Covariates for persistence
G <- X
# G <- matrix(cbind(1,rnorm(T)),T,2)
theta <- matrix(c(-0.25,0.25),2,1)
phi <- pnorm(G%*%theta)
summary(phi)

# Covariates for detectability (site and visit)
W <- array(NA,dim=c(sum(T),2,J))
W[,1,] <- 1
W[,2,] <- rnorm(sum(T)*J)
alpha <- matrix(c(-1,-0.075),2,1)
p <- apply(W,3,function(x) pnorm(x%*%alpha))
summary(p)

psi <- 0.5 #Probability of occupancy at t=0

# Simulate occupancy latent state
z <- numeric(sum(T)+N)
z[z_0.idx] <- rbinom(N,1,psi)

for(i in 1:N){
	for(j in z_1.idx[i]:z_T.idx[i]){ #Simulation z_1:Z_T
		# print(c(j,j-1,j-i))
		if(z[j-1]==1) z[j] <- rbinom(1,1,phi[j-i])
		if(z[j-1]==0) z[j] <- rbinom(1,1,gamma[j-i])
	}
}

# Simulate observed data
Y <- sapply(1:J,function(x) rbinom(sum(T),1,z[-z_0.idx]*p[,x]))

# Fit occupancy model to one data source
source("/Users/brost/Dropbox/belugas/dynamic_occupancy/occ.dyn.1.mcmc.R")
out1 <- occ.dyn.1.mcmc(site,Y,W,X,G,1000,z.true=NULL)
out1 <- occ.dyn.1.mcmc(Y,W,X,G,1000,z.true=z)

mod <- out1
mod <- out2
matplot(mod$beta,type="l");abline(h=beta,col=1:2,lty=2)
matplot(mod$alpha,type="l");abline(h=alpha,col=1:2,lty=2)
matplot(mod$theta,type="l");abline(h=theta,col=1:2,lty=2)
matplot(mod$psi,type="l");abline(h=psi,col=1,lty=2)
plot(z[-z_0.idx],mod$z.mean)
apply(mod$beta,2,mean);beta
apply(mod$alpha,2,mean);alpha
apply(mod$theta,2,mean);theta
apply(mod$beta,2,quantile,c(0.025,0.975))
apply(mod$alpha,2,quantile,c(0.025,0.975))
apply(mod$theta,2,quantile,c(0.025,0.975))







# # # # #########################################################################
# # # # #########################################################################
# # # # ##### Simulate dynamic occupancy
# # # # #########################################################################
# # # # #########################################################################

# # # # #####
# # # # #####  Simulate one data source with dynamic occupancy
# # # # #####  Fully heterogeneous model (covariates on occupancy and detection)
# # # # #####

# # # # T <- 1000 #Number of primary periods
# # # # J <- 5 #Number of visits per primary period

# # # # # Covariates for colonization
# # # # X <- matrix(cbind(1,rnorm(T)),T,2)
# # # # beta <- matrix(c(-0.75,0.05),2,1)
# # # # gamma <- pnorm(X%*%beta)
# # # # summary(gamma)

# # # # # Covariates for persistence
# # # # G <- X
# # # # # G <- matrix(cbind(1,rnorm(T)),T,2)
# # # # theta <- matrix(c(-0.25,0.25),2,1)
# # # # phi <- pnorm(G%*%theta)
# # # # summary(phi)

# # # # # Covariates for detectability (site and visit)
# # # # W <- array(NA,dim=c(T,2,J))
# # # # W[,1,] <- 1
# # # # for(i in 1:J){
	# # # # W[,2,i] <- rnorm(nrow(W[,,i]))
# # # # }
# # # # alpha <- matrix(c(-1,-0.075),2,1)
# # # # p <- apply(W,3,function(x) pnorm(x%*%alpha))
# # # # summary(p)

# # # # psi <- 0.5 #Probability of occupancy at t=0

# # # # # Simulate occupancy latent state
# # # # z <- numeric(T+1)
# # # # z[1] <- rbinom(1,1,psi) #Occupancy at t=0
# # # # for(i in 2:(T+1)){
	# # # # if(z[i-1]==1) z[i] <- rbinom(1,1,phi[i-1])
	# # # # if(z[i-1]==0) z[i] <- rbinom(1,1,gamma[i-1])
# # # # }

# # # # # Simulate observed data
# # # # Y <- sapply(1:J,function(x) rbinom(T,1,z[-1]*p[,x]))

# # # # # Fit occupancy model to one data source
# # # # source("dynamic_occupancy/occ.dyn.1.mcmc.R")
# # # # out1 <- occ.dyn.1.mcmc(Y,W,X,G,500,z.true=z)
# # # # matplot(out1$beta,type="l");abline(h=beta,col=1:2,lty=2)
# # # # matplot(out1$alpha,type="l");abline(h=alpha,col=1:2,lty=2)
# # # # matplot(out1$theta,type="l");abline(h=theta,col=1:2,lty=2)
# # # # plot(out1$psi,type="l");abline(h=psi,col=1,lty=2)
# # # # plot(z,out1$z.mean)
# # # # apply(out1$beta,2,mean);beta
# # # # apply(out1$alpha,2,mean);alpha
# # # # apply(out1$theta,2,mean);theta