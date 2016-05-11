occ.probit.2.mcmc <- function(Y1,Y2,W1,W2,X,n.mcmc){

	#
	#  	Brian Brost (11 JUN 2014)
	#
	#  	Y: n x J matrix containing full detection history
	# 	W: covariates for detection probability (p)
	#  	X: covariates for occupancy probability (psi)
	#
	#
	
	####
	####  Libraries and Subroutines
	####
	
	library(mvtnorm)
	
	expit <- function(logit){
		 exp(logit)/(1+exp(logit)) 
	}

	logit <- function(expit){
		log(expit/(1-expit))
	}
	
	truncnormsamp <- function(mu,sig2,low,high,nsamp){
		flow <- pnorm(low,mu,sqrt(sig2)) 
		fhigh <- pnorm(high,mu,sqrt(sig2)) 
		u <- runif(nsamp) 
		tmp <- flow+u*(fhigh-flow)
		x <- qnorm(tmp,mu,sqrt(sig2))
		x
	}
	
	####
	####  Create Variables 
	####
# browser()
	n <- nrow(Y1)
	qX <- ncol(X)
	qW1 <- ncol(W1)
	qW2 <- ncol(W2)
	J1 <- apply(Y1,1,function(x) sum(!is.na(x)))
	J2 <- apply(Y2,1,function(x) sum(!is.na(x)))

	beta.save <- matrix(0,n.mcmc,qX)
	alpha1.save <- matrix(0,n.mcmc,qW1)
	alpha2.save <- matrix(0,n.mcmc,qW2)
	z.mean <- numeric(n)
	N.save <- numeric(n.mcmc)
	v <- numeric(n)
	u1 <- numeric(n*dim(W1)[3]) # matrix(0,n,dim(W1)[3])	
	u2 <- numeric(n*dim(W2)[3]) # matrix(0,n,dim(W2)[3])	
	y1 <- apply(Y1,1,sum,na.rm=TRUE)
	y2 <- apply(Y2,1,sum,na.rm=TRUE)
	
	####
	####  Priors and Starting Values 
	####
	
	beta.mn <- matrix(0,qX,1)
	alpha1.mn <- matrix(0,qW1,1)
	alpha2.mn <- matrix(0,qW2,1)
	beta.var <- diag(qX)*10
	alpha1.var <- diag(qW1)*10
	alpha2.var <- diag(qW2)*10
	z <- ifelse(y1>0|y2>0,1,0)

# browser()
	beta <- as.vector(glm(z ~ 0+X,family=binomial())$coefficients)
	alpha1 <- as.vector(glm(cbind(y1[y1>0],J1[y1>0]-y1[y1>0]) ~ 
		0+W1[y1>0,,1],family=binomial())$coefficients)
	alpha2 <- as.vector(glm(cbind(y2[y2>0],J2[y2>0]-y2[y2>0]) ~ 
		0+W2[y2>0,,1],family=binomial())$coefficients)
# browser()
	Y1 <- c(Y1)
	Y2 <- c(Y2)
	W1 <- apply(W1,2,I)
	W2 <- apply(W2,2,I)
	
		
	####
	####  Begin MCMC Loop 
	####
	
	for(k in 1:n.mcmc){
		if(k%%100==0) cat(k,""); flush.console()
		
		####
		####  Sample v (auxilliary variable for z) 
		####
# browser()		
		z0 <- z%in%0 #Excludes NA values
		z1 <- z%in%1 #Excludes NA values		
		# z0 <- z==0 #Includes NA values
		# z1 <- z==1 #Includes NA values
		v[z1] <- truncnormsamp(matrix(X[z1,],,qX)%*%beta,1,0,Inf,sum(z1))
		v[z0] <- truncnormsamp(matrix(X[z0,],,qX)%*%beta,1,-Inf,0,sum(z0))

		# library(msm)	
		# v[z1] <- rtnorm(sum(z1),(X%*%beta)[z1],lower=0)
		# v[z0] <- rtnorm(sum(z0),(X%*%beta)[z0],upper=0)		


		####
		####  Sample psi (beta) 
		####
# browser()				
		A <- solve(t(X)%*%X+solve(beta.var))
		b <- t(X)%*%v+solve(beta.var)%*%beta.mn
		beta <- A%*%b+t(chol(A))%*%matrix(rnorm(qX),qX,1)
		# beta <- t(rmvnorm(1,A%*%b,A))
		psi <- pnorm(X%*%beta)
		
				
		####
		####  Sample u1 (auxilliary variable for p1, first data source) 
		####
# browser()		
		idx1 <- which(Y1==1&z1)
		idx0 <- which(Y1==0&z1)	
		u1[idx1] <- truncnormsamp((matrix(W1[idx1,],,qW1)%*%alpha1),1,0,Inf,length(idx1))			
		u1[idx0] <- truncnormsamp((matrix(W1[idx0,],,qW1)%*%alpha1),1,-Inf,0,length(idx0))


		####
		####  Sample p1 (alpha1) 
		####
# browser()	
		u.tmp <- u1[z1]
		W.tmp <- W1[z1,]

		# Excluce records where z=1 but is.na(W)
		idx <- which(!is.na(W.tmp))
		u.tmp <- u.tmp[idx]		
		W.tmp <- W.tmp[idx]
		
		A <- solve(t(W.tmp)%*%W.tmp+solve(alpha1.var))
		b <- u.tmp%*%W.tmp+t(alpha1.mn)%*%solve(alpha1.var)
		alpha1 <- A%*%t(b)+t(chol(A))%*%matrix(rnorm(qW1),qW1,1)		
		p1 <- matrix(apply(W1,1,function(x) pnorm(x%*%alpha1)),,max(J1))


		####
		####  Sample u2 (auxilliary variable for p2, first data source) 
		####
# browser()		
		idx1 <- which(Y2==1&z1)
		idx0 <- which(Y2==0&z1)	
		u2[idx1] <- truncnormsamp((matrix(W2[idx1,],,qW2)%*%alpha2),1,0,Inf,length(idx1))			
		u2[idx0] <- truncnormsamp((matrix(W2[idx0,],,qW2)%*%alpha2),1,-Inf,0,length(idx0))


		####
		####  Sample p1 (alpha1) 
		####
# browser()	
		u.tmp <- u2[z1]
		W.tmp <- W2[z1,]		
		
		# Excluce records where z=1 but is.na(W)
		idx <- which(!is.na(W.tmp))
		u.tmp <- u.tmp[idx]		
		W.tmp <- W.tmp[idx]
				
		A <- solve(t(W.tmp)%*%W.tmp+solve(alpha2.var))
		b <- u.tmp%*%W.tmp+t(alpha2.mn)%*%solve(alpha2.var)
		alpha2 <- A%*%t(b)+t(chol(A))%*%matrix(rnorm(qW2),qW2,1)		
		p2 <- matrix(apply(W2,1,function(x) pnorm(x%*%alpha2)),,max(J2))
		
		
		####
		####  Sample z 
		####
# browser()		
		# p1.tmp <- matrix(p1^Y1*(1-p1)^(1-Y1),,max(J1))
		# p1.tmp <- apply(p1.tmp,1,prod)

		idx <- which(y1==0&y2==0)		
		p1.tmp <- apply(1-p1,1,prod,na.rm=TRUE)
		p2.tmp <- apply(1-p2,1,prod,na.rm=TRUE)		
		num.tmp <- psi*p1.tmp*p2.tmp 
		psi.tmp <- num.tmp/(num.tmp+(1-psi))
		z[idx] <- rbinom(length(idx),1,psi.tmp[idx])

		# p <- p1
		# p.tmp <- apply(1-p,1,prod)
		# # p.tmp <- 1-p
		# num.tmp <- psi*p.tmp
		# psi.tmp <- num.tmp/(num.tmp+(1-psi))
		# z[y==0] <- rbinom(sum(y==0),1,psi.tmp[y==0])
		
		####
		####  Save Samples 
		####
		
		beta.save[k,] <- beta
		alpha1.save[k,] <- alpha1
		alpha2.save[k,] <- alpha2
		z.mean <- z.mean+z/n.mcmc
		N.save[k] <- sum(z)
	
		}
	cat("\n")
	
	####
	####  Write Output 
	####
	
	list(beta=beta.save,alpha1=alpha1.save,alpha2=alpha2.save,N.save=N.save,z.mean=z.mean,n.mcmc=n.mcmc)
	
}








# occ.2.mcmc <- function(Y1,Y2,J1,J2,W1,W2,X,n.mcmc){

	# #
	# #  	Brian Brost (11 JUN 2014)
	# #
	# #  	Y: n x J matrix containing full detection history
	# # 	W: covariates for detection probability (p)
	# #  	X: covariates for occupancy probability (psi)
	# #
	# #
	
	# ####
	# ####  Libraries and Subroutines
	# ####
	
	# library(mvtnorm)
	
	# expit <- function(logit){
		 # exp(logit)/(1+exp(logit)) 
	# }

	# logit <- function(expit){
		# log(expit/(1-expit))
	# }
	
	# truncnormsamp <- function(mu,sig2,low,high,nsamp){
		# flow <- pnorm(low,mu,sqrt(sig2)) 
		# fhigh <- pnorm(high,mu,sqrt(sig2)) 
		# u <- runif(nsamp) 
		# tmp <- flow+u*(fhigh-flow)
		# x <- qnorm(tmp,mu,sqrt(sig2))
		# x
	# }
	
	# ####
	# ####  Create Variables 
	# ####
# # browser()
	# n <- nrow(Y1)
	# qX <- ncol(X)
	# qW1 <- ncol(W1)
	# qW2 <- ncol(W2)
	# beta.save <- matrix(0,n.mcmc,qX)
	# alpha1.save <- matrix(0,n.mcmc,qW1)
	# alpha2.save <- matrix(0,n.mcmc,qW2)
	# z.mean <- numeric(n)
	# N.save <- numeric(n.mcmc)
	# v <- numeric(n)
	# u1 <- numeric(n*dim(W1)[3]) # matrix(0,n,dim(W1)[3])	
	# u2 <- numeric(n*dim(W2)[3]) # matrix(0,n,dim(W2)[3])	
	# y1 <- apply(Y1,1,sum)
	# y2 <- apply(Y2,1,sum)
	
	# ####
	# ####  Priors and Starting Values 
	# ####
	
	# beta.mn <- matrix(0,qX,1)
	# alpha1.mn <- matrix(0,qW1,1)
	# alpha2.mn <- matrix(0,qW2,1)
	# beta.var <- diag(qX)*10
	# alpha1.var <- diag(qW1)*10
	# alpha2.var <- diag(qW2)*10
	# z <- ifelse(y1>0|y2>0,1,0)

# # browser()
	# beta <- as.vector(glm(z ~ 0+X,family=binomial())$coefficients)
	# alpha1 <- as.vector(glm(cbind(y1[y1>0],J1[y1>0]-y1[y1>0]) ~ 
		# 0+W1[y1>0,,1],family=binomial())$coefficients)
	# alpha2 <- as.vector(glm(cbind(y2[y2>0],J2[y2>0]-y2[y2>0]) ~ 
		# 0+W2[y2>0,,1],family=binomial())$coefficients)

	# Y1 <- c(Y1)
	# Y2 <- c(Y2)
	# W1 <- apply(W1,2,I)
	# W2 <- apply(W2,2,I)

	
	# ####
	# ####  Begin MCMC Loop 
	# ####
	
	# for(k in 1:n.mcmc){
		# if(k%%100==0) cat(k," "); flush.console()
		
		# ####
		# ####  Sample v (auxilliary variable for z) 
		# ####
# # browser()		
		# z0 <- z%in%0 #Excludes NA values
		# z1 <- z%in%1 #Excludes NA values		
		# # z0 <- z==0 #Includes NA values
		# # z1 <- z==1 #Includes NA values
		# v[z1] <- truncnormsamp(matrix(X[z1,],,qX)%*%beta,1,0,Inf,sum(z1))
		# v[z0] <- truncnormsamp(matrix(X[z0,],,qX)%*%beta,1,-Inf,0,sum(z0))

		# # library(msm)	
		# # v[z1] <- rtnorm(sum(z1),(X%*%beta)[z1],lower=0)
		# # v[z0] <- rtnorm(sum(z0),(X%*%beta)[z0],upper=0)		


		# ####
		# ####  Sample psi (beta) 
		# ####
# # browser()				
		# A <- solve(t(X)%*%X+solve(beta.var))
		# b <- t(X)%*%v+solve(beta.var)%*%beta.mn
		# beta <- A%*%b+t(chol(A))%*%matrix(rnorm(qX),qX,1)
		# # beta <- t(rmvnorm(1,A%*%b,A))
		# psi <- pnorm(X%*%beta)
		
				
		# ####
		# ####  Sample u1 (auxilliary variable for p1, first data source) 
		# ####
# # browser()		
		# idx1 <- which(Y1==1&z1)
		# idx0 <- which(Y1==0&z1)	
		# u1[idx1] <- truncnormsamp((matrix(W1[idx1,],,qW1)%*%alpha1),1,0,Inf,length(idx1))			
		# u1[idx0] <- truncnormsamp((matrix(W1[idx0,],,qW1)%*%alpha1),1,-Inf,0,length(idx0))


		# ####
		# ####  Sample p1 (alpha1) 
		# ####
# # browser()	
		# u.tmp <- u1[z1]
		# W.tmp <- W1[z1,]		
		# A <- solve(t(W.tmp)%*%W.tmp+solve(alpha1.var))
		# b <- u.tmp%*%W.tmp+t(alpha1.mn)%*%solve(alpha1.var)
		# alpha1 <- A%*%t(b)+t(chol(A))%*%matrix(rnorm(qW1),qW1,1)		
		# p1 <- matrix(apply(W1,1,function(x) pnorm(x%*%alpha1)),,max(J1))


		# ####
		# ####  Sample u2 (auxilliary variable for p2, first data source) 
		# ####
# # browser()		
		# idx1 <- which(Y2==1&z1)
		# idx0 <- which(Y2==0&z1)	
		# u2[idx1] <- truncnormsamp((matrix(W2[idx1,],,qW2)%*%alpha2),1,0,Inf,length(idx1))			
		# u2[idx0] <- truncnormsamp((matrix(W2[idx0,],,qW2)%*%alpha2),1,-Inf,0,length(idx0))


		# ####
		# ####  Sample p1 (alpha1) 
		# ####
# # browser()	
		# u.tmp <- u2[z1]
		# W.tmp <- W2[z1,]		
		# A <- solve(t(W.tmp)%*%W.tmp+solve(alpha2.var))
		# b <- u.tmp%*%W.tmp+t(alpha2.mn)%*%solve(alpha2.var)
		# alpha2 <- A%*%t(b)+t(chol(A))%*%matrix(rnorm(qW2),qW2,1)		
		# p2 <- matrix(apply(W2,1,function(x) pnorm(x%*%alpha2)),,max(J2))

		
		# ####
		# ####  Sample z 
		# ####
# # browser()		
		# # p1.tmp <- matrix(p1^Y1*(1-p1)^(1-Y1),,max(J1))
		# # p1.tmp <- apply(p1.tmp,1,prod)

		# idx <- which(y1==0&y2==0)		
		# p1.tmp <- apply(1-p1,1,prod)
		# p2.tmp <- apply(1-p2,1,prod)		
		# num.tmp <- psi*p1.tmp*p2.tmp 
		# psi.tmp <- num.tmp/(num.tmp+(1-psi))
		# z[idx] <- rbinom(length(idx),1,psi.tmp[idx])

		# # p <- p1
		# # p.tmp <- apply(1-p,1,prod)
		# # # p.tmp <- 1-p
		# # num.tmp <- psi*p.tmp
		# # psi.tmp <- num.tmp/(num.tmp+(1-psi))
		# # z[y==0] <- rbinom(sum(y==0),1,psi.tmp[y==0])
		
		# ####
		# ####  Save Samples 
		# ####
		
		# beta.save[k,] <- beta
		# alpha1.save[k,] <- alpha1
		# alpha2.save[k,] <- alpha2
		# z.mean <- z.mean+z/n.mcmc
		# N.save[k] <- sum(z)
	
		# }
	# cat("\n")
	
	# ####
	# ####  Write Output 
	# ####
	
	# list(beta=beta.save,alpha1=alpha1.save,alpha2=alpha2.save,N.save=N.save,z.mean=z.mean,n.mcmc=n.mcmc)
	
# }






# occ.mcmc <- function(Y,J,W,X,n.mcmc,alpha.tune=NULL,beta.tune=NULL){

	# #
	# #  Mevin Hooten (20111031), Last Updated: 20131029
	# #
	# #  W: covariates for detection probability (p)
	# #  X: covariates for occupancy probability (psi)
	# #
	# #
	
	# ####
	# ####  Libraries and Subroutines
	# ####
	
	# library(mvtnorm)
	
	# expit <- function(logit){
		 # exp(logit)/(1+exp(logit)) 
	# }

	# logit <- function(expit){
		# log(expit/(1-expit))
	# }
	
	# truncnormsamp <- function(mu,sig2,low,high,nsamp){
		# flow <- pnorm(low,mu,sqrt(sig2)) 
		# fhigh <- pnorm(high,mu,sqrt(sig2)) 
		# u <- runif(nsamp) 
		# tmp <- flow+u*(fhigh-flow)
		# x <- qnorm(tmp,mu,sqrt(sig2))
		# x
	# }
	
	# ####
	# ####  Create Variables 
	# ####
# # browser()
	# n <- nrow(Y)
	# qX <- ncol(X)
	# qW <- ncol(W)
	# beta.save <- matrix(0,n.mcmc,qX)
	# alpha.save <- matrix(0,n.mcmc,qW)
	# z.mean <- numeric(n)
	# N.save <- numeric(n.mcmc)
	# v <- numeric(n)
	# u <- matrix(0,n,dim(W)[3])	
	# y <- apply(Y,1,sum)
	
	
	# ####
	# ####  Priors and Starting Values 
	# ####
	
	# beta.mn <- matrix(0,qX,1)
	# alpha.mn <- matrix(0,qW,1)
	# beta.var <- diag(qX)*10
	# alpha.var <- diag(qW)*10
# # alpha.var <- rep(10,qW)
	# z <- ifelse(y>0,1,0)

# # browser()
	# beta <- as.vector(glm(z ~ 0+X,family=binomial())$coefficients)
	# if(is.null(beta.tune)){beta.tune <- 0.1*abs(beta)}
	
	# alpha <- as.vector(glm(cbind(y[y>0],J[y>0]-y[y>0]) ~ 
		# 0+W[y>0,,1],family=binomial())$coefficients)
	# if(is.null(alpha.tune)){alpha.tune <- 0.1*abs(alpha)}
	
	# ####
	# ####  Begin MCMC Loop 
	# ####
	
	# for(k in 1:n.mcmc){
		# if(k%%1000==0) cat(k," "); flush.console()
		
		# ####
		# ####  Sample v (auxilliary variable for z) 
		# ####
# # browser()		
		# z0 <- z==0
		# z1 <- z==1
		# v[z1] <- truncnormsamp(X[z1,]%*%beta,1,0,Inf,sum(z1))
		# v[z0] <- truncnormsamp(X[z0,]%*%beta,1,-Inf,0,sum(z0))

		# # library(msm)	
		# # v[z1] <- rtnorm(sum(z1),(X%*%beta)[z1],lower=0)
		# # v[z0] <- rtnorm(sum(z0),(X%*%beta)[z0],upper=0)		


		# ####
		# ####  Sample psi (beta) 
		# ####
# # browser()				
		# A <- solve(t(X)%*%X+solve(beta.var))
		# b <- t(X)%*%v+solve(beta.var)%*%beta.mn
		# beta <- A%*%b+t(chol(A))%*%matrix(rnorm(qX),qX,1)
		# # beta <- t(rmvnorm(1,A%*%b,A))
		# psi <- pnorm(X%*%beta)
		
				
		# ####
		# ####  Sample u (auxilliary variable for p) 
		# ####
# # browser()		
		# # y0 <- y==0&z==1
		# # y1 <- y>0&z==1
		# # u[y1] <- truncnormsamp(W[y1,,1]%*%alpha,1,0,Inf,sum(y1))
		# # u[y0] <- truncnormsamp(W[y0,,1]%*%alpha,1,-Inf,0,sum(y0))		
		
		# for(i in 1:ncol(Y)){
			# # i <- 1
			# y1 <- Y[,i]==1&z1		
			# y0 <- Y[,i]==0&z1		
			# u[y1,i] <- truncnormsamp((W[y1,,i]%*%alpha),1,0,Inf,sum(y1))
			# u[y0,i] <- truncnormsamp((W[y0,,i]%*%alpha),1,-Inf,0,sum(y0))
		# }
	
		# ####
		# ####  Sample p (alpha) 
		# ####
# # browser()	

		# # u.tmp <- c(u[z1,])
		# # W.tmp <- apply(W[z1,,],2,I)		
		# # A <- solve(t(W.tmp)%*%W.tmp+solve(alpha.var))
		# # b <- t(W.tmp)%*%u.tmp+solve(alpha.var)%*%alpha.mn
		# # alpha <- A%*%b+t(chol(A))%*%matrix(rnorm(qW),qW,1)
		# # # beta <- t(rmvnorm(1,A%*%b,A))
		# # p <- pnorm(W[,,1]%*%alpha)
	
		# u.tmp <- c(u[z1,])
		# W.tmp <- apply(W[z1,,],2,I)		
		# A <- solve(t(W.tmp)%*%W.tmp+solve(alpha.var))
		# b <- u.tmp%*%W.tmp+t(alpha.mn)%*%solve(alpha.var)
		# alpha <- A%*%t(b)+t(chol(A))%*%matrix(rnorm(qW),qW,1)		
		# p <- apply(W,3,function(x) pnorm(x%*%alpha))
	
	
		# ####
		# ####  Sample z 
		# ####
		
		# p.tmp <- apply(1-p,1,prod)
		# # p.tmp <- 1-p
		# num.tmp <- psi*p.tmp^J 
		# psi.tmp <- num.tmp/(num.tmp+(1-psi))
		# z[y==0] <- rbinom(sum(y==0),1,psi.tmp[y==0])
		
		# ####
		# ####  Save Samples 
		# ####
		
		# beta.save[k,] <- beta
		# alpha.save[k,] <- alpha
		# z.mean <- z.mean+z/n.mcmc
		# N.save[k] <- sum(z)
	
		# }
	# cat("\n")
	
	# ####
	# ####  Write Output 
	# ####
	
	# list(beta.save=beta.save,alpha.save=alpha.save,N.save=N.save,z.mean=z.mean,n.mcmc=n.mcmc)
	
# }





# occ.mcmc <- function(y,J,W,X,n.mcmc,alpha.tune=NULL,beta.tune=NULL){

	# #
	# #  Mevin Hooten (20111031), Last Updated: 20131029
	# #
	# #  W: covariates for detection probability (p)
	# #  X: covariates for occupancy probability (psi)
	# #
	# #
	
	# ####
	# ####  Libraries and Subroutines
	# ####
	
	# library(mvtnorm)
	
	# expit <- function(logit){
		 # exp(logit)/(1+exp(logit)) 
	# }

	# logit <- function(expit){
		# log(expit/(1-expit))
	# }
	
	# truncnormsamp <- function(mu,sig2,low,high,nsamp){
		# flow <- pnorm(low,mu,sqrt(sig2)) 
		# fhigh <- pnorm(high,mu,sqrt(sig2)) 
		# u <- runif(nsamp) 
		# tmp <- flow+u*(fhigh-flow)
		# x <- qnorm(tmp,mu,sqrt(sig2))
		# x
	# }
	
	# ####
	# ####  Create Variables 
	# ####
	
	# n <- nrow(y)
	# qX <- ncol(X)
	# qW <- ncol(W)
	# beta.save <- matrix(0,n.mcmc,qX)
	# alpha.save <- matrix(0,n.mcmc,qW)
	# z.mean <- numeric(n)
	# N.save <- numeric(n.mcmc)
	# v <- numeric(n)
	# u <- matrix(0,n,J)	
	# y.idot <- apply(y,1,sum)
	
	
	# ####
	# ####  Priors and Starting Values 
	# ####
	
	# beta.mn <- matrix(0,qX,1)
	# alpha.mn <- rep(0,qW)
	# beta.var <- diag(qX)*10
	# alpha.var <- diag(qW)*10
# alpha.var <- rep(10,qW)
	# z <- ifelse(y.idot>0,1,0)
	
# # browser()
	# beta <- as.vector(glm(z ~ 0+X,family=binomial())$coefficients)
	# if(is.null(beta.tune)){beta.tune <- 0.1*abs(beta)}
	
	# alpha <- as.vector(glm(cbind(y.idot[y.idot>1],J[y.idot>1]-y.idot[y.idot>1]) ~ 
		# 0+W[y.idot>1,],family=binomial())$coefficients)
	# if(is.null(alpha.tune)){alpha.tune <- 0.1*abs(alpha)}
	
	# ####
	# ####  Begin MCMC Loop 
	# ####
	
	# for(k in 1:n.mcmc){
		# if(k%%1000==0) cat(k," "); flush.console()
		
		# ####
		# ####  Sample v (auxilliary variable for z) 
		# ####
# # browser()		
		# z0 <- z==0
		# z1 <- z==1
		# v[z1] <- truncnormsamp((X%*%beta)[z1],1,0,Inf,sum(z1))
		# v[z0] <- truncnormsamp((X%*%beta)[z0],1,-Inf,0,sum(z0))

		# # library(msm)	
		# # v[z1] <- rtnorm(sum(z1),(X%*%beta)[z1],lower=0)
		# # v[z0] <- rtnorm(sum(z0),(X%*%beta)[z0],upper=0)		


		# ####
		# ####  Sample psi (beta) 
		# ####
# # browser()				
		# A <- solve(t(X)%*%X+solve(beta.var))
		# b <- t(X)%*%v+solve(beta.var)%*%beta.mn
		# beta <- A%*%b+t(chol(A))%*%matrix(rnorm(qX),qX,1)
		# # beta <- t(rmvnorm(1,A%*%b,A))
		# psi <- pnorm(X%*%beta)
		
		# # beta.star <- rnorm(qX,beta,beta.tune)
		# # mh1 <- sum(dbinom(z,1,expit(X%*%beta.star),log=TRUE)) +
			# # sum(dnorm(beta.star,beta.mn,sqrt(beta.var),log=TRUE))
		# # mh2 <- sum(dbinom(z,1,expit(X%*%beta),log=TRUE)) +
			# # sum(dnorm(beta,beta.mn,sqrt(beta.var),log=TRUE))
		# # if(exp(mh1-mh2) > runif(1)){
			# # beta <- beta.star
		# # }
		# # psi <- expit(X%*%beta)
		

		# ####
		# ####  Sample p (alpha) 
		# ####
	
		
		
	
	
		# alpha.star <- rnorm(qW,alpha,alpha.tune)
		# mh1 <- sum(dbinom(y.idot[z==1],J[z==1],expit(W[z==1,]%*%alpha.star),log=TRUE)) +
			# sum(dnorm(alpha.star,alpha.mn,sqrt(alpha.var),log=TRUE))
		# mh2 <- sum(dbinom(y.idot[z==1],J[z==1],expit(W[z==1,]%*%alpha),log=TRUE)) +
			# sum(dnorm(alpha,alpha.mn,sqrt(alpha.var),log=TRUE))
		# if(exp(mh1-mh2) > runif(1)){
			# alpha <- alpha.star
		# }
		# p <- expit(W%*%alpha)
		 
		# ####
		# ####  Sample z 
		# ####
		
		# num.tmp <- psi*(1-p)^J 
		# psi.tmp <- num.tmp/(num.tmp+(1-psi))
		# z[y.idot==0] <- rbinom(sum(y.idot==0),1,psi.tmp[y.idot==0])
		
		# ####
		# ####  Save Samples 
		# ####
		
		# beta.save[k,] <- beta
		# alpha.save[k,] <- alpha
		# z.mean <- z.mean+z/n.mcmc
		# N.save[k] <- sum(z)
	
		# }
	# cat("\n")
	
	# ####
	# ####  Write Output 
	# ####
	
	# list(beta.save=beta.save,alpha.save=alpha.save,N.save=N.save,z.mean=z.mean,n.mcmc=n.mcmc)
	
# }
