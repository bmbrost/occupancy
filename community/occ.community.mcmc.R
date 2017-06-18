occ.community.mcmc <- function(Y,J,W=NULL,X,priors,start,n.mcmc,adapt=TRUE){

	#
	#  	Brian Brost (11 JUN 2014)
	#	
	#	Fully heterogeneous occupany model using auxilliary variables
	#
	#  	Y: n x J matrix containing full detection history
	# 	W: covariates for detection probability (p)
	#  	X: covariates for occupancy probability (psi)
	#
	#
	
	
	###
	###  Libraries and Subroutines
	###
	
	library(mvtnorm)
	
	truncnormsamp <- function(mu,sig2,low,high,nsamp){
		flow <- pnorm(low,mu,sqrt(sig2)) 
		fhigh <- pnorm(high,mu,sqrt(sig2)) 
		u <- runif(nsamp) 
		tmp <- flow+u*(fhigh-flow)
		x <- qnorm(tmp,mu,sqrt(sig2))
		x
	}
	

	###
	###  Create Variables 
	###
# browser()
	n <- nrow(Y)
	R <- ncol(Y)
	qX <- ncol(X)
	qW <- ifelse(is.null(W),1,ncol(W))
	# J <- apply(Y,1,function(x) sum(!is.na(x)))
	# y <- apply(Y,1,sum,na.rm=TRUE)
	# z <- ifelse(y>0,1,0)	
	# Y.long <- c(Y)
	# W.long <- apply(W,2,I)
	v <- matrix(0,n,R)
	u <- numeric(n)

	# v <- numeric(n)
	# u <- numeric(n*dim(W)[3])	
	
	###
	###  Priors
	###
	
	mu.beta <- matrix(priors$mu.beta,qX,1)	
	mu.alpha <- matrix(priors$mu.alpha,qW,1)	
	sigma.beta <- priors$sigma.beta
	sigma.alpha <- priors$sigma.alpha

	Sigma.beta <- diag(qX)*sigma.beta^2
	Sigma.alpha <- diag(qW)*sigma.alpha^2	
	
	
	###
	###  Starting Values 
	###

	beta <- start$beta
	alpha <- start$alpha
	
	
	###
	###  Create receptacles for output
	###

	beta.save <- array(0,c(n,qX,n.mcmc))
	alpha.save <- array(0,c(n,qW,n.mcmc))	
	z.mean <- matrix(0,n,R)
	N.save <- numeric(R)

	# beta.save <- matrix(0,n.mcmc,qX)
	# alpha.save <- matrix(0,n.mcmc,qW)
	# z.mean <- numeric(n)
	# N.save <- numeric(n.mcmc)

	
	###
	###  Begin MCMC Loop 
	###
	
	for(k in 1:n.mcmc){
		if(k%%1000==0) cat(k," "); flush.console()

browser()		
		
		###
		###  Sample v (auxilliary variable for z) 
		###

		z0 <- z==0
		z1 <- z==1

		for(i in 1:R){
			i <- 1
			z0 <- z[,i]==0
			z1 <- z[,i]==1
			v[z0,i] <- truncnormsamp(matrix(X[i,],,qX)%*%beta[,z0],1,0,Inf,sum(z0))

		}

		v[z1] <- truncnormsamp(matrix(X[z1,],,qX)%*%beta,1,0,Inf,sum(z1))
		v[z0] <- truncnormsamp(matrix(X[z0,],,qX)%*%beta,1,-Inf,0,sum(z0))

		# library(msm)	
		# v[z1] <- rtnorm(sum(z1),(X%*%beta)[z1],lower=0)
		# v[z0] <- rtnorm(sum(z0),(X%*%beta)[z0],upper=0)		


		###
		###  Sample psi
		###

		A <- solve(t(X)%*%X+solve(Sigma.beta))
		b <- t(X)%*%v+solve(Sigma.beta)%*%mu.beta
		beta <- A%*%b+t(chol(A))%*%matrix(rnorm(qX),qX,1)
		# beta <- t(rmvnorm(1,A%*%b,A))
		psi <- pnorm(X%*%beta)
		
				
		###
		###  Sample u (auxilliary variable for p) 
		###
		
		idx1 <- which(Y.long==1&z1)
		idx0 <- which(Y.long==0&z1)	
		u[idx1] <- truncnormsamp((matrix(W.long[idx1,],,qW)%*%alpha),1,0,Inf,length(idx1))			
		u[idx0] <- truncnormsamp((matrix(W.long[idx0,],,qW)%*%alpha),1,-Inf,0,length(idx0))

		
		###
		###  Sample p 
		###

		u.tmp <- u[z1]
		W.tmp <- W.long[z1,]		
		A <- solve(t(W.tmp)%*%W.tmp+solve(Sigma.alpha))
		b <- u.tmp%*%W.tmp+t(mu.alpha)%*%solve(Sigma.alpha)
		alpha <- A%*%t(b)+t(chol(A))%*%matrix(rnorm(qW),qW,1)		
		p <- matrix(apply(W.long,1,function(x) pnorm(x%*%alpha)),,max(J))

		
		###
		###  Sample z 
		###
		
		p.tmp <- apply(1-p,1,prod)
		num.tmp <- psi*p.tmp
		psi.tmp <- num.tmp/(num.tmp+(1-psi))
		z[y==0] <- rbinom(sum(y==0),1,psi.tmp[y==0])
		
		
		###
		###  Save samples 
		###

		beta.save[k,] <- beta
		alpha.save[k,] <- alpha
		z.mean <- z.mean+z/n.mcmc
		N.save[k] <- sum(z)
	
		}
	cat("\n")
	
	###
	###  Write output 
	###
	
	list(beta=beta.save,alpha=alpha.save,N=N.save,z.mean=z.mean,
	Y=Y,W=W,X=X,priors=priors,start=start,n.mcmc=n.mcmc)
}





# occ.probit.1.mcmc <- function(Y,W,X,n.mcmc){

	# #
	# #  	Brian Brost (11 JUN 2014)
	# #	
	# #	Fully heterogeneous occupany model using auxilliary variables
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
	# n <- nrow(Y)
	# qX <- ncol(X)
	# qW <- ncol(W)
	# beta.save <- matrix(0,n.mcmc,qX)
	# alpha.save <- matrix(0,n.mcmc,qW)
	# z.mean <- numeric(n)
	# N.save <- numeric(n.mcmc)
	# v <- numeric(n)
	# # u <- matrix(0,n,dim(W)[3])	
	# y <- apply(Y,1,sum,na.rm=TRUE)
	# u <- numeric(n*dim(W)[3])	
	# J <- apply(Y,1,function(x) sum(!is.na(x)))
	
	# ####
	# ####  Priors and Starting Values 
	# ####
	
	# beta.mn <- matrix(0,qX,1)
	# alpha.mn <- matrix(0,qW,1)
	# beta.var <- diag(qX)*10
	# alpha.var <- diag(qW)*10
	# z <- ifelse(y>0,1,0)

# # browser()
	# beta <- as.vector(glm(z ~ 0+X,family=binomial())$coefficients)
	# alpha <- as.vector(glm(cbind(y[y>0],J[y>0]-y[y>0]) ~ 
		# 0+W[y>0,,1],family=binomial())$coefficients)

# Y <- c(Y)
# W <- apply(W,2,I)

	
	# ####
	# ####  Begin MCMC Loop 
	# ####
	
	# for(k in 1:n.mcmc){
		# if(k%%100==0) cat(k," "); flush.console()
		
		# ####
		# ####  Sample v (auxilliary variable for z) 
		# ####
# # browser()		
		# z0 <- z==0
		# z1 <- z==1
		# v[z1] <- truncnormsamp(matrix(X[z1,],,qX)%*%beta,1,0,Inf,sum(z1))
		# v[z0] <- truncnormsamp(matrix(X[z0,],,qX)%*%beta,1,-Inf,0,sum(z0))

		# # v[z1] <- truncnormsamp(X[z1,]%*%beta,1,0,Inf,sum(z1))
		# # v[z0] <- truncnormsamp(X[z0,]%*%beta,1,-Inf,0,sum(z0))

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
		
		# idx1 <- which(Y==1&z1)
		# idx0 <- which(Y==0&z1)	
		# u[idx1] <- truncnormsamp((matrix(W[idx1,],,qW)%*%alpha),1,0,Inf,length(idx1))			
		# u[idx0] <- truncnormsamp((matrix(W[idx0,],,qW)%*%alpha),1,-Inf,0,length(idx0))

		# # u[idx1] <- truncnormsamp((W[idx1,]%*%alpha),1,0,Inf,length(idx1))			
		# # u[idx0] <- truncnormsamp((W[idx0,]%*%alpha),1,-Inf,0,length(idx0))
		
		# # for(i in 1:ncol(Y)){
			# # # i <- 1
			# # y1 <- Y[,i]==1&z1		
			# # y0 <- Y[,i]==0&z1		
			# # u[y1,i] <- truncnormsamp((W[y1,,i]%*%alpha),1,0,Inf,sum(y1))
			# # u[y0,i] <- truncnormsamp((W[y0,,i]%*%alpha),1,-Inf,0,sum(y0))
		# # }
	
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
	
		# u.tmp <- u[z1]
		# W.tmp <- W[z1,]		
		# A <- solve(t(W.tmp)%*%W.tmp+solve(alpha.var))
		# b <- u.tmp%*%W.tmp+t(alpha.mn)%*%solve(alpha.var)
		# alpha <- A%*%t(b)+t(chol(A))%*%matrix(rnorm(qW),qW,1)		
		# p <- matrix(apply(W,1,function(x) pnorm(x%*%alpha)),,max(J))

		# ####
		# ####  Sample z 
		# ####
		
		# p.tmp <- apply(1-p,1,prod)
		# # p.tmp <- 1-p
		# num.tmp <- psi*p.tmp#^J 
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
