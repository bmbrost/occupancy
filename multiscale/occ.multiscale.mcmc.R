#
#  Mevin Hooten (20111031), Last Updated: 20131029
#
#  W: covariates for detection probability (p)
#  X: covariates for occupancy probability (psi)
#
#


occ.multiscale.mcmc <- function(y,group.idx,W,U,X,priors,start,tune,n.mcmc,adapt=TRUE){

	###
	###  Libraries and subroutines
	###
	
	expit <- function(logit){
		exp(logit)/(1+exp(logit)) 
	}
	
	logit <- function(expit){
		log(expit/(1-expit))
	}

	get.tune <- function(tune,keep,k,target=0.44){  # adaptive tuning
		a <- min(0.01,1/sqrt(k))
		# a <- min(0.025,1/sqrt(k))
		exp(ifelse(keep<target,log(tune)-a,log(tune)+a))
	}

# y.lik <- function(y,p,phi,log=FALSE){
	# tmp <- (1-phi)*p^y*(1-p)^(1-y)+phi*y
	# if(log) tmp <- log(tmp)
	# tmp
# }

	###
	###  Create variables 
	###
browser()	
	n <- length(unique(group.idx[,1]))
	J <- tapply(group.idx[,2],group.idx[,1],function(x) length(unique(x)))
	K <- tapply(group.idx[,2],group.idx[,1],table)
	K.long <- unlist(K)
	# K <- tapply(group.idx[,2],list(group.idx[,1],group.idx[,2]),length,simplify=FALSE)

# n <- dim(Y)[3]
# J <- apply(Y,3,ncol)
# K <- apply(Y,c(3,2), function(x) sum(!is.na(x)))

	qX <- ncol(X)
	qU <- ncol(U)
	qW <- ncol(W)

# U.mat <- apply(U,2,I)
# W.mat <- apply(U,2,I)

# y <- apply(Y,c(3,2),sum,na.rm=TRUE)
# y0 <- which(y==0)
# n.y0 <- length(y0)


	y0 <- ifelse(y==1,0,1)		
	
	###
	###  Priors
	###
	
	mu.beta <- matrix(priors$mu.beta,qX,1)
	mu.gamma <- matrix(priors$mu.gamma,qX,1) 	
	mu.alpha <- matrix(priors$mu.alpha,qW,1)
	sigma.beta <- priors$sigma.beta
	sigma.gamma <- priors$sigma.gamma
	sigma.alpha <- priors$sigma.alpha

	
	###
	###  Starting values 
	###
# browser()	
	beta <- as.vector(start$beta)
	gamma <- as.vector(start$gamma)
	alpha <- as.vector(start$alpha)
	z <- start$z
	
	# Create an indicator vector to map this z to z.long without computation; same for a.long
	z.long <- c(sapply(1:n,function(x) rep(z[x],J[x])))

	a <- start$a
	a.long <- c(sapply(1:sum(J),function(x) rep(a[x],K.long[x])))

# A <- start$A
	psi <- expit(X%*%beta)  # occupancy probability
	p <- expit(W%*%alpha)  # detection probability
# p <- apply(W,c(3,4),function(x) expit(x%*%alpha))  # detection probability
# Q <- start$Q	
	
	###
	###  Create receptacles for output
	###
	
	beta.save <- matrix(0,n.mcmc,qX)
	gamma.save <- matrix(0,n.mcmc,qU)
	alpha.save <- matrix(0,n.mcmc,qW)
	z.mean <- numeric(n)
	a.mean <- a*0
# A.mean <- matrix(0,n,max(J))
	N.save <- numeric(n.mcmc)
# Q.mean <- matrix(0,n,ncol(Y))

	keep <- list(beta=0,gamma=0,alpha=0)
	keep.tmp <- keep
	Tb <- 50  # frequency of adaptive tuning
	
		
	###
	###  Begin MCMC loop 
	###
	
	for(k in 1:n.mcmc){
		if(k%%1000==0) cat(k," "); flush.console()
	
		###
		### Adaptive tuning
		###
		
		if(adapt==TRUE & k%%Tb==0) {  # Adaptive tuning
			keep.tmp <- lapply(keep.tmp,function(x) x/Tb)
			tune$beta <- get.tune(tune$beta,keep.tmp$beta,k)
			tune$gamma <- get.tune(tune$gamma,keep.tmp$gamma,k)
			tune$alpha <- get.tune(tune$alpha,keep.tmp$alpha,k)
			keep.tmp <- lapply(keep.tmp,function(x) x*0)
	   	} 	

	
	  	###
		###  Sample beta (psi)
	  	###
	
	  	beta.star <- rnorm(qX,beta,tune$beta)
		psi.star <- expit(X%*%beta.star)
			mh.star <- sum(dbinom(z,1,psi.star,log=TRUE))+
			sum(dnorm(beta.star,mu.beta,sigma.beta,log=TRUE))
		mh.0 <- sum(dbinom(z,1,psi,log=TRUE))+
			sum(dnorm(beta,mu.beta,sigma.beta,log=TRUE))
		if(exp(mh.star-mh.0) > runif(1)){
			beta <- beta.star
			psi <- psi.star
			keep$beta <- keep$beta+1
			keep.tmp$beta <- keep.tmp$beta+1
		}


	  	###
		###  Sample gamma (theta)
	  	###
# browser()
		z1 <- z.long==1	
	  	gamma.star <- rnorm(qU,gamma,tune$gamma)
		theta.star <- expit(U%*%gamma.star)
		mh.star <- dbinom(a[z1,],1,theta.star[z1,],log=TRUE)+
			sum(dnorm(gamma.star,mu.gamma,sigma.gamma,log=TRUE))		
		mh.0 <- dbinom(a[z1,],1,theta[z1,],log=TRUE)*z.long+
			sum(dnorm(gamma,mu.gamma,sigma.gamma,log=TRUE))			
		if(exp(mh.star-mh.0) > runif(1)){
			gamma <- gamma.star
			theta <- theta.star
			keep$theta <- keep$theta+1
			keep.tmp$theta <- keep.tmp$theta+1
		}

	
		###
	  	###  Sample alpha (p)
	  	###
# browser()
		a1 <- a.long==1
	  	alpha.star <- rnorm(qW,alpha,tune$alpha)
	  	p.star <- expit(W%*%alpha.star)
	  	mh.star <- dbinom(y[a1],1,p.star[a1])+
	 		sum(dnorm(alpha.star,mu.alpha,sigma.alpha,log=TRUE))
	 	mh.0 <-	dbinom(y[a1],1,p[a1])+
	 		sum(dnorm(alpha,mu.alpha,sigma.alpha,log=TRUE))
	  	# mh.star <- sum(log((dbinom(Y[z1,],1,p.star[z1,])^(1-Q[z1,]))))+
	 		# sum(dnorm(alpha.star,mu.alpha,sigma.alpha,log=TRUE))
	 	# mh.0 <-	sum(log((dbinom(Y[z1,],1,p[z1,])^(1-Q[z1,]))))+
	 		# sum(dnorm(alpha,mu.alpha,sigma.alpha,log=TRUE))
		if(exp(mh.star-mh.0) > runif(1)){
			alpha <- alpha.star
			p <- p.star
			keep$alpha <- keep$alpha+1
			keep.tmp$alpha <- keep.tmp$alpha+1
	  	}

	 
		###
	  	###  Sample z 
	  	###
# browser()
		j.idx <- c(sapply(1:n,function(x) rep(x,J[x])))
		z.tmp <- tapply(a,j.idx,sum)	
		z <- ifelse(z.tmp>0,1,0)
		z.long <- c(sapply(1:n,function(x) rep(z[x],J[x])))

		# tmp <- tapply(dbinom(a,1,theta),idx,prod)
		# p1 <- psi*tmp
		# tmp <- a
		# a0 <- ifelse(a==1,0,1)
		# tmp <- tapply(a0^(1-z.long))
		# z.long
		# p0 <- 
		
		# p1 <- psi*apply(Y^Q*p^Y*((1-p)^(1-Y))^(1-Q),1,prod)
		# p0 <- (1-psi)*apply(Y0^(1-Q)*Y^Q,1,prod)
		# psi.tmp <- p1/(p1+p0)	
		# z <- rbinom(n,1,psi.tmp)


		###
	  	###  Sample a 
	  	###
# browser()
		k.idx <- c(sapply(1:sum(J),function(x) rep(x,K.long[x])))
		theta1 <- c(z.long*theta)*tapply(dbinom(y,1,p),k.idx,prod)
		theta0 <- c(1-z.long*theta)*tapply(y0,k.idx,prod)
		theta.tmp <- theta1/(theta1+theta0)	
		a <- rbinom(n,1,theta.tmp)


		a.long <- c(sapply(1:sum(J),function(x) rep(a[x],K.long[x])))						


		###
	  	###  Save samples 
	  	###
	
	  	beta.save[k,] <- beta
	  	gamma.save[k,] <- alpha
	  	alpha.save[k,] <- alpha
	  	A.mean <- A.mean+A/n.mcmc
	  	z.mean <- z.mean+z/n.mcmc
	  	N.save[k] <- sum(z)
	  	Q.mean <- Q.mean+Q
	
	}
	cat("\n")
	
	###
	###  Write output 
	###

	z.mean <- z.mean/n.mcmc
	Q.mean <- Q.mean/n.mcmc
	
	keep <- lapply(keep,function(x) x/n.mcmc)
	end <- list(beta=beta,gamma=gamma,alpha=alpha,z=z,A=A)  # Q=Q,phi=phi ending values
	
	list(beta=beta.save,gamma=gamma,alpha.save=alpha.save,
		A.mean=A.mean,N.save=N.save,z.mean=z.mean,#Q.mean=Q.mean,
		keep=keep,end=end,Y=Y,X=X,U=U,W=W,priors=priors,start=start,tune=tune,n.mcmc=n.mcmc)
}
