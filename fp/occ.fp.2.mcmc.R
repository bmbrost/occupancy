occ.fp.2.mcmc <- function(Y,W,X,priors,start,tune,n.mcmc,adapt=TRUE){

#
#  Mevin Hooten (20111031), Last Updated: 20131029
#
#  W: covariates for detection probability (p)
#  X: covariates for occupancy probability (psi)
#
#

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
		# a <- min(0.01,1/sqrt(k))
		a <- min(0.025,1/sqrt(k))
		exp(ifelse(keep<target,log(tune)-a,log(tune)+a))
	}

	y.loglik <- function(y,p,phi){
		(1-phi)*p^y*(1-p)^(1-y)+phi*y
	}
	
	###
	###  Create variables 
	###
# browser()	
	n <- nrow(Y)
	qX <- ncol(X)
	qW <- ncol(W)
	J <- apply(Y,1,function(x) sum(!is.na(x)))
	y <- apply(Y,1,sum,na.rm=TRUE)
	y0 <- which(y==0)
	n.y0 <- length(y0)

	
	###
	###  Priors
	###
	
	mu.beta <- matrix(priors$mu.beta,qX,1)
	mu.alpha <- matrix(priors$mu.alpha,qW,1)
	sigma.beta <- priors$sigma.beta
	sigma.alpha <- priors$sigma.alpha

	
	###
	###  Starting values 
	###
	
	beta <- as.vector(start$beta)
	alpha <- as.vector(start$alpha)
	z <- start$z
	psi <- expit(X%*%beta)
	p <- apply(W,3,function(x) expit(x%*%alpha))
	
	
	###
	###  Create receptacles for output
	###
	
	beta.save <- matrix(0,n.mcmc,qX)
	alpha.save <- matrix(0,n.mcmc,qW)
	z.mean <- numeric(n)
	N.save <- numeric(n.mcmc)

	keep <- list(beta=0,alpha=0)
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
		mh.0 <- sum(dbinom(z,1,psi,log=TRUE))	+sum(dnorm(beta,mu.beta,sigma.beta,log=TRUE))
		if(exp(mh.star-mh.0) > runif(1)){
			beta <- beta.star
			psi <- psi.star
			keep$beta <- keep$beta+1
			keep.tmp$beta <- keep.tmp$beta+1
		}

	
		###
	  	###  Sample p (alpha)
	  	###
# browser()
		z1 <- z==1	
	  	alpha.star <- rnorm(qW,alpha,tune$alpha)
	  	p.star <- apply(W,3,function(x) expit(x%*%alpha.star))	 	
	 	mh.star <- sum(log((1-phi)*dbinom(Y[z1,],1,p.star[z1,])+phi*Y[z1,]))
	 		+sum(dnorm(alpha.star,mu.alpha,sigma.alpha,log=TRUE))
	 		# sum(log((1-phi)*p.star[z1,]^Y[z1,]*(1-p.star[z1,])^(1-Y[z1,])
	 			# +phi*Y[z1,]))
	 		# sum(dbinom(Y[z1,],1,p.star[z1,],log=TRUE))
	 	mh.0 <- sum(log((1-phi)*dbinom(Y[z1,],1,p[z1,])+phi*Y[z1,]))
	 		+sum(dnorm(alpha,mu.alpha,sigma.alpha,log=TRUE))
	 		# sum(log((1-phi)*p[z1,]^Y[z1,]*(1-p[z1,])^(1-Y[z1,])+phi*Y[z1,]))
	 		# sum(dbinom(Y[z1,],1,p[z1,],log=TRUE))
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
		p1 <- psi*apply((1-phi)*p^Y*(1-p)^(1-Y)+phi*Y,1,prod)
		p0 <- (1-psi)*apply(phi^Y*(1-phi)^(1-Y),1,prod)
		psi.tmp <- p1/(p1+p0)	
		z <- rbinom(n,1,psi.tmp)

		# p0.tmp <- psi*apply(p^Y*(1-p)^(1-Y),1,prod)
		# psi.tmp <- p0.tmp/(p0.tmp+(1-psi))	
		# z[y0] <- rbinom(n.y0,1,psi.tmp[y0])
	
	
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

	keep <- lapply(keep,function(x) x/n.mcmc)
	end <- list(beta=beta,alpha=alpha,z=z)  # starting values
	
	list(beta=beta.save,alpha.save=alpha.save,N.save=N.save,z.mean=z.mean,keep=keep,end=end,
		Y=Y,X=X,W=W,priors=priors,start=start,tune=tune,n.mcmc=n.mcmc)
}
