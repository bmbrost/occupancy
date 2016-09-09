occ.fp.mcmc <- function(Y,W,X,priors,start,tune,n.mcmc,adapt=TRUE){

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
# Y.long <- c(Y)
# W.long <- apply(W,2,I)

	
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
	Q <- start$Q
	psi <- expit(X%*%beta)
	p <- apply(W,3,function(x) expit(x%*%alpha))
	pi <- start$pi
	
	###
	###  Create receptacles for output
	###
	
	beta.save <- matrix(0,n.mcmc,qX)
	alpha.save <- matrix(0,n.mcmc,qW)
	z.mean <- numeric(n)
	N.save <- numeric(n.mcmc)
	Q.mean <- matrix(n,qW)

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
		###  Sample pi
	  	###
	  	
	  	
	  	###
		###  Sample beta (psi)
	  	###
	
	  	beta.star <- rnorm(qX,beta,tune$beta)
		psi.star <- expit(X%*%beta.star)
		mh.star <- sum(dbinom(z,1,psi.star,log=TRUE))+
			sum(dnorm(beta.star,mu.beta,sigma.beta,log=TRUE))
		mh.0 <- sum(dbinom(z,1,psi,log=TRUE))
			+sum(dnorm(beta,mu.beta,sigma.beta,log=TRUE))
		if(exp(mh.star-mh.0) > runif(1)){
			beta <- beta.star
			psi <- psi.star
			keep$beta <- keep$beta+1
			keep.tmp$beta <- keep.tmp$beta+1
		}

		 	
		###
	  	###  Sample z 
	  	###
browser()
		
		# Samples for instances where z==1
		z1 <- which(z==1)
		p.tmp <- p[z1,]
		Y.tmp <- Y[z1,]
		Q.tmp <- Q[z1,]
		psi.tmp <- psi[z1,]
	 	p1 <- (psi.tmp*p.tmp^Y.tmp*(1-p.tmp)^(1-Y.tmp))^(1-Q.tmp)*
	 		(pi^Y.tmp*(1-pi)^(1-Y.tmp))^(Q.tmp)
		p0 <- ((1-psi.tmp)*pi^Y.tmp*(1-pi)^(1-Y.tmp))
		psi.tmp <- p1/(p1+p0)	
		psi.tmp <- apply(psi.tmp,1,prod)
		z[z1] <- rbinom(length(z1),1,psi.tmp)

		# Samples for instances where z==0
		z0 <- which(z==0)
		p.tmp <- p[z0,]
		Y.tmp <- Y[z0,]
		psi.tmp <- psi[z0,]
		p0 <- (1-psi.tmp)*p.tmp^Y.tmp*(1-p.tmp)^(1-Y.tmp)
		psi.tmp <- psi.tmp/(psi.tmp+p0)
		psi.tmp <- apply(psi.tmp,1,prod)
		z[z0] <- rbinom(length(z0),1,psi.tmp)


	 	###
	 	### Sample Q
	 	###
# browser()	 	
		Q <- matrix(0,n,qW) 
		z0 <- which(z==0)
		z1 <- which(z==1)
	 	p1 <- pi*pi^Y*(1-pi)^(1-Y)
		p0 <- (1-pi)*pi^Y*(1-pi)^(1-Y)
	 	pi.tmp <- p1/(p1+p0)
	 	Q <- rbinom(n*qW,1,pi.tmp)
	 	
			
	
		###
	  	###  Sample p (alpha)
	  	###
# browser()
		# z1 <- z==1	
	  	# alpha.star <- rnorm(qW,alpha,tune$alpha)
	  	# p.star <- apply(W,3,function(x) expit(x%*%alpha.star))	 	
	 	# mh.star <- sum(dbinom(Y[z1,],1,p.star[z1,],log=TRUE))
	 		# +sum(dnorm(alpha.star,mu.alpha,sigma.alpha,log=TRUE))
	 	# mh.0 <- sum(dbinom(Y[z1,],1,p[z1,],log=TRUE))
	 		# +sum(dnorm(alpha,mu.alpha,sigma.alpha,log=TRUE))
		# if(exp(mh.star-mh.0) > runif(1)){
			# alpha <- alpha.star
			# p <- p.star
			# keep$alpha <- keep$alpha+1
			# keep.tmp$alpha <- keep.tmp$alpha+1
	  	# }


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
