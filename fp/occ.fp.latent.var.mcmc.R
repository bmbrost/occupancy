occ.fp.latent.var.mcmc <- function(Y,ctrl,W,X,priors,start,tune,n.mcmc,adapt=TRUE){

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
		a <- min(0.01,1/sqrt(k))
		# a <- min(0.025,1/sqrt(k))
		exp(ifelse(keep<target,log(tune)-a,log(tune)+a))
	}

	y.lik <- function(y,p,phi,log=FALSE){
		p.tmp <- p+phi-p*phi
		tmp <- p.tmp^y*(1-p.tmp)^(1-y)
		# tmp <- (1-phi)*p^y*(1-p)^(1-y)+phi*y
		if(log) tmp <- log(tmp)
		tmp
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
	# y0 <- which(y==0)
	# n.y0 <- length(y0)
	Y0 <- ifelse(Y==1,0,1)		
	
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
# browser()	
	beta <- as.vector(start$beta)
	alpha <- as.vector(start$alpha)
	z <- start$z
	psi <- expit(X%*%beta)
	p <- apply(W,3,function(x) expit(x%*%alpha))
	Q <- start$Q	
	
	###
	###  Create receptacles for output
	###
	
	beta.save <- matrix(0,n.mcmc,qX)
	alpha.save <- matrix(0,n.mcmc,qW)
	phi.save <- numeric(n.mcmc)
	z.mean <- numeric(n)
	N.save <- numeric(n.mcmc)
	Q.mean <- matrix(0,n,ncol(Y))

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
		###  Sample phi
	  	###

		# Update of phi using proper full-conditional distribution that includes Y		
		phi <- rbeta(1,sum(Q)+ctrl$v+priors$a,sum(1-Q)+ctrl$M-ctrl$v+priors$b)
		
		# Update of phi using 'cut' function to prevent feedback of Y in the model (see Plummer 2015)	
		# phi <- rbeta(1,ctrl$v+priors$a,ctrl$M-ctrl$v+priors$b)


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
	 	mh.star <- sum(log((dbinom(Y[z1,],1,p.star[z1,])^(1-Q[z1,]))))+
	 		sum(dnorm(alpha.star,mu.alpha,sigma.alpha,log=TRUE))
	 	mh.0 <-	sum(log((dbinom(Y[z1,],1,p[z1,])^(1-Q[z1,]))))+
	 		sum(dnorm(alpha,mu.alpha,sigma.alpha,log=TRUE))
	 	# mh.star <- sum(log((dbinom(Y,1,p.star)^(1-Q))^z))+
	 		# sum(dnorm(alpha.star,mu.alpha,sigma.alpha,log=TRUE))
	 	# mh.0 <-	sum(log((dbinom(Y,1,p)^(1-Q))^z))+
	 		# sum(dnorm(alpha,mu.alpha,sigma.alpha,log=TRUE))
		if(exp(mh.star-mh.0) > runif(1)){
			alpha <- alpha.star
			p <- p.star
			keep$alpha <- keep$alpha+1
			keep.tmp$alpha <- keep.tmp$alpha+1
	  	}

	 
		###
		### Sample Q
		###

# browser()
		p1 <- phi*Y^z*Y^(1-z)
		p0 <- (1-phi)*(p^Y*(1-p)^(1-Y))^z*Y0^(1-z)
		phi.tmp <- p1/(p1+p0)
		Q <- apply(phi.tmp,2,function(x) rbinom(n,1,x))
		# matrix(rbinom(n*ncol(Y),1,phi.tmp),n)
		
		###
	  	###  Sample z 
	  	###
# browser()
		p1 <- psi*apply(Y^Q*p^Y*((1-p)^(1-Y))^(1-Q),1,prod)
		p0 <- (1-psi)*apply(Y0^(1-Q)*Y^Q,1,prod)
		psi.tmp <- p1/(p1+p0)	
		z <- rbinom(n,1,psi.tmp)

			
		###
	  	###  Save samples 
	  	###
	
	  	beta.save[k,] <- beta
	  	alpha.save[k,] <- alpha
	  	phi.save[k] <- phi
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
	end <- list(beta=beta,alpha=alpha,z=z,Q=Q,phi=phi)  # starting values
	
	list(beta=beta.save,alpha=alpha.save,phi=phi.save,N=N.save,
		z.mean=z.mean,Q.mean=Q.mean,keep=keep,end=end,Y=Y,X=X,W=W,
		priors=priors,start=start,tune=tune,n.mcmc=n.mcmc)
}
