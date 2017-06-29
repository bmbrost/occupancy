occ.community.mcmc <- function(Y,J,W,X,priors,start,tune,n.mcmc,adapt=TRUE){

	#
	#  	Brian Brost (28 JUN 2017)
	#	
	#	Community occupancy model with independent detection and occupancy intercepts
	#
	#  	Y: 
	#	J:
	# 	W: 
	#  	X: 
	#	priors:
	#	start:
	#	tune:
	#	n.mcmc:
	#	adapt:
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

	get.hill <- function(psi,q){  # calculate Hill numbers a la Broms et al. (2015)
		sum.psi <- sum(psi)
		relative.psi <- psi/sum.psi
		if(q!=1){
			out <- (sum(relative.psi^q))^(1/(1-q))
		}
		if(q==1){
			out <- exp(-1*sum(relative.psi*log(relative.psi)))
		}
		out
	}


	###
	###  Create Variables 
	###
# browser()
	n <- ncol(Y)  # number of species in community
	R <- nrow(Y)  # number of sites and sampling periods
	qX <- ncol(X)
	qW <- ncol(W)
	
	Y.inv <- ifelse(Y==0,1,0)  # for updating occupancy state (z)
	v <- matrix(0,nrow=R,ncol=n)  # auxilliary variable for updating beta


	###
	###  Starting Values 
	###

	alpha <- matrix(start$alpha,qW,n)  # detection coefficients (p)
	beta <- matrix(start$beta,qX,n)  # occupancy coefficients (psi)
	mu.alpha <- matrix(start$mu.alpha,qW,1)  # mean for alpha	
	mu.beta <- matrix(start$mu.beta,qX,1)  # mean for beta			
	p <- expit(W%*%alpha)  # detection probability
	z <- start$z  # latent occupancy state
	
	Sigma.beta <- diag(qX)*start$sigma.beta^2  # variance-covariance of occupancy coefficients
	Sigma.beta.inv <- solve(Sigma.beta)	
	
	sigma.alpha <- start$sigma.alpha  # standard deviation of detection coefficients
	Sigma.alpha <- diag(qW)*sigma.alpha^2   # variance-covariance of detection coefficients	
	Sigma.alpha.inv <- solve(Sigma.alpha)

	
	###
	###  Priors
	###
		
	Sigma.mu.alpha <- diag(qW)*priors$sigma.mu.alpha^2  # prior covariance for mu.alpha
	Sigma.mu.alpha.inv <- solve(Sigma.mu.alpha)
	Sigma.mu.beta <- diag(qX)*priors$sigma.mu.beta^2  # prior covariance for mu.beta
	Sigma.mu.beta.inv <- solve(Sigma.mu.beta)	
	
	
	###
	###  Create receptacles for output
	###

	
	alpha.save <- array(0,c(n.mcmc,qW,n))
	beta.save <- array(0,c(n.mcmc,qX,n))
	z.mean <- matrix(0,R,n)
	mu.alpha.save <- matrix(0,n.mcmc,qW)
	mu.beta.save <- matrix(0,n.mcmc,qX)
	sigma.alpha.save <- numeric(n.mcmc)
	sigma.beta.save <- numeric(n.mcmc)
	richness.save <- matrix(0,n.mcmc,R)
	hill0.save <- matrix(0,n.mcmc,R)
	hill1.save <- matrix(0,n.mcmc,R)
	hill2.save <- matrix(0,n.mcmc,R)

	keep <- list(alpha=0)
	keep.tmp <- keep
	Tb <- 50  # frequency of adaptive tuning


	###
	###  Begin MCMC Loop 
	###
	
	for(k in 1:n.mcmc){
		if(k%%1000==0) cat(k," "); flush.console()		

		###
		### Adaptive tuning
		###
		
		if(adapt==TRUE & k%%Tb==0) {  # Adaptive tuning
			keep.tmp <- lapply(keep.tmp,function(x) x/Tb)
			keep.tmp$alpha <- keep.tmp$alpha/n
			tune$alpha <- get.tune(tune$alpha,keep.tmp$alpha,k)
			keep.tmp <- lapply(keep.tmp,function(x) x*0)
	   	} 	

		for (i in 1:n){  # loop through species
# browser()
			###
			###  Sample v (auxilliary variable for z) 
			###
		
			z0 <- z[,i]==0
			z1 <- z[,i]==1
			v[z0,i] <- truncnormsamp(matrix(X[z0,],,qX)%*%beta[,i],1,-Inf,0,sum(z0))
			v[z1,i] <- truncnormsamp(matrix(X[z1,],,qX)%*%beta[,i],1,0,Inf,sum(z1))
				

			###
	  		###  Sample alpha (p)
		  	###

			alpha.star <- rnorm(qW,alpha[,i],tune$alpha)  # proposal
			p.star <- expit(W%*%alpha.star)  # detection probability
		 	mh.star <- sum(dbinom(Y[z1,i],J[z1],p.star[z1],log=TRUE))+
				dnorm(alpha.star,mu.alpha,sigma.alpha,log=TRUE)
		 	mh.0 <- sum(dbinom(Y[z1,i],J[z1],p[z1,i],log=TRUE))+
		 		dnorm(alpha[,i],mu.alpha,sigma.alpha,log=TRUE)
			if(exp(mh.star-mh.0) > runif(1)){
				alpha[,i] <- alpha.star
				p[,i] <- p.star
				keep$alpha <- keep$alpha+1
				keep.tmp$alpha <- keep.tmp$alpha+1
		  	}			
	
	
			###
			###  Sample beta (psi)
			###

			A.inv <- solve(t(X)%*%X+Sigma.beta.inv)  # for update of beta
			b <- t(X)%*%v[,i]+Sigma.beta.inv%*%mu.beta			
			beta[,i] <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qX),qX,1)
			psi[,i] <- pnorm(X%*%beta[,i])
	
	
			###
			###  Sample z 
			###

			p1 <- psi[,i]*p[,i]^Y[,i]*(1-p[,i])^(J-Y[,i])
			p0 <- (1-psi[,i])*Y.inv[,i]
			psi.tmp <- p1/(p1+p0)
			z[,i] <- rbinom(R,1,psi.tmp)					

		}  # end loop through species
		

		###
		###  Sample mu.alpha
		###

		A.inv <- solve(n*Sigma.alpha.inv+Sigma.mu.alpha.inv)
		b <- Sigma.alpha.inv%*%rowSums(alpha)
	    mu.alpha <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qW),qW,1)
	    # mu.alpha <- rnorm(1,A.inv*b,A.inv)


		###
		###  Sample mu.beta
		###

		A.inv <- solve(n*Sigma.beta.inv+Sigma.mu.beta.inv)
		b <- Sigma.beta.inv%*%rowSums(beta)
	    mu.beta <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qX),qX,1)
		# mu.beta <- t(rmvnorm(1,A.inv%*%b,A.inv))


		###
		###  Sample sigma.alpha
		###

		tmp <- sum(apply(alpha,2,function(x) x-mu.alpha)^2)
		r.tmp <- 1/(tmp/2+1/priors$r)
		q.tmp <- (qW*n)/2+priors$q
		sigma.alpha <- 1/rgamma(1,q.tmp,,r.tmp)
		Sigma.alpha <- diag(qW)*sigma.alpha
		Sigma.alpha.inv <- solve(Sigma.alpha)
		sigma.alpha <- sqrt(sigma.alpha)


		###
		###  Sample sigma.beta 
		###	

		tmp <- sum(apply(beta,2,function(x) x-mu.beta)^2)
		r.tmp <- 1/(tmp/2+1/priors$r)
		q.tmp <- (qX*n)/2+priors$q
		sigma.beta <- 1/rgamma(1,q.tmp,,r.tmp)
		Sigma.beta <- diag(qX)*sigma.beta
		Sigma.beta.inv <- solve(Sigma.beta)				
		sigma.beta <- sqrt(sigma.beta)


		###
		### Calculate Hill numbers a la Broms et al. (2015)
		###
		
		hill0 <- apply(psi,1,function(x) get.hill(x,q=0))  # richness
		hill1 <- apply(psi,1,function(x) get.hill(x,q=1))  # Shannon diversity
		hill2 <- apply(psi,1,function(x) get.hill(x,q=2))  # Simpson diversity
		
		
		###
		###  Save samples 
		###

		alpha.save[k,,] <- alpha
		beta.save[k,,] <- beta
		z.mean <- z.mean+z/n.mcmc
		mu.alpha.save[k,] <- mu.alpha
		mu.beta.save[k,] <- mu.beta
		sigma.beta.save[k] <- sigma.beta
		sigma.alpha.save[k] <- sigma.alpha
		richness.save[k,] <- rowSums(z)
		hill0.save[k,] <- hill0
		hill1.save[k,] <- hill1
		hill2.save[k,] <- hill2		
		
	}
	cat("\n")
	
	###
	###  Write output 
	###
	
	keep$alpha <- keep$alpha/n
	keep <- lapply(keep,function(x) x/n.mcmc)

	end <- list(alpha=alpha,beta=beta,z=z,mu.alpha=mu.alpha,mu.beta=mu.beta,
		sigma.alpha=sigma.alpha,sigma.beta=sigma.beta)  # ending values

	list(alpha=alpha.save,beta=beta.save,z.mean=z.mean,mu.alpha=mu.alpha.save,mu.beta=mu.beta.save,
		sigma.alpha=sigma.alpha.save,sigma.beta=sigma.beta.save,
		richness=richness.save,hill0=hill0.save,hill1=hill1.save,hill2=hill2.save,keep=keep,end=end,
		Y=Y,W=W,X=X,priors=priors,start=start,tune=tune,n.mcmc=n.mcmc,adapt=adapt)
}

