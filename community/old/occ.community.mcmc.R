occ.community.mcmc <- function(Y,J,W,X,priors,start,tune,n.mcmc,adapt=TRUE){

	#
	#  	Brian Brost (20 JUN 2017)
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


	###
	###  Create Variables 
	###
# browser()
	n <- nrow(Y)
	R <- ncol(Y)
	qX <- ncol(X)
	qW <- ncol(W)
	v <- matrix(0,n,R)

	
	###
	###  Priors
	###
		
	Sigma.mu.alpha <- diag(qW)*priors$sigma.mu.alpha^2
	Sigma.mu.alpha.inv <- solve(Sigma.mu.alpha)
	Sigma.mu.beta <- diag(qX)*priors$sigma.mu.beta^2
	Sigma.mu.beta.inv <- solve(Sigma.mu.beta)	

	
	###
	###  Starting Values 
	###

	alpha <- start$alpha
	beta <- start$beta
	mu.alpha <- matrix(start$mu.alpha,qW,1)	
	mu.beta <- matrix(start$mu.beta,qX,1)		
	p <- expit(alpha)
	z <- start$z
	
	Sigma.beta <- diag(qX)*start$sigma.beta^2
	Sigma.beta.inv <- solve(Sigma.beta)	
	
	sigma.alpha <- start$sigma.alpha
	Sigma.alpha <- diag(qW)*sigma.alpha^2	
	Sigma.alpha.inv <- solve(Sigma.alpha)
	
	
	###
	###  Create receptacles for output
	###

	
	alpha.save <- array(0,c(n.mcmc,qW,n))
	beta.save <- array(0,c(n.mcmc,qX,n))
	z.mean <- matrix(0,n,R)
	mu.alpha.save <- numeric(n.mcmc)
	mu.beta.save <- matrix(0,n.mcmc,qX)
	sigma.alpha.save <- numeric(n.mcmc)
	sigma.beta.save <- numeric(n.mcmc)
	N.save <- matrix(0,n.mcmc,R)

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
# browser()	   		   			
		alpha.star <- rnorm(n,alpha,tune$alpha)  # proposals for alpha_i
		p.star <- expit(alpha.star)  # proposals for p_i
		A <- solve(t(X)%*%X+Sigma.beta.inv)  # for update of beta

		for (i in 1:n){  # loop through species
# i <- 1
			###
			###  Sample v (auxilliary variable for z) 
			###
			
			z0 <- z[i,]==0
			z1 <- z[i,]==1
			v[i,z0] <- truncnormsamp(matrix(X[z0,],,qX)%*%beta[,i],1,-Inf,0,sum(z0))
			v[i,z1] <- truncnormsamp(matrix(X[z1,],,qX)%*%beta[,i],1,0,Inf,sum(z1))
				

			###
	  		###  Sample alpha_i and updata p_i
		  	###
# browser()
			z1 <- z[i,]==1	
		 	mh.star <- sum(dbinom(Y[i,z1],J[z1],p.star[i],log=TRUE))+
				dnorm(alpha.star[i],mu.alpha,sigma.alpha,log=TRUE)
		 	mh.0 <- sum(dbinom(Y[i,z1],J[z1],p[i],log=TRUE))+
		 		dnorm(alpha[i],mu.alpha,sigma.alpha,log=TRUE)
			if(exp(mh.star-mh.0) > runif(1)){
				alpha[i] <- alpha.star[i]
				p[i] <- p.star[i]
				keep$alpha <- keep$alpha+1
				keep.tmp$alpha <- keep.tmp$alpha+1
		  	}			
	
	
			###
			###  Sample beta_i and update psi_i
			###
# browser()	
			b <- t(X)%*%v[i,]+Sigma.beta.inv%*%mu.beta			
			beta[,i] <- A%*%b+t(chol(A))%*%matrix(rnorm(qX),qX,1)
			psi[i,] <- pnorm(X%*%beta[,i])
	
	
			###
			###  Sample z 
			###
# browser()			
			idx <- Y[i,]==0
			p1 <- psi[i,]*p[i]^Y[i,]*(1-p[i])^(J-Y[i,])
			p0 <- 1-psi[i,]
			psi.tmp <- p1/(p1+p0)
			z[i,idx] <- rbinom(sum(idx),1,psi.tmp[idx])					

		}  # end loop through species
		

		###
		###  Sample mu.alpha
		###

		A.inv <- solve(n*Sigma.alpha.inv+Sigma.mu.alpha.inv)
		b <- Sigma.alpha.inv%*%sum(alpha)
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
# browser()		
		tmp <- sum(apply(alpha,2,function(x) x-mu.alpha)^2)
		r.tmp <- 1/(tmp/2+1/priors$r)
		q.tmp <- n/2+priors$q
		sigma.alpha <- 1/rgamma(1,q.tmp,,r.tmp)
		Sigma.alpha <- diag(qW)*sigma.alpha
		Sigma.alpha.inv <- solve(Sigma.alpha)
		sigma.alpha <- sqrt(sigma.alpha)


		###
		###  Sample sigma.beta 
		###	
# browser()
		tmp <- sum(apply(beta,2,function(x) x-mu.beta)^2)
		r.tmp <- 1/(tmp/2+1/priors$r)
		q.tmp <- (qX*n)/2+priors$q
		sigma.beta <- 1/rgamma(1,q.tmp,,r.tmp)
		Sigma.beta <- diag(qX)*sigma.beta
		Sigma.beta.inv <- solve(Sigma.beta)				
		sigma.beta <- sqrt(sigma.beta)

		
		###
		###  Save samples 
		###

		alpha.save[k,,] <- alpha
		beta.save[k,,] <- beta
		z.mean <- z.mean+z/n.mcmc
		mu.alpha.save[k] <- mu.alpha
		mu.beta.save[k,] <- mu.beta
		sigma.beta.save[k] <- sigma.beta
		sigma.alpha.save[k] <- sigma.alpha
		N.save[k,] <- colSums(z)
	
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
		sigma.alpha=sigma.alpha.save,sigma.beta=sigma.beta.save,N=N.save,keep=keep,
		Y=Y,W=W,X=X,priors=priors,start=start,tune=tune,n.mcmc=n.mcmc,adapt=adapt)
}

