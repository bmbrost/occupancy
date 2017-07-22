occ.community.multiscale.mcmc <- function(Y,J,groups,W,U,X,priors,start,tune,n.mcmc,n.thin=1,adapt=TRUE){

	#
	#  	Brian Brost (28 JUN 2017)
	#	
	#	Community occupancy model with independent detection and occupancy intercepts
	#
	#  	Y: 
	#	J:
	#	groups:
	# 	W: 
	#	U:
	#  	X: 
	#	priors:
	#	start:
	#	tune:
	#	n.mcmc:
	#	adapt:
	#

	###
	###  Libraries and subroutines
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
			tmp <- ifelse(relative.psi==0,0,log(relative.psi))  # per email from K. Broms on 20170711
			out <- exp(-1*sum(relative.psi*tmp))
		}
		out
	}


	###
	###  Create variables 
	###
# browser()	
	n <- ncol(Y)  # number of species in community
	qW <- ncol(W) 
	qU <- ncol(U)
	qX <- ncol(X)

	Y.inv <- ifelse(Y==0,1,0)  # for updating use state (a)
	v.beta <- matrix(0,nrow(X),n)  # auxilliary variable for updating of beta
	v.gamma <- matrix(0,nrow(U),n)  # auxilliary variable for updating of gamma	
				
	# Create indicator variable that maps latent occupancy state (z) to use state (a)
	z.map <- match(paste(groups$U$unit,groups$U$time),paste(groups$X$unit,groups$X$time))


	###
	###  Starting values 
	###
	
	alpha <- matrix(start$alpha,qW,n)  # detection coefficients (p)
	gamma <- matrix(start$gamma,qU,n)  # use coefficientes (theta) 
	beta <- matrix(start$beta,qX,n)  # occupancy coefficients (psi)
	mu.alpha <- matrix(start$mu.alpha,qW,1)	 # mean for alpha
	mu.gamma <- matrix(start$mu.gamma,qU,1)	 # mean for gamma
	mu.beta <- matrix(start$mu.beta,qX,1)  # mean for beta
	p <- expit(W%*%alpha)  # detection probability
	theta <- t(apply(U,1,function(x) expit(x%*%gamma)))  # probability of use
	psi <- t(apply(X,1,function(x) pnorm(x%*%beta)))  # occupancy probability
	a <- start$a  # use state
	z <- start$z  # latent occupancy state	

	sigma.alpha <- start$sigma.alpha  # standard deviation of detection coefficients
	Sigma.alpha <- diag(qW)*sigma.alpha
	Sigma.alpha.inv <- solve(Sigma.alpha)	
	sigma.gamma <- start$sigma.gamma  # standard deviation of use coefficients
	Sigma.gamma <- diag(qU)*sigma.gamma^2
	Sigma.gamma.inv <- solve(Sigma.gamma)	
	sigma.beta <- start$sigma.beta  # standard deviation of occupancy coefficients
	Sigma.beta <- diag(qX)*sigma.beta^2
	Sigma.beta.inv <- solve(Sigma.beta)	

	
	###
	###  Priors
	###
	
	Sigma.mu.alpha <- diag(qW)*priors$sigma.mu.alpha^2  # prior variance-covariance for mu.alpha
	Sigma.mu.alpha.inv <- solve(Sigma.mu.alpha)	
	Sigma.mu.gamma <- diag(qU)*priors$sigma.mu.gamma^2  # prior variance-covariance for mu.gamma
	Sigma.mu.gamma.inv <- solve(Sigma.mu.gamma)	
	Sigma.mu.beta <- diag(qX)*priors$sigma.mu.beta^2  # prior variance-covariance for mu.beta
	Sigma.mu.beta.inv <- solve(Sigma.mu.beta)	

	
	###
	###  Create receptacles for output
	###
	
	alpha.save <- array(0,c(n.mcmc/n.thin,qW,n))	
	gamma.save <- array(0,c(n.mcmc/n.thin,qU,n))
	beta.save <- array(0,c(n.mcmc/n.thin,qX,n))
	a.mean <- matrix(0,nrow(U),n)
	z.mean <- matrix(0,nrow(X),n)
	mu.alpha.save <- matrix(0,n.mcmc/n.thin,qW)
	mu.gamma.save <- matrix(0,n.mcmc/n.thin,qU)
	mu.beta.save <- matrix(0,n.mcmc/n.thin,qX)
	sigma.alpha.save <- numeric(n.mcmc/n.thin)
	sigma.gamma.save <- numeric(n.mcmc/n.thin)
	sigma.beta.save <- numeric(n.mcmc/n.thin)
	richness.save <- matrix(0,n.mcmc/n.thin,nrow(z))
	hill0.save <- matrix(0,n.mcmc/n.thin,nrow(X))
	hill1.save <- matrix(0,n.mcmc/n.thin,nrow(X))
	hill2.save <- matrix(0,n.mcmc/n.thin,nrow(X))
	
	keep <- list(alpha=0)  # number of MH proposals accepted
	keep.tmp <- keep  # for adaptive tuning
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
			tune$alpha <- get.tune(tune$alpha,keep.tmp$alpha/n,k)
			keep.tmp <- lapply(keep.tmp,function(x) x*0)
	   	} 	
		
		for(i in 1:n){  # loop through species

			###
		  	###  Sample a (use state)
		  	###
# browser()
			z.tmp <- z[z.map,i]
			p1 <- (z.tmp*theta[,i])*(p[,i]^Y[,i])*((1-p[,i])^(J-Y[,i]))
			p0 <- (1-z.tmp*theta[,i])*Y.inv[,i]
			theta.tmp <- p1/(p1+p0)
			a[,i] <- rbinom(J,1,theta.tmp)
			# boxplot(p1~Y[,i]>0)
			# boxplot(p0~Y[,i]>0)
			# boxplot(theta.tmp~Y[,i]>0)
			# cbind(Y[,i],theta.tmp,p1,p0)
			
			
			###
	  		###  Sample alpha (detection)
		  	###

			alpha.star <- rnorm(qW,alpha[,i],tune$alpha)  # proposal
			p.star <- expit(W%*%alpha.star)  # detection probability
			a1 <- a[,i]==1  # update depends only on sites that are 'used'
		 	mh.star <- sum(dbinom(Y[a1,i],J[a1],p.star[a1],log=TRUE))+
				dnorm(alpha.star,mu.alpha,sigma.alpha,log=TRUE)
		 	mh.0 <- sum(dbinom(Y[a1,i],J[a1],p[a1,i],log=TRUE))+
				dnorm(alpha[,i],mu.alpha,sigma.alpha,log=TRUE)
			if(exp(mh.star-mh.0) > runif(1)){
				alpha[,i] <- alpha.star
				p[,i] <- p.star
				keep$alpha <- keep$alpha+1
				keep.tmp$alpha <- keep.tmp$alpha+1
		  	}			


		  	###
			###  Sample v.gamma (auxilliary variable for z) 
			###

			z.tmp <- z.tmp==1
			a1 <- which(a[,i]==1&z.tmp)
			a0 <- which(a[,i]==0&z.tmp)
			v.gamma[a0,i] <- truncnormsamp(U[a0,]%*%gamma[,i],1,-Inf,0,length(a0))
			v.gamma[a1,i] <- truncnormsamp(U[a1,]%*%gamma[,i],1,0,Inf,length(a1))


			###
			###  Sample gamma (theta)
			###

			v.gamma.tmp <- v.gamma[z.tmp,i]
			U.tmp <- U[z.tmp,]
			A.inv <- solve(t(U.tmp)%*%U.tmp+Sigma.gamma.inv)
			b <- t(U.tmp)%*%v.gamma.tmp+Sigma.gamma.inv%*%mu.gamma			
			gamma[,i] <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qU),qU,1)
			theta[,i] <- pnorm(U%*%gamma[,i])			
			

			###
			###  Sample v.beta (auxilliary variable for z) 
			###

			z0 <- z[,i]==0
			z1 <- z[,i]==1
			v.beta[z0,i] <- truncnormsamp(X[z0,]%*%beta[,i],1,-Inf,0,sum(z0))
			v.beta[z1,i] <- truncnormsamp(X[z1,]%*%beta[,i],1,0,Inf,sum(z1))

			
			###
			###  Sample beta (psi)
			###

			A.inv <- solve(t(X)%*%X+Sigma.beta.inv)  # for update of beta		
			b <- t(X)%*%v.beta[,i]+Sigma.beta.inv%*%mu.beta			
			beta[,i] <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qX),qX,1)
			psi[,i] <- pnorm(X%*%beta[,i])			


		 	###
			###  Sample gamma (theta)
		  	###
	
			# gamma.star <- rnorm(qU,gamma[,i],tune$gamma)  # proposal
			# theta.star <- expit(U%*%gamma.star)  # probability of use
			# z1 <- z.tmp==1  # update depends only on sites within occupied units
			# mh.star <- sum(dbinom(a[z1,i],1,theta.star[z1],log=TRUE))+
				# sum(dnorm(gamma.star,mu.gamma,sigma.gamma,log=TRUE))		
			# mh.0 <- sum(dbinom(a[z1,i],1,theta[z1,i],log=TRUE))+
				# sum(dnorm(gamma[,i],mu.gamma,sigma.gamma,log=TRUE))			
			# if(exp(mh.star-mh.0) > runif(1)){
				# gamma[,i] <- gamma.star
				# theta[,i] <- theta.star
				# keep$gamma <- keep$gamma+1
				# keep.tmp$gamma <- keep.tmp$gamma+1
			# }


		  	###
			###  Sample beta (psi)
		  	###

			# beta.star <- rnorm(qX,beta[,i],tune$beta)  # proposal
			# psi.star <- expit(X%*%beta.star)  # occupancy probability
			# mh.star <- sum(dbinom(z[,i],1,psi.star,log=TRUE))+
				# sum(dnorm(beta.star,mu.beta,sigma.beta,log=TRUE))
			# mh.0 <- sum(dbinom(z[,i],1,psi[,i],log=TRUE))+
				# sum(dnorm(beta[,i],mu.beta,sigma.beta,log=TRUE))
			# if(exp(mh.star-mh.0) > runif(1)){
				# beta[,i] <- beta.star
				# psi[,i] <- psi.star
				# keep$beta <- keep$beta+1
				# keep.tmp$beta <- keep.tmp$beta+1
			# }

		}  # end loop through species  
		

		###
	  	###  Sample z (occupancy state)
	  	###

		a.inv <- ifelse(a==1,0,1)  
		p1 <- psi*apply((theta^a)*((1-theta)^(1-a)),2,function(x) tapply(x,z.map,prod))
		p0 <- (1-psi)*apply(a.inv,2,function(x) tapply(x,z.map,prod))
		psi.tmp <- p1/(p1+p0)	
		z <- t(apply(psi.tmp,1,function(x) rbinom(n,1,x)))
		# boxplot(p0~apply(a,2,function(x) tapply(x,z.map,sum))>0)
		# boxplot(psi.tmp~apply(a,2,function(x) tapply(x,z.map,sum))>0)
		
		
		###
		###  Sample mu.alpha
		###

		alpha.sum <- rowSums(alpha)
		A.inv <- solve(n*Sigma.alpha.inv+Sigma.mu.alpha.inv)
		b <- Sigma.alpha.inv%*%alpha.sum
	    mu.alpha <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qW),qW,1)

		
		###
		###  Sample mu.gamma
		###

		gamma.sum <- rowSums(gamma)
		A.inv <- solve(n*Sigma.gamma.inv+Sigma.mu.gamma.inv)
		b <- Sigma.gamma.inv%*%gamma.sum
	    mu.gamma <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qU),qU,1)
		
		
		###
		###  Sample mu.beta
		###

		beta.sum <- rowSums(beta)
		A.inv <- solve(n*Sigma.beta.inv+Sigma.mu.beta.inv)
		b <- Sigma.beta.inv%*%beta.sum
	    mu.beta <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qX),qX,1)
		# mu.beta <- t(rmvnorm(1,A.inv%*%b,A.inv))

		
		###
		###  Sample sigma.alpha
		###	

		tmp <- sum(apply(alpha,2,function(x) x-mu.alpha)^2)
		r.tmp <- 1/(tmp/2+1/priors$r)
		q.tmp <- ((qW)*n)/2+priors$q
		sigma.alpha <- 1/rgamma(1,q.tmp,,r.tmp)
		Sigma.alpha <- diag(qW)*sigma.alpha
		Sigma.alpha.inv <- solve(Sigma.alpha)				
		sigma.alpha <- sqrt(sigma.alpha)


		###
		###  Sample sigma.gamma
		###	

		tmp <- sum(apply(gamma,2,function(x) x-mu.gamma)^2)
		r.tmp <- 1/(tmp/2+1/priors$r)
		q.tmp <- ((qU)*n)/2+priors$q
		sigma.gamma <- 1/rgamma(1,q.tmp,,r.tmp)
		Sigma.gamma <- diag(qU)*sigma.gamma
		Sigma.gamma.inv <- solve(Sigma.gamma)				
		sigma.gamma <- sqrt(sigma.gamma)


		###
		###  Sample sigma.beta 
		###	

		tmp <- sum(apply(beta,2,function(x) x-mu.beta)^2)
		r.tmp <- 1/(tmp/2+1/priors$r)
		q.tmp <- ((qX)*n)/2+priors$q
		sigma.beta <- 1/rgamma(1,q.tmp,,r.tmp)
		Sigma.beta <- diag(qX)*sigma.beta
		Sigma.beta.inv <- solve(Sigma.beta)				
		sigma.beta <- sqrt(sigma.beta)


		###
		### Calculate unit-level diversity using Hill numbers a la Broms et al. (2015)
		###
		
		hill0 <- rowSums(z)
		# hill0 <- apply(psi*z,1,function(x) get.hill(x[x!=0],q=0))  # richness
		hill1 <- apply(psi*z,1,function(x) get.hill(x,q=1))  # Shannon diversity
		hill2 <- apply(psi*z,1,function(x) get.hill(x,q=2))  # Simpson diversity
		
		
		###
	  	###  Save samples 
	  	###

		if(k%%n.thin==0){
			k.tmp <- k/n.thin
		  	alpha.save[k.tmp,,] <- alpha
		  	gamma.save[k.tmp,,] <- gamma
		  	beta.save[k.tmp,,] <- beta
		  	a.mean <- a.mean+a
		  	z.mean <- z.mean+z
		  	mu.alpha.save[k.tmp,] <- mu.alpha
		  	mu.gamma.save[k.tmp,] <- mu.gamma
		  	mu.beta.save[k.tmp,] <- mu.beta
			sigma.alpha.save[k.tmp] <- sigma.alpha
			sigma.gamma.save[k.tmp] <- sigma.gamma
			sigma.beta.save[k.tmp] <- sigma.beta
			richness.save[k.tmp,] <- rowSums(z)
			hill0.save[k.tmp,] <- hill0
			hill1.save[k.tmp,] <- hill1
			hill2.save[k.tmp,] <- hill2		
		}
				
	}
	cat("\n")
	
	###
	###  Write output 
	###

	z.mean <- z.mean/(n.mcmc/n.thin)
  	a.mean <- a.mean/(n.mcmc/n.thin)
	
	keep <- lapply(keep,function(x) x/(n.mcmc*n))

	end <- list(beta=beta,gamma=gamma,alpha=alpha,z=z,a=a,mu.alpha=mu.alpha,  # ending values
		mu.gamma=mu.gamma,mu.beta=mu.beta,
		sigma.alpha=sigma.alpha,sigma.gamma=sigma.gamma,sigma.beta=sigma.beta)
	
	list(alpha=alpha.save,gamma=gamma.save,beta=beta.save,a.mean=a.mean,z.mean=z.mean,
		mu.alpha=mu.alpha.save,mu.gamma=mu.gamma.save,mu.beta=mu.beta.save,
		sigma.alpha=sigma.alpha.save,sigma.gamma=sigma.gamma.save,sigma.beta=sigma.beta.save,
		richness=richness.save,hill0=hill0.save,hill1=hill1.save,hill2=hill2.save,keep=keep,end=end,
		Y=Y,J=J,groups=groups,W=W,U=U,X=X,priors=priors,start=start,tune=tune,n.mcmc=n.mcmc,n.thin=n.thin,
		adapt=adapt)
}
