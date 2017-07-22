occ.community.correlated.mcmc <- function(Y,J,W,X,priors,start,tune,n.mcmc,n.thin=1,adapt=TRUE){

	#
	#  	Brian Brost (28 JUN 2017)
	#	
	#	Community occupancy model where covariance between detection and occupancy intercepts is modeled
	#
	#  	Y: 
	#	n: number of species in Y
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
			tmp <- ifelse(relative.psi==0,0,log(relative.psi))  # per email from K. Broms on 20170711
			out <- exp(-1*sum(relative.psi*tmp))
		}
		out
	}


	###
	###  Create Variables 
	###
# browser()
	n <- ncol(Y)  # number of species in community
	R <- nrow(Y)  # number of units and sampling periods
	qW <- ncol(W)
	qX <- ncol(X)
		
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
	psi <- t(apply(X,1,function(x) pnorm(x%*%beta)))  # occupancy probabilities
	z <- start$z  # latent occupancy state
	mu.0 <- matrix(c(mu.alpha[1],mu.beta[1]),2,1)  # mean of detection and occupancy intercepts
	
	Sigma <- start$Sigma  # variance-covariance of detection and occupancy intercepts
	Sigma.inv <- solve(Sigma)
	Sigma.beta <- diag(qX-1)*start$sigma.beta^2  # variance-covariance of occupancy coefficients
	Sigma.beta.inv <- solve(Sigma.beta)	


	###
	###  Priors
	###
		
	S0 <- priors$S0  # Wishart scale matrix
	nu <- priors$nu  # Wishart df
	Sigma.0 <- diag(2)*priors$sigma.mu.0^2  # prior covariance for mu.0
	Sigma.0.inv <- solve(Sigma.0)
	Sigma.mu.beta <- diag(qX-1)*priors$sigma.mu.beta^2  # prior covariance for mu.beta
	Sigma.mu.beta.inv <- solve(Sigma.mu.beta)	
		
	
	###
	###  Create receptacles for output
	###

	alpha.save <- array(0,c(n.mcmc/n.thin,qW,n))	
	beta.save <- array(0,c(n.mcmc/n.thin,qX,n))
	z.mean <- matrix(0,R,n)
	mu.alpha.save <- matrix(0,n.mcmc/n.thin,qW)
	mu.beta.save <- matrix(0,n.mcmc/n.thin,qX)
	Sigma.save <- array(0,dim=c(2,2,n.mcmc/n.thin))
	sigma.beta.save <- numeric(n.mcmc/n.thin)
	hill0.save <- matrix(0,n.mcmc/n.thin,R)
	hill1.save <- matrix(0,n.mcmc/n.thin,R)
	hill2.save <- matrix(0,n.mcmc/n.thin,R)

	keep <- list(alpha=0,beta=0)
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
			keep.tmp <- lapply(keep.tmp,function(x) x/(Tb))
			tune$alpha <- get.tune(tune$alpha,keep.tmp$alpha/n,k)
			tune$beta <- get.tune(tune$beta,keep.tmp$beta/n,k)
			keep.tmp <- lapply(keep.tmp,function(x) x*0)
	   	} 		
		
		for (i in 1:n){  # loop through species

			###
			###  Sample v (auxilliary variable for z) 
			###
# browser()			
			z0 <- z[,i]==0
			z1 <- z[,i]==1
			v[z0,i] <- truncnormsamp(X[z0,]%*%beta[,i],1,-Inf,0,sum(z0))
			v[z1,i] <- truncnormsamp(X[z1,]%*%beta[,i],1,0,Inf,sum(z1))


			###
	  		###  Sample alpha_0 (detection intercept)
		  	###

			alpha.star <- alpha[,i]
			alpha.star[1] <- rnorm(1,alpha[1,i],tune$alpha)  # proposal
			p.star <- expit(W%*%alpha.star)  # detection probability
		 	mh.star <- sum(dbinom(Y[z1,i],J[z1],p.star[z1],log=TRUE))+
				dmvnorm(c(alpha.star[1],beta[1,i]),mu.0,Sigma,log=TRUE)
		 	mh.0 <- sum(dbinom(Y[z1,i],J[z1],p[z1,i],log=TRUE))+
				dmvnorm(c(alpha[1,i],beta[1,i]),mu.0,Sigma,log=TRUE)
			if(exp(mh.star-mh.0) > runif(1)){
				alpha[,i] <- alpha.star
				p[,i] <- p.star
				keep$alpha <- keep$alpha+1
				keep.tmp$alpha <- keep.tmp$alpha+1
		  	}			
			

			###
			### Sample beta_0 (occupancy intercept)
			###

			beta.star <- beta[,i]
			beta.star[1] <- rnorm(1,beta[1,i],tune$beta)  # proposal
			mh.star <- sum(dnorm(v[,i],X%*%beta.star,1,log=TRUE),na.rm=TRUE)+
				dmvnorm(c(alpha[1,i],beta.star[1]),mu.0,Sigma,log=TRUE)
			mh.0 <- sum(dnorm(v[,i],X%*%beta[,i],1,log=TRUE),na.rm=TRUE)+
				dmvnorm(c(alpha[1,i],beta[1,i]),mu.0,Sigma,log=TRUE)
			if(exp(mh.star-mh.0) > runif(1)){
				beta[,i] <- beta.star
				keep$beta <- keep$beta+1
				keep.tmp$beta <- keep.tmp$beta+1
		  	}			

			
			###
			###  Sample beta_i (psi)
			###

			A.inv <- solve(t(X[,-1])%*%X[,-1]+Sigma.beta.inv)
			b <- t(X[,-1])%*%c(v[,i]-beta[1,i])+Sigma.beta.inv%*%mu.beta[-1]
			beta[-1,i] <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qX-1),qX-1,1)
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
		###  Sample mu.0 (alpha_0 and beta_0)
		###

		beta.sum <- rowSums(beta)
		alpha.sum <- rowSums(alpha)
		A.inv <- solve(n*Sigma.inv+Sigma.0.inv)
		b <- Sigma.inv%*%c(alpha.sum[1],beta.sum[1])
	    mu.0 <- matrix(A.inv%*%(b)+t(chol(A.inv))%*%matrix(rnorm(2),2,1),2,1)
		mu.alpha[1,1] <- mu.0[1,1]	
		mu.beta[1,1] <- mu.0[2,1]


		###
		###  Sample mu.beta
		###

		A.inv <- solve(n*Sigma.beta.inv+Sigma.mu.beta.inv)
		b <- Sigma.beta.inv%*%beta.sum[-1]
	    mu.beta[-1,1] <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qX-1),qX-1,1)
		# mu.beta <- t(rmvnorm(1,A.inv%*%b,A.inv))
	
	
		###
		### Sample Sigma
		###

		tmp <- apply(cbind(alpha[1,],beta[1,]),1,function(x) crossprod(t(x-mu.0)))
		Sn <- S0+matrix(rowSums(tmp),2,2)
		Sigma <- solve(rWishart(1,nu+n,solve(Sn))[,,1])
		Sigma.inv <- solve(Sigma)
		

		###
		###  Sample sigma.beta 
		###	

		tmp <- apply(beta,2,function(x) x-mu.beta)^2
		tmp <- sum(tmp[-1,])
		r.tmp <- 1/(tmp/2+1/priors$r)
		q.tmp <- ((qX-1)*n)/2+priors$q
		sigma.beta <- 1/rgamma(1,q.tmp,,r.tmp)
		Sigma.beta <- diag(qX-1)*sigma.beta
		Sigma.beta.inv <- solve(Sigma.beta)				
		sigma.beta <- sqrt(sigma.beta)

		# if(qW>1){  # Add updates if detection varies by more than just species
# browser()
			# ###
			# ###  Sample mu.alpha
			# ###


			# ###
			# ###  Sample sigma.alpha
			# ###	
			
		# }

		
		###
		### Calculate Hill numbers a la Broms et al. (2015)
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
			beta.save[k.tmp,,] <- beta
			z.mean <- z.mean+z
			mu.alpha.save[k.tmp,] <- mu.alpha
			mu.beta.save[k.tmp,] <- mu.beta
			Sigma.save[,,k.tmp] <- Sigma
			sigma.beta.save[k.tmp] <- sigma.beta
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

	keep$alpha <- keep$alpha/n
	keep$beta <- keep$beta/n
	keep <- lapply(keep,function(x) x/n.mcmc)

	end <- list(alpha=alpha,beta=beta,z=z,mu.alpha=mu.alpha,mu.beta=mu.beta,  # ending values
		Sigma=Sigma,sigma.beta=sigma.beta)

	list(alpha=alpha.save,beta=beta.save,z.mean=z.mean,mu.alpha=mu.alpha.save,mu.beta=mu.beta.save,
		Sigma=Sigma.save,sigma.beta=sigma.beta.save,
		hill0=hill0.save,hill1=hill1.save,hill2=hill2.save,keep=keep,end=end,
		Y=Y,W=W,X=X,priors=priors,start=start,tune=tune,n.mcmc=n.mcmc,n.thin=n.thin,adapt=adapt)
}
