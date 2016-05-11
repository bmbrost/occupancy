occ.dyn.1.mcmc <- function(site,Y,W,X,G,n.mcmc,z.true=NULL){

	#
	#  	Brian Brost (20 SEP 2014)
	#	
	#	Fully heterogeneous dynamic occupany model using auxilliary variables
	#	Note: this model does not include random effects to allow site-level variability
	#	See occ.dyn.1.sim.R script to simulate data according to this model
	#
	# 	Arguments to function:
	#		[N=number of sites, i indexes site, t indexes primary period, j indexes sampling event]
	#	  	Y: (T x N) x J matrix containing full detection history
	# 		W: (T x N) x qW x J matrix of covariates for detection probability (p)
	#	  	X: (T x N) x qX covariates for colonization (gamma)
	#		G: (T x N) x qG covariates for persistence (phi)
		
	
	####
	####  Libraries and Subroutines
	####
	
	library(mvtnorm)
  #library(abind)
	
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
	y <- apply(Y,1,sum,na.rm=TRUE)
	row.idx <- t(sapply(unique(site),function(x) range(which(site==x))))
	T <- as.numeric(table(site))
	J <- apply(Y,1,function(x) sum(!is.na(x)))
	J.max <- max(J,na.rm=TRUE)
	Y.vec <- c(Y)
	W.mat <- apply(W,2,I) #Flatten W array
	N <- length(T) #Number of sites
	qX <- ncol(X)
	qW <- ncol(W)
	qG <- ncol(G)

	beta.save <- matrix(0,n.mcmc,qX)
	theta.save <- matrix(0,n.mcmc,qG)
	alpha.save <- matrix(0,n.mcmc,qW)
	psi.save <- matrix(0,n.mcmc,N)
	z.mean <- numeric(sum(T))
	# N.save <- numeric(n.mcmc)
	v <- numeric(sum(T))
	q <- numeric(sum(T))
	# u <- matrix(0,n,dim(W)[3])	
	u <- numeric(sum(T)*dim(W)[3])	

	# Create indices for z_0, z_1, and z_T, and y_1 and y_T
	y_1.idx <- c(t(sapply(unique(site),function(x) min(which(site==x)))))
	y_T.idx <- c(t(sapply(unique(site),function(x) max(which(site==x)))))
	z_0.idx <- y_1.idx + 0:(N-1) 
	z_1.idx <- z_0.idx + 1
	z_T.idx <- y_T.idx + 1:N
	
	
	####
	####  Priors Specifications 
	####
	
	beta.mn <- matrix(0,qX,1)
	theta.mn <- matrix(0,qG,1)
	alpha.mn <- matrix(0,qW,1)
	beta.var <- diag(qX)*10
	theta.var <- diag(qG)*10
	alpha.var <- diag(qW)*10

	alpha.psi <- 0.25 #scale1 of beta distribution for psi
	beta.psi <- 1 #scale2 of beta distribution for psi
	
	####
	####  Starting Values 
	####

	# Starting values for z_0 (occupancy at t=0) and z
	z <- numeric(sum(T)+N)
	z[-z_0.idx] <- ifelse(y>0,1,0)
	z[z_0.idx] <- ifelse(y[z_1.idx]>0,1,0)
	
	# Starting values for persistence (phi)
	# z0.idx <- z[-(T+1)]==1&z[-1]==0 #idx for z_t-1=1 & z_t=0
	# z1.idx <- z[-(T+1)]==1&z[-1]==1 #idx for z_t-1=1 & z_t=1
	# tmp <- rep(NA,T)
	# tmp[z0.idx] <- 0
	# tmp[z1.idx] <- 1	
	# theta <- matrix(c(glm(tmp~0+G,family=binomial)$coefficients),qG,1)
	# theta <- matrix(c(-0.25,0.25),2,1)
  theta <- matrix(0,qG,1)
	phi <- pnorm(G%*%theta)
	

	# Starting values for colonization (gamma)
	# z0.idx <- z[-(T+1)]==0&z[-1]==0 #idx for z_t-1=0 & z_t=0
	# z1.idx <- z[-(T+1)]==0&z[-1]==1 #idx for z_t-1=0 & z_t=1
	# tmp <- rep(NA,T)
	# tmp[z0.idx] <- 0
	# tmp[z1.idx] <- 1	
	# beta <- matrix(c(glm(tmp~0+X,family=binomial)$coefficients),qX,1)
	# beta <- matrix(c(-0.75,0.05),2,1)
  beta <- matrix(0,qX,1)
	gamma <- pnorm(X%*%beta)

	# Starting values for detection (p)
	# idx <- which(y[-1]>0)
	# alpha <- matrix(c(glm(Y.vec[idx]~0+W.mat[idx,],family=binomial)$coefficients),qW,1)
	# alpha <- matrix(c(-1,-0.075),2,1)
  alpha <- matrix(0,qW,1)
	p <- apply(W,3,function(x) pnorm(x%*%alpha))

	# Starting values for persistence (psi)
	# psi <- tapply(z[-z_0.idx],site,sum)/T
	psi <- rep(0.5,N)
		
	####
	####  Begin MCMC Loop 
	####
	
	for(k in 1:n.mcmc){
		if(k%%100==0) cat(k," "); flush.console()
		
		####
		#### Sample parameters related to colonization (gamma), i.e., when z_t-1=0 for t=1:T
		####

		# Sample v (auxilliary variable for gamma)
		z.test <- rep(NA,length(z))
		z.test[-z_0.idx] <- z[-z_T.idx]==0&z[-z_0.idx]==0 #idx for z_t-1=0 & z_t=0
		z.test[-z_0.idx] <- z[-z_T.idx]==0&z[-z_0.idx]==1 #idx for z_t-1=0 & z_t=1
		# data.frame(t=unlist(sapply(1:N,function(x) 0:T[x])),z,z.test)
		
		z0.idx <- z[-z_T.idx]==0&z[-z_0.idx]==0 #idx for z_t-1=0 & z_t=0
		z1.idx <- z[-z_T.idx]==0&z[-z_0.idx]==1 #idx for z_t-1=0 & z_t=1

# which(z1.idx)
		v[z0.idx] <- truncnormsamp(matrix(X[z0.idx,],,qX)%*%beta,1,-Inf,0,sum(z0.idx))
		v[z1.idx] <- truncnormsamp(matrix(X[z1.idx,],,qX)%*%beta,1,0,Inf,sum(z1.idx))		

		# Sample beta, calculate gamma
		idx <- z0.idx|z1.idx
		X.tmp <- X[idx,]
		v.tmp <- v[idx]
		A <- solve(t(X.tmp)%*%X.tmp+solve(beta.var))
		b <- t(X.tmp)%*%v.tmp+solve(beta.var)%*%beta.mn
		beta <- A%*%b+t(chol(A))%*%matrix(rnorm(qX),qX,1)
		# beta <- t(rmvnorm(1,A%*%b,A))
		gamma <- pnorm(X%*%beta)


		####
		#### Sample parameters related to persistence (phi), i.e., when z_t-1=1 for t=1:T
		####
			
		# Sample q (auxilliary variable for phi)
		z0.idx <- z[-z_T.idx]==1&z[-z_0.idx]==0 #idx for z_t-1=1 & z_t=0
		z1.idx <- z[-z_T.idx]==1&z[-z_0.idx]==1 #idx for z_t-1=1 & z_t=1
		q[z0.idx] <- truncnormsamp(matrix(G[z0.idx,],,qG)%*%theta,1,-Inf,0,sum(z0.idx))
		q[z1.idx] <- truncnormsamp(matrix(G[z1.idx,],,qG)%*%theta,1,0,Inf,sum(z1.idx))		

		# Sample theta, calculate phi
		idx <- z0.idx|z1.idx
		G.tmp <- G[idx,]
		q.tmp <- q[idx]
		A <- solve(t(G.tmp)%*%G.tmp+solve(theta.var))
		b <- t(G.tmp)%*%q.tmp+solve(theta.var)%*%theta.mn
		theta <- A%*%b+t(chol(A))%*%matrix(rnorm(qG),qG,1)
		# phi <- t(rmvnorm(1,A%*%b,A))
		phi <- pnorm(G%*%theta)


		####
		#### Sample parameters related to detection probability (p), i.e., when z_t=1 for t=1:T
		####
# browser()
			
		# Sample u (auxilliary variable for p)
		z1.idx <- z[-z_0.idx]==1
		z1.y0.idx <- which(z1.idx&Y.vec==0)
		z1.y1.idx <- which(z1.idx&Y.vec==1)
		u[z1.y0.idx] <- truncnormsamp((matrix(W.mat[z1.y0.idx,],,qW)%*%alpha),1,-Inf,0,length(z1.y0.idx))
		u[z1.y1.idx] <- truncnormsamp((matrix(W.mat[z1.y1.idx,],,qW)%*%alpha),1,0,Inf,length(z1.y1.idx))		
			
		# Sample alpha
		W.tmp <- W.mat[z1.idx,]		
		u.tmp <- u[z1.idx]
		A <- solve(t(W.tmp)%*%W.tmp+solve(alpha.var))
		b <- u.tmp%*%W.tmp+t(alpha.mn)%*%solve(alpha.var)
		alpha <- A%*%t(b)+t(chol(A))%*%matrix(rnorm(qW),qW,1)		
		p <- matrix(apply(W.mat,1,function(x) pnorm(x%*%alpha)),,J.max)


		####
		####  Sample psi, occupancy probability for z_0 (i.e., t=0) 
		####
# browser()				
		psi <- rbeta(1,(sum(z[z_0.idx])+alpha.psi),(sum(1-z[z_0.idx])+beta.psi))

		####
		####  Sample z_t, t=0:T, latent state variable for occupancy
		####
# browser()
		# Sample z_0 (i.e., for t=0) 

		z_1 <- z[z_0.idx+1] #z_1
		phi.tmp <- phi[y_1.idx]
		gamma.tmp <- gamma[y_1.idx]
		A <- (psi)*(phi.tmp^z_1)*((1-phi.tmp)^(1-z_1))
		B <- (1-psi)*(gamma.tmp^z_1)*((1-gamma.tmp)^(1-z_1))
		psi.tilde <- A/(A+B)
		z[z_0.idx] <- rbinom(N,1,psi.tilde)

# z[z_0.idx] <- z.true[z_0.idx]
# browser()
		y0.idx <- y==0
		y0.idx[y_T.idx] <- FALSE #set y_T to FALSE
		y0.idx <- which(y0.idx)
		z.idx <- y0.idx + site[y0.idx] #Matching index for z

# browser()			
		for(i in 1:length(y0.idx)){ #loop through primary sampling periods, t=1:(T-1)
			# i <- 1		
			# i <- 13
			# i <- y0.idx[1]
			# i <- y0.idx[length(y0.idx)]
			y.idx.tmp <- y0.idx[i]
			z.idx.tmp <- z.idx[i]
# print(c(y.idx.tmp,z.idx.tmp))
			z_tm1 <- z[z.idx.tmp-1] #z_t-1
			z_tp1 <- z[z.idx.tmp+1] #z_t+1
			phi.tmp <- phi[y.idx.tmp]
			gamma.tmp <- gamma[y.idx.tmp]
			p.tmp <- p[y.idx.tmp,]
			Y.tmp <- Y[y.idx.tmp,]
			p.tmp <- prod((p.tmp^Y.tmp)*(1-p.tmp)^(1-Y.tmp),na.rm=TRUE)
			A <- (p.tmp)*(phi.tmp^z_tm1)*(phi.tmp^z_tp1)*((1-phi.tmp)^(1-z_tp1))*(gamma.tmp^(1-z_tm1))
			B <- ((1-phi.tmp)^z_tm1)*(gamma.tmp^z_tp1)*((1-gamma.tmp)^(1-z_tm1))*((1-gamma.tmp)^(1-z_tp1))
			psi.tilde <- A/(A+B)
# print(psi.tilde)
			z[z.idx.tmp] <- rbinom(1,1,psi.tilde)
		}

		#  Sample z_T (t=T)
# browser()
		z_Tm1 <- z[z_T.idx-1] #z_T-1
		p.tmp <- matrix(p[y_T.idx,],N,J.max)
		Y.tmp <- matrix(Y[y_T.idx,],N,J.max)
		p.tmp <- sapply(1:N,function(x) prod((p.tmp[x,]^Y.tmp[x,])*(1-p.tmp[x,])^(1-Y.tmp[x,]),na.rm=TRUE))
		# p.tmp <- prod((p.tmp^Y.tmp)*(1-p.tmp)^(1-Y.tmp),na.rm=TRUE)
		phi.tmp <- phi[y_T.idx]
		gamma.tmp <- gamma[y_T.idx]
		A <- (p.tmp)*(phi.tmp^z_Tm1)*(gamma.tmp^(1-z_Tm1))
		B <- ((1-phi.tmp)^z_Tm1)*((1-gamma.tmp)^(1-z_Tm1))
		psi.tilde <- A/(A+B)
		z[z_T.idx] <- rbinom(N,1,psi.tilde)

# z[z_T.idx] <- z.true[z_T.idx]
# z <- z.true

		####
		####  Save Samples 
		####
	
		beta.save[k,] <- beta
		alpha.save[k,] <- alpha
		theta.save[k,] <- theta
		psi.save[k,] <- psi
		z.mean <- z.mean+z[-z_0.idx]/n.mcmc

		# N.save[k] <- sum(z)
	
		}
	cat("\n")
	
	####
	####  Write Output 
	####
	
	list(beta=beta.save,theta=theta.save,alpha=alpha.save,psi=psi.save,
		z.mean=z.mean,n.mcmc=n.mcmc)#N.save=N.save
}




### This model doesn't allow for multiple sites, and it appears
### something is up with its ability to estimate alpha, beta, and theta

# occ.dyn.1.mcmc <- function(Y,W,X,G,n.mcmc,z.true=NULL){

	# #
	# #  	Brian Brost (11 JUN 2014)
	# #	
	# #	Fully heterogeneous dynamic occupany model using auxilliary variables
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
	# library(abind)
	
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
	# y <- c(NA,apply(Y,1,sum,na.rm=TRUE))
	# J <- apply(Y,1,function(x) sum(!is.na(x)))
	# J.max <- max(J,na.rm=TRUE)
	# Y.vec <- c(Y)
	# W.mat <- apply(W,2,I)
	# T <- nrow(Y) #Number of primary periods
	# N <- 1 #Number of sites
	# qX <- ncol(X)
	# qW <- ncol(W)
	# qG <- ncol(G)

	# ###
	# ### Pad inputs for z_0
	# ###

	# # Y <- rbind(NA,Y)
	# # W <- abind(array(NA,dim=c(1,2,3)),W,along=1)
	# # X <- rbind(NA,X)
			
	# beta.save <- matrix(0,n.mcmc,qX)
	# theta.save <- matrix(0,n.mcmc,qX)
	# alpha.save <- matrix(0,n.mcmc,qW)
	# psi.save <- numeric(n.mcmc)
	# z.mean <- numeric(sum(T+1))
	# # N.save <- numeric(n.mcmc)
	# v <- numeric(sum(T))
	# q <- numeric(sum(T))
	# # u <- matrix(0,n,dim(W)[3])	
	# u <- numeric(T*dim(W)[3])	

	# ####
	# ####  Priors Specifications 
	# ####
	
	# beta.mn <- matrix(0,qX,1)
	# theta.mn <- matrix(0,qX,1)
	# alpha.mn <- matrix(0,qW,1)
	# beta.var <- diag(qX)*10
	# theta.var <- diag(qX)*10
	# alpha.var <- diag(qW)*10

	# alpha.psi <- 0.25 #scale1 of beta distribution for psi
	# beta.psi <- 1 #scale2 of beta distribution for psi
	
	# ####
	# ####  Starting Values 
	# ####

	# # Starting values for z_0 (occupancy at t=0) and z
	# z <- ifelse(y>0,1,0)
	# z[1] <- ifelse(z[2]==1,1,0)
	
	# # Starting values for persistence (phi)
	# z0.idx <- z[-(T+1)]==1&z[-1]==0 #idx for z_t-1=1 & z_t=0
	# z1.idx <- z[-(T+1)]==1&z[-1]==1 #idx for z_t-1=1 & z_t=1
	# tmp <- rep(NA,T)
	# tmp[z0.idx] <- 0
	# tmp[z1.idx] <- 1	
	# theta <- matrix(c(glm(tmp~0+G,family=binomial)$coefficients),qG,1)
	# phi <- pnorm(G%*%theta)

	# # Starting values for colonization (gamma)
	# z0.idx <- z[-(T+1)]==0&z[-1]==0 #idx for z_t-1=0 & z_t=0
	# z1.idx <- z[-(T+1)]==0&z[-1]==1 #idx for z_t-1=0 & z_t=1
	# tmp <- rep(NA,T)
	# tmp[z0.idx] <- 0
	# tmp[z1.idx] <- 1	
	# beta <- matrix(c(glm(tmp~0+X,family=binomial)$coefficients),qX,1)
	# gamma <- pnorm(X%*%beta)

	# # Starting values for detection (p)
	# idx <- which(y[-1]>0)
	# alpha <- matrix(c(glm(Y.vec[idx]~0+W.mat[idx,],family=binomial)$coefficients),qW,1)
	# p <- apply(W,3,function(x) pnorm(x%*%alpha))

	# # Starting values for persistence (psi)
	# psi <- sum(z[-1])/T
	
	
	# ####
	# ####  Begin MCMC Loop 
	# ####
	
	# for(k in 1:n.mcmc){
		# if(k%%100==0) cat(k," "); flush.console()
		
		# ####
		# #### Sample parameters related to colonization (gamma), i.e., when z_t-1=0 for t=1:T
		# ####
			
		# # Sample v (auxilliary variable for gamma)
		# z0.idx <- z[-(T+1)]==0&z[-1]==0 #idx for z_t-1=0 & z_t=0
		# z1.idx <- z[-(T+1)]==0&z[-1]==1 #idx for z_t-1=0 & z_t=1
		# # cbind(z,c(NA,z0.idx),c(NA,z1.idx))
		# v[z0.idx] <- truncnormsamp(matrix(X[z0.idx,],,qX)%*%beta,1,-Inf,0,sum(z0.idx))
		# v[z1.idx] <- truncnormsamp(matrix(X[z1.idx,],,qX)%*%beta,1,0,Inf,sum(z1.idx))		

		# # Sample beta, calculate gamma
		# idx <- z0.idx|z1.idx
		# X.tmp <- X[idx,]
		# v.tmp <- v[idx]
		# A <- solve(t(X.tmp)%*%X.tmp+solve(beta.var))
		# b <- t(X.tmp)%*%v.tmp+solve(beta.var)%*%beta.mn
		# beta <- A%*%b+t(chol(A))%*%matrix(rnorm(qX),qX,1)
		# # beta <- t(rmvnorm(1,A%*%b,A))
		# gamma <- pnorm(X%*%beta)


		# ####
		# #### Sample parameters related to persistence (phi), i.e., when z_t-1=1 for t=1:T
		# ####
			
		# # Sample q (auxilliary variable for phi)
		# z0.idx <- z[-(T+1)]==1&z[-1]==0 #idx for z_t-1=1 & z_t=0
		# z1.idx <- z[-(T+1)]==1&z[-1]==1 #idx for z_t-1=1 & z_t=1
		# # cbind(z,c(NA,z0.idx),c(NA,z1.idx))
		# q[z0.idx] <- truncnormsamp(matrix(G[z0.idx,],,qG)%*%theta,1,-Inf,0,sum(z0.idx))
		# q[z1.idx] <- truncnormsamp(matrix(G[z1.idx,],,qG)%*%theta,1,0,Inf,sum(z1.idx))		

		# # Sample theta, calculate phi
		# idx <- z0.idx|z1.idx
		# G.tmp <- G[idx,]
		# q.tmp <- q[idx]
		# A <- solve(t(G.tmp)%*%G.tmp+solve(theta.var))
		# b <- t(G.tmp)%*%q.tmp+solve(theta.var)%*%theta.mn
		# theta <- A%*%b+t(chol(A))%*%matrix(rnorm(qG),qG,1)
		# # phi <- t(rmvnorm(1,A%*%b,A))
		# phi <- pnorm(G%*%theta)


		# ####
		# #### Sample parameters related to detection probability (p), i.e., when z_t=1 for t=1:T
		# ####
# # browser()
			
		# # Sample u (auxilliary variable for p)
		# z1.idx <- z[-1]==1
		# z1.y0.idx <- which(z1.idx&Y.vec==0)
		# z1.y1.idx <- which(z1.idx&Y.vec==1)
		# u[z1.y0.idx] <- truncnormsamp((matrix(W.mat[z1.y0.idx,],,qW)%*%alpha),1,-Inf,0,length(z1.y0.idx))
		# u[z1.y1.idx] <- truncnormsamp((matrix(W.mat[z1.y1.idx,],,qW)%*%alpha),1,0,Inf,length(z1.y1.idx))		
			
		# # Sample alpha
		# W.tmp <- W.mat[z1.idx,]		
		# u.tmp <- u[z1.idx]
		# A <- solve(t(W.tmp)%*%W.tmp+solve(alpha.var))
		# b <- u.tmp%*%W.tmp+t(alpha.mn)%*%solve(alpha.var)
		# alpha <- A%*%t(b)+t(chol(A))%*%matrix(rnorm(qW),qW,1)		
		# p <- matrix(apply(W.mat,1,function(x) pnorm(x%*%alpha)),,J.max)


		# ####
		# ####  Sample psi, occupancy probability for z_0 (i.e., t=0) 
		# ####
# # browser()				
		# psi <- rbeta(1,z[1]+alpha.psi,(1-z[1]+beta.psi))


		# ####
		# ####  Sample z_t, t=0:T, latent state variable for occupancy
		# ####
# # browser()
		# # Sample z_0 (i.e., for t=0) 

		# z_1 <- z[2] #z_T-1
		# phi.tmp <- phi[1]
		# gamma.tmp <- gamma[1]
		# A <- (psi)*(phi.tmp^z_1)*((1-phi.tmp)^(1-z_1))
		# B <- (1-psi)*(gamma.tmp^z_1)*((1-gamma.tmp)^(1-z_1))
		# psi.tilde <- A/(A+B)
		# z[1] <- rbinom(1,1,psi.tilde)
# # z[1] <- z.true[1]

		# # Sample z_t when y_t==0 for (t=1:T-1)
		# y0.idx <- which(y[-(T+1)]==0)
		# for(i in y0.idx){ #loop through primary sampling periods, t=1:(T-1)
			# # i <- y0.idx[1]
			# # i <- y0.idx[length(y0.idx)]
			# z_tm1 <- z[i-1] #z_t-1
			# z_tp1 <- z[i+1] #z_t+1
			# phi.tmp <- phi[i-1]
			# gamma.tmp <- gamma[i-1]
			# p.tmp <- p[i-1,]
			# Y.tmp <- Y[i-1,]
			# p.tmp <- prod((p.tmp^Y.tmp)*(1-p.tmp)^(1-Y.tmp),na.rm=TRUE)
			# A <- (p.tmp)*(phi.tmp^z_tm1)*(phi.tmp^z_tp1)*((1-phi.tmp)^(1-z_tp1))*(gamma.tmp^(1-z_tm1))
			# B <- ((1-phi.tmp)^z_tm1)*(gamma.tmp^z_tp1)*((1-gamma.tmp)^(1-z_tm1))*((1-gamma.tmp)^(1-z_tp1))
			# psi.tilde <- A/(A+B)
			# z[i] <- rbinom(1,1,psi.tilde)
		# }
# # z <- z.true

		# #  Sample z_T (t=T)
		# z_Tm1 <- z[T] #z_T-1
		# p.tmp <- prod(p[T,]^Y[T,]*(1-p[T,])^(1-Y[T,]),na.rm=TRUE)
		# phi.tmp <- phi[T]
		# gamma.tmp <- gamma[T]
		# A <- (p.tmp)*(phi.tmp^z_Tm1)*(gamma.tmp^(1-z_Tm1))
		# B <- ((1-phi.tmp)^z_Tm1)*((1-gamma.tmp)^(1-z_Tm1))
		# psi.tilde <- A/(A+B)
		# z[T+1] <- rbinom(1,1,psi.tilde)
# # z[T+1] <- z.true[T+1]
# # z <- z.true
		# ####
		# ####  Save Samples 
		# ####
		
		# beta.save[k,] <- beta
		# alpha.save[k,] <- alpha
		# theta.save[k,] <- theta
		# psi.save[k] <- psi
		# z.mean <- z.mean+z/n.mcmc

		# # N.save[k] <- sum(z)
	
		# }
	# cat("\n")
	
	# ####
	# ####  Write Output 
	# ####
	
	# list(beta=beta.save,theta=theta.save,alpha=alpha.save,psi=psi.save,
		# z.mean=z.mean,n.mcmc=n.mcmc)#N.save=N.save
# }







### padded y with NA for z_0, working on z_t update
# occ.dyn.1.mcmc <- function(Y,W,X,n.mcmc){

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
	# y <- c(NA,apply(Y,1,sum,na.rm=TRUE))
	# T <- nrow(Y) #Number of primary periods
	# N <- 1 #Number of sites
	# J <- apply(Y,1,function(x) sum(!is.na(x)))
	# qX <- ncol(X)
	# qW <- ncol(W)
		
	# beta.save <- matrix(0,n.mcmc,qX)
	# theta.save <- matrix(0,n.mcmc,qX)
	# alpha.save <- matrix(0,n.mcmc,qW)
	# z.mean <- numeric(sum(T))
	# # N.save <- numeric(n.mcmc)
	# v <- numeric(sum(T))
	# q <- numeric(sum(T))
	# # u <- matrix(0,n,dim(W)[3])	
	# u <- numeric(T*dim(W)[3])	

	# ####
	# ####  Priors Specifications 
	# ####
	
	# beta.mn <- matrix(0,qX,1)
	# theta.mn <- matrix(0,qX,1)
	# alpha.mn <- matrix(0,qW,1)
	# beta.var <- diag(qX)*10
	# theta.var <- diag(qX)*10
	# alpha.var <- diag(qW)*10

	# ####
	# ####  Starting Values 
	# ####

	# # Starting values for persistence (phi)
	# theta <- matrix(c(0.75,0.1),2,1)
	# phi <- pnorm(X%*%theta)

	# # Starting values for colonization (gamma)
	# beta <- matrix(c(-0.5,0.1),2,1)
	# gamma <- pnorm(X%*%beta)

	# # Starting values for detection (p)
	# alpha <- matrix(c(0,-0.1),2,1)
	# p <- apply(W,3,function(x) pnorm(x%*%alpha))

	# # Starting values for persistence (psi)
	# psi <- 0.5
	
	# # Starting values for z_0 (occupancy at t=0) and z
	# z_0 <- 0
	# z <- c(z_0,ifelse(y>0,1,0))

# # browser()
	# # beta <- as.vector(glm(z ~ 0+X,family=binomial())$coefficients)
	# # alpha <- as.vector(glm(cbind(y[y>0],J[y>0]-y[y>0]) ~ 
		# # 0+W[y>0,,1],family=binomial())$coefficients)

# Y.vec <- c(Y)
# W.mat <- apply(W,2,I)


	# ####
	# ####  Begin MCMC Loop 
	# ####
	
	# for(k in 1:n.mcmc){
		# if(k%%1000==0) cat(k," "); flush.console()
		
		# ####
		# #### Sample parameters related to colonization (gamma), i.e., when z_t-1=0
		# ####
			
		# # Sample v (auxilliary variable for gamma)
# v <- numeric(T)
		# z0.idx <- z[-T]==0&z[-1]==0 #idx for z_t-1=0 & z_t=0
		# z1.idx <- z[-T]==0&z[-1]==1 #idx for z_t-1=0 & z_t=1
		# # cbind(z,c(NA,z0.idx),c(NA,z1.idx))
		# v[z0.idx] <- truncnormsamp(matrix(X[z0.idx,],,qX)%*%beta,1,-Inf,0,sum(z0.idx))
		# v[z1.idx] <- truncnormsamp(matrix(X[z1.idx,],,qX)%*%beta,1,0,Inf,sum(z1.idx))		

		# # Sample beta, calculate gamma
		# idx <- z0.idx|z1.idx
		# X.tmp <- X[idx,]
		# v.tmp <- v[idx]
		# A <- solve(t(X.tmp)%*%X.tmp+solve(beta.var))
		# b <- t(X.tmp)%*%v.tmp+solve(beta.var)%*%beta.mn
		# beta <- A%*%b+t(chol(A))%*%matrix(rnorm(qX),qX,1)
		# # beta <- t(rmvnorm(1,A%*%b,A))
		# gamma <- pnorm(X%*%beta)


		# ####
		# #### Sample parameters related to persistence (phi), i.e., z_t-1=1
		# ####
			
		# # Sample q (auxilliary variable for phi)
# q <- numeric(T)
		# z0.idx <- z[-T]==1&z[-1]==0 #idx for z_t-1=1 & z_t=0
		# z1.idx <- z[-T]==1&z[-1]==1 #idx for z_t-1=1 & z_t=1
		# # cbind(z,c(NA,z0.idx),c(NA,z1.idx))
		# q[z0.idx] <- truncnormsamp(matrix(X[z0.idx,],,qX)%*%beta,1,-Inf,0,sum(z0.idx))
		# q[z1.idx] <- truncnormsamp(matrix(X[z1.idx,],,qX)%*%beta,1,0,Inf,sum(z1.idx))		

		# # Sample theta, calculate phi
		# idx <- z0.idx|z1.idx
		# X.tmp <- X[idx,]
		# q.tmp <- q[idx]
		# A <- solve(t(X.tmp)%*%X.tmp+solve(theta.var))
		# b <- t(X.tmp)%*%q.tmp+solve(theta.var)%*%theta.mn
		# theta <- A%*%b+t(chol(A))%*%matrix(rnorm(qX),qX,1)
		# # phi <- t(rmvnorm(1,A%*%b,A))
		# phi <- pnorm(X%*%theta)


		# ####
		# #### Sample parameters related to detection probability (p)
		# ####
			
		# # Sample u (auxilliary variable for p)
# u <- numeric(T*dim(W)[3])	
		# z1.idx <- z[-1]==1
		# z1.y0.idx <- which(z1.idx&Y.vec==0)
		# z1.y1.idx <- which(z1.idx&Y.vec==1)
		# u[z1.y0.idx] <- truncnormsamp((matrix(W.mat[z1.y0.idx,],,qW)%*%alpha),1,-Inf,0,length(z1.y0.idx))
		# u[z1.y1.idx] <- truncnormsamp((matrix(W.mat[z1.y1.idx,],,qW)%*%alpha),1,0,Inf,length(z1.y1.idx))			
		# # Sample alpha
		# W.tmp <- W.mat[z1.idx,]		
		# u.tmp <- u[z1.idx]
		# A <- solve(t(W.tmp)%*%W.tmp+solve(alpha.var))
		# b <- u.tmp%*%W.tmp+t(alpha.mn)%*%solve(alpha.var)
		# alpha <- A%*%t(b)+t(chol(A))%*%matrix(rnorm(qW),qW,1)		
		# p <- matrix(apply(W.mat,1,function(x) pnorm(x%*%alpha)),,max(J))


		# ####
		# ####  Sample psi (beta) 
		# ####
# # browser()				
		# # A <- solve(t(X)%*%X+solve(beta.var))
		# # b <- t(X)%*%v+solve(beta.var)%*%beta.mn
		# # beta <- A%*%b+t(chol(A))%*%matrix(rnorm(qX),qX,1)
		# # # beta <- t(rmvnorm(1,A%*%b,A))
		# # psi <- pnorm(X%*%beta)
			

		# ####
		# ####  Sample z_0 (t=0) 
		# ####
# browser()
		# z[1] <- 0
		
		# ####
		# ####  Sample z_t when y_t==0, for t=1:(T-1)
		# ####
		
		# y0.idx <- which(y==0)

		# for(i in y0.idx){ #loop through primary sampling periods, t=1:(T-1)
		
			# i <- y0.idx[1]
		
			# z_tm1 <- z[i-1] #z_t-1
			# z_tp1 <- z[i+1] #z_t+1

			# p.tmp <- prod(1-p[i,],na.rm=TRUE)
			# phi.tmp <- phi[i]
			# gamma.tmp <- gamma[i]
			# A <- (p.tmp)*(phi.tmp^z_tm1)*(phi.tmp^z_tp1)*((1-phi.tmp)^(1-z_tp1))*(gamma.tmp^(1-z_tm1))
			# B <- ((1-phi.tmp)^z_tm1)*(gamma.tmp^z_tp1)*((1-gamma.tmp)^(1-z_tm1))*((1-gamma.tmp)^(1-z_tp1))
			# psi.tilde <- A/(A+B)
			
			# y0.idx <- which(y[-T]==0)
			# z[i] <- rbinom(length(y0.idx),1,psi.tilde[y0.idx])
		# }


		# ####
		# ####  Sample z_T (t=T)
		# ####

		# z[T+1] <- 0	
		
		
		# ####
		# ####  Save Samples 
		# ####
		
		# beta.save[k,] <- beta
		# alpha.save[k,] <- alpha
		# theta.save[k,] <- theta
		# z.mean <- z.mean+z/n.mcmc
		# # N.save[k] <- sum(z)
	
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
