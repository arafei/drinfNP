
################################################################
#                                                              #
#                        SIMULATION STUDY                      #
#                                                              #
################################################################
# Author: Ali Rafei (arafei {AT} umich {DOT} edu)
# Date:   02/01/2021

# The following R/stan codes are prepared for making robust Bayesian inference
# for non-probability samples. The function performs joint modeling of PS and
# outcome using linear-in-weight prediction (LWP) and Gaussian Process of PS 
# (GPPP) approaches.

# cal_bayes(psmp, nsmp, method=c("LWP", "GPPP"), spec=c(T, T), family="gaussian", nMCMC=1000L, nMC=200L, cores=1L)
# Predicts the outcome variable for units of the reference survey, and draw point and variance estimates using Rubin's combining rules
# Arguments:
# psmp: Probability sample with auxiliary variables and sampling weights as the last column
# nsmp: Nonprobability sample with auxiliary variables and the outcome variable as the last column
# method: calibrating method, either LWP or GPPP
# spec: specification status of the underlying models: c(T, T) means both PS and outcome models are correctly specified, respecively.
# family: determine the type of the outcome: currently available options are c("gaussian", "binomial", "nb")
# nMCMC: number of total MCMC draws. Note that an equal number of MCMC draws will be set aside as burn-in period.
# M: The size of the random sample selected from the MCMC draws to calculate the adjusted mean and variance
# cores: The number of cores is for parallel processing. It is also equivalent to the number of MCMC chains generated.
# Output of the function: Method, Mean, SE, and 95%CIs (LCL, UCL)

# Required packages:
library(MASS)
library(survey)
library(rstan)

logit <- function(x){
	log(x/(1-x))
}

# Weighted variance (TSL)
weighted.var <- function(y0, w0){
	n <- length(y0)
	x <- sum(y0 * w0)
	z <- sum(w0)
	(var(y0*w0) + ((x/z)^2)*var(w0) - 2*(x/z)*cov(y0*w0, w0))*n/(z^2)
}

#weighted variance conditional on \pi^B_i
cw.var <- function(y0, w0, id0=NULL){
	w0 <- w0[!is.na(y0)]
	y0 <- y0[!is.na(y0)]
	M <- 1
	if(!is.null(id0)){
		y0 <- as.numeric(rowsum(y0, id0))
		w0 <- w0[!duplicated(id0)]
		M <- mean(table(id0))
	}	
	sum(w0^2)*var(y0)/((M^2)*sum(w0)^2)
}

#Weighted variance conditional on y_i
cy.var <- function(y0, w0){
	n <- length(y0)
	x <- sum(y0 * w0)
	z <- sum(w0)
	y1 <- sum(y0)
	y2 <- sum(y0^2)
	(y2 + n*((x/z)^2) - 2*(x/z)*y1)*var(w0)/(z^2)
}

# Weighted se for clustered data
weighted.se <- function(y0, w0, id0=NULL){
	w0 <- w0[!is.na(y0)]
	y0 <- y0[!is.na(y0)]
	if(!is.null(id0)){
		wy <- as.numeric(rowsum(w0*y0, id0))
		w0 <- as.numeric(rowsum(w0, id0))
	}else{
		wy <- w0*y0
	}
	n <- length(w0)
	x <- sum(wy)
	z <- sum(w0)
	sqrt((var(wy) + ((x/z)^2)*var(w0) - 2*(x/z)*cov(wy, w0))*n/(z^2))
}

#Example of input for cal_bayes 
#psmp=smp0[, c("X", "X0", "XD", "Z0", "D", "D2", "wght", "Y")]; nsmp=smp1[, c("X", "X0", "XD", "Z0", "D", "D2", "wght", "Y")]; method=c("GPPP", "LWP"); spec=c(T, T); family="gaussian"; nMCMC=1000L; skip=500L; M=500; cores=1L; knot=10; V=3; C=1.25; hpar1=c(3, 3, 3, 3, 3, 3, 1, 2); hpar2=c(3, 3, 3, 3, 3)

#Set stan's options
rstan_options(auto_write = TRUE);
options(mc.cores = parallel::detectCores());

#call stan codes 
GP_C <- stan_model("GPPP_k_normal_C_1.stan")
#GP_B <- stan_model("GPPP_k_normal_B_1.stan")
#GP_NB <- stan_model("GPPP_k_normal_NB_1.stan")
LWP_C <- stan_model("LWP_normal_C_1.stan")
#LWP_B <- stan_model("LWP_normal_B_1.stan")
#LWP_NB <- stan_model("LWP_normal_NB_1.stan")
#GP_C <- stan_model("GPPP_k_beta_C_1.stan")
#GP_B <- stan_model("GPPP_k_beta_B.stan")
#GP_NB <- stan_model("GPPP_k_beta_NB.stan")
#LWP_C <- stan_model("LWP_beta_C_1.stan")
#LWP_B <- stan_model("LWP_beta_B.stan")
#LWP_NB <- stan_model("LWP_beta_NB.stan")

# SMP structure: X, W, Y, Z
maGPPP <- function(psmp, nsmp, N, method=c("LWP", "GPPP"), spec=c(T, T), family="gaussian", nMCMC=1000L, skip=500L, M=200L, cores=1L, knot=20, V=3, C=1.25, hpar1, hpar2){
	n0 <- nrow(psmp)
	n1 <- nrow(nsmp)
	ncov <- ncol(nsmp)-2
	y <- nsmp[, ncov+2]
	y0 <- psmp[, ncov+2]
	W0 <- psmp[, ncov+1]
	W1 <- nsmp[, ncov+1]
	W <- c(W1, W0)
	z <- rep(1:0, c(n1, n0))
	smp <- rbind(nsmp[, 1:ncov, drop=FALSE], psmp[, 1:ncov, drop=FALSE])
	ml <- sample(1:(nMCMC-skip), M, replace=F)
	chain1 <- cores
	if(spec[1]){
		PS_f <- c("X0")
	}else{
		PS_f <- c("Z0")
	}
	if(spec[2]){
		PM_f <- c("X", "D", "XD")
	}else{
		PM_f <- c("Z0", "D2")
	}
	P_f <- union(PS_f, PM_f)
	rwght <- Mu <- SE <- nms <- param <- gof <- LL <- UL <- c()
	if("GPPP"%in%method){
		if(family=="gaussian"){
			data <- list(Npop=N, N0=n0, N1=n1, Y=y, Z=z, X=as.matrix(smp[, P_f]), W=W0, ind1=1:n1, ind0=(n1+1):(n0+n1), Vp=match(PS_f, P_f), Vy=match(PM_f, P_f), Kp=length(PS_f), Ky=length(PM_f), K=length(P_f), V=V, M=knot, C=C, offset=NULL, sig=hpar1)
			fit2 <- sampling(GP_C, chains=cores, iter=nMCMC/cores, warmup=skip, seed=2312021, control = list(adapt_delta=0.95), data=data, cores=cores)
		}else if(family=="binomial"){
			data <- list(Npop=N, N0=n0, N1=n1, Y=y, Z=z, X=as.matrix(smp[, P_f]), W=W0, ind1=1:n1, ind0=(n1+1):(n0+n1), Vp=match(PS_f, P_f), Vy=match(PM_f, P_f), Kp=length(PS_f), Ky=length(PM_f), K=length(P_f), V=V, M=knot, C=C, offset=NULL, sig=hpar1)
			#fit2 <- sampling(GP_B, chains=cores, iter=nMCMC/cores, warmup=skip, seed=2312021, control = list(adapt_delta=0.95), data=data, cores=cores)
		}else if(family=="nb"){
			data <- list(Npop=N, N0=n0, N1=n1, Y=y, Z=z, X=as.matrix(smp[, P_f]), W=W0, ind1=1:n1, ind0=(n1+1):(n0+n1), Vp=match(PS_f, P_f), Vy=match(PM_f, P_f), Kp=length(PS_f), Ky=length(PM_f), K=length(P_f), V=V, M=knot, C=C, offset=log(10000+100*smp$D), sig=hpar1)
			#fit2 <- sampling(GP_NB, chains=cores, iter=nMCMC/cores, warmup=skip, seed=2312021, control = list(adapt_delta=0.95), data=data, cores=cores)
		}else{return("The family you specified does not exist. The currently available options are gaussian, binomial and nb")}
		fit20 <- extract(fit2)
		param <- c(param, mean(fit20$c0), mean(fit20$c1), mean(fit20$theta), mean(fit20$a0), mean(fit20$a1), mean(fit20$b0), apply(fit20$b1, 2, mean), apply(fit20$b2, 2, mean), mean(fit20$gp0), mean(fit20$sigma), mean(fit20$alpha), mean(fit20$rho), mean(fit20$L))		#
		yh0 <- t(fit20$y_predict[ml, ])
		gof <- cbind(gof, c(cor(y0, apply(yh0, 1, mean)), mean((y0 - apply(yh0, 1, mean))^2)))
		mu <- fit20$ybar
		Mu <- c(Mu, mean(mu, na.rm=T))
		SE <- c(SE,   sd(mu, na.rm=T))
		LL <- c(LL, quantile(mu, probs=0.025))
		UL <- c(UL, quantile(mu, probs=0.975))
		nms <- c(nms, "PAPP_GPPP")
	}
	if("LWP"%in%method){
		if(family=="gaussian"){
			data <- list(Npop=N, N0=n0, N1=n1, Y=y, Z=z, X=as.matrix(smp[, P_f]), W=W0, ind1=1:n1, ind0=(n1+1):(n0+n1), Vp=match(PS_f, P_f), Vy=match(PM_f, P_f), Kp=length(PS_f), Ky=length(PM_f), K=length(P_f), offset=NULL, sig=hpar2)
			fit2 <- sampling(LWP_C, chains=cores, iter=nMCMC/cores, warmup=skip, seed=2312021, control = list(adapt_delta=0.95), data=data, cores=cores)
		}else if(family=="binomial"){
			data <- list(Npop=N, N0=n0, N1=n1, Y=y, Z=z, X=as.matrix(smp[, P_f]), W=W0, ind1=1:n1, ind0=(n1+1):(n0+n1), Vp=match(PS_f, P_f), Vy=match(PM_f, P_f), Kp=length(PS_f), Ky=length(PM_f), K=length(P_f), offset=NULL, sig=hpar2)
			#fit2 <- sampling(LWP_B, chains=cores, iter=nMCMC/cores, warmup=skip, seed=2312021, control = list(adapt_delta=0.95), data=data, cores=cores)
		}else if(family=="nb"){
			data <- list(Npop=N, N0=n0, N1=n1, Y=y, Z=z, X=as.matrix(smp[, P_f]), W=W0, ind1=1:n1, ind0=(n1+1):(n0+n1), Vp=match(PS_f, P_f), Vy=match(PM_f, P_f), Kp=length(PS_f), Ky=length(PM_f), K=length(P_f), offset=log(10000+100*smp$D), sig=hpar2)
			#fit2 <- sampling(LWP_NB, chains=cores, iter=nMCMC/cores, warmup=skip, seed=2312021, control = list(adapt_delta=0.95), data=data, cores=cores)
		}else{return("The family you specified does not exist. The currently available options are gaussian, binomial and nb")}
		fit20 <- extract(fit2)
		param <- c(param, mean(fit20$c0), mean(fit20$c1), mean(fit20$theta), mean(fit20$a0), mean(fit20$a1), mean(fit20$b0), apply(fit20$b1, 2, mean), mean(fit20$b2), mean(fit20$sigma))
		yh0 <- t(fit20$y_predict[ml, ])
		gof <- cbind(gof, c(COR=cor(y0, apply(yh0, 1, mean)), MSE=mean((y0 - apply(yh0, 1, mean))^2)))
		mu <- fit20$ybar
		Mu <- c(Mu, mean(mu, na.rm=T))
		SE <- c(SE,   sd(mu, na.rm=T))
		LL <- c(LL, quantile(mu, probs=0.025))
		UL <- c(UL, quantile(mu, probs=0.975))
		nms <- c(nms, "PAPP_LPW")
	}
	colnames(gof) <- nms
	rownames(gof) <- c("COR", "MSE")
	output <- list(Method=nms, Mu=Mu, LCL=LL, UCL=UL, SE=SE, GOF=gof, PAR=param)
}
