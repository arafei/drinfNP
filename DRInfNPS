
################################################################
#                                                              #
#     Doubly Robust Inference for Non-probability Samples      #
#                                                              #
################################################################
# Author: Ali Rafei (arafei {AT} umich {DOT} edu)
# Date:   01/14/2018

# The following codes are prepared to evaluate the performance of
# quasi-randomization approaches including pseudo-weighting and 
# predictive mean matching. Several models including GLM, BART,
# CART and Random Forest are utilized.

# calibrate(psmp, nsmp, method=c("PSDO", "PMM"), model=c("GLM1", "GLM2", "BART", "CART", "RF"), q=0.01, JRR=FALSE, JKn=c(100, 100))
# Generates a set of pseudo-weights for the non-probability sample and a set of replicate pseudo-weights based on a modified Jackknife method
# Arguments:
# psmp: Probability sample with auxiliary variables and weights in the last column
# nsmp: Nonprobability sample with auxiliary variables
# method: calibrating method, either PSDO or PMM
# model: The model should be used for calibration: GLM, BART, CART or Random Forest
# q: The number of subclasses in PMM: 0.01 means to use percentiles of the distribution
# JRR: logical value to determine if Jackknife replicate method is required
# JKn0: favorite number of replicates for the probability sample
# JKn1: favorite number of replicates for the nonprobability sample


# Required packages:
library(MASS)
library(survey)
library(mgcv)
library(mvtnorm)
require(rpart)
library(kernlab)
library(BART)
library(BayesTree)
#library(gpe)
library(rstanarm)

# Function for trimming sampling weights and replication weights
# The upper threshold if defined for each replication weight independently



# Function for trimming sampling weights and replication weights
# The upper threshold if defined for each replication weight independently


trimw <- function(w, m='t1', c=10, max=Inf){
	if(m=='t1'){
		kn <- sqrt(c*mean(w^2, na.rm=T))
		i <- 0
		while(sum(w>kn & !is.na(w))>=1 & i<max){
			s <- sum(w, na.rm=T)
			w[w>kn & !is.na(w)] <- kn
			w[!is.na(w)] <- w[!is.na(w)]*s/sum(w, na.rm=T)
			kn <- sqrt(c*sum(w^2, na.rm=T)/sum(!is.na(w), na.rm=T))
			i <- i+1
		}
	} else if(m=='t2'){
		kn <- median(w, na.rm=T)+c*IQR(w, na.rm=T)
		i <- 0
		while(sum(w>kn, na.rm=T)>=1 & i<max){
			s <- sum(w, na.rm=T)
			w[w>kn & !is.na(w)] <- kn
			w[!is.na(w)] <- w[!is.na(w)]*s/sum(w, na.rm=T)
			kn <- median(w, na.rm=T)+c*IQR(w, na.rm=T)
			i <- i+1
		}
	} else if(m=='t3'){
		a <- 2+mean(w, na.rm=T)*(sum(w, na.rm=T)-1)/((length(w[!is.na(w)])-1)*var(w, na.rm=T))
		b <- (sum(w, na.rm=T)-1)*(1+mean(w, na.rm=T)*(sum(w, na.rm=T)-1)/((length(w[!is.na(w)])-1)*var(w, na.rm=T)))
		kn <- 1/(length(w[!is.na(w)])*qbeta(0.99, a, b))
		i <- 0
		while(sum(w>kn, na.rm=T)>=1 & i<max){
			s <- sum(w, na.rm=T)
			w[w>kn & !is.na(w)] <- kn
			w[!is.na(w)] <- w[!is.na(w)]*s/sum(w, na.rm=T)
			a <- 2+mean(w, na.rm=T)*(sum(w, na.rm=T)-1)/((length(w[!is.na(w)])-1)*var(w, na.rm=T))
			b <- (sum(w, na.rm=T)-1)*(1+mean(w, na.rm=T)*(sum(w, na.rm=T)-1)/((length(w[!is.na(w)])-1)*var(w, na.rm=T)))
			kn <- 1/(length(w[!is.na(w)])*qbeta(0.99, a, b))
			i <- i+1
		}
	}
	w	
}


#BS GLM: psmp=smp0[, c("X1", "X2", "X12", "X13", "D", "X1X2", "wght")]; nsmp=smp1[, c("X1", "X2", "X12", "X13", "D", "X1X2", "Y1", "wght0")]; method=c("PSDO0", "PSDO1", "IPSW", "IPSWC", "OR", "AIPW", "PSPP"); spec=c(T, T); model="GLM"; Npop=N; BS=TRUE; nBS=nBS; trim=c("no")
#psmp=smp0[, c("X1", "X12", "X13", "X2", "X1X2", "D", "wght")]; nsmp=smp1[, c("X1", "X12", "X13", "X2", "X1X2", "D", "Y1", "wght0")]; method=c("PSDO0", "PSDO1", "OR", "AIPW"); model="GLM"; spec=c(T, T); Npop=N; m0=nmc1[, i]; m1=nmc2[, i]; m2=nmc3[, i]; nMCMC=1000; trim=c("no"); C=10; MAX=1
#psmp1=psmp; nsmp1=nsmp; method1=method; model1=model; Npop1=Npop; nMCMC1=1000; trim=c("no"); trim1=trim;
#X=smp; Y=z; di=c(rep(1, n1), W); threshold = 1e-10; max_iter = 100
#X=pop$X; Y=z1; di=rep(1, N-n1); threshold = 1e-10; max_iter = 100


chen_glm <- function(X, Y, di, threshold = 1e-10, max_iter = 100){
	X <- as.matrix(X)
	n <- nrow(X)
	np <- ncol(X) + 1
	X <- cbind(rep(1, n), X)
	calc_p <- function(X, beta){
		beta <- as.vector(beta)
		return(1 / (1+ exp(-X%*%beta)))
	}  
	beta <- rep(0, np)
	diff <- 10000 
	iter_count <- 0
	while(diff > threshold){
		pi <- as.vector(calc_p(X, beta))
		U <- t(X)%*%Y-t(di[Y==0]*X[Y==0, ])%*%pi[Y==0]
		WW <- pi[Y==0]*(1-pi[Y==0])*di[Y==0]
		H <- t(X[Y==0, ])%*%(WW*X[Y==0,])
		beta_change <- ginv(H) %*% U
		beta <- beta + beta_change
		diff <- sum(beta_change^2)
		iter_count <- iter_count + 1
		if(iter_count > max_iter) {
			stop("This isn't converging, mate.")
		}
	} 
	list(coef=as.numeric(beta), predict=as.numeric(calc_p(X, beta)))
}

wglm <- function(X, Y, di, threshold = 1e-10, max_iter = 100){
	X <- as.matrix(X)
	n <- nrow(X)
	np <- ncol(X) + 1
	X <- cbind(rep(1, n), X)
	calc_p <- function(X, beta){
		beta <- as.vector(beta)
		return(1 / (1+ exp(-X%*%beta)))
	}
	beta <- rep(0, np)
	diff <- 10000 
	iter_count <- 0
	while(diff > threshold){
		pi <- as.vector(calc_p(X, beta))
		U <- t(di*X)%*%(Y - pi)
		H <- t(X)%*%(di*pi*(1-pi)*X)
		beta_change = ginv(H) %*% U
		beta <- beta + beta_change
		diff <- sum(beta_change^2)
		iter_count <- iter_count + 1
		if(iter_count > max_iter) {
			stop("This isn't converging, mate.")
		}
	} 
	list(coef=as.numeric(beta), predict=as.numeric(calc_p(X, beta)))
}


predict1.gp <- function (object, newdata = NULL, nMCMC=1000, ...){

    if (is.null(newdata)) {
        newdata <- object$data$training_data
    }
    projector <- getProjector(object$posterior)

    res <- projector(object$posterior, newdata, variance = "matrix")
    res$ppd <- mvrnorm(n=nMCMC, mu=res$mu, Sigma=res$var)

    link <- object$likelihood$link
    fit <- link(res$mu)
    ppd <- apply(res$ppd, c(1, 2), link)

    ans <- list(fit = fit, ppd=ppd)

    return(ans)
}
#environment(predict1.gp) <- environment(gp)


# SMP structure: X, W, Y, Z
calsmp <- function(psmp1, nsmp1, method1=c("PSDO0", "PSDO1", "IPS", "IPSWC", "IPSW", "PSPP"), model1=c("GLM", "CART", "GP"), PS_f, PM_f, nMCMC1, Npop1=100000, trim1=c("no", "t1", "t2"), C1=10, MAX1=10){
	n0 <- nrow(psmp1)
	n1 <- nrow(nsmp1)
	ncov <- ncol(nsmp1)-2
	y <- nsmp1[, ncov+1]
	ob <- length(unique(y))==2
	W <- c(nsmp1[, ncov+2], psmp1[, ncov+1])
	z <- rep(1:0, c(n1, n0))
	smp <- rbind(nsmp1[, 1:ncov, drop=FALSE], psmp1[, 1:ncov, drop=FALSE])
	X_smp <- makeind(smp)
	Wi_log <- log((1/W[z==0])/(1-1/W[z==0]))
	weight <- c()
	rwght <- c()
	out_m <- c()
	out_v <- c()
	if("GLM"%in%model1){
		if("PSDO0"%in%method1){
			fit2 <- glm(as.formula(paste("z", paste(c(PS_f, "D"), collapse = " + "), sep = " ~ ")), data=smp, family=binomial(link="logit"))
			PZi <- predict(fit2, type="response")
			weight <- cbind(weight, W*(1-PZi)/PZi)
		}
		if("PSDO1"%in%method1){
			fit1 <- lm(Wi_log~X1+X2, data=smp[z==0, , drop=FALSE])
			Wt <- predict(fit1, newdata=smp[, , drop=FALSE])
			Wt <- (1+exp(Wt))/exp(Wt)
			fit2 <- glm(as.formula(paste("z", paste(PS_f, collapse = " + "), sep = " ~ ")), data=smp, family=binomial(link="logit"))
			PZi <- predict(fit2, type="response")
			weight <- cbind(weight, Wt*(1-PZi)/PZi)
		}
		if("IPS"%in%method1){
			#fit2 <- glm(as.formula(paste("z", paste(PS_f, collapse = " + "), sep = " ~ ")), data=smp, family=binomial(link="logit"))
			#PZi <- predict(fit2, type="response")
			weight <- cbind(weight, 1/PZi)
		}
		if("IPSW"%in%method1){
			PZi <- wglm(X=smp[, PS_f], Y=z, di=c(rep(1, n1), W[z==0]))$predict
			weight <- cbind(weight, 1/PZi)
		}
		if("IPSWC"%in%method1){
			PZi <- chen_glm(X=smp[, PS_f], Y=z, di=c(rep(1, n1), W[z==0]))$predict
			weight <- cbind(weight, 1/PZi)
		}
	}else if("BART"%in%model1){
		if("PSDO0"%in%method1){
			fit2 <- pbart(x.train=X_smp[, c(PS_f, "D")], y.train=z, ntree=50L, nskip=100L, ndpost=nMCMC1)
			PZi <- fit2$prob.train.mean
			weight <- cbind(weight, W*(1-PZi)/PZi)
		}
		if("PSDO1"%in%method1){
			fit1 <- wbart(x.train=X_smp[z==0, c("X1", "X2"), drop=F], y.train=Wi_log, x.test=X_smp[, c("X1", "X2"), drop=F], ntree=50L, nskip=100L, ndpost=nMCMC1)
			Wt <- fit1$yhat.test.mean
			Wt <- (1+exp(Wt))/exp(Wt)
			fit2 <- pbart(x.train=X_smp[, PS_f, drop=F], y.train=z, ntree=50L, nskip=100L, ndpost=nMCMC1)
			PZi <- fit2$prob.train.mean
			weight <- cbind(weight, Wt*(1-PZi)/PZi)
		}
		if("IPS"%in%method1){
			#fit2 <- glm(as.formula(paste("z", paste(PS_f, collapse = " + "), sep = " ~ ")), data=smp, family=binomial(link="logit"))
			#PZi <- predict(fit2, type="response")
			weight <- cbind(weight, 1/PZi)
		}
		if("IPSW"%in%method1){
			PZi <- wglm(X=smp[, PS_f], Y=z, di=c(rep(1, n1), W[z==0]))$predict
			weight <- cbind(weight, 1/PZi)
		}
		if("IPSWC"%in%method1){
			PZi <- chen_glm(X=smp[, PS_f], Y=z, di=c(rep(1, n1), W[z==0]))$predict
			weight <- cbind(weight, 1/PZi)
		}
	}else if("CART"%in%model1){
		if("PSDO0"%in%method1){
			fit2 <- rpart(as.formula(paste("factor(z)", paste(c(PS_f, "D"), collapse = " + "), sep = " ~ ")), smp, control=rpart.control(minbucket=50, cp=0.01, xval=0),  method="class")
			PZi <- as.numeric(predict(fit2, type="prob")[, 2])
			weight <- cbind(weight, W*(1-PZi)/PZi)
		}
		if("PSDO1"%in%method1){
			fit1 <- rpart(Wi_log~X1+X2, data=smp[z==0, , drop=FALSE], control=rpart.control(minbucket=50, cp=0.01,   xval=0), method="anova")
			Wt <- predict(fit1, newdata=smp[, , drop=FALSE], type="vector")
			Wt <- (1+exp(Wt))/exp(Wt)
			fit2 <- rpart(as.formula(paste("factor(z)", paste(PS_f, collapse = " + "), sep = " ~ ")), smp, control=rpart.control(minbucket=50, cp=0.01, xval=0),  method="class")
			PZi <- as.numeric(predict(fit2, type="prob")[, 2])
			weight <- cbind(weight, Wt*(1-PZi)/PZi)
		}
	}else if("GP"%in%model1){
		if("PSDO0"%in%method1){
			fit2 <- gausspr(as.formula(paste("factor(z)", paste(c(PS_f, "D"), collapse = " + "), sep = " ~ ")), data=smp, type="classification")
			PZi <- as.numeric(predict(fit2, newdata=smp[z==1, PS_f, drop=FALSE], type="probabilities")[, 2])
			weight <- cbind(weight, W*(1-PZi)/PZi)
		}
		if("PSDO1"%in%method1){
			fit1 <- gausspr(Wi_log~X1+X2, data=smp[z==0, , drop=FALSE], type="regression")
			Wt <- predict(fit1, newdata=smp[, , drop=FALSE], type="response")
			Wt <- (1+exp(Wt))/exp(Wt)
			fit2 <- gausspr(as.formula(paste("factor(z)", paste(PS_f, collapse = " + "), sep = " ~ ")), data=smp, type="classification")
			PZi <- as.numeric(predict(fit2, newdata=smp[z==1, PS_f, drop=FALSE], type="probabilities")[, 2])
			weight <- cbind(weight, Wt*(1-PZi)/PZi)
		}
	}
	# normalizing weights
	weight <- t(N*t(weight)/apply(weight, 2, sum))
	# triming weights
	if("PSDO0"%in%method1 | "PSDO1"%in%method1 | "IPS"%in%method1 | "IPSW"%in%method1 | "IPSWC"%in%method1){
		if("no"%in%trim1){
			rwght <- cbind(rwght, weight)
		}
		if("t1"%in%trim1){
			rwght_t <- apply(weight, 2, function(x)trimw(x, m="t1", c=C1, max=MAX1))
			rwght <- cbind(rwght, rwght_t)
		}
		if("t2"%in%trim1){
			rwght_t <- apply(weight, 2, function(x)trimw(x, m="t2", c=C1, max=MAX1))
			rwght <- cbind(rwght, rwght_t)
		}
		#rwght <- t(Npop1*t(rwght)/apply(rwght, 2, function(x)sum(x, na.rm=T)))
		out_m <- c(out_m, apply(rwght[z==1, ], 2, function(x)weighted.mean(y, w=x, na.rm=T)))
		out_v <- c(out_v, apply(rwght[z==1, ], 2, function(x)weighted.var(y0=y, w0=x)))
	}
	if("OR"%in%method1 | "AIPW"%in%method1){
		if("GLM"%in%model1){
			if(ob){
				fit3 <- glm(as.formula(paste("y", paste(PM_f, collapse = " + "), sep = " ~ ")), data=smp[z==1, , drop=FALSE], family=binomial(link="logit"))
				yhi <- predict(fit3, newdata=smp[, PM_f, drop=FALSE], type="response")	
			}else{
				fit3 <- lm(as.formula(paste("y", paste(PM_f, collapse = " + "), sep = " ~ ")), data=smp[z==1, , drop=FALSE])
				yhi <- predict(fit3, newdata=smp[, PM_f, drop=FALSE])			
			}
		}else if("BART"%in%model1){
			if(ob){
				fit3 <- pbart(x.train=X_smp[z==1, PM_f, drop=F], y.train=y, x.test=X_smp[, PM_f, drop=F], ntree=50L, nskip=100L, ndpost=nMCMC1)
				yhi <- fit3$prob.test.mean
			}else{
				fit3 <- wbart(x.train=X_smp[z==1, PM_f, drop=F], y.train=y, x.test=X_smp[, PM_f, drop=F], ntree=50L, nskip=100L, ndpost=nMCMC1)
				yhi <- fit3$yhat.test.mean
			}
		}else if("CART"%in%model1){
			if(ob){
				fit3 <- rpart(as.formula(paste("y", paste(PM_f, collapse = " + "), sep = " ~ ")), data=smp[z==1, , drop=FALSE], control=rpart.control(minbucket=50, cp=0.01,   xval=0), method="anova")
				yhi <- predict(fit3, newdata=smp[, PM_f, drop=FALSE], type="vector")
			}else{
				fit3 <- rpart(as.formula(paste("factor(y)", paste(PM_f, collapse = " + "), sep = " ~ ")), data=smp[z==1, , drop=FALSE], control=rpart.control(minbucket=50, cp=0.01, xval=0),  method="class")
				yhi <- as.numeric(predict(fit3, newdata=smp[, PM_f, drop=FALSE], type="prob")[, 2])	
			}
		}else if("GP"%in%model1){
			if(ob){
				fit3 <- gausspr(as.formula(paste("factor(y)", paste(PM_f, collapse = " + "), sep = " ~ ")), data=smp[z==1, , drop=FALSE], type="classification")
				yhi <- as.numeric(predict(fit3, newdata=smp[, PM_f, drop=FALSE], type="probabilities")[, 2])			
			}else{
				fit3 <- gausspr(as.formula(paste("y", paste(PM_f, collapse = " + "), sep = " ~ ")), data=smp[z==1, , drop=FALSE], type="regression")
				yhi <- predict(fit3, newdata=smp[, PM_f, drop=FALSE], type="response")			
			}
		}
		out_m <- c(out_m, weighted.mean(yhi[z==0], w=W[z==0], na.rm=T))
		out_v <- c(out_v, weighted.var(y0=yhi[z==0], w0=W[z==0]))
	}
	if("AIPW"%in%method1){
		ei <- y-yhi[z==1]
		out_m <- c(out_m, apply(rwght[z==1, ], 2, function(x){weighted.mean(ei, w=x, na.rm=T) + weighted.mean(yhi[z==0], w=W[z==0], na.rm=T)}))
		out_v <- c(out_v, apply(rwght[z==1, ], 2, function(x){weighted.var(y0=ei, w0=x) + weighted.var(y0=yhi[z==0], w0=W[z==0])}))
	}
	if("PSPP"%in%method1){
		if("GLM"%in%model1){
			bs <- "re"; nknots <- 20; d2 <- 2;
			if(ob){
				yhi <- apply(rwght, 2, function(xx){fit4 <- gam(as.formula(paste("y", paste(c("s(PS, bs=bs, k=nknots, m=d2)", PM_f), collapse = " + "), sep = " ~ ")), family=binomial(link="logit"), data=data.frame(smp[z==1, ], PS=log((1/xx[z==1])/(1-1/xx[z==1])))); predict(fit4, newdata=data.frame(smp[z==0, ], PS=log((1/xx[z==0])/(1-1/xx[z==0]))), type="response")})
			}else{
				yhi <- apply(rwght, 2, function(xx){fit4 <- gam(as.formula(paste("y", paste(c("s(PS, bs=bs, k=nknots, m=d2)", PM_f), collapse = " + "), sep = " ~ ")), family=gaussian(), data=data.frame(smp[z==1, ], PS=log((1/xx[z==1])/(1-1/xx[z==1])))); predict(fit4, newdata=data.frame(smp[z==0, ], PS=log((1/xx[z==0])/(1-1/xx[z==0]))), type="response")})
			}	
		}else if("BART"%in%model1){
			if(ob){
				yhi <- apply(rwght, 2, function(xx){X_smp1 <- data.frame(X_smp, PS=xx); fit4 <- pbart(x.train=X_smp1[z==1, c(PM_f, "PS"), drop=F], y.train=y, x.test=X_smp1[z==0, c(PM_f, "PS"), drop=F], ntree=50L, nskip=100L, ndpost=nMCMC1); fit4$prob.test.mean})
			}else{
				yhi <- apply(rwght, 2, function(xx){X_smp1 <- data.frame(X_smp, PS=xx); fit4 <- wbart(x.train=X_smp1[z==1, c(PM_f, "PS"), drop=F], y.train=y, x.test=X_smp1[z==0, c(PM_f, "PS"), drop=F], ntree=50L, nskip=100L, ndpost=nMCMC1); fit4$yhat.test.mean})
			}
		}
		out_m <- c(out_m, apply(yhi, 2, function(xx)weighted.mean(xx, w=W[z==0], na.rm=T)))		
		out_v <- c(out_v, apply(yhi, 2, function(x)weighted.var(y0=x, w0=W[z==0])))
	}
	list(mu=out_m, var=out_v, RW=rwght[z==1, ])
}


cal_jk <- function(psmp, nsmp, method=c("PSDO0", "PSDO1", "IPS", "IPSW", "PSPP", "IPSWC"), model=c("GLM", "BART", "CART", "GP"), spec=c(T, T), Npop=10000, nMCMC=1000, JK=FALSE, JKn=c(100, 100), trim=c("no", "t1", "t2"), C=NULL, MAX=NULL){
	n0 <- nrow(psmp)
	n1 <- nrow(nsmp)
	if(spec[1]){
		PS_f0 <- c("X1", "X12", "X2", "X1X2")
	}else{
		PS_f0 <- c("X1", "X2")
	}
	if(spec[2]){
		PM_f0 <- c("X12", "X13", "X2", "X1X2", "D")
	}else{
		PM_f0 <- c("X1", "X2")
	}
	OT <- calsmp(psmp1=psmp, nsmp1=nsmp, method1=method, model1=model, PS_f=PS_f0, PM_f=PM_f0, nMCMC1=nMCMC, Npop1=Npop, trim1=trim, C1=C, MAX1=MAX)
	Mu <- OT$mu
	SE <- sqrt(OT$var)
	LL <- Mu - qnorm(0.975)*SE
	UL <- Mu + qnorm(0.975)*SE
	rwght <- OT$RW
	nms <- c(method[1:5], paste0("AIPW_", method[1:4]), paste0("PSPP_", method[1:4]))
	if("BART"%in%model){
		nms <- c(method[1:3], paste0("AIPW_", method[1:2]), paste0("PSPP_", method[1:4]))
	}
	#nms <- paste(rep(nms, length(trim)), rep(trim, rep(length(nms), length(trim))), sep="_")
	colnames(rwght) <- nms[1:4]
	output <- list(Model=nms, Mu=Mu, SE=SE, LCL=LL, UCL=UL, RW=rwght)
	if(JK){
		jk_id <- c(rep(1:JKn[1], rep(n0%/%JKn[1], JKn[1])), rep(JKn[1], n0%%JKn[1]))
		M1_i <- sapply(1:JKn[1], function(j)calsmp(psmp1=psmp[jk_id!=j, ], nsmp1=nsmp, method1=method, model1=model, PS_f=PS_f0, PM_f=PM_f0, Npop1=Npop, nMCMC1=nMCMC, trim1=trim, C1=C, MAX1=MAX)$mu)
		jk_id <- c(rep((JKn[1]+1):(JKn[1]+JKn[2]), rep(n1%/%JKn[2], JKn[2])), rep((JKn[1]+JKn[2]), n1%%JKn [2]))	
		M2_i <- sapply((JKn[1]+1):(JKn[1]+JKn[2]), function(j)calsmp(psmp1=psmp, nsmp1=nsmp[jk_id!=j, ], method1=method, model1=model, PS_f=PS_f0, PM_f=PM_f0, Npop1=Npop, nMCMC1=nMCMC, trim1=trim, C1=C, MAX1=MAX)$mu)
		M_i <- cbind(M1_i, M2_i)
		SE <- apply(M_i, 1, function(x)sqrt(sum((x-mean(x, na.rm=T))^2, na.rm=T)))
		LL <- Mu - qnorm(0.975)*SE
		UL <- Mu + qnorm(0.975)*SE
		output <- list(Model=nms, Mu=Mu, LCL=LL, UCL=UL, SE=SE, RW=rwght)
	}
	output
}


cal_bs <- function(psmp, nsmp, method=c("PSDO0", "PSDO1", "IPS", "IPSW", "PSPP", "IPSWC"), model=c("GLM", "CART", "BART", "GP"), spec=c(T, T), Npop=10000, BS=FALSE, nBS=c(100, 100), nMCMC=1000, trim=c("no", "t1", "t2"), C=NULL, MAX=NULL){
	n0 <- nrow(psmp)
	n1 <- nrow(nsmp)
	if(spec[1]){
		PS_f0 <- c("X1", "X12", "X2", "X1X2")
	}else{
		PS_f0 <- c("X1", "X2")
	}
	if(spec[2]){
		PM_f0 <- c("X12", "X13", "X2", "X1X2", "D")
	}else{
		PM_f0 <- c("X1", "X2")
	}
	OT <- calsmp(psmp1=psmp, nsmp1=nsmp, method1=method, model1=model, PS_f=PS_f0, PM_f=PM_f0, nMCMC1=nMCMC, Npop1=Npop, trim1=trim, C1=C, MAX1=MAX)
	Mu <- OT$mu
	SE <- sqrt(OT$var)
	LL <- Mu - qnorm(0.975)*SE
	UL <- Mu + qnorm(0.975)*SE
	rwght <- OT$RW
	nms <- c(method[1:3], paste0("AIPW_", method[1:2]), paste0("PSPP_", method[1:2]))
	#if("BART"%in%model){
	#	nms <- c(method[1:3], paste0("AIPW_", method[1:2]), paste0("PSPP_", method[1:2]))
	#}
	#nms <- paste(rep(nms, length(trim)), rep(trim, rep(length(nms), length(trim))), sep="_")
	colnames(rwght) <- nms[1:2]
	output <- list(Model=nms, Mu=Mu, SE=SE, LCL=LL, UCL=UL, RW=rwght)
	if(BS){
		svyd <- svydesign(ids=~1, strata=NULL, weights=~wght, data=psmp, nest=F)
		svyd.r <- as.svrepdesign(design=svyd, type="subbootstrap", replicates=nBS)
		repwt0 <- as.matrix(svyd.r$repweights)
		repwt0[repwt0==0] <- NA
		nsmp.id <- sapply(1:nBS, FUN=function(i)sample(1:n1, n1, replace=T))
		mu <- va <- matrix(NA, length(nms), nBS)
		for(j in 1:nBS){
			nsmp_bs <- nsmp[nsmp.id[, j], ]
			psmp_bs <- psmp
			psmp_bs$repwt0 <- repwt0[, j]
			psmp_bs <- na.omit(psmp_bs)
			psmp_bs <- psmp_bs[rep(1:nrow(psmp_bs), round(psmp_bs$repwt0)), ]
			psmp_bs$repwt0 <- NULL
			psmp_bs[, ncol(psmp_bs)] <- psmp_bs[, ncol(psmp_bs)]*n0/(n0-1)
			fit_bs <- calsmp(psmp1=psmp_bs, nsmp1=nsmp_bs, method1=method, model1=model, PS_f=PS_f0, PM_f=PM_f0, Npop1=Npop, nMCMC1=nMCMC, trim1=trim, C1=C, MAX1=MAX)
			mu[, j] <- fit_bs$mu
			va[, j] <- fit_bs$var
		}
		Mu <- apply(mu, 1, mean)
		SE <- sqrt(apply(mu, 1, var)) #(1+1/nBS)* + apply(va, 1, mean)
		LL <- Mu - qnorm(0.975)*SE
		UL <- Mu + qnorm(0.975)*SE
		output <- list(Model=nms, Mu=Mu, LCL=LL, UCL=UL, SE=SE, RW=rwght)
	}
	output
}

