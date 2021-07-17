
library(mgcv)
library(survey)
library(MASS)
library(polyapost)

logit <- function(x){
	log(x/(1-x))
}


#thetas = theta; design.x = cov.m; ICCy=0.2; a=1000; na=100; family="gaussian"
#thetas = theta; design.x = cov.m; ICCy=0.2; a=1000; na=100; family="gaussian"
getY <- function(theta, X, ICC, a, na, sigma=1, family="gaussian"){
	dim(theta) <- c(length(theta), 1) #make it matrix of ncol=1
	N <- nrow(X); 
	X <- as.matrix(cbind(1, X))
	if(family=="binomial"){
		p <- exp(X %*% theta)/(1+exp(X %*% theta)) #Pr(y=1|D,E)
		if (ICC==0){
			y <- rbinom(N, 1, p)
		}else{
			ei0 <- rep(rnorm(a), each=na)
			eij <- rnorm(N)
			Uij <- rbinom(N, 1, sqrt(ICC))
			thetaij <- qnorm(p); 
			threshold <- Uij*ei0 + (1-Uij)*eij
			y <- ifelse(threshold <= thetaij, 1, 0)
		}
	}else{
		mu <- X %*% theta
		if (ICC==0){
			y <- rnorm(N, mu, sigma)
		}else{
			eij <- rnorm(N, 0, sigma)
			su <- ICC*(sigma^2)/(1-ICC)
			ei0 <- rep(rnorm(a, 0, sqrt(su)), each=na)
			y <- mu + ei0 + eij
		}
	}
	return(as.numeric(y))
}



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


# Weighted variance
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

# Weighted variance
weighted.var <- function(y0, w0, id0=NULL){
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
	(var(wy) + ((x/z)^2)*var(w0) - 2*(x/z)*cov(wy, w0))*n/(z^2)
}

weighted.var1 <- function(y0, w0, id0){
	w0 <- w0[!is.na(y0)]
	y0 <- y0[!is.na(y0)]
	wy <- as.numeric(rowsum(w0*y0, id0))
	w <- as.numeric(rowsum(w0, id0))
	n <- length(w)
	x <- sum(wy)
	z <- sum(w)
	c(dx2=n*var(wy),
	     dz2=n*var(w),
	     dxz=n*cov(wy, w))
}

weighted.seh <- function(y0, w0, id0=NULL, st0){
	out <- sapply(unique(st0), FUN=function(h)weighted.var1(y0[st0==h], w0[st0==h], id0=id0[st0==h]))
	yw <- weighted.mean(y0, w0)
	sqrt((sum(out[1, ]) + (yw^2)*sum(out[2, ]) - 2*yw*sum(out[3, ]))/(sum(w0)^2))
}

#psmp=smp0[, c("X", "X0", "XD", "Z0", "D", "D2", "strat", "psu", "wght", "Y")]; nsmp=smp1[, c("X", "X0", "XD", "Z0", "D", "D2", "strat", "psu", "wght", "Y")]; method=c("GPPP", "PSPP", "LWP", "AIPW"); N=N; family="gaussian"; Bt1=20; Bt2=10; spec=c(T, F);

maGPPP <- function(psmp, nsmp, N, method=c("GPPP", "PSPP", "LWP", "APIW"), spec=c("T", "T"), family="gaussian", Bt1=5, Bt2=5){
	n0 <- nrow(psmp)
	n1 <- nrow(nsmp)
	n <- n0 + n1
	ncov <- ncol(nsmp)-2
	names(nsmp)[ncov+2] <- "Y"
	names(psmp)[ncov+2] <- "Y"
	psmp$W <- psmp[, ncov+1]
	nsmp$W <- nsmp[, ncov+1]
	if(spec[1]){
		PS_f <- c("X0")
	}else{
		PS_f <- c("Z0")
	}
	if(spec[2]){
		PM_f <- c("X", "D", "strat", "XD")
	}else{
		PM_f <- c("Z0")
	}
	P_f <- c(PM_f, setdiff(PS_f, PM_f), "W", "Y")
	FW_M1 <- GPPP_M1 <- PSPP_M1 <- LWP_M1 <- AIPW_M1 <- rep(NA, Bt1)
	nsmp_id <- cbind(1:n1, sapply(1:Bt1, function(x)sample(1:n1, n1, replace=T)))

	svyd <- svydesign(ids=~psu, strata=~strat, weights=~W, data=psmp, nest=T)
	svyd.r <- as.svrepdesign(design=svyd, type="subbootstrap", replicates=nBt1)
	repwt0 <- as.matrix(svyd.r$repweights)
	for(j in 1:Bt1){
		st.BB0 <- psmp					#[rep(1:n0, repwt0[, j]), ]
		st.BB0$W <- st.BB0$W*repwt0[, j]
		st.BB0 <- st.BB0[st.BB0$W!=0, ]
		#N <- sum(st.BB0$W)
		np0 <- nrow(st.BB0)
		ns <- N - np0
		Samwt0 <- N*st.BB0$W/sum(st.BB0$W)

		st.BB1 <- nsmp[nsmp_id[, j], ]
		st.BB <- rbind(st.BB1[, P_f], st.BB0[, P_f])
		st.BB$Z <- rep(1:0, c(n1, np0))
		FW_M2 <- GPPP_M2 <- PSPP_M2 <- LWP_M2 <- AIPW_M2 <- rep(NA, Bt2)
		for(boot in 1:Bt2){
			sm <- wtpolyap(1:np0, Samwt0, ns)
			st.BB$wf <- c(rep(NA, n1), table(factor(sm, levels=1:np0)))
			SPX <- rbind(st.BB1[, P_f], st.BB0[sm, P_f])
			SPX$Z <- rep(1:0, c(n1, N))
			fit1 <- glm(as.formula(paste("Z", paste(PS_f, collapse = " + "), sep = " ~ ")), data=SPX, family=binomial(link="logit"))
			st.BB$PS1_log <- predict(fit1, newdata=st.BB, type="link")
			st.BB$PS1 <- exp(st.BB$PS1_log)
			st.BB$PW1 <- 1/st.BB$PS1
			FW_M2[boot] <- weighted.mean(st.BB$Y[st.BB$Z==0], st.BB$wf[st.BB$Z==0], na.rm=T)
			if("GPPP"%in%method){
				bs <- "gp"; d2 <- c(3, -1, 1); knot=10;
				if(family=="binomial"){
					fit2 <- gam(as.formula(paste("Y", paste(c("s(PS1_log, bs=bs, m=d2, k=knot)", PM_f), collapse = " + "), sep = " ~ ")), family=binomial(link="logit"), data=st.BB[st.BB$Z==1, ])		#, method="REML"	
				}else if(family=="gaussian"){
					fit2 <- gam(as.formula(paste("Y", paste(c("s(PS1_log, bs=bs, m=d2, k=knot)", PM_f), collapse = " + "), sep = " ~ ")), family=gaussian(), data=st.BB[st.BB$Z==1, ])				#, method="REML"
				}
				Yh <- predict(fit2, newdata=st.BB[st.BB$Z==0, ], exclude=NULL, type="response")
				GPPP_M2[boot] <- weighted.mean(Yh, st.BB$wf[st.BB$Z==0], na.rm=T) + mean(st.BB$Y[st.BB$Z==1] - fit2$fitted.values, na.rm=T)
			}
			if("PSPP"%in%method){
				bs <- "ps"; d2 <- c(2, 0); knot=10;
				if(family=="binomial"){
					fit2 <- gam(as.formula(paste("Y", paste(c("s(PS1_log, bs=bs, m=d2, k=knot)", PM_f), collapse = " + "), sep = " ~ ")), family=binomial(link="logit"), data=st.BB[st.BB$Z==1, ])		#, method="REML"
				}else if(family=="gaussian"){
					fit2 <- gam(as.formula(paste("Y", paste(c("s(PS1_log, bs=bs, m=d2, k=knot)", PM_f), collapse = " + "), sep = " ~ ")), family=gaussian(), data=st.BB[st.BB$Z==1, ])				#, method="REML"
				}
				Yh <- predict(fit2, newdata=st.BB[st.BB$Z==0, ], exclude=NULL, type="response")
				PSPP_M2[boot] <- weighted.mean(Yh, st.BB$wf[st.BB$Z==0], na.rm=T) + mean(st.BB$Y[st.BB$Z==1] - fit2$fitted.values, na.rm=T)
			}
			if("LWP"%in%method){
				if(family=="binomial"){
					fit2 <- glm(as.formula(paste("Y", paste(c("PW1", PM_f), collapse = " + "), sep = " ~ ")), family=binomial(link="logit"), data=st.BB[st.BB$Z==1, ])
				}else if(family=="gaussian"){
					fit2 <- lm(as.formula(paste("Y", paste(c("PW1", PM_f), collapse = " + "), sep = " ~ ")), data=st.BB[st.BB$Z==1, ])
				}
				Yh <- predict(fit2, newdata=st.BB[st.BB$Z==0, ], exclude=NULL, type="response")
				LWP_M2[boot] <- weighted.mean(Yh, st.BB$wf[st.BB$Z==0], na.rm=T) + mean(st.BB$Y[st.BB$Z==1] - fit2$fitted.values, na.rm=T)
			}
			if("AIPW"%in%method){
				if(family=="binomial"){
					fit2 <- glm(as.formula(paste("Y", paste(PM_f, collapse = " + "), sep = " ~ ")), family=binomial(link=logit), data=st.BB[st.BB$Z==1, ])
				}else if(family=="gaussian"){		
					fit2 <- lm(as.formula(paste("Y", paste(PM_f, collapse = " + "), sep = " ~ ")), data=st.BB[st.BB$Z==1, ])
				}
				Yh <- predict(fit2, newdata=st.BB[st.BB$Z==0, ], type="response")
				AIPW_M2[boot] <- weighted.mean(Yh, st.BB$wf[st.BB$Z==0], na.rm=T) + weighted.mean(st.BB$Y[st.BB$Z==1] - fit2$fitted.values, w=st.BB$PW1[st.BB$Z==1])
			}
		}
		FW_M1[j] <- mean(FW_M2)
		GPPP_M1[j] <- mean(GPPP_M2)
		PSPP_M1[j] <- mean(PSPP_M2)
		LWP_M1[j]  <- mean(LWP_M2)
		AIPW_M1[j] <- mean(AIPW_M2)
	}
	Mu <- c(mean(GPPP_M1), mean(PSPP_M1), mean(LWP_M1), mean(AIPW_M1))
	SE <- sqrt((1 + 1/Bt1)*c(var(GPPP_M1), var(PSPP_M1), var(LWP_M1), var(AIPW_M1)))

	LL <- Mu - qnorm(0.975)*SE
	UL <- Mu + qnorm(0.975)*SE
	list(Method=c("GPPP", "PSPP", "LPW", "AIPW"), Mu=Mu, LCL=LL, UCL=UL, SE=SE, FW_Mu=mean(FW_M1), FW_SE=sqrt((1 + 1/Bt1)*var(GPPP_M1)))
}
