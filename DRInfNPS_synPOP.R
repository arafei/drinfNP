
######################################################################################
#
#     Doubly Robust Inference for Non-probability Samples Using a 
#     Synthetic Population based on Finite Population Bayesian Bootstrapping 
#
######################################################################################


library(nlme)
library(mgcv)
library(polyapost)
library(survey)
library(BART)
library(BayesTree)

FPBB <- function(psmp, nsmp, k=1, method=c("PSPP", "APIW"), model=c("GLM", "BART"), spec=c("T", "T"), N=100000, Bt1=5, Bt2=20, nknots=20, nMCMC=1000, ntree=50, nskip=100, cores=10){
	n0 <- nrow(psmp)
	n1 <- nrow(nsmp)
	n <- n0 + n1
	ncov <- ncol(nsmp)-1
	names(nsmp)[ncov+1] <- "Y"
	ob <- length(unique(nsmp$Y))==2
	psmp$W <- (N/k-n1)*psmp[, ncov+1]/sum(psmp[, ncov+1])
	if(spec[1]){
		PS_f <- c("X1", "X12", "X2", "X1X2")
	}else{
		PS_f <- c("X1", "X2")
	}
	if(spec[2]){
		PM_f <- c("X12", "X13", "X2", "X1X2", "D")
	}else{
		PM_f <- c("X1", "X2")
	}
	svyd <- svydesign(ids=~1, strata=NULL, weights=~W, data=psmp, nest=F)
	svyd.r <- as.svrepdesign(design=svyd, type="bootstrap", replicates=Bt1)
	repwt0 <- as.matrix(svyd.r$repweights)
	repwt0[repwt0==0] <- NA
	svyd <- svydesign(ids=~1, strata=NULL, weights=~1, data=nsmp, nest=F)
	svyd.r <- as.svrepdesign(design=svyd, type="bootstrap", replicates=Bt1)
	repwt1 <- as.matrix(svyd.r$repweights)
	PSPP_M1 <- AIPW_M1 <- PSPP_V1 <- AIPW_V1 <- rep(NA, Bt1)
	for(j in 1:Bt1){
		st.bb0 <- cbind(psmp, repwt0[, j])
		st.BB0 <- na.omit(st.bb0)
		st.BB1 <- nsmp[rep(1:n1, repwt1[, j]), ]
		np0 <- nrow(st.BB0)
		Samwt0 <- st.BB0$W*st.BB0[, "repwt0[, j]"]
		Samwt0 <- Samwt0*(N/k-n1)/sum(Samwt0)
		ns <- N/k-n1-np0
		ids <- 1:np0
		PSPP_M2 <- AIPW_M2 <- rep(NA, Bt2)
		for(boot in 1: Bt2){
			sm <- wtpolyap(ids, Samwt0, ns)
			l <- sapply(1:np0, function(i)length(sm[sm==i]))
			SPX <- rbind(st.BB1[, 1:ncov], st.BB0[rep(ids, l), 1:ncov])
			SPX$Z <- c(rep(1, n1), rep(0, N/k-n1))
			SPX$Y <- NA; SPX$Y[SPX$Z==1] <- st.BB1$Y
			if(model=="GLM"){
				fit1 <- glm(as.formula(paste("Z", paste(PS_f, collapse = " + "), sep = " ~ ")), data=SPX, family=binomial(link="logit"))
				SPX$PS_logit <- log(fit1$fitted.values/(1-fit1$fitted.values))
				if("PSPP"%in%method){
					bs <- "re"; d2 <- 2;
					if(ob){
						fit2 <- gam(as.formula(paste("Y", paste(c("s(PS_logit, bs=bs, k=nknots, m=d2)", PM_f), collapse = " + "), sep = " ~ ")), family=binomial(link="logit"), data=SPX[SPX$Z==1, ],method="REML")
					}else{
						fit2 <- gam(as.formula(paste("Y", paste(c("s(PS_logit, bs=bs, k=nknots, m=d2)", PM_f), collapse = " + "), sep = " ~ ")), data=SPX[SPX$Z==1, ],method="REML")
					}
					SPX$Y[SPX$Z==0] <- predict(fit2, newdata=SPX[SPX$Z==0, ], type="response")
					PSPP_M2[boot] <- mean(SPX$Y)
				}
				if("AIPW"%in%method){
					if(ob){
						fit2 <- glm(as.formula(paste("Y", paste(PM_f, collapse = " + "), sep = " ~ ")), family=binomial(link=logit), data=SPX[SPX$Z==1, ])
					}else{		
						fit2 <- lm(as.formula(paste("Y", paste(PM_f, collapse = " + "), sep = " ~ ")), data=SPX[SPX$Z==1, ])
					}
					SPX$Y[SPX$Z==0] <- predict(fit2, newdata=SPX[SPX$Z==0, ], type="response")
					pwt <- 1/fit1$fitted.values[SPX$Z==1]
					AIPW_M2[boot] <- weighted.mean(fit2$residuals, w=pwt) + mean(SPX$Y)
				}
			}else if(model=="BART"){
				fit1 <- pbart(x.train=makeind(SPX[, PS_f, drop=F]), y.train=SPX$Z, ntree=ntree, nskip=nskip, ndpost=nMCMC)
				SPX$PS_logit <- log(fit1$prob.train.mean/(1-fit1$prob.train.mean))
				if("PSPP"%in%method){
					if(ob){
						fit2 <- pbart(x.train=makeind(SPX[SPX$Z==1, c("PS_logit", PM_f), drop=F]), y.train=SPX$Y[SPX$Z==1], x.test=makeind(SPX[SPX$Z==0, c("PS_logit", PM_f), drop=F]), ntree=50L, nskip=100L, ndpost=nMCMC)
						SPX$Y[SPX$Z==0] <- fit2$prob.test.mean
					}else{
						fit2 <- wbart(x.train=makeind(SPX[SPX$Z==1, c("PS_logit", PM_f), drop=F]), y.train=SPX$Y[SPX$Z==1], x.test=makeind(SPX[SPX$Z==0, c("PS_logit", PM_f), drop=F]), ntree=50L, nskip=100L, ndpost=nMCMC)
						SPX$Y[SPX$Z==0] <- fit2$yhat.test.mean
					}
					PSPP_M2[boot] <- mean(SPX$Y)
				}
				if("AIPW"%in%method){
					if(ob){
						fit2 <- pbart(x.train=makeind(SPX[SPX$Z==1, PM_f, drop=F]), y.train=SPX$Y[SPX$Z==1], x.test=makeind(SPX[SPX$Z==0, PM_f, drop=F]), ntree=ntree, nskip=nskip, ndpost=nMCMC)
						SPX$Y[SPX$Z==0] <- fit2$prob.test.mean
						resid <- SPX$Y[SPX$Z==1] - fit2$prob.train.mean
					}else{		
						fit2 <- wbart(x.train=makeind(SPX[SPX$Z==1, PM_f, drop=F]), y.train=SPX$Y[SPX$Z==1], x.test=makeind(SPX[SPX$Z==0, PM_f, drop=F]), ntree=ntree, nskip=nskip, ndpost=nMCMC)
						SPX$Y[SPX$Z==0] <- fit2$yhat.test.mean
						resid <- SPX$Y[SPX$Z==1] - fit2$yhat.train.mean
					}
					pwt <- 1/fit1$prob.train.mean[SPX$Z==1]
					AIPW_M2[boot] <- weighted.mean(resid, w=pwt) + mean(SPX$Y)
				}
			}
		}
		PSPP_M1[j] <- mean(PSPP_M2)
		AIPW_M1[j] <- mean(AIPW_M2)
		PSPP_V1[j] <- var(PSPP_M2)
		AIPW_V1[j] <- var(AIPW_M2)
	}
	Mu <- c(mean(PSPP_M1), mean(AIPW_M1))
	SE <- sqrt(c(var(PSPP_M1), var(AIPW_M1)))  #mean(PSPP_V1) + (1+1/Bt1)*    mean(AIPW_V1) + (1+1/Bt1)*

	LL <- Mu - qnorm(0.975)*SE
	UL <- Mu + qnorm(0.975)*SE
	list(Model=c("PSPP", "AIPW"), Mu=Mu, LCL=LL, UCL=UL, SE=SE)
}
