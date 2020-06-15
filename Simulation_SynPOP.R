

rm(list=ls())
setwd("C:\\Users\\arafei\\Desktop\\paper3_PSPP_sim")

library(pps)
library(survey)
library(betareg)
library(MASS)

source("func_FPBB_PSPP3N-n.R")

set.seed(06042020)

rho=0.8; sigma=1;
nSim <- 1000
N <- 10000
n0 <- 200
n1 <- 100
nBt1 <- 20
nBt2 <- 5

# generating covariates X_1 and X_2 in the population
#cov <- matrix(c(1, rho, rho, 1), 2, 2)
#mu <- c(1, 0)

#pop <- as.data.frame(mvrnorm(n = N, mu=mu, Sigma=cov, tol = 1e-6, empirical = FALSE, EISPACK = FALSE))
#names(pop) <- c("X1", "D")
#pop$X2 <- rbinom(N, 1, 0.5)
#pop$X12 <- pop$X1^2
#pop$X13 <- pop$X1^3
#pop$X1X2 <- pop$X1*pop$X2

# generating errors
#pop$e_y <- rnorm(N, 0, 1)

#sigma <- 1
#pop$Y1 <- 2+0.4*pop$X12+0.2*pop$X13-pop$X2-0.1*pop$X1X2-pop$D+sigma*pop$e_y
#pop$Y2 <- rbinom(N, 1, exp(-1+0.2*pop$X12+0.1*pop$X13-pop$X2-0.5*pop$X1X2-pop$D)/(1+exp(-1+0.2*pop$X12+0.1*pop$X13-pop$X2-0.5*pop$X1X2-pop$D)))

#cor(pop$Y1, 2+0.4*pop$X12+0.2*pop$X13-pop$X2-0.1*pop$X1X2-pop$D)
#summary(pop$Y1)
#summary(pop$Y2)

#rtexp0 <- function(x){
#uniroot(function(t){mean(exp(t+x)/(1+exp(t+x)))-n0/N}, interval=c(-20, 20))$root
#}

#rtexp1 <- function(x){
#uniroot(function(t){mean(exp(t+x)/(1+exp(t+x)))-n1/N}, interval=c(-20, 20))$root
#}

#gamma0 <- rtexp0(0.5*pop$D)
#gamma1 <- rtexp1(+0.4*pop$X1-0.2*pop$X12+0.6*pop$X2+0.1*pop$X1X2)

#pop$pi0 <- exp(gamma0+0.5*pop$D)/(1+exp(gamma0+0.5*pop$D))
#pop$pi1 <- exp(gamma1+0.4*pop$X1-0.2*pop$X12+0.6*pop$X2+0.1*pop$X1X2)/(1+exp(gamma1+0.4*pop$X1-0.2*pop$X12+0.6*pop$X2+0.1*pop$X1X2))
#summary(pop$pi0)
#summary(pop$pi1)

#cor(pop$Y1, 0.5*pop$D)
#cor(pop$Y1, -0.1*pop$X12+0.2*pop$X2+0.3*pop$X1X2)

#write.csv(pop, "sim_pop.csv", col.names=T, sep=",")
pop <- read.csv("sim_pop.csv", header=T, sep=",")
pop$X <- NULL

mu1 <- mean(pop$Y1)
mu2 <- mean(pop$Y2)

#ZZ0 <- list()
#ZZ1 <- list()
#for(i in 1:nSim){
#	ZZ0[[i]] <- (1:N)[as.logical(rbinom(N, 1, pop$pi0))]
#	ZZ1[[i]] <- (1:N)[as.logical(rbinom(N, 1, pop$pi1))]
#}

#saveRDS(ZZ0, file="smp_idx0.rds")
#saveRDS(ZZ1, file="smp_idx1.rds")
ZZ0 <- readRDS(file="smp_idx0.rds")
ZZ1 <- readRDS(file="smp_idx1.rds")

b <- 2
p <- 4*b+4
res1 <- as.data.frame(matrix(NA, nSim, 7*p))
res2 <- as.data.frame(matrix(NA, nSim, 7*p))

names(res1) <- c(paste0("mu", 1:p), paste0("se", 1:p), paste0("bias", 1:p), paste0("rmse", 1:p), paste0("ci_ll", 1:p), paste0("ci_ul", 1:p), paste0("cov", 1:p))
names(res2) <- c(paste0("mu", 1:p), paste0("se", 1:p), paste0("bias", 1:p), paste0("rmse", 1:p), paste0("ci_ll", 1:p), paste0("ci_ul", 1:p), paste0("cov", 1:p))

for(i in 566:nSim){
	smp0 <- pop[ZZ0[[i]], ]
	smp1 <- pop[ZZ1[[i]], ]

	smp0$wght <- 1/smp0$pi0
	smp1$wght <- 1/smp1$pi1
	smp1$wght0 <- 1/smp1$pi0

	out10 <- FPBB(psmp=smp0[, c("X1", "X2", "X12", "X13", "X1X2", "D", "wght")], nsmp=smp1[, c("X1", "X2", "X12", "X13", "X1X2", "D", "Y1")], k=1, method=c("PSPP", "AIPW"), model="GLM", spec=c(T, T), N=10000, Bt1=nBt1, Bt2=nBt2, nknots=20, nMCMC=200, ntree=50, nskip=100, cores=10)
	out11 <- FPBB(psmp=smp0[, c("X1", "X2", "X12", "X13", "X1X2", "D", "wght")], nsmp=smp1[, c("X1", "X2", "X12", "X13", "X1X2", "D", "Y1")], k=1, method=c("PSPP", "AIPW"), model="GLM", spec=c(T, F), N=10000, Bt1=nBt1, Bt2=nBt2, nknots=20, nMCMC=200, ntree=50, nskip=100, cores=10)
	out12 <- FPBB(psmp=smp0[, c("X1", "X2", "X12", "X13", "X1X2", "D", "wght")], nsmp=smp1[, c("X1", "X2", "X12", "X13", "X1X2", "D", "Y1")], k=1, method=c("PSPP", "AIPW"), model="GLM", spec=c(F, T), N=10000, Bt1=nBt1, Bt2=nBt2, nknots=20, nMCMC=200, ntree=50, nskip=100, cores=10)
	out13 <- FPBB(psmp=smp0[, c("X1", "X2", "X12", "X13", "X1X2", "D", "wght")], nsmp=smp1[, c("X1", "X2", "X12", "X13", "X1X2", "D", "Y1")], k=1, method=c("PSPP", "AIPW"), model="GLM", spec=c(F, F), N=10000, Bt1=nBt1, Bt2=nBt2, nknots=20, nMCMC=200, ntree=50, nskip=100, cores=10)
	#out14 <- FPBB(psmp=smp0[, c("X1", "X2", "X12", "X13", "X1X2", "D", "wght")], nsmp=smp1[, c("X1", "X2", "X12", "X13", "X1X2", "D", "Y1")], k=1, method=c("PSPP", "AIPW"), model="BART", spec=c(F, F), N=10000, Bt1=nBt1, Bt2=nBt2, nknots=20, nMCMC=200, ntree=50, nskip=100, cores=10)
	#out20 <- FPBB(psmp=smp0[, c("X1", "X2", "X12", "X13", "X1X2", "D", "wght")], nsmp=smp1[, c("X1", "X2", "X12", "X13", "X1X2", "D", "Y2")], k=1, method=c("PSPP", "AIPW"), model="GLM", spec=c(T, T), N=10000, Bt1=nBt1, Bt2=nBt2, nknots=20, nMCMC=200, ntree=50, nskip=100, cores=10)
	#out21 <- FPBB(psmp=smp0[, c("X1", "X2", "X12", "X13", "X1X2", "D", "wght")], nsmp=smp1[, c("X1", "X2", "X12", "X13", "X1X2", "D", "Y2")], k=1, method=c("PSPP", "AIPW"), model="GLM", spec=c(T, F), N=10000, Bt1=nBt1, Bt2=nBt2, nknots=20, nMCMC=200, ntree=50, nskip=100, cores=10)
	#out22 <- FPBB(psmp=smp0[, c("X1", "X2", "X12", "X13", "X1X2", "D", "wght")], nsmp=smp1[, c("X1", "X2", "X12", "X13", "X1X2", "D", "Y2")], k=1, method=c("PSPP", "AIPW"), model="GLM", spec=c(F, T), N=10000, Bt1=nBt1, Bt2=nBt2, nknots=20, nMCMC=200, ntree=50, nskip=100, cores=10)
	#out23 <- FPBB(psmp=smp0[, c("X1", "X2", "X12", "X13", "X1X2", "D", "wght")], nsmp=smp1[, c("X1", "X2", "X12", "X13", "X1X2", "D", "Y2")], k=1, method=c("PSPP", "AIPW"), model="GLM", spec=c(F, F), N=10000, Bt1=nBt1, Bt2=nBt2, nknots=20, nMCMC=200, ntree=50, nskip=100, cores=10)
	#out24 <- FPBB(psmp=smp0[, c("X1", "X2", "X12", "X13", "X1X2", "D", "wght")], nsmp=smp1[, c("X1", "X2", "X12", "X13", "X1X2", "D", "Y2")], k=1, method=c("PSPP", "AIPW"), model="BART", spec=c(F, F), N=10000, Bt1=nBt1, Bt2=nBt2, nknots=20, nMCMC=200, ntree=50, nskip=100, cores=10)

	# Estimating unweighted mean and 95%CI	
	svyd_uw0 <- svydesign(ids=~1, weights=~1, data=smp0)
	res1$mu1[i] <- svymean(~Y1, design=svyd_uw0)[1]
	res1$se1[i] <- SE(svymean(~Y1, design=svyd_uw0))[1]
	res1$ci_ll1[i] <- confint(svymean(~Y1, design=svyd_uw0), level=0.95)[1]
	res1$ci_ul1[i] <- confint(svymean(~Y1, design=svyd_uw0), level=0.95)[2]

	# Estimating fully weighted mean and 95%CI
	svyd_tw0 <- svydesign(ids=~1, weights=~wght, data=smp0)
	res1$mu2[i] <- svymean(~Y1, design=svyd_tw0)[1]
	res1$se2[i] <- SE(svymean(~Y1, design=svyd_tw0))
	res1$ci_ll2[i] <- confint(svymean(~Y1, design=svyd_tw0), level=0.95)[1]
	res1$ci_ul2[i] <- confint(svymean(~Y1, design=svyd_tw0), level=0.95)[2]

	# Estimating unweighted mean and 95%CI	
	svyd_uw1 <- svydesign(ids=~1, weights=~1, data=smp1)
	res1$mu3[i] <- svymean(~Y1, design=svyd_uw1)[1]
	res1$se3[i] <- SE(svymean(~Y1, design=svyd_uw1))[1]
	res1$ci_ll3[i] <- confint(svymean(~Y1, design=svyd_uw1), level=0.95)[1]
	res1$ci_ul3[i] <- confint(svymean(~Y1, design=svyd_uw1), level=0.95)[2]

	# Estimating fully weighted mean and 95%CI
	svyd_tw1 <- svydesign(ids=~1, weights=~wght, data=smp1)
	res1$mu4[i] <- svymean(~Y1, design=svyd_tw1)[1]
	res1$se4[i] <- SE(svymean(~Y1, design=svyd_tw1))
	res1$ci_ll4[i] <- confint(svymean(~Y1, design=svyd_tw1), level=0.95)[1]
	res1$ci_ul4[i] <- confint(svymean(~Y1, design=svyd_tw1), level=0.95)[2]

	res1[i, paste0("mu", 5:p)] <- c(out10$Mu, out11$Mu, out12$Mu, out13$Mu)
	res1[i, paste0("se", 5:p)] <- c(out10$SE, out11$SE, out12$SE, out13$SE)

	# Estimating unweighted mean and 95%CI	
	res2$mu1[i] <- svymean(~Y2, design=svyd_uw0)[1]
	res2$se1[i] <- SE(svymean(~Y2, design=svyd_uw0))[1]
	res2$ci_ll1[i] <- confint(svymean(~Y2, design=svyd_uw0), level=0.95)[1]
	res2$ci_ul1[i] <- confint(svymean(~Y2, design=svyd_uw0), level=0.95)[2]

	# Estimating fully weighted mean and 95%CI
	res2$mu2[i] <- svymean(~Y2, design=svyd_tw0)[1]
	res2$se2[i] <- SE(svymean(~Y2, design=svyd_tw0))
	res2$ci_ll2[i] <- confint(svymean(~Y2, design=svyd_tw0), level=0.95)[1]
	res2$ci_ul2[i] <- confint(svymean(~Y2, design=svyd_tw0), level=0.95)[2]

	# Estimating unweighted mean and 95%CI	
	res2$mu3[i] <- svymean(~Y2, design=svyd_uw1)[1]
	res2$se3[i] <- SE(svymean(~Y2, design=svyd_uw1))[1]
	res2$ci_ll3[i] <- confint(svymean(~Y2, design=svyd_uw1), level=0.95)[1]
	res2$ci_ul3[i] <- confint(svymean(~Y2, design=svyd_uw1), level=0.95)[2]

	# Estimating fully weighted mean and 95%CI
	res2$mu4[i] <- svymean(~Y2, design=svyd_tw1)[1]
	res2$se4[i] <- SE(svymean(~Y2, design=svyd_tw1))
	res2$ci_ll4[i] <- confint(svymean(~Y2, design=svyd_tw1), level=0.95)[1]
	res2$ci_ul4[i] <- confint(svymean(~Y2, design=svyd_tw1), level=0.95)[2]

	#res2[i, paste0("mu", 5:p)] <- c(out20$Mu, out21$Mu, out22$Mu, out23$Mu)
	#res2[i, paste0("se", 5:p)] <- c(out20$SE, out21$SE, out22$SE, out23$SE)
}


res1[, paste0("bias", 1:p)] <- res1[, paste0("mu", 1:p)]-mu1
res1[, paste0("rmse", 1:p)] <- res1[, paste0("bias", 1:p)]^2
res1[, paste0("ci_ll", 5:p)] <- res1[, paste0("mu", 5:p)]-qt(0.975, n1-10)*res1[, paste0("se", 5:p)]
res1[, paste0("ci_ul", 5:p)] <- res1[, paste0("mu", 5:p)]+qt(0.975, n1-10)*res1[, paste0("se", 5:p)]
res1[, paste0("cov", 1:p)] <- res1[, paste0("ci_ll", 1:p)]<mu1 & mu1<res1[, paste0("ci_ul", 1:p)]

res2[, paste0("bias", 1:p)] <- res2[, paste0("mu", 1:p)]-mu2
res2[, paste0("rmse", 1:p)] <- res2[, paste0("bias", 1:p)]^2
res2[, paste0("ci_ll", 5:p)] <- res2[, paste0("mu", 5:p)]-qt(0.975, n1-10)*res2[, paste0("se", 5:p)]
res2[, paste0("ci_ul", 5:p)] <- res2[, paste0("mu", 5:p)]+qt(0.975, n1-10)*res2[, paste0("se", 5:p)]
res2[, paste0("cov", 1:p)] <- res2[, paste0("ci_ll", 1:p)]<mu2 & mu2<res2[, paste0("ci_ul", 1:p)]

write.csv(res1, "res1_TT_glm_JK_d2_logitp.csv", col.names=T, sep=",")
write.csv(res2, "res2_TT_glm_JK_d2_logitp.csv", col.names=T, sep=",")

tb <- as.data.frame(matrix(NA, p, 19))
names(tb) <- c("method", "bias_y1", "rMSE_y1", "cov_rate_y1", "se_ratio_y1", "MeanSE_y1", "TrueSE_y1", "Q5_y1", "Q95_y1", "biasSE_y1", "bias_y2", "rMSE_y2", "cov_rate_y2", "se_ratio_y2", "MeanSE_y2", "TrueSE_y2", "Q5_y2", "Q95_y2", "biasSE_y2")
tb$method <- c("UW_PS", "TW_PS", "UW_NS", "TW_NS", paste(out10$Model, c("TT"), sep="_"), paste(out10$Model, c("TF"), sep="_"), paste(out10$Model, c("FT"), sep="_"), paste(out10$Model, c("FF"), sep="_"))

tb[, 2] <- 100*apply(res1[, paste0("bias", 1:p)], 2, function(x)mean(x, na.rm=T))/mu1
tb[, 3] <- sqrt(apply(res1[, paste0("rmse", 1:p)], 2, function(x)mean(x, na.rm=T)))
tb[, 4] <- 100*apply(res1[, paste0("cov", 1:p)], 2, function(x)mean(x, na.rm=T))
tb[, 5] <- apply(res1[, paste0("se", 1:p)], 2, function(x)mean(x, na.rm=T))/apply(res1[, paste0("mu", 1:p)], 2, function(x)sd(x, na.rm=T))
tb[, 6] <- apply(res1[, paste0("se", 1:p)], 2, function(x)mean(x, na.rm=T))
tb[, 7] <- apply(res1[, paste0("mu", 1:p)], 2, function(x)sd(x, na.rm=T))
tb[, 8] <- 100*apply(res1[, paste0("bias", 1:p)], 2, function(x)quantile(x, probs=0.05, na.rm=T))/mu1
tb[, 9] <- 100*apply(res1[, paste0("bias", 1:p)], 2, function(x)quantile(x, probs=0.95, na.rm=T))/mu1
tb[, 10] <- 100*apply(res1[, paste0("mu", 1:p)], 2, function(x)sd(x, na.rm=T))/mu1

tb[, 11] <- 100*apply(res2[, paste0("bias", 1:p)], 2, function(x)mean(x, na.rm=T))/mu2
tb[, 12] <- sqrt(apply(res2[, paste0("rmse", 1:p)], 2, function(x)mean(x, na.rm=T)))
tb[, 13] <- 100*apply(res2[, paste0("cov", 1:p)], 2, function(x)mean(x, na.rm=T))
tb[, 14] <- apply(res2[, paste0("se", 1:p)], 2, function(x)mean(x, na.rm=T))/apply(res2[, paste0("mu", 1:p)], 2, function(x)sd(x, na.rm=T))
tb[, 15] <- apply(res2[, paste0("se", 1:p)], 2, function(x)mean(x, na.rm=T))
tb[, 16] <- apply(res2[, paste0("mu", 1:p)], 2, function(x)sd(x, na.rm=T))
tb[, 17] <- 100*apply(res2[, paste0("bias", 1:p)], 2, function(x)quantile(x, probs=0.05, na.rm=T))/mu2
tb[, 18] <- 100*apply(res2[, paste0("bias", 1:p)], 2, function(x)quantile(x, probs=0.95, na.rm=T))/mu2
tb[, 19] <- 100*apply(res2[, paste0("mu", 1:p)], 2, function(x)sd(x, na.rm=T))/mu2

tb

write.csv(tb, "tb_TT_glm_JK_d2_logitp.csv", col.names=T, sep=",")
