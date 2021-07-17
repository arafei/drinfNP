
#rm(list=ls())
setwd("/home/arafei/MyJobs/Thesis/Chapter IV/Simulation/Bayesian/Own4/Gaussian/FB1")

library(survey)
library(MASS)
library(mvtnorm)
library(foreach)
library(doParallel)

source("maGPPP_func.R")

for(n0 in 500){
for(n1 in c(500, 1000)){
for(lambda in c(0.3, 0.6)){

set.seed(01012021)

V1=4; V2=1; rho=0.5; Cov=rho*sqrt(V1*V2); sigma=0.8; #lambda=0.6;
nSim <- 216
N <- 100000
#n1 <- 500
#n0 <- 1000
nMC <- 1000L

# generating covariates X and D in the population
cov <- matrix(c(V1, Cov, Cov, V2), 2, 2)
mu <- c(0, 0)

pop <- as.data.frame(mvrnorm(n = N, mu=mu, Sigma=cov, tol = 1e-6, empirical = FALSE, EISPACK = FALSE))
names(pop) <- c("X0", "D")

pop$X1 <- pop$X0
pop$X2 <- (pop$X0/3)^3
pop$X3 <- exp(pop$X0/2)/5
pop$X4 <- 5*sin(pi*pop$X0/3)
pop$Z0 <- pop$X0^2
pop$XD <- pop$X0*pop$D
pop$D2 <- pop$D^2
pop$e_y <- rnorm(N, 0, 1)

summary(pop)

#par(mfrow=c(2, 2))
#plot(pop$X0, pop$X1)
#plot(pop$X0, pop$X2)
#plot(pop$X0, pop$X3)
#plot(pop$X0, pop$X4)

rtcor <- function(x, e, r){
uniroot(function(t){cor(x+t*e, x) - r}, interval=c(0, 20))$root
}

s1 <- rtcor(3 + pop$X1 + pop$D + 0.2*pop$XD, pop$e_y, r=sigma)
pop$Y1 <- 3 + pop$X1 + pop$D + 0.2*pop$XD + s1*pop$e_y
s2 <- rtcor(3 + pop$X2 + pop$D + 0.2*pop$XD, pop$e_y, r=sigma)
pop$Y2 <- 3 + pop$X2 + pop$D + 0.2*pop$XD + s2*pop$e_y
s3 <- rtcor(3 + pop$X3 + pop$D + 0.2*pop$XD, pop$e_y, r=sigma)
pop$Y3 <- 3 + pop$X3 + pop$D + 0.2*pop$XD + s3*pop$e_y
s3 <- rtcor(3 + pop$X4 + pop$D + 0.2*pop$XD, pop$e_y, r=sigma)
pop$Y4 <- 3 + pop$X4 + pop$D + 0.2*pop$XD + s3*pop$e_y

summary(pop$Y1)
summary(pop$Y2)
summary(pop$Y3)
summary(pop$Y4)

cor(pop[, c("X1", "X2", "X3" ,"X4", "D2", "XD")])

cor(pop$Y1, 3 + pop$X1 + pop$D + 0.2*pop$XD)
cor(pop$Y2, 3 + pop$X2 + pop$D + 0.2*pop$XD)
cor(pop$Y3, 3 + pop$X3 + pop$D + 0.2*pop$XD)
cor(pop$Y4, 3 + pop$X4 + pop$D + 0.2*pop$XD)

rtexp0 <- function(x){
	uniroot(function(t){mean(exp(t+x)/(1+exp(t+x)))-n0/N}, interval=c(-20, 20))$root
}

px0 <- -0.4*pop$D
summary(px0)
gamma0 <- rtexp0(px0)
pop$pi0 <- exp(gamma0+px0)/(1+exp(gamma0+px0))
sum(pop$pi0)
summary(1/pop$pi0)

rtexp1 <- function(x){
	uniroot(function(t){mean(exp(t+x)/(1+exp(t+x)))-n1/N}, interval=c(-20, 20))$root
}

px1 <- lambda*pop$X0
summary(px1)

gamma1 <- rtexp1(px1)

pop$pi1 <- exp(gamma1+px1)/(1+exp(gamma1+px1))

sum(pop$pi1)

summary(1/pop$pi1)

cor(pop$Y1, px0)
cor(pop$Y2, px0)
cor(pop$Y3, px0)
cor(pop$Y4, px0)

cor(pop$Y1, px1)
cor(pop$Y2, px1)
cor(pop$Y3, px1)
cor(pop$Y4, px1)

#par(mfrow=c(2, 4))
#plot(log(pop$pi1), pop$Y1, pch=".", xlab="log(p_i)", ylab="y_1")
#plot(log(pop$pi1), pop$Y2, pch=".", xlab="log(p_i)", ylab="y_2")
#plot(log(pop$pi1), pop$Y3, pch=".", xlab="log(p_i)", ylab="y_3")
#plot(log(pop$pi1), pop$Y4, pch=".", xlab="log(p_i)", ylab="y_4")
#plot(pop$pi1, pop$Y1, pch=".")
#plot(pop$pi1, pop$Y2, pch=".")
#plot(pop$pi1, pop$Y3, pch=".")
#plot(pop$pi1, pop$Y4, pch=".")
#plot(1/pop$pi1, pop$Y1, pch=".", xlab="w_i", ylab="y_1")
#plot(1/pop$pi1, pop$Y2, pch=".", xlab="w_i", ylab="y_2")
#plot(1/pop$pi1, pop$Y3, pch=".", xlab="w_i", ylab="y_3")
#plot(1/pop$pi1, pop$Y4, pch=".", xlab="w_i", ylab="y_4")

pop$Y <- pop$Y2
pop$X <- pop$X2
mu1 <- mean(pop$Y1)
mu  <- mean(pop$Y2)
mu3 <- mean(pop$Y3)
mu4 <- mean(pop$Y4)

ZZ0 <- ZZ1 <- list()
for(i in 1:nSim){
	ZZ0[[i]] <- (1:N)[as.logical(rbinom(N, 1, pop$pi0))]
	ZZ1[[i]] <- (1:N)[as.logical(rbinom(N, 1, pop$pi1))]	
}

mean(sapply(ZZ0, length))
mean(sapply(ZZ1, length))

#####################################################

b <- 2
s <- 4*34
p <- 4+ 4*b

cores <- as.numeric(Sys.getenv('SLURM_NTASKS_PER_NODE', unset='8'))
cl <- makeCluster(cores)
registerDoParallel(cl)

combine <- function(x, ...) {  
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
}

start_time <- Sys.time()
res <- foreach(i=1:nSim, .combine=combine, .multicombine = TRUE, .packages = c("MASS", "survey", "rstan")) %dopar% {
	b <- 2
	s <- 4*34
	p <- 4+4*b
	res1 <- as.data.frame(matrix(NA, 1, 10*p + s))
	names(res1) <- c(paste0("mu", 1:p), paste0("se", 1:p), paste0("bias", 1:p), paste0("rmse", 1:p), paste0("ci_ll", 1:p), paste0("ci_ul", 1:p), paste0("cov", 1:p), paste0("cle", 1:p), paste0("cor", 1:p), paste0("mse", 1:p), paste0("par", 1:s))

	smp0 <- pop[ZZ0[[i]], ]
	smp1 <- pop[ZZ1[[i]], ]

	smp0$wght <- 1/smp0$pi0
	smp1$wght <- 1/smp1$pi0
	smp1$wght0 <- 1/smp1$pi1
	#smp0$wght <- N*smp0$wght/sum(smp0$wght)
	#smp1$wght0<- N*smp1$wght0/sum(smp1$wght0) 

	#start_time <- Sys.time()
	 out10 <- maGPPP(psmp=smp0[, c("X", "X0", "XD", "Z0", "D", "D2", "wght", "Y")], nsmp=smp1[, c("X", "X0", "XD", "Z0", "D", "D2", "wght", "Y")], N=N, method=c("GPPP", "LWP"), spec=c(T, T), family="gaussian", nMCMC=nMC, skip=nMC/2, M=nMC/2, cores=1L, knot=10, V=3, C=1.25, hpar1=c(3, 3, 3, 3, 3, 3, 1, 2), hpar2=c(3, 3, 3, 3, 3))
	 out11 <- maGPPP(psmp=smp0[, c("X", "X0", "XD", "Z0", "D", "D2", "wght", "Y")], nsmp=smp1[, c("X", "X0", "XD", "Z0", "D", "D2", "wght", "Y")], N=N, method=c("GPPP", "LWP"), spec=c(F, T), family="gaussian", nMCMC=nMC, skip=nMC/2, M=nMC/2, cores=1L, knot=10, V=3, C=1.25, hpar1=c(3, 3, 3, 3, 3, 3, 1, 2), hpar2=c(3, 3, 3, 3, 3))
	 out12 <- maGPPP(psmp=smp0[, c("X", "X0", "XD", "Z0", "D", "D2", "wght", "Y")], nsmp=smp1[, c("X", "X0", "XD", "Z0", "D", "D2", "wght", "Y")], N=N, method=c("GPPP", "LWP"), spec=c(T, F), family="gaussian", nMCMC=nMC, skip=nMC/2, M=nMC/2, cores=1L, knot=10, V=3, C=1.25, hpar1=c(3, 3, 3, 3, 3, 3, 1, 2), hpar2=c(3, 3, 3, 3, 3))
	 out13 <- maGPPP(psmp=smp0[, c("X", "X0", "XD", "Z0", "D", "D2", "wght", "Y")], nsmp=smp1[, c("X", "X0", "XD", "Z0", "D", "D2", "wght", "Y")], N=N, method=c("GPPP", "LWP"), spec=c(F, F), family="gaussian", nMCMC=nMC, skip=nMC/2, M=nMC/2, cores=1L, knot=10, V=3, C=1.25, hpar1=c(3, 3, 3, 3, 3, 3, 1, 2), hpar2=c(3, 3, 3, 3, 3))
	#out13 <- out10 <- out11 <- out12

	# Estimating unweighted mean and 95%CI	
	res1$mu1[1] <- mean(smp0$Y, na.rm=T)
	res1$se1[1] <- sd(smp0$Y, na.rm=T)/sqrt(n0)
	res1$ci_ll1[1] <- res1$mu1[1] - qnorm(0.975)*res1$se1[1]
	res1$ci_ul1[1] <- res1$mu1[1] + qnorm(0.975)*res1$se1[1]

	# Estimating fully weighted mean and 95%CI
	res1$mu2[1] <- weighted.mean(smp0$Y, w=smp0$wght, na.rm=T)
	res1$se2[1] <- weighted.se(smp0$Y, smp0$wght)
	res1$ci_ll2[1] <- res1$mu2[1] - qnorm(0.975)*res1$se2[1]
	res1$ci_ul2[1] <- res1$mu2[1] + qnorm(0.975)*res1$se2[1]

	# Estimating unweighted mean and 95%CI	
	res1$mu3[1] <- mean(smp1$Y, na.rm=T)
	res1$se3[1] <- sd(smp1$Y, na.rm=T)/sqrt(n1)
	res1$ci_ll3[1] <- res1$mu3[1] - qnorm(0.975)*res1$se3[1]
	res1$ci_ul3[1] <- res1$mu3[1] + qnorm(0.975)*res1$se3[1]

	# Estimating fully weighted mean and 95%CI
	res1$mu4[1] <- weighted.mean(smp1$Y, w=smp1$wght0, na.rm=T)
	res1$se4[1] <- weighted.se(smp1$Y, smp1$wght0)
	res1$ci_ll4[1] <- res1$mu4[1] - qnorm(0.975)*res1$se4[1]
	res1$ci_ul4[1] <- res1$mu4[1] + qnorm(0.975)*res1$se4[1]

	res1[1, paste0("mu", 5:p)] <- c(out10$Mu, out11$Mu, out12$Mu, out13$Mu)
	res1[1, paste0("se", 5:p)] <- c(out10$SE, out11$SE, out12$SE, out13$SE)
	res1[1, paste0("ci_ll", 5:p)] <- c(out10$LCL, out11$LCL, out12$LCL, out13$LCL)
	res1[1, paste0("ci_ul", 5:p)] <- c(out10$UCL, out11$UCL, out12$UCL, out13$UCL)
	res1[1, paste0("cor", 5:p)] <- c(out10$GOF[1, ], out11$GOF[1, ], out12$GOF[1, ], out13$GOF[1, ])
	res1[1, paste0("mse", 5:p)] <- c(out10$GOF[2, ], out11$GOF[2, ], out12$GOF[2, ], out13$GOF[2, ])
	res1[1, paste0("par", 1:s)] <- c(out10$PAR, out11$PAR, out12$PAR, out13$PAR)

	res1
}
end_time <- Sys.time()
end_time - start_time

stopCluster(cl)

res1 <- as.data.frame(res)

res1[, paste0("bias", 1:p)] <- res1[, paste0("mu", 1:p)]-mu
res1[, paste0("rmse", 1:p)] <- res1[, paste0("bias", 1:p)]^2
res1[, paste0("cov", 1:p)] <- res1[, paste0("ci_ll", 1:p)]<mu & mu<res1[, paste0("ci_ul", 1:p)]
res1[, paste0("cle", 1:p)] <- res1[, paste0("ci_ul", 1:p)] - res1[, paste0("ci_ll", 1:p)]

write.csv(res1, paste0("res1_TT_glm_BS_", n0, "_", n1, "_", rho, "_", lambda, "_2.csv"), col.names=T, sep=T)

tb <- as.data.frame(matrix(NA, p, 13))
names(tb) <- c("method", "bias", "rMSE", "cov_rate", "cov_len", "se_ratio", "MeanSE", "TrueSE", "Q05", "Q95", "biasSE", "cor", "mse")
tb$method <- c("UW_PS", "TW_PS", "UW_NS", "TW_NS", paste(c("GPPP", "LWP"), c("TT"), sep="_"), paste(c("GPPP", "LWP"), c("FT"), sep="_"), paste(c("GPPP", "LWP"), c("TF"), sep="_"), paste(c("GPPP", "LWP"), c("FF"), sep="_"))

tb[, 2] <- 100*apply(res1[, paste0("bias", 1:p)], 2, function(x)mean(x, na.rm=T))/mu
tb[, 3] <- 100*sqrt(apply(res1[, paste0("rmse", 1:p)], 2, function(x)mean(x, na.rm=T)))/mu
tb[, 4] <- 100*apply(res1[, paste0("cov", 1:p)], 2, function(x)mean(x, na.rm=T))
tb[, 5] <- apply(res1[, paste0("cle", 1:p)], 2, function(x)mean(x, na.rm=T))
tb[, 6] <- apply(res1[, paste0("se", 1:p)], 2, function(x)mean(x, na.rm=T))/apply(res1[, paste0("mu", 1:p)], 2, function(x)sd(x, na.rm=T))
tb[, 7] <- apply(res1[, paste0("se", 1:p)], 2, function(x)mean(x, na.rm=T))
tb[, 8] <- apply(res1[, paste0("mu", 1:p)], 2, function(x)sd(x, na.rm=T))
tb[, 9] <- 100*apply(res1[, paste0("bias", 1:p)], 2, function(x)quantile(x, probs=0.05, na.rm=T))/mu
tb[, 10] <- 100*apply(res1[, paste0("bias", 1:p)], 2, function(x)quantile(x, probs=0.95, na.rm=T))/mu
tb[, 11] <- 100*apply(res1[, paste0("mu", 1:p)], 2, function(x)sd(x, na.rm=T))/mu
tb[, 12] <- apply(res1[, paste0("cor", 1:p)], 2, mean)
tb[, 13] <- apply(res1[, paste0("mse", 1:p)], 2, mean)

tb

apply(res1[, paste0("par", 1:s)], 2 ,mean)

write.csv(tb, paste0("tb_TT_glm_BS_", n0, "_", n1, "_", rho, "_", lambda, "_2.csv"), col.names=T, sep=",")
}}}