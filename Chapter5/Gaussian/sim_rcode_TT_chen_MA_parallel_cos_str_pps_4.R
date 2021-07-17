
#rm(list=ls())
setwd("/home/arafei/MyJobs/Thesis/Chapter V/Simulation/BS/Own4/Gaussian/Design 2stage")

library(survey)
library(MASS)
library(mvtnorm)
library(foreach)
library(doParallel)
library(purrr)
library(pps)


source("func_MB_GPPP_str1.R")

set.seed(01012021)

V1=4; V2=1; rho=0.5; sigma=0.8; lambda=0.6;
nSim <- 1008
N <- 50000
N_na <- 50
N_a <- N/N_na
Nh_a <- 20
H <- N_a/Nh_a

n1 <- 500
n0 <- 500
n0_na <- 5
n0_a <- n0/n0_na
nBt1 <- 50
nBt2 <- 5

D0 <- rexp(N_a, 1)
pop <- data.frame(D0=D0)
psu <- rep(1:N_a, each=N_na)
pop <- pop[psu, , drop=FALSE]
pop$psu <- psu; psu <- NULL
pop$strat <- rep(1:H, each=N/H)-1

#pop$D1 <- rchisq(N, 10)
pop$D1 <- runif(N, 1, 5)
pop$pi02 <- pop$D1/ave(pop$D1, pop$psu, FUN=sum)

rtrat <- function(x){
uniroot(function(t){max(t+x)/min(t+x)-30}, interval=c(0, 20))$root
}
c <- rtrat(pop$D0)
px0 <- c+pop$D0
max(px0)/min(px0)
#pop$pi0 <- n0*px0/sum(px0)
#pop$pi01 <- px0 <- ave(px0, pop$strat, FUN=function(x)x/sum(x))*n0/H
pop$pi01 <- px0 <- ave(px0, pop$strat, FUN=function(x)x/sum(unique(x)))
sum(pop$pi01)
summary(1/pop$pi01)

pop$pi0 <- pop$pi01*pop$pi02
sum(pop$pi0)
summary(1/pop$pi0)

pop$wght <- 1/pop$pi0
pop$D <- log(pop$wght)

pop$e_x <- rnorm(N, 0, sqrt(V1))

rtcor <- function(x, e, r){
	uniroot(function(t){cor(x+t*e, x) - r}, interval=c(0, 20))$root
}

s1 <- rtcor(pop$D, pop$e_x, r=rho)
pop$X0 <- -7 + pop$D + s1*pop$e_x
cor(pop$X0, pop$D)

pop$X1 <- pop$X0
pop$X2 <- (pop$X0/3)^3
pop$X3 <- exp(pop$X0/2)/5
pop$X4 <- 5*sin(pi*pop$X0/3)
pop$Z0 <- pop$X0^2
pop$XD <- pop$X0*pop$D
pop$D2 <- pop$D^2

summary(pop)

#par(mfrow=c(2, 2))
#plot(pop$X0, pop$X1)
#plot(pop$X0, pop$X2)
#plot(pop$X0, pop$X3)
#plot(pop$X0, pop$X4)

pop$Y1 <- getY(theta=c(1, 1, 1, -0.1, 0.2), X=pop[, c("X1", "D", "strat", "XD")], ICC=0.2, a=N_a, na=N_na, sigma=2, family="gaussian")
pop$Y2 <- getY(theta=c(1, 1, 1, -0.1, 0.2), X=pop[, c("X2", "D", "strat", "XD")], ICC=0.2, a=N_a, na=N_na, sigma=2, family="gaussian")
pop$Y3 <- getY(theta=c(1, 1, 1, -0.1, 0.2), X=pop[, c("X3", "D", "strat", "XD")], ICC=0.2, a=N_a, na=N_na, sigma=2, family="gaussian")
pop$Y4 <- getY(theta=c(1, 1, 1, -0.1, 0.2), X=pop[, c("X4", "D", "strat", "XD")], ICC=0.2, a=N_a, na=N_na, sigma=2, family="gaussian")

summary(pop$Y1)
summary(pop$Y2)
summary(pop$Y3)
summary(pop$Y4)

cor(pop$Y1, 1 + pop$X1 + pop$D - 0.1*pop$strat + 0.2*pop$XD)		# + 0.5*pop$strat
cor(pop$Y2, 1 + pop$X2 + pop$D - 0.1*pop$strat + 0.2*pop$XD) 		# + 0.5*pop$strat
cor(pop$Y3, 1 + pop$X3 + pop$D - 0.1*pop$strat + 0.2*pop$XD) 		# + 0.5*pop$strat
cor(pop$Y4, 1 + pop$X4 + pop$D - 0.1*pop$strat + 0.2*pop$XD) 		# + 0.5*pop$strat
cor(pop[, c("X1", "D", "strat", "XD")])

#rtexp0 <- function(x){
#	uniroot(function(t){mean(exp(t+x)/(1+exp(t+x)))-n0/N}, interval=c(-20, 20))$root
#}
#px0 <- -0.4*pop$D			
#summary(px0)
#gamma0 <- rtexp0(px0)
#pop$pi0 <- exp(gamma0+px0)/(1+exp(gamma0+px0))
#sum(pop$pi0)
#summary(1/pop$pi0)


rtexp1 <- function(x){
	uniroot(function(t){mean(exp(t+x)/(1+exp(t+x)))-n1/N}, interval=c(-20, 20))$root
}
px1 <- lambda*pop$X0
summary(px1)
gamma1 <- rtexp1(px1)
pop$pi1 <- exp(gamma1+px1)/(1+exp(gamma1+px1))
sum(pop$pi1)
summary(1/pop$pi1)
pop$wght0 <- 1/pop$pi1

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

pop$Y <- pop$Y4
pop$X <- pop$X4
mu1 <- mean(pop$Y1)
mu2 <- mean(pop$Y2)
mu3 <- mean(pop$Y3)
mu  <- mean(pop$Y4)

#p0 <- as.numeric(rowsum(pop$pi01, pop$psu))/n0_na
p0 <- pop$pi01[!duplicated(pop$psu)]
h <- rep(1:H, each=N_a/H)

ZZ0 <- ZZ1 <- list()
for(i in 1:nSim){
	tZZ0 <- as.numeric(sapply(1:H, FUN=function(j)((1:N_a)[h==j])[sample(1:Nh_a, n0_a/H, replace=F, prob=p0[h==j]/sum(p0[h==j]))]))
	ZZ0[[i]] <- rep(N_na*(tZZ0-1), each=n0_na) + as.numeric(sapply(tZZ0, function(x)sample(N_na, n0_na, replace=F, prob=pop$pi02[pop$psu==x])))
	ZZ1[[i]] <- (1:N)[as.logical(rbinom(N, 1, pop$pi1))]
}

mean(sapply(ZZ0, length))
mean(sapply(ZZ1, length))

#####################################################

b <- 4
p <- 4 + 4*b

cores <- as.numeric(Sys.getenv('SLURM_NTASKS_PER_NODE', unset='8'))
cl <- makeCluster(cores)
registerDoParallel(cl)

combine <- function(x, ...) {  
  mapply(rbind, x, ..., SIMPLIFY = FALSE)
}

outx <- list(Method=c("GPPP", "PSPP", "LPW", "AIPW"), Mu=rep(NA, 4), LCL=rep(NA, 4), UCL=rep(NA, 4), SE=rep(NA, 4), FW_Mu=NA, FW_SE=NA)

#start_time <- Sys.time()
res <- foreach(i=1:nSim, .combine=combine, .multicombine = TRUE, .packages = c("MASS", "survey", "rstan", "polyapost", "mgcv")) %dopar% {
	b <- 4
	p <- 4 + 4*b
	res1 <- as.data.frame(matrix(NA, 1, 10*p))
	names(res1) <- c(paste0("mu", 1:p), paste0("se", 1:p), paste0("bias", 1:p), paste0("rmse", 1:p), paste0("ci_ll", 1:p), paste0("ci_ul", 1:p), paste0("cov", 1:p), paste0("cle", 1:p))

	smp0 <- pop[ZZ0[[i]], ]
	smp1 <- pop[ZZ1[[i]], ]

	#start_time <- Sys.time()
	 out10 <- tryCatch(maGPPP(psmp=smp0[, c("X", "X0", "XD", "Z0", "D", "D2", "strat", "psu", "wght", "Y")], nsmp=smp1[, c("X", "X0", "XD", "Z0", "D", "D2", "strat", "psu", "wght", "Y")], N=N, method=c("GPPP", "PSPP", "LWP", "AIPW"), spec=c(T, T), family="gaussian", Bt1=nBt1, Bt2=nBt2), error = function(e){outx})
	 out11 <- tryCatch(maGPPP(psmp=smp0[, c("X", "X0", "XD", "Z0", "D", "D2", "strat", "psu", "wght", "Y")], nsmp=smp1[, c("X", "X0", "XD", "Z0", "D", "D2", "strat", "psu", "wght", "Y")], N=N, method=c("GPPP", "PSPP", "LWP", "AIPW"), spec=c(F, T), family="gaussian", Bt1=nBt1, Bt2=nBt2), error = function(e){outx})
	 out12 <- tryCatch(maGPPP(psmp=smp0[, c("X", "X0", "XD", "Z0", "D", "D2", "strat", "psu", "wght", "Y")], nsmp=smp1[, c("X", "X0", "XD", "Z0", "D", "D2", "strat", "psu", "wght", "Y")], N=N, method=c("GPPP", "PSPP", "LWP", "AIPW"), spec=c(T, F), family="gaussian", Bt1=nBt1, Bt2=nBt2), error = function(e){outx})
	 out13 <- tryCatch(maGPPP(psmp=smp0[, c("X", "X0", "XD", "Z0", "D", "D2", "strat", "psu", "wght", "Y")], nsmp=smp1[, c("X", "X0", "XD", "Z0", "D", "D2", "strat", "psu", "wght", "Y")], N=N, method=c("GPPP", "PSPP", "LWP", "AIPW"), spec=c(F, F), family="gaussian", Bt1=nBt1, Bt2=nBt2), error = function(e){outx})
	#end_time <- Sys.time()

	# Estimating unweighted mean and 95%CI	
	res1$mu1[1] <- mean(smp0$Y, na.rm=T)
	res1$se1[1] <- sd(smp0$Y, na.rm=T)/sqrt(n0)
	res1$ci_ll1[1] <- res1$mu1[1] - qnorm(0.975)*res1$se1[1]
	res1$ci_ul1[1] <- res1$mu1[1] + qnorm(0.975)*res1$se1[1]

	# Estimating fully weighted mean and 95%CI
	res1$mu2[1] <- out10$FW_Mu		#weighted.mean(smp0$Y, w=smp0$wght, na.rm=T)
	res1$se2[1] <- out10$FW_SE		#weighted.seh(smp0$Y, smp0$wght, smp0$psu, smp0$strat)
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

	res1
}
#end_time <- Sys.time()
#end_time - start_time

stopCluster(cl)

res1 <- as.data.frame(res)

res1[, paste0("bias", 1:p)] <- res1[, paste0("mu", 1:p)]-mu
res1[, paste0("rmse", 1:p)] <- res1[, paste0("bias", 1:p)]^2
res1[, paste0("cov", 1:p)] <- res1[, paste0("ci_ll", 1:p)]<mu & mu<res1[, paste0("ci_ul", 1:p)]
res1[, paste0("cle", 1:p)] <- res1[, paste0("ci_ul", 1:p)] - res1[, paste0("ci_ll", 1:p)]

write.csv(res1, paste0("res1_TT_glm_BS_", n0, "_", n1, "_", rho, "_", lambda, "_4.csv"), col.names=T, sep=T)

mns <- c("GPPP", "PSPP", "LWP", "AIPW")
tb <- as.data.frame(matrix(NA, p, 11))
names(tb) <- c("method", "bias", "rMSE", "cov_rate", "cov_len", "se_ratio", "MeanSE", "TrueSE", "Q05", "Q95", "biasSE")
tb$method <- c("UW_PS", "TW_PS", "UW_NS", "TW_NS", paste(mns, "TT", sep="_"), paste(mns, "FT", sep="_"), paste(mns, "TF", sep="_"), paste(mns, "FF", sep="_"))

tb[, 2] <- 100*apply(res1[, paste0("bias", 1:p)], 2, function(x)mean(x, na.rm=T))/mu
tb[, 3] <- 100*sqrt(apply(res1[, paste0("rmse", 1:p)], 2, function(x)mean(x, na.rm=T)))/mu
tb[, 4] <- 100*apply(res1[, paste0("cov", 1:p)], 2, function(x)mean(x, na.rm=T))
tb[, 5] <- 100*apply(res1[, paste0("cle", 1:p)], 2, function(x)mean(x, na.rm=T))/mu
tb[, 6] <- apply(res1[, paste0("se", 1:p)], 2, function(x)mean(x, na.rm=T))/apply(res1[, paste0("mu", 1:p)], 2, function(x)sd(x, na.rm=T))
tb[, 7] <- apply(res1[, paste0("se", 1:p)], 2, function(x)mean(x, na.rm=T))
tb[, 8] <- apply(res1[, paste0("mu", 1:p)], 2, function(x)sd(x, na.rm=T))
tb[, 9] <- 100*apply(res1[, paste0("bias", 1:p)], 2, function(x)quantile(x, probs=0.05, na.rm=T))/mu
tb[, 10] <- 100*apply(res1[, paste0("bias", 1:p)], 2, function(x)quantile(x, probs=0.95, na.rm=T))/mu
tb[, 11] <- 100*apply(res1[, paste0("mu", 1:p)], 2, function(x)sd(x, na.rm=T))/mu

tb

write.csv(tb, paste0("tb_TT_glm_BS_", n0, "_", n1, "_", rho, "_", lambda, "_4.csv"), col.names=T, sep=",")
