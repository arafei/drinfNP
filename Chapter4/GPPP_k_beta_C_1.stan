//stan function running join PS and outcome modeling using Gaussian Process of Propensity Prediction (GPPP)
//note 1: This set of stan codes is for continuous gaussian outcomes
//note 2: \pi^R_i is estimated through a linear regression with a logit-transformed outcome
functions {
	int[] rep_int(int x, int n){
		int out_int[n];
		for(i in 1:n) out_int[i] = x;
		return out_int;
	}
	//Generalized inverse gaussian distribution
	real generalized_inverse_gaussian_lpdf(real x, int p, real a, real b) {
    		return  p * 0.5 * log(a / b) 
			- log(2 * modified_bessel_second_kind(p, sqrt(a * b)))
			+ (p - 1) * log(x) - (a * x + b / x) * 0.5;
	}
	//Spectral density function
	real lambda(real L, int m) {
		real lam;
		lam = ((m*pi())/(2*L))^2;		
		return lam;
	}
	real spd(real alpha, real rho, real w, int v) {
		real S;
		if(v==9){
			S = (alpha^2)*sqrt(5*pi()/4)*rho*((sqrt(5)*rho*fabs(w))^2)*modified_bessel_second_kind(2, sqrt(5)*rho*fabs(w))/tgamma(2.5);	//QR cov for v=3/2
		}
		if(v==7){
			S = 3*sqrt(pi())*(alpha^2)*(rho^2)*fabs(w)*modified_bessel_second_kind(1, sqrt(3)*rho*fabs(w))/tgamma(1.5);	//QR cov for v=3/2
		}
		if(v==5){
			S = (4*(alpha^2)*sqrt(pi()*(5^5))/(tgamma(2.5)*rho^5))*1.0/(((5/(rho^2)) + 4*(pi()^2)*w^2)^3);	//Matern cov for v=3/2
		}
		if(v==3){
			  S = 4*alpha^2 * (sqrt(3)/rho)^3 * 1.0/((sqrt(3)/rho)^2 + w^2)^2;
			//S = (2*(alpha^2)*sqrt(pi()*(3^3))/(tgamma(1.5)*rho^3))*1.0/(((3/(rho^2)) + 4*(pi()^2)*w^2)^2);	//Matern cov for v=3/2
		}
		if(v==2){
			S = (alpha^2) * sqrt(2*pi()) * rho * exp(-0.5*(rho^2)*(w^2));					//squared exponential cov
		}
		if(v==1){
			S = (2*(alpha^2)*sqrt(pi())/(tgamma(0.5)*rho))*1.0/((1/rho)^2 + 4*(pi()^2)*w^2);			//absolute exponential cov (Matern v=1/2)
		}
		return S;
	}
	vector phi(real L, int m, vector x) {
		vector[rows(x)] fi;
		fi = 1/sqrt(L) * sin(m*pi()/(2*L) * (x+L));			
		return fi;
	}
	vector agp(vector X, vector beta, real rho, real alpha, real L, int M, int N, int V){
		matrix[N, M] PHI;
		vector[M] diagSPD;
		vector[N] out;
		for(m in 1:M){
			PHI[:,m] = phi(L, m, X);
			diagSPD[m] =  sqrt(spd(alpha, rho, sqrt(lambda(L, m)), V)); 
		}
		out = PHI[,] * (diagSPD .* beta);
		return out;
	}
}
data {
	real C;						//boundary condition factor
	int  V;						//choice of kernel
	int Npop;					//Population true size
	int<lower=1> M;					//number of basis functions		
	int<lower=1> N0;				//sample size of probability sample
	int<lower=1> N1;				//sample size of non-probability sample
	int<lower=1> Kp;				//number of covariates for propensity model
	int<lower=1> Ky;				//number of covariates for outcome model
	int<lower=1> K;					//total number of auxiliary variables
	int<lower=1> Vp;				//indicator of covariates for propensity modeling
	int<lower=1> Vy[Ky];				//indicator of covariates for outcome modeling
	int<lower=1> ind1[N1];				//indices of training observations
	int<lower=1> ind0[N0];				//indices of test observations
	matrix[N0+N1, K] X;				//matrix of auxiliary variables for the combined sample
	vector[N1] Y;					//outcome variable in the non-probability sample
	int<lower=0,upper=1> Z[N0+N1];			//indicator of being in non-probability sample
	vector<lower=1>[N0] W;				//sampling weights of the reference survey
	real sig[8];
}
transformed data {
	int<lower=1> N=N0+N1;
	vector<lower=0,upper=1>[N0] Pi = inv(W);
	real mean_Y;
	real sd_Y;
	vector[N1] Yc;
	matrix[N, K] Xc;  				// centered version of X without an intercept
	vector[K] means_X1;  				// column means of X before centering
	for (i in 1:K) {
		Xc[, i] = (X[, i] - mean(X[ind1, i]));		//sd(X[ind1, i])
	}
	mean_Y = mean(Y);
	sd_Y = sd(Y);
	Yc = (Y - mean_Y);
}
parameters {
	real b0;			// temporary intercept for centered predictors
	vector[Ky] b1;			// population-level effects
	vector[M]  b2;
	real gp0;
	real a0;
	real a1;
	real c0;
	real c1;
	real<lower=0> theta;
	real<lower=0> rho;
	real<lower=0> sigma;
	real<lower=0> alpha;
	simplex[N0] psi;
}
transformed parameters{
	real L;
	vector[N] f;
	vector<lower=0,upper=1>[N0] psi_s;
	vector<lower=0,upper=1>[N] Pi_hat = inv_logit(c0 + Xc[:,Vp]*c1);
	vector[N] Z_hat = a0 + Xc[:,Vp]*a1;
	vector[N] PS_log = Z_hat + log(Pi_hat);
	PS_log = PS_log - mean(PS_log[ind1]);											//sd(PS_log[ind1])			
	L = C*max(fabs(PS_log));
	f = PS_log * gp0 + agp(PS_log, b2, rho, alpha, L, M, N, V) + b0 + Xc[,Vy] * b1;						//
	psi_s = (psi .* (1-Pi) ./ Pi)/sum(psi .* (1-Pi) ./ Pi);
}
model{
	// priors
	target += student_t_lpdf(c0 | sig[1], 0, 1);									//normal_lpdf(c0 | 0, sig[1]);
	target += student_t_lpdf(c1 | sig[1], 0, 1);									//normal_lpdf(c1 | 0, sig[1]);
	target += student_t_lpdf(a0 | sig[2], 0, 1);									//normal_lpdf(a0 | 0, sig[2]);
	target += student_t_lpdf(a1 | sig[2], 0, 1);									//normal_lpdf(a1 | 0, sig[2]);
	target += student_t_lpdf(b0 | sig[3], 0, 1);									//normal_lpdf(b0 | 0, sig[3]);
	target += student_t_lpdf(b1 | sig[3], 0, 1);									//normal_lpdf(b1 | 0, sig[3]);
	target += student_t_lpdf(gp0 | sig[3], 0, 1);									//normal_lpdf(gp0 | 0, sig[3]);
	target += normal_lpdf(b2 | 0, 1);										//student_t_lpdf(b2 | sig[3], 0, 1);
	target += student_t_lpdf(theta | sig[4], 0, 1) - 1 * student_t_lccdf(0 | sig[4], 0, 1);				//inv_gamma_lpdf(theta | sig[4], sig[4]-1);
	target += student_t_lpdf(sigma | sig[5], 0, 1) - 1 * student_t_lccdf(0 | sig[5], 0, 1);				//inv_gamma_lpdf(sigma | sig[5], sig[5]-1);
	target += student_t_lpdf(alpha | sig[6], 0, 1) - 1 * student_t_lccdf(0 | sig[6], 0, 1);				//inv_gamma_lpdf(alpha  | sig[6], sig[6]-1);		
	target += generalized_inverse_gaussian_lpdf(rho | 0, sig[7], sig[8]);						//inv_gamma_lpdf(rho   | sig[7], sig[8]);
	target += dirichlet_lpdf(psi | rep_vector(1.0, N0));
	// likelihood
	target += beta_lpdf(Pi | Pi_hat[ind0]*theta, (1.0 - Pi_hat[ind0])*theta);
	target += bernoulli_logit_lpmf(Z | Z_hat);									//Binary regression				
	target += normal_lpdf(Yc | f[ind1], sigma);									//Partall linear HSGP
	target += multinomial_lpmf(rep_int(1, N0) | psi);
}
generated quantities{
	//predicting the outcome for units of reference survey
	vector[N0] y_predict = to_vector(normal_rng(f[ind0], sigma)) + mean_Y;
	vector<lower=0>[N0] w_predict = to_vector(multinomial_rng(psi_s, Npop - N0)) + 1.0;
	real ybar = sum(w_predict .* y_predict)/sum(w_predict);
}
