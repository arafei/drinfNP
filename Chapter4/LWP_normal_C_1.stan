//stan function runing join PS and outcome modeling for linear-in-weight prediction (LWP)
//note 1: This set of stan codes is for continuous gaussian outcomes
//note 2: \pi^R_i is estimated through a outcome linear regression with a logit-transformed
functions{
	int[] rep_int(int x, int n){
		int out_int[n];
		for(i in 1:n) out_int[i] = x;
		return out_int;
	}
}
data {
	int Npop;					//Population true size
	int<lower=1> N0;				//sample size in S_R, n_R
	int<lower=1> N1;				//sample size in S_A, n_A
	int<lower=1> Kp;				//number of covariates associated with PS model
	int<lower=1> Ky;				//number of covariates associated with outcome model
	int<lower=1> K;					//total number of covariates
	int<lower=1> Vp;				//variable indicator for covariates of PS model
	int<lower=1> Vy[Ky];				//variable indicator for covariates of outcome model
	int<lower=1> ind0[N0];				//indices of being in S_R
	int<lower=1> ind1[N1];				//indices of being in S_A
	matrix[N0+N1, K] X;				//matrix of total auxiliary variables
	vector[N1] Y;					//vector of outcome variable in S_A
	int<lower=0,upper=1> Z[N0+N1];			//indicator of being in S_A
	vector[N0] W;					//sampling weights in S_R
	vector[5] sig;
}
transformed data {
	int<lower=1> N=N0+N1;				//total sample size
	vector<lower=0,upper=1>[N0] Pi = inv(W);	//selection probabilities in S_R
	vector[N0] Pi_logit = logit(Pi);		//logit of selection probabilities in S_R
	real mean_Y;
	vector[N1] Yc;
	matrix[N, K] Xc;  				// centered version of X without an intercept
	vector[K] means_X1;  				// column means of X before centering
	for (i in 1:K) {
		Xc[, i] = X[, i] - mean(X[ind1, i]);
	}
	mean_Y = mean(Y);
	Yc = Y - mean_Y; 
}
parameters {
	real b0;
	vector[Ky] b1;  				//slopes for predictors of the PM
	real b2;  					//slope associated with pseudo-weights in the PM
	real<lower=0> sigma;				//error variance of the PM
	real a0;
	real a1;
	real c0;
	real c1;
	real<lower=0> theta;
	simplex[N0] psi;
}
transformed parameters{
	vector[N] f_train;
	vector[N0] psi_s;
	vector[N] P_hat = c0 + Xc[:,Vp]*c1;
	vector[N] Z_hat = a0 + Xc[:,Vp]*a1;
	vector<lower=0,upper=1>[N] PS = inv_logit(P_hat);
	vector<lower=0>[N] PW = exp(-Z_hat) .* inv(PS);
	f_train= b0 + Xc[:,Vy] * b1 + PW * b2;
	psi_s = psi .* (1-Pi) ./ Pi;
	psi_s = psi_s/sum(psi_s);
}
model{
	// priors
	target += student_t_lpdf(c0 | sig[1], 0, 1);									//normal_lpdf(c0 | 0, sig[1]);
	target += student_t_lpdf(c1 | sig[1], 0, 1);									//normal_lpdf(c1 | 0, sig[1]);
	target += student_t_lpdf(a0 | sig[2], 0, 1);									//normal_lpdf(a0 | 0, sig[2]);
	target += student_t_lpdf(a1 | sig[2], 0, 1);									//normal_lpdf(a1 | 0, sig[2]);
	target += student_t_lpdf(b0 | sig[3], 0, 1);									//normal_lpdf(b0 | 0, sig[3]);
	target += student_t_lpdf(b1 | sig[3], 0, 1);									//normal_lpdf(b1 | 0, sig[3]);
	target += student_t_lpdf(b2 | sig[3], 0, 1);									//normal_lpdf(gp0 | 0, sig[3]);
	target += student_t_lpdf(theta | sig[4], 0, 1) - 1 * student_t_lccdf(0 | sig[4], 0, 1);		//inv_gamma_lpdf(theta | sig[4], sig[4]-1);						//
	target += student_t_lpdf(sigma | sig[5], 0, 1) - 1 * student_t_lccdf(0 | sig[5], 0, 1);		//inv_gamma_lpdf(sigma | sig[5], sig[5]-1);						//	
	target += dirichlet_lpdf(psi | rep_vector(1.0, N0));
	// likelihood
	target += normal_lpdf(Pi_logit | P_hat[ind0], theta);
	target += bernoulli_logit_lpmf(Z | Z_hat);
	target += normal_lpdf(Yc | f_train[ind1], sigma);
	target += multinomial_lpmf(rep_int(1, N0) | psi);
}
generated quantities{
	//predicting the outcome for units of reference survey
	vector[N0] y_predict = to_vector(normal_rng(f_train[ind0], sigma)) + mean_Y;
	vector[N0] w_predict = to_vector(multinomial_rng(psi_s, Npop - N0)) + 1.0;
	real ybar = (sum(w_predict .* y_predict) + sum(Yc - f_train[ind1]))/Npop;
}
