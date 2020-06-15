# drinfNP
Doubly robust inference for non-probability samples

This package contains a set of functions that adjust for the selection bias in non-probability samples when a relavant probability survey with a set of overlapping covariates is available. To weaken the modeling assumpltions, this package employs two classes of doubly robust adjustment methods, i.e. Augmented Inverse Propensity Weighting and Penalized Spline of Propensity Prediction. 

Furthermore, there are three options for users to generate the pseudo-selection probabilities for units of the non-probability sample:
(1) Propensity Adjusted Probability Prediction (PAPP) (Elliott & Valliant, 2017)
(2) Pseudo-Maximum Likelihood-based estimation (PMLE) (Chen et al 2019)
(3) Pseudo-Maximum Likelihood-based estimation (PMLE) (Dever & Valliant 2010)

In addition to Generalized Linear Models, this package allows users to employ Bayesian Additive Regression Trees (BART), which is a strong flexible non-parametric tool, for predicting both the outcome variable for units of the reference survey and the pseudo-selection probabilities for units of the non-probability sample.

Variance of the doubly robust estimator of the population mean can be obtained through both Bootstrap and Jackknife Repeated Replication techniques.
