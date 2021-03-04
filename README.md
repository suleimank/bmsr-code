# Bayesian multi-source regression

## Contents
The repository contains the following:

* BMSR: Bayesian multi-source regression
	* STAN code for Bayesian multi-source regression (bmsr.stan)
	* R script for training and predicting from the BMSR method (bmsr.stan)
	* Demo code to run the method on simulated data (demo_bmsr.stan)
 
* BMSMTR: Bayesian multi-source multi-task regression
	* STAN code for Bayesian multi-source multi-task regression (bmsmtr.stan)
	* R script for training and predicting from the BMSMTR method (bmsr.stan)
	* Demo code to run the method on simulated data (demo_bmsmtr.stan)

## Description
The bmsr package implements joint regression from multiple data sources in a Bayesian framework. The package provides implementation for both single-task and multi-task regression. The model is implemented using STAN and interface is provided using R programming language. Options for training the model using both NUTS sampler and variational inference are provided. The package is structured for ease of use and the included demo shows the model execution on real-life as well as simulated datasets.

## R-package
For an Installable package see [R-package](https://github.com/suleimank/bmsr) 

## Citation
Cite as: To Appear <citation information comes here>
