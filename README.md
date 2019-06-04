# BMLCimpute: Bayesian Multilevel Latent Class Models for the Multiple Imputation of Nested Categorical Data   
An R package for the multiple imputation of single-level and nested categorical data by means of Bayesian Multilevel Latent Class models. 

## Author
Davide Vidotto <d.vidotto@uvt.nl> 

## Description
```BMLCimpute``` allows researchers and users of categorical datasets with missing data to perform Multiple Imputation via Bayesian latent class models. 
    Data can be either single- or multi-level. Model estimation and imputations are implemented via a Gibbs sampler run with the Rcpp package interface. 
    The function ```multilevelLCMI``` performs the imputations. Prior to the imputation step, data should be processed with the function ```convData```; the 
    resulting list is then passed as input to the ```multilevelLCMI```. Complete datasets are obtained via the ```compData``` function. Check package
	documentation in ```inst/doc``` for further information. 

## Functions

* ```multilevelLCMI``` for the imputations and model estimation (internally calls Rcpp code) 
* ```convData``` for data preparation (preprocessing) 
* ```compData``` for dataset completion

## Install
devtools::install_github("davidevdt/BMLCimpute")

## Version
0.0.1

## Depends 
R (>= 3.3.3)

## License 
GPL-2



