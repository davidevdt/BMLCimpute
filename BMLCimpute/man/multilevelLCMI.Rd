\name{multilevelLCMI}
\alias{multilevelLCMI}
\title{
	multilevelLCMI
}
\description{
Perform single- and multi-level multiple imputation of categorical data through single/multi level Bayesian Latent Class analysis.
}
\usage{
multilevelLCMI(convData,L,K,it1,it2,it3,it.print,v,I=5,pri2=1.0,pri1=1.0,priresp=1.0,priresp2=1.0,random=TRUE,estimates=TRUE,count=FALSE,plot.loglik=FALSE,prec=4,scale=1.0)
}
\arguments{
  \item{convData}{
	Dataset produced as output by the 'convData' function. 
}
 \item{L}{
	Number of higher-level mixture components. When L=1, single-level Latent Class multiple imputation is performed.
}
  \item{K}{
	Number of Latent Classes at the lower-level. 
}
  \item{it1}{
	Number of Gibbs sampler iterations for the burn-in (must be larger than 0).
}
  \item{it2}{
	Number of Gibbs sampler iterations for the imputations. 
}
  \item{it3}{
	Every it3 iterations, the sampler stores new parameter estimates for model estimation. Meaningful only when estimates=TRUE. 
}
  \item{it.print}{
	Every it.print iterations, the state of the Gibbs sampler is screen-printed.
}
  \item{v}{
	The Gibbs sampler will produce the first set of imputations at the iteration number (it1+V), where V~Unif(1,v). Subsequent imputations are automatically spaced from each others across
	the remaining iterations, so that the last imputation (imputation I) occurs at the iteration number (it1+it2).
}
  \item{I}{
	Number of imputations to be performed. 
}
  \item{pri2}{
	Hyperparameter value for the higher-level mixture probabilities. Default to 1.  
}
  \item{pri1}{
	Hyperparameter value for the lower-level mixture probabilities. Default to 1.  
}
  \item{priresp}{
	Hyperparameter value for the lower-level conditional response probabilities. Default to 1. 
}
  \item{priresp2}{
	Hyperparameter value for the higher-level conditional response probabilities. Default to 1. 
}
  \item{random}{
	Logical. Should the model parameters be initialized at random values? If TRUE, parameters are initialized through draws from uniform Dirichlet distributions. If FALSE, parameters are initialized to be equal to 1/D, with D the number of categories in the (observed/latent) variable of interest. 
}
  \item{estimates}{
	Logical. If TRUE, the function returns the posterior means of the model parameters. Default to TRUE.
}
  \item{count}{
	Logical. Should the output include the posterior distribution of L and K? (only K if L=1)
}
  \item{plot.loglik}{
	Logical. Should the output include a traceplot of the log-likelihood ratios obtained through the Gibbs sampler iterations? Helpful for assessing convergence. 
}
  \item{prec}{
	When estimates=TRUE, prec defines the number of digits with which the estimates are returned. 
}
 \item{scale}{
	Re-scale the log-likelihood value by a factor equal to 'scale'; this parameter is useful to avoid underflow in the calculation of the log-likelihood (which can occur in large datasets) and consequently to prevent error messages for the visualization of the log-likelihood traceplot. This parameter is meaningful only when 'plot.loglik' is set to TRUE. Default value equal to 1.0.
}
}
\details{
Function for performing Multiple Imputation with the Bayesian Multilevel Latent Class Model. The model takes as input the list produced by the 'convData' function, in which the dataset converted and prepared for the imputations is present, along with other parameters specified by the user (e.g., number of latent classes and specification of the prior distribution hyperparameters). The function can also offer (when the corresponding boolean parameter is activated) a graphical representation of the posterior distribution of the number of occupied classes during the Gibbs sampler iterations. In this way, 'multilevelLCMI' can also perform model selection in a pre-imputation stage. Symmetric Dirichlet priors are used. 
}
\value{
\item{imp}{Set of imputations for the level-1 variables and units.}
\item{imp2}{Set of imputations for the level-2 variables and units.}
\item{piL}{osterior means of the level-2 class probabilities. Calculated only if estimates = TRUE.}
\item{piLses}{Posterior standard deviations of the level-2 class probabilities. Calculated only if estimates = TRUE.}
\item{piK}{Posterior means of the level-1 class probabilities. Calculated only if estimates = TRUE.}
\item{piKses}{Posterior standard deviations of the level-2 class probabilities. Calculated only if estimates = TRUE.}
\item{picondlev1}{Posterior means of the level-1 conditional probabilities. Calculated only if estimates = TRUE.}
\item{picondlev1ses}{Posterior standard deviations of the level-1 conditional probabilities. Calculated only if estimates = TRUE.}
\item{picondlev2}{Posterior means of the level-2 conditional probabilities. Calculated only if estimates = TRUE.}
\item{picondlev2ses}{Posterior standard deviations of the level-2 conditional probabilities. Calculated only if estimates = TRUE.}
\item{DIC}{DIC index for the BMLC model. Calculated only if estimates = TRUE.}
\item{freqL}{Posterior distribution of the number of latent classes at level-2. Calculated only if count is set to TRUE.}
\item{freqK}{Posterior distribution of the number of latent classes at level-1. Calculated only if count is set to TRUE.}
\item{time}{Running time of the Gibbs sampler iterations.}
}
\references{
Vidotto, Vermunt, van Deun. Multilevel Latent Class models for the Multiple imputation of Nested Categorical Data.
}
\author{
	Davide Vidotto (Tilburg University)
}
\seealso{
	Rcpp.
}

\keyword{ multilevel analysis }
\keyword{ multiple imputation }
