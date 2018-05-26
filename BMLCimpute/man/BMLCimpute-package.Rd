\name{BMLCimpute-package}
\alias{BMLCimpute-package}
\alias{BMLCimpute}
\docType{package}
\title{
	BMLCimpute: Multilevel Latent Class Models for the Multiple Imputation of Nested Categorical Data
}
\description{
	BMLCimpute performs multiple imputation of nested categorical data by means of Bayesian Multilevel Latent Class models. 
}
\details{
	Through the function multilevelLCMI the routine performs multiple imputation of single- or multi-level
	categorical data via single- or multi-level latent class (mixture) models. Data should be first processed with the function convData; the resulting list is then passed as input to the functions multilevelLCMI in order to perform the imputations. Complete datasets are obtained via the function compData.
}
\author{
Davide Vidotto

Maintainer: Davide Vidotto <d.vidotto@uvt.nl>
}
\references{
	Vidotto, Vermunt, van Deun (in press). Multilevel Latent Class models for the Multiple imputation of Nested Categorical Data.
}
\keyword{ multilevel mixture models }
\keyword{ multilevel analysis }
\keyword{ missing data }
\keyword{ Bayesian multiple imputation }
\seealso{
  Rcpp package
}
\examples{
library(BMLCimpute)
data(simul_incomplete)

# PreProcess the Data
cd <- convData(simul_incomplete, GID = 1, UID = 2, var2 = 8:12)

# Model Selection
mmLC<- multilevelLCMI( convData = cd, L = 10, K = 15, it1 = 1000, it2 = 3000, it3 = 100, it.print = 250, v = 10, I = 0, pri2 = 1 / 10, pri1 = 1 / 15, priresp = 0.01, priresp2 = 0.01, random = TRUE, estimates = FALSE, count = TRUE, plot.loglik = FALSE, prec = 3, scale = 1.0)

# Select posterior maxima of the number of classes for the imputations (Other alternatives are possible, such as posterior modes or posterior quantiles)
L = max(which(mmLC[[12]] != 0))
K = max(apply(mmLC[[13]], 1, function(x) max(which( x != 0))), na.rm = TRUE)

# Impute
mmLC<- multilevelLCMI( convData = cd, L = L, K = K, it1 = 2000, it2 = 4000, it3 = 100, it.print = 250, v = 10, I = 5, pri2 = 500, pri1 = 50, priresp = 0.01, priresp2 = 0.01, random = TRUE, estimates = FALSE, count = TRUE, plot.loglik = TRUE, prec = 4, scale = 1.0)

# Obtain the dataset completed with the first set of imputations (ind = 1)
complete_data = compData( convData = cd, implev1 = mmLC[[1]], implev2 = mmLC[[2]], ind = 1 )
}



















