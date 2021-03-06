% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multilevelLCMI.R
\name{multilevelLCMI}
\alias{multilevelLCMI}
\alias{multilevelLCMI.default}
\alias{print.multilevelLCMI}
\title{Multilevel Latent Class models for the Multiple Imputation of Categorical Data (multilevelLCMI).}
\usage{
multilevelLCMI(convData, L, K, it1, it2, it3, it.print, v, I = 5, pri2 = 1, 
   pri1 = 1, priresp = 1, priresp2 = 1, random = TRUE, estimates = TRUE,
   count = FALSE, plot.loglik = FALSE, prec = 3, scale = 1)

\method{multilevelLCMI}{default}(convData, L, K, it1, it2, it3, it.print, v, I = 5, pri2 = 1,
     pri1 = 1, priresp = 1, priresp2 = 1, random = TRUE, estimates = TRUE,
     count = FALSE, plot.loglik = FALSE, prec = 3, scale = 1)

\method{print}{multilevelLCMI}(x, ... )
}
\arguments{
\item{convData}{Dataset produced as output by the \code{convData} function.}

\item{L}{Number of higher-level mixture components. When L=1, single-level Latent Class multiple imputation is performed.}

\item{K}{Number of Latent Classes at the lower-level.}

\item{it1}{Number of Gibbs sampler iterations for the burn-in (must be larger than 0).}

\item{it2}{Number of Gibbs sampler iterations for the imputations.}

\item{it3}{Every \code{it3} iterations, the sampler stores new parameter estimates for the calculation of psoterior estimates. Meaningful only when
\code{estimates=TRUE}.}

\item{it.print}{Every \code{it.print} iterations, the state of the Gibbs sampler is screen-printed.}

\item{v}{The Gibbs sampler will produce the first set of imputations at the iteration number \code{(it1+V)}, where V~Unif(1,v). Subsequent imputations are
automatically spaced from each others across the remaining iterations, so that the last imputation (imputation I) occurs at the iteration number 
\code{(it1+it2)}.}

\item{I}{Number of imputations to be performed.}

\item{pri2}{Hyperparameter value for the higher-level mixture probabilities. Default to 1.}

\item{pri1}{Hyperparameter value for the lower-level mixture probabilities. Default to 1.}

\item{priresp}{Hyperparameter value for the lower-level conditional response probabilities. Default to 1.}

\item{priresp2}{Hyperparameter value for the higher-level conditional response probabilities. Default to 1.}

\item{random}{Logical. Should the model parameters be initialized at random values? If \code{TRUE}, parameters are initialized through draws from 
uniform Dirichlet distributions. If \code{FALSE}, parameters are initialized to be equal to 1/D, with D the number of categories in the (observed/latent) 
variable of interest.}

\item{estimates}{Logical. If \code{TRUE}, the function returns the posterior means of the model parameters. Default to \code{TRUE}.}

\item{count}{Logical. Should the output include the posterior distribution of L and K? (only K if L=1)}

\item{plot.loglik}{Logical. Should the output include a traceplot of the log-likelihood ratios obtained through the Gibbs sampler iterations?
Helpful for assessing convergence.}

\item{prec}{When \code{estimates=TRUE}, \code{prec} defines the number of digits with which the estimates are returned.}

\item{scale}{Re-scale the log-likelihood value by a factor equal to \code{scale}; this parameter is useful to avoid underflow in the calculation
 of the log-likelihood (which can occur in large datasets) and consequently to prevent error messages for the visualization of the log-likelihood traceplot. 
This parameter is meaningful only when \code{plot.loglik} is set to TRUE. Default value equal to 1.0.}

\item{x}{A \code{multilevelLCMI} object (\code{print} method).}

\item{...}{Not used.}
}
\value{
A \code{multilevelLCMI} object, a list containing: 
\item{imp}{Set of imputations for the level-1 variables and units.}
\item{imp2}{Set of imputations for the level-2 variables and units.}
\item{piL}{posterior means of the level-2 class probabilities. Calculated only if \code{estimates = TRUE}.}
\item{piLses}{Posterior standard deviations of the level-2 class probabilities. Calculated only if \code{estimates = TRUE}.}
\item{piK}{Posterior means of the level-1 class probabilities. Calculated only if \code{estimates = TRUE}.}
\item{piKses}{Posterior standard deviations of the level-2 class probabilities. Calculated only if \code{estimates = TRUE}.}
\item{picondlev1}{Posterior means of the level-1 conditional probabilities. Calculated only if \code{estimates = TRUE}.}
\item{picondlev1ses}{Posterior standard deviations of the level-1 conditional probabilities. Calculated only if \code{estimates = TRUE}.}
\item{picondlev2}{Posterior means of the level-2 conditional probabilities. Calculated only if \code{estimates = TRUE}.}
\item{picondlev1ses}{Posterior standard deviations of the level-2 conditional probabilities. Calculated only if \code{estimates = TRUE}.}
\item{DIC}{DIC index for the BMLC model. Calculated only if \code{estimates = TRUE.}}
\item{freqL}{Posterior distribution of the number of latent classes at level-2. Calculated only if \code{count} is set to \code{TRUE}.}
\item{freqK}{Posterior distribution of the number of latent classes at level-1. Calculated only if \code{count} is set to \code{TRUE}.}
\item{time}{Running time of the Gibbs sampler iterations.}
}
\description{
Perform single- and multi-level multiple imputation of categorical data through single/multi level Bayesian Latent Class models.
}
\details{
Function for performing Multiple Imputation with the Bayesian Multilevel Latent Class Model. The model takes the list produced by the 'convData' 
    function as input,  in which the dataset converted and prepared for the imputations is present, along with other parameters specified by the user (e.g., 
    number of latent classes and specification of the prior distribution hyperparameters). The function can also offer (when the corresponding boolean parameter
    is activated) a graphical representation of the posterior distribution of the number of occupied classes during the Gibbs sampler iterations. In this way, 
    multilevelLCMI' can also perform model selection in a pre-imputation stage. For model selection, set \code{count=TRUE}. Symmetric Dirichlet priors are used.
}
\examples{
\dontrun{

library(BMLCimpute)

# Load data 
data(simul_incomplete)

# Preprocess the Data
cd <- convData(simul_incomplete, GID = 1, UID = 2, var2 = 8:12)

# Model Selection
set.seed(1)
mmLC <- multilevelLCMI( convData = cd, L = 10, K = 10, it1 = 1000, it2 = 3000, it3 = 100,
    it.print = 250, v = 10, I = 0, pri2 = 1 / 10, pri1 = 1 / 15, priresp =  0.01,
    priresp2 = 0.01, random = TRUE, estimates = FALSE, count = TRUE, plot.loglik = FALSE,
    prec = 3, scale = 1.0)

# Select posterior maxima of the number of classes for the imputations 
# (Other alternatives are possible, such as posterior modes or posterior quantiles)
L = max(which(mmLC[[12]] != 0))
K = max(apply(mmLC[[13]], 1, function(x) max(which( x != 0))), na.rm = TRUE)

# Perform 5 imutations on the dataset 
mmLC <- multilevelLCMI( convData = cd, L = L, K = K, it1 = 2000, it2 = 4000, it3 = 100, 
     it.print = 250, v = 10, I = 5, pri2 = 500, pri1 = 50, priresp = 0.01, priresp2 = 0.01,
     random = TRUE, estimates = FALSE, count = TRUE, plot.loglik = TRUE, prec = 4, scale = 1.0)

# Obtain the dataset completed with the first set of imputations (ind = 1)
complete_data = compData( convData = cd, implev1 = mmLC[[1]], implev2 = mmLC[[2]], ind = 1 )

} 

}
\author{
D. Vidotto <d.vidotto@uvt.nl>
}
\references{
\itemize{

\item[1] Vidotto D., Vermunt J.K., Van Deun K. (2018). 'Bayesian Multilevel Latent Class Models for the Multiple Imputation of Nested Categorical Data'. Journal 
of Educational and Beahvioral Statistics 43(5), 511-539. 

}
}

