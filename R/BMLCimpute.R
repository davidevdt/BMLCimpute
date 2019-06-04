#' BMLCimpute : Bayesian Multilevel Latent Class Models for the Multiple Imputation of Nested Categorical Data 
#'
#' @description  A package for the multiple imputation of single-level and nested categorical data by means of Bayesian Multilevel Latent Class models.
#'
#' @details `BMLCimpute` allows researchers and users of categorical datasets with missing data to perform Multiple Imputation via Bayesian latent class models. 
#'     Data can be either single- or multi-level. Model estimation and imputations are implemented via a Gibbs sampler run with the Rcpp package interface. 
#'     The function \code{multilevelLCMI} performs the imputations. Prior to the imputation step, data should be processed with the function \code{convData}; the 
#'     resulting list is then passed as input to the \code{multilevelLCMI}. Complete datasets are obtained via the \code{compData} function.
#'
#' @section Functions: 
#' \itemize{
#'     \item \code{multilevelLCMI} for the imputations and model estimation (internally calls Rcpp code); 
#'     \item \code{convData} for data preparation (preprocessing); 
#'     \item \code{compData} for dataset completion.
#' }
#'
#' @author D. Vidotto <d.vidotto@@uvt.nl> 
#'
#' @docType package 
#' @name BMLCimpute 
#' @aliases BMLCimpute 
#' 
#' @useDynLib BMLCimpute, .registration = TRUE
#'
#' @importFrom Rcpp evalCpp
#' @importFrom grDevices colors 
#' @importFrom graphics axis legend matplot par plot 
#' @importFrom stats rgamma runif 
#' 
#' @references
#' \itemize{
#' 
#' \item [1] Vidotto D., Vermunt J.K., Van Deun K. (2018). 'Bayesian Multilevel Latent Class Models for the Multiple Imputation of Nested Categorical Data'. Journal 
#' of Educational and Beahvioral Statistics 43(5), 511-539. 
#'
#' }
#'
#' @examples 
#' \dontrun{
#' 
#' library(BMLCimpute)
#'
#' # Load data 
#' data(simul_incomplete)
#' 
#' # Preprocess the Data
#' cd <- convData(simul_incomplete, GID = 1, UID = 2, var2 = 8:12)
#' 
#' # Model Selection
#' set.seed(1)
#' mmLC <- multilevelLCMI( convData = cd, L = 10, K = 10, it1 = 1000, it2 = 3000, it3 = 100,
#'     it.print = 250, v = 10, I = 0, pri2 = 1 / 10, pri1 = 1 / 15, priresp =  0.01,
#'     priresp2 = 0.01, random = TRUE, estimates = FALSE, count = TRUE, plot.loglik = FALSE,
#'     prec = 3, scale = 1.0)
#' 
#' # Select posterior maxima of the number of classes for the imputations 
#' # (Other alternatives are possible, such as posterior modes or posterior quantiles)
#' L = max(which(mmLC[[12]] != 0))
#' K = max(apply(mmLC[[13]], 1, function(x) max(which( x != 0))), na.rm = TRUE)
#' 
#' # Perform 5 imutations on the dataset 
#' mmLC <- multilevelLCMI( convData = cd, L = L, K = K, it1 = 2000, it2 = 4000, it3 = 100, 
#'      it.print = 250, v = 10, I = 5, pri2 = 500, pri1 = 50, priresp = 0.01, priresp2 = 0.01,
#'      random = TRUE, estimates = FALSE, count = TRUE, plot.loglik = TRUE, prec = 4, scale = 1.0)
#' 
#' # Obtain the dataset completed with the first set of imputations (ind = 1)
#' complete_data = compData( convData = cd, implev1 = mmLC[[1]], implev2 = mmLC[[2]], ind = 1 )
#' 
#' } 
NULL
