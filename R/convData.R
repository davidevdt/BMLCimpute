#' @title Prepare data for imputations with BMLCimpute (convData)
#'
#' @description This function takes a categorical dataset as input (categories can be denoted by numbers)
#'     and returns a list of objects that will be used by the 'multilevelLCMI' function to perform the imputations. 
#'
#' @details 
#' Convert a raw categorical dataset with missing data into one ready to be imputed with the multilevelLCMI function.
#' In particular, the function will transform factor variables into numeric ones, where numbers denote a different category. A 
#' coding list is returned along with the converted dataset.
#' 
#' @param dat Raw (categorical) data frame with missing data. It can also be a data matrix. The \code{GID} and \code{UID} 
#'     arguments, if passed to the function, must be in the first two columns of the dataset.
#'
#' @param GID Group (level-2 unit) indicator (expressed as column number corresponding to the group ID in the dataset).
#'     It can be omitted in single-level datasets.    
#' 
#' @param UID Lower-level unit indicator (expressed as column number corresponding to the unit ID in the dataset). Optional.
#' 
#' @param var2 Higher-level (group-specific) variables (expressed as a vector of column numbers in the dataset corresponding
#'     to the variables measured at the higher levels). Optional.  
#'
#' 
#' @return A \code{convData} object, a list containing the following items: 
#' \item{convDat}{The converted dataset}
#' \item{codLev1}{List containing the new (and original) scores which will be used for the imputations
#' (Level-1 variables).}         
#' \item{codLev1}{Vector containing the number of categories observed for each variable (Level-1 variables).}
#' \item{nCatLev1}{Vector containing the number of categories observed for each variable (Level-1 variables).}
#' \item{codLev2}{List containing the new (and original) scores which will be used for the imputations (Level-2 variables).}
#' \item{nCatLev2}{List containing the new (and original) scores which will be used for the imputations (Level-2 variables).}
#' \item{GroupIDs}{Matrix containing original and new Group ID's.}
#' \item{GID}{The column Group ID number (as entered by the user).}
#' \item{UID}{The column Unit ID number (as entered by the user).}
#' \item{var2}{The column numbers for level-2 variables (as entered in the input).}
#' \item{doVar2}{Boolean. Shall the BMLC model impute variables at level-2? (Result of \code{!is.null(var2))}.}
#' \item{namesLev1}{Vector of variable names (level-1 variables).}
#' \item{namesLev2}{Vector of variable names (level-2 variables).}
#' \item{GroupName}{Group ID variable name.}
#' \item{CaseName}{Unit ID variable name.}
#' \item{caseID}{Unit ID vector (re-permuted).}
#' \item{sort_}{Vector containing the original permutation of the dataset rows.}
#'
#' @author D. Vidotto <d.vidotto@@uvt.nl>
#'
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


#' @export 
convData <- function(dat, GID = NULL, UID = NULL, var2 = NULL) {
    
    # If the Group ID is not in the first 2 columns, stop the program.
    if (!is.null(GID)) {
        if (GID > 2) {
            stop("Group ID must be in either the first two columns.", call. = FALSE)
        }
    }
    
    # If the Unit ID is not in the first 2 columns, stop the program.
    if (!is.null(UID)) {
        if (UID > 2) {
            stop("Case ID must be in either the first two columns.", call. = FALSE)
        }
    }
    
    # Are level-2 variables present in the dataset?
    doVar2 <- TRUE
    if (is.null(var2)) {
        doVar2 <- FALSE
    }
    
    # For single-level datasets : create a fictitious Group-ID (with only 1 group)
    if (!is.null(GID)) {
        GN <- colnames(dat)[GID]
    } else {
        GID <- 1
        GroupID <- rep(1, nrow(dat))
        GN <- "GroupID"
        dat <- data.frame(GroupID, dat)
        colnames(dat)[1] <- GN
    }
    
    
    
    # Create a copy of the dataset, converting categories into numbers; create vector of variable names
    Y <- as.data.frame(dat)
    if (is.null(UID)) {
        Y <- apply(Y, 2, function(x) as.numeric(as.factor(x)))
    } else {
        Y <- apply(Y[, -UID], 2, function(x) as.numeric(as.factor(x)))
    }
    n <- nrow(Y)
    nms_tmp <- colnames(dat)
    
    # Column Name for Unit ID
    CN <- NULL
    CID <- NULL
    if (!is.null(UID)) {
        CN <- colnames(dat[UID])
        CID <- dat[, UID]
    }
    
    # Separate level-2 and level-1 variable names
    nms2 <- NULL
    if (doVar2) {
        if (!is.null(UID)) {
            nms2 <- nms_tmp[var2]
            nms <- nms_tmp[-c(GID, UID, var2)]
        } else {
            nms2 <- nms_tmp[var2]
            nms <- nms_tmp[-c(GID, var2)]
        }
    } else {
        if (!is.null(UID)) {
            nms <- nms_tmp[-c(GID, UID)]
        } else {
            nms <- nms_tmp[-GID]
        }
    }
    rm(nms_tmp)
    
    # New codes for Group ID
    dtGID <- data.frame(dat[, GID], Y[, 1])
    dtGID <- unique(dtGID)
    dtGID <- t(dtGID[order(dtGID[, 2]), ])
    rownames(dtGID) <- cbind("Original GroupID", "New GroupID")
    colnames(dtGID) <- NULL
    
    # Sort rows by Group ID
    sortCase <- 1:n
    sortCase <- sortCase[order(Y[, 1])]
    Y <- Y[order(Y[, 1]), ]
    CID <- CID[order(Y[, 1])]
    
    # Create level-2 variables dataset (and list with new and original codes)
    ncat2 <- rep(0, 1)
    cod2 <- 0
    ncat2 <- 0
    
    if (doVar2) {
        ncat2 <- rep(0, length(var2))
        cod2 <- vector("list", length(var2))
        for (tt in 1:length(var2)) {
            ncat2[tt] <- length(table(Y[, var2[tt]]))
            cod2[[tt]] <- matrix(0, 2, ncat2[tt])
            names(cod2)[tt] <- nms2[tt]
            cod2[[tt]][1, ] <- names(table(dat[, nms2[tt]]))
            cod2[[tt]][2, ] <- names(table(Y[, var2[tt]]))
            rownames(cod2[[tt]]) <- c("Original Values", "Coded Values")
        }
    }
    
    # Original and new codes for level-1 variables
    var1 <- which(colnames(Y) %in% nms)
    ncat <- rep(0, length(var1))
    cod <- vector("list", length(var1))
    
    for (tt in 1:length(var1)) {
        ncat[tt] <- length(table(Y[, var1[tt]]))
        cod[[tt]] <- matrix(0, 2, ncat[tt])
        names(cod)[tt] <- nms[tt]
        cod[[tt]][1, ] <- names(table(dat[, nms[tt]]))
        cod[[tt]][2, ] <- names(table(Y[, var1[tt]]))
        rownames(cod[[tt]]) <- c("Original Values", "Coded Values")
    }
    
    
    
    ret <- list(convDat = Y, codLev1 = cod, nCatLev1 = ncat, codLev2 = cod2, nCatLev2 = ncat2, GroupIDs = dtGID, GID = GID, UID = UID, var2 = var2, 
        doVar2 = doVar2, namesLev1 = nms, namesLev2 = nms2, groupName = GN, caseName = CN, caseID = CID, sort_ = sortCase)
    
    class(ret) <- "convData"
    
    cat("Dataset converted.\n")
    invisible(ret)
    
}
