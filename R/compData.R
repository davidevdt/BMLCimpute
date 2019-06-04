#' @title Create a completed dataset with the imputated values of BMLCimpute (compData)  
#'  
#' 
#'
#' @description Plug the imputations obtained with multilevelLCMI into the original dataset, in order to obtain a 
#' completed dataset.  
#'
#'
#' @details This function takes a 'convData' list, the imputations provided by 'multilevelLCMI'
#' and the imputation index (ind in {1,..., M} where M is the number of imputations) and returns
#' the completed dataset. 
#'
#' 
#' 
#' @param convData Ouptut list produced by the 'convData' function 
#' @param implev1 The set of imputations for the level-1 variables provided by the 'multilevelLCMI' function.
#'        It corresponds to the first element of the list returned by 'multilevelLCMI'.
#'
#' @param implev2 The set of imputations for the level-2 variables (when present) provided by the
#'        'multilevelLCMI' function. It corresponds to the second element of the list returned by 'multilevelLCMI'.
#'
#' @param ind The imputation index; an integer value that ranges from 1 to M, where M is the number of imputations
#'        computed by \code{BMLCimpute}. 
#'
#' @return The imputed dataset    
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
compData <- function(convData, implev1, implev2 = NULL, ind) {
    
    if (class(convData) != "convData") {
        stop("convData object not provided.")
    }
    
    if (class(implev1) != "impmlcmi") {
        stop("<implev1> must be a set of imputations coming from a mlcmi object.")
    }
    
    if ((!is.null(implev2)) & (class(implev2) != "impmlcmi")) {
        stop("<implev1> must be a set of imputations coming from a mlcmi object.")
    }
    
    doVar2 <- convData$doVar2  # Boolean; are level-2 variables present in the dataset? 
    GID_ <- convData$GID  # Group ID column index 
    UID <- convData$UID  # Unit ID column index 
    var2 <- convData$var2  # level-2 variable indices; it can be NULL
    if (!is.null(convData$UID)) {
        var2 <- var2 - 1
    }
    dat <- convData$convDat  # The  dataset converted by convData 
    GN <- convData$groupName  # Group ID variable name
    CN <- convData$caseName  # Unit ID variable name
    CID <- convData$caseID  # Case ID vector 
    SC <- convData$sort_  # Original permutation of the data rows
    GIDC <- convData$GroupIDs[1, ]  # Group ID vector 
    
    
    Y <- as.data.frame(dat)  # Dataset 
    n <- nrow(dat)  # (Total) Sample size
    GID <- 1
    J <- length(table(Y[, GID]))
    
    
    # Preparation objects for level-2 variables imputation
    TT2 <- 1  # Number of level-2 variables 
    ncat2 <- rep(0, 1)  # Number of categories for each level-2 variable  
    cod2 <- 0  # Codes for variables at level-2 
    nms2 <- 0  # Names of level-2 variables 
    Y2 <- matrix(0, 1, 1)  # Level-2 dataset 
    R2 <- matrix(0, 1, 1)  # Level-2 missingess indicator 
    YY2 <- matrix(0, 1, 1)  # Auxiliary matrix 
    
    GU <- convData$convDat[, 1]  # Group ID recoded 
    
    # If variables at level-2 are present, initialize the objects just created
    if (doVar2) {
        Y2 <- cbind(GU, Y[, var2])
        Y2 <- unique(Y2)
        Y2 <- Y2[order(Y2[, 1]), ]
        Y2b <- Y2
        R2 <- ifelse(is.na(Y2[, -1]), 1, 0)
        nms2 <- convData$namesLev2
        TT2 <- ncol(Y2[, -1])
        ncat2 <- convData$nCatLev2
        cod2 <- convData$codLev2
        piresp2 <- vector("list", TT2)
        nlat2 <- vector("list", TT2)
        YY2 <- as.matrix(Y2[, -1])
    }
    
    # Level-2 unit-specific sample size
    if (J > 1) {
        nj <- table(GU)
    } else {
        nj <- rep(0, 1)
        nj[1] <- n
    }
    
    # Isolate level-1 data from ID and level-2 variables
    if (doVar2 == FALSE) {
        YY <- Y[, -c(GID)]
        colnames(YY) <- colnames(Y[, -c(GID)])
    } else {
        YY <- Y[, -c(GID, var2)]
        colnames(YY) <- colnames(Y[, -c(GID, var2)])
    }
    nms <- colnames(YY)
    
    # Level-1 imputation : Initializations
    
    TT <- ncol(YY)  # Total number of level-1 variables 
    ncat <- convData$nCatLev1  # Number of categories for each variable  
    cod <- convData$codLev1  # Codes for the categories (original/new) 
    Y <- cbind(GU, YY)  # Level-1 dataset 
    
    datJ <- vector("list", J)  # Level-1 group-specific data 
    R <- vector("list", J)  # Level-1 group-specific missingness indicator 
    
    imp <- implev1  # Level-1 imputations 
    imp2 <- implev2  # Level-2 imputations 
    rm(implev1)
    rm(implev2)
    rm(YY)
    
    for (j in 1:J) {
        datJ[[j]] <- matrix(0, nj[j], TT)
        R[[j]] <- matrix(0, nj[j], TT)
        datJ[[j]] <- Y[GU == j, 2:(TT + 1)]
        colnames(datJ[[j]]) <- nms
        R[[j]] <- ifelse(is.na(datJ[[j]]), 1, 0)
        datJ[[j]] <- as.matrix(datJ[[j]])
        R[[j]] <- as.matrix(R[[j]])
    }
    
    
    
    # Fill the missing entries of the dataset with the imputations
    
    if (doVar2) {
        
        cmpletedata <- NULL
        
        YR2 <- as.matrix(Y2b[, -1])  # Auxiliary level-2 matrix 
        
        # Fill the missing entries of YR2
        for (t2 in 1:TT2) {
            for (m2 in 1:ncat2[t2]) {
                YR2[which(Y2b[, t2 + 1] == cod[[t2]][2, m2]), t2] <- as.numeric(noquote(cod[[t2]][1, m2]))
            }
            
            
            if (sum(R2[, t2]) > 0) {
                YR2[as.numeric(row.names(imp2[[t2]])), t2] <- imp2[[t2]][, ind]
            }
        }
        
        YR <- vector("list", J)  # Auxiliary level-1 matrices 
        
        # Fill the missing entries of YR
        for (j in 1:J) {
            YR[[j]] <- as.matrix(Y[GU == j, -1])
            
            for (tt in 1:TT) {
                for (m in 1:ncat[tt]) {
                  YR[[j]][which(Y[GU == j, tt + 1] == cod[[tt]][2, m]), tt] <- as.numeric(noquote(cod[[tt]][1, m]))
                }
                if (sum(R[[j]][, tt]) > 0) {
                  YR[[j]][as.numeric(row.names(imp[[j]][[tt]])), tt] <- imp[[j]][[tt]][, ind]
                }
            }
            
            if (is.null(UID)) {
                YR[[j]] <- cbind(rep(j, nj[j]), YR[[j]])
            } else {
                YR[[j]] <- cbind(rep(j, nj[j]), CID[GU == j], YR[[j]])
            }
            cmpletedata <- rbind(cmpletedata, YR[[j]])
            
        }
        
        cmpletedata <- as.data.frame(cmpletedata)
        YR2 <- as.data.frame(YR2)
        
        # Assemble level-1 and level-2 variables
        cmpletedata <- data.frame(cmpletedata, matrix(NA, nrow(cmpletedata), TT2))
        for (t2 in 1:TT2) {
            for (j in 1:J) {
                if (is.null(UID)) {
                  cmpletedata[which(cmpletedata[, 1] == j), t2 + (TT + 1)] <- YR2[j, t2]
                } else {
                  cmpletedata[which(cmpletedata[, 1] == j), t2 + (TT + 2)] <- YR2[j, t2]
                }
            }
        }
        
        cmpletedata <- as.data.frame(cmpletedata)
        
        if (is.null(UID)) {
            cmpletedata[, 1] <- rep(GIDC, nj)
            colnames(cmpletedata) <- c(GN, nms, nms2)
        } else {
            colnames(cmpletedata) <- c(GN, CN, nms, nms2)
            if (UID < GID_) {
                tmp_ <- cmpletedata[, 1]
                cmpletedata[, 1] <- cmpletedata[, 2]
                cmpletedata[, 2] <- tmp_
                colnames(cmpletedata) <- c(CN, GN, nms, nms2)
                rm(tmp_)
            }
        }
        
        rownames(cmpletedata) <- SC
        SC_ = sort(SC)
        cmpletedata <- cmpletedata[order(SC_), ]
        
        if (is.null(UID)) {
            index_1 = which((!((1:(TT + TT2 + 1)) %in% (var2))))[-1]
            index_2 = which((1:(TT + TT2 + 1) %in% (var2)))
            nm_tmp_1 = colnames(cmpletedata)[2:(1 + TT)]
            nm_tmp_2 = colnames(cmpletedata)[(TT + 2):(dim(cmpletedata)[2])]
            tmpVar1 = cmpletedata[, 2:(1 + TT)]
            tmpVar2 = cmpletedata[, (TT + 2):(dim(cmpletedata)[2])]
            cmpletedata[, index_1] <- tmpVar1
            cmpletedata[, index_2] <- tmpVar2
            colnames(cmpletedata)[index_1] <- nm_tmp_1
            colnames(cmpletedata)[index_2] <- nm_tmp_2
            rm(tmpVar1)
            rm(tmpVar2)
            rm(index_1)
            rm(index_2)
            rm(nm_tmp_1)
            rm(nm_tmp_2)
        } else {
            index_1 = which((!((1:(TT + TT2 + 2)) %in% (var2 + 1))))[-c(1, 2)]
            index_2 = which((1:(TT + TT2 + 2) %in% (var2 + 1)))
            nm_tmp_1 = colnames(cmpletedata)[3:(2 + TT)]
            nm_tmp_2 = colnames(cmpletedata)[(TT + 3):(dim(cmpletedata)[2])]
            tmpVar1 = cmpletedata[, 3:(2 + TT)]
            tmpVar2 = cmpletedata[, (TT + 3):(dim(cmpletedata)[2])]
            cmpletedata[, index_1] <- tmpVar1
            cmpletedata[, index_2] <- tmpVar2
            colnames(cmpletedata)[index_1] <- nm_tmp_1
            colnames(cmpletedata)[index_2] <- nm_tmp_2
            rm(tmpVar1)
            rm(tmpVar2)
            rm(index_1)
            rm(index_2)
            rm(nm_tmp_1)
            rm(nm_tmp_2)
        }
    } else {
        # Fill missing entries only on level-1 variables
        cmpletedata <- NULL
        
        YR <- vector("list", J)
        for (j in 1:J) {
            
            YR[[j]] <- as.matrix(Y[GU == j, -1])
            
            for (tt in 1:TT) {
                for (m in 1:ncat[tt]) {
                  YR[[j]][which(Y[GU == j, tt + 1] == cod[[tt]][2, m]), tt] <- as.numeric(noquote(cod[[tt]][1, m]))
                }
                
                if (sum(R[[j]][, tt]) > 0) {
                  YR[[j]][as.numeric(row.names(imp[[j]][[tt]])), tt] <- imp[[j]][[tt]][, ind]
                }
            }
            
            if (J > 1) {
                if (is.null(UID)) {
                  YR[[j]] <- cbind(rep(j, nj[j]), YR[[j]])
                } else {
                  YR[[j]] <- cbind(rep(j, nj[j]), CID[GU == j], YR[[j]])
                }
                cmpletedata <- rbind(cmpletedata, YR[[j]])
            } else {
                if (is.null(UID)) {
                  cmpletedata <- YR[[j]]
                } else {
                  cmpletedata <- cbind(CID, YR[[j]])
                }
            }
        }
        
        cmpletedata <- as.data.frame(cmpletedata)
        
        if (J > 1) {
            cmpletedata[, 1] <- rep(GIDC, nj)
            if (is.null(UID)) {
                colnames(cmpletedata) <- c(GN, nms)
            } else {
                colnames(cmpletedata) <- c(GN, CN, nms)
                if (UID < GID_) {
                  tmp_ <- cmpletedata[, 1]
                  cmpletedata[, 1] <- cmpletedata[, 2]
                  cmpletedata[, 2] <- tmp_
                  colnames(cmpletedata) <- c(CN, GN, nms)
                  rm(tmp_)
                }
            }
        } else {
            if (is.null(UID)) {
                colnames(cmpletedata) <- nms
            } else {
                colnames(cmpletedata) <- c(CN, nms)
            }
            rownames(cmpletedata) <- SC
            SC_ = sort(SC)
            cmpletedata <- cmpletedata[order(SC_), ]
        }
    }
    
    
	dataset <- cmpletedata
    
	cat("Dataset filled with imputing values.\n")
    return(invisible( dataset ))
    
}


















