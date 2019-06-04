#' @title Multilevel Latent Class models for the Multiple Imputation of Categorical Data (multilevelLCMI).
#'
#' @aliases multilevelLCMI multilevelLCMI.default print.multilevelLCMI 
#'
#' @usage 
#'
#' multilevelLCMI(convData, L, K, it1, it2, it3, it.print, v, I = 5, pri2 = 1, 
#'    pri1 = 1, priresp = 1, priresp2 = 1, random = TRUE, estimates = TRUE,
#'    count = FALSE, plot.loglik = FALSE, prec = 3, scale = 1)
#'
#' \method{multilevelLCMI}{default}(convData, L, K, it1, it2, it3, it.print, v, I = 5, pri2 = 1,
#'      pri1 = 1, priresp = 1, priresp2 = 1, random = TRUE, estimates = TRUE,
#'      count = FALSE, plot.loglik = FALSE, prec = 3, scale = 1)
#'
#' \method{print}{multilevelLCMI}(x, ... )
#' 
#'
#' @description Perform single- and multi-level multiple imputation of categorical data through single/multi level Bayesian Latent Class models. 
#' 
#' 
#' @details Function for performing Multiple Imputation with the Bayesian Multilevel Latent Class Model. The model takes the list produced by the 'convData' 
#'     function as input,  in which the dataset converted and prepared for the imputations is present, along with other parameters specified by the user (e.g., 
#'     number of latent classes and specification of the prior distribution hyperparameters). The function can also offer (when the corresponding boolean parameter
#'     is activated) a graphical representation of the posterior distribution of the number of occupied classes during the Gibbs sampler iterations. In this way, 
#'     multilevelLCMI' can also perform model selection in a pre-imputation stage. For model selection, set \code{count=TRUE}. Symmetric Dirichlet priors are used.
#'
#' @param convData Dataset produced as output by the \code{convData} function. 
#' 
#' @param L Number of higher-level mixture components. When L=1, single-level Latent Class multiple imputation is performed.
#'
#' @param K Number of Latent Classes at the lower-level. 
#' 
#' @param it1 Number of Gibbs sampler iterations for the burn-in (must be larger than 0).
#'
#' @param it2 Number of Gibbs sampler iterations for the imputations.
#' 
#' @param it3 Every \code{it3} iterations, the sampler stores new parameter estimates for the calculation of psoterior estimates. Meaningful only when  
#'    \code{estimates=TRUE}. 
#'
#' @param it.print Every \code{it.print} iterations, the state of the Gibbs sampler is screen-printed.    
#'
#' @param v The Gibbs sampler will produce the first set of imputations at the iteration number \code{(it1+V)}, where V~Unif(1,v). Subsequent imputations are
#'     automatically spaced from each others across the remaining iterations, so that the last imputation (imputation I) occurs at the iteration number 
#'     \code{(it1+it2)}.
#'
#' @param I Number of imputations to be performed. 
#'
#' @param pri2 Hyperparameter value for the higher-level mixture probabilities. Default to 1. 
#'
#' @param pri1 Hyperparameter value for the lower-level mixture probabilities. Default to 1.
#'
#' @param priresp Hyperparameter value for the lower-level conditional response probabilities. Default to 1.
#'
#' @param priresp2 Hyperparameter value for the higher-level conditional response probabilities. Default to 1. 
#'
#' @param random Logical. Should the model parameters be initialized at random values? If \code{TRUE}, parameters are initialized through draws from 
#'     uniform Dirichlet distributions. If \code{FALSE}, parameters are initialized to be equal to 1/D, with D the number of categories in the (observed/latent) 
#'     variable of interest.
#'
#' @param estimates Logical. If \code{TRUE}, the function returns the posterior means of the model parameters. Default to \code{TRUE}.
#'
#' @param count Logical. Should the output include the posterior distribution of L and K? (only K if L=1)
#'
#' @param plot.loglik Logical. Should the output include a traceplot of the log-likelihood ratios obtained through the Gibbs sampler iterations?
#'     Helpful for assessing convergence.
#'
#' @param prec When \code{estimates=TRUE}, \code{prec} defines the number of digits with which the estimates are returned.
#'
#' @param scale Re-scale the log-likelihood value by a factor equal to \code{scale}; this parameter is useful to avoid underflow in the calculation
#'     of the log-likelihood (which can occur in large datasets) and consequently to prevent error messages for the visualization of the log-likelihood traceplot. 
#'    This parameter is meaningful only when \code{plot.loglik} is set to TRUE. Default value equal to 1.0.
#' 
#' @param x A \code{multilevelLCMI} object (\code{print} method). 
#'
#' @param ... Not used. 
#' 
#' @return A \code{multilevelLCMI} object, a list containing: 
#' \item{imp}{Set of imputations for the level-1 variables and units.}
#' \item{imp2}{Set of imputations for the level-2 variables and units.}
#' \item{piL}{posterior means of the level-2 class probabilities. Calculated only if \code{estimates = TRUE}.}
#' \item{piLses}{Posterior standard deviations of the level-2 class probabilities. Calculated only if \code{estimates = TRUE}.}
#' \item{piK}{Posterior means of the level-1 class probabilities. Calculated only if \code{estimates = TRUE}.}
#' \item{piKses}{Posterior standard deviations of the level-2 class probabilities. Calculated only if \code{estimates = TRUE}.}
#' \item{picondlev1}{Posterior means of the level-1 conditional probabilities. Calculated only if \code{estimates = TRUE}.}
#' \item{picondlev1ses}{Posterior standard deviations of the level-1 conditional probabilities. Calculated only if \code{estimates = TRUE}.}
#' \item{picondlev2}{Posterior means of the level-2 conditional probabilities. Calculated only if \code{estimates = TRUE}.}
#' \item{picondlev1ses}{Posterior standard deviations of the level-2 conditional probabilities. Calculated only if \code{estimates = TRUE}.}
#' \item{DIC}{DIC index for the BMLC model. Calculated only if \code{estimates = TRUE.}}
#' \item{freqL}{Posterior distribution of the number of latent classes at level-2. Calculated only if \code{count} is set to \code{TRUE}.}
#' \item{freqK}{Posterior distribution of the number of latent classes at level-1. Calculated only if \code{count} is set to \code{TRUE}.}
#' \item{time}{Running time of the Gibbs sampler iterations.}
#'
#' @references
#' \itemize{
#' 
#' \item[1] Vidotto D., Vermunt J.K., Van Deun K. (2018). 'Bayesian Multilevel Latent Class Models for the Multiple Imputation of Nested Categorical Data'. Journal 
#' of Educational and Beahvioral Statistics 43(5), 511-539. 
#'
#' }
#'
#' @author D. Vidotto <d.vidotto@@uvt.nl> 
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
#'
 

#' @export  
multilevelLCMI <- function(convData, L, K, it1, it2, it3, it.print, v, I = 5, pri2 = 1, pri1 = 1, priresp = 1, priresp2 = 1, random = TRUE, 
    estimates = TRUE, count = FALSE, plot.loglik = FALSE, prec = 3, scale = 1) {
    
    UseMethod("multilevelLCMI")
    
}



#' @export 
multilevelLCMI.default <- function(convData, L, K, it1, it2, it3, it.print, v, I = 5, pri2 = 1, pri1 = 1, priresp = 1, priresp2 = 1, random = TRUE, 
    estimates = TRUE, count = FALSE, plot.loglik = FALSE, prec = 3, scale = 1) {
    
    
    if (class(convData) != "convData") {
        stop("convData object not provided.")
    }
    
    
    # If the hyperparameter values are too low, posterior estimates can become unstable
    if ((pri2 < 0.01 | pri1 < 0.01 | priresp < 0.01 | (priresp2 < 0.01 & !(convData$doVar2))) & estimates == TRUE) {
        warning("some hyperparameters value is very low (<.01). This might cause some NaN in the final estimates.", call. = FALSE, immediate. = TRUE)
    }
    
    # L and K must be at least equal to 1
    if (L < 1 | K < 1) {
        stop("L or K are <1.", call. = FALSE)
    }
    
    
    # The number of burn-in iterations must be larger than 0
    if (it1 == 0) {
        stop("Set the number of iterations <<it1>>.", call. = FALSE)
    }
    
    doVar2 <- convData$doVar2  # Boolean; are level-2 variables present in the dataset? 
    GID <- 1  # Group ID column index 
    var2 <- convData$var2  # level-2 variable indices; it can be NULL
    if (!is.null(convData$UID)) {
        var2 <- var2 - 1
    }
    dat <- convData$convDat  # The  dataset converted by convData 
    
    
    # Select the iterations in which the imputations will occur
    sel <- floor(seq(floor(runif(1, 1, v)), it2, length.out = I))
    
    
    Y <- as.data.frame(dat)  # Original dataset 
    n <- nrow(dat)  # (Total) Sample size 
    GID <- 1
    J <- length(table(Y[, GID]))  # Number of groups 
    
    
    
    # Vector that stores log-likelihood values
    pllk <- rep(0, 1)
    if (plot.loglik) {
        pllk <- rep(0, (it1 + it2))
    }
    
    # Preparation objects for level-2 variables imputation
    TT2 <- 1  # Number of level-2 variables
    ncat2 <- rep(0, 1)  # Number of categories for each level-2 variable  
    cod2 <- 0  # Codes for variables at level-2
    nms2 <- 0  # Names of level-2 variables 
    Y2 <- matrix(0, 1, 1)  # Level-2 dataset 
    R2 <- matrix(0, 1, 1)  # Level-2 missingess indicator 
    YY2 <- matrix(0, 1, 1)  # Auxiliary matrix 
    nlat2 <- vector("list", 1)  # Level-2 conditional counts
    piresp2 <- vector("list", 1)  # Level-2 conditional probabilities
    f = it3/it2  # Adjustment factor for posterior means 
    
    
    
    GU <- as.numeric(as.factor(Y[, GID]))  # Group ID recoded 
    
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
        for (tt in 1:TT2) {
            nlat2[[tt]] <- matrix(0, L, ncat2[tt])
            piresp2[[tt]] <- matrix(1/ncat2[tt], L, ncat2[tt])
            if (random) {
                for (l in 1:L) {
                  piresp2[[tt]][l, ] <- rgamma(ncat2[[tt]], 1, 1)
                  piresp2[[tt]][l, ] <- piresp2[[tt]][l, ]/sum(piresp2[[tt]][l, ])
                }
            }
        }
        YY2 <- as.matrix(Y2[, -1])
    }
    
    # Group-specific sample size
    if (J > 1) {
        nj <- table(GU)
    } else {
        nj <- rep(0, 1)
        nj[1] <- n
    }
    
    # Isolate level-1 data from ID and level-2 variables
    if (doVar2 == FALSE) {
        YY <- Y[, -c(GID)]
        colnames(YY) <- colnames(Y[, -GID])
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
    rm(YY)
    
    datJ <- vector("list", J)  # Level-1 group-specific data 
    R <- vector("list", J)  # Level-1 group-specific missingness indicator 
    npattern <- vector("list", J)  # Group-specific number of patterns observed at level-1 
    npattern2 <- rep(0, J)  # Number of patterns observed at level-2 
    pat2 <- rep(-1, J)  # Patterns observed at level-2
    pat <- vector("list", J)  # Group-specific patterns observed at level-2 
    w <- rep(-1, J)  # Level-2 latent class indicator 
    z <- vector("list", J)  # Group-specific level-1 latent class indicator  
    nresp <- vector("list", TT)  # Level-1 conditional probabilities counts  
    piresp <- vector("list", TT)  # Level-1 conditional probabilities 
    
    pclasslev1 <- vector("list", J)  # Group-specific level-1 posterior memberships
    jointlev2 <- matrix(0, L, J)  # Level-2 posterior memberships 
    grouploglik <- rep(0, J)  # Group-specific log-likelihood 
    loglik <- 0  # Initialization log-likelihood 
    nlatw <- rep(0, L)  # Level-2 latent class counts   
    nlatk <- matrix(0, L, K)  # Level-1 latent class counts 
    piK <- matrix(1/K, L, K)  # Level-1 latent class probabilities  
    piW <- rep(1/L, L)  # Level-2 latent class probabilities 
    
    
    if (random == TRUE) {
        piW <- rgamma(L, 1, 1)
        piW <- piW/sum(piW)
        for (l in 1:L) {
            nlatk[l, ] <- rep(0, K)
            piK[l, ] <- rgamma(K, 1, 1)
            piK[l, ] <- piK[l, ]/(sum(piK[l, ]))
        }
        for (tt in 1:TT) {
            piresp[[tt]] <- vector("list", L)
            nresp[[tt]] <- vector("list", L)
            for (l in 1:L) {
                nresp[[tt]][[l]] <- matrix(0, K, ncat[[tt]])
                piresp[[tt]][[l]] <- matrix(0, K, ncat[tt])
                for (k in 1:K) {
                  piresp[[tt]][[l]][k, ] <- rgamma(ncat[[tt]], 1, 1)
                  piresp[[tt]][[l]][k, ] <- piresp[[tt]][[l]][k, ]/sum(piresp[[tt]][[l]][k, ])
                }
            }
        }
    } else {
        for (tt in 1:TT) {
            piresp[[tt]] <- vector("list", L)
            nresp[[tt]] <- vector("list", L)
            for (l in 1:L) {
                piresp[[tt]][[l]] <- matrix(1/ncat[[tt]], K, ncat[[tt]])
                nresp[[tt]][[l]] <- matrix(0, K, ncat[[tt]])
            }
        }
    }
    
    
    imp <- vector("list", J)  # Group-specific level-1 imputation set
    nimp <- 0  # Imputation counter 
    imp2 <- vector("list", 1)  # Level-2 imputation set
    
    if (doVar2) {
        imp2 <- vector("list", TT2)
        for (tt2 in 1:TT2) {
            if (sum(R2[, tt2] > 0)) {
                if (I > 0) {
                  imp2[[tt2]] <- matrix(NA, sum(R2[, tt2]), I)
                }
            }
        }
    }
    
    condp1 <- matrix(0, L, J)  # Level-1 conditional probabilities given the level-2 latent classes
    
    
    for (j in 1:J) {
        datJ[[j]] <- matrix(0, nj[j], TT)
        R[[j]] <- matrix(0, nj[j], TT)
        datJ[[j]] <- Y[GU == j, 2:(TT + 1)]
        colnames(datJ[[j]]) <- nms
        R[[j]] <- ifelse(is.na(datJ[[j]]), 1, 0)
        datJ[[j]] <- as.matrix(datJ[[j]])
        R[[j]] <- as.matrix(R[[j]])
        npattern[[j]] <- rep(0, nj[j])
        pat[[j]] <- rep(-1, nj[j])
        z[[j]] <- rep(-1, nj[j])
        pclasslev1[[j]] <- vector("list", L)
        for (l in 1:L) {
            pclasslev1[[j]][[l]] <- matrix(0, K, nj[j])
        }
        imp[[j]] <- vector("list", TT)
        for (tt in 1:TT) {
            if (sum(R[[j]][, tt] > 0)) {
                if (I > 0) {
                  imp[[j]][[tt]] <- matrix(NA, sum(R[[j]][, tt]), I)
                  colnames(imp[[j]][[tt]]) <- 1:I
                  rownames(imp[[j]][[tt]]) <- which(R[[j]][, tt] == 1)
                }
            }
        }
    }
    
    
    # Posterior means and standard deviations estimates
    piWm <- rep(0, L)  # Level-2 latent classes means 
    piWs <- rep(0, L)  # Level-2 latent classes standard deviations 
    piKm <- matrix(0, L, K)  # Level-1 latent classes means 
    piKs <- matrix(0, L, K)  # Level-1 latent classes standard deviations 
    pirespm <- vector("list", TT)  # Level-1 conditional probabilities means 
    piresps <- vector("list", TT)  # Level-1 conditional probabilities standard deviations 
    piresp2m <- vector("list", TT2)  # Level-2 conditional probabilities means 
    piresp2s <- vector("list", TT2)  # Level-2 conditional probabilities standard deviations
    
    
    
    for (tt in 1:TT) {
        pirespm[[tt]] <- piresps[[tt]] <- vector("list", L)
        for (l in 1:L) {
            pirespm[[tt]][[l]] <- matrix(0, K, ncat[tt])
            piresps[[tt]][[l]] <- matrix(0, K, ncat[tt])
        }
    }
    
    
    if (doVar2) {
        for (tt in 1:TT2) {
            piresp2m[[tt]] <- matrix(0, L, ncat2[tt])
            piresp2s[[tt]] <- matrix(0, L, ncat2[tt])
        }
    }
    
    # Running time
    ptm <- proc.time()
    
    # Run the Gibbs Sampler
    mc <- cppCycle(L, K, J, it1, it2, it3, it.print, pri2, pri1, priresp, priresp2, estimates, doVar2, sel, TT2, ncat2, YY2, R2, nlat2, piresp2, 
        nj, TT, ncat, datJ, R, npattern, pat, w, z, nresp, piresp, pclasslev1, jointlev2, grouploglik, loglik, nlatw, nlatk, piK, piW, imp, 
        imp2, nimp, count, piWm, piKm, pirespm, piresp2m, piWs, piKs, piresps, piresp2s, f, plot.loglik, pllk, condp1, npattern2, pat2, scale)
    
    ptmf <- proc.time() - ptm
    
    # Store the output of the Gibbs sampler
    piWm <- mc$piWm
    piKm <- mc$piKm
    pirespm <- mc$pirespm
    DIC <- mc$DIC
    imp <- mc$imp
    piWs <- mc$piWs
    piKs <- mc$piKs
    piresps <- mc$piresps
    countIT <- mc$countIT
    
    
    if (doVar2) {
        imp2 <- mc$imp2
        piresp2m <- mc$piresp2m
        piresp2s <- mc$piresp2s
    } else {
        imp2 <- NULL
        piresp2m <- piresp2s <- NULL
    }
    
    # Posterior means and standard deviations
    if (!estimates) {
        DIC <- piWm <- piKm <- pirespm <- piWs <- piKs <- piresps <- piresp2m <- piresp2s <- NULL
    } else {
        names(piWm) <- names(piWs) <- noquote(paste("W=", 1:L, sep = ""))
        rownames(piKm) <- rownames(piKs) <- noquote(paste("W=", 1:L, sep = ""))
        colnames(piKm) <- colnames(piKs) <- noquote(paste("X=", 1:K, sep = ""))
        piWm <- round(piWm, prec)
        piKm <- round(piKm, prec)
        piWs <- round(piWs, prec)
        piKs <- round(piKs, prec)
        
        for (tt in 1:TT) {
            names(pirespm)[tt] <- names(piresps)[tt] <- nms[tt]
            for (l in 1:L) {
                names(pirespm[[tt]])[l] <- names(piresps[[tt]])[l] <- noquote(paste("W=", l, sep = ""))
                rownames(pirespm[[tt]][[l]]) <- rownames(piresps[[tt]][[l]]) <- noquote(paste("X=", 1:K, sep = ""))
                colnames(pirespm[[tt]][[l]]) <- colnames(piresps[[tt]][[l]]) <- noquote(cod[[tt]][1, ])
            }
        }
        
        pirespm <- lapply(pirespm, function(x) lapply(x, function(y) round(y, prec)))
        piresps <- lapply(piresps, function(x) lapply(x, function(y) round(y, prec)))
        pirespm <- lapply(pirespm, function(x) lapply(x, function(y) t(y)))
        piresps <- lapply(piresps, function(x) lapply(x, function(y) t(y)))
        
        if (doVar2) {
            for (tt in 1:TT2) {
                names(piresp2m)[tt] <- names(piresp2s)[tt] <- nms2[tt]
                rownames(piresp2m[[tt]]) <- rownames(piresp2s[[tt]]) <- noquote(paste("W=", 1:L, sep = ""))
                colnames(piresp2m[[tt]]) <- colnames(piresp2s[[tt]]) <- noquote(cod2[[tt]][1, ])
            }
            piresp2m <- lapply(piresp2m, function(x) round(x, prec))
            piresp2s <- lapply(piresp2s, function(x) round(x, prec))
            piresp2m <- lapply(piresp2m, function(y) t(y))
            piresp2s <- lapply(piresp2s, function(y) t(y))
        }
    }
    
    # Create imputation sets (level-2 and level-1 )
    if (doVar2) {
        
        if (it2 > 0 & I > 0) {
            
            imp2_ <- vector("list", TT2)
            for (tt2 in 1:TT2) {
                if (sum(R2[, tt2]) > 0) {
                  imp2_[[tt2]] <- matrix(NA, sum(R2[, tt2]), I)
                  for (i in 1:I) {
                    for (m2 in 1:ncat2[tt2]) {
                      imp2_[[tt2]][which(imp2[[tt2]][, i] == cod2[[tt2]][2, m2]), i] <- as.numeric(noquote(cod2[[tt2]][1, m2]))
                    }
                  }
                }
            }
            
            imp2 <- imp2_
            rm(imp2_)
            names(imp2) <- nms2
            for (tt2 in 1:TT2) {
                if (sum(R2[, tt2] > 0)) {
                  if (I > 0) {
                    colnames(imp2[[tt2]]) <- 1:I
                    rownames(imp2[[tt2]]) <- which(R2[, tt2] == 1)
                  }
                }
            }
            
            imp22 <- vector("list", J)
            for (j in 1:J) {
                imp22[[j]] <- vector("list", TT)
                for (tt in 1:TT) {
                  if (sum(R[[j]][, tt]) > 0) {
                    
                    imp22[[j]][[tt]] <- matrix(NA, sum(R[[j]][, tt]), I)
                    
                    for (i in 1:I) {
                      for (m in 1:ncat[tt]) {
                        imp22[[j]][[tt]][which(imp[[j]][[tt]][, i] == cod[[tt]][2, m]), i] <- as.numeric(noquote(cod[[tt]][1, m]))
                      }
                    }
                  }
                }
            }
            imp <- imp22
            rm(imp22)
            
            for (j in 1:J) {
                names(imp)[j] <- noquote(paste("Group", j, sep = " "))
                for (tt in 1:TT) {
                  names(imp[[j]])[tt] <- nms[tt]
                  if (sum(R[[j]][, tt]) > 0) {
                    colnames(imp[[j]][[tt]])[1:I] <- 1:I
                    rownames(imp[[j]][[tt]]) <- which(R[[j]][, tt] == 1)
                  }
                }
            }
        }
        # Create imputation sets (only level-1 )
    } else {
        
        if (it2 > 0 && I > 0) {
            
            imp22 <- vector("list", J)
            
            for (j in 1:J) {
                imp22[[j]] <- vector("list", TT)
                
                for (tt in 1:TT) {
                  if (sum(R[[j]][, tt]) > 0) {
                    
                    imp22[[j]][[tt]] <- matrix(NA, sum(R[[j]][, tt]), I)
                    
                    
                    for (i in 1:I) {
                      for (m in 1:ncat[tt]) {
                        imp22[[j]][[tt]][which(imp[[j]][[tt]][, i] == cod[[tt]][2, m]), i] <- as.numeric(noquote(cod[[tt]][1, m]))
                      }
                    }
                  }
                }
            }
            
            imp <- imp22
            rm(imp22)
            
            for (j in 1:J) {
                names(imp)[j] <- noquote(paste("Group", j, sep = " "))
                for (tt in 1:TT) {
                  names(imp[[j]])[tt] <- nms[tt]
                  if (sum(R[[j]][, tt]) > 0) {
                    colnames(imp[[j]][[tt]])[1:I] <- 1:I
                    rownames(imp[[j]][[tt]]) <- which(R[[j]][, tt] == 1)
                  }
                }
            }
        }
    }
    
    if (it2 == 0 || I == 0) {
        imp = "None"
        imp2 = "None"
    }
    
    
    # Plots of log-likelihoods and classes distribution
    
    if (plot.loglik & !count) {
        maxlik <- mc$maxlik
        pllk <- mc$pllk
        pllk <- maxlik - pllk
        plot(pllk, type = "l", lwd = 2, xlab = "t", ylab = "ll(max)-ll", main = "log-likelihood ratios", ylim = c(-0.1, max(pllk, na.rm = TRUE)))
        countK <- NULL
        countL <- NULL
    } else if (!plot.loglik & count) {
        countK <- mc$countK
        countIT[countIT == 0] <- 1
        countK <- countK/countIT
        countL <- mc$countL
        countL <- countL/it2
        names(countL) <- paste("L=", 1:L, sep = "")
        rownames(countK) <- paste("W=", 1:L, sep = "")
        colnames(countK) <- paste("K=", 1:K, sep = "")
        if (L > 1) {
            par(mfrow = c(2, 1))
            plot(countL, type = "h", lwd = "3", col = ifelse(countL == max(countL), "red", "black"), xlab = "Level 2 Class", ylab = "Freq.", 
                main = "Distribution of L", xaxt = "n")
            axis(side = 1, at = 1:L, labels = paste("L=", 1:L, sep = ""))
            cols <- sample(1:657, L, replace = F)
            matplot(t(countK), type = "s", lwd = 3, ylim = c(0, 1.5), lty = 1, xlab = "Level 1 Class", ylab = "Freq.", main = "Distribution of K|W", 
                yaxt = "n", col = colors()[cols], xaxt = "n")
            if (L > 2) {
                legend(1, 1.5, legend = paste("W = ", 1:L, sep = ""), fill = colors()[cols], ncol = floor(L/3))
            } else {
                legend(1, 1.5, legend = paste("W = ", 1:L, sep = ""), fill = colors()[cols])
            }
            axis(side = 1, at = 1:K, labels = paste("K = ", 1:K, sep = ""))
            axis(side = 2, at = seq(0, 1, by = 0.1), labels = seq(0, 1, by = 0.1))
            par(mfrow = c(1, 1))
        } else {
            plot(1:K, countK, type = "h", lwd = "3", col = ifelse(countK == max(countK), "red", "black"), xlab = "Level 1 Class", ylab = "Freq.", 
                main = "Distribution of K", xaxt = "n", ylim = c(0, 1))
            axis(side = 1, at = 1:K, labels = paste("K=", 1:K, sep = ""))
        }
    } else if (plot.loglik & count) {
        maxlik <- mc$maxlik
        pllk <- mc$pllk
        pllk <- maxlik - pllk
        countK <- mc$countK
        countIT[countIT == 0] <- 1
        countK <- countK/countIT
        rownames(countK) <- paste("W=", 1:L, sep = "")
        colnames(countK) <- paste("K=", 1:K, sep = "")
        countL <- mc$countL
        countL <- countL/it2
        names(countL) <- paste("L=", 1:L, sep = "")
        if (L > 1) {
            par(mfrow = c(3, 1))
            plot(pllk, type = "l", lwd = 2, xlab = "t", ylab = "ll(max)-ll", main = "log-likelihood ratios", ylim = c(-0.1, max(pllk, na.rm = TRUE)))
            plot(countL, type = "h", lwd = "3", col = ifelse(countL == max(countL), "red", "black"), xlab = "Level 2 Class", ylab = "Freq.", 
                main = "Distribution of L", xaxt = "n")
            axis(side = 1, at = 1:L, labels = paste("L=", 1:L, sep = ""))
            cols <- sample(1:657, L, replace = F)
            matplot(t(countK), type = "s", lwd = 3, ylim = c(0, 1.5), lty = 1, xlab = "Level 1 Class", ylab = "Freq.", main = "Distribution of K|W", 
                yaxt = "n", col = colors()[cols], xaxt = "n")
            if (L > 2) {
                legend(1, 1.5, legend = paste("W = ", 1:L, sep = ""), fill = colors()[cols], ncol = floor(L/3))
            } else {
                legend(1, 1.5, legend = paste("W = ", 1:L, sep = ""), fill = colors()[cols])
            }
            axis(side = 1, at = 1:K, labels = paste("K = ", 1:K, sep = ""))
            axis(side = 2, at = seq(0, 1, by = 0.1), labels = seq(0, 1, by = 0.1))
            par(mfrow = c(1, 1))
        } else {
            par(mfrow = c(2, 1))
            plot(pllk, type = "l", lwd = 2, xlab = "t", ylab = "ll(max)-ll", main = "log-likelihood ratios", ylim = c(-0.1, max(pllk, na.rm = TRUE)))
            plot(1:K, countK, type = "h", lwd = "3", col = ifelse(countK == max(countK), "red", "black"), xlab = "Level 1 Class", ylab = "Freq.", 
                main = "Distribution of K", xaxt = "n", ylim = c(0, 1))
            axis(side = 1, at = 1:K, labels = paste("K=", 1:K, sep = ""))
            par(mfrow = c(1, 1))
        }
    } else {
        countK <- NULL
        countL <- NULL
        countG <- NULL
    }
    
    
    
    ret <- list(implev1 = imp, implev2 = imp2, piL = piWm, piLses = piWs, piK = piKm, piKses = piKs, picondlev1 = pirespm, picondlev1ses = piresps, 
        picondlev2 = piresp2m, picondlev2ses = piresp2s, DIC = DIC, freqL = countL, freqK = countK, time = ptmf)
    
    ret$call <- match.call()
    
    class(ret) <- "multilevelLCMI"
    class(ret$implev1) <- "impmlcmi"
    class(ret$implev2) <- "impmlcmi"

	ret 

}





#' @export 
print.multilevelLCMI <- function(x, ... ){
	
    cat("Call:\n")
    print(x$call)
    cat("\n")
    cat("Elpsed time: ", x$time[[3]], "\n")	

}



