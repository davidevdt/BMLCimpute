# Package BMLCimpute
# Function multilevelLCMI
# 
# Function for performing Multiple Imputation with the Bayesian Multilevel Latent Class Model. The model takes as input the list produced by the 'convData' function, in which the dataset converted and prepared for the imputations is present, along with other parameters specified by the user (e.g., number of latent classes and specification of the prior distribution hyperparameters). The function can also offer (when the corresponding boolean parameter is activated) a graphical representation of the posterior distribution of the number of occupied classes during the Gibbs sampler iterations. In this way, 'multilevelLCMI' can also perform model selection in a pre-imputation stage. Symmetric Dirichlet priors are used. 
#
# Inputs :   
# 
# @param convData Dataset produced as output by the 'convData' function 
# @param L Number of level-2 latent classes 
# @param K Number of level-1 latent classes 
# @param it1 Number of burn-in iterations for the Gibbs Sampler 
# @param it2 Number of Gibbs sampler iterations for the imputations, and more in general for the computation of the posterior distribution (e.g., of the number of allocated classes) 
# @param it3 If estimates = TRUE, the sampled estimates of the model parameters (i.e., class and conditional probabilities; NOT the number of allocated classes!) are stored every 'it3' iterations for the computation of the posterior means
# @param it.print Output the current status of the algorithm (iteration number; produced log-likelihood) every 'it.print' iterations 
# @param v Perform the first imputation after it1 + v_ iterations, where v_ ~ floor(Unif(1, v))
# @param I Total number of imputations to be performed (integer) (default 5)
# @param pri2 Hyperparameter of the (symmetric) Dirichlet distribution for the level-2 latent class probabilities (default 1.0)
# @param pri1 Hyperparameter of the (symmetric) Dirichlet distribution for the level-1 latent class probabilities (default 1.0)
# @param priresp Hyperparameter of the (symmetric) Dirichlet distribution for the level-1 conditional probabilities (default 1.0)
# @param priresp2 Hyperparameter of the (symmetric) Dirichlet distribution for the level-2 conditional probabilities (default 1.0)
# @param random Boolean. If TRUE, initialize the model parameters by drawing them from a Dirichlet distribution with all hyperparameters set equal to 1. If FALSE, parameters are set to 1/D, where D is the number of categories in the (observed/latent) variable of interest. Default to TRUE 
# @param estimates Boolean. If TRUE, the function returns the posterior means of the model parameters. Default to TRUE
# @param count Boolean. If TRUE, both graphical and numerical representations of the posterior distribution for the number of classes allocated (at both levels) during the Gibbs sampler iterations are returned. Default to FALSE
# @param plot.loglik Boolean. If TRUE, a traceplot of the log-likelihood calculated at each iteration of the Gibbs sampler (including burn-in) is returned. Useful to check convergence. Defaut to FALSE
# @param prec Number of digits for the decimal part of the returned posterior probabilities; only meaningful when 'estimates' is set to TRUE. Equal to 3 by default
# @param scale Re-scale the log-likelihood value by a factor equal to 'scale'; this parameter is useful to avoid underflow in the calculation of the log-likelihood (which can occur in large datasets) and consequently to prevent error messages for the visualization of the log-likelihood traceplot. This parameter is meaningful only when 'plot.loglik' is set to TRUE. Default value equal to 1.0   
#
# Outputs : 
# 
# @return implev1 Set of imputations for the level-1 variables and units
# @return implev2 Set of imputations for the level-2 variables and units
# @return piL Posterior means of the level-2 class probabilities. Calculated only if estimates = TRUE
# @return piLses Posterior standard deviations of the level-2 class probabilities. Calculated only if estimates = TRUE
# @return piK Posterior means of the level-1 class probabilities. Calculated only if estimates = TRUE
# @return piKses Posterior standard deviations of the level-2 class probabilities. Calculated only if estimates = TRUE
# @return picondlev1 Posterior means of the level-1 conditional probabilities. Calculated only if estimates = TRUE
# @return picondlev1ses Posterior standard deviations of the level-1 conditional probabilities. Calculated only if estimates = TRUE
# @return picondlev2 Posterior means of the level-2 conditional probabilities. Calculated only if estimates = TRUE
# @return picondlev2ses Posterior standard deviations of the level-2 conditional probabilities. Calculated only if estimates = TRUE
# @return DIC DIC index for the BMLC model. Calculated only if estimates = TRUE
# @return freqL Posterior distribution of the number of latent classes at level-2. Calculated only if count is set to TRUE
# @return freqK Posterior distribution of the number of latent classes at level-1. Calculated only if count is set to TRUE
# @return time Running time of the Gibbs sampler iterations       


multilevelLCMI <- 
function( convData, L, K, it1, it2, it3, it.print, v, I = 5, pri2 = 1.0, pri1 = 1.0, priresp = 1.0, priresp2 = 1.0, random = TRUE, estimates = TRUE, count = FALSE, plot.loglik = FALSE, prec = 3, scale = 1.0){

	
	# If the hyperparameter values are too low, posterior estimates can become unstable 
	if( ( pri2 < 0.01 | pri1 < 0.01 | priresp < 0.01 | (priresp2 < 0.01 & !(convData$doVar2)) ) & estimates == TRUE){
		warning("some hyperparameters value is very low (<.01). This might cause some NaN in the final estimates.", call. = FALSE, immediate. = TRUE)
	}

	# L and K must be at least equal to 1 
	if( L < 1 | K < 1 ){
		stop("L or K are <1.", call. = FALSE)
	}

	
	# The number of burn-in iterations must be larger than 0 
	if( it1 == 0 ){
		stop("Set the number of iterations <<it1>>.", call. = FALSE)
	}

	doVar2 <- convData$doVar2				# Boolean; are level-2 variables present in the dataset? 	
	GID <- 1								# Group ID column index 
	var2 <- convData$var2					# level-2 variable indices; it can be NULL
	if( !is.null(convData$UID) ){
		var2 <- var2 - 1
	}
	dat <- convData$convDat					# The  dataset converted by convData 

	
	# Select the iterations in which the imputations will occur 
	sel <- floor(seq(floor(runif(1, 1, v)), it2, length.out = I))

	
	Y <- as.data.frame(dat)								# Original dataset 
	n <- nrow(dat)										# (Total) Sample size 
	GID <- 1
	J <- length(table(Y[ ,GID]))						# Number of groups 
	
	

	# Vector that stores log-likelihood values 
	pllk <- rep(0, 1)										
	if( plot.loglik ){
		pllk <- rep(0, (it1 + it2))
	}

	# Preparation objects for level-2 variables imputation 
	TT2<-1												# Number of level-2 variables
	ncat2<-rep(0,1)										# Number of categories for each level-2 variable  
	cod2<-0												# Codes for variables at level-2
	nms2<-0												# Names of level-2 variables 
	Y2<-matrix(0,1,1)									# Level-2 dataset 
	R2<-matrix(0,1,1)									# Level-2 missingess indicator 
	YY2<-matrix(0,1,1)									# Auxiliary matrix 
	nlat2<-vector("list",1)								# Level-2 conditional counts
	piresp2<-vector("list",1)							# Level-2 conditional probabilities
	f = it3/it2											# Adjustment factor for posterior means 



	GU<-as.numeric(as.factor(Y[,GID]))					# Group ID recoded 

	# If variables at level-2 are present, initialize the objects just created  
	if( doVar2 ){
		Y2 <- cbind(GU, Y[ ,var2])
		Y2 <- unique(Y2)
		Y2 <- Y2[order(Y2[ ,1]),]
		Y2b <- Y2
		R2 <- ifelse(is.na(Y2[ ,-1]), 1, 0)
		nms2 <- convData$namesLev2
		TT2 <- ncol(Y2[ ,-1])
		ncat2 <- convData$nCatLev2
		cod2 <- convData$codLev2
		piresp2 <- vector("list", TT2)
		nlat2 <- vector("list", TT2)
		for( tt in 1:TT2 ){
			nlat2[[tt]] <- matrix(0, L, ncat2[tt])
			piresp2[[tt]] <- matrix(1 / ncat2[tt], L, ncat2[tt])
			if( random ){
				for( l in 1:L ){
					piresp2[[tt]][l, ] <- rgamma(ncat2[[tt]], 1, 1)
					piresp2[[tt]][l, ] <- piresp2[[tt]][l, ]/sum(piresp2[[tt]][l, ])
				}
			}
		}
		YY2<-as.matrix(Y2[,-1])
	}

	# Group-specific sample size 
	if(J>1){
		nj<-table(GU)
	} else {
		nj <- rep(0,1)
		nj[1] <- n
	}

	# Isolate level-1 data from ID and level-2 variables  
	if( doVar2 == FALSE ){
		YY <- Y[ ,-c(GID)]
		colnames(YY) <- colnames(Y[ ,-GID])
	} else {
		YY <- Y[ ,-c(GID, var2)]
		colnames(YY) <- colnames(Y[ ,-c(GID, var2)])
	}
	nms <- colnames(YY)
	
	# Level-1 imputation : Initializations 
	
	TT <- ncol(YY)										# Total number of level-1 variables 
	ncat <- convData$nCatLev1							# Number of categories for each variable  
	cod <- convData$codLev1								# Codes for the categories (original/new) 
	Y <- cbind(GU, YY)									# Level-1 dataset
	rm(YY)

	datJ <- vector("list", J)							# Level-1 group-specific data 
	R <- vector("list", J)								# Level-1 group-specific missingness indicator 
	npattern <- vector("list", J)						# Group-specific number of patterns observed at level-1 
	npattern2 <- rep(0, J)								# Number of patterns observed at level-2 
	pat2 <- rep(-1, J) 									# Patterns observed at level-2
	pat <- vector("list", J)							# Group-specific patterns observed at level-2 
	w <- rep(-1, J)										# Level-2 latent class indicator 
	z <- vector("list", J)								# Group-specific level-1 latent class indicator  
	nresp <- vector("list", TT)							# Level-1 conditional probabilities counts  
	piresp <- vector("list", TT)						# Level-1 conditional probabilities 

	pclasslev1 <- vector("list", J)						# Group-specific level-1 posterior memberships
	jointlev2 <- matrix(0.0, L, J)						# Level-2 posterior memberships 
	grouploglik <- rep(0.0, J)							# Group-specific log-likelihood 
	loglik <- 0.0										# Initialization log-likelihood 
	nlatw <- rep(0, L)									# Level-2 latent class counts   
	nlatk <- matrix(0, L, K)							# Level-1 latent class counts 
	piK <- matrix(1/K, L, K)							# Level-1 latent class probabilities  
	piW <- rep(1/L,L)									# Level-2 latent class probabilities 	

	
	if( random == TRUE ){
		piW <- rgamma(L, 1, 1)
		piW <- piW / sum(piW)
		for( l in 1:L ){
			nlatk[l, ] <- rep(0, K)
			piK[l, ] <- rgamma(K, 1, 1)
			piK[l, ] <- piK[l, ] / (sum(piK[l, ]))
		}
		for( tt in 1:TT ){
			piresp[[tt]] <- vector("list", L)
			nresp[[tt]] <- vector("list", L)
			for( l in 1:L ){
				nresp[[tt]][[l]] <- matrix(0, K, ncat[[tt]])
				piresp[[tt]][[l]] <- matrix(0, K, ncat[tt])
				for( k in 1:K ){
					piresp[[tt]][[l]][k, ] <- rgamma(ncat[[tt]], 1, 1)
					piresp[[tt]][[l]][k, ] <- piresp[[tt]][[l]][k, ]/sum(piresp[[tt]][[l]][k, ])
				}
			}
		}
	} else {
		for( tt in 1:TT ){
			piresp[[tt]] <- vector("list", L)
			nresp[[tt]] <- vector("list", L)
			for( l in 1:L ){
				piresp[[tt]][[l]] <- matrix(1 / ncat[[tt]], K, ncat[[tt]])
				nresp[[tt]][[l]] <- matrix(0, K, ncat[[tt]])
			}
		}
	}
	
	
	imp <- vector("list", J)							# Group-specific level-1 imputation set		
	nimp <- 0											# Imputation counter 
	imp2 <- vector("list", 1)							# Level-2 imputation set	

	if( doVar2 ){
		imp2 <- vector("list", TT2)
		for( tt2 in 1:TT2 ){
			if( sum(R2[ ,tt2] > 0) ){
				if( I > 0 ){
					imp2[[tt2]] <- matrix(NA, sum(R2[ ,tt2]), I)
				}
			}
		}
	}

	condp1 <- matrix(0.0, L, J)							# Level-1 conditional probabilities given the level-2 latent classes		

	
	for( j in 1:J ){
		datJ[[j]] <- matrix(0, nj[j], TT)
		R[[j]] <- matrix(0, nj[j], TT)
		datJ[[j]] <- Y[GU == j, 2:(TT+1)]
		colnames(datJ[[j]]) <- nms
		R[[j]] <- ifelse(is.na(datJ[[j]]), 1, 0)
		datJ[[j]] <- as.matrix(datJ[[j]])
		R[[j]] <- as.matrix(R[[j]])
		npattern[[j]] <- rep(0, nj[j])
		pat[[j]] <- rep(-1, nj[j])
		z[[j]] <- rep(-1, nj[j])
		pclasslev1[[j]] <- vector("list", L)
		for( l in 1:L ){
			pclasslev1[[j]][[l]] <- matrix(0, K, nj[j])
		}
		imp[[j]] <- vector("list", TT)
		for( tt in 1:TT ){
			if( sum(R[[j]][ ,tt] > 0) ){
				if( I > 0 ){
					imp[[j]][[tt]] <- matrix(NA, sum(R[[j]][ ,tt]), I)
					colnames(imp[[j]][[tt]]) <- 1:I
					rownames(imp[[j]][[tt]]) <- which(R[[j]][ ,tt] == 1)
				}
			}
		}
	}

	
	# Posterior means and standard deviations estimates  
	piWm <- rep(0, L)								# Level-2 latent classes means 							
	piWs <- rep(0, L)								# Level-2 latent classes standard deviations 
	piKm <- matrix(0, L, K)							# Level-1 latent classes means 
	piKs <- matrix(0, L, K)							# Level-1 latent classes standard deviations 
	pirespm <- vector("list", TT)					# Level-1 conditional probabilities means 
	piresps <- vector("list", TT)					# Level-1 conditional probabilities standard deviations 
	piresp2m <- vector("list", TT2)					# Level-2 conditional probabilities means 
	piresp2s <- vector("list", TT2)					# Level-2 conditional probabilities standard deviations

	

	for( tt in 1:TT ){
		pirespm[[tt]] <- piresps[[tt]] <- vector("list", L)
		for( l in 1:L ){
			pirespm[[tt]][[l]] <- matrix(0, K, ncat[tt])
			piresps[[tt]][[l]] <- matrix(0, K, ncat[tt])
		}
	}


	if(doVar2){
		for( tt in 1:TT2 ){
			piresp2m[[tt]] <- matrix(0, L, ncat2[tt])
			piresp2s[[tt]] <- matrix(0, L, ncat2[tt])
		}
	}
	
	# Running time 
	ptm <- proc.time()
	
	# Run the Gibbs Sampler 
	mc<-cppCycle(L, K, J, it1, it2, it3, it.print, pri2, pri1, priresp, priresp2, estimates, doVar2, sel, TT2, ncat2,YY2, R2, nlat2, piresp2, nj, TT, ncat, datJ, R, npattern, pat, w, z, nresp, piresp, pclasslev1, jointlev2, grouploglik, loglik, nlatw, nlatk, piK, piW, imp, imp2, nimp, count, piWm, piKm, pirespm, piresp2m, piWs, piKs, piresps, piresp2s, f, plot.loglik, pllk, condp1, npattern2, pat2, scale)

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


	if(doVar2){
		imp2 <- mc$imp2
		piresp2m <- mc$piresp2m
		piresp2s <- mc$piresp2s
	}else{
		imp2 <- NULL
		piresp2m <- piresp2s <- NULL
	}
	
	# Posterior means and standard deviations 
	if( !estimates ){
		DIC <- piWm <- piKm <- pirespm <- piWs <- piKs <- piresps <- piresp2m <- piresp2s <- NULL
	} else {
		names(piWm) <- names(piWs) <- noquote(paste("W=", 1:L, sep = ""))
		rownames(piKm) <- rownames(piKs) <- noquote(paste("W=", 1:L, sep = ""))
		colnames(piKm) <- colnames(piKs) <- noquote(paste("X=", 1:K, sep = ""))
		piWm <- round(piWm, prec)
		piKm <- round(piKm, prec)
		piWs <- round(piWs, prec)
		piKs <- round(piKs, prec)
		
		for( tt in 1:TT ){
			names(pirespm)[tt] <- names(piresps)[tt] <- nms[tt]
			for( l in 1:L ){
				names(pirespm[[tt]])[l] <- names(piresps[[tt]])[l] <- noquote(paste("W=", l, sep = ""))
				rownames(pirespm[[tt]][[l]]) <- rownames(piresps[[tt]][[l]]) <- noquote(paste("X=", 1:K, sep=""))
				colnames(pirespm[[tt]][[l]]) <- colnames(piresps[[tt]][[l]]) <- noquote(cod[[tt]][1, ])
			}
		}
		
		pirespm <- lapply(pirespm, function(x) lapply(x, function(y) round(y, prec)))
		piresps <- lapply(piresps, function(x) lapply(x, function(y) round(y, prec)))
		pirespm <- lapply(pirespm, function(x) lapply(x, function(y) t(y)))
		piresps <- lapply(piresps, function(x) lapply(x, function(y) t(y)))
		
		if( doVar2 ){
			for(tt in 1:TT2){
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
	if( doVar2 ){

		if( it2 > 0 & I > 0){		
					
			imp2_ <- vector("list", TT2)
			for( tt2 in 1:TT2 ){
				if( sum(R2[ ,tt2]) > 0 ){
					imp2_[[tt2]] <- matrix(NA, sum(R2[ ,tt2]), I)
					for( i in 1:I ){
						for( m2 in 1:ncat2[tt2] ){
							imp2_[[tt2]][which(imp2[[tt2]][ ,i] == cod2[[tt2]][2, m2]), i]<-as.numeric(noquote(cod2[[tt2]][1, m2]))
						}
					}
				}
			}

			imp2 <- imp2_
			rm(imp2_)
			names(imp2) <- nms2
			for( tt2 in 1:TT2 ){
				if(sum(R2[ ,tt2] > 0)){
					if( I > 0 ){
						colnames(imp2[[tt2]]) <- 1:I
						rownames(imp2[[tt2]]) <- which(R2[ ,tt2] == 1)
					}
				}
			}
			
			imp22 <- vector("list", J)
			for(j in 1:J){
				imp22[[j]]<-vector("list",TT)
				for( tt in 1:TT ){
					if( sum(R[[j]][ ,tt]) > 0 ){

						imp22[[j]][[tt]] <- matrix(NA, sum(R[[j]][ ,tt]), I)
						
						for( i in 1:I ){
							for( m in 1:ncat[tt] ){
								imp22[[j]][[tt]][which(imp[[j]][[tt]][ ,i] == cod[[tt]][2, m]), i] <- as.numeric(noquote(cod[[tt]][1, m]))
							}
						}
					}
				}
			}
			imp <- imp22
			rm(imp22)	

			for( j in 1:J ){
				names(imp)[j] <- noquote(paste("Group", j, sep=" "))
				for(tt in 1:TT){
					names(imp[[j]])[tt] <- nms[tt]
					if( sum(R[[j]][ ,tt]) > 0 ){
						colnames(imp[[j]][[tt]])[1:I] <- 1:I		
						rownames(imp[[j]][[tt]]) <- which(R[[j]][ ,tt] == 1)
					}
				}
			}		
		}
	# Create imputation sets (only level-1 )
	} else {

		if( it2 > 0 && I > 0){
			
			imp22 <- vector("list", J)
			
			for( j in 1:J ){
				imp22[[j]] <- vector("list", TT)
				
				for( tt in 1:TT ){
					if( sum(R[[j]][,tt]) > 0 ){

						imp22[[j]][[tt]] <- matrix(NA, sum(R[[j]][ ,tt]), I)

						
						for( i in 1:I ){
							for( m in 1:ncat[tt] ){
								imp22[[j]][[tt]][which(imp[[j]][[tt]][ ,i] == cod[[tt]][2, m]), i] <- as.numeric(noquote(cod[[tt]][1, m]))
							}
						}					
					}
				}
			}
			
			imp <- imp22
			rm(imp22)
			
			for( j in 1:J ){
				names(imp)[j] <- noquote(paste("Group", j, sep = " "))
				for( tt in 1:TT ){
					names(imp[[j]])[tt] <- nms[tt]
					if(sum(R[[j]][ ,tt]) > 0){
						colnames(imp[[j]][[tt]])[1:I] <- 1:I
						rownames(imp[[j]][[tt]]) <- which(R[[j]][ ,tt] == 1)
					}
				}
			}
		}
	}

	if( it2 == 0 || I == 0 ){
		imp = NULL
		imp2 = NULL
	}


	# Plots of log-likelihoods and classes distribution 

	if( plot.loglik & !count ){
		maxlik <- mc$maxlik
		pllk <- mc$pllk
		pllk <- maxlik - pllk
		plot(pllk, type = "l", lwd = 2, xlab = "t", ylab = "ll(max)-ll", main = "log-likelihood ratios", ylim = c(-.1, max(pllk, na.rm = TRUE)))
		countK <- NULL
		countL <- NULL
	} else if( !plot.loglik & count ){
		countK <- mc$countK
		countIT[countIT == 0] <- 1
		countK <- countK / countIT
		countL <- mc$countL
		countL <- countL/it2
		names(countL) <- paste("L=", 1:L, sep = "")
		rownames(countK) <- paste("W=", 1:L, sep = "")
		colnames(countK) <- paste("K=", 1:K , sep = "")
		if(L>1){
			par(mfrow = c(2,1))
			plot(countL, type = "h", lwd = "3", col = ifelse(countL == max(countL), "red", "black"), xlab = "Level 2 Class", ylab = "Freq.", main = "Distribution of L", xaxt = "n")
			axis(side = 1, at = 1:L, labels = paste("L=", 1:L, sep = ""))
			cols <- sample(1:657, L, rep = F)
			matplot(t(countK), type = "s", lwd = 3, ylim = c(0,1.5), lty = 1, xlab = "Level 1 Class", ylab = "Freq.", main = "Distribution of K|W", yaxt = "n", col = colors()[cols], xaxt = "n")
			if( L > 2 ){
				legend(1, 1.5, legend = paste("W = ", 1:L, sep = ""), fill = colors()[cols], ncol = floor(L / 3))
			} else {
				legend(1, 1.5, legend = paste("W = ", 1:L, sep = ""), fill = colors()[cols])
			}
			axis(side = 1, at = 1:K, labels = paste("K = ", 1:K, sep = ""))
			axis(side = 2, at = seq(0, 1, by = 0.1), labels = seq(0, 1, by = 0.1))
			par(mfrow = c(1, 1))
		} else {
			plot(1:K, countK, type = "h", lwd = "3", col = ifelse(countK == max(countK), "red", "black"), xlab = "Level 1 Class", ylab = "Freq.", main = "Distribution of K", xaxt = "n", ylim = c(0, 1))
			axis(side = 1, at = 1:K, labels = paste("K=", 1:K, sep = ""))
		}
	} else if( plot.loglik & count ){
		maxlik <- mc$maxlik
		pllk <- mc$pllk
		pllk <- maxlik-pllk
		countK <- mc$countK
		countIT[countIT == 0] <- 1
		countK <- countK / countIT
		rownames(countK) <- paste("W=", 1:L, sep = "")
		colnames(countK) <- paste("K=", 1:K, sep = "")
		countL <- mc$countL
		countL <- countL / it2
		names(countL) <- paste("L=", 1:L, sep = "")
		if( L > 1 ){
			par(mfrow = c(3, 1))
			plot(pllk, type = "l", lwd = 2, xlab = "t",ylab = "ll(max)-ll", main = "log-likelihood ratios", ylim = c(-.1, max(pllk, na.rm = TRUE)))
			plot(countL, type = "h", lwd = "3", col = ifelse(countL == max(countL), "red", "black"), xlab = "Level 2 Class", ylab = "Freq.", main = "Distribution of L", xaxt = "n")
			axis(side = 1, at = 1:L, labels = paste("L=", 1:L, sep = ""))
			cols <- sample(1:657, L, rep = F)
			matplot(t(countK), type = "s", lwd = 3, ylim = c(0, 1.5), lty = 1, xlab = "Level 1 Class", ylab = "Freq.", main = "Distribution of K|W", yaxt = "n", col = colors()[cols], xaxt = "n")
			if( L > 2 ){
				legend(1, 1.5, legend = paste("W = ", 1:L, sep = ""), fill = colors()[cols], ncol = floor(L / 3))
			} else {
				legend(1, 1.5, legend = paste("W = ", 1:L, sep = ""), fill = colors()[cols])
			}
			axis(side = 1, at = 1:K, labels = paste("K = ", 1:K, sep = ""))
			axis(side = 2, at = seq(0, 1, by = 0.1), labels = seq(0, 1, by = 0.1))
			par(mfrow = c(1, 1))
		}else{
			par(mfrow = c(2, 1))
			plot(pllk, type = "l", lwd = 2, xlab = "t", ylab = "ll(max)-ll", main = "log-likelihood ratios", ylim = c(-.1, max(pllk, na.rm = TRUE)))
			plot(1:K, countK, type = "h", lwd = "3", col = ifelse(countK == max(countK), "red", "black"), xlab = "Level 1 Class", ylab = "Freq.", main = "Distribution of K", xaxt = "n", ylim = c(0, 1))
			axis(side = 1, at = 1:K, labels = paste("K=", 1:K, sep = ""))
			par(mfrow = c(1, 1))
		}
	} else {
		countK <- NULL
		countL <- NULL
		countG <- NULL
	}

	

	return(list(
		implev1 = imp,
		implev2 = imp2,
		piL = piWm,
		piLses = piWs,
		piK = piKm,
		piKses = piKs,
		picondlev1 = pirespm,
		picondlev1ses = piresps,
		picondlev2 = piresp2m,
		picondlev2ses = piresp2s,
		DIC = DIC,
		freqL = countL,
		freqK = countK,
		time = ptmf
	))
}
