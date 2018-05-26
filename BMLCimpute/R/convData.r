# Package BMLCimpute
# Function convData 
# 
# This function takes a categorical dataset as input (categories can be denoted by numbers) and returns a list of objects that will be used by the 'multilevelLCMI' function to perform the imputations. 
#
# Inputs : 
# 
# @param dat The categorical dataset
# @param GID Integer. Column number denoting the Group ID (optional for single-level datasets; mandatory for multilevel datasets). If present, it should be in either column 1 (GID = 1) or column 2 (GID = 2)
# @param UID Integer. Column number denoting Unit ID (optional). If present, it should be in either column 1 (UID = 1) or column 2 (UID = 2)
# @param var2 Integer. Vector of column number(s) denoting the variables observed at the higher-units level. Mandatory only if such variables are present in the dataset
#
# Output : 
# 
# @return convData The converted dataset
# @return codLev1 List containing the new (and original) scores which will be used for the imputations (Level-1 variables)
# @return nCatLev1 Vector containing the number of categories observed for each variable (Level-1 variables) 
# @return codLev2 List containing the new (and original) scores which will be used for the imputations (Level-2 variables)
# @return nCatLev2 Vector containing the number of categories observed for each variable (Level-2 variables)
# @return GroupIDs Matrix containing original and new Group ID's 
# @return GID The column Group ID number (as entered in the input)
# @return UID The column Unit ID number (as entered in the input)
# @return var2 The column numbers for level-2 variables (as entered in the input)
# @return doVar2 Boolean. Shall the BMLC model impute variables at level-2? (Result of !is.null(var2))
# @return namesLev1 Vector of variable names (level-1 variables)
# @return namesLev2 Vector of variable names (level-2 variables)
# @return groupName Group ID variable name
# @return caseName Unit ID variable name
# @return caseID Unit ID vector (re-permuted) 
# @return sort_ Vector containing the original permutation of the dataset rows


convData <- function( dat, GID = NULL, UID = NULL, var2 = NULL ){

	# If the Group ID is not in the first 2 columns, stop the program. 
	if(! is.null(GID) ){
		if( GID > 2 ){
			stop("Group ID must be in either the first two columns.", call. = FALSE)
		}
	}

	# If the Unit ID is not in the first 2 columns, stop the program. 
	if(! is.null(UID) ){
		if( UID > 2 ){
			stop("Case ID must be in either the first two columns.", call. = FALSE)
		}
	}

	# Are level-2 variables present in the dataset? 
	doVar2 <- TRUE
	if( is.null(var2) ){
		doVar2 <- FALSE
	}

	# For single-level datasets : create a fictitious Group-ID (with only 1 group)
	if(! is.null(GID) ){
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
	if( is.null(UID) ){
		Y <- apply(Y, 2, function(x) as.numeric(as.factor(x)))
	} else {
		Y <- apply(Y[ ,-UID], 2, function(x) as.numeric(as.factor(x)))
	}
	n <- nrow(Y)
	nms_tmp <- colnames(dat)

	# Column Name for Unit ID
	CN <- NULL
	CID <- NULL
	if(! is.null(UID) ){
		CN <- colnames(dat[UID])
		CID <- dat[ ,UID]
	}

	# Separate level-2 and level-1 variable names
	nms2 <- NULL
	if( doVar2 ){
		if( !is.null(UID) ){
			nms2 <- nms_tmp[var2]
			nms <- nms_tmp[-c(GID, UID, var2)]
		}else{
			nms2 <- nms_tmp[var2]
			nms <- nms_tmp[-c(GID, var2)]
		}
	} else {
		if( !is.null(UID) ){
			nms <- nms_tmp[-c(GID, UID)]
		}else{
			nms <- nms_tmp[-GID]
		}
	}
	rm(nms_tmp)

	# New codes for Group ID
	dtGID <- data.frame(dat[,GID], Y[,1])
	dtGID <- unique(dtGID)
	dtGID <- t(dtGID[order(dtGID[,2]),])
	rownames(dtGID) <- cbind("Original GroupID", "New GroupID")
	colnames(dtGID) <- NULL

	# Sort rows by Group ID 
	sortCase <- 1:n
	sortCase <- sortCase[order(Y[,1])]
	Y <- Y[order(Y[ ,1]), ]
	CID <- CID[order(Y[ ,1])]

	# Create level-2 variables dataset (and list with new and original codes)
	ncat2<-rep(0,1)
	cod2<-0
	ncat2<-0

	if( doVar2 ){
		ncat2 <- rep(0, length(var2))
		cod2 <- vector("list", length(var2))
		for( tt in 1:length(var2) ){
			ncat2[tt] <- length(table(Y[ ,var2[tt]]))
			cod2[[tt]] <- matrix(0,2,ncat2[tt])
			names(cod2)[tt] <- nms2[tt]
			cod2[[tt]][1,] <- names(table(dat[ ,nms2[tt]]))
			cod2[[tt]][2,] <- names(table(Y[,var2[tt]]))
			rownames(cod2[[tt]]) <- c("Original Values", "Coded Values")
		}
	}

	# Original and new codes for level-1 variables  
	var1 <- which(colnames(Y) %in% nms)
	ncat <- rep(0, length(var1))
	cod <- vector("list", length(var1))

	for( tt in 1:length(var1) ){
		ncat[tt] <- length(table(Y[ ,var1[tt]]))
		cod[[tt]] <- matrix(0, 2, ncat[tt])
		names(cod)[tt] <- nms[tt]
		cod[[tt]][1,] <- names(table(dat[ ,nms[tt]]))
		cod[[tt]][2,] <- names(table(Y[ ,var1[tt]]))
		rownames(cod[[tt]]) <- c("Original Values","Coded Values")
	}



	return( list(
		convDat = Y,
		codLev1 = cod,
		nCatLev1 = ncat,
		codLev2 = cod2,
		nCatLev2 = ncat2,
		GroupIDs = dtGID,
		GID = GID,
		UID = UID,
		var2 = var2,
		doVar2 = doVar2,
		namesLev1 = nms,
		namesLev2 = nms2,
		groupName = GN,
		caseName = CN,
		caseID = CID,
		sort_ = sortCase
	))

}
