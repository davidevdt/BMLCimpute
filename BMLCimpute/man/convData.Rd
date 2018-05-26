\name{convData}
\alias{convData}
\title{
	convData
}
\description{
This function takes a categorical dataset as input (categories can be denoted by numbers) and returns a list of objects that will be used by the 'multilevelLCMI' function to perform the imputations. 
}
\usage{
# Dataset called dat the group-ID is in the first column and the level-2 variables are in the 5th - 10th columns:
convData(dat, GID=1, UID=NULL, var2=5:10)
}
\arguments{
\item{dat}{
	Raw (categorical) dataset with missing data. The GID and UID variables, if passed to the function, must be in the first two columns of the dataset. 
}
 \item{GID}{
	Group Indicator (expressed as column number corresponding to the group ID in the dataset). It can be omitted in single-level datasets.  
}
  \item{UID}{
	Lower-level unit indicator (expressed as column number corresponding to the unit ID in the dataset). Optional. 
}
  \item{var2}{
	Higher-level variables (expressed as a vector of column numbers in the dataset corresponding to the variables measured at the higher levels). Optional. 
}


}
\details{
Convert a raw categorical dataset into one ready to be imputed with the multilevelLCMI function. In particular, the function will transform factor variables into numeric ones, in which numbers denote a different category. A coding list is returned along with the converted dataset. 
}


\value{
\item{convDat}{The converted dataset.}
\item{codLev1}{List containing the new (and original) scores which will be used for the imputations (Level-1 variables).}
\item{nCatLev1}{Vector containing the number of categories observed for each variable (Level-1 variables).}
\item{codLev2}{List containing the new (and original) scores which will be used for the imputations (Level-2 variables).}
\item{nCatLev2}{Vector containing the number of categories observed for each variable (Level-2 variables).}
\item{GroupIDs}{Matrix containing original and new Group ID's.}
\item{GID}{The column Group ID number (as entered in the input).}
\item{UID}{The column Unit ID number (as entered in the input). }
\item{var2}{The column numbers for level-2 variables (as entered in the input).}
\item{doVar2}{Boolean. Shall the BMLC model impute variables at level-2? (Result of !is.null(var2)).}
\item{namesLev1}{Vector of variable names (level-1 variables).}
\item{namesLev2}{Vector of variable names (level-2 variables).}
\item{GroupName}{Group ID variable name.}
\item{CaseName}{Unit ID variable name.}
\item{caseID}{Unit ID vector (re-permuted).}
\item{sort_}{Vector containing the original permutation of the dataset rows.}
}
\keyword{ preprocessing }
