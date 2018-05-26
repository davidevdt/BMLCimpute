\name{compData}
\alias{compData}
\title{
	compData
}
\description{
Plug the imputations obtained with multilevelLCMI into the original dataset, in order to obtain a complete dataset.  
}
\usage{
compData(convData,implev1,implev2,ind)
}
\arguments{
\item{convData}{
	Ouptut list produced by the 'convData' function.
}
 \item{implev1}{The set of imputations for the level-1 variables provided by the 'multilevelLCMI' function. It corresponds to the first element of the list returned by 'multilevelLCMI'.}
 \item{implev2}{The set of imputations for the level-2 variables (when present) provided by the 'multilevelLCMI' function. It corresponds to the second element of the list returned by 'multilevelLCMI'. Default to NULL.}
 \item{ind}{ The imputation index; an integer value that ranges from 1 to M, where M is the number of imputations (specified as I in the multilevelLCMI input.}


}
\details{
This function takes a 'convData' list, the imputations provided by 'multilevelLCMI' and the imputation index (ind is an integer which can ragne from 1 to M, where M is the number of imputation sets) and returns the completed dataset. 
}


\value{
\item{dataset}{The imputed dataset.}
}


\keyword{ postprocessing }
