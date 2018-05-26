\name{simul_incomplete}
\docType{data}
\alias{simul_incomplete}
\title{Toy multilevel categorical dataset with missing entries}
\description{
The dataset contains 1000 level-1 artificial observations grouped into 200 level-2 groups. The dataset consists of: one identifier for the group, one identifier for the units, six level-1 binary variables (5 predictors and 1 response), five level-2 binary variables. The level-1 variables X1, X2, X5 and the level-2 variables Z1 and Z3 have missing observations, generated through a MAR mechanism. The dataset was generated following the simulation study in Vidotto, Vermunt, van Deun (in press) : Multilevel Latent Class models for the Multiple imputation of Nested Categorical Data.
}
\usage{data(simul)}
\format{A data frame with 13 columns and 1,000 rows from 200 groups.
\tabular{llll}{
[,1] \tab GroupID     \tab numeric  \tab Level-2 unit Identifier\cr
[,2] \tab UnitID   \tab numeric  \tab Area Identifier\cr
[,3] \tab X1,...,X5 \tab binary  \tab Level-1 predictors \cr
[,4] \tab Z1,...,Z5 \tab binary  \tab Level-2 predictors \cr
[,5] \tab Y   \tab binary  \tab Response Variable\cr
}
}
\references{
Vidotto, Vermunt, van Deun (in press). Multilevel Latent Class models for the Multiple imputation of Nested Categorical Data.
}
\keyword{data}