\name{set1_scaled_diff}
\alias{set1_scaled_diff}
\docType{data}
\title{
The standardized F-statistics and ECF-statistics for the comparison 
between simulated data1 and the control data
}
\description{

This data is a list containing two components, the first
component is a vector of length 500 which contains the 
F-statistics measuring differential expression of single 
genes, and the second component is a 500 by 500 matrix containing
the ECF-statistics measuring the differential correlation of gene pairs

}
\usage{data(set1_scaled_diff)}
\format{

  The format is:

List of 2

 $ : num [1:500] -0.416 -0.392 -0.326 -0.335 0.237 ...

 $ : num [1:500, 1:500] -0.986 -0.352 -0.828 -0.79 -0.208 ...
}
\examples{
data(set1_scaled_diff)
}
\keyword{datasets}
