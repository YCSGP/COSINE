\name{get_quantiles}
\alias{get_quantiles}
\title{Get the five quantiles of the weight parameter lambda}
\description{Use random sampling to get a large number 
of subnetworks and then calculate the distribution of 
the ratio between the edge-score term and node-score term}
\usage{
get_quantiles(diff_expr, diff_coex, klist, pop_size)
}
\arguments{
  \item{diff_expr}{The vector storing the F-statistics 
    measuring the differential expression of each gene}
  \item{diff_coex}{The matrix storing the ECF-statistics 
   measuring the differential correlation of each gene pair}
  \item{klist}{A vector of the sizes (number of genes) 
               of random subnetworks to be sampled}
  \item{pop_size}{The number of random subnetworks to be sampled}
}
\value{
    A list containing two components:
    \item{ratio}{The five quantiles of the log-ratios 
                 between the edge-score term and node-score term}
    \item{lambda}{The five quantiles of lambdas calculated 
                  based on the ratios}
}
\author{Haisu Ma}
\examples{

data(set1_scaled_diff)

klist<-c(20,25)

set1_quantile<-get_quantiles(diff_expr=set1_scaled_diff[[1]],
diff_coex=set1_scaled_diff[[2]],klist,pop_size=20)

}