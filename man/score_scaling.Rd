\name{score_scaling}
\alias{score_scaling}
\title{To get the normalzied F-statistics and ECF-statistics}
\description{Perform standardization of the node score and edge score}
\usage{
score_scaling(diff_expr, diff_coex)
}
\arguments{
  \item{diff_expr}{The vector storing the F-statistics measuring the differential
                   expression of each gene}
  \item{diff_coex}{The matrix storing the ECF-statistics measuring the differential
                   correlation of each gene pair}
}
\value{
  A list containing two components:
  \item{diff_expr}{A vector of the standardized F-statistics}
  \item{diff_coex}{A matrix of the standardized ECF-statistics}
}

\author{Haisu Ma}
\examples{

data(set1_unscaled_diff)

scaled_diff_set1 <- score_scaling(diff_expr=set1_unscaled_diff[[1]],
diff_coex=set1_unscaled_diff[[2]])

}