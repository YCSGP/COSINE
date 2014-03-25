\name{diff_gen}
\alias{diff_gen}
\title{Calculate the F-statistics and ECF-statistics}
\description{The "diff_gen" function calculates the 
F-statistics which measures the differential expression 
of each gene and the ECF-statistcs which measures the 
differential correlation of each gene pair between two datasets
}

\usage{
diff_gen(data1, data2)
}

\arguments{
  \item{data1}{one of the two gene expression datasets}
  \item{data2}{a second gene expression dataset}
}

\value{
  A list containing:
  \item{diff_expr}{a vector of the F-statistics for each gene}
  \item{diff_coex}{a square matrix storing the 
                   ECF-statistics for each gene pair}
}

\author{Haisu Ma}

\examples{

#Load two of the simualted datasets

data(simulated_data)

set1_data <- simulated_data[[1]]
control_data <- simulated_data[[7]]

#Compute the F-statistics and ECF-statistics for the first 10 genes

diff_gen_test <- diff_gen(set1_data[,1:10],control_data[,1:10])

}