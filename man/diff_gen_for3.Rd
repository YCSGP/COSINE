\name{diff_gen_for3}
\alias{diff_gen_for3}
\title{Generate the F-statistics and ECF-statistics for the 
comparison of three datasets}
\description{
It works very similarly to "diff_gen" except that it 
performs the calculation across three rather than two datasets}
\usage{
diff_gen_for3(data1, data2, data3)
}
\arguments{
  \item{data1}{The first dataset of gene expression profiles}
  \item{data2}{The second dataset of gene expression profiles}
  \item{data3}{The third dataset of gene expression profiles}
}
\value{
  A list containing:
  \item{diff_expr}{a vector of the F-statistics for each gene}
  \item{diff_coex}{a square matrix storing the ECF-statistics 
                   for each gene pair}
}
\author{
Haisu Ma
}
\examples{

#Load the simulated datasets

data(simulated_data)

set1_data<-simulated_data[[1]]

set2_data<-simulated_data[[2]]

control_data<-simulated_data[[7]]

#Calculate the F-statistics and ECF-statistics 
#for the first five genes

diff_gen_for3_test <- diff_gen_for3(set1_data[,1:5],
set2_data[,1:5],control_data[,1:5])

}