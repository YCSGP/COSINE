\name{f.test}
\alias{f.test}
\title{To get the F-statistics for each gene}
\description{
Calculate the F-statistics measuring the differential 
expression of each gene
}
\usage{
f.test(data, type)
}
\arguments{
  \item{data}{A vector containing the expression values of 
              a gene across two datasets}
  \item{type}{A vector indicating the type of each sample 
             (coming  from dataset 1 or dataset 2) }
}
\value{The F-statistics of a specific gene}
\author{Haisu Ma}
\examples{

#Load the simulated datasets

data(simulated_data)
data1 <- simulated_data[[1]]
data2 <- simulated_data[[7]]

#Calculate the F-statistics for genes 1~10

num_sample <- dim(data1)[1]
diff_expr <- rep(0,10)
type <- c(rep(0,num_sample),rep(1,num_sample))

for(i in 1:10){
       data <- c(data1[,i],data2[,i])
       diff_expr[i] <- f.test(data,type)
    }

}