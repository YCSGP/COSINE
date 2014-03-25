\name{cond.fyx}

\alias{cond.fyx}

\title{Compute the ECF-statistics measuring the differential 
       correlation of gene pairs}

\description{A function to calculate the expected conditional 
F-statistics as a measure of differential gene co-expression patterns.}

\usage{
cond.fyx(data.y, data.x, type)
}

\arguments{
  \item{data.y}{A vector containing the expression values 
                of one gene across two datasets}
  \item{data.x}{A vector containing the expression values 
                of another gene across two datasets}
  \item{type}{A vector indicating the type of each sample, 
        whose length is the sum of the sample sizes of data.y and data.x}
}

\value{The ECF-statistics of a specific gene}

\references{
http://bioinformatics.oxfordjournals.org/content/20/17/3146.long
}

\author{Yinglei Lai}

\examples{

#load two of the simulated datasets

data(simulated_data)
set1_data<-simulated_data[[1]]
control_data<-simulated_data[[7]]

num_sample <- dim(set1_data)[1]
num_gene <- dim(set1_data)[2]

type <- c(rep(0,num_sample),rep(1,num_sample))

#Compute the ECF-statistic for the gene pair of gene 1 and gene 2

i=1
j=2
data.x <- c(set1_data[,i],control_data[,i])
data.y <- c(set1_data[,j],control_data[,j])

ecf <- ( cond.fyx(data.y,data.x,type) + 
        cond.fyx(data.x,data.y,type) )/2

}