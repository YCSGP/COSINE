\name{DataSimu}
\alias{DataSimu}
\title{Simulation of the six datasets and the case dataset}
\description{This function simulates six datasets with various 
differential expression and differential correlation patterns.}
\usage{
DataSimu()
}
\value{
  A list containing:
  \item{set1.data,set2.data,...,set6.data}{ case datasets 
        with different expression patterns to be 
        compared with the control datasets }
  \item{control.data}{ the control dataset }
}
\author{ Haisu Ma }
\examples{

set.seed(666)
simulated.data<-DataSimu()

}

\keyword{ datagen }
