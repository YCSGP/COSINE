\name{get_quantiles_PPI}
\alias{get_quantiles_PPI}
\title{Get the five quantile values of lambda for analysis 
of gene expression and PPI network data}
\description{Perform random sampling a large number of times
to get the distribution of node score term and edge score term
and calculate the quantiles of lambda}
\usage{
get_quantiles_PPI(scaled_node_score, scaled_edge_score, PPI, klist, pop_size)
}
\arguments{
  \item{scaled_node_score}{The scaled F-statistics measuring the 
differential expression of each gene}
  \item{scaled_edge_score}{The scaled ECF-statistics measuring the
differential correlation of each gene pair}
  \item{PPI}{The two column matrix containing the protein interaction pairs}
  \item{klist}{A list of size of random networks to be sampled}
  \item{pop_size}{The number of random networks to be sampled}
}
\value{
  A list containing two components:
  \item{ratio }{The quantiles of the ratio of the edge 
                score term and node score term}
  \item{lambda }{The quantiles of lambda}
}
\author{Haisu Ma}
\examples{
data(scaled_node_score)
data(scaled_edge_score)
data(PPI)
quantiles<-get_quantiles_PPI(scaled_node_score,scaled_edge_score,
PPI,klist=seq(50,60,by=5),pop_size=10)
}