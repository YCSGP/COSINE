\name{random_network_sampling_PPI}
\alias{random_network_sampling_PPI}
\title{To sample random sub-network from the PPI data}
\description{Randomly sample a set of nodes from the gene pool,
and check the number of edges contained, if there are edges among
the nodes, return the random sub-network}
\usage{
random_network_sampling_PPI(size, PPI, all_genes)
}
\arguments{
  \item{size}{number of nodes to be sampled}
  \item{PPI}{The PPI network data}
  \item{all_genes}{The names of all genes}
}
\value{
  \item{current }{The names of selected genes}
}
\author{Haisu Ma}
\examples{
data(PPI)
all_genes<-union(PPI[,1],PPI[,2])
ran_net<-random_network_sampling_PPI(size=30,PPI,all_genes)
}