\name{GA_search_PPI}
\alias{GA_search_PPI}
\title{Run genetic algorithm to search for the PPI sub-network}
\description{
This function performs the stochastic search using genetic 
algorithm to find the globally optimal subnetwork which gives 
rise to the highest score defined by a scoring function, which
measures the extent of the differential expression of the PPI
subnetwork across several datasets.
}
\usage{
GA_search_PPI(lambda, scaled_node_score, scaled_edge_score, PPI, 
num_iter = 1000, muCh = 0.05, zToR = 10, minsize = 10)}
\arguments{
  \item{lambda}{One of the five quantiles of the weight parameter lambda}
  \item{scaled_node_score}{A vector storing the F-statistics measuring 
                  the differential expression of each gene, 
                  which length equals the number of genes N}
  \item{scaled_edge_score}{A vector storing the ECF-statistics measuring 
                  the differential correlation of each gene pair}
  \item{PPI}{A two-column matrix containing the protein interaction pairs}
  \item{num_iter}{The number of iterations to be performed 
                  by the genetic algorithm}
  \item{muCh}{the mutation chance used by genetic algorithm}
  \item{zToR}{zero to one ratio}
  \item{minsize}{The minimal size of selected sub-network}
}
\value{
    A list containing the following components:
    \item{Subnet_size}{A vector containing the size of the 
                       subnetwork identified using each lambda}
    \item{Best_Scores}{A vector containing the best scores 
                       of the subnetworks}
    \item{Subnet}{A list containing the extracted subnetworks 
                 (a list of genes) for each of the five lambda values}
    \item{GA_obj}{A list of the returned objects of the 
                  genetic algorithm function}
}
\author{Haisu Ma}
\examples{

data(scaled_node_score)
data(scaled_edge_score)
data(PPI)

GA_result<-GA_search_PPI(lambda=0.5,scaled_node_score,scaled_edge_score,PPI,
num_iter=1, muCh=0.05, zToR=10, minsize=50)

}