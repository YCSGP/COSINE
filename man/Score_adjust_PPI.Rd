\name{Score_adjust_PPI}
\alias{Score_adjust_PPI}
\title{To adjust the score of the selected PPI 
sub-network using random sampling}
\description{Randomly sample a large number of node-sets 
and edge-sets from the background PPI network
to derive the null distribution of the scores for 
subnetworks with certain size, and to compute the
adjusted scores for the selected sub-network}
\usage{
Score_adjust_PPI(scaled_node_score, scaled_edge_score, 
PPI, lam, subnet, num_random_sampling, best_score)
}
\arguments{
  \item{scaled_node_score}{The scaled F-statistics of each node(gene) in the network}
  \item{scaled_edge_score}{The scaled ECF-statistics of each edge in the network}
  \item{PPI}{A matrix with two columns containing the interacting gene pairs}
  \item{lam}{The weight parameter lambda used for the selection of this sub-network}
  \item{subnet}{A vector of the index of selected genes}
  \item{num_random_sampling}{Number of random subnetworks to be sampled}
  \item{best_score}{The original score of selected sub-network}
}
\value{The adjusted score of the selected sub-network}
\author{Haisu Ma}
\examples{
data(scaled_node_score)
data(scaled_edge_score)
data(PPI)
data(set1_GA)
adj_Score<-Score_adjust_PPI(scaled_node_score,scaled_edge_score,
PPI,lam=0.1,subnet=set1_GA$Subnet[[1]],num_random_sampling=2,
best_score=set1_GA$Best_Scores[1])

}