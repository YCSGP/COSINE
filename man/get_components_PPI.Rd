\name{get_components_PPI}
\alias{get_components_PPI}
\title{Get all the components (connected clusters) of the sub-network}
\description{Map the edges in the selected sub-network 
to the background PPI network and get all the clusters
with size larger than the minimumset by the user}
\usage{
get_components_PPI(gene_names, vector, PPI, minsize)
}
\arguments{
  \item{gene_names}{The gene names of all the nodes}
  \item{vector}{A binary vector indicating whether each node is selected or not}
  \item{PPI}{A two column matrix including the protein interaction data}
  \item{minsize}{The minimal size of clusters}
}
\value{A list with each element corresponding to one cluster in the selected sub-network}
\author{Haisu Ma}
\examples{

data(scaled_node_score)
data(scaled_edge_score)
data(PPI)
gene_names<-names(scaled_node_score)
vector<-rep(0,length(scaled_node_score))
vector[sample(1:length(scaled_node_score),length(scaled_node_score)/3)]<-1
components<-get_components_PPI(gene_names,vector,PPI,minsize=3)

}