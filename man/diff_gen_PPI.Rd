\name{diff_gen_PPI}
\alias{diff_gen_PPI}
\title{Generate the scaled node score and scaled
edge score for nodes and edges in the background network}
\description{Compute the F-statistic and ECF-statistic 
and then standardize them}
\usage{
diff_gen_PPI(data1, data2, PPI)
}
\arguments{
  \item{data1}{The first gene expression dataset(with rows corresponding to samples)}
  \item{data2}{The second gene expression dataset}
  \item{PPI}{A matrix with two columns containing the protein interaction pairs}
}
\value{
  A list containing:
  \item{scaled_node_score }{The standardized F-statistic measuring 
                            the differential expression of each gene}
  \item{scaled_edge_score }{The standardized ECF-statistic measuring
                            the differential correlation of each gene pair}
}
\author{Haisu Ma}
\examples{

data(simulated_data)
data(PPI)
data1 <- simulated_data[[1]]
data2 <- simulated_data[[7]]
colnames(data1)<-colnames(data2)<-as.character(1:500)
test <- diff_gen_PPI(data1[,1:20],data2[,1:20],PPI)
}