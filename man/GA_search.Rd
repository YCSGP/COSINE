\name{GA_search}
\alias{GA_search}
\title{
Use genetic algorithm to search for the globally optimal subnetwork
}
\description{
This function performs the stochastic search using genetic 
algorithm to find the globally optimal subnetwork which gives 
rise to the highest score defined by a scoring function, which
measures the extent of the differential expression of the 
subnetwork across several datasets.
}
\usage{
GA_search(lambda, diff_expr, diff_coex, num_iter = 1000, 
muCh = 0.05, zToR = 10)
}
\arguments{
  \item{lambda}{A vector containing the five quantiles of the 
               weight parameter lambda}
  \item{diff_expr}{A vector storing the F-statistics measuring 
                  the differential expression of each gene, 
                  which length equals the number of genes N}
  \item{diff_coex}{An N by N matrix with entry (i,j) corresponding 
              to the ECF-statistics of gene pair (i,j), which 
              measures the differential correlation between genes i and j}
  \item{num_iter}{The number of iterations to be performed 
                  by the genetic algorithm}
  \item{muCh}{the mutation chance used by genetic algorithm}
  \item{zToR}{zero to one ratio}
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

\references{
http://cran.r-project.org/web/packages/genalg/index.html
}

\author{
Haisu Ma
}

\examples{

# Load the scaled F-statistics and ECF-statistics 
# for the simulated datasets

data(set1_scaled_diff)

# Get the quantiles of lambda

klist<-c(25,30)
set1_quantile<-get_quantiles(diff_expr=set1_scaled_diff[[1]], 
diff_coex=set1_scaled_diff[[2]],klist,pop_size=10)
lambda<-set1_quantile[[2]]

#Perform genetic algorithm to search-just show the first iteration here

set1_GA<-GA_search(lambda[1:2],diff_expr=set1_scaled_diff[[1]],
diff_coex=set1_scaled_diff[[2]], num_iter=1, muCh=0.05, zToR=50)

}




