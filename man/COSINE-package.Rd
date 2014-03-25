\name{COSINE-package}
\alias{COSINE-package}
\alias{COSINE}
\docType{package}
\title{
COndition SpecIfic subNEtwork identification
}
\description{This is a package to identify the single globally optimal
subnetwork which differs the most between two or more datasets}
\details{
\tabular{ll}{
Package: \tab COSINE\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2014-03-25\cr
License: \tab GPL (version 2 or later)\cr
LazyLoad: \tab yes\cr
}

install.packages("COSINE")

}
\author{
Haisu Ma
Maintainer: Haisu Ma <haisu.ma.pku.2008@gmail.com>
}
\examples{

set.seed(666)

simulated_data <- DataSimu()

set1_data<-simulated_data[[1]][,1:100]

control_data<-simulated_data[[7]][,1:100]

diff_set1 <- diff_gen(set1_data,control_data)

scaled_diff_set1 <- score_scaling(diff_expr=diff_set1[[1]],
diff_coex=diff_set1[[2]])

klist<-c(25,30)

set1_quantile<-get_quantiles(diff_expr=scaled_diff_set1[[1]],
diff_coex=scaled_diff_set1[[2]],klist,pop_size=10)

lambda<-set1_quantile[[2]]

set1_GA<-GA_search(lambda,diff_expr=scaled_diff_set1[[1]],
diff_coex=scaled_diff_set1[[2]], num_iter=1, muCh=0.05, zToR=10)

set1_choose_lambda <- choose_lambda(diff_expr=scaled_diff_set1[[1]],
diff_coex=scaled_diff_set1[[2]],lambda,
subnet_size=set1_GA$Subnet_size,num_random_sampling=10,
best_score=set1_GA$Best_Scores)

}
