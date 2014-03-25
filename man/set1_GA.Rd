\name{set1_GA}
\alias{set1_GA}
\docType{data}
\title{
Result of genetic algorithm search for simulated data set1
}
\description{
This data set contains the result of the subnetwork extraction using 
genetic algorithm applied to the analysis of the differential expression
pattern between simulated dataset1 and the control dataset
}
\usage{data(set1_GA)}
\format{

  The format is:

List of 4

 $ Subnet_size: num [1:5] 23 58 61 75 99

 $ Best_Scores: num [1:5] 16.7 24.7 26.6 29.3 77.7

 $ Subnet     :List of 5

  ..$ : int [1:23] 36 40 41 112 121 148 163 184 185 206 ...
  ..$ : int [1:58] 25 36 38 39 40 41 61 71 78 79 ...
  ..$ : int [1:61] 25 36 38 39 40 41 48 61 78 79 ...
  ..$ : int [1:75] 8 25 36 38 39 40 41 61 71 78 ...
  ..$ : int [1:99] 3 5 8 11 25 29 36 38 39 40 ...

 $ GA_obj     :List of 5

  ..$ :List of 11
  .. ..$ type          : chr "binary chromosome"
  .. ..$ size          : int 500
  .. ..$ popSize       : num 200
  .. ..$ iters         : num 1000
  .. ..$ suggestions   : NULL
  .. ..$ population    : num [1:200, 1:500] 0 0 0 0 0 0 0 0 0 0 ...
  .. ..$ elitism       : num 40
  .. ..$ mutationChance: num 0.05
  .. ..$ evaluations   : num [1:200] -16.7 -16.7 -16.7 -16.7 -16.7 ...
  .. ..$ best          : num [1:1000] -2.85 -2.85 -3.68 -4.49 -4.59 ...
  .. ..$ mean          : num [1:1000] 0.0192 -0.743 -1.4033 -2.1104 -2.7789 ...
  .. ..- attr(*, "class")= chr "rbga"

  ..$ :List of 11
  .. ..$ type          : chr "binary chromosome"
  .. ..$ size          : int 500
  .. ..$ popSize       : num 200
  .. ..$ iters         : num 1000
  .. ..$ suggestions   : NULL
  .. ..$ population    : num [1:200, 1:500] 0 0 0 0 0 0 0 0 0 0 ...
  .. ..$ elitism       : num 40
  .. ..$ mutationChance: num 0.05
  .. ..$ evaluations   : num [1:200] -24.7 -24.7 -24.7 -24.7 -24.7 ...
  .. ..$ best          : num [1:1000] -4.27 -5.88 -6.64 -6.73 -7.95 ...
  .. ..$ mean          : num [1:1000] 0.00387 -1.27523 -2.59219 -3.73991 -4.7593 ...
  .. ..- attr(*, "class")= chr "rbga"

  ..$ :List of 11
  .. ..$ type          : chr "binary chromosome"
  .. ..$ size          : int 500
  .. ..$ popSize       : num 200
  .. ..$ iters         : num 1000
  .. ..$ suggestions   : NULL
  .. ..$ population    : num [1:200, 1:500] 0 0 0 0 0 0 0 0 0 0 ...
  .. ..$ elitism       : num 40
  .. ..$ mutationChance: num 0.05
  .. ..$ evaluations   : num [1:200] -26.6 -26.6 -26.6 -26.6 -26.6 ...
  .. ..$ best          : num [1:1000] -5.13 -6.63 -6.84 -7.4 -8.47 ...
  .. ..$ mean          : num [1:1000] 0.0412 -1.3408 -2.6099 -3.9691 -5.1945 ...
  .. ..- attr(*, "class")= chr "rbga"

  ..$ :List of 11
  .. ..$ type          : chr "binary chromosome"
  .. ..$ size          : int 500
  .. ..$ popSize       : num 200
  .. ..$ iters         : num 1000
  .. ..$ suggestions   : NULL
  .. ..$ population    : num [1:200, 1:500] 0 0 0 0 0 0 0 0 0 0 ...
  .. ..$ elitism       : num 40
  .. ..$ mutationChance: num 0.05
  .. ..$ evaluations   : num [1:200] -29.3 -29.2 -29.2 -29.2 -29.2 ...
  .. ..$ best          : num [1:1000] -4.43 -4.71 -6.21 -6.92 -8.31 ...
  .. ..$ mean          : num [1:1000] -0.126 -1.546 -2.71 -3.8 -4.788 ...
  .. ..- attr(*, "class")= chr "rbga"

  ..$ :List of 11
  .. ..$ type          : chr "binary chromosome"
  .. ..$ size          : int 500
  .. ..$ popSize       : num 200
  .. ..$ iters         : num 1000
  .. ..$ suggestions   : NULL
  .. ..$ population    : num [1:200, 1:500] 0 0 0 0 0 0 0 0 0 0 ...
  .. ..$ elitism       : num 40
  .. ..$ mutationChance: num 0.05
  .. ..$ evaluations   : num [1:200] -77.7 -77.6 -77.3 -77.3 -77.3 ...
  .. ..$ best          : num [1:1000] -11.7 -15.2 -15.2 -20.6 -21.4 ...
  .. ..$ mean          : num [1:1000] -0.0381 -2.8258 -5.5397 -8.3003 -10.7343 ...
  .. ..- attr(*, "class")= chr "rbga"

}

\details{
   This dataset is a list containing the following components:
   Subnet_size: A vector of length 5 showing the size of the 
                 selected subnetwork using five different lambdas.
   Best_Scores: The scores of the selected subnetworks 
                      corresponding to five lambdas.
   Subnet: The selected subnetworks (gene indices) for five lambdas.
   GA_obj: The objected returned by the function "rbga.bin", 
           which stores the results of the genetic algorithm.
}

\examples{
data(set1_GA)
}

\keyword{datasets}