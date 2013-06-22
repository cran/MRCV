\name{marginal.table}
\alias{marginal.table}
\title{
Create a Marginal Table
}
\description{
The \code{marginal.table} function is used to summarize data arising from two or three MRCVs.  For the two MRCV case, a cross-tabulation of the positive responses for each (Wi, Yj) pair is presented as a table.  For the three MRCV case, a cross-tabulation of the positive responses for each (Wi, Yj) pair is presented conditional on the positive responses for each Zk.
}
\usage{
marginal.table(data, I, J, K = NULL)
}
\arguments{
  \item{data}{
A data frame containing the raw data where rows correspond to the individual item response vectors, and columns correspond to the binary items W1, \ldots, WI, Y1, \ldots, YJ, and Z1, \ldots, ZK (in this order).  The third set of items is only necessary when analyzing the relationship among three MRCVs.
}
  \item{I}{
The number of items corresponding to row variable W.
}
  \item{J}{
The number of items corresponding to column variable Y.
}
  \item{K}{
The number of items corresponding to strata variable Z.
}
}
\value{
The \code{marginal.table} function uses the \code{\link{tabular}} function (package \pkg{tables}) to produce a table containing the IxJ positive counts (or K tables containing the IxJ positive counts in the case of three MRCVs).
}
\seealso{
The \code{\link{item.response.table}} function for creating an item response table or data frame that summarizes both the positive and negative responses for each (Wi, Yj) pair (conditional on the response for each Zk in the case of three MRCVs).
}
\examples{
# Create a marginal table for 2 MRCVs
farmer.mtable <- marginal.table(data = farmer2, I = 3, J = 4) 
farmer.mtable

# Create a marginal table for 3 MRCVs
farmer.mtable.three <- marginal.table(data = farmer3, I = 3, J = 4, K = 5) 
farmer.mtable.three
}