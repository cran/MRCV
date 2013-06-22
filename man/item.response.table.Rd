\name{item.response.table}
\alias{item.response.table}
\title{
Create an Item Response Table or Data Frame
}
\description{
The \code{item.response.table} function is used to summarize data arising from two or three MRCVs.  For the two MRCV case, a cross-tabulation of the positive and negative responses for each (Wi, Yj) pair is presented as a table or data frame.  For the three MRCV case, a cross-tabulation of the positive and negative responses for each (Wi, Yj) pair is presented conditional on the response for each Zk.
}
\usage{
item.response.table(data, I, J, K = NULL, create.dataframe = FALSE)
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
  \item{create.dataframe}{
A logical value indicating whether the results should be presented as a data frame instead of a table.
}
}
\value{
For \code{create.dataframe = FALSE}, \code{item.response.table} uses the \code{\link{tabular}} function (package \pkg{tables}) to produce a 2Ix2J table of marginal counts (or multiple 2Ix2J tables in the case of three MRCVs).

For \code{create.dataframe = TRUE}, \code{item.response.table} returns the same information as above but presents it as a data frame.  For the two MRCV case, the data frame contains 2Ix2J rows and 5 columns generically named \code{W}, \code{Y}, \code{wi}, \code{yj}, and \code{count}.  For the three MRCV case, the data frame contains 2Ix2Jx2K rows and 7 columns generically named \code{W}, \code{Y}, \code{Z}, \code{wi}, \code{yj}, \code{zk}, and \code{count}.
}
\seealso{
The \code{\link{marginal.table}} function for creating a marginal table that summarizes only the positive responses for each pair.
}
\examples{
# Create an item response table for 2 MRCVs
farmer.irtable <- item.response.table(data = farmer2, I = 3, J = 4) 
farmer.irtable

# Create an item response data frame for 2 MRCVs
farmer.irdataframe <- item.response.table(data = farmer2, I = 3, J = 4, 
    create.dataframe = TRUE) 
farmer.irdataframe

# Create an item response table for 3 MRCVs
farmer.irtable.three <- item.response.table(data = farmer3, I = 3, J = 4, K = 5) 
farmer.irtable.three
}