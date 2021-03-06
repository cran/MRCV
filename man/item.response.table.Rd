\name{item.response.table}
\alias{item.response.table}
\title{
Create an Item Response Table or Data Frame
}
\description{
The \code{item.response.table} function is used to summarize data arising from one, two, or three MRCVs.  For the one and two MRCV cases, a cross-tabulation of the positive and negative responses for each (Wi, Yj) pair is presented as a table or data frame (where Wi = W for the one MRCV case).  For the three MRCV case, a cross-tabulation of the positive and negative responses for each (Wi, Yj) pair is presented conditional on the response for each Zk.
}
\usage{
item.response.table(data, I, J, K = NULL, create.dataframe = FALSE)
}
\arguments{
  \item{data}{
A data frame containing the raw data where rows correspond to the individual item response vectors, and columns correspond to the items W1, \ldots, WI, Y1, \ldots, YJ, and Z1, \ldots, ZK (in this order).
}
  \item{I}{
The number of items corresponding to row variable W.  I = 1 for the one MRCV case.
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
For \code{create.dataframe = FALSE}, \code{item.response.table} uses the \code{\link{tabular}} function (package \pkg{tables}) to produce tables of marginal counts.

For \code{create.dataframe = TRUE}, \code{item.response.table} returns the same information as above but presents it as a data frame.  For the one MRCV case, the data frame contains rx2J rows and 4 columns generically named \code{W}, \code{Y}, \code{yj}, and \code{count}.  For the two MRCV case, the data frame contains 2Ix2J rows and 5 columns named \code{W}, \code{Y}, \code{wi}, \code{yj}, and \code{count}.  For the three MRCV case, the data frame contains 2Ix2Jx2K rows and 7 columns named \code{W}, \code{Y}, \code{Z}, \code{wi}, \code{yj}, \code{zk}, and \code{count}.
}
\seealso{
The \code{\link{marginal.table}} function for creating a marginal table that summarizes only the positive responses for each pair.
}
\examples{
# Create an item response table for 1 SRCV and 1 MRCV
farmer.irtable.one <- item.response.table(data = farmer1, I = 1, J = 5) 
farmer.irtable.one

# Create an item response data frame for 1 SRCV and 1 MRCV
farmer.irdataframe.one <- item.response.table(data = farmer1, I = 1, J = 5, 
    create.dataframe = TRUE) 
farmer.irdataframe.one

# Create an item response table for 2 MRCVs
farmer.irtable.two <- item.response.table(data = farmer2, I = 3, J = 4) 
farmer.irtable.two

# Create an item response table for 3 MRCVs
farmer.irtable.three <- item.response.table(data = farmer3, I = 3, J = 4, K = 5) 
farmer.irtable.three
}