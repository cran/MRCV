\name{MRCV-internal}
\alias{MRCV-internal}
\alias{check.margins}
\alias{check.min}
\alias{check.zero}
\alias{data.format}
\alias{est.jack}
\alias{ipf.genloglin}
\title{
Internal MRCV Package Functions
}
\description{
Hidden internal MRCV package functions not meant to be called by the user.
}
\usage{
check.margins(data, I, J, K, nvars, model.vars, item.names)

check.min(value)

check.zero(value, add.constant = .5)

data.format(data, I, J, K, nvars, add.constant = .5, model.vars = NULL, 
    predict.func = FALSE)
    
est.jack(mu.hat, i, I, J, K)

ipf.genloglin(data, I, J, K, nvars, p, p.theta.2, p.theta.3, x.theta.2, 
    x.theta.3)
}
\arguments{
  \item{data}{
For \code{check.margins}, the marginal counts for a bootstrap resample.

For \code{data.format} and observed data, a data frame containing the raw MRCV data. For bootstrap resamples, a data frame containing the joint table cell counts.

For \code{ipf.genloglin}, a data frame containing the raw MRCV data.
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
  \item{nvars}{
The number of MRCVs.
}
  \item{model.vars}{
For \code{check.margins} and the two MRCV case, a data frame containing 2Ix2J rows and 4 columns generically named \code{W}, \code{Y}, \code{wi}, and \code{yj}.  For the three MRCV case, a data frame containing 2Ix2Jx2K rows and 6 columns generically named \code{W}, \code{Y}, \code{Z}, \code{wi}, \code{yj}, and \code{zk}.

For \code{data.format} and observed data: \code{model.vars = NULL}. For bootstrap resamples and the two MRCV case, a data frame with 2^(I+J) rows representing all possible combinations of (Wi,Yj). For the three MRCV case, a data frame with 2^(I+J+K) rows reprsenting all possible combinations of (Wi, Yj, Zk).
}
  \item{item.names}{
The names of the MRCV items as labeled and ordered in the data frame containing the raw MRCV data.
}
  \item{value}{
A numeric vector of length = 1.
}
  \item{add.constant}{
A positive constant to be added to all zero marginal cell counts.
}
  \item{predict.func}{
A logical value indicating whether \code{\link{predict.genloglin}} is being called.
}
  \item{mu.hat}{
A data frame containing the estimated counts from a model based on \emph{n}-1 subjects.
}
  \item{i}{
The current index value for the \code{\link{predict.genloglin}} \code{for()} loop.
}
  \item{p}{
A data frame with 2^(I+J) (or 2^(I+J+K)) rows and column variables \code{W1},\ldots,\code{WI}, \code{Y1},\ldots,\code{YJ}, \code{Z1},\ldots,\code{ZK}, and \code{p} (in this order).  The third set of items is only necessary when analyzing the relationship among 3 MRCVs. \code{p} is the (m-1)th iterative estimate of the multinomial probability generated under the null hypothesis.
}
  \item{p.theta.2}{
For the two MRCV case. See Gange (1995, p. 136) for details.  Note that \code{ipf.genloglin} uses cell counts instead of probabilities.
}
  \item{p.theta.3}{
For the three MRCV case. See Gange (1995, p. 136) for details.  Note that \code{ipf.genloglin} uses cell counts instead of probabilities.
}
  \item{x.theta.2}{
For the two MRCV case. See Gange (1995, p. 136) for details.  Note that \code{ipf.genloglin} uses cell counts instead of probabilities.
}
  \item{x.theta.3}{
For the three MRCV case. See Gange (1995, p. 136) for details.  Note that \code{ipf.genloglin} uses cell counts instead of probabilities.
}
}
\details{
The \code{check.margins} function is used to determine whether a bootstrap resample can be validly used by the \code{\link{genloglin}} function.  A resample is valid only if all binary items have at least one positive and one negative response.

The \code{check.min} function is called by \code{\link{MI.test}} to compute the minimum fora Bonferroni adjustment that results in a p-value equal to the minimum of 1 and the adjusted p-value.

The \code{check.zero} function is used to add a small constant to zero cell counts.  This is necessary to prevent zero or undefined odds ratio estimates, and to ensure correct convergence of the iterative proportional fitting algorithm described in Gange (1995).

The \code{data.format} function is used to convert raw or bootstrap data for two or three MRCVs into a data frame containing the marginal counts.

The \code{est.jack} function calculates model-predicted odds ratios.

The \code{ipf.genloglin} function carries out the mth iteration of the iterative proportional fitting algorithm described by Gange (1995, p. 136).  This function is called iteratively by \code{\link{genloglin}} until the convergence criterion is reached.
}
\value{
The \code{check.margins} function returns a logical value indicating whether the resample is valid.

The \code{check.min} function returns a numeric object equal to 1 if \code{value} > 1; otherwise, the original value is returned.

The \code{check.zero} function returns a numeric object equal to \code{add.constant} if \code{value} = 0; otherwise, the original value is returned.

The \code{data.format} function returns a summary data frame containing the marginal counts.

The \code{est.jack} function returns the model-predicted odds ratio using data from \emph{n}-1 subjects.

The \code{ipf.genloglin} function returns an updated version of \code{p} containing the mth iterative estimate of the multinomial probabilities.
}
\references{
Gange, S. (1995)  Generating multivariate categorical variables using the iterative proportional fitting algorithm.  \emph{The American Statistician}, \bold{49}, 134--138.
}
\keyword{internal}