\name{SPMI.test}
\alias{SPMI.test}
\alias{SPMI.stat}
\title{
Test for Marginal Independence Between Two MRCVs
}
\description{
The \code{SPMI.test} function offers three approaches for testing for simultaneous pairwise marginal independence (SPMI) between two MRCVs.
}
\usage{
SPMI.test(data, I, J, type = "all", B = 1999, B.max = B, summary.data = FALSE, 
    add.constant = 0.5, plot.hist = FALSE, print.status = TRUE)

SPMI.stat(data, I, J, summary.data = FALSE, add.constant = 0.5)
}
\arguments{
  \item{data}{
For \code{summary.data = FALSE}: a data frame containing the raw data where rows correspond to the individual item response vectors, and columns correspond to the binary items, W1, \ldots, WI and Y1, \ldots, YJ (in this order).

For \code{summary.data = TRUE}: a data frame containing 2Ix2J rows and 5 columns generically named \code{W}, \code{Y}, \code{wi}, \code{yj}, and \code{count}.
}
  \item{I}{
The number of items corresponding to row variable W.
}
  \item{J}{
The number of items corresponding to column variable Y.
}
  \item{type}{
A character string specifying one of the following approaches for testing for SPMI: \code{"boot"} specifies a nonparametric bootstrap procedure; \code{"rs2"} specifies a Rao-Scott second-order adjustment; \code{"bon"} specifies a Bonferroni adjustment; \code{"all"} specifies all three approaches.
}
  \item{B}{
The desired number of bootstrap resamples.
}
  \item{B.max}{
The maximum number of bootstrap resamples.  Resamples for which at least one item has all positive or negative responses are thrown out; \code{SPMI.test} uses the first \code{B} valid resamples or all valid resamples if that number is less than \code{B}.
}
  \item{summary.data}{
A logical value indicating whether \code{data} is a summary file containing the item response data instead of the raw data.  Only \code{type = "bon"} is available for \code{summary.data = TRUE}.
}
  \item{add.constant}{
A positive constant to be added to all zero marginal cell counts.
}
  \item{plot.hist}{
A logical value indicating whether plots of the emprical bootstrap sampling distributions should be provided.
}
  \item{print.status}{
A logical value indicating whether bootstrap progress updates should be provided.
}
}
\details{
The \code{SPMI.test} function calls \code{SPMI.stat}, which in turn calculates a modified Pearson statistic defined in Bilder and Loughin (2004).  The modified statistic is the sum of the IxJ Pearson statistics used to test for independence of each (Wi,Yj) pair.  The asymptotic distribution of the modified statistic is a linear combination of independent chi-square(1) random variables, so traditional methods for analyzing the association between categorical variables W and Y are inappropriate.  For this reason, the \code{SPMI.stat} function is not normally called directly.

The \code{SPMI.test} function offers three testing methods.  The nonparametric bootstrap method empirically estimates the sampling distribution of the modified statistic.  Resampling under the null hypothesis is done by independently sampling the W and Y vectors with replacement for s = 1, \ldots, \emph{n} where \emph{n} is the sample size.  In addition to bootstrapping the Pearson statistic, bootstrap p-value combination methods are available in which the minimum or the product of the IxJ p-values (each using a chi-square(1) approximation for the corresponding Pearson statistic) is calculated for each resample.

The Rao-Scott approach applies a second-order adjustment to the modified Pearson statistic and its sampling distribution.  Formulas are provided in Appendix A of Bilder and Loughin (2004).

Using the Bonferroni adjustment, SPMI is rejected if at least one of the IxJ p-values is less than alpha/(IxJ).  The \code{SPMI.test} function provides the adjusted p-values so SPMI is rejected if at least one of the p-values is less than alpha.  An overall Bonferroni adjusted p-value is given as the minimum of the IxJ p-values multiplied by IxJ, where SPMI is rejected if this value is less than alpha.  Note that the Bonferroni adjustment tends to produce an overly conservative test when I and J are large.
}
\value{
--- For \code{SPMI.test}, a list containing at least \code{general}, a list containing the following objects:
\itemize{
\item{\code{data}:} {The original data frame supplied to the \code{data} argument.}
\item{\code{I}:} {The original value supplied to the \code{I} argument.}
\item{\code{J}:} {The original value supplied to the \code{J} argument.}
\item{\code{summary.data}:} {The original value supplied to the \code{summary.data} argument.}
\item{\code{X.sq.S}:} {The modified Pearson statistic; NA if not all Wi and Yj have at least one positive and one negative response.}
\item{\code{X.sq.S.ij}:} {A matrix containing the IxJ Pearson statistics.}
}

--- For \code{type = "boot"}, the primary list additionally includes \code{boot}, a list containing the following objects:
\itemize{
\item{\code{B.use}:} {The number of bootstrap resamples used.}
\item{\code{B.discard}:} {The number of bootstrap resamples discarded due to having at least one item with all positive or negative responses.}
\item{\code{p.value.boot}:} {The bootstrap p-value for the test of SPMI.}
\item{\code{p.combo.min.boot}:} {The bootstrap p-value for the minimum p-value combination method.}
\item{\code{p.combo.prod.boot}:} {The bootstrap p-value for the product p-value combination method.}
\item{\code{X.sq.S.star}:} {A numeric vector containing the modified Pearson statistics calculated for each resample.}
\item{\code{X.sq.S.ij.star}:} {A matrix containing the IxJ Pearson statistics calculated for each resample.}
\item{\code{p.combo.min.star}:} {A numeric vector containing the minimum p-value calculated for each resample.}
\item{\code{p.combo.prod.star}:} {A numeric vector containing the product p-value calculated for each resample.}
}

--- For \code{type = "rs2"}, the primary list additionally includes \code{rs2}, a list containing the following objects:
\itemize{
\item{\code{X.sq.S.rs2}:} {The Rao-Scott second-order adjusted Pearson statistic.}
\item{\code{df.rs2}:} {The degrees of freedom for testing the second-order Rao-Scott adjusted Pearson statistic.}
\item{\code{p.value.rs2}:} {The p-value based on the Rao-Scott second-order adjustment.}
}

--- For \code{type = "bon"}, the primary list additionally includes \code{bon}, a list containing the following objects:
\itemize{
\item{\code{p.value.bon}:} {The Bonferroni adjusted p-value for the test of SPMI.}
\item{\code{X.sq.S.ij.p.bon}:} {An IxJ matrix containing the Bonferroni adjusted p-values for the IxJ Pearson statistics.}
}

--- For \code{type = "all"}, the list includes all of the above objects.

--- \code{SPMI.stat} returns a list containing the following objects:
\itemize{
\item{\code{X.sq.S}:} {Defined above.}
\item{\code{X.sq.S.ij}:} {Defined above.}
\item{\code{valid.margins}:} {The number of (Wi, Yj) response tables with at least one positive and one negative response for Wi and Yj.}
}
}
\references{
Bilder, C. and Loughin, T. (2004)  Testing for marginal independence between two categorical variables with multiple responses.  \emph{Biometrics}, \bold{36}, 433--451.
}
\seealso{
The \code{\link{genloglin}} function offers a generalized loglinear modeling approach for testing the relationship among two or three MRCVs.
}
\examples{
# Test for SPMI using the second-order Rao-Scott adjustment
test.stat.rs2 <- SPMI.test(data = farmer2, I = 3, J = 4, type = "rs2")
test.stat.rs2

# Test for SPMI using all three approaches
# A small B is used for demonstration purposes; normally, a larger B should be used
\dontrun{
test.stat.all <- SPMI.test(data = farmer2, I = 3, J = 4, type = "all", B = 99, 
    plot.hist = TRUE)
test.stat.all}

# Use SPMI.test() with summary data
# Convert raw data file to summary file for this example 
farmer.irdframe <- item.response.table(data = farmer2, I = 3, J = 4, create.dataframe = 
    TRUE)

# Test for SPMI using the Bonferroni adjustment
test.stat.bon <- SPMI.test(data = farmer.irdframe, I = 3, J = 4, type = "bon", 
    summary.data = TRUE)
test.stat.bon
}