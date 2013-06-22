\name{MRCV-package}
\alias{MRCV-package}
\alias{MRCV}
\docType{package}
\title{
Statistical Methods for Analyzing the Association Among Two or Three MRCVs
}
\description{
The MRCV package provides functions for analyzing the association between two or three multiple response categorical variables (MRCVs).  A modified Pearson chi-square statistic can be used to test for marginal independence betwen two MRCVs, or a more general loglinear modeling approach can be used to compare various other structures of association among two or three MRCVs.  Bootstrap- and asymptotic-based standardized residuals and model-predicted odds ratios are available, in addition to other descriptive information.
}
\details{
\tabular{ll}{
Package: \tab MRCV\cr
Version: \tab 0.1-0\cr
Date: \tab 2013-06-22\cr
Depends: \tab tables, R (>= 3.0.0)\cr
LazyData: \tab TRUE\cr
License: \tab GPL (>= 3)\cr
}

\strong{Notation:}\cr  
Define row variable, W, column variable, Y, and strata variable, Z, as MRCVs with binary items (i.e., categories) Wi for i = 1, \ldots, I, j = 1, \ldots, J, and k = 1, \ldots, K, respectively.  Also, define a marginal count as the number of subjects who responded (Wi=a, Yj=b, Zk=c) for a, b, and c belonging to the set \{0, 1\}.

\strong{Format of Data Frame:}\cr
Many of the functions require a data frame containing the raw data structured such that \emph{n} rows correspond to the individual item response vectors, and the columns correspond to the binary items, W1, \ldots, WI, Y1, \ldots, YJ, and Z1, \ldots, ZK (in this order).  The third set of items is only necessary when analyzing the relationship among three MRCVs.  Some of the functions use a summary version of the raw data frame (converted automatically without need for user action) formatted to have 2Ix2J rows (or 2Ix2Jx2K rows in the three MRCV case) and 5 columns generically named \code{W}, \code{Y}, \code{wi}, \code{yj}, and \code{count} (or 7 columns generically named \code{W}, \code{Y}, \code{Z}, \code{wi}, \code{yj}, \code{zk}, and \code{count} in the three MRCV case). The column named \code{count} contains the marginal counts defined above.

\strong{Descriptive Functions:}\cr
Users can call the \code{\link{item.response.table}} function to obtain a cross-tabulation of the positive and negative responses for each combination of items, or the \code{\link{marginal.table}} function to obtain a cross-tabulation of only the positive responses.

\strong{Functions to Test for Marginal Independence:}\cr
Methods proposed by Agresti and Liu (1999), Bilder and Loughin (2004), and Thomas and Decady (2004) are implemented using the \code{\link{SPMI.test}} function.  This function calculates a modified Pearson chi-square statistic that can be used to test for simultaneous pairwise marginal independence (SPMI) between two MRCVs.  SPMI is a test of whether each Wi is pairwise independent of each Yj.  The modified statistic is the sum of the IxJ Pearson statistics used to test for independence of each (Wi, Yj) pair.  The asymptotic distribution of the modified statistic is a linear combination of independent chi-square(1) random variables, so traditional methods for analyzing the association between categorical variables W and Y are inappropriate.  The \code{\link{SPMI.test}} function offers three methods, a nonparametric bootstrap approach, a Rao-Scott second-order adjustment, and a Bonferroni adjustment, that can be used in conjunction with the modified statistic to construct an appropriate test for independence.

\strong{Functions for Performing Regression Modeling:}\cr
Regression modeling methods described by Bilder and Loughin (2007) are implemented using  \code{\link{genloglin}} and methods \code{\link{summary.genloglin}}, \code{\link{residuals.genloglin}}, \code{\link{anova.genloglin}}, and \code{\link{predict.genloglin}}.  The \code{\link{genloglin}} function provides parameter estimates and Rao-Scott adjusted standard errors for models involving two or three MRCVs.  The \code{\link{anova.genloglin}} function offers second-order Rao-Scott and bootstrap adjusted model comparison and goodness-of-fit (Pearson and LRT) statistics.   The \code{\link{residuals.genloglin}} and \code{\link{predict.genloglin}} functions provide bootstrap- and asymptotic-based standardized Pearson residuals and model-based odds ratios, respectively.
}
\author{
Natalie Koziol, Chris Bilder

Maintainer: Natalie Koziol <nak371@gmail.com>
}
\references{
Agresti, A. and Liu, I.-M. (1999)  Modeling a categorical variable allowing arbitrarily many category choices.  \emph{Biometrics}, \bold{55}, 936--943.

Bilder, C. and Loughin, T. (2004)  Testing for marginal independence between two categorical variables with multiple responses.  \emph{Biometrics}, \bold{36}, 433--451.

Bilder, C. and Loughin, T. (2007)  Modeling association between two or more categorical variables that allow for multiple category choices.  \emph{Communications in Statistics--Theory and Methods}, \bold{36}, 433--451.

Thomas, D. and Decady, Y. (2004)  Testing for association using multiple response survey data: Approximate procedures based on the Rao-Scott approach.  \emph{International Journal of Testing}, \bold{4}, 43--59.
}