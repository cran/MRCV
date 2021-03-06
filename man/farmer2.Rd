\name{farmer2}
\alias{farmer2}
\docType{data}
\title{
Data for Two MRCVs from the Kansas Farmer Survey
}
\description{
Responses for two MRCVs from a survey of Kansas farmers.  This data was first examined by Bilder and Loughin (2007).
}
\usage{farmer2}
\format{
The data frame contains the following 7 columns:

Columns 1-3 correspond to the response categories for the question "Which of the following do you test your swine waste for?"  Binary responses (1 = Yes, 0 = No) are provided for each category.
\itemize{
\item{\code{W1}:} {Nitrogen}
\item{\code{W2}:} {Phosphorous}
\item{\code{W3}:} {Salt}
}
Columns 4-7 correspond to the response categories for the question "What swine waste disposal methods do you use?"  Binary responses (1 = Yes, 0 = No) are provided for each category.
\itemize{
\item{\code{Y1}:} {Lagoon}
\item{\code{Y2}:} {Pit}
\item{\code{Y3}:} {Natural drainage}
\item{\code{Y4}:} {Holding tank}
}
}
\source{
Bilder, C. and Loughin, T. (2007)  Modeling association between two or more categorical variables that allow for multiple category choices.  \emph{Communications in Statistics--Theory and Methods}, \bold{36}, 433--451.
}