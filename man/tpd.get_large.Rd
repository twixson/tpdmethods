% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tpd.R
\name{tpd.get_large}
\alias{tpd.get_large}
\title{Get large points from a bivariate \code{matrix}}
\usage{
tpd.get_large(data, thresh = 0.975, r0 = 0)
}
\arguments{
\item{data}{a bivariate \code{matrix} of points}

\item{thresh}{the radial quantile used to distinguish large points from
points in the bulk of the distribution. (default is 0.975)}

\item{r0}{instead of setting a quantile the user may set a radial value. If
\code{r0} is greater than zero the \code{thresh} is ignored. (default is 0).}
}
\value{
a bivariate \code{matrix} of large points.
}
\description{
This function extracts the large points from a bivariate \code{matrix}
}
\examples{
myData <- matrix(evd::rfrechet(1000, shape = 2), ncol = 2)
out <- tpd.get_large(myData)
}
