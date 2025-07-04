% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tlets.R
\name{gen_maq}
\alias{gen_maq}
\title{Generate a TL-MA(\code{q}) time series}
\usage{
gen_maq(n, thetas)
}
\arguments{
\item{n}{the length of the desired time series}

\item{thetas}{the \code{vector} of MA(\code{q}) parameters}
}
\value{
a \code{vector}. The length-\code{n} times series with MA(\code{q}) parameters
\code{thetas}.
}
\description{
This function generates a TL-MA(\code{q}) time series of length \code{n} where \code{q}
is determined by \code{length(thetas)}.
}
\examples{
out <- gen_maq(n = 1200, thetas = c(0.8, 0.2, 0.3, 0.1))
}
\references{
\insertRef{mhatre2024arma}{tpdmethods}
}
