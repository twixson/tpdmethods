% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pca.R
\name{plot_eigen}
\alias{plot_eigen}
\title{Create a plot of several eigenvectors.}
\usage{
plot_eigen(eigen_vecs, num_vecs = 5, var_names = NULL)
}
\arguments{
\item{eigen_vecs}{The \code{matrix} of eigenvectors.}

\item{num_vecs}{the eigenvectors you want to plot. If you input a single
value (n) then eigenvectors 1:n will be plotted. If you input a vector
then eigenvectors in that vector will be plotted. (default is 5)}

\item{var_names}{The names of the variables, used for making nice x-axis
tick mark labels.}
}
\value{
a \code{ggplot2} plot object
}
\description{
This function is a \code{ggplot2} wrapper which plots several eigenvectors weights
on a common binned scale which is automatically computed.
}
\examples{
myData <- matrix(evd::rfrechet(625, shape = 2), ncol = 5)
myTPDM <- tpd(myData)
myEigen <- get_eigen(myTPDM)
out <- plot_eigen(myEigen$vectors, num_vecs = 3, var_names = 1:5)
}
