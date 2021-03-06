% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/decomp.r
\name{movemed}
\alias{movemed}
\title{\code{movemed} smoothes a time series by moving medians}
\usage{
movemed(y, q)
}
\arguments{
\item{y}{the series, a vector or a time series}

\item{q}{scalar,  span of moving median}
}
\value{
g  vector,  smooth component
}
\description{
\code{movemed} smoothes a time series by moving medians
}
\examples{
data(BIP)
g <- movemed(GDP,12) 
\donttest{ plot(GDP) ; t <- seq(from = 1970, to = 2009.5,by=0.25) ; lines(t,g,col="red") } 
}
