% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/frequdom.r
\name{lagwintu}
\alias{lagwintu}
\title{\code{lagwintu} Tukey's Lag-window for indirect spectrum estimation}
\usage{
lagwintu(NL)
}
\arguments{
\item{NL}{number of lags used for estimation}
}
\value{
win vector, one-sided weights
}
\description{
\code{lagwintu} Tukey's Lag-window for indirect spectrum estimation
}
\examples{
win <- lagwintu(5)   
}
