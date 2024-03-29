% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/WritePin.R
\name{WritePin}
\alias{WritePin}
\title{Write .pin file}
\usage{
WritePin(ExName, ADMBDir, Vals)
}
\arguments{
\item{ExName}{An object of class \code{character} that specifies the name of
the ADMB exectuable (without extension).}

\item{ADMBDir}{An object of class \code{character} that specifies the
location of the ADMB executable.}

\item{Vals}{An object of class \code{numeric} of length 3 that specifies the
starting values.}
}
\value{
\code{NULL} if successfully exectuted and .pin file written. Otherwise error message.
}
\description{
A function to write the .pin value (initials values) for ADMB.
}
\examples{
\dontrun{

}
}
\author{
Adrian Hordyk
}
\seealso{
\code{\link{}}
}

