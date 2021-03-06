% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_misc.R
\name{genW}
\alias{genW}
\alias{rowStdrt}
\alias{constructW}
\title{Standard spatial contiguity matrixes}
\usage{
genW(n, type = "rook", seed = NULL)

rowStdrt(W)

constructW(coords, labels)
}
\arguments{
\item{n}{a number of objects with spatial interaction to be arranged.See 'Details' for objects arranging principle}

\item{type}{an optional type of spatial interaction. Currently 'rook' and ''queen' values are supported, to produce Rook and Queen Contiguity matrix. See references for more info. 
By default set to rook.}

\item{seed}{an optional random number generator seed for random matrices}

\item{W}{a spatial contiguity matrix to be standatised}

\item{coords}{a matrix of two columns, where every row is a longitude-latitude pair of object coordinates}

\item{labels}{a vector of object lables to mark rows and columns of the resulting contiguity matrix}
}
\description{
\code{genW} generates an spatial contiguity matrix (rook or queen)

\code{rowStdrt} standartizes spatial contiguity matrix by rows

\code{constructW} contructs a spatial contiguity matrix using object longitude and latitude coordinates
}
\details{
To generate spatial interaction between \code{n} objects the function arranges them on a chess board. 
A number of columns is calculated as a square root of \code{n}, rounded to the top. The last row contains empty cells, if n is not quadratic

The function divides every element in an argument matrix by the sum of elements in its row. Some spatial estimation requires this standartisation (generally - for faster calculations)

The function contructs a spatial contiguity matrix 
using object longitude and latitude coordinates. Eucledean distance is currently used.
}
\examples{
# Completely filled 10x10 rook contiguity matrix
rookW <- genW(100)
rookW

# Partly filled 10x10 rook contiguity matrix
rookW <- genW(90)
rookW

# Completely filled 10x10 queen contiguity matrix
queenW <- genW(100, type="queen")
queenW

# Completely filled 10x10 queen contiguity matrix
queenW <- genW(100, type="queen")
queenW

# Standartisation
stQueenW <- rowStdrt(queenW)
stQueenW


data(airports)

W <- constructW(cbind(airports$lon, airports$lat),airports$ICAO_code)

}
\references{
Anselin, L. (1988). Spatial Econometrics: Methods and Models. Kluwer Academic Publishers, Dordrecht, The Netherlands.
}
\keyword{frontier}
\keyword{spatial}
\keyword{stochastic}
