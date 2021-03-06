% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doGRAS.R
\name{doGRAS.long}
\alias{doGRAS.long}
\title{Do the GRAS algorithm for a "long" data.table}
\usage{
doGRAS.long(dt, rowcol, colcol, rowsum, colsum, ...)
}
\arguments{
\item{dt}{the "long" data.table}

\item{rowcol}{the column which holds the row-dimension of the reshaped matrix}

\item{colcol}{the column which holds the column-dimension of the reshaped matrix}

\item{rowsum}{vector with the row totals}

\item{colsum}{vector with the column totals}

\item{...}{parameter to be passed to doGRAS}
}
\value{
the updated data.table in the same format as before
}
\description{
This function calculates an updated "long" data.table, based on a data.table dt,
which meets the given row and columns totals. Note: dt can contain
negative elements.
}
\references{
Junius T. and J. Oosterhaven (2003), The solution of
    updating or regionalizing a matrix with both positive and
    negative entries, Economic Systems Research, 15, pp. 87-96.

    Lenzen M., R. Wood and B. Gallego (2007), Some comments on the
    GRAS method, Economic Systems Research, 19, pp. 461-465.

    Temurshoev, U., R.E. Miller and M.C. Bouwmeester (2013), A note
    on the GRAS method, Economic Systems Research, 25, pp. 361-367.
}
\author{
Oliver Reiter
}
\keyword{gras,}
\keyword{matrix}
\keyword{updating}
