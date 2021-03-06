% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ibs.R
\name{ibs}
\alias{ibs}
\title{Integrals of B-Splines}
\usage{
ibs(
  x,
  df = NULL,
  knots = NULL,
  degree = 3,
  intercept = FALSE,
  Boundary.knots = NULL,
  ...
)
}
\arguments{
\item{x}{The predictor variable.  Missing values are allowed and will be
returned as they are.}

\item{df}{Degree of freedom that equals to the column number of the returned
matrix.  One can specify \code{df} rather than \code{knots}, then the
function chooses \code{df - degree - as.integer(intercept)} internal
knots at suitable quantiles of \code{x} ignoring missing values and
those \code{x} outside of the boundary.  If internal knots are specified
via \code{knots}, the specified \code{df} will be ignored.}

\item{knots}{The internal breakpoints that define the splines.  The default
is \code{NULL}, which results in a basis for ordinary polynomial
regression.  Typical values are the mean or median for one knot,
quantiles for more knots.}

\item{degree}{A nonnegative integer specifying the degree of the piecewise
polynomial. The default value is \code{3} for cubic splines. Zero degree
is allowed for piecewise constant basis functions.}

\item{intercept}{If \code{TRUE}, the complete basis matrix will be returned.
Otherwise, the first basis will be excluded from the output.}

\item{Boundary.knots}{Boundary points at which to anchor the splines.  By
default, they are the range of \code{x} excluding \code{NA}.  If both
\code{knots} and \code{Boundary.knots} are supplied, the basis
parameters do not depend on \code{x}. Data can extend beyond
\code{Boundary.knots}.}

\item{...}{Optional arguments that are not used.}
}
\value{
A numeric matrix of \code{length(x)} rows and \code{df} columns if
    \code{df} is specified or \code{length(knots) + degree +
    as.integer(intercept)} columns if \code{knots} are specified instead.
    Attributes that correspond to the arguments specified are returned
    mainly for other functions in this package.
}
\description{
Generates basis matrix for integrals of B-splines.
}
\details{
The implementation is based on the close form recursion formula.
}
\examples{
library(splines2)

x <- seq.int(0, 1, 0.01)
knots <- c(0.2, 0.4, 0.7, 0.9)
ibsMat <- ibs(x, knots = knots, degree = 1, intercept = TRUE)

## get the corresponding B-splines by bSpline()
bsMat0 <- bSpline(x, knots = knots, degree = 1, intercept = TRUE)
## or by the deriv() method
bsMat <- deriv(ibsMat)
stopifnot(all.equal(bsMat0, bsMat, check.attributes = FALSE))

## plot B-spline basis with their corresponding integrals
op <- par(mfrow = c(1, 2))
matplot(x, bsMat, type = "l", ylab = "B-spline basis")
abline(v = knots, lty = 2, col = "gray")
matplot(x, ibsMat, type = "l", ylab = "Integral of B-spline basis")
abline(v = knots, lty = 2, col = "gray")

## reset to previous plotting settings
par(op)
}
\references{
De Boor, Carl. (1978). \emph{A practical guide to splines}.
Vol. 27. New York: Springer-Verlag.
}
\seealso{
\code{\link{bSpline}} for B-splines;
\code{\link{dbs}} for derivatives of B-splines;
}
