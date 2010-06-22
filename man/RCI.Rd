\name{RCI}
\Rdversion{1.1}
\alias{RCI}
\alias{plot.hci}
\alias{confint.hci}
\alias{default.hci}
\alias{print.hci}
\alias{summary.hci}
\alias{concentrationIndex}
\alias{varianceConcentrationIndex}
\title{
A function to calculate the relative concentration index
}
\description{
This function calculates the relative concentration index due to Kakwani \emph{et al}.
It returns an object of class \code{hci} for which confidence intervals, summaries and 
plots are defined.
}
\usage{
RCI(x, y, wt = NULL)
concentrationIndex(object)
\method{plot}{hci}(x, ...)
\method{print}{hci}(x, ...)
\method{summary}{hci}(object, ...)
\method{confint}{hci}(object, parm = NULL, level = 0.95, ...)
varianceConcentrationIndex(object)
}
\arguments{
  \item{x}{
For \code{RCI}, a numeric variable containing the health measure. For the \code{print} and \code{plot} methods: an object of class \code{hci}.
}
  \item{y}{
a numeric variable containing the wealth measure.
}
  \item{wt}{
a numeric variable containing the sample weights. Defaults to \code{NULL}, which is equal to a constant sample weight for all cases.
}
  \item{object}{
an object of class \code{hci} - the result of a call to \code{RCI}.
}
  \item{parm}{
always set to \code{NULL}; provided for consistency with the generic \code{confint}.  
}

  \item{level}{
confidence level at which confidence intervals should be calculated. Defaults to 0.95.
}
  \item{...}{
arguments to pass on to \code{plot}, \code{print} and \code{summary}.
}
}
\details{
\code{RCI} returns the relative concentration index due to Kakwani \emph{et al} in an object of class \code{hci}.
A number of extractor functions have been defined. The \code{plot} function returns a plot of the Lorenz curve.
The relative concentration index is defined as twice the surface between this curve and the diagonal. \code{RCI} does
not check whether the health measure has the same sign for all cases. If this is not true the relative concentration
index is not bounded between [-1,+1]. It is up to the user to handle this, e.g. using the \code{\link{correctSign}} function.
}
\value{
For \code{RCI}, an object of class \code{hci}; this is a list containing the following components:
\item{concentrationIndex}{relative concentration index}
\item{varC}{variance of the concentration index}
\item{fractionalRank}{weighted fractional rank for each case}
\item{outcome}{health measure for each case}
\item{n}{number of cases}
\item{call}{the original call to \code{RCI}}
For \code{plot.hci}, a concentration plot.
}
\references{
Kakwani \emph{et al.}
}
\author{
Peter Konings and Niko Speybroeck
}

\seealso{
\code{\link{correctSign}}
}
\examples{
data(kenya)
example <- RCI(x = kenya$wealth, y = kenya$stunting)
plot(example)
summary(example)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univar}
\keyword{hplot}
\keyword{methods}
\keyword{print}
\keyword{misc}
