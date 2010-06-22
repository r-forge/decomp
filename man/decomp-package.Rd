\name{decomp-package}
\Rdversion{1.1}
\alias{decomp-package}
\alias{decomp}
\docType{package}
\title{
A package to calculate and decompose health inequalities
}
\description{
This package provides various ways to calculate and decompose health inequalities.
}
\details{
\tabular{ll}{
Package: \tab decomp\cr
Type: \tab Package\cr
Version: \tab 0.3\cr
Date: \tab 2009-08-26\cr
License: \tab GPL 2.0 or later\cr
LazyLoad: \tab yes\cr
Suggests: \tab gplots \cr
}
}
\author{Peter Konings <peter.konings@esat.kuleuven.be> and Niko Speybroeck <nico.speybroeck@itg.be>
}
\references{
Konings \emph{et al.}, 2009
Speybroeck \emph{et al.}, 2009
}

\keyword{ package }

\examples{
data(kenya)
kenya.rci <- RCI(x = kenya$stunting, y = kenya$wealth)
summary(kenya.rci)
plot(kenya.rci)
simple.model <- lm(stunting ~ wealth + sex + urban + meducY, data =  kenya)
simple.contrib <- contribution(simple.model, kenya$wealth)
simple.contrib
plot(simple.contrib)
}
