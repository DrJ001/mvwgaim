\name{summary.mvwgaim}
\alias{summary.mvwgaim}
\title{Summary and print methods for the class \code{"mvwgaim"}.
}
\description{
Print a QTL summary from the \code{"mvwgaim"} object in presentable
format.
}
\usage{
\method{summary}{mvwgaim}(object, genObj, LOD = TRUE, \ldots)
}
\arguments{
  \item{object}{an object of class \code{"mvwgaim"} (see Details)}
  \item{genObj}{a genetic data structure inheriting the class \code{"cross"}}
  \item{LOD}{logical value. If TRUE LOD scores for QTL are calculated, defaults to \code{TRUE}}
  \item{\ldots}{further arguments passed to or from other methods}
}
\details{
  It is important that the \code{genObj} is not missing in
  \code{summary.mvwgaim()} as it contains vital summary information about
  each of the QTL detected.
}
\value{A summary of the QTL component of the \code{"mvwgaim"} object is
  printed to the screen. For each QTL detected, if an \code{"interval"}
  analysis was performed then \code{summary.mvwgaim()} prints the
  chromosome, name and distance of each flanking marker, size,
  p-value and LOD score if desired. If a \code{"marker"} analysis was
  performed then the chromosome, name and distance of the associated
  marker, size, p-value and LOD score are printed. 
}
\author{Julian Taylor
}
\seealso{\code{\link{mvwgaim.asreml}}}
\examples{
\dontrun{
# convert "cross" object to "interval" object

mapi <- cross2int(map, impute = "Martinez", id = "Genotype")

# baseline multi-environment model

hwt.fa2 <- asreml(hectlitwt_500 ~ 1, random = ~ fa(Expt, 2):Genotype + at(Expt):Range
                  + at(Expt):Row, residual = ~ dsum(~ ar1(Row):ar1(Range) | Expt),
                  G.param = hwt.fa1$G.param, R.param = hwt.fa1$R.param,
                  data = heath, maxit = 30)

# Run whole genome mvwgaim

hwt.wg <- mvwgaim(hwt.asdo, hwt.fa2, mapi, Trait = "Expt", merge.by =
          "Genotype", method = "fixed", gen.type = "interval", n.fa = 2,
          na.action = na.method(x = "include")) 

# summarise

summary(hwt.wg, mapi)
}
}
\keyword{regression}
\keyword{methods}
