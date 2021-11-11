\name{mvwgaim-package}
\alias{mvwgaim-package}
\docType{package}
\title{Whole genome multi-section analysis using ASReml-R V4}
\description{
This package provides a flexible computational implementation of
multi-section QTL analysis using new technological advances in ASReml-R
V4 functionality. Here, "multi-section" is a universal term encapsulating
three possible different types, multi-environment, multi-variate or multi-treatment.  
}
\details{
\tabular{ll}{
Package: \tab mvgaim \cr
Type: \tab Package\cr
Version: \tab 0.99\cr
Date: \tab 2021-11-10\cr
License: \tab GPL 2\cr
}
Welcome to version 0.99 of \pkg{mvwgaim}! The documentation given in this help file is
only brief and users should consult individual help files for the main functions of the
package. The package requires users have an installed and fully licensed
version of ASReml-R V4 (see
\url{https://www.vsni.co.uk/software/asreml-r} for more information).

It should be noted the package contains no functions for preparing your
genetic marker data for analysis. This can be done externally in many
ways but the genetic data object \code{genObj} to be passed into
\code{\link{mvwgaim.asreml}} should be a \pkg{qtl} \code{"cross"}
object. Within each of the chromosomes of \code{genObj} there at least needs to be an
\code{imputed.data} or \code{interval.data} element containing complete
(no missing values) set of markers. It may be helpful to consult the
\code{cross2int()} function in the \pkg{wgaim} package.

Multi-section QTL analysis is conducted using the function
\code{\link{mvwgaim.asreml}} which, as its first two arguments, requires
multi-section \code{asreml} base models (see the help dcumentation of
\code{mvwgaim.asreml} for more information on these models. For a whole genome set of genetic
markers in \code{genObj}, the function builds on the supplied base model
by including an additive relationship matrix into the variance
of a random multi-section by marker term in the model. The mvwgaim algorithm
then uses forward selection to choose significant markers to place in the
fixed component of the model. This forward selection prociedure avoids
computationally cumbersome marker scanning methods and subsequent thresholding
techniques. 

A summary method is available for the returned \code{"mvwgaim"}
object and provide users with a detailed report on the putative QTL
found by the algorithm including their size, the closest linked markers,
significance (including LOD score if desired). Outlier statistics of
each iteration can be graphically displayed using the \code{outStat}
function in the \pkg{wgaim} R package and putative QTL are overlayed at
the selected peaks.

It is expected this function will be integrated into a wider QTL/GWAS/GS
type package in the future. Stay tuned!
}
\author{
Author: Julian Taylor & Ari Verbyla
Maintainer: Julian Taylor <julian.taylor@adelaide.edu.au>
}
\keyword{package}
\seealso{
\code{\link{qtl-package}} 
}
