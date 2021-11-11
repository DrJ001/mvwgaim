\name{mvwgaim.asreml}
\alias{mvwgaim.asreml}
\alias{mvwgaim}
\title{Genome wide association analysis of single and multi-section
  experiments using ASReml-R V4
}
\description{
A comprehensive marker scanning function for genome wide association
analysis of single and multi-section (environment, trait, treatment) experiments using ASReml-R V4
}
\usage{
\method{mvwgaim}{asreml}(baseDiag, baseModel, genObj, merge.by = NULL,
       fix.lines = TRUE, Trait = NULL, gen.type = "interval", method =
       "fixed", selection = "interval", n.fa = 0, exclusion.window = 20,
       breakout = -1, TypeI = 0.05, trace = TRUE, verboseLev = 0, \ldots)
}
\arguments{
\item{baseDiag}{a linear mixed model object of class \code{"asreml"} representing
  a base multi-environment \code{asreml()} model containing a Genotype
  by Environment term with a diagonal variance structure for the environment effects.  
}  
\item{baseModel}{a linear mixed model object of class \code{"asreml"} representing
  a base multi-environment \code{asreml()} model containing a Genotype
  by Environment term with a variance/covariance structure for the environment effects.  
}
\item{genObj}{a list object containing the genotypic data in the form of
  a \pkg{qtl} \code{"cross"} object (see Details).
}
 \item{merge.by}{a character string or name of the column(s) in
   \code{baseModel} phenotypic data \code{genObj} to merge the
   phenotypic and genotypic data sets.
 }
\item{Trait}{a character string column name in the phenotypic
       data used in \code{baseModel} to identify to the multiple
       sections of the data (see Details).
     }
\item{fix.lines}{a logical value. If \code{TRUE} then lines existing in
  the phenotype data that do not exist in \code{genObj} are fixed and placed in
  the \code{fixed} component of the \code{asreml()} models (see
  Details). It is recommended to set this to TRUE. Defaults to \code{TRUE}.
}
\item{gen.type}{a character string determining the type of genetic data to
  be used in the analysis. Possibilities are \code{"interval"} and
  \code{"markers"}. The default is \code{"interval"}.
}
\item{method}{a character string determining the type of algorithm to
  be used in the analysis. Possibilities are \code{"random"} and
  \code{"fixed"}. The default is \code{"random"}. (see Details).
}
\item{selection}{a character string determining the type of selection
  method that is used to select QTL in the analysis. Possibilities are
  \code{"interval"} and \code{"chromosome"}. The default is
  \code{"interval"}. (see Details).
}
\item{n.fa}{number of factors for the Factor Analytic whole genome random marker by
  \code{Trait} interaction term. Default = 0 (see Details).
}
\item{exclusion.window}{For each QTL, the distance in centimorgans on
  the left and right side of each QTL that is excluded from further analysis.
}
\item{breakout}{a numerical integer equivalent to the iteration
 where the algorithm breaks out. The default is
  -1 which ensures the algorithm scans all markers within the specified
  chromosomes. (see Details)
}
\item{TypeI}{a numerical value determining the familywise alpha level of significance
    for detecting a QTL. The default is 0.05.
  }
\item{trace}{a automatic tracing facility. If \code{trace = TRUE} then
    all \code{asreml} output is piped to the screen during the analysis.
  If \code{trace = "file.txt"}, then output from all asreml models is
  piped to "\code{file.txt}". Both trace mechanisms will display a
  message if a QTL is detected.
}
\item{verboseLev}{Numerical value, either 0 or 1, determining the level of tracing
   outputted during execution of the algorithm A 0 value will produce
   the standard model fitting output from the fitted ASReml models
   involved in the forward selection. A value of 1 will add a table of
   chromosome and interval outlier statistics for each iteration.
}
\item{\dots}{
Any extra arguments to be passed to the \code{asreml}
  calls. within the function (see \code{?asreml} and \code{?asreml.options} for more
  information).
}
}
\details{
In the initial call to \code{mvwgaim.asreml()}, the marker or interval
information is collated from \code{genObj}. If \code{gen.type =
  "interval"} then interval marker data is collated from the
\code{"interval.data"} components of the chromosomes in
\code{genObj}. If \code{gen.type = "markers"} then markers are
collated from the \code{"imputed.data"} components of the chromosomes in
\code{genObj}.

The \code{Trait} argument is used as a placeholder for potentially
different types of sectioning in the data. For example, \code{Trait} may
be the column factor describing different experiments or
environments. It may also be the column factor partitioning multiple
measurements of the same experimental unit. 

It is recommended to set \code{fix.lines = TRUE} to ensure additive and
non-additive genetic variances are estimated from lines in the
\code{merge.by} component of the phenotypic data that have genetic
marker data in \code{genObj}. Lines in the phenotype
\code{merge.by} factor not existing in \code{genObj} will be placed
as a fixed factor (called \code{Gomit}) in the \code{asreml} model.
Note, if there are others factors in the model that have some
potential confounding with \code{Gomit} then \code{asreml} will
indicate this with a simple message \code{'Terms with zero df listed
in attribute 'zerodf' of the wald table'} at the end of its iterative
maximisation. This confounding will have no effect on the outcome and
can be safely ignored. If \code{fix.lines = FALSE} is set then all
available lines in the \code{merge.by} component of the phenotypic data
will be used to estimate the non-additive genetic variance.
In this instance, users also need to be aware that \code{asreml} will
output a large number of warnings due to an inherent mismatch in the levels
of the lines contained in the phenotype data compared to the lines in
\code{genObj}.

The \code{n.fa} argument provides flexibility in specifying the Factor
Analytic model for the whole genome random marker by \code{Trait}
interaction effects. For \code{n.fa = 0} a \code{diag(Trait}) prefix term is
used and separate additive genetic variances are assumed for all levels specified
in \code{Trait}. For \code{n.fa = n} then a \code{fa(Trait, n)} prefix
term is used then a Factor Analytic model of order n is used as a
parsimonious approximation to the interaction effects. To ensure model
stability the Factor Analytic models are built sequentially. For
example, if \code{n.fa = 2} is set then the models for \code{n.fa = 0}
and \code{n.fa = 1} are fitted in sequence before the model with
\code{n.fa = 2}.

The mvwgaim algorithm is a forward selection procedure that sequentially
chooses significant intervals/markers linked to the response and places them
additively as QTL main and QTL by Trait interaction terms in the fixed
component of the model. (Note: these \code{exclusion.window} arguments
has a low numerical default to match a \code{genObj}
with centimorgan map units. These would need to be dramatically
increased if map units, for example, were base-pairs from a physical map.)
Once all significant intervals/markers have been
selected the QTL by Trait interaction term for each selected
interval/marker is tested and if significant the main QTL term is
removed. If the QTL by Trait interaction term is not
signficant the term is removed and the main QTL term is retained. 

Users can break out of the algorithm using the \code{breakout}
argument. If a numerical value greater than zero is given then the
scanning algorithm will breakout after the specified number of
iterations. This may be useful for quick investigation of a particular
markers within a chromosome. This feature will be enhanced in the future
versions of the function.

It is recommended that \code{trace = "file.txt"} be used to pipe the
sometimes invasive tracing of \code{asreml} licensing and fitting
numerics for each model to a file. Errors, warnings and messages will
still appear on screen during this process. Note some warnings that
appear may be passed through from an asreml call and are outputted upon
exit. These may be ignored as they are handled during the execution of
the function.
}
\value{An object of class \code{"mvwgaim"} which also inherits the class
  \code{"asreml"} by default. The object returned is actually an \code{asreml}
  object (see \code{asreml.object}) with the addition of components from
  the QTL detection listed below.

\item{QTL}{A list of components from the significant QTL detected
    including a character vector of the significant QTL. The complete
    set of Wald statistics are also returned for graphically displaying
    using \code{plotProfile()}.
  }
 }
 \references{
Verbyla, A. P., & Cullis, B. R. (2012) Multivariate whole genome average
interval mapping: QTL analysis for multiple traits and/or
    environments. \emph{Theoretical and applied genetics}, \bold{125},
    933-953.     
}
\author{
Julian Taylor
}
\seealso{
\code{\link{summary.mvwgaim}}
}
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
