\name{wgaim-internal}
\alias{getQTL}
\alias{waldTest}
\alias{waldTest.asreml}
\alias{mvwgaim.default}
\alias{fa.modify}
\alias{favar.init}
\alias{pchisq.mixture}
\alias{qchisq.mixture}
\alias{qtlMSelect}
\title{Internal wgaim functions}
\description{
These are internal \code{gwasreml} functions that are not intended to be
directly called by the user.
}
\usage{
getQTL(object, genObj)
waldTest(object, \ldots)
\method{waldTest}{asreml}(object, cc, keep.fac = TRUE, ...)
\method{mvwgaim}{default}(baseDiag, \ldots)
fa.modify(model, Trait, by)
favar.init(char1, char2, G.param)
pchisq.mixture(x, ntrait=2)
qchisq.mixture(prob, ntrait=2, trace=F, maxiter=10)
qtlMSelect(asm, phenoData, genObj, gen.type, selection, n.fa, Trait,exclusion.window, state, verboseLev)
igrep(patterns, x, \ldots)
}
\author{Julian Taylor}
\keyword{internal}