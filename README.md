
### R/mvwgaim: An R package for efficient multi-(environment/variate/treatment) whole genome QTL analysis

**Authors**: Julian Taylor & Ari Verbyla

This is the public facing GitHub repository version of the R package mvwgaim.

**R/wgaim** is a multi-(environment/variate/treatment) whole genome average
interval mapping R package that implements the multi-environment wgaim algorithm
derived in Verbyla et al. (2007, 2012). The packages main QTL analysis function uses ASReml-R V4 for its core linear mixed modelling. To use full functionality of the package users will require a valid license for ASReml-R V4 and this can be obtained from [https://www.vsni.co.uk/software/asreml-r](https://www.vsni.co.uk/software/asreml-r). 

To install the package from GitHub you will need to do the following: 

1. Install the [devtools](https://cran.r-project.org/package=devtools) package. Do this by invoking R and then typing


```r
install.packages("devtools")
```

2. Install wgaim using 


```r
devtools::install_github("DrJ001/mvwgaim")
```

#### Getting Started

For a quick but complete introduction of the functionality of the package please
visit the help files of the package.

#### References

Verbyla, A.P., Cullis, B.R. & Thompson, R. (2007) The analysis of QTL by simultaneous use of the of the full linkage map. *Theoretical and Applied Genetics*, **116**, 95-111.

Verbyla, A. P., & Cullis, B. R. (2012) Multivariate whole genome average
interval mapping: QTL analysis for multiple traits and/or environments. *Theoretical and applied genetics* **125**, 933-953.
