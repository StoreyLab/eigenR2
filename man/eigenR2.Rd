% --- Source file: man/eigenR2.Rd ---
\name{eigenR2}
\alias{eigenR2}
\title{Eigen-R2 for Dissecting Variation in High-dimensional Studies}
\description{
  Eigen-R2 is a high-dimensional version of the classic R-square
  statistic.  It can be applied when one wants to determine the
  aggregate R-square value for many related response variables
  according to a common set of independent variables.
}
\usage{
  eigenR2(dat, model, null.model=NULL, adjust=FALSE, eigen.sig=NULL, mod.fit.func=NULL)
}

\arguments{
  \item{dat}{The matrix of related response variables, where each row
  represents a distinct response variable.}
  \item{model}{A model matrix describing the fit of each response
  variable to the independent variables (see \code{\link{model.matrix}}). }
  \item{null.model}{If different than the simple mean, an optional null
  model to be used in calculating the denominator of the R-square statistic. }
  \item{adjust}{If \emph{TRUE}, eigen-R2 makes the analogous small
  sample adjustment avaiable for the R-square statistic. }
  \item{eigen.sig}{An optional significance threshold to select the
  right eigen-vectors explaining more variation than expected by
  chance for use in the eigen-R2 calculation. If not \emph{NULL},
  the value of the p-value threshold to apply. }
  \item{mod.fit.func}{An optional function to estimate R-square. By default, the
  sum of squares will be minimized with respect to the model.}
}

\details{
  This function can be used for determining the aggregate R-square value for many related response variables
  according to a common set of independent variables.

  This function can also be used to compute conditional eigen-R2,
  by specifying the \code{null.model}. See examples.

  By specifying the significance threshold \code{eigen.sig}, only
  the statistically significant right eigenvectors will be used to compute eigen-R2. 

  By default, \code{mod.fit.func} is \emph{NULL} and the function will compute 
  eigen-R2s with simple linear regression. If \code{mod.fit.func} is not
  \emph{NULL}, an R-square estimation function should be provided. The 
  input of the function should be the response variables and the function
  returns the estimated R2 for the response variables.

}
\value{
  The function returns a list containing:
  \item{eigenR2}{The eigen-R2 estimate for the model. }
  \item{weights}{The proportion of total variation of the response
    variables that each eigenvector captures. }
  \item{eg.R2s}{The R-square for each significant eigenvectors. }
  \item{p.eg}{The p-value for each eigenvector, testing whether it
    explains more variation among the response variables than expected
    by chance. }

}
\references{
  Chen LS and Storey JD (2008) Eigen-R2 for dissecting variation 
  in high-dimensional studies.
}
\author{Lin S. Chen \email{lschen@princeton.edu} and John D. Storey \email{jstorey@princeton.edu}}
\seealso{\code{\link{plot.eigenR2}}}

\examples{
  \dontrun{
  ## Load the simulated data.
  ## The data set contains age, genotype and IDs for 50 arrays.
  ## The expression matrix is a 200 genes by 50 array matrix
  data(eigSdat)
  exp <- eigSdat$exp
  exp <- t(apply(exp, 1, function(x) x-mean(x)))
  varList <- eigSdat$varList
  age <- varList$age
  geno <- varList$geno
  ID <- varList$ID

  ## the eigen-R2 for variable "age"
  mod1 <- model.matrix(~1+age)
  eigenR2.age <- eigenR2(dat = exp, model = mod1) 

  ## the eigen-R2 for variable "age"
  ## adjust for small sample size
  eigenR2.age.adj <- eigenR2(dat=exp, model=mod1, adjust=TRUE)
  ## Or specifying the significance threshold of eigen-genes
  eigenR2.age.adj2 <- eigenR2(dat=exp, model=mod1, eigen.sig=0.01)

  ## the conditional eigen-R2 for "genotype" given "age"
  mod2 <- model.matrix(~1+age+as.factor(geno))
  eigenR2.g <- eigenR2(dat=exp, model=mod2, null.model=mod1)

  ## use "mod.fit.func" to specify the estimation function
  func1 <- function(y) {summary(lm(y~mod1))$r.squared}
  eigenR2.age2 <- eigenR2(dat=exp, mod.fit.func=func1)

  ## use linear mixed effect model to estimate R2s
  library(nlme)
  func2 <- function(y) {
      m1 <- lme(y~1+age+as.factor(geno), random=~1|ID)
      r2 <- 1-sum(resid(m1)^2)/sum((y-mean(y))^2)
      return(r2)
  }
  eigenR2.ag <- eigenR2(dat=exp, mod.fit.func=func2)

  ## Plot the information on the eigen-genes
  plot.eigenR2(eigenR2.ag)
  }
}
\keyword{misc}