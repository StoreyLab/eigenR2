% --- Source file: man/eigenR2plot.Rd ---
\name{eigenR2plot}
\alias{eigenR2plot}
\alias{plot.eigenR2}
\title{Graphical display of eigenR2 objects}
\description{
  Graphical display of the information captured the eigenvectors.
}

\usage{
  eigenR2plot(eigenR2obj)
  plot.eigenR2(x,...)
}

\arguments{
  \item{eigenR2obj, x}{An eigenR2 object.}
  \item{...}{Plotting arguments. }
}

\details{
  The function graphically displays the information captured by the eigenvectors, 
  including:
  \enumerate{
    \item A plot that shows the proportion of total variation each
    eigen-vector captures.
    \item A plot of the R-square value for each significant eigenvectors.
    \item A plot of the p-values for each eigenvector, if \code{eigen.sig}
    is specified.
  }
  
}

\value{
  None.
}

\references{
  Chen LS and Storey JD (2008) Eigen-R2 for dissecting the variation 
  in high-dimensional studies.
}

\author{Lin S. Chen \email{lschen@princeton.edu} and John D. Storey \email{jstorey@princeton.edu}}
\seealso{\code{\link{eigenR2}}}
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

  ## the eigen R-square for variable "age"
  mod2 <- model.matrix(~1+age+as.factor(geno))
  eigenR2.g <- eigenR2(dat=exp, model=mod2)
  
  plot.eigenR2(eigenR2.g)
 }
}
\keyword{misc}
