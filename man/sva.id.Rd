% --- Source file: man/sva.id.Rd ---
\name{sva.id}
\alias{sva.id}
\title{Estimate the number of significant eigenvectors to include in analysis}
\description{
  Estimate the number of significant eigenvectors from the residuals of
  a high-dimensional data set given the model fit. 
}
\usage{
  sva.id(dat, mod, B=20, eigen.sig=0.10, seed=NULL)
}
\arguments{
  \item{dat}{A m response variables by n samples matrix of data}
  \item{mod}{A n by k model matrix corresponding to the primary model fit, 
  if mod is NULL, we estimate the number of significant eigenvectors in
  the expression data. (see model.matrix)}
  \item{B}{The number of null iterations to perform.}
  \item{eigen.sig}{The significance cutoff for eigenvectors.}
  \item{seed}{A numeric seed for reproducible results. Optional.} 
}
\details{
  Note that this function is a modified function from the package 
  \code{sva} by Leek JT and Storey JD. 
  The model matrix should include a column for an intercept. sva.id 
  estimates the number of surrogate variables to include in the analysis
  as described in Leek and Storey (2007).
}


\value{
  A list containing:
  \item{n.sv}{The number of significant surrogate variables}
  \item{p.sv}{The p-values for each eigenvector} 

}
\references{
   Leek JT and Storey JD (2007) Capturing heterogeneity in gene expression   
   studies by "surrogate variable analysis." \emph{PLoS Genetics}, 
   \bold{3}:e161. 
}

\author{Jeff T. Leek \email{jtleek@gmail.com} and John D. Storey \email{jstorey@princeton.edu}}

\seealso{\code{\link{eigenR2}}}

\keyword{misc}