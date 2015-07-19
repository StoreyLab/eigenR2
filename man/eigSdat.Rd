% --- Source file: man/eigSdat.Rd ---
\name{eigSdat} 
\docType{data} 
\alias{eigSdat}
\title{A simulated data set with independent variables of interest and a gene expression matrix}
\description{This simulated data set consists of two parts. The first
  part "\code{varList}" contains three variables: age, genotype and ID. The second part "\code{exp}" is an expression matrix, in which all the gene expression levels are affected by age and genotype with varying effect sizes.}

\usage{data(eigSdat)}

\format{"\code{varList}" is a data set consisting of 50 rows and 3 columns, 
 each row is one observation, and "\code{exp}" is a 200 x 50 gene 
 expression matrix with genes in rows and arrays in columns. }

\references{
  Chen LS and Storey JD (2008) Eigen-R2 for dissecting variation in high-dimensional studies.
}
\keyword{datasets}