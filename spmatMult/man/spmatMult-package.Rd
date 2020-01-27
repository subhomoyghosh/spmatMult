\name{spmatMult-package}
\alias{spmatMult-package}
\alias{spmatMult}
\docType{package}
\title{
Efficient sparse matrix multiplication using multicores
}
\description{
 Current R's Blas or Openblas doesn't provide multithreaded sparse-matrix multiplication. This is a sparse matrix multiplication package in RcppArmadillo and R using OpenMp.
}
\details{

This is an adaptation of Yadav et. al. 's \code{\link[https://www.geosci-model-dev-discuss.net/gmd-2016-204/]{technical note}}. Generally, this code is an adaptation to FG Gustavson's sparse matrix multiplication in RcppArmadillo so that R users can leverage this too!

Usage of the functions in this package from R assumes Matrix library is loaded.

Note: You need to have RcppArmadillo configured with OpenMp. For Mac user, Coatless has an awesome repo  \code{\link[https://github.com/coatless/r-macos-rtools/] here}.

}

\references{
Technical Note: Improving the computational efficiency of sparse matrix multiplication in linear atmospheric inverse problems, Yadav et al. 
}

}

\author{

Subhomoy Ghosh

Maintainer: Subhomoy Ghosh
}



\examples{

## Generate two random matrices

nrA<- 1000; nrB<- 10000
ncA<- 10000; ncB<- 5000
p1<- .001; p2<- .05

A<- rsparsematrix(nrA,ncA,p1);
B<- rsparsematrix(nrB,ncB,p2)

## Multiply A and B
# start with something close but small, function will automatically 
# adjust sparsity density

isProdSym<- 0;
dens<- .01;  
nthreads<- 4

tic<- Sys.time()
C1<- spsp_prodSp_openmp(A,B,isProdSym,dens,nthreads)
toc<- Sys.time()
toc-tic

tic<- Sys.time()
C2<- A \%*\% B
toc<- Sys.time()
toc-tic

all.equal(C1,C2)


}
