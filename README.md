# spmatMult

This is a package that provides efficient sparse sparse matrix multiplication using parallel architecture. It uses 
Rcpp and RcppArmadillo. Sparse sparse matrix multiplication algorithm used here is inspired by [Yadav et al.'s](https://www.geosci-model-dev-discuss.net/gmd-2016-204/gmd-2016-204.pdf) approach  
to implement Gustavson (1978) algorithm.

There are four functions that user can use. One of the main functions is spsp_prodSp_openmp(). It has three extra arguments 
on top of the two input matrices.

1) Input matrix A
2) Input matrix B
3) isProdSym: Is the product symmetric? (0/1)
4) dens: What is the projected density of the product? (try to be less conservative here, it will be automatically adjusted) 
5) nthreads: How many threads to use? 

These arguments basically sum up the entire package description. An example is added on the package description page! 
