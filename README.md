# spmatMult

This is a R package that provides efficient sparse sparse matrix multiplication using parallel architecture. It uses 
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

### Install instruction
devtools::install_github('subhomoyghosh/spmatMult/spmatMult')


### Minimal example 

```@{r}
#### Generate two random matrices

nrA<- 1000; nrB<- 10000
ncA<- 10000; ncB<- 5000
p1<- .001; p2<- .05

A<- rsparsematrix(nrA,ncA,p1);
B<- rsparsematrix(nrB,ncB,p2)

#### Multiply A and B
#### start with something close but small, function will automatically 
#### adjust sparsity density

isProdSym<- 0;
dens<- .01;  
nthreads<- 4

tic<- Sys.time()
C1<- spsp_prodSp_openmp(A,B,isProdSym,dens,nthreads)
toc<- Sys.time()
toc-tic

tic<- Sys.time()
C2<- A %*% B
toc<- Sys.time()
toc-tic

all.equal(C1,C2)
```
