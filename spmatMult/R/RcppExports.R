# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

rcpparma_hello_world <- function() {
    .Call(`_spmatMult_rcpparma_hello_world`)
}

rcpparma_outerproduct <- function(x) {
    .Call(`_spmatMult_rcpparma_outerproduct`, x)
}

rcpparma_innerproduct <- function(x) {
    .Call(`_spmatMult_rcpparma_innerproduct`, x)
}

rcpparma_bothproducts <- function(x) {
    .Call(`_spmatMult_rcpparma_bothproducts`, x)
}

spsp_prodSp_serial <- function(A, B, isProdSym, p) {
    .Call(`_spmatMult_spsp_prodSp_serial`, A, B, isProdSym, p)
}

spsp_prodSp_openmp <- function(A, B, isProdSym, p, nthreads) {
    .Call(`_spmatMult_spsp_prodSp_openmp`, A, B, isProdSym, p, nthreads)
}

spsp_prodSp_fast_openmp <- function(A, B, isProdSym, nthreads) {
    .Call(`_spmatMult_spsp_prodSp_fast_openmp`, A, B, isProdSym, nthreads)
}

spsp_prodMat_openmp <- function(A, B, isProdSym, nthreads) {
    .Call(`_spmatMult_spsp_prodMat_openmp`, A, B, isProdSym, nthreads)
}

