

#include <RcppArmadillo.h>
#include<omp.h>



using namespace Rcpp ;
using namespace std;
using namespace arma;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp,cpp11)]]

//// Additional function to sort row indices by columns in a CSC matrix
sp_mat sort_row_indices(uvec &ri, uvec &cptr, vec &val, uword nrow, uword ncols){
  
  // No. of non-zero vals per column from the col-ptr
  uvec nc = diff(cptr);
  field<uvec> ris(ncols);
  uvec rr;
  
  for(int it=0; it< ncols; it++) { 
    uword ncl = nc(it);
    if(ncl > 0){ 
      uword nvals_uptoPrev= cptr[it];
      uvec ricol= ri.subvec(nvals_uptoPrev,nvals_uptoPrev+ncl-1);
      bool checkr = ricol.is_sorted();
      if(checkr==0){ 
        rr = sort_index(ricol);
        ri.subvec(nvals_uptoPrev,nvals_uptoPrev+ncl-1) = ricol.elem(rr);
        vec valcol = val.subvec(nvals_uptoPrev,nvals_uptoPrev+ncl-1);
        val.subvec(nvals_uptoPrev,nvals_uptoPrev+ncl-1) = valcol.elem(rr);
      }
    }
  } 
  sp_mat X(ri,cptr,val,nrow,ncols);
  
  return X;
}



//// add an element every element of an unsigned vector
uvec element_add(uvec &A,int &x){

  for(int i=0;i<A.n_elem;i++){
    A[i]+=x;
  }

  return A;
}




/// column bind multiple sparse matrices
// [[Rcpp::export]]
sp_mat cbind_field_spmat(field<sp_mat> A){

  int n = A.n_rows;
  int nmax = 0;
  int ncols = 0;
  int nrow = A(0,0).n_rows;


  for(int i=0;i<n;i++){
    nmax += A(i,0).n_nonzero;
    ncols += A(i,0).n_cols;
  }

  uvec ir(nmax);
  vec vals(nmax);
  uvec cptr(ncols+1,fill::zeros);
  int idx = 0;
  int idc =1;

  for(int i=0;i<n;i++){
    int nval = A(i,0).n_nonzero;
    int nc = A(i,0).n_cols; 
    
     if(nval>0){
       uvec Air(A(i,0).row_indices,nval);
       ir.subvec(idx,idx+nval-1)= Air ;

       vec Aval(A(i,0).values,nval);
       vals.subvec(idx,idx+nval-1)= Aval;
      
     }
      uvec cptrA(A(i,0).col_ptrs,nc+1);
      cptrA= element_add(cptrA,idx);
      cptr.subvec(idc,idc+nc-1) =cptrA.subvec(1,nc);
      idx+=nval;
      idc +=nc;
  }

  // Put in the sp_mat container
  sp_mat X(ir,cptr,vals,nrow,ncols);

  return X;

}



// [[Rcpp::export]]
uvec resizethis(uvec &A){
  int len=A.n_elem;
  A.resize(3*len);
  return A;
}








////////////////////////////////////////////////////////////
////////////////// SPARSE MATRIX MULTIPLICATION SUITE (V2) *** 
///////////////////////////////////////////////////////////


/////////////////////////// SP-SERIAL


// SERIAL -- This function evaluates A * B where both A & B are sparse and the resultant
// product is also sparse 


// [[Rcpp::export]]
sp_mat spsp_to_spgen_serial(const arma::sp_mat &A, const arma::sp_mat &B, double p){
  
   
  // Define matrix sizes
  const uword mA= A.n_rows;
  const uword nB= B.n_cols;
  const uword nnzA = A.n_nonzero;
  const uword nnzB=B.n_nonzero;
  sp_mat C(mA,nB);
  
  if(nnzA>0 && nnzB>0){
    
    // approximate initial number of non-zeros in the resultant matrix
     uword nnzC = ceil(mA * nB * p);
  
   // Initialize colptr, row_index and value vectors for the resultant sparse matrix
     uvec colptrC(nB+1);
     colptrC.zeros();
  
  
    uvec rowvalC(nnzC);
    rowvalC.zeros();
  
  
    colvec nzvalC(nnzC);
    nzvalC.zeros();
  
  
   //setenv("OMP_STACKSIZE","500M",1);
  
   // Counters and other scratch variables
    uword i, jp, j, kp, k; 
    uword nzmax = 0;
    double nzB, nzA; 
    ivec xb(mA);
    xb.fill(-1);
    const int increasefactor =2;
  
  
   //uvec nzmax_vec(nB,fill::zeros);
  
   // Loop Logic :: CSC format requires traversing by entries of column to optimize efficiency // 
  
   //// Outer loop is over columns of B
  
   ////// For each column of B, mid loop is over non-zero elements of that columns of B
  
   ///////// For each nnz. element of that column of B, inner loop is over nnz. elements of corresponding
   ///////// column (that matches with row number of that element of B) of A 
  
   ////// So, as we traverse down elements of a column of B, we travers right over corr. columns of A
   ///// & then sum accordingly
  
   for(i=0; i< nB; i++) { 
    
    
     for ( jp = B.col_ptrs[i]; jp < B.col_ptrs[i+1]; jp++) {
      
       j = B.row_indices[jp];
       nzB = B.values[jp];
      
      
       for ( kp = A.col_ptrs[j]; kp < A.col_ptrs[j+1]; kp++ ){
        
         k = A.row_indices[kp];
         nzA = A.values[kp];
        
        
         if (xb(k) != -1){
           nzvalC(xb(k))+= nzA * nzB;
         } else {
           rowvalC(nzmax) = k;
           xb(k) = nzmax;
           nzvalC(nzmax) = nzA * nzB;
           nzmax +=1;
          
         }
       }
     }
    
    
     // fill it back with -1s (to save time: fill only where it changed, rest are -1 anyway) 
     for(k = colptrC(i); k < nzmax; ++k){
       xb(rowvalC(k)) = -1;
     }
    
    colptrC(i+1) = nzmax;
    
     //* Check if more memory is required
     if ( (nnzC - nzmax)< mA) {
       rowvalC.resize(increasefactor*nnzC);
       nzvalC.resize(increasefactor*nnzC);
       nnzC *= increasefactor;
     } 
   }
  
   // Put in the sp_mat container: it is already ordered! 
  uvec rit = rowvalC.subvec(0,nzmax-1);
  vec nzvalt = nzvalC.subvec(0,nzmax-1);
  C = sort_row_indices(rit,colptrC,nzvalt,mA,nB);
  
  } 
  
  return C;
}  




// This function evaluates A * B where both A & B are sparse and the resultant
// product is also sparse & symmetric (compute lower traingular only)

// [[Rcpp::export]]
sp_mat spsp_to_spsym_serial(const arma::sp_mat &A, const arma::sp_mat &B, double p){

  // Define matrix sizes
  const uword mA= A.n_rows;
  const uword nB= B.n_cols;
  const uword nnzA = A.n_nonzero;
  const uword nnzB=B.n_nonzero;
  sp_mat C(mA,nB);
  
  if(nnzA>0 && nnzB>0){
    

  // approximate initial number of non-zeros in the resultant matrix
  uword nnzC = ceil(mA * nB * p);

  // Initialize colptr, row_index and value vectors for the resultant sparse matrix
  uvec colptrC(nB+1);
  colptrC.zeros();

  uvec rowvalC(nnzC);
  rowvalC.zeros();

  colvec nzvalC(nnzC);
  nzvalC.zeros();

  //setenv("OMP_STACKSIZE","500M",1);

  // Counters and other scratch variables
  uword i, jp, j, kp, k;
  uword nzmax = 0;
  double nzB, nzA;
  ivec xb(mA);
  xb.fill(-1);
  const int increasefactor =2;


  // Loop Logic :: CSC format requires traversing by entries of column to optimize efficiency //

  //// Outer loop is over columns of B

  ////// For each column of B, mid loop is over non-zero elements of that columns of B

  ///////// For each nnz. element of that column of B, inner loop is over nnz. elements of corresponding
  ///////// column (that matches with row number of that element of B) of A

  ////// So, as we traverse down elements of a column of B, we travers right over corr. columns of A
  ///// & then sum accordingly

  for(i=0; i< nB; i++) {

    for ( jp = B.col_ptrs[i]; jp < B.col_ptrs[i+1]; jp++) {

      j = B.row_indices[jp];
      nzB = B.values[jp];

      for ( kp = A.col_ptrs[j]; kp < A.col_ptrs[j+1]; kp++ ){

        k = A.row_indices[kp];

        if(i<=k) {

          nzA = A.values[kp];

          if (xb(k) != -1){
            nzvalC(xb(k))+= nzA * nzB;
          } else {
            rowvalC(nzmax) = k;
            xb(k) = nzmax;
            nzvalC(nzmax) = nzA * nzB;
            nzmax +=1;
          }
        }
      }
    }


    // fill it back with -1s (to save time: fill only where it changed, rest are -1 anyway)
    for(k = colptrC(i); k < nzmax; ++k){
      xb(rowvalC(k)) = -1;
    }

    colptrC(i+1) = nzmax;

    //* Check if more memory is required
    if ( (nnzC - nzmax)< mA) {
      rowvalC.resize(increasefactor*nnzC);
      nzvalC.resize(increasefactor*nnzC);
      nnzC *= increasefactor;
    }

  }

  // Put in the sp_mat container: it is already ordered!
   uvec rit = rowvalC.subvec(0,nzmax-1);
   vec nzvalt = nzvalC.subvec(0,nzmax-1);
   C = sort_row_indices(rit,colptrC,nzvalt,mA,nB);
  }

  return C;
}

//////////////////////////////////////////

///////////////////////////// SP - OMP

// openmp version of the spsp_to_spgen_serial. It in fact uses spsp_to_spgen_serial()

// [[Rcpp::export]]
sp_mat spsp_to_spgen_openmp(const arma::sp_mat &A, const arma::sp_mat &B, double p,
                            int nthreads){

  
  int nchunks = nthreads;
  int ncolB= B.n_cols;
  int ncols_perChunk= ncolB/nchunks;
  int d = ncols_perChunk -1;
  sp_mat C;
  field<sp_mat> Cf(nchunks);
  int i;

  if(ncolB<nchunks){
    C = spsp_to_spgen_serial(A,B,p);
  } else {

  // set stacksize, not really necessary unless you are running into memory issues
    setenv("OMP_STACKSIZE","100M",1);
    
  omp_set_num_threads(nthreads);
 #pragma omp parallel for shared(Cf,A,B,d,p,ncols_perChunk,ncolB,nchunks) private(i) default(none) schedule(auto)
    for(i=0;i<nchunks;i++){
      //Rcout << "i= "<< i << std::endl;
      if(i < (nchunks-1)) {
        sp_mat Bchunk=B.cols(ncols_perChunk*i,ncols_perChunk*i+d);
        Cf(i,0)= spsp_to_spgen_serial(A,Bchunk,p);
      } else {
        sp_mat Bchunk=B.cols(ncols_perChunk*i,ncolB-1);
        Cf(i,0)= spsp_to_spgen_serial(A,Bchunk,p);
      }
    }

    C= cbind_field_spmat(Cf);
  }

  return C;
}



// openmp version of the spsp_to_spgen_serial. It in fact uses spsp_to_spgen_serial()


sp_mat spsp_to_spsym_openmp(const arma::sp_mat &A, const arma::sp_mat &B, double p,
                            int nthreads){

  omp_set_num_threads(nthreads);
  int nchunks = nthreads;
  int ncolB= B.n_cols;
  int ncols_perChunk= ncolB/nchunks;
  int d = ncols_perChunk -1;
  sp_mat C;
  field<sp_mat> Cf(nchunks);
  int i;

  if(ncolB<nchunks){
    C = spsp_to_spsym_serial(A,B,p);
  } else {

    // set stacksize, not really necessary unless you are running into memory issues
    setenv("OMP_STACKSIZE","100M",1);
    
    omp_set_num_threads(nthreads);
#pragma omp parallel for shared(Cf,A,B,d,p,ncols_perChunk,ncolB,nchunks) private(i) default(none) schedule(auto)
    for(i=0;i<nchunks;i++){
      if(i < (nchunks-1)) {
        Cf(i,0)= spsp_to_spsym_serial(A,B.cols(ncols_perChunk*i,ncols_perChunk*i+d),p);
      } else {
        Cf(i,0)= spsp_to_spsym_serial(A,B.cols(ncols_perChunk*i,ncolB-1),p);
      }
    }

    C= cbind_field_spmat(Cf);
  }

  return C;
}

///////////////////////////////////////////





//////////////////////////// SP - FAST (SMALL MATRIX)

// This function evaluates A * B where both A & B are sparse and the resultant
// product is also sparse
// It's basically a hack! = we first compute it as dense but then use Armadillo
// sp_mat() container to map it back to sparse!!
/// -- THIS IS FASTEST BUT DON'T USE IT WHERE RESULTING MATRIX IS LARGE**

sp_mat spsp_to_spgen_fast_openmp(const arma::sp_mat &A, const arma::sp_mat &B, int nthreads){

  // define matrix sizes
  const uword mA= A.n_rows;
  const uword nB= B.n_cols;
  const uword nnzA = A.n_nonzero;
  const uword nnzB=B.n_nonzero;
  sp_mat Cab(mA,nB);

  if(nnzA>0 && nnzB>0){


  // initialize the resultant dense matrix
  mat C(mA,nB,fill::zeros);
  double* Cptr = C.begin();


  // counters and other variables
  uword i, jp, j, kp, k;

  omp_set_num_threads(nthreads);
#pragma omp parallel for shared(Cptr,A,B) private(i,j,k,jp,kp) default(none) schedule(static)
  for(i=0; i< nB; i++) {

    uword nelem_bfr_ith_col = i*mA;

    for ( jp = B.col_ptrs[i]; jp < B.col_ptrs[i+1]; jp++) {

      j = B.row_indices[jp];
      double nzB = B.values[jp];

      for ( kp = A.col_ptrs[j]; kp < A.col_ptrs[j+1]; kp++ ){

        k = A.row_indices[kp];
        double nzA = A.values[kp];
        Cptr[nelem_bfr_ith_col+k] += nzA * nzB;

       }
     }
   }
   Cab=conv_to<arma::sp_mat>::from(C);
 }

  return Cab;
}



// sparse * sparse = sparse in fast method
sp_mat spsp_to_spsym_fast_openmp(const arma::sp_mat &A, const arma::sp_mat &B,int nthreads){

  // define matrix sizes
  const uword mA= A.n_rows;
  const uword nB= B.n_cols;
  const uword nnzA = A.n_nonzero;
  const uword nnzB=B.n_nonzero;
  sp_mat Cab(mA,nB);

  if(nnzA>0 && nnzB>0){


  // initialize the resultant dense matrix
  mat C(mA,nB,fill::zeros);
  double* Cptr = C.begin();


  // counters and other variables
  uword i, jp, j, kp, k;

  omp_set_num_threads(nthreads);
#pragma omp parallel for shared(Cptr,A,B) private(i,j,k,jp,kp) default(none) schedule(static)
  for(i=0; i< nB; i++) {

    uword nelem_bfr_ith_col = i*mA;

    for ( jp = B.col_ptrs[i]; jp < B.col_ptrs[i+1]; jp++) {

      j = B.row_indices[jp];

      if(i<=k)  {

        double nzB = B.values[jp];
        for ( kp = A.col_ptrs[j]; kp < A.col_ptrs[j+1]; kp++ ){
          k = A.row_indices[kp];
          double nzA = A.values[kp];
          Cptr[nelem_bfr_ith_col+k] += nzA * nzB;
         }
       }
     }
   }
   Cab = conv_to<arma::sp_mat>::from(C);
  }

  return Cab;
}

//////////////////////////////////////



///////////////////////////// SP - DENSE
// This function evaluates A * B where both A & B are sparse and the resultant
// product is dense

mat spsp_to_dnsgen_openmp(const arma::sp_mat &A, const arma::sp_mat &B, int nthreads){

  // define matrix sizes
  const uword mA= A.n_rows;
  const uword nB= B.n_cols;
  const uword nnzA = A.n_nonzero;
  const uword nnzB=B.n_nonzero;

  // initialize the resultant dense matrix
  mat C(mA,nB,fill::zeros);


  if(nnzA>0 && nnzB>0){

  // counters and other variables
  double* Cptr = C.begin();
  uword i, jp, j, kp, k;

  omp_set_num_threads(nthreads);
#pragma omp parallel for shared(Cptr,A,B) private(i,j,k,jp,kp) default(none) schedule(static)
  for(i=0; i< nB; i++) {

    uword nelem_bfr_ith_col = i*mA;

    for ( jp = B.col_ptrs[i]; jp < B.col_ptrs[i+1]; jp++) {

      j = B.row_indices[jp];
      double nzB = B.values[jp];

      for ( kp = A.col_ptrs[j]; kp < A.col_ptrs[j+1]; kp++ ){

        k = A.row_indices[kp];
        double nzA = A.values[kp];
        Cptr[nelem_bfr_ith_col+k] += nzA * nzB;

       }
     }
    }
  }
  return C;
}



// This function evaluates A * B where both A & B are sparse and the resultant
// product is dense-symmetric


mat spsp_to_dnssym_openmp(const arma::sp_mat &A, const arma::sp_mat &B, int nthreads){

  // define matrix sizes
  const uword mA= A.n_rows;
  const uword nB= B.n_cols;
  const uword nnzA = A.n_nonzero;
  const uword nnzB=B.n_nonzero;

  // initialize the resultant dense matrix
  mat C(mA,nB,fill::zeros);


  if(nnzA>0 && nnzB>0){

    // counters and other variables
    double* Cptr = C.begin();
    uword i, jp, j, kp, k;

  omp_set_num_threads(nthreads);
#pragma omp parallel for shared(Cptr,A,B) private(i,j,k,jp,kp) default(none) schedule(static)
  for(i=0; i< nB; i++) {

    uword nelem_bfr_ith_col = i*mA;

    for ( jp = B.col_ptrs[i]; jp < B.col_ptrs[i+1]; jp++) {

      j = B.row_indices[jp];
      double nzB = B.values[jp];

      for ( kp = A.col_ptrs[j]; kp < A.col_ptrs[j+1]; kp++ ){

        k = A.row_indices[kp];

        if(i <= k){
          double nzA = A.values[kp];
          Cptr[nelem_bfr_ith_col+k] += nzA * nzB;
         }
       }
     }
   }
  }
  return C;
}

///////////////////////////////////////




//////////////////////////// SP - SCALAR
// This function evaluates A * B where both A & B are sparse but either A or B is a scalar

sp_mat spspscalar_to_sp(const arma::sp_mat &A, const arma::sp_mat &B){

  // define matrix sizes
  const uword mA= A.n_rows;
  const uword nA = A.n_cols;
  const uword mB = B.n_rows;
  const uword nB= B.n_cols;


  if(mA==1 & nA==1){

    return A.values[0] * B ;

  } else if (mB==1 & nB==1){

    return B.values[0] * A;

  }


}
//////////////////////////////


//////////////////////////////////////////////////////////////////
// This function is the MAIN function that calls the above
//////////////////////////////////////////////////////////////////



// [[Rcpp::export]]
sp_mat spsp_prodSp_serial(const sp_mat &A, const sp_mat &B, const uword &isProdSym, double p){

  // choose from prior expert knowledge a non-zero fraction for AB
  //double p =.01;


  // define matrix sizes
  uvec dims(4);

  dims.at(0) = A.n_rows;
  dims.at(1) = A.n_cols;
  dims.at(2) = B.n_rows;
  dims.at(3) = B.n_cols;
  int dA= dims(0)*dims(1);
  int dB=dims(2)*dims(3);

  if( dA ==1 | dB ==1){

    return spspscalar_to_sp(A,B);

  } else if(isProdSym==0){

    return spsp_to_spgen_serial(A,B,p);

  } else if(isProdSym==1){

    return spsp_to_spsym_serial(A,B,p);

  }

}



// [[Rcpp::export]]
sp_mat spsp_prodSp_openmp(const sp_mat &A, const sp_mat &B, const uword &isProdSym, double p, int nthreadsUser){
  
  //// MATRIX SIZES
  int ncA=A.n_cols;
  int nrB=B.n_rows;
  int nrAB=A.n_rows;
  int ncAB=B.n_cols;
  int dA= nrAB*ncA;
  int dB=nrB*ncAB;
  //int dAB= nrAB*ncAB;
  ////
  
  
  
  
  ///// THREADS MANAGEMENT
  // Obtain environment containing function
  Rcpp::Environment base("package:parallel"); 
  
  // Make function callable from C++
  Rcpp::Function detectCores = base["detectCores"];
  bool logical=FALSE;
  bool alltest=FALSE;
  int ncores= Rcpp::as<int >(detectCores(alltest,logical));
  
  //int nthreadsMax = omp_get_max_threads();
  int nthreads=1;
  if(nthreadsUser>ncores){
    nthreads = ncores;
  }else if(nthreadsUser==ncores){
    nthreads = ncores-1;
  } else if(nthreadsUser<0){
    nthreads = 1;
  } else {
    nthreads=nthreadsUser;
  }
  //////
  
  
  //// SPARSITY DENSITY
  int nnzA=A.n_nonzero;
  int nnzB=B.n_nonzero;
  double pA = nnzA/(double)dA;
  double pB = nnzB/(double)dB;
  double pAB=.1;
  //double pAB = std::max(p,std::min(pA+pB,1.0));
  
  // choose from prior expert knowledge a non-zero fraction for AB
  if(pA<.05 && pB<.05){
    pAB = pA + pB;
  }
  
  if(p > .1){
    pAB=p;
  }
  
  // check if p*#nrowA > nrowA if not set it s.t. it covers at least one col
  double ncABperCore=ncAB/nthreads;
  int nrpAB = ceil(p*nrAB*ncABperCore);
  if(nrpAB<=nrAB){
    pAB = std::max(.1,2*((nrAB)/(double)(nrAB*ncABperCore))); 
  }
  /////
  
  
  
  
  ////  EMPTY CASE
  if(nnzA==0 || nnzB==0) {
    sp_mat AB(nrAB,ncAB);
    return AB;
  }


  //// NON-EMPTY CASE
  if( dA ==1 || dB ==1){

    return spspscalar_to_sp(A,B);

  } else if(isProdSym==0){

    return spsp_to_spgen_openmp(A,B,pAB,nthreads);

  } else if(isProdSym==1){

    return spsp_to_spsym_openmp(A,B,pAB,nthreads);

  }
  //
  
  
}




// [[Rcpp::export]]
sp_mat spsp_prodSp_fast_openmp(const sp_mat &A, const sp_mat &B, const uword &isProdSym, int nthreadsUser){
  
  
  // define matrix sizes
  uvec dims(4);
  
  dims.at(0) = A.n_rows;
  dims.at(1) = A.n_cols;
  dims.at(2) = B.n_rows;
  dims.at(3) = B.n_cols;
  int dA= dims(0)*dims(1);
  int dB=dims(2)*dims(3);
  
  
    ///// THREADS MANAGEMENT
  // Obtain environment containing function
  Rcpp::Environment base("package:parallel"); 
  
  // Make function callable from C++
  Rcpp::Function detectCores = base["detectCores"];
  bool logical=FALSE;
  bool alltest=FALSE;
  int ncores= Rcpp::as<int >(detectCores(alltest,logical));
  
  //int nthreadsMax = omp_get_max_threads();
  int nthreads=1;
  if(nthreadsUser>ncores){
    nthreads = ncores;
  }else if(nthreadsUser==ncores){
    nthreads = ncores-1;
  } else if(nthreadsUser<0){
    nthreads = 1;
  } else {
    nthreads=nthreadsUser;
  }
  //////
   
  
  
  
  if( dA ==1 | dB ==1){
    
    return spspscalar_to_sp(A,B);
    
  } else if(isProdSym==0){
    
    return spsp_to_spgen_fast_openmp(A,B,nthreads);
    
  } else if(isProdSym==1){
    
    return spsp_to_spsym_fast_openmp(A,B,nthreads);
    
  }
  
}


// [[Rcpp::export]]
mat spsp_prodMat_openmp(const sp_mat &A, const sp_mat &B,const uword &isProdSym, int nthreadsUser){
  
  
  ///// THREADS MANAGEMENT
  // Obtain environment containing function
  Rcpp::Environment base("package:parallel"); 
  
  // Make function callable from C++
  Rcpp::Function detectCores = base["detectCores"];
  bool logical=FALSE;
  bool alltest=FALSE;
  int ncores= Rcpp::as<int >(detectCores(alltest,logical));
  
  //int nthreadsMax = omp_get_max_threads();
  int nthreads=1;
  if(nthreadsUser>ncores){
    nthreads = ncores;
  }else if(nthreadsUser==ncores){
    nthreads = ncores-1;
  } else if(nthreadsUser<0){
    nthreads = 1;
  } else {
    nthreads=nthreadsUser;
  }
  //////
  

  
  
  if(isProdSym==1) {
    
    return spsp_to_dnssym_openmp(A,B,nthreads);
    
  } else {
    
    return spsp_to_dnsgen_openmp(A,B,nthreads);
    
  }
  
}

////////////////////////////////////////////////
////////////////////////////////////////////////

