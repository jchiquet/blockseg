/***
 * @file Lars2D.h
 *
 */
#ifndef __seg2D_HPP
#define __seg2D_HPP

#define ARMA_NO_DEBUG
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS
#define ARMA_64BIT_WORD

#include "activeSet.h"

#ifndef ARMA_HAVE_GETTIMEOFDAY
#define ARMA_HAVE_GETTIMEOFDAY
#endif

#ifndef ARMA_DEFAULT_OSTREAM
#define ARMA_DEFAULT_OSTREAM Rcpp::Rcout
#endif


#define ZERO 2e-16 // practical zero

using namespace arma;

class Lars2D {
public: 
  Lars2D(const mat &Y,
	 const uword maxBreaks,
	 const uword maxVar   ,
	 const bool verbose = true,
	 const bool Beta    = false);
  
  // public method to launch the LARS algorithm
  void run() ; 

  //! Access the set of coefficients along the regularization path
  const arma::sp_mat        & getBeta      () const { return BetaPath   ; }
  const std::vector<double> & getLambda    () const { return LambdaPath ; }
  const std::vector<ivec>   & getActions   () const { return ActionsPath; }
  const std::vector<uvec>   & getRowBreaks () const { return BreaksRow  ; }
  const std::vector<uvec>   & getColBreaks () const { return BreaksCol  ; }
  arma::vec XProd(const arma::mat& X) ;


private:
  uword n         ; // Problem dimension (in matrix form)  
  uword m         ; // Problem dimension (once vectorized)  
  uword maxBreaks ; // maximal number of activated segment
  uword maxVar    ; // maximal number of activated segment
  bool verbose    ; // verbose mode status
  bool Beta       ; // should the Beta matrix be sent back?
  mat Xtr         ; // vector of current correlations
  activeSet A     ; // active set class
  mat CholFact    ; // Cholesky factorization (triangular matrix)

  // VARIABLES FOR  HANDLING THE SOLUTION PATH 
  // Solution path sparsely encoded
  arma::sp_mat BetaPath          ; // list of sparse matrices  
  std::vector<double> LambdaPath ; // vector of values of lambda along the path
  std::vector<ivec>  ActionsPath ; // vector keeping trace of what comes in and out
  std::vector<uvec>  BreaksRow   ; // vector keeping trace of row breaks 
  std::vector<uvec>  BreaksCol   ; // vector keeping trace of column breaks 

  void CholeskyInsertLastActive();
  void CholeskyDeleteIndActive(uword j);

  // METHODS FOR FAST COMPUTATION
  arma::vec XtProd(const arma::mat& X) ;
  arma::vec XtXAProd(const arma::vec& wA) ;

};

#endif

