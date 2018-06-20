// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "Lars2D.h"

using namespace arma;
using namespace Rcpp;

#if !defined(ARMA_DEFAULT_OSTREAM)
#define ARMA_DEFAULT_OSTREAM Rcpp::Rcout
#endif

// [[Rcpp::export]]
List doLARS2D(SEXP R_Y, SEXP R_maxBreaks, SEXP R_maxVar, SEXP R_verbose, SEXP R_Beta){
  // disable messages being printed to the err2 stream
  std::ostream nullstream(0);
  set_stream_err2(nullstream);

  Rcpp::traits::input_parameter< const arma::mat& >::type Y(R_Y);
  Rcpp::traits::input_parameter< const uword& >::type maxBreaks(R_maxBreaks);
  Rcpp::traits::input_parameter< const uword& >::type maxVar(R_maxVar);
  Rcpp::traits::input_parameter< const bool& >::type verbose(R_verbose);
  Rcpp::traits::input_parameter< const bool& >::type Beta(R_Beta);

  // Initializing the Lars object
  Lars2D seg2D(Y, maxBreaks, maxVar, verbose, Beta);

  // Run for something
  seg2D.run() ;
  return List::create(Named("Beta")      = wrap(seg2D.getBeta())     ,
		      Named("Lambda")    = wrap(seg2D.getLambda())   ,
		      Named("RowBreaks") = wrap(seg2D.getRowBreaks()),
		      Named("ColBreaks") = wrap(seg2D.getColBreaks()),
		      Named("Actions")   = wrap(seg2D.getActions()) );
}

// [[Rcpp::export]]
NumericVector Xprod(SEXP R_Y){
  Rcpp::traits::input_parameter< const arma::mat& >::type Y(R_Y);

  Lars2D seg2D(Y, 1, 1, false, false);
  
  return(wrap(seg2D.XProd(Y)));
}

