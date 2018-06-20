/*
 * Author: Julien CHIQUET
 *         Statistique et Génome
 */
#ifndef _ACTIVE_SET_H
#define _ACTIVE_SET_H

#define ARMA_NO_DEBUG
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS

#include <RcppArmadillo.h>

#ifndef ARMA_HAVE_GETTIMEOFDAY
#define ARMA_HAVE_GETTIMEOFDAY
#endif

#ifndef ARMA_DEFAULT_OSTREAM
#define ARMA_DEFAULT_OSTREAM Rcpp::Rcout
#endif


#define ZERO 2e-16 // practical zero

using namespace arma;

class activeSet {

public:
  // CONSTRUCTORS
  activeSet();
  activeSet(const uword size);
  
  // ACCESSORS
  const uword & getnActive( ) const { return nActive   ; }
  const uvec  & getActive ( ) const { return whoActive ; }
  uvec getDormant( ) const { return find(amIActive < 1) ; }
  const vec   & getValues ( ) const { return valActive ; }
  void setValues (vec newValues) { valActive = newValues ; }
  bool isActive(uword i) { return (amIActive(i) == 1) ; }

  // METHODS FOR ACTIVE SET UPDATES
  void deletion(uword indVarOut);
  void addition(uword varIn    );

private:
  // VARIABLES FOR HANDLING THE ACTIVE SET //
  uword nActive   ; // number of active variables
  vec valActive   ; // current values of active variables
  uvec whoActive  ; // indices of active variables
  uvec whoDormant ; // indices of dormant (inactive) variables
  uvec amIActive  ; // activity status
};

#endif
