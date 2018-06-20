#include "activeSet.h"

// ______________________________________________________
// ACTIVE SET CLASS
//
// handle vector of active variable, with method for adding/removing
// elements

activeSet::activeSet() {}

activeSet::activeSet(const uword size) {
  // No active variable at start-up
  valActive.zeros(size) ;
  amIActive.zeros(size) ;
  nActive    = 0        ;
  whoDormant = linspace<uvec>(0,size-1,size) ;
}

// PRIVATE METHODS FOR ACTIVE SET UPDATES
void activeSet::deletion(uword indVarOut) {
  valActive.shed_row(indVarOut)       ; // update the vector of active parameters
  amIActive[whoActive[indVarOut]] = 0 ; // update the active set
  whoActive.shed_row(indVarOut)       ;
  nActive--                           ; // update the number of active variable
  // actions.resize(actions.n_elem+1)    ; // update current actions
  // actions[actions.n_elem -1] = -whoActive[indVarOut] ;
}

void activeSet::addition(uword varIn) {
  valActive.resize(nActive+1) ; // update the vector of active parameters
  valActive(nActive) = 0.0    ;
  whoActive.resize(nActive+1) ; // update the active set
  whoActive[nActive] = varIn  ;
  nActive++                   ; // update the number of active variable
  amIActive[varIn] = 1        ;
}
