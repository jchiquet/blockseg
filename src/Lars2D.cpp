#include "RcppArmadillo.h"
#include "Lars2D.h"

Lars2D::Lars2D(const mat &Y,
               const uword maxBreaks,
               const uword maxVar,
               const bool verbose,
               const bool Beta) :
    maxBreaks(maxBreaks) , // max number of change point per dimension
    maxVar   (maxVar)    , // max number of activated variables
    verbose  (verbose)   , // verbose mode 
    Beta     (Beta)        // should Beta be sent back?
{
  
  n   = Y.n_cols   ; // problem dimension (matrix form)
  Xtr = XtProd(Y)  ; // current correlation 
  m   = Xtr.n_elem ; // problem dimension (once vectorized)
  
  // ACTIVE-SET
  A = activeSet(m) ;
  
  // CHOLESKY FACTORIZATION
  CholFact = mat()  ;
}

void Lars2D::run() {
    
  // "THE" MAIN LOOP
  uword k = 0     ; // current iterations step
  vec  BetaValues ; // set of vectors to handle the sparse encoding of
  uvec BetaRowInd ; // the solution path
  uvec BetaColInd ; 
  uvec RowBreaks  ;
  uvec ColBreaks  ;
  bool drop = false ;

  while(A.getnActive() < maxVar && ( RowBreaks.n_elem < maxBreaks || ColBreaks.n_elem < maxBreaks)) {
      ivec actions   ;  // reset the current actions
      k++ ;

    // Identify the largest nonactive gradient (also the next lambda)
    double lambda = norm(Xtr.elem(A.getDormant()), "inf") ;

    // --------------------------------------------------------
    // S1. ADDING NEW VARIABLE(S)
    // Find which variable(s) to let in among dormant (inactive) variables
    uvec dormant = A.getDormant() ;
    uvec toAdd = dormant.elem(find(abs(Xtr.elem(dormant))>=lambda-ZERO)) ;

    // Update active set and Cholesky factorization
    for (uword i=0; i != (toAdd.n_elem); i++) {
      A.addition(toAdd[i]) ;
      CholeskyInsertLastActive() ;
      if(verbose) {
          Rprintf("Step %u:\t variable %u \tadded\n", k, toAdd[i]) ;
      }
      // update current actions
      actions.resize(actions.n_elem+1)    ; 
      actions[actions.n_elem -1] = toAdd[i]+1 ;
    }
    // S1. END ------------------------------------------------

    // --------------------------------------------------------
    // S2. COMPUTE THE DIRECTION OF DESCENT
    // 
    // wA_tilde (resp. wa) is the unormalized (resp. normalized) direction
    vec sA  = sign(Xtr.elem(A.getActive())) ; // sA are the sign of the active correlation
    vec wA_tilde = solve(trimatu(CholFact), solve( trimatl(strans(CholFact)), sA)) ; 
    double B  = 1/sqrt(sum(wA_tilde % sA)) ; vec wA  = wA_tilde * B ; 
    // X wA    gives equiangulor vectors
    // X'XA wA gives the corresponding  direction for y.hat,
    vec alpha  = XtXAProd(wA) ; 
    // S2. END ------------------------------------------------

    // --------------------------------------------------------
    // S3. COMPUTE THE DIRECTION STEP 
    // 
    // How far we go  before the next competitor arrives
    double gamhat = lambda/B ;
    if (A.getnActive() < m) {
      for (uword i = 0; i < m; i++) {
	if (A.isActive(i))
	  continue;
	double gam1 = (lambda - Xtr[i]) / (B - alpha[i]);
	double gam2 = (lambda + Xtr[i]) / (B + alpha[i]);
	if ((gam1 > ZERO) && (gam1 < gamhat))
            gamhat = gam1;
	if ((gam2 > ZERO) && (gam2 < gamhat))
            gamhat = gam2;
      }
    }
    // S3. END ------------------------------------------------
    // --------------------------------------------------------
    // S4. DROP VARIABLES CROSSING THE ZERO LINE
    double gam0 = 2*fabs(gamhat);
    uvec indToDrop ;
    for (uword i = 0; i < A.getnActive(); i++) {
      double val = -A.getValues()[i] / wA[i];
      if ((val > 0) && (val < gam0)) {
	gam0 = val;
        indToDrop.resize(indToDrop.n_elem+1)    ; 
        indToDrop[indToDrop.n_elem -1] = i ;
      }
    }
    if (gam0 < gamhat) {
        gamhat = gam0;      
	drop = true ;
    } else {
        drop = false ;
    }
    
    // Moving beta with the right amount of gamma
    A.setValues(A.getValues() + gamhat * wA) ;
    // Updating the gradient according to the updated residuals
    Xtr = Xtr - gamhat * alpha ;
        
    if (drop) {
        indToDrop = sort(indToDrop,"descend");
        for(uword i=0; i < indToDrop.n_elem; i++) {
            CholeskyDeleteIndActive(indToDrop[i]) ;
            if(verbose) {
                Rprintf("Step %u:\t variable %u \tdeleted\n", k, A.getActive()[indToDrop[i]]) ;
            } 
            // update current actions
            actions.resize(actions.n_elem+1)    ; 
            actions[actions.n_elem -1] = -(A.getActive()[indToDrop[i]]+1) ;
            A.deletion(indToDrop[i]) ;
        }
    }
    // S4. END ------------------------------------------------

    // Record variables related to the regularization path
    BetaValues = join_cols(BetaValues, A.getValues()) ;
    BetaRowInd = join_cols(BetaRowInd, conv_to<uvec>::from((k-1)*ones(A.getnActive()))) ;
    BetaColInd = join_cols(BetaColInd, A.getActive()) ;
    uvec div = floor((A.getActive())/n);
    RowBreaks = 1 + A.getActive() - div*n ;
    ColBreaks = 1+ div ;
    BreaksRow.push_back(RowBreaks) ;
    BreaksCol.push_back(ColBreaks) ;
    LambdaPath.push_back(lambda)   ;
    ActionsPath.push_back(actions) ;
    RowBreaks = unique(RowBreaks)  ; // not just keep the unique breakpoints
    ColBreaks = unique(ColBreaks)  ;

  } // END MAIN WHILE LOOP  

  if (Beta) { // if requested, then back the matrices of Beta values
      umat loc= join_rows(BetaRowInd,BetaColInd).t();
      BetaPath = sp_mat(loc, BetaValues, k, m) ;
  }
}

// Private functions.

// METHODS FOR CHOLESKY UPDATES
void Lars2D::CholeskyInsertLastActive() {
  
  // Building the currently added vector xtxAj
  uvec div = floor((A.getActive())/n);
  uvec qA  = 1 + div ;
  uvec rA  = 1 + A.getActive() - div*n ;
  for (uword i=0;i<A.getnActive();i++) {
    if (qA(A.getnActive()-1)>qA(i)) {
      qA(i) = qA(A.getnActive()-1) ;
    } 
    if (rA(A.getnActive()-1)>rA(i)) {
      rA(i) = rA(A.getnActive()-1) ;
    }
  }
  
  vec xtx_j = conv_to<vec>::from(n-qA+1) % conv_to<vec>::from(n-rA+1) ;
  // Updating the Cholesky Factorization 
  if (A.getnActive() == 1) {
    CholFact = sqrt(as_scalar(xtx_j));
  } else {
    colvec rp  = arma::zeros<arma::colvec>(A.getnActive(),1);
    rp.subvec(0,A.getnActive()-2) = solve (trimatl(CholFact.t()), xtx_j.subvec(0,A.getnActive()-2));
    rp(A.getnActive()-1) = sqrt(xtx_j(A.getnActive()-1)-dot(rp,rp));
    CholFact = join_rows( join_cols(CholFact, zeros<mat>(1,A.getnActive()-1)), rp);
  }
  
}

void Lars2D::CholeskyDeleteIndActive(uword j) {
  
  vec x = zeros<vec>(2,1); // currently removed column
  mat R = zeros<mat>(2,2); // rotation matrix
  
  CholFact.shed_col(j);
  int nActive = CholFact.n_cols;
  double r;
  for (int k=j; k<nActive; k++) {
    x = CholFact.submat(k,k,k+1,k);
    if (x[1] != 0) {
      r = norm(x,2);
      R(0,0) = x(0) ; R(0,1) = x(1) ; R(1,0) = -x(1) ; R(1,1) = x(0) ; 
      R = R / r;
      x(0) = r; x(1) = 0;
    } else {
      R = arma::eye(2,2);
    }
    CholFact.submat(k,k,k+1,k) = x;
    if (k + 1 < nActive) {
      CholFact.submat(k,k+1,k+1,nActive-1) = R * CholFact.submat(k,k+1,k+1,nActive-1);
    }
  }
  CholFact.shed_row(nActive);
}

arma::vec Lars2D::XtProd(const arma::mat& X) {
  return(vectorise(flipud(fliplr(cumsum(cumsum(flipud(fliplr(X)),0),1)))));
}

arma::vec Lars2D::XProd(const arma::mat& X) {
    return(vectorise(cumsum(arma::cumsum(X,0),1))) ;
}

arma::vec Lars2D::XtXAProd(const arma::vec& wA) {
  vec w = zeros(m) ;
  w.elem(A.getActive()) = wA ;
  mat X = reshape(w, n, n) ;
  
  return(vectorise(flipud(fliplr(cumsum(cumsum(flipud(fliplr(cumsum(cumsum(X,0),1))),0),1))))) ;
} 
