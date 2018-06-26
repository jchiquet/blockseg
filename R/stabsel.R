##' \code{stab.blockSeg} algorithm
##'
##' Model selection for the \code{blockSeg} algorithm.
##'
##' @param Y matrix of observations.
##' @param nsimu a positive integer.
##' @param max.break a positive integer less than number of columns divided by 2 and number of rows divided by 2.
##' @param max.var a positive integer less than number of columns times number of rows.
##' By default, \code{ncol(Y)**2/8}.
##' @param random.break logical. To change the position of the first row (resp. column); the rows
##' before this position are moved to the end. By default TRUE.
##' @param sym.break logical. In the case of symmetric matrices, it is possible to accumulate breaks
##' in row and columns to improve the quality of the estimation. By default FALSE.
##' Warning: a check is made on the dimensions of the matrix but not on the fact that it is symmetrical
##' or not; this choice was made for the case where the user would like to have symmetrical breaks
##' even if the matrix is not (not recommended by the authors of the package).
##' @param mc.cores a positive integer giving the number of cores used. If you use windows, the parallelization is impossible.
##' By default, 2
##' @param verbose logical. To display each step. By default TRUE.
##'
##' @rdname stab.blockSeg
##'
##' @examples
##'  ## model parameters
##' n <- 100
##' K <- 5
##' mu <- suppressWarnings(matrix(rep(c(1,0),ceiling(K**2/2)), K,K))
##' Y <- rblockdata(n,mu,sigma=.5)$Y
##' res <- stab.blockSeg(Y, 100, 20)
##'
##' @export stab.blockSeg
stab.blockSeg <- function(Y, nsimu, max.break,max.var = floor(ncol(Y)**2/8),
                          random.break=TRUE,sym.break=FALSE,mc.cores=2, verbose=TRUE){
  if (!(is.matrix(Y)||(class(Y)=="dgeMatrix"))){
    stop("Y must be the observations data (or a transformation)")
  }

  if (!is.numeric(max.break)){
    stop("max.break must be an integer between 1 and n")
  } else if ((max.break<=0)||(max.break>(min(dim(Y))/2))||(length(max.break)!=1)||(floor(max.break)!=max.break)){
    stop("max.break must be an integer between 1 and n/2")
  }
  if (!is.numeric(max.var)){
    stop("max.var must be an integer between 1 and n1 times n2")
  } else if ((max.var<=0)||(max.var>(length(Y)/4))||(length(max.var)!=1)||(floor(max.var)!=max.var)){
    stop("max.var must be an integer between 1 and n1/2 times n2/2")
  }
  if (!is.logical(random.break)){
    stop("random.break must be logical : TRUE if you want that the frist row (resp. column) changes
         randomly")
  }
  if (!is.logical(sym.break)){
    stop("sym.break must be logical : TRUE if you want that the same breaks for the row and columns")
  }
  if (sym.break){
    if (nrow(Y)!=ncol(Y)){
      warning("If sym.break==TRUE, the matrix Y has to have the same number of rows and columns.
              Enforcing to FALSE.")
      sym.break=FALSE
    }
  }
  if (!is.logical(verbose)){
    stop("verbose must be logical : TRUE if you want the details of the procedure")
  }
  if (!is.numeric(mc.cores)){
    stop("max.break must be an integer between 1 and n")
  } else if ((mc.cores<=0)||(length(mc.cores)!=1)||(floor(mc.cores)!=mc.cores)){
    stop("max.break must be a positive integer")
  }

  ## =============================================================
  ## INITIALIZATION & PARAMETERS RECOVERY
  if (Sys.info()[['sysname']] == "Windows") {
    warning("\nWindows does not support fork, enforcing mc.cores to '1'.")
    mc.cores <- 1
  }
  if (random.break){
    stabsel.bloc <- function(i) {
      if (verbose) {setTxtProgressBar(pb, i)}
      ind <- sort(sample(nrow(Y),nrow(Y)/2))
      out <- blockSeg(Y[ind,ind], max.break, max.var, FALSE, FALSE)
      nstep <- length(out@Lambda)
      return(cbind(tabulate(ind[sort(out@RowBreaks[[nstep]])], nrow(Y)),
                   tabulate(ind[sort(out@ColBreaks[[nstep]])], nrow(Y))))
    }
  }else{
    stabsel.bloc <- function(i) {
      if (verbose) {
        setTxtProgressBar(pb, i)
      }
      ind <- sort(sample(nrow(Y), nrow(Y)/2))
      cas <- floor(stats::runif(1)*nrow(Y)/2)
      if (cas!=0){
        ind=ind[c((cas+1):(nrow(Y)/2),1:cas)]
      }
      out <- blockSeg(Y[ind, ind], max.break, max.var, FALSE,
                      FALSE)
      nstep <- length(out@Lambda)
      return(cbind(tabulate(ind[sort(out@RowBreaks[[nstep]][order(ind)])],
                            nrow(Y)), tabulate(ind[sort(out@ColBreaks[[nstep]][order(ind)])],
                                               nrow(Y))))
    }
  }

  if (verbose) {
    cat("\n\nSTABILITY SELECTION for the blockSeg procedure")
    cat("\nRunning",nsimu," subsampling parallely on",mc.cores,"core(s)\n\n")
    pb <- txtProgressBar(min = 0, max = nsimu, style = 3)
  }

  if (mc.cores>1){
    res.cum <- Reduce("+", parallel::mclapply(1:nsimu, stabsel.bloc, mc.cores=mc.cores))
  } else {
    res.cum <- Reduce("+", lapply(1:nsimu, stabsel.bloc))
  }
  cat("\n")
  if (sym.break){
    res.cum=matrix(rep(rowSums(res.cum),2),ncol=2)
  }
  new(Class="stab.blockSeg", RowBreaks=res.cum[,1], ColBreaks=res.cum[,2])
}
