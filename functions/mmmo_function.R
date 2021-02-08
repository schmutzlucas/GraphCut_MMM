mmmo <- function(ref_present, ref_future, var_present, var_future) {
  width <- ncol(ref_present) ### Number of longitudes
  height <- nrow(ref_present) ### Number of latitudes
  nlabs <- dim(var_present)[3] ### Number of longitudes
  
  dim(ref_present) <- c(height * width, 1) 
  dim(var_present) <- c(height * width, nlabs) 

  C <- cbind(rep(1,nlabs), diag(nlabs))
  b <- c(1,rep(0,nlabs))
  d <- t(ref_present) %*% var_present  
  Dmat <- (t(var_present) %*% var_present)
  sc <- norm(Dmat,"2")
  weights <- quadprog::solve.QP(Dmat = Dmat/sc, dvec=d/sc, Amat=C, bvec=b, meq=0, factorized=FALSE)$solution

  var_mmmo <- apply(var_future, 1:2, weighted.mean, w = weights)
  bias_mmmo <- if(is.null(ref_future)) NULL else (var_mmmo - ref_future)
  bias_var_mmmo <- list("Var" = var_mmmo, "Bias" = bias_mmmo)
  return(bias_var_mmmo)
}
