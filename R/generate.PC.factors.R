generate.PC.factors <- function(data, no.factors = 3){
  
  data <- data[, apply(data, 2, sd, na.rm = TRUE) > 0]
  X <- scale(data, center = TRUE)
  XtX <- t(X) %*% X
  
  N <- dim(X)[2]
  eigen_decomposition <- eigen(XtX)
  eigenvalues <- eigen_decomposition$values
  eigenvectors <- eigen_decomposition$vectors
  
  loadings <- eigenvectors[,1:no.factors]*N
  factors <- (X %*% loadings)/N
  colnames(factors) <- paste0("M",1:no.factors)
  macro_factors <- ts(factors,start=time(data)[1],end=rev(time(data))[1],frequency=12)
  
  return(macro_factors)
}