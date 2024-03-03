dyn_multipliers <- function(VAR, h){
  
  c_case <- VAR$c_case
  p <- VAR$p 
  n <- VAR$n
  A <- VAR$A
  
  AR <- t(A[(c_case+1):nrow(A), ]) # AR coefficient [A1,A2,...,Ap]
  AR_3d <- array(NA, dim=c(n,n,p))
  for (pp in 1:p){
    AR_3d[,,pp] <- AR[,((pp-1)*n+1):(pp*n)]
  }
  
  C <- array(NA, dim=c(n,n,h+p))
  for (pp in 1:p){
    C[,,pp] <- matrix(0, n, n)
  }
  C[,,p+1] <- diag(n)
  for (pp in (p+2):(p+h)){
    acc <- 0
    for (j in 1:p){
      acc <- acc + AR_3d[,,j]%*%C[,,pp-j]
    }
    C[,,pp] <- acc
  }
  
  # reindicize
  return(C[,,(p+1):(dim(C)[3])])
  
}





