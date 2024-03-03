estimateVAR <- function(data, p, c_case, exdata=NULL){
  
  # Prep
  t <- nrow(data)
  n <- ncol(data)
  
  # Y matrix
  Y <- unname(as.matrix(data))
  
  # X matrix: Generate lags
  X <- NULL
  for (pp in 1:p){
    X <- cbind(X, rbind(matrix(NA, nrow=pp, ncol=n), Y[1:(nrow(data)-pp),]))
  }
  # lags introduce some NA's => truncate
  X <- X[(p+1):nrow(X), ]
  Y <- Y[(p+1):nrow(Y), ]
  # add constant and/or trend to X matrix
  if (c_case == 0){
    c <- NULL
  } else if (c_case == 1){
    c <- matrix(1, nrow=nrow(X), ncol=1)
  } else if (c_case == 2){
    c <- cbind(matrix(1, nrow=nrow(X), ncol=1), matrix(1:nrow(X), nrow=nrow(X), ncol=1))
  } else {
    print('c_case variable needs to be set to 0, 1 or 2.')
  }
  # add exogenous variables
  Xex <- NULL
  if (!is.null(exdata)){
    for (pp in 1:p){
      Xex <- cbind(Xex, rbind(matrix(NA, nrow=pp, ncol=ncol(exdata)), matrix(exdata[1:(nrow(data)-pp),], ncol=ncol(exdata))))
    }
    Xex <- Xex[(p+1):nrow(Xex), ]
  }
  X <- cbind(c, X, Xex)
  
  # estimation
  if (qr(t(X) %*% X)$rank == 1){
    Sys.sleep(1)
    stop('System is computationally singular.')
  }
  A <- solve(t(X) %*% X) %*% (t(X) %*% Y) # OLS: beta = inv(x'x) * (x'y)
  u <- Y - X%*%A # reduced-form residuals
  Omega <- cov(u)  # variance-covariance matrix
  if (sum(is.na(Omega)) > 0){
    Sys.sleep(1)
    stop('Variance-covariance matrix not positive definitive')
  } else if (sum(eigen(Omega)$values > 0) < n){
    Sys.sleep(1)
    stop('Variance-covariance matrix not positive definitive')
  } else {
    S <- t(chol(Omega))  # lower-trinagular matrix (Cholesky decomposition)
  }
  
  # Housekeeping
  return(VAR = list('data' = data,
                    'exdata' = exdata,
                    'X' = X[, (c_case+1):(c_case+n*p)],
                    'Xex' = Xex,
                    'c_case' = c_case,
                    'p' = p,
                    't' = t,
                    'n' = n,
                    'A' = A,
                    'u' = u,
                    'Omega' = Omega,
                    'S' = S))
}
