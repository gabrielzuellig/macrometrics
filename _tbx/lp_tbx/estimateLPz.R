estimateLPz <- function(data, z, h, p, c_case, exdata, alpha=90, NWSE=T){
  
  # Prep parameters/containers
  t <- nrow(data)
  n <- ncol(data)
  gamma <- matrix(NA, nrow=h, ncol=n)
  gammabands <- array(NA, dim=c(h,n,2))

  zscore <- qnorm(1-(1-alpha/100)/2)
  
  # Y matrix
  Y <- unname(as.matrix(data))
  
  # X matrix: Generate lags
  X <- NULL
  for (pp in 1:p){
    X <- cbind(X, rbind(matrix(NA, nrow=pp, ncol=n), Y[1:(nrow(data)-pp),]))
  }
  
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
      Xex <- cbind(Xex, rbind(matrix(NA, nrow=pp, ncol=ncol(Xex)), exdata[1:(nrow(data)-pp),]))
    }
    Xex <- Xex[(p+1):nrow(Xex), ]
  }
  X <- cbind(c, X, Xex)
  
  ## Estimation
  print('Progress...')
  counter <- 1
  pb <- txtProgressBar(min = 0, max = n*h, initial = 0, style=3) 
  for (nn in 1:n){
    for (hh in 0:(h-1)){
      
      # prepare LHS variable (changes in each loop while RHS stays constant)
      lhs = as.matrix(c(Y[(hh+1):nrow(Y),nn], rep(NA, hh)), ncol=1)
      
      # drop observations with any NA's
      regdata <- data.frame(cbind(lhs, z, X))
      nona <- rowSums(is.na(regdata)) == 0
      regdata <- regdata[nona, ]
      
      # reduced-form estimation (OLS with Newey-West standard errors)
      reg <- lm(X1 ~ 0 + ., data=regdata)
      if (NWSE){
        vcovNW <- unname(NeweyWest(reg))
        se <- sqrt(diag(vcovNW))
      } else {
        se <- summary(reg)$coefficients[,2]
      }
      
      gamma[hh+1,nn] <- reg$coefficients[1]
      gammabands[hh+1,nn,1] <- reg$coefficients[1] - zscore*se[1]
      gammabands[hh+1,nn,2] <- reg$coefficients[1] + zscore*se[1]
      
      # show progress report
      setTxtProgressBar(pb,counter)
      counter <- counter + 1
      
    }
  }
  
  # Housekeeping
  return(LP = list('data' = data,
                   'exdata' = exdata,
                   'X' = X,
                   'Xex' = Xex,
                   'c_case' = c_case,
                   'p' = p,
                   't' = t,
                   'n' = n,
                   'gamma' = gamma,
                   'gammabands' = gammabands,
                   'h' = h))
  
}