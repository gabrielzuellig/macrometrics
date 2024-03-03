testIRFsign <- function(irf_vec, sign_to_check, rev_cond=F){
  
  # this functions checks whether the sign restriction on the kk variable is 
  # verified 
  # INPUT: 
  # - irf_vec = hx1 vector of impulse response to the ith shock of the kth
  #           variable (h = EndMAT(kk,ii)-StartMAT(kk,ii)); 
  # - sign_to_check = +1 Positive 
  #                   -1 Negative 
  # - rev_cond = 1 (checks that conditions are verified as shown) 
  #            = 0 (check whether the reverse are verified)
  # 
  # OUTPUT: 
  # - a = 1 when sign restriction is verified (for the all h periods) 
  #       0 otherwise
  
  a <- 0
  
  if (rev_cond){
    if (sign_to_check == -1){
      ggg <- irf_vec <= 0
    } else if (sign_to_check == 1){
      ggg <- irf_vec >= 0
    }
  } else {
    if (sign_to_check == 1){
      ggg <- irf_vec >= 0
    } else if (sign_to_check == -1){
      ggg <- irf_vec <= 0
    }
  }
  
  if (sum(ggg) == length(irf_vec)){
    a <- 1
  }
  return(a)
  
}