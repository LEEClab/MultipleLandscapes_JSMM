loglik <- function(THETAH,data){
  
  SP <- data$SP
  ns <- length(SP) # number of species
  nl <- dim(THETAH)[3] # number of landscapes
  
  like <- matrix(0, ns, nl)
  
  for (ll in 1:nl){ # for each landscape
    
    # bird movement data
    all.tracks <- data$tracks [[ll]] 
    SPS <- all.tracks[,1]
    
    # #environmental covariates
    # semiopen <- data$semiopen [[ll]]
    # forest <- data$forest [[ll]]
    intercept <- 1
    
    for (k in 1: ns) {
      thetah <-  THETAH[k, ,ll] # theta is a 3 x 1 vector with the species-specific parameters for species k
      tracks <- all.tracks[SPS == SP[k],]
      
      beta_i <- thetah[1] # typical perching time of species k
      # beta_s <- thetah[2] # perching time in semi-open habitat
      # beta_f <- thetah[3] # perching time in forest
      sd <- exp(thetah[2])
      
      pt <- tracks[,3]
      # loc <- tracks[,2]
      
      pr <- numeric(length(pt))
      L <- numeric(length(pt))
      li <- 0
      
      if(length(pt) == 0 || is.na(pt)){
        
        li <- 0
        
      }else{
        
        for(ss in 1: length(pt)){ # for each one of the time lapses
          
          # x <- loc[ss]
          
          L[ss] <- sum(c(intercept) * # , forest[x], semiopen[x]) * 
                         c(beta_i)) # , beta_f, beta_s))
          
          pr[ss] <- dnorm(log(pt[ss]), L[ss], sd, log = TRUE)
          
          if(is.na(pr[ss])) pr[ss] <- 0
        }
        
        
        li <- li + sum(pr)
        
      }
      
      like[k, ll] <- li
    }
  }
  
  return(like)
  
}