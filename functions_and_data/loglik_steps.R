loglik <- function(THETAH,data){
  
  SP <- data$SP
  ns <- length(SP) # number of species
  nl <- dim(THETAH)[3] # number of landscapes
  
  like <- matrix(0, ns, nl)
  
  for (ll in 1:nl){ # for each landscape
    
    # bird movement data
    all.tracks <- data$tracks [[ll]] 
    SPS <- all.tracks[,1]
    
    #environmental covariates
    semiopen <- data$semiopen [[ll]]
    forest <- data$forest [[ll]]
    d <- data$d
    
    for (k in 1: ns) {
      thetah <-  THETAH[k, ,ll] # theta is a 3 x 1 vector with the species-specific parameters for species k
      tracks <- all.tracks[SPS == SP[k],]
      
      alpha <- exp(thetah[1]) # typical step length made by species k
      beta_s <- thetah[2] # affinity of the species to the corridor habitat compared to the open areas
      beta_f <- thetah[3] # affinity of the species to the forest habitat compared to the open areas
      li <- 0
      
      idt <- tracks[ ,2]
      uid <- unique(idt)
      ta <- tracks[ ,3]
      
      if(length(uid) == 0 || is.na(uid)){
        
        li <- 0
        
      }else{
        
        for(ss in 1: length(uid)){ #for each one of the tracks
          lta <- ta[idt == uid[ss]] #data for the focal track
          
          for(i in 1: (length(lta)-1)) { #for each movement step
            x <- lta[i]
            pr <- exp(-d[x,]/alpha) * exp(beta_s * semiopen) * exp(beta_f * forest)
            pr <- pr / sum(pr)
            pr <- pr + 10^(-10)
            pr <- pr / sum(pr)
            li <- li + log(pr[lta[i + 1]])
          }
        }
        
        like[k, ll] <- li
      }
    }
  }
  return(like)
}