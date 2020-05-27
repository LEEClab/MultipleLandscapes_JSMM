jsmm.mcmc <- function(loglikelihood = loglik, data = data, n.iter = n.iter, n.adapt.iter = iter, n.thin = 1, rotate = TRUE){
  
  # The following objects should be changed to use different prior distributions: mu_z, PSI_S, PSI_V, nu_V, rhos, priorRho
  
  tpm0 <- proc.time() #initial processing time
  
  ## First, we separate some variables from the "data" list
  ns <- data$ns # number of studied species
  np <- data$np # number of species-specific parameters
  nl <- data$nl # number of landscapes
  nc <- data$nc # number of landscapes covariates (+ intercetp)
  
  SP <- data$SP # sequence of studied species
  includePhylogeny <- data$includePhylogeny
  includeTraits <- data$includeTraits
  CC <- data$CC # Phylogenetic correlation matrix
  TT <- data$TT # Species traits matrix
  b <- data$b # landscape covariates matrix
  
  ## And, then, we set up some variables and prior distributions.
  
  ## To include the intercept into the model, we set TT[,1] <- 1 for all species. In the absence of trait information, only the intercept is included.
  if(!includeTraits){
    TT <- matrix(1, nrow = ns)
    nt <- 1
    colnames(TT) <- "intercept"
  }
  
  if(includeTraits){
    TT <- TT
    nt <- ncol(TT)
    # traits <- colnames(TT)
  }
  #########################  
  ## Prior distributions:
  
  ## Zeta
  mu_z <- matrix(0, nrow = np * nc * nt)
  S_z  <- 100 * diag(np * nc * nt) ##### pra que 100 *?
  
  ## Sigma
  PSI_S <- diag(nc*np)
  nu_S  <- nc*np
  
  ## Upsilon
  PSI_V <- diag(np)
  nu_V  <- np
  
  ## Rho
  ## Pre-compute inverse and determinant of the D-matrix
  
  if(includePhylogeny){
    
    rhos <- seq(0, 0.99, 0.01)
    nr   <- length(rhos)
    priorRho  <- rep(1/nr, nr)
    # priorRho <- rep(0.5/(nr - 1), nr)
    # priorRho[1] <- 0.5
    lpriorRho <- log(priorRho)
    iWs <- array(dim = c(ns, ns, nr))
    ldetWs <- numeric(nr)
    
    for(i in 1: nr){
      W <- rhos[i] * CC + (1 - rhos[i]) * diag(ns)
      iWs[ , , i] <- solve(W)
      ldetWs[i]   <- log(det(as.matrix(W)))
    }
  }
  
  ## Before running the MCMC we set initial values for the
  ## parameters to be estimated: Theta, Thetah, Sigma, Upsilon, Zeta, and rho.
  
  if (is.null(data$THETAINIT)){
    XTheta <- numeric(ns * np * nc)
  } else {
    XTheta <- data$THETAINIT # this ns x np matrix should be added to the data list in order to use initial values different from 0
  }
  
  Thetah <- array(0,c(ns, np, nl))
  
  Sigma  <- diag(nc * np)  # identity matrix
  iSI <- solve(Sigma) # inverse of matrix Sigma
  
  V <- diag(np)
  iV <- solve(V)
  
  Z <- matrix(0, nrow = nt, ncol = np * nc)
  z <- as.vector(Z)
  
  # X   <- TT %x% diag(np)
  iVII <- solve(V %x% diag(nl) %x% diag(ns))
  
  Q <- diag(np) %x% diag(nc) %x% TT
  B <- diag(np) %x% b %x% diag(ns) ## The matrix B can be written as B=I_(n_p )*b*I_(n_s ), where the n_l×n_c  matrix b describes the landscape covariates for each landscape.
  
  M   <- B %*% XTheta   # mu_k; expected values of THETAh
  
  if(includePhylogeny){
    RI <- round(nr/2)
    rho <- rhos[RI]
    iW <- iWs[ , ,RI]
    ldetW <- ldetWs[RI]
  }
  
  if(!includePhylogeny){
    rho <- 1
    iW <- diag(ns)
    # ldetW <- log(det(diag(ns)))
  }
  
  iXSI <- iW %x% iSI
  
  ## Initial likelihoods for Thetah
  
  li1 <- loglikelihood(Thetah, data)
  
  RES = as.vector(B %*% XTheta) - as.vector(aperm(Thetah, c(1,3,2)))
  li2 <- -(1/2) * RES %*% iVII %*% RES
  
  print("inital log-likelihood: ")
  print(li1)
  print(li2)
  print("sampling starts")
  
  ## Acceptance rates
  
  ac <- array(0, c(ns, np, 2, nl))
  kk <- array(1, c(ns, np, nl)) # sd for the proposal distributions
  acr <- array(NA, c(ns, np, nl))
  
  ns1s2 <- array(0, c(ns, nl))
  s1 <- array(0, c(ns, np, nl))
  s2 <- array(0, c(ns, np, np, nl))
  la <- array(1, c(ns, np, nl))
  vect <- array(0, c(ns, np, np, nl))
  
  for (ll in 1 : nl){
    for (k in 1 : ns){
      vect[k,,,ll] <- diag(np)
    }
  }
  
  ## Posteriors to be stored
  Post_Thetah <- array(NA, c(n.iter, ns, np, nl)) # array of thetah 
  Post_Theta <- array(NA, c(n.iter, ns * np * nc)) # vector of theta
  Post_Z <- array(NA, c(n.iter, nt * np * nc)) # vector of zeta
  Post_Sigma <- array(NA, c(n.iter, np * nc, np * nc)) # matrix of sigma
  Post_V <- array(NA, c(n.iter, np, np)) # matrix of upsilon
  Post_rho <- numeric(n.iter)
  Post_LIKE <- array(NA, c(n.iter, ns, nl))
  
  ##  MCMC sampling scheme
  
  for (i in 1:(n.iter + n.adapt.iter)){
    
    print(i)
    
    for (ii in 1:n.thin){
      
      ## UPDATE THETAH
      
      for (ll in 1: nl){ # update thetah for each landscape
        
        for (l in 1:np){ # update thetah for each parameter
          
          NTHETAH <- Thetah #_ls
          
          for(k in 1:ns){ # update thetah for each species
            
            nTHETAH <- Thetah[k, ,ll]
            mult <- rnorm(1, mean=0, sd = (kk[k,l,ll]*sqrt(la[k,l,ll])))
            
            for (l2 in 1:np){
              
              nTHETAH[l2] <- nTHETAH[l2] + mult*vect[k,l,l2,ll]
              
            }
            
            NTHETAH[k, ,ll] <- nTHETAH
          }
          
          nli1 <- loglikelihood(NTHETAH, data) # likelihood table ns x 1, for only landscape ll
          
          for (k in 1:ns){
            
            N2THETAH <- Thetah #_ls 
            N2THETAH[k, ,ll] <- NTHETAH[k, ,ll]
            
            RES = as.vector(B %*% XTheta) - as.vector(aperm(N2THETAH, c(1,3,2))) 
            nli2 = -(1/2) * RES %*% iVII %*% RES
            
            ac[k,l,1,ll] <- ac[k,l,1,ll] + 1
            
            if(is.finite(nli1[k,ll]) & is.finite(nli2)){
              
              if(runif(1) < exp(nli1[k,ll] - li1[k,ll] + nli2 - li2)){
                
                Thetah[k,,ll] <- NTHETAH[k,,ll]
                li1[k,ll] <- nli1[k,ll]
                li2 <- nli2
                ac[k,l,2,ll] <- ac[k,l,2,ll] + 1
                
              }
            }
          }
        }
      }
      
      XThetah <- as.vector(aperm(Thetah, c(1,3,2))) # transform to vector; thetah: ns x nl x np
      
      ## The other parameters are sampled directly from their full conditional distribution 
      
      ## UPDATE THETA
      
      Vt <- solve(iXSI + t(B) %*% iVII %*% B)
      Vt <- (Vt + t(Vt))/2
      musTHE = Vt %*% (iXSI %*% (Q %*% z) + t(B) %*% iVII %*% XThetah)
      # set.seed(42)
      XTheta <- as.vector(rmvnorm(1, mean = musTHE, sigma = Vt))

      ## UPDATE Z
      
      Vs <- solve(solve(S_z) + t(Q) %*% iXSI %*% Q)
      Vs <- (Vs + t(Vs))/2
      mus <- Vs %*% (solve(S_z) %*% mu_z + t(Q) %*% iXSI %*% XTheta)
      # set.seed(42)
      z <- as.vector(rmvnorm(1,mean = mus, sigma = Vs))

      ## UPDATE Sigma
      
      RES <- XTheta - as.vector(Q %*% z)
      RES <- matrix(RES, ncol = nc*np, byrow = TRUE)
      A <- t(RES) %*% iW %*% RES
      PSIA_SI <- PSI_S + A
      PSIA_SI <- (PSIA_SI + t(PSIA_SI))/2
      ##set.seed(42)
      Sigma <- riwish((nu_S + ns), PSIA_SI) # Inverse Wishart Matrix Distribution; np*nc x np*nc matrix
      Sigma <- (Sigma + t(Sigma))/2
      iSI = solve(Sigma)
      
      ## UPDATE Upsilon
      
      RES <- XThetah - as.vector(B %*% XTheta)
      RES <- matrix(RES, ncol = np)
      A <- t(RES) %*% solve(diag(nl)%x%diag(ns)) %*% RES
      PSIA_V <- PSI_V + A
      PSIA_V <- (PSIA_V + t(PSIA_V))/2
      ##set.seed(42)
      V <- riwish((nu_V + nl*ns), PSIA_V) # Inverse Wishart Matrix Distribution; np x np matrix
      V <- (V + t(V))/2
      iV = solve(V)
      iVII <- solve(V %x% diag(nl) %x% diag(ns)) ## V %x% solve(diag(nl)) %x% solve(diag(ns)) ???????
      
      ## UPDATE rho
      
      if(includePhylogeny){
        RES <- as.vector(Q %*% z) - XTheta
        likeRho <- numeric(nr)
        
        for(iii in 1:nr){
          likeRho[iii] <- (-1/2) * (np * ldetWs[iii] + RES %*% (iWs[ , , iii] %x% iSI) %*% RES) ### t(RES) ??????? ver supporting info
        }
        
        postRho <- lpriorRho + likeRho
        pr <- exp(postRho) / sum(exp(postRho))
        RI <- sample(seq(1 : nr), size = 1, prob = pr)
        iW <- iWs[ , , RI]
        ldetW <- ldetWs[RI]
        rho <- rhos[RI]
        
        RES <- XThetah - as.vector(B %*% XTheta)
        li2 <- -(1/2) * RES %*% iVII %*% RES
        
      }
      
      iXSI <- iW %x% iSI
      
      ## Adaptation
      
      if (i <= n.adapt.iter){
        
        for(ll in 1:nl){
          for(k in 1:ns){
            
            q <- 1 + exp(-(i*n.thin)/500)
            w <- 1 - 0.1*exp(-(i*n.thin)/500)
            
            acr[k,,ll] <- ac[k,,2,ll]/ac[k,,1,ll]
            kk[k,,ll] <- sapply(kk[k,,ll] * q^(acr[k,,ll] - 0.44), trunca)
            s1[k,,ll] <- s1[k,,ll] + w * Thetah[k,,ll]
            s2[k,,,ll] <- s2[k,,,ll] + w * (Thetah[k,,ll] %*% t(Thetah[k,,ll]))
            ns1s2[k,ll] <- ns1s2[k,ll] + w
            
            if (rotate && ((i * n.thin) > 50)) {
              cov <- (s2[k, , ,ll] - (s1[k, ,ll] %*% t(s1[k, ,ll]) / ns1s2[k, ll])) / (ns1s2[k, ll] - 1)
              met <- cov + 10 ^ (-5) * diag(np)
              lavect <- eigen(met)
              la[k,,ll] <- abs(lavect$values)
              vect[k,,,ll] <- lavect$vectors
            }
            
          }
        }
        
        ac = ac*w
        
      }
      
    } # END OF THINNING LOOP
    
    ## Store the posteriors
    
    if(i > n.adapt.iter) {
      
      Post_Thetah[i-n.adapt.iter,,,] <- Thetah
      Post_Theta[i-n.adapt.iter,] <- XTheta
      Post_Z[i-n.adapt.iter,] <- z
      Post_Sigma[i-n.adapt.iter,,] <- Sigma
      Post_V[i-n.adapt.iter,,] <- V
      Post_rho[i-n.adapt.iter] <- rho
      Post_LIKE[i-n.adapt.iter,,] <- li1
      
    }
    
  }
  
  tpm1 <- proc.time() #final time
  
  ## Organizing outputs
  
  dtpm <- (tpm1[3] - tpm0[3])/60 # time elapsed
  
  rnames_summary <- c("Time elapsed (minutes):", "Total number of iterations:",
                      "Number of adaptation iterations:","Length of thinned posterior",
                      "Number of species", "Number of parameters", "Number of traits")
  results_summary <- matrix(c(as.numeric(dtpm), (n.iter+n.adapt.iter)*n.thin, n.adapt.iter*n.thin, n.iter, ns, np, nt),
                            dimnames = list(rnames_summary,"Model summary"))
  
  posterior <- list(THETAH = Post_Thetah, THETA = Post_Theta, UPSILON = Post_V, ZETA = Post_Z, SIGMA = Post_Sigma, RHO = Post_rho, LIKE = Post_LIKE)
  
  return(list(results_summary = results_summary, posterior = posterior, la = la, vect = vect, ac = ac, kk = kk))
  
}
