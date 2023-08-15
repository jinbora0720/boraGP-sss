create_BF <- function(coords.T, coords.S, nn.indx.T, sig_sq, phi, nu = NULL,
                      base_cov = "exponential", tol = 0, cores = 2) {
  # coords.T: coords for observed locations
  # coords.S: coords for reference locations
  # nn.indx.T: matrix of nearest neighbor index for coords.T
  # sig_sq: spatial variance
  # phi: decay parameter
  # nu: smoothness parameter, applicable only for base_cov = "matern"
  # base_cov = "exponential" or "matern"
  # tol: tolerance for solve
  # cores: number of cores for parallelization
  
  # set the number of cores
  doParallel::registerDoParallel(cores = cores)
  
  # number of observed locations
  n <- nrow(coords.T)
  
  # number of reference locations
  k <- nrow(coords.S)
  
  # create B and F
  if (base_cov == "exponential") {
    BFlist <- foreach::foreach (i = 1:n) %dopar% {
      coords_tmp <- rbind(coords.T[i,], coords.S[nn.indx.T[i,],])
      Cor <- exp(-phi*as.matrix(dist(coords_tmp)))
      if (tol > 0) {
        CinvC <- Cor[1,-1]%*%solve(Cor[-1,-1] + diag(tol, nrow(Cor)-1))
      } else {
        CinvC <- Cor[1,-1]%*%solve(Cor[-1,-1])
      }
      B <- CinvC
      Fmat <- sig_sq*(Cor[1,1]-CinvC%*%Cor[-1,1])
      return(list(B = B, Fmat = Fmat))
    }
  } else {
    BFlist <- foreach::foreach (i = 1:n) %dopar% {
      coords_tmp <- rbind(coords.T[i,], coords.S[nn.indx.T[i,],])
      Cor <- boraGP::Cov_matern(as.matrix(dist(coords_tmp)), sigmasq = 1,
                                phi = phi, nu = nu)
      if (tol > 0) {
        CinvC <- Cor[1,-1]%*%solve(Cor[-1,-1] + diag(tol, nrow(Cor)-1))
      } else {
        CinvC <- Cor[1,-1]%*%solve(Cor[-1,-1])
      }
      B <- CinvC
      Fmat <- sig_sq*(Cor[1,1]-CinvC%*%Cor[-1,1])
      return(list(B = B, Fmat = Fmat))
    }
  }
  
  Felem <- sapply(BFlist, function(x) as.numeric(x$Fmat))
  Fmat <- diag(Felem)
  B <- matrix(0, n, k)
  for (i in 1:n) {
    B[i, nn.indx.T[i,]] <- BFlist[[i]]$B
  }
  
  return(list(B = B, Fmat = Fmat)) 
}

spNNGP_freeS <- function(y, X = NULL, coords.T, 
                         X.S = NULL, coords.S, X.0 = NULL, coords.0 = NULL,
                         starting, 
                         tuning = list("phi" = 0.1, "sigma.sq" = 0.1, 
                                       "tau.sq" = 0.1, "nu" = 0), 
                         priors, 
                         cov.model = "exponential",
                         sub.sample = list(start = 1, end = 1000, thin = 1),
                         cores, neighbor.info, tol = 0,
                         RAMH = TRUE, verbose = TRUE, predict = TRUE,
                         debug = list(fix_cov_params = FALSE)) {
  
  # spNNGP with user-chosen S 
  # family = "gaussian"
  # method = "response"
  # y: response at observed locations (coords.T)
  # X: covariates at observed locations (coords.T)
  # coords.T: observed locations
  # X.S: covariates for reference locations
  # coords.S: reference locations
  # X.0: covariates for new locations
  # coords.0: new locations for prediction (neither S nor T)
  # neighbor.info: list including nn.indx.0, ord, and n.indx
  # cores: number of cores for parallelization
  # if cov.model = "matern", nu must be fixed
  # RAMH = TRUE: Robust Adaptive Metropolis Hastings (Vihola, 2012)
  # predict = TRUE: return p.y.0
  # tol: tolerance for solve, applies only when predict = TRUE
  
  ##################
  # pre-processing #
  ##################
  # set the number of cores
  doParallel::registerDoParallel(cores = cores)
  
  # error 
  if (is.null(coords.0) & predict) {
    error("coords.0 should be provided for prediction.")
  }
  
  # dimension
  y <- as.vector(y)
  n <- length(y)
  if (is.null(X)) {
    X <- matrix(1, nrow = n, ncol = 1)
    x.names <- "(Intercept)"
    cat("Run intercept-only model.\n")
  }
  x.names <- colnames(X)
  p <- ncol(X)
  
  n.0 <- nrow(coords.0)
  if (is.null(X.0) & predict) {
    X.0 <- matrix(1, nrow = n.0, ncol = 1)
  }
  k <- nrow(coords.S)
  if (is.null(X.S) & predict) {
    X.S <- matrix(1, nrow = k, ncol = 1)
  }
  
  # distinguish nn.indx.0 for coords.T and coords.0
  nn.indx.T <- neighbor.info$nn.indx.0[1:n,]
  nn.indx.0 <- neighbor.info$nn.indx.0[-c(1:n),]
  
  col.names <- c("sigma.sq", "tau.sq", "phi")
  if(cov.model == "matern"){
    col.names <- c(col.names, "nu")
  }
  
  # prior 
  names(priors) <- tolower(names(priors))
  sigma.sq.IG <- priors[["sigma.sq.ig"]]
  tau.sq.IG <- priors[["tau.sq.ig"]]
  phi.Unif <- priors[["phi.unif"]]
  nu.Unif <- priors[["nu.unif"]]
  
  # starting value
  names(starting) <- tolower(names(starting))   
  beta <- as.vector(coefficients(lm(y~X-1)))
  sigma.sq <- starting[["sigma.sq"]][1]
  tau.sq <- starting[["tau.sq"]][1]
  phi <- starting[["phi"]][1]
  nu <- starting[["nu"]][1]
  if (cov.model == "exponential") {
    nu <- 0.5
  }
  
  # Sigma
  Ctilde <- boraGP:::create_Ctilde(coords = coords.S, 
                                   neighbor.info = neighbor.info, 
                                   sig_sq = sigma.sq, phi = phi, nu = nu,
                                   base_cov = cov.model, tol = 0, cores = cores)
  BtFt <- create_BF(coords.T = coords.T, coords.S = coords.S, 
                    nn.indx.T = nn.indx.T, 
                    sig_sq = sigma.sq, phi = phi, nu = nu,
                    base_cov = cov.model, tol = 0, cores = cores)
  Bt <- BtFt$B
  Ft <- BtFt$Fmat
  Sigma <- tcrossprod(Bt%*%Ctilde, Bt) + Ft + diag(tau.sq, n)
  Sigma_inv <- solve(Sigma)
  
  # tuning (nu is fixed)
  names(tuning) <- tolower(names(tuning))
  sigma.sq.tuning <- tuning[["sigma.sq"]][1]
  tau.sq.tuning <- tuning[["tau.sq"]][1]
  phi.tuning <- tuning[["phi"]][1]
  
  if (RAMH) {
    # tuning parameters for Robust Adaptive MH
    Sn <- diag(1,3)
  } else {
    # acceptance rate 
    acc_rate <- 0
  }
  
  # mcmc 
  start <- sub.sample$start
  end <- sub.sample$end
  thin <- sub.sample$thin
  s.indx <- seq(start, end, by = thin)
  save <- length(s.indx)
  
  ############
  # sampling #
  ############
  out <- list()
  out$p.theta.samples <- matrix(0, nrow = save, ncol = length(col.names))
  out$p.beta.samples <- matrix(0, nrow = save, ncol = p)
  if (predict) {
    out$p.y.0 <- matrix(0, nrow = n.0, ncol = save)
  }
  
  ptm <- proc.time()
  
  if (verbose) {
    cat('burnin...\n')
    pb = txtProgressBar(style=3,width=50)
    pb2 = txtProgressBar(style=3,width=50)
  }
  
  for (s in 1:end) {
    # update beta
    V_b <- solve(crossprod(X, Sigma_inv%*%X))
    L_b <- t(chol(V_b))
    mu_b <- crossprod(X, Sigma_inv%*%y)
    beta <- V_b%*%mu_b + L_b%*%rnorm(p)
    Xb <- X%*%beta
    
    # update (tau.sq, sigma.sq, phi) with Robust Adaptive MH (Vihola)
    if (!debug$fix_cov_params) {
      ## propose
      if (RAMH) {
        U <- rnorm(3)
        logstar <- c(log(tau.sq),
                     log(sigma.sq),
                     log(phi - phi.Unif[1]) - log(phi.Unif[2] - phi)) + Sn%*%U
      } else {
        logstar <- c(log(tau.sq),
                     log(sigma.sq),
                     log(phi - phi.Unif[1]) - log(phi.Unif[2] - phi)) +
          diag(c(sqrt(tau.sq.tuning), sqrt(sigma.sq.tuning), sqrt(phi.tuning)))%*%rnorm(3)
      }
      tau.sq.star <- exp(logstar[1])
      sigma.sq.star <- exp(logstar[2])
      phi.star <- (exp(logstar[3])*phi.Unif[2] + phi.Unif[1])/(1+exp(logstar[3]))
      
      ## log likelihood
      Ctilde.star <- boraGP:::create_Ctilde(coords = coords.S, 
                                            neighbor.info = neighbor.info, 
                                            sig_sq = sigma.sq.star, 
                                            phi = phi.star, nu = nu,
                                            base_cov = cov.model, tol = 0, 
                                            cores = cores)
      BtFt.star <- create_BF(coords.T = coords.T, coords.S = coords.S, 
                             nn.indx.T = nn.indx.T, 
                             sig_sq = sigma.sq.star, phi = phi.star, nu = nu,
                             base_cov = cov.model, tol = 0, cores = cores)
      Bt.star <- BtFt.star$B
      Ft.star <- BtFt.star$Fmat
      Sigma.star <- tcrossprod(Bt.star%*%Ctilde.star, Bt.star) + 
        Ft.star + diag(tau.sq.star, n)
      Sigma.star_inv <- solve(Sigma.star)
      
      Ll <- -0.5*det(Sigma, log = TRUE)-0.5*crossprod(y-Xb, Sigma_inv%*%(y-Xb))
      Ll.star <- -0.5*det(Sigma.star, log = TRUE)-
        0.5*crossprod(y-Xb, Sigma.star_inv%*%(y-Xb))
      
      ## compute 
      logr <- (invgamma::dinvgamma(tau.sq.star, shape = tau.sq.IG[1], 
                                   rate = tau.sq.IG[2], log = TRUE) + 
                 invgamma::dinvgamma(sigma.sq.star, shape = sigma.sq.IG[1], 
                                     rate = sigma.sq.IG[2], log = TRUE) + 
                 Ll.star +
                 log(tau.sq.star) + log(sigma.sq.star) + 
                 log(phi.star - phi.Unif[1]) + log(phi.Unif[2] - phi.star)) - 
        (invgamma::dinvgamma(tau.sq, shape = tau.sq.IG[1], 
                             rate = tau.sq.IG[2], log = TRUE) + 
           invgamma::dinvgamma(sigma.sq, shape = sigma.sq.IG[1], 
                               rate = sigma.sq.IG[2], log = TRUE) + 
           Ll +
           log(tau.sq) + log(sigma.sq) + 
           log(phi - phi.Unif[1]) + log(phi.Unif[2] - phi))
      
      ## accept or reject 
      if (log(runif(1)) < logr) {
        tau.sq <- tau.sq.star
        sigma.sq <- sigma.sq.star
        phi <- phi.star 
        Sigma <- Sigma.star
        Sigma_inv <- Sigma.star_inv
        Bt <- Bt.star
        Ctilde <- Ctilde.star
        if (!RAMH) {
          acc_rate <- acc_rate + 1/n.samples
        }
      }
      
      ## update Sn
      if (RAMH) {
        etan <- min(1, 3*s^(-2/3))
        an <- min(1, exp(logr))
        U_l2norm <- sum(U^2)
        Sn <- t(chol(tcrossprod(Sn%*%(diag(1,3) + tcrossprod(U)*etan*(an - 0.234)/U_l2norm), Sn)))
      }
    }
    
    ########
    # save #
    ########
    if ((s > burn) & (s-burn) %% thin == 0) {
      if(cov.model == "matern"){
        out$p.theta.samples[(s-burn)/thin,] <- c(sigma.sq, tau.sq, phi, nu)
      } else {
        out$p.theta.samples[(s-burn)/thin,] <- c(sigma.sq, tau.sq, phi)
      }
      out$p.beta.samples[(s-burn)/thin,] <- beta
      
      if (predict) {
        # for coords.0
        if (cov.model == "exponential") {
          BFlist <- foreach::foreach (i = 1:n.0) %dopar% {
            coords_tmp <- rbind(coords.0[i,], coords.S[nn.indx.0[i,],])
            Cor <- exp(-phi*as.matrix(dist(coords_tmp)))
            if (tol > 0) {
              CinvC <- Cor[1,-1]%*%solve(Cor[-1,-1] + diag(tol, nrow(Cor)-1))
            } else {
              CinvC <- Cor[1,-1]%*%solve(Cor[-1,-1])
            }
            B <- CinvC
            Fmat <- sigma.sq*(Cor[1,1]-CinvC%*%Cor[-1,1])
            v0 <- as.numeric(tcrossprod(B%*%Ctilde[nn.indx.0[i,],nn.indx.0[i,]], B) +
                               Fmat + tau.sq)
            V0y <- tcrossprod(B%*%Ctilde[nn.indx.0[i,],], Bt)
            V0y_inv <- V0y%*%Sigma_inv
            v0gy <- v0 - tcrossprod(V0y_inv, V0y)
            return(as.numeric(X.0[i,]%*%beta + V0y_inv%*%(y-Xb) + 
                                sqrt(v0gy)*rnorm(1)))
          }
        } else {
          BFlist <- foreach::foreach (i = 1:n.0) %dopar% {
            coords_tmp <- rbind(coords.0[i,], coords.S[nn.indx.0[i,],])
            Cor <- boraGP::Cov_matern(as.matrix(dist(coords_tmp)), sigmasq = 1,
                                      phi = phi, nu = nu)
            if (tol > 0) {
              CinvC <- Cor[1,-1]%*%solve(Cor[-1,-1] + diag(tol, nrow(Cor)-1))
            } else {
              CinvC <- Cor[1,-1]%*%solve(Cor[-1,-1])
            }
            B <- CinvC
            Fmat <- sigma.sq*(Cor[1,1]-CinvC%*%Cor[-1,1])
            v0 <- as.numeric(tcrossprod(B%*%Ctilde[nn.indx.0[i,],nn.indx.0[i,]], B) +
                               Fmat + tau.sq)
            V0y <- tcrossprod(B%*%Ctilde[nn.indx.0[i,],], Bt)
            V0y_inv <- V0y%*%Sigma_inv
            v0gy <- v0 - tcrossprod(V0y_inv, V0y)
            return(as.numeric(X.0[i,]%*%beta + V0y_inv%*%(y-Xb) + 
                                sqrt(v0gy)*rnorm(1)))
          }
        }
        out$p.y.0[,(s-burn)/thin] <- unlist(BFlist)
      }
    }
    
    if (verbose) {
      setTxtProgressBar(pb, s/burn)
      if (s > burn) {
        if (s == burn + 1) {
          close(pb)
          cat('saving...\n')
          setTxtProgressBar(pb2,(s-burn)/(end-burn))
        } else setTxtProgressBar(pb2,(s-burn)/(end-burn))
      }
    }
  }
  
  out$run.time <- proc.time() - ptm
  out$p.theta.samples <- mcmc(out$p.theta.samples)
  out$p.beta.samples <- mcmc(out$p.beta.samples)
  colnames(out$p.theta.samples) <- col.names
  colnames(out$p.beta.samples) <- x.names
  if (!RAMH) {
    out$acc_rate <- acc_rate
  }
  
  return(out)
}
