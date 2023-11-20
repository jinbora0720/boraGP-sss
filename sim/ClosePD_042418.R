#this draft moves cross-validation subsetting of gamma matrices from semivarioNE
#to KrigeNE. Also objects are removed from KrigeNE to make it more efficient and
#take up less RAM

##Finally added the function krigeCV to perform Kriging and produce Cross-validation metrics
#while making sure Kriging objects don't take up too much space during simulations!



variogNE<-function (geodata, coords = geodata$coords, data = geodata$data, 
                    uvec = "default", breaks = "default", trend = "cte", lambda = 1, 
                    option = c("bin", "cloud", "smooth"), estimator.type = c("classical", 
                                                                             "modulus"), nugget.tolerance, max.dist, pairs.min = 2, 
                    bin.cloud = TRUE, direction = "omnidirectional", tolerance = pi/8,  #BD defaulted bin.cloud to "TRUE", will make life easier...... but also make the object bigger
                    unit.angle = c("radians", "degrees"), angles = FALSE, messages,
                    NonEuc = FALSE, NonEucDist,                                                 #BD My Non-Euclidean Option and distance vector... perhaps I can create a new object with a third element in the list?
                    ...) 
{
  if (missing(geodata)) 
    geodata <- list(coords = coords, data = data)
  call.fc <- match.call()
  if (missing(messages)) 
    messages.screen <- as.logical(ifelse(is.null(getOption("geoR.messages")), 
                                         TRUE, getOption("geoR.messages")))
  else messages.screen <- messages
  keep <- list(...)
  if (is.null(keep$keep.NA)) 
    keep.NA <- FALSE
  else keep.NA <- keep$keep.NA
  unit.angle <- match.arg(unit.angle)
  
  if(NonEuc == TRUE && missing(NonEucDist))                           #BD
    stop("Non-Euclidean distance values must be provided")      #BD
  
  if(NonEuc == TRUE && class(NonEucDist) != "dist")            #BD
    stop("Non-Euclidean distance must be an object of class 'dist'")    #BD
  
  #BD will have to decide later if we want to be able to compute directional variograms
  if (NonEuc == TRUE && direction != "omnidirectional")                        #BD
    stop("directional variograms not allowed with non-Euclidean distances") #BD
  
  #BD this section is all directional information 
  if (mode(direction) == "numeric") {
    if (length(direction) > 1) 
      stop("only one direction is allowed")
    if (length(tolerance) > 1) 
      stop("only one tolerance value is allowed")
    if (unit.angle == "degrees") {
      ang.deg <- direction
      ang.rad <- (ang.deg * pi)/180
      tol.deg <- tolerance
      tol.rad <- (tol.deg * pi)/180
    }
    else {
      ang.rad <- direction
      ang.deg <- (ang.rad * 180)/pi
      tol.rad <- tolerance
      tol.deg <- (tol.rad * 180)/pi
    }
    if (ang.rad > pi | ang.rad < 0) 
      stop("direction must be an angle in the interval [0,pi[ radians")
    if (tol.rad > pi/2 | tol.rad < 0) 
      stop("tolerance must be an angle in the interval [0,pi/2] radians")
    if (tol.deg >= 90) {
      direction <- "omnidirectional"
      cat("variog: computing omnidirectional variogram\n")
    }
    else {
      if (messages.screen) {
        cat(paste("variog: computing variogram for direction = ", 
                  round(ang.deg, digits = 3), " degrees (", round(ang.rad, 
                                                                  digits = 3), " radians)\n", sep = ""))
        cat(paste("        tolerance angle = ", round(tol.deg, 
                                                      digits = 3), " degrees (", round(tol.rad, digits = 3), 
                  " radians)\n", sep = ""))
      }
    }
  }
  
  else if (messages.screen) 
    cat("variog: computing omnidirectional variogram\n")
  
  coords <- as.matrix(coords)
  data <- as.matrix(data)
  if (nrow(coords) != nrow(data)) 
    stop("coords and data have incompatible dimensions")
  
  data.var <- apply(data, 2, var)
  n.data <- nrow(coords)
  n.datasets <- ncol(data)
  data <- drop(data)
  option <- match.arg(option)
  estimator.type <- match.arg(estimator.type)
  
  if (abs(lambda - 1) > 1e-04) {
    if (abs(lambda) < 1e-04) 
      data <- log(data)
    else data <- ((data^lambda) - 1)/lambda
  }
  
  #BD probably can't deal with trends right now, will block the function from using them
  if (NonEuc == TRUE && trend != "cte")                        #BD
    stop("cannot currently specify trend in mean model with non-Euclidean distances (Only Ordinary Kriging allowed at this time), instead input residuals from a spatially-indexed 'spatially independent' model.")              #BD
  
  xmat <- unclass(trend.spatial(trend = trend, geodata = geodata))
  if (nrow(xmat) != n.data) 
    stop("coords and trend have incompatible sizes")
  if (trend != "cte") {
    if (is.vector(data)) {
      temp.fit <- lm(data ~ xmat + 0)
      beta.ols <- temp.fit$coeff
      data <- temp.fit$residuals
      temp.fit <- NULL
      names(data) <- NULL
    }
    else {      #BD I'm guessing this is if you have more than one realization of the data?
      only.res <- function(y, x) lm(y ~ xmat + 0)$residuals
      data <- apply(data, 2, only.res, x = xmat)
      only.beta <- function(y, x) lm(y ~ xmat + 0)$coef
      beta.ols <- apply(data, 2, only.beta, x = xmat)
    }
  }
  
  else beta.ols <- colMeans(as.matrix(data))
  if (NonEuc == FALSE){                            #BD keeping old code in this if statement
    u <- as.vector(dist(as.matrix(coords))) #BD old code moved
  }
  else u <- as.vector(NonEucDist)                       #BD attaching NonEuc distance instead
  if (missing(nugget.tolerance) || nugget.tolerance < 1e-11) {
    nugget.tolerance <- 1e-12
    nt.ind <- FALSE
  }
  else {
    if (mode(nugget.tolerance) != "numeric") 
      stop("nugget.tolerance must be numeric")
    nt.ind <- TRUE
  }
  if (any(is.na(u)))                                                      #BD
    warning("some distance values are NA, check that this is correct") #BD
  min.dist <- min(u, na.rm=TRUE)                      #BD allows for NA, needed for 3rd semivariogram test
  if (min.dist < nugget.tolerance) 
    nt.ind <- TRUE
  
  
  if (direction != "omnidirectional" | angles) {
    u.ang <- .C("tgangle", as.double(as.vector(coords[, 1])), 
                as.double(as.vector(coords[, 2])), as.integer(dim(coords)[1]), 
                res = as.double(rep(0, length(u))), PACKAGE = "geoR")$res
    if (any(is.na(u.ang))) 
      stop("NA returned in angle calculations maybe due to co-located data")
    u.ang <- atan(u.ang)
    u.ang[u.ang < 0] <- u.ang[u.ang < 0] + pi
  }
  
  
  
  if (option == "bin" && bin.cloud == FALSE && direction == 
      "omnidirectional") {    # BD very specific situation, the else chunk below is mainly for when bin.cloud == TRUE, probably best to just say don't let this happen for mine, since it's nice to have the bin.cloud anyway!
    if(NonEuc == TRUE){                                                         #BD
      stop("to compute Non-Euclidean variograms please compute cloud  #BD
                             values for each class (i.e.  bin.cloud==TRUE)")            #BD
    }                                                               #BD
    
    else{   
      if (missing(max.dist)) 
        umax <- max(u, na.rm=TRUE)                 #BD allowing for NA values in distance
      else umax <- max(u[u < max.dist], na.rm=TRUE)   #BD same here
      dbins <- geoR:::.define.bins(max.dist = umax, uvec = uvec, breaks = breaks,  
                                   nugget.tolerance = nugget.tolerance) #BD added ":::" to access hidden function
      uvec <- dbins$uvec
      bins.lim <- dbins$bins.lim
      nbins <- length(bins.lim) - 1
      if (missing(max.dist)) 
        max.dist <- max(bins.lim, na.rm=TRUE)           #BD still allowing for NA values in distance
      if (bins.lim[1] < 1e-16) 
        bins.lim[1] <- -1
      bin.f <- function(data) {
        cbin <- vbin <- sdbin <- rep(0, nbins)
        #BD This is C code that calculates variogram values using Euclidean distance. This is the option if just doing a Euclidean distance.
        .C("binit", as.integer(n.data), as.double(as.vector(coords[,   
                                                                   1])), as.double(as.vector(coords[, 2])), as.double(as.vector(data)), 
           as.integer(nbins), as.double(as.vector(bins.lim)), 
           as.integer(estimator.type == "modulus"), as.double(max.dist), 
           cbin = as.integer(cbin), vbin = as.double(vbin), 
           as.integer(TRUE), sdbin = as.double(sdbin), PACKAGE = "geoR")[c("vbin", 
                                                                           "cbin", "sdbin")]
      }
      
      ###BD what I want to do is extract distance from the vectorized distance matrix (object 'u'), for a given two points, calculate the semivariance and then bin them for each variance,
      #BD This is already done below actually...
      
      result <- array(unlist(lapply(as.data.frame(data), bin.f)), 
                      dim = c(nbins, 3, n.datasets))
      indp <- (result[, 2, 1] >= pairs.min)
      result[!indp, 1, ] <- NA
      if (bins.lim[1] < 0) 
        bins.lim[1] <- 0
      if (!nt.ind) {
        uvec <- uvec[-1]
        indp <- indp[-1]
        bins.lim <- bins.lim[-1]
        result <- result[-1, , , drop = FALSE]
      }
      if (keep.NA) 
        result <- list(u = uvec, v = result[, 1, ], n = result[, 
                                                               2, 1], sd = result[, 3, ], bins.lim = bins.lim, 
                       ind.bin = indp)
      else result <- list(u = uvec[indp], v = result[indp, 
                                                     1, ], n = result[indp, 2, 1], sd = result[indp, 3, 
                                                     ], bins.lim = bins.lim, ind.bin = indp)
    }
  }
  
  
  #BD now Non-Euclidean should run just normally in here... just adding semivariance values for all comparisons.
  
  else {
    data <- as.matrix(data)
    v <- matrix(0, nrow = length(u), ncol = n.datasets)
    for (i in 1:n.datasets) {
      v[, i] <- as.vector(dist(data[, i])) #BD taking absolute difference between data points
      if (estimator.type == "modulus") 
        v[, i] <- v[, i, drop = FALSE]^(0.5) #BD this is calculating just the first part of the estimator, in hidden function '.rfm.bin' the rest of the equation is carried out.
      else v[, i] <- (v[, i, drop = FALSE]^2)/2 
    }
    if (!missing(max.dist)) {
      v <- v[u <= max.dist, , drop = FALSE]
      if (direction != "omnidirectional") 
        u.ang <- u.ang[u <= max.dist]
      u <- u[u <= max.dist]
    }
    if (direction != "omnidirectional") {
      ang.lower <- ang.rad - tol.rad
      ang.upper <- ang.rad + tol.rad
      if (ang.lower >= 0 & ang.upper < pi) 
        ang.ind <- (!is.na(u.ang) & ((u.ang >= ang.lower) & 
                                       (u.ang <= ang.upper)))
      if (ang.lower < 0) 
        ang.ind <- (!is.na(u.ang) & ((u.ang < ang.upper) | 
                                       (u.ang > (pi + ang.lower))))
      if (ang.upper >= pi) 
        ang.ind <- (!is.na(u.ang) & ((u.ang > ang.lower) | 
                                       (u.ang < (ang.upper - pi))))
      v <- v[ang.ind, , drop = FALSE]
      u <- u[ang.ind]
    }
    data <- drop(data)
    v <- drop(v)
    if (option == "cloud") {
      result <- list(u = u, v = v)
      if (angles) 
        result$angles <- u.ang
    }
    if (option == "bin") {
      if (missing(max.dist)) 
        umax <- max(u, na.rm=T)
      else umax <- max(u[u < max.dist], na.rm=T)
      if (bin.cloud == "diff"){    #BD Is this code obsolete? just in case
        if (NonEuc == TRUE){  #BD
          stop("Non-Euclidean distances cannot be used to #BD
                                             calculate differences in bins at this time")   #BD
        }                                   #BD
        else dd <- diffpairs(coords, data)$diff #BD this probably doesn't matter... but could give wrong variogram 'v' values since 'diffpars' is calculating Euclidean distance
      }                   #BD
      else dd <- 0
      result <- geoR:::.rfm.bin(cloud = list(u = u, v = v, d = dd),  #BD calls out to a hidden function, but nothing needs to be changed there, just averaging values over bins, still based on u matrix... added ':::' to call from the package namespace
                                estimator.type = estimator.type, uvec = uvec, 
                                breaks = breaks, nugget.tolerance = nugget.tolerance, 
                                bin.cloud = bin.cloud, max.dist = umax, keep.NA = keep.NA)
      if (keep.NA) {
        if (pairs.min > 0) {
          indp <- (result$n < pairs.min)
          if (!nt.ind) {   #BD nugget.tolerance indicator
            for (i in 1:5) result[[i]] <- result[[i]][-1]
            indp <- indp[-1]
          }
          if (is.matrix(result$v)) {
            result$v[indp, ] <- result$sd[indp, ] <- NA
          }
          else {
            result$v[indp] <- result$sd[indp] <- NA
          }
        }
        result$ind.bin <- indp
      }
      else {
        if (pairs.min > 0) {
          if (!nt.ind) {
            for (i in 1:5) result[[i]] <- result[[i]][-1]
          }
          indp <- (result$n >= pairs.min)
          if (is.matrix(result$v)) {
            result$v <- result$v[indp, ]
            result$sd <- result$sd[indp, ]
          }
          else {
            result$v <- result$v[indp]
            result$sd <- result$sd[indp]
          }
          result$u <- result$u[indp]
          result$n <- result$n[indp]
        }
        result$ind.bin <- indp
      }
    }
    if (option == "smooth") {
      if (is.matrix(v)) 
        stop("smooth not yet available for more than one data-set")
      temp <- ksmooth(u, v, ...)
      result <- list(u = temp[[1]], v = temp[[2]])
    }
    if (missing(max.dist)) 
      max.dist <- max(u)
  }
  if (nt.ind) {
    if (!exists(".variog4.nomessage", where = 1)) 
      cat("variog: co-locatted data found, adding one bin at the origin\n")
    if (all(result$u[1:2] < 1e-11)) 
      result$u[2] <- sum(result$bins.lim[2:3])/2
  }
  result <- c(result, list(var.mark = data.var, beta.ols = beta.ols, 
                           output.type = option, max.dist = max.dist, estimator.type = estimator.type, 
                           n.data = n.data, lambda = lambda, trend = trend, pairs.min = pairs.min))
  result$nugget.tolerance <- nugget.tolerance
  if (direction != "omnidirectional") 
    result$direction <- ang.rad
  else result$direction <- "omnidirectional"
  if (direction != "omnidirectional") 
    result$tolerance <- tol.rad
  else result$tolerance <- "none"
  result$uvec <- uvec
  result$call <- call.fc
  oldClass(result) <- "variogram"
  return(result)
}

















Embed<-function(dmat){
  n<-dim(dmat)[1]  #takes row count of distance matrix
  J<-matrix(1,n,n) #makes a nxn matrix of ones
  H<-diag(n)-(1/n)*J #turns 100 into a diagnoal matrix of 1s (100x100) and subtracts it by J divided by 100 (a full matrix os 0.01) ... Hat Matrix?
  A<-dmat^2 #squares values of the distance matrix
  B<-(-1/2)*H%*%A%*%H #Matrix multiplication divided by -0.5  ... not sure why, is this Euclidian calculation? or is this multidimensional scaling
  return(eigen(B)) #returns the 100 eigen values... and vectors?
} 














#added function for wave
#added an option to determine how many eigenvalues can be used for distance of sqrt(gamma)


#semivarioNE <- function(dist, vario, model="exponential", MDSgamma=FALSE, eig=kth, ith=1){ 
semivarioNE <- function(dist, vario, model="exponential", MDSgamma=FALSE, ith=1, lam=kth){   #sample, cv,
  if (class(dist) != "dist")
    stop ("distance object must be of class dist")      
  else if (all(class(vario)!="variofit"))
    stop("covariance parameters must be from a variofit object")     
  else {
    distm <- as.matrix(dist)
    gamma <- matrix(nrow=dim(distm)[1], ncol=dim(distm)[2], 
                    data=rep(0,dim(distm)[1]*dim(distm)[2]))
  } 
  if (model =="spherical"){                
    stop ("Multidimensional scaling of Non-Euclidean distances 
                      not valid with spherical semivariogram model")       
  }
  
  else if (model=="linear"){
    gamma <- vario[[1]] + vario[[2]][1]*distm  
  }
  else if (model=="gaussian"){
    gamma <- vario[[1]] + vario[[2]][1]*(1-exp(-((distm/vario[[2]][2])^2))) 
  }     
  else if (model=="power"){
    gamma <- vario[[1]] + vario[[2]][1]*(distm^vario[[2]][2])   
  }   
  else if (model=="matern"){
    gamma_pre <- ((distm/vario[[2]][2]) ^ vario[[4]]) * (besselK((distm/vario[[2]][2]), nu=vario[[4]]))  
    gamma_pre[gamma_pre==Inf] <- .Machine$double.xmax
    gamma <- vario[[1]] + vario[[2]][1] *
      (1 - (1/(2^(vario[[4]]-1)*gamma(vario[[4]]) ) )* gamma_pre )  
  }  
  else if (model=="exponential") {  
    gamma<- vario[[1]] + vario[[2]][1]*
      (1-exp(-distm/vario[[2]][2])) 
  }
  else if (model=="powered.exponential"){
    gamma <- vario[[1]] + vario[[2]][1] * 
      (1-exp(-1*(distm/vario[[2]][[2]])^vario[[4]]))
  }
  else if (model=="cauchy"){
    gamma <- vario[[1]] + vario[[2]][1]*
      (1- (1+(distm/vario[[2]][[2]])^2)^-vario[[4]] )
  }
  else stop ("Unknown semivariogram model")    
  
  diag(gamma) <- 0
  gamma[gamma<0] <- vario[[1]]   #to fix for matern!
  
  if (MDSgamma==TRUE) {
    cat("MDS on semivariogram ")
    
    NR <- nrow(as.matrix(gamma))
    gamma.5 <- sqrt(gamma)
    MDSgamma.5 <- cmdscale(gamma.5, k=NR-1, eig=T) #remember will only calculate points (vectors) connected to positive eigenvalues
    #MDSgamma.5 <- Embed(gamma.5)
    
    ###Old kth dimenion
    #for (i in 1:(NR-1)){    
    #        #if (round(MDSgamma.5$values[i],5) <= 0) {#      
    #        if (round(MDSgamma.5$eig[i],5) <= 0) {
    #                kth <- i-1
    #                break
    #        }
    #        else kth <- NR-1
    #}
    
    ### finding ideal dimension to estimate to when log10(eigen[i] - eigen[i+1]) < -1.5 
    kpre <- NULL
    for (i in 1:(NR-1)){ 
      kpre[i] <- MDSgamma.5$eig[i] - MDSgamma.5$eig[i+1]
    }
    for(i in 1:length(kpre)){
      if (log10(kpre[i]) < -1.5){
        kth <- i
        break
      }
    }
    
    kth
    #gamma.5New<-dist(MDSgamma.5$vectors[,1:kth] %*% diag(sqrt(MDSgamma.5$values[1:kth])))
    #gamma.5New<-dist(MDSgamma.5$points[,1:eig])
    #gamma.5New<-dist(MDSgamma.5$points[,1:10]) #reducing the number of dimensions used
    gamma.5New<-dist(MDSgamma.5$points[,1:lam])  # creating the diagnostic function across all possible dimensions
    #gammaNew<-(gamma.5New)^2
    #### new version since raising to the power of 2 is not a valid semivariogram
    ## Power version
    #gammaNew<-(gamma.5New)^1.9
    
    
    ## Gaussian version
    ### trying to fit it as well
    VarioGam <- list(u=as.vector(as.dist(gamma.5New)),
                     
                     ##for v6JHMI
                     #v=as.vector(as.dist(ith[[1]])), 
                     #var.mark=max(ith[[1]]), beta.ols=NA,
                     #For v6
                     #v=as.vector(as.dist(semiVNE[[ith]][[1]])), 
                     #   var.mark=max(semiVNE[[ith]][[1]]), beta.ols=NA, 
                     v=as.vector(as.dist(semiVNE[[1]])), 
                     var.mark=max(semiVNE[[1]]), beta.ols=NA, 
                     output.type="cloud",
                     max.dist=max(gamma.5New), estimator.type="classical", 
                     n.data=300,
                     lambda=1, trend="cte", pairs.min=2, nugget.tolerance=1e-12, 
                     direction="omnidirectional", tolerance="none", uvec="default",
                     call=NA)
    class(VarioGam) <- "variogram"
    #vtest <- variofit(VarioGam, cov.model="gaussian")
    vtest <- variofit(VarioGam, ini.cov.pars=c(100, 10), cov.model="gaussian")
    
    gammaNew<- vtest[[1]] + vtest[[2]][1]*(1-exp(-((gamma.5New/vtest[[2]][2])^2))) 
    
    #gammaNew<- 0 + 4560*(1-exp(-((gamma.5New/68)^2)))
    #gammaNew<-(gamma.5New)^2
    
    
    
    
    
    #this is needed if we included non-sampled locations
    #gamma_samp <- as.matrix(gammaNew)[sample,sample]
    
    #gamma_train<-as.matrix(gamma_samp)[-cv,-cv]  
    #gamma_predict<-as.matrix(gamma_samp)[-cv, cv, drop=FALSE]
    
    #semiV <- list(as.matrix(gammaNew),kth,round(MDSgamma.5$values,5))
    semiV <- list(as.matrix(gammaNew),kth,round(MDSgamma.5$eig,5), MDSgamma.5$points)
  }
  else{
    cat("No MDS performed")
    
    #gamma
    #gamma_samp <- gamma[sample,sample]
    #gamma_train<-as.matrix(gamma_samp)[-cv,-cv]  
    #gamma_predict<-as.matrix(gamma_samp)[-cv, cv, drop=FALSE]
    
    semiV <- list(gamma)
  }
  
  class(semiV) <- "semivarianceNE"
  return(semiV)
}































#included new option for subsetting by rownames.

#now semivariance must just be a matrix...! No longer calling out the object...

#will write the function so if missing semivariance and missing covariance matrix stop
#Also stop if both are filled in, can only have one. 

krigeNE <- function (geodata, coords = geodata$coords, data = geodata$data, 
                     locations, krige, output, semivariance, covariance, cv, cv.type="subset", 
                     tolerance=.Machine$double.eps) {   
  #now cv will subset from this matrix.
  #cv.type will determine how to subset...
  
  #BD ... probably don't need borders... this means I need to determine prediction locations and their distances BEFORE I run this function. I should also know what points I want (created from borders)
  
  #BD also distance epsilon will not be used, most distances will be far enough not to matter though.
  
  if (missing(geodata)) 
    stop ("geodata object missing")                                         #BD
  #BD removed krige1D lines of code
  call.fc <- match.call()
  base.env <- sys.frame(sys.nframe())
  if (missing(krige)) 
    krige <- krige.control()
  else {
    if (length(class(krige)) == 0 || class(krige) != "krige.geoR") {
      if (!is.list(krige)) 
        stop("krigeNE: the argument krige only takes a list or an output of the function krige.control")
      else {
        krige.names <- c("type.krige", "trend.d", "trend.l", 
                         "obj.model", "beta", "cov.model", "cov.pars", 
                         "kappa", "nugget", "micro.scale", "dist.epsilon", 
                         "lambda", "aniso.pars")
        krige.user <- krige
        krige <- list()                 
        if (length(krige.user) > 0) {
          for (i in 1:length(krige.user)) {
            n.match <- match.arg(names(krige.user)[i], 
                                 krige.names)
            krige[[n.match]] <- krige.user[[i]]       
          }
        }
        if (is.null(krige$type.krige)) 
          krige$type.krige <- "ok"
        if (is.null(krige$trend.d)) 
          krige$trend.d <- "cte"
        if (is.null(krige$trend.l)) 
          krige$trend.l <- "cte"
        if (is.null(krige$obj.model)) 
          krige$obj.model <- NULL
        if (is.null(krige$beta)) 
          krige$beta <- NULL
        if (is.null(krige$cov.model)) 
          krige$cov.model <- "matern"
        if (is.null(krige$cov.pars)) 
          stop("covariance parameters (sigmasq and phi) should be provided in cov.pars")
        if (is.null(krige$kappa)) 
          krige$kappa <- 0.5
        if (is.null(krige$nugget)) 
          krige$nugget <- 0
        if (is.null(krige$micro.scale)) 
          krige$micro.scale <- 0
        if (is.null(krige$dist.epsilon)) 
          krige$dist.epsilon <- 1e-10
        if (is.null(krige$aniso.pars)) 
          krige$aniso.pars <- NULL
        if (is.null(krige$lambda)) 
          krige$lambda <- 1
        krige <- krige.control(type.krige = krige$type.krige, 
                               trend.d = krige$trend.d, trend.l = krige$trend.l, 
                               obj.model = krige$obj.model, beta = krige$beta, 
                               cov.model = krige$cov.model, cov.pars = krige$cov.pars, 
                               kappa = krige$kappa, nugget = krige$nugget, 
                               micro.scale = krige$micro.scale, dist.epsilon = krige$dist.epsilon, 
                               aniso.pars = krige$aniso.pars, lambda = krige$lambda)
      }
    }
  }                                                                  #BD I don't think I need to remove the aniso.pars, we're probably never going to use it anyway.
  cov.model <- krige$cov.model
  kappa <- krige$kappa
  lambda <- krige$lambda
  beta <- krige$beta
  cov.pars <- krige$cov.pars
  nugget <- krige$nugget
  micro.scale <- krige$micro.scale
  aniso.pars <- krige$aniso.pars
  if (missing(output)) 
    output <- output.control()
  else {
    if (length(class(krige)) == 0 || class(output) != "output.geoR") {
      if (!is.list(output)) 
        stop("krigeNE: the argument output can take only a list or an output of the function output.control")
      else {
        output.names <- c("n.posterior", "n.predictive", 
                          "moments", "n.back.moments", "simulations.predictive", 
                          "mean.var", "quantile", "threshold", "signal", 
                          "messages.screen")                            
        output.user <- output
        output <- list()
        if (length(output.user) > 0) {
          for (i in 1:length(output.user)) {
            n.match <- match.arg(names(output.user)[i], 
                                 output.names)
            output[[n.match]] <- output.user[[i]]
          }
        }
        if (is.null(output$n.posterior)) 
          output$n.posterior <- 1000
        if (is.null(output$n.predictive)) 
          output$n.predictive <- NULL
        if (is.null(output$moments)) 
          output$moments <- TRUE
        if (is.null(output$n.back.moments)) 
          output$n.back.moments <- 1000
        if (is.null(output$simulations.predictive)) {
          if (is.null(output$n.predictive)) 
            output$simulations.predictive <- NULL
          else output$simulations.predictive <- ifelse(output$n.predictive > 
                                                         0, TRUE, FALSE)
        }
        if (is.null(output$mean.var)) 
          output$mean.var <- NULL
        if (is.null(output$quantile)) 
          output$quantile <- NULL
        if (is.null(output$threshold)) 
          output$threshold <- NULL
        if (is.null(output$sim.means)) 
          output$sim.means <- NULL
        if (is.null(output$sim.vars)) 
          output$sim.vars <- NULL
        if (is.null(output$signal)) 
          output$signal <- NULL
        if (is.null(output$messages.screen)) 
          output$messages.screen <- TRUE
        output <- output.control(n.posterior = output$n.posterior, 
                                 n.predictive = output$n.predictive, moments = output$moments, 
                                 n.back.moments = output$n.back.moments, simulations.predictive = output$simulations.predictive, 
                                 mean.var = output$mean.var, quantile = output$quantile, 
                                 threshold = output$threshold, sim.means = output$sim.means, 
                                 sim.vars = output$sim.vars, signal = output$signal, 
                                 messages = output$messages.screen)
      }
    }
  }
  signal <- ifelse(is.null(output$signal), FALSE, output$signal)
  messages.screen <- output$messages.screen
  n.predictive <- output$n.predictive
  n.back.moments <- output$n.back.moments
  n.predictive <- ifelse(is.null(n.predictive), 0, n.predictive)
  simulations.predictive <- ifelse(is.null(output$simulations.predictive), 
                                   FALSE, TRUE)
  mean.estimator <- output$mean.estimator
  sim.means <- output$sim.means
  if (is.null(sim.means)) 
    sim.means <- ifelse(simulations.predictive, TRUE, FALSE)
  sim.vars <- output$sim.vars
  if (is.null(sim.vars)) 
    sim.vars <- FALSE
  if (is.null(mean.estimator) & simulations.predictive) 
    mean.estimator <- TRUE
  quantile.estimator <- output$quantile.estimator
  probability.estimator <- output$probability.estimator
  if (!is.null(probability.estimator)) {
    if (length(probability.estimator) > 1 & length(probability.estimator) != 
        nrow(locations)) 
      stop("krigeNE: probability.estimator must either have length 1, or have length = nrow(locations)\n")
  }
  if (simulations.predictive & n.predictive == 0) #BD simulations.predictive is set to false, so n.predictive remains at zero...               
    stop("cannot currently run simulations in non-Euclidean version of kriging")    #BD will prevent simulations for now...
  
  
  
  if (krige$type.krige == "ok") 
    beta.prior <- "flat"
  if (krige$type.krige == "sk") 
    stop ("KrigeNE only performs ordinary kriging")                     #BD
  
  #BD removed all code related to coords and locations
  
  if (krige$trend.d != "cte" || krige$trend.l !="cte")                         #BD
    stop("Non-Euclidean Kriging only works with a constant mean for trend (Ordinary Kriging)") #BD
  
  if (messages.screen) {                                                      #BD removed non-constant trend code
    cat(as.character(krige$trend.d)[1], cte = "krigeNE: model with constant mean")
    cat("\n")
  }
  if (class(krige$trend.d) == "trend.spatial") 
    trend.d <- unclass(krige$trend.d)
  else trend.d <- unclass(trend.spatial(trend = krige$trend.d, 
                                        geodata = geodata))                                     #BD this will just create a matrix (one column) of 1s if using "cte"
  if (nrow(trend.d) != nrow(geodata$coords))                          #replaced coords object
    stop("coords and trend.d have incompatible sizes")
  beta.size <- ncol(trend.d)
  if (class(krige$trend.l) == "trend.spatial") 
    trend.l <- unclass(krige$trend.l)
  else trend.l <- unclass(trend.spatial(trend = krige$trend.l, 
                                        geodata = list(coords = locations)))                         
  if (beta.size > 1) 
    beta.names <- paste("beta", (0:(beta.size - 1)), sep = "")
  else beta.names <- "beta"
  
  #BD removed anisotropy code
  
  if (!isTRUE(all.equal(lambda, 1))) {
    if (messages.screen) 
      cat("krigeNE: performing the Box-Cox data transformation\n")
    data <- BCtransform(x = data, lambda = lambda)$data
    warning ("krigeNE not formally tested for box-cox transformations")     #BD just to remind myself I need to make sure my new code works the same for Box-Cox transformations, it may very well not, especially if Box-Cox changes the lambda weights.
  }
  if (is.vector(cov.pars)) {
    sigmasq <- cov.pars[1]
    phi <- cov.pars[2]
    cpars <- c(1, phi)
  }
  else {
    stop("current version of krigeNE does not accept nested covariance models\n")   #BD
  }
  sill.partial <- sum(sigmasq)                                                #BD But why is the sum of an element being taken here???
  if (sill.partial < 1e-16) {
    tausq.rel <- 0
    tausq.rel.micro <- 0
  }
  else {
    tausq.rel <- nugget/sum(sigmasq)                                        #BD I guess relative nugget?
    tausq.rel.micro <- micro.scale/sum(sigmasq)
  }
  n <- length(data)
  ni <- nrow(trend.l)
  kc <- list()
  
  nug.factor <- ifelse(signal, tausq.rel.micro, tausq.rel)
  
  ### Kriging section, variance/covariance input here
  
  if(missing(semivariance) & missing(covariance)){
    stop(" A semivariance or covariance matrix must be supplied")
  }
  else if(!missing(semivariance) & !missing(covariance)){
    stop("Only one semivariance or covariance matrix, not both, must be supplied")
  }
  
  else if(missing(covariance)){
    cat("krigeNE: performing Kriging with semivariance matrix \n")
    
    ##### beforehand I need to have calculated initial semivariane matrix from WLS parameter estimates, 
    ### re-estimated using MDS on square-root of matrix elements
    ### Maybe I can create a function that includes these as a list.... 
    ### let's say I've done that, first part of the list is cap-Gamma and second part is small gamma.
    
    # if (class(semivariance) != "semivarianceNE"){
    #                warning("KrigeNE: object semivariance must be output from function 'semivarianceNE'")
    #       }
    
    
    
    
    if(cv.type=="subset"){
      SVtrain <- semivariance[-cv,-cv]                   #square matrix of sampled locations (observed X bserved)
      SVpredict <-semivariance[-cv,cv,drop=FALSE]        #non-square matrix of unsampled locations to be predicted at (observed X prediction)
    }
    
    else if (cv.type=="row_names"){
      SVtrain <- semivariance[!rownames(semivariance[[1]]) %in% 
                                paste("V",cv, sep=""), 
                              !rownames(semivariance[[1]]) %in% 
                                paste("V",cv, sep="") ]  #square matrix of sampled locations (observed X bserved)
      SVpredict <-semivariance[!rownames(semivariance[[1]]) %in% 
                                 paste("V",cv, sep=""),
                               paste("V",cv, sep=""),
                               drop=FALSE]        #non-square matrix of unsampled locations to be predicted at (observed X prediction)
    }
    
    else if (cv.type=="row_names_euc"){
      SVtrain <- semivariance[!rownames(semivariance[[1]]) %in% cv, 
                              !rownames(semivariance[[1]]) %in% cv]  #square matrix of sampled locations (observed X bserved)
      SVpredict <-semivariance[!rownames(semivariance[[1]]) %in% cv,
                               cv, drop=FALSE]        #non-square matrix of unsampled locations to be predicted at (observed X prediction)
    }
    
    
    Capg <- cbind(rbind(SVtrain, 1), 
                  c(rep(1,nrow(SVtrain)),0))
    remove(SVtrain)
    s0g <- rbind(SVpredict,1)
    #    v0 <- loccoords(coords = coords, locations = locations)     #BD this computes Euclidean distances for the 's0g' matrix above (semivariances between sampled and unsampled locations).
    remove(SVpredict)
    lambda0 <- crossprod(solve(Capg, tol=tolerance), s0g)
    remove(Capg)
    pred<-colSums(data * lambda0[-nrow(lambda0),,drop=FALSE])  
    kcvar <- diag(crossprod(lambda0,s0g))
    remove(lambda0, s0g)
    
    #BD
  }                                                                   #BD
  
  
  else{
    cat("krigeNE: performing Kriging with covariance matrix \n")
    
    if(cv.type=="subset"){
      Vcov <- as.matrix(covariance)[-cv,-cv] #putting in new covariance matrix of measured points (train) from Non-Euclidean distance instead
      v0 <- covariance[-cv,cv, drop=FALSE]
    }
    
    else if (cv.type=="row_names"){
      Vcov <- covariance[!rownames(covariance) %in% 
                           paste("V",cv, sep=""), 
                         !rownames(covariance) %in% 
                           paste("V",cv, sep="") ]  #square matrix of sampled locations (observed X observed)
      v0 <-covariance[!rownames(covariance) %in% 
                        paste("V",cv, sep=""),
                      paste("V",cv, sep=""),
                      drop=FALSE]        #non-square matrix of unsampled locations to be predicted at (observed X prediction)
    }
    
    else if (cv.type=="row_names_euc"){
      Vcov <- covariance[!rownames(covariance) %in% cv, 
                         !rownames(covariance) %in% cv]  #square matrix of sampled locations (observed X observed)
      v0 <-covariance[!rownames(covariance) %in% cv,
                      cv, drop=FALSE]        #non-square matrix of unsampled locations to be predicted at (observed X prediction)
    }
    
    
    
    CapSig <- cbind(rbind(Vcov, 1), 
                    c(rep(1,nrow(Vcov)),0))
    
    c_u <- rbind(v0,1)
    
    lambda0 <- crossprod(solve(CapSig, tol=tolerance), c_u)
    remove(CapSig, c_u)
    
    pred<-colSums(data * lambda0[-nrow(lambda0),,drop=FALSE])
    
    ###variance setup
    
    #faster way to calcualte third term of variance function, 
    #only works if predicting at one location
    
    if (ncol(lambda0)==1){
      lambdai <- matrix(,nrow=nrow(Vcov),ncol=ncol(Vcov))
      lambdaj <- matrix(,nrow=nrow(Vcov),ncol=ncol(Vcov))
      for(g in 1:nrow(Vcov)){
        lambdai[g,] <-lambda0[-nrow(lambda0),,drop=FALSE]
      }
      for(g in 1:ncol(Vcov)){
        lambdaj[,g] <-lambda0[-nrow(lambda0),,drop=FALSE]
      }
      
      var3rdterm <- sum(lambdai*lambdaj*Vcov)
      
      remove(lambdai, lambdaj)
      
    }
    
    #Longer way if more than one prediction location
    
    else{
      var3rdterm   <-     rowSums( mapply( function(z){
        rowSums(mapply(function(y) {
          return(lambda0[z,]*lambda0[y,]*Vcov[y,z])
        },  y=1:nrow(Vcov)))
      }, z=1:ncol(Vcov)) )
    }
    
    
    kcvar <- max(as.matrix(covariance))  - 
      2*colSums(v0* lambda0[-nrow(lambda0),,drop=FALSE]) + 
      var3rdterm
    
    #### old 3rd term of variance function       
    #        sum( mapply( function(z){
    #                sum(mapply(function(y) {
    #                        return(lambda0[z]*lambda0[y]*Vcov[z,y])
    #                },  y=1:ncol(Vcov)))
    #        }, z=1:nrow(Vcov)) )
    
    
    remove(lambda0, var3rdterm, Vcov, v0)
    
    
  }
  
  kc$predict <- pred
  kc$krige.var <- kcvar
  
  kc$beta.est <- NA  
  
  remove(pred, kcvar)
  
  if (messages.screen) {                                           #BD 
    cat("krigeNE Reminder: betas cannot currently be estimated in this function")   #BD
    cat("\n")     
  }
  
  
  ### BD Old Code below, will keep here as reference as  I may want to go back and use this once covariance matrices are calculated. ###
  #           Vcov <- varcov.spatial(coords = coords, cov.model = cov.model,
  #                                  kappa = kappa, nugget = tausq.rel, 
  #                                  cov.pars = cpars)$varcov     #BD This is where we input the variance/covariance matrix
  
  
  # Vcov <- as.matrix(covariance)[-cv,-cv] #putting in new covariance matrix of measured points (train) from Non-Euclidean distance instead     
  
  #ivtt <- solve(Vcov, trend.d)                #BD Comes out as a one column matrix of values multiplying by a column of ones?
  #ttivtt <- crossprod(ivtt, trend.d)          #and this leads to a single value... have to find out from Frank what is going on here...
  #beta.flat <- drop(solve(ttivtt, crossprod(ivtt, as.vector(data))))   #BD still left with one value, this is beta-hat, cannot currently be estimated using semivariances
  #remove("ivtt")
  
  
  
  #BD removed 'n.predictive' code (coinciding data)
  
  ### BD more cold code I may want to keep when I bring back covariance matrices ###
  #    v0 <- loccoords(coords = coords, locations = locations)     #BD this computes Euclidean distances for the 's0g' matrix above (semivariances between sampled and unsampled locations).
  #    ind.v0 <- which(v0 < krige$dist.epsilon) #BD takes out those small values as a vector
  #    v0 <- cov.spatial(obj = v0, cov.model = cov.model, kappa = kappa, 
  #        cov.pars = cpars)  #BD Now makes a covariance matrix based on those distances....
  #    v0[ind.v0] <- 1 + nug.factor        #BD just give coincidental data the nugget + 1
  
  
  #SVpredict <-semivariance[-cv,cv,drop=FALSE]        #non-square matrix of unsampled locations to be predicted at (observed X prediction)
  #v0 <- covariance[-cv,cv, drop=FALSE]        #BD again just putting in covarinace matrix instead
  
  #ivv0 <- solve(Vcov, v0)             #BD This is where the real math happens... but am I doing this right?
  #tv0ivv0 <- colSums(v0 * ivv0)
  #remove("Vcov") #, "ind.v0")
  #b <- crossprod(cbind(data, trend.d), ivv0) #ends up as a 2x1 matrix
  #if (n.predictive == 0) 
  #remove("v0",  "ivv0")
  #tv0ivdata <- drop(b[1, ]) #takes first row of b
  #b2 <- t(trend.l) - b[-1, , drop = FALSE] #renamed so I can look at b...
  #takes second row of b (or all rows but the first)
  
  #BD removed simple kriging code
  
  # pred<-colSums(data * lambda0[-nrow(lambda0),,drop=FALSE])  
  #kc$predict <- tv0ivdata + drop(crossprod(b2, beta.flat))  
  #bitb <- colSums(b2 * solve(ttivtt, b2))
  #kc$krige.var <- sill.partial * drop(1 + nug.factor - tv0ivv0 + bitb)   #maybe this isn't right either...
  
  #kc$krige.var <- tv0ivv0/35623.12 
  
  #kc$krige.var <-   36767.96 * 
  #        drop(1 + 0.04941725 - tv0ivv0/36767.96 + bitb/36767.96)
  
  #nugget?   10.0434
  
  
  #kc$krige.var <-  (max(as.matrix(covariance)) - min(as.matrix(covariance)) ) * 
  #        drop(1 + 0 - ((tv0ivv0 + bitb) /max(as.matrix(covariance))) )
  
  #(max(as.matrix(covariance)) - min(as.matrix(covariance)) ) * 
  #        drop(1 + min(as.matrix(covariance)) - tv0ivv0 + bitb) 
  #(solve(ttivtt, crossprod(ivtt, as.vector(data))))
  #crossprod(ivv0, v0)
  #solve(ttivtt, crossprod(ivtt, as.vector(data))) %*% t(v0)
  
  #solve(ttivtt, crossprod(ivtt, as.vector(data)))
  
  
  #kc$krige.var <-  (max(as.matrix(covariance)) - min(as.matrix(covariance)) )
  
  #BD moved kc$beta.est
  # remove("bitb", "b", "tv0ivv0", "tv0ivdata")         #BD slight alteration/simplifcation
  #        kc$beta.est <- beta.flat              #BD used to be 'beta.flat' object, but can bring that back later once covariance matrices are used
  #       names(kc$beta.est) <- beta.names
  
  
  
  kc$distribution <- "normal"
  #if (any(kc$krige.var < 0)) 
  #        cat("krige.conv: negative kriging variance found! Investigate why this is happening.\n")                                           #BD I think this should happen first before imputing zeros...
  #kc$krige.var[kc$krige.var < 1e-08] <- 0
  #BD removed n.predictive code for simulations.
  if (!isTRUE(all.equal(lambda, 1))) {
    if (messages.screen) {
      cat("krigeNE: back-transforming the predicted mean and variance\n")
      if (!isTRUE(all.equal(lambda, 0)) & !isTRUE(all.equal(lambda, 
                                                            0.5))) 
        cat("krigeNE: back-transforming by simulating from the predictive.\n           (run the function a few times and check stability of the results.\n")
    }
    kc[c("predict", "krige.var")] <- backtransform.moments(lambda = lambda, 
                                                           mean = kc$predict, variance = kc$krige.var, distribution = "normal", 
                                                           n.simul = n.back.moments)[c("mean", "variance")]
  }
  message <- "krigeNE: Kriging performed using global neighbourhood"
  if (messages.screen) 
    cat(paste(message, "\n"))
  kc$message <- message
  kc$call <- call.fc
  attr(kc, "sp.dim") <- "2d"                  #BD just 2d option now
  attr(kc, "prediction.locations") <- call.fc$locations
  attr(kc, "parent.env") <- parent.frame()
  if (!is.null(call.fc$coords)) 
    attr(kc, "data.locations") <- call.fc$coords
  else attr(kc, "data.locations") <- substitute(a$coords, list(a = substitute(geodata)))
  #   if (!is.null(call.fc$borders)) 
  #        attr(kc, "borders") <- call.fc$borders
  oldClass(kc) <- "kriging"
  return(kc)
}













































































krigeCV <- function(train.geo, pred, semiV, predloc, CVtype, CV, KrigeC, Messages=TRUE){
  if (class(train.geo) !="geodata")
    stop ("train.geo data must be class 'geodata'")
  if(length(train.geo$data) != dim(predloc)[1])
    stop ("train.geo and predloc must be have same number of slots")
  
  if(Messages==FALSE)
    OC <- output.control(messages=F)
  else OC <- output.control()
  
  
  
  
  if (CVtype == "row_names"){
    KrigeEuc <- krigeNE(geodata=train.geo, locations=predloc, 
                        semivariance=semiV[[1]], cv.type="row_names_euc",  #keep cv.type and cv set since they don't have V in front
                        cv=CV, krige=KrigeC, output=OC)  
  }
  
  else{
    KrigeEuc <- krigeNE(geodata=train.geo, locations=predloc, 
                        semivariance=semiV[[1]], cv.type="CVtype",
                        cv=CV, krige=KrigeC, output=OC)  
  }
  
  
  KrigeMDSdist <- krigeNE(geodata=train.geo, locations=predloc, 
                          semivariance=semiV[[2]], cv.type=CVtype, 
                          cv=CV, krige=KrigeC, output=OC)
  
  KrigeNEdist <- krigeNE(geodata=train.geo, locations=predloc, 
                         semivariance=semiV[[3]], cv.type=CVtype, 
                         cv=CV, krige=KrigeC, output=OC)
  
  KrigeMDSgamma <- krigeNE(geodata=train.geo, locations=predloc, 
                           semivariance=semiV[[4]], cv.type=CVtype, 
                           cv=CV, krige=KrigeC, OC)
  
  
  #and here goes cross-validation code...
  ### don't forget afterwards to remove kriging objects since they are so large!
  
  ### Changing CI here to a vector...
  
  EucSE <- sqrt(KrigeEuc$krige.var)
  EucRMSE <- sqrt(mean((KrigeEuc$pred - pred$idata)^2))
  EucRMSEstd <- sqrt(mean(((KrigeEuc$pred - pred$idata)/EucSE)^2))
  EucCI95 <- ifelse(pred$idata >= KrigeEuc$pred - 
                      1.96*sqrt(KrigeEuc$krige.var) &
                      pred$idata <= KrigeEuc$pred + 
                      1.96*sqrt(KrigeEuc$krige.var), 1, 0)
  EucUncert <- mean(EucSE[is.finite(EucSE)])
  EucNegVar <- any(KrigeEuc$krige.var<0)
  remove(KrigeEuc, EucSE)
  
  MDSdistSE <- sqrt(KrigeMDSdist$krige.var)
  MDSdistRMSE <- sqrt(mean((KrigeMDSdist$pred - pred$idata)^2))
  MDSdistRMSEstd <- sqrt(mean(((KrigeMDSdist$pred - pred$idata)/MDSdistSE)^2))
  MDSdistCI95 <- ifelse(pred$idata >= KrigeMDSdist$pred - 
                          1.96*sqrt(KrigeMDSdist$krige.var) &
                          pred$idata <= KrigeMDSdist$pred + 
                          1.96*sqrt(KrigeMDSdist$krige.var), 1, 0)
  MDSdistUncert <- mean(MDSdistSE[is.finite(MDSdistSE)])
  MDSdistNegVar <- any(KrigeMDSdist$krige.var<0)
  remove(KrigeMDSdist, MDSdistSE)
  
  NEdistSE <- sqrt(KrigeNEdist$krige.var)
  NEdistRMSE <- sqrt(mean((KrigeNEdist$pred - pred$idata)^2))
  NEdistRMSEstdpre <- ((KrigeNEdist$pred - pred$idata)/NEdistSE)^2
  NEdistRMSEstd <- sqrt(mean(NEdistRMSEstdpre[is.finite(NEdistRMSEstdpre)]))
  NEdistCI95 <- ifelse(pred$idata >= KrigeNEdist$pred - 
                         1.96*sqrt(KrigeNEdist$krige.var) &
                         pred$idata <= KrigeNEdist$pred + 
                         1.96*sqrt(KrigeNEdist$krige.var), 1, 0)
  NEdistUncert <- mean(NEdistSE[is.finite(NEdistSE)])
  NEdistNegVar <- any(KrigeNEdist$krige.var<0)
  remove(KrigeNEdist, NEdistSE)
  
  MDSgammaSE <- sqrt(KrigeMDSgamma$krige.var)
  MDSgammaRMSE <- sqrt(mean((KrigeMDSgamma$pred - pred$idata)^2))
  MDSgammaRMSEstd <- sqrt(mean(((KrigeMDSgamma$pred - pred$idata)/MDSgammaSE)^2))
  MDSgammaCI95 <- ifelse(pred$idata >= KrigeMDSgamma$pred - 
                           1.96*sqrt(KrigeMDSgamma$krige.var) &
                           pred$idata <= KrigeMDSgamma$pred + 
                           1.96*sqrt(KrigeMDSgamma$krige.var), 1, 0)
  MDSgammaUncert <- mean(MDSgammaSE[is.finite(MDSgammaSE)])
  MDSgammaNegVar <- any(KrigeMDSgamma$krige.var<0)
  remove(KrigeMDSgamma, MDSgammaSE)
  
  return(c(EucRMSE,MDSdistRMSE,NEdistRMSE,MDSgammaRMSE, EucRMSEstd, 
           MDSdistRMSEstd,NEdistRMSEstd,MDSgammaRMSEstd, EucUncert, MDSdistUncert,
           NEdistUncert,MDSgammaUncert, EucNegVar,MDSdistNegVar,
           NEdistNegVar,MDSgammaNegVar, EucCI95,
           MDSdistCI95,NEdistCI95,MDSgammaCI95))
}




























grfNE<-function (n, grid = "irreg", nx, ny, xlims = c(0, 1), ylims = c(0, 1), 
                 borders, nsim = 1, cov.model = "matern", 
                 cov.pars = stop("missing covariance parameters sigmasq and phi"), 
                 NEdist = stop("missing Non-Euclidean distances"), #BD
                 kappa = 0.5, nugget = 0, lambda = 1, aniso.pars = NULL, mean = 0, 
                 method, RF = FALSE, messages)       #BD changed RF default to FALSE as I have not corrected that one for NEdistances.
{
  if(class(NEdist)!="dist")                               #BD
    stop("Non-Euclidean Distances must be of class dist") #BD
  if(length(NEdist) != (nrow(grid)^2 - nrow(grid))/2 )    #BD
    stop("dimensions of grid do not match dimensions of NEdistance object")   #BD
  #should make sure distance object matches number of coords...
  call.fc <- match.call()
  if (missing(messages)) 
    messages.screen <- as.logical(ifelse(is.null(getOption("geoR.messages")), 
                                         TRUE, getOption("geoR.messages")))
  else messages.screen <- messages
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    warning(".Random.seed not initialised. Creating it with by calling runif(1)")
    runif(1)
  }
  rseed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  cov.model <- match.arg(cov.model, choices = geoRCovModels)
  if (cov.model == "stable") 
    cov.model <- "powered.exponential"
  if (cov.model == "matern" && kappa == 0.5) 
    cov.model <- "exponential"
  tausq <- nugget
  if (is.vector(cov.pars)) {
    sigmasq <- cov.pars[1]
    phi <- cov.pars[2]
    nst <- 1
  }
  else {
    sigmasq <- cov.pars[, 1]
    phi <- cov.pars[, 2]
    nst <- nrow(cov.pars)
  }
  sill.total <- tausq + sum(sigmasq)
  messa <- geoR:::.grf.aux1(nst, nugget, sigmasq, phi, kappa, cov.model) #BD calling hidden geoR function
  results <- list(coords = NULL, data = NULL)
  if ((!missing(nx) && nx == 1) | (!missing(ny) && ny == 1) | 
      diff(xlims) == 0 | diff(ylims) == 0) {
    sim1d <- TRUE
    if (messages.screen) 
      cat("simulations in 1D\n")
  }
  else sim1d <- FALSE
  if (mode(grid) == "character") {
    grid <- match.arg(grid, choices = c("irreg", "reg"))
    attr(results, "grid") <- grid
    if (!missing(borders) & !sim1d) {
      results$borders <- borders
      if (grid == "irreg") 
        results$coords <- splancs::csr(poly = borders, #BD generating random points in a polygon... don't think I want that for now... maybe later.
                                       npoints = n)
      else {
        if (!missing(nx) && !missing(ny)) {
          bb <- bbox(borders)
          results$coords <- splancs::gridpts(poly = borders, 
                                             xs = diff(bb[1, ])/(nx - 1), ys = diff(bb[1, 
                                             ])/(nx - 1))
        }
        else results$coords <- splancs::gridpts(poly = borders, 
                                                npts = n)
        xgrid <- round(sort(unique(results$coords[, 1])), 
                       digits = 12)
        ygrid <- round(sort(unique(results$coords[, 2])), 
                       digits = 12)
        attr(results, "xgrid") <- c(range(xgrid), unique(diff(xgrid)))
        attr(results, "ygrid") <- c(range(ygrid), unique(diff(ygrid)))
        names(attr(results, "xgrid")) <- c("xmin", "xmax", 
                                           "xstep")
        names(attr(results, "ygrid")) <- c("ymin", "ymax", 
                                           "ystep")
      }
    }
  }
  else {
    results$coords <- as.matrix(grid)
    x1vals <- sort(unique(round(results$coords[, 1], digits = 12)))
    x2vals <- sort(unique(round(results$coords[, 2], digits = 12)))
    if (length(unique(diff(x1vals)) == 1) && length(unique(diff(x2vals)) == 
                                                    1)) 
      attr(results, "grid") <- "reg"
    if (messages.screen) 
      cat("grf: simulation on a set of locations provided by the user\n")
  }
  if (!is.matrix(results$coords) & !is.data.frame(results$coords)) {
    if (missing(nx)) {
      if (sim1d) 
        nx <- ifelse(diff(xlims) == 0, 1, n)
      else nx <- ifelse((mode(grid) == "character" && grid == 
                           "reg"), round(sqrt(n)), n)
    }
    if (missing(ny)) {
      if (sim1d) 
        ny <- ifelse(diff(ylims) == 0, 1, n)
      else ny <- ifelse((mode(grid) == "character" && grid == 
                           "reg"), round(sqrt(n)), n)
    }
    if (mode(grid) == "character" && grid == "irreg") {
      results$coords <- cbind(x = runif(nx, xlims[1], xlims[2]), 
                              y = runif(ny, ylims[1], ylims[2]))
      if (messages.screen) 
        cat(paste("grf: simulation(s) on randomly chosen locations with ", 
                  n, " points\n"))
    }
    else {
      xpts <- seq(xlims[1], xlims[2], length = nx)
      ypts <- seq(ylims[1], ylims[2], length = ny)
      xspacing <- ifelse(length(xpts) == 1, 0, diff(xpts[1:2]))
      yspacing <- ifelse(length(ypts) == 1, 0, diff(ypts[1:2]))
      results$coords <- as.matrix(expand.grid(x = xpts, 
                                              y = ypts))
      equal.spacing <- ifelse(abs(xspacing - yspacing) < 
                                1e-12, TRUE, FALSE)
      if (messages.screen) 
        cat(paste("grf: generating grid ", nx, " * ", 
                  ny, " with ", (nx * ny), " points\n"))
      attr(results, "xgrid") <- c(xmin = xlims[1], xmax = xlims[2], 
                                  xstep = xspacing)
      attr(results, "ygrid") <- c(ymin = ylims[1], ymax = ylims[2], 
                                  ystep = yspacing)
    }
    if (!sim1d & missing(borders)) {
      lbor <- as.matrix(expand.grid(xlims, ylims))
      results$borders <- lbor[chull(lbor), ]
    }
  }
  
  
  n <- nrow(results$coords)
  
  
  if (length(unique(round(results$coords[, 1], digits = 12))) == 
      1 | length(unique(round(results$coords[, 2], digits = 12))) == 
      1) 
    sim1d <- TRUE
  else sim1d <- FALSE
  if (!RF && !is.null(aniso.pars)) {
    if (length(aniso.pars) != 2 | mode(aniso.pars) != "numeric") 
      stop("anisotropy parameters must be provided as a numeric vector with two elements: the rotation angle (in radians) and the anisotropy ratio (a number greater than 1)")
    if (messages.screen) 
      cat("grf: transforming to the isotropic space \n")
    results$coords <- coords.aniso(coords = results$coords, 
                                   aniso.pars = aniso.pars)
  }
  if (missing(method)) {
    method <- "cholesky"
    if (n > 500 && RF) 
      method <- "RF"
  }
  method <- match.arg(method, choices = c("cholesky", "svd", 
                                          "eigen", "RF", "circular.embedding"))
  if (messages.screen) {
    cat(messa$nst)
    cat(messa$nugget)
    cat(messa$cov.structures)
    if (method == "RF") 
      cat("grf: simulation using the function GaussRF from package RandomFields \n")
    else cat(paste("grf: decomposition algorithm used is: ", 
                   method, "\n"))
  }
  if (all(phi == 0)) 
    results$data <- matrix(rnorm((n * nsim), mean = 0, sd = sqrt(sill.total)), 
                           nrow = n, ncol = nsim)
  else {
    if (method == "RF") {                   #if I start dealing with random fields...
      RandomFields::RFoldstyle(old = TRUE)
      assign("setRF", geoR2RF(cov.model = cov.model, cov.pars = cov.pars, 
                              nugget = nugget, kappa = kappa, aniso.pars = aniso.pars), 
             pos = 1)
      if (!exists("xpts") || is.null(xpts)) {
        results$data <- RandomFields::GaussRF(x = results$coords[, 
                                                                 1], y = results$coords[, 2], model = get("setRF", 
                                                                                                          pos = 1), grid = FALSE, n = nsim)
      }
      else {
        results$data <- drop(matrix(RandomFields::GaussRF(x = xpts, 
                                                          y = ypts, model = get("setRF", pos = 1), grid = TRUE, 
                                                          n = nsim), ncol = nsim))
      }
    }
    
    else results$data <- drop(crossprod(varcov.spatial( 
      dists.lowertri=as.vector(NEdist), #BD This is where the Non-Euclidean distances are added in 
      #coords = results$coords, #BD no longer want this
      cov.model = cov.model, kappa = kappa, nugget = nugget, 
      cov.pars = cov.pars, only.decomposition = TRUE, func.inv = method)$sqrt.varcov, 
      matrix(rnorm((n * nsim)), nrow = n, ncol = nsim)))
  }
  if (length(mean) != 1 & length(mean) != dim(as.matrix(results$data))[1] & 
      length(mean) != length(results$data)) 
    stop("the mean must be a scalar or a vector of the same size as the data")
  results$data <- results$data + mean
  if (lambda != 1) {
    if (lambda != 0) 
      results$data <- (results$data * lambda + 1)^(1/lambda)
    else results$data <- exp(results$data)
    messa$transformation <- paste("grf: Data transformed (Box-Cox), for lambda =", 
                                  lambda)
    if (messages.screen) 
      cat(messa$transformation)
    cat("\n")
  }
  if (!RF && !is.null(aniso.pars)) {
    if (messages.screen) 
      cat("grf: back-transforming to the anisotropic space \n")
    results$coords <- coords.aniso(coords = results$coords, 
                                   aniso.pars = aniso.pars, reverse = TRUE)
  }
  else {
    aniso.pars <- "no anisotropy parameters provided/used"
  }
  if (messages.screen) 
    cat(paste("grf: End of simulation procedure. Number of realizations:", 
              nsim, "\n"))
  results[c("cov.model", "nugget", "cov.pars", "kappa", "lambda", 
            "aniso.pars", "method", ".Random.seed", "messages", "call")] <- list(cov.model, 
                                                                                 nugget, cov.pars, kappa = kappa, lambda = lambda, aniso.pars, 
                                                                                 method, rseed, messa, call.fc)
  attr(results, "borders") <- call.fc$borders
  attr(results, "sp.dim") <- ifelse(sim1d, "1d", "2d")
  oldClass(results) <- c("grf", "geodata", "variomodel")
  return(results)
}


testF <- function(coV){
  
  Krige <- t(mcmapply(krigeNE, cv=1:117, geodata=WQtrain.geo, locations=WQpredloc,
                      MoreArgs = list(covariance=as.matrix(coV), 
                                      krige=KC,output=OC), tolerance=1e-30, 
                      mc.cores=1))[,1:2]
  
  dev <- SE <- E <- rep(NA, 117)
  
  E <- unlist(Krige[,1]) - unlist(lapply(WQpred, function(x) x$data))  
  SE <- sqrt(unlist(Krige[,2]))
  
  dev <- unlist(lapply(WQpred, function(x) x$data))  - 
    mean(unlist(lapply(WQpred, function(x) x$data)))
  
  
  #for log-transformed
  
  #E <- expm1(unlist(Krige[,1])) - expm1(unlist(lapply(WQpred, function(x) x$data)))  
  #SE <- expm1(sqrt(unlist(Krige[,2])))
  
  #dev <- expm1(unlist(lapply(WQpred, function(x) x$data)))  - 
  #        mean(expm1(unlist(lapply(WQpred, function(x) x$data))))
  
  
  #expm1
  
  ME <- mean(E)
  MSE <- mean((E/SE)[is.finite(E/SE)])
  MAE <- mean(abs(E))
  RMSE <- sqrt(mean(E^2))
  RMSEstd <- sqrt(mean(((E/SE)[is.finite(E/SE)]^2)))
  Uncert <- mean(SE[is.finite(SE)])
  R2 <- 1- ( sum(E^2) / sum(dev^2) )
  NegVar <- sum(as.numeric(Krige[,2]<0))
  
  CI95 <- mapply( function(z){
    ifelse(WQpred[[z]]$data >= unlist(Krige[z,1]) - 1.96*SE[z] &
             WQpred[[1]]$data <= unlist(Krige[z,1]) + 
             1.96*SE[z], 1, 0)
  }, z=1:117)
  
  OUT <- list(Krige,E,SE, dev, ME, MSE, MAE,  RMSE,RMSEstd,Uncert,R2, NegVar, CI95)
  names(OUT) <- c("Krige","E","SE", "dev", "ME", "MSE", "MAE", "RMSE","RMSEstd",
                  "Uncert", "R2", "NegVar", "CI95")
  
  return(OUT)
  
}
