# Modified Ramsay's horseshoe
rm(list = ls())

# dependencies
library(tidyverse)
theme_set(theme_bw())
library(mgcv)
library(fields)
library(rgeos)
library(sf)
library(boraGP)
library(INLA)

# path 
path <- "~/boraGP-sss/"

# plot extensions
extension <- c(".pdf", ".png")

# horseshoe polygon
fsb <- fs.boundary()
p <- Polygon(cbind(fsb$x, fsb$y))
p <- Polygons(list(p), ID = "none")
poly <- SpatialPolygons(list(p))

# complement polygon
bb <- Polygon(cbind(c(-0.9, -0.9, 3.4, 3.4, -0.9), 
                    c(0.9, -0.9, -0.9, 0.9, 0.9)), hole=FALSE)
bbox <- SpatialPolygons(list(Polygons(list(bb), ID = "none")))
rest <- gDifference(bbox, poly)

rest_sf <- st_as_sf(rest)
hs_sf <- st_as_sf(poly)

#------#
# data #
#------#
xm <- 125 # 250
yn <- 50 # 100
x <- seq(-1, 4, length = xm)
y <- seq(-1, 1, length = yn)
xx <- rep(x, yn)
yy <- rep(y, rep(xm, yn))
z <- mgcv::fs.test(xx, yy)
zz <- matrix(z, xm, yn)

# is.na(w) = T outside the horseshoe
data <- data.frame(easting = xx, northing = yy, w = as.numeric(zz)) %>% 
  na.omit()

# create y 
tau_sq <- 0.1^2
set.seed(123)
data_hs <- data %>%
  mutate(y = w + sqrt(tau_sq)*rnorm(nrow(data)))

# remove 4 points that intersect with barriers 
coords_hs_sf <- st_as_sf(data_hs %>% select(easting, northing), 
                         coords = c("easting", "northing"))
dlt_idx <- unlist(st_intersects(rest_sf, coords_hs_sf$geometry))
data_hs <- data_hs[-dlt_idx,]
n <- nrow(data_hs)
rm(coords_hs_sf, dlt_idx)

# priors
nu <- 1
maxdist <- max(dist(as.matrix(data_hs[,c("easting", "northing")])))
d.low <- 0.25*maxdist
d.high <- 0.75*maxdist
dist_mat <- matrix(c(0, d.low, d.high, d.low, 0, 0, d.high, 0, 0), 
                   nrow = 3, ncol = 3, byrow = T)
phi.cand <- seq(1, 5, by=0.01)
cov.low = cov.high <- rep(0, length(phi.cand))
for (i in 1:length(phi.cand)) {
  cov <- Cov_matern(dist = dist_mat, sigmasq = 1, phi = phi.cand[i], nu = nu)
  cov.low[i] <- cov[1,2]
  cov.high[i] <- cov[1,3]
}
phi.high <- phi.cand[which.min(abs(cov.low - 0.05))]
phi.low <- phi.cand[which.min(abs(cov.high - 0.05))]
phi <- mean(c(phi.high, phi.low))

sigma.sq <- round(var(data_hs$y))
tau.sq <- 0.01
starting <- list("phi" = phi, "sigma.sq" = sigma.sq, 
                 "tau.sq" = tau.sq, "nu" = nu)
tuning <- list("phi" = 0.5, "sigma.sq" = 1, "tau.sq" = 0.005, "nu" = 0)         # nu is fixed 
priors <- list("phi.Unif" = c(phi.low, phi.high),                               # approx. 3/(0.55*maxdist), 3/(0.19*maxdist)
               "sigma.sq.IG" = c(2, sigma.sq), "tau.sq.IG" = c(2, tau.sq), 
               "nu.Unif" = c(nu-0.5,nu+0.5))

n.samples <- 10000
burn <- 5000

# simulation variables 
mlist <- c(10, 15, 20)
ntrlist <- c(300, 600, 1000)
rplc <- 30
set.seed(120)
seedsave <- sample.int(1e3, size = rplc*length(mlist)*length(ntrlist))
seedmat <- data.frame(M = rep(mlist, each = rplc*length(ntrlist)), 
                      N_tr = rep(rep(ntrlist, each = rplc), 
                                 times = length(mlist)), 
                      seed = seedsave)

# #------#
# # NNGP #
# #------#
# # save results
# predres = paramres = timeres <- list()
# 
# for (i in 1:length(ntrlist)) {
#   n_tr <- ntrlist[i]
#   n_tt <- n - n_tr
# 
#   for (j in 1:length(mlist)) {
#     m <- mlist[j]
#     listname <- paste0("n=",n_tr, ", m=", m)
# 
#     # save results
#     timeres_tmp <- matrix(0, nrow = rplc, ncol = 3)
#     colnames(timeres_tmp) <- c("user", "system", "elapsed")
#     predres_tmp = paramres_tmp <- matrix(0, nrow = rplc, ncol = 4)
#     colnames(predres_tmp) <- c("rmspe", "mape", "coverage", "meanwidth")
#     colnames(paramres_tmp) <- c("beta", "sigsq", "tausq", "phi")
#     seedlist <- seedmat %>% filter(M == m, N_tr == n_tr) %>% 
#       select(seed) %>% unlist() %>% as.numeric()
# 
#     cat(listname,'\n')
#     pb = txtProgressBar(style=3,width=50)
#     for (s in 1:rplc) {
#       seed <- seedlist[s]
#       set.seed(seed)
#       idx_tr <- sample.int(n, n_tr)
# 
#       # training data
#       coords_tr <- data_hs[idx_tr, c("easting", "northing")]
#       rownames(coords_tr) <- NULL
#       y_tr <- data_hs[idx_tr,"y"]
# 
#       # ordering (default)
#       ord <- order(coords_tr[,1])
# 
#       # test data
#       coords_tt <- data_hs[-idx_tr, c("easting", "northing")]
#       rownames(coords_tt) <- NULL
#       y_tt <- data_hs[-idx_tr,"y"]
# 
#       m.s <- spNNGP(y_tr ~ 1, coords = coords_tr, starting = starting,
#                     method = "response", n.neighbors = m,
#                     tuning = tuning, priors = priors, cov.model = "matern",
#                     n.samples = n.samples, n.omp.threads = 10, 
#                     ord = ord, verbose = F)
#       p.s <- predict(m.s, X.0 = matrix(1, nrow = n_tt, ncol = 1),
#                      sub.sample = 
#                        list(start = burn+1, end = n.samples, thin = 1),
#                      coords.0 = as.matrix(coords_tt), 
#                      n.omp.threads = 10, verbose = F)
# 
#       # time
#       timeres_tmp[s,] <- as.numeric(m.s$run.time)[1:3] + 
#         as.numeric(p.s$run.time)[1:3]
# 
#       # summary results
#       ystarNNGP <- rowMeans(p.s$p.y.0)
#       yquantNNGP <- apply(p.s$p.y.0, 1, 
#                           function(x) quantile(x, probs = c(0.025, 0.975)))
#       predres_tmp[s,] <- data.frame(y_tt = y_tt, yhat = ystarNNGP,
#                                     lower = yquantNNGP[1,], 
#                                     upper = yquantNNGP[2,]) %>%
#         mutate(error = y_tt - yhat,
#                width = upper-lower,
#                cover = ifelse((y_tt > lower & y_tt < upper), 1, 0)) %>%
#         summarise(rmspe = sqrt(mean(error^2)),
#                   mape = mean(abs(error)),
#                   coverage = mean(cover),
#                   meanwidth = mean(width)) %>%
#         as.numeric()
# 
#       ## beta
#       paramres_tmp[s,1] <- mean(m.s$p.beta.samples)
# 
#       ## theta
#       paramres_tmp[s,2:4] <- colMeans(m.s$p.theta.samples[,1:3])
# 
#       setTxtProgressBar(pb, s/rplc)
#     }
#     close(pb)
# 
#     timeres[[listname]] <- timeres_tmp
#     predres[[listname]] <- predres_tmp
#     paramres[[listname]] <- paramres_tmp
# 
#   }
# }
# 
# saveRDS(list(timeres = timeres, predres = predres,
#              paramres = paramres, seedmat = seedmat),
#         paste0(path, "sim/horseshoe_NNGP.RDS"))

# #-------------#
# # Barrier SGF #
# #-------------#
# # save results 
# predres = paramres = timeres <- list()
# 
# # for barrierSGF, varying m is meaningless 
# # use seed for each n at m = 10
# m <- 10
# 
# for (i in 1:length(ntrlist)) {
#   n_tr <- ntrlist[i]
#   n_tt <- n - n_tr
#   listname <- paste0("n=",n_tr)
#   
#   # save results
#   timeres_tmp <- matrix(0, nrow = rplc, ncol = 3) 
#   colnames(timeres_tmp) <- c("user", "system", "elapsed")
#   predres_tmp = paramres_tmp <- matrix(0, nrow = rplc, ncol = 4) 
#   colnames(predres_tmp) <- c("rmspe", "mape", "coverage", "meanwidth")
#   colnames(paramres_tmp) <- c("beta", "sigsq", "tausq", "phi")
#   seedlist <- seedmat %>% filter(M == m, N_tr == n_tr) %>% 
#     select(seed) %>% unlist() %>% as.numeric()
#   
#   cat(listname,'\n')
#   pb = txtProgressBar(style=3,width=50)
#   for (s in 1:rplc) {
#     seed <- seedlist[s]
#     set.seed(seed)
#     idx_tr <- sample.int(n, n_tr)
#     
#     # training data
#     coords_tr <- data_hs[idx_tr, c("easting", "northing")]
#     rownames(coords_tr) <- NULL
#     y_tr <- data_hs[idx_tr,"y"]
#     
#     # test data
#     coords_tt <- data_hs[-idx_tr, c("easting", "northing")]
#     rownames(coords_tt) <- NULL
#     y_tt <- data_hs[-idx_tr,"y"]
#     
#     # mesh 
#     max.edge <- 0.2
#     mesh <- inla.mesh.2d(boundary = poly,
#                          loc = coords_tr,
#                          max.edge = c(1,5)*max.edge,
#                          cutoff = 0.04,
#                          offset = c(max.edge, 1.5))
#     
#     # barrier model
#     tl <- length(mesh$graph$tv[,1])
#     posTri <- matrix(0, tl, 2)
#     for (t in 1:tl){
#       temp <- mesh$loc[mesh$graph$tv[t, ], ]
#       posTri[t,] <- colMeans(temp)[c(1,2)]
#     }
#     posTri <- SpatialPoints(posTri)
#     
#     normal <- over(poly, SpatialPoints(posTri), returnList = T)
#     barrier.triangles <- setdiff(1:tl, unlist(normal))
#     poly.barrier <- inla.barrier.polygon(mesh, barrier.triangles)
#     
#     # connect observations to mesh nodes
#     A.obs <- inla.spde.make.A(mesh, loc = as.matrix(coords_tr))
#     stk.obs <- inla.stack(data = list(y = y_tr),
#                           effects = list(s = 1:mesh$n,                          # spatial random effects
#                                          data.frame(int = rep(1,n_tr))),
#                           A = list(A.obs, 1),
#                           remove.unused = FALSE, tag = "obs")
#     
#     # same for prediction
#     proj.pred <- inla.mesh.projector(mesh, loc = as.matrix(coords_tt))
#     A.pred <- inla.spde.make.A(mesh, loc = proj.pred$loc)
#     stk.pred <- inla.stack(data = list(y = NA),
#                            A = list(A.pred, 1),
#                            effects = list(s = 1:mesh$n,
#                                           data.frame(int = rep(1,n_tt))),
#                            tag = "pred")
#     stk <- inla.stack(stk.obs, stk.pred)
#     
#     barrier.model <- inla.barrier.pcmatern(mesh, 
#                                            barrier.triangles = barrier.triangles,
#                                            prior.range = c(1, .5),              # P(range < 1) = 0.5
#                                            prior.sigma = c(3, 0.01))            # P(sigma > 3) = 0.01
#     
#     formula <- y~ -1 + int + f(s, model=barrier.model) 
#     time_inla <- system.time({
#       res <- inla(formula,
#                   data=inla.stack.data(stk),
#                   control.predictor=list(A=inla.stack.A(stk), compute = TRUE),
#                   control.compute=list(return.marginals.predictor=TRUE),
#                   family="gaussian", 
#                   control.inla= list(int.strategy = "eb"), num.threads = 10)
#     })
#     
#     # time 
#     timeres_tmp[s,] <- as.numeric(time_inla)[1:3]
#     
#     # summary results
#     index.pred <- c(inla.stack.index(stk, "pred")$data)
#     y_inla <- res$summary.fitted.values[index.pred, "mean"]
#     inla_sd <- sqrt(res$summary.fitted.values[index.pred,"sd"]^2                # variability from the coefficient
#                     + 1/res$summary.hyperpar[1,"mean"])                         # tausq
#     y_inla_low <- y_inla + qnorm(.025)*inla_sd
#     y_inla_high <- y_inla + qnorm(.975)*inla_sd
#     predres_tmp[s,] <- data.frame(y_tt = y_tt, yhatI = y_inla, 
#                                   lowerI = y_inla_low, upperI = y_inla_high) %>% 
#       mutate(errorI = y_tt - yhatI,
#              widthI = upperI-lowerI,
#              coverI = ifelse((y_tt > lowerI & y_tt < upperI), 1, 0)) %>% 
#       summarise(rmspe = sqrt(mean(errorI^2)),
#                 mape = mean(abs(errorI)),
#                 coverage = mean(coverI),
#                 meanwidth = mean(widthI)) %>% 
#       as.numeric()
#     
#     ## beta 
#     paramres_tmp[s,1] <- as.numeric(res$summary.fixed[1])
#     
#     ## theta
#     # sigma = exp(theta1)
#     paramres_tmp[s,2] <- res$internal.marginals.hyperpar[2] %>%
#       lapply(function(m) {
#         inla.tmarginal(function(x) exp(x)^2, m)
#       }) %>%
#       sapply(function(m)
#         unlist(inla.zmarginal(m, silent = TRUE))[1]) %>% as.numeric()
#     
#     # tausq
#     paramres_tmp[s,3] <- res$internal.marginals.hyperpar[1] %>%
#       lapply(function(m) {
#         inla.tmarginal(function(x) 1/exp(x), m)
#       }) %>%
#       sapply(function(m)
#         unlist(inla.zmarginal(m, silent = TRUE))[1]) %>% as.numeric()
#     
#     # phi = sqrt(8)/r where spatial range r = exp(theta2)
#     paramres_tmp[s,4] <- res$internal.marginals.hyperpar[3] %>%
#       lapply(function(m) {
#         inla.tmarginal(function(x) sqrt(8)/exp(x), m) 
#       }) %>%
#       sapply(function(m)
#         unlist(inla.zmarginal(m, silent = TRUE))[1]) %>% as.numeric()
#     
#     setTxtProgressBar(pb, s/rplc) 
#   }
#   close(pb)
#   
#   timeres[[listname]] <- timeres_tmp
#   predres[[listname]] <- predres_tmp
#   paramres[[listname]] <- paramres_tmp
# }
# 
# saveRDS(list(timeres = timeres, predres = predres, 
#              paramres = paramres, seedmat = seedmat),
#         paste0(path, "sim/horseshoe_barrierSGF.RDS"))

# #---------#
# # BORA-GP #
# #---------#
# # save results 
# timenb = predres = paramres = timeres <- list()
# 
# for (i in 1:length(ntrlist)) {
#   n_tr <- ntrlist[i]
#   n_tt <- n - n_tr
#   
#   for (j in 1:length(mlist)) {
#     m <- mlist[j]
#     listname <- paste0("n=",n_tr, ", m=", m)
#     
#     # save results
#     timenb_tmp = timeres_tmp <- matrix(0, nrow = rplc, ncol = 3) 
#     colnames(timenb_tmp) = colnames(timeres_tmp) <- 
#       c("user", "system", "elapsed")
#     predres_tmp = paramres_tmp <- matrix(0, nrow = rplc, ncol = 4) 
#     colnames(predres_tmp) <- c("rmspe", "mape", "coverage", "meanwidth")
#     colnames(paramres_tmp) <- c("beta", "sigsq", "tausq", "phi")
#     seedlist <- seedmat %>% filter(M == m, N_tr == n_tr) %>% 
#       select(seed) %>% unlist() %>% as.numeric()
#     
#     for (s in 1:rplc) {
#       seed <- seedlist[s]
#       set.seed(seed)
#       idx_tr <- sample.int(n, n_tr)
#       
#       # training data
#       coords_tr <- data_hs[idx_tr, c("easting", "northing")]
#       rownames(coords_tr) <- NULL
#       y_tr <- data_hs[idx_tr,"y"]
#       
#       # ordering (default)
#       ord <- order(coords_tr[,1])
#       
#       # test data
#       coords_tt <- data_hs[-idx_tr, c("easting", "northing")]
#       rownames(coords_tt) <- NULL
#       y_tt <- data_hs[-idx_tr,"y"]
#       
#       coords_all <- rbind(coords_tr, coords_tt)
#       coords_sf <- st_as_sf(coords_all, coords = c("easting", "northing"))
#       
#       time_nb <- system.time(
#         barrier_nninfo_all <- barrier_neighbor(coords = coords_tr,
#                                                coords.0 = coords_tt,
#                                                ord = ord,
#                                                n.neighbors = m,
#                                                barrier = rest_sf,
#                                                cores = 20,
#                                                verbose = F,
#                                                debug = list(barrier_n.indx = NULL,
#                                                             barrier_dist = NULL,
#                                                             barrier_nn.indx.0_list = NULL,
#                                                             barrier_dist0 = NULL,
#                                                             ref_window = NULL,
#                                                             nonref_window = NULL,
#                                                             ref_fill = TRUE,
#                                                             nonref_fill = TRUE))
#       )
#       barrier_nn.indx.0_list <- barrier_nninfo_all$barrier_nn.indx.0_list
#       barrier_nn.indx.0 <- do.call("rbind", barrier_nn.indx.0_list)
#       barrier_n.indx <- barrier_nninfo_all$barrier_n.indx
#       barrier_nn.indx <- as.integer(unlist(barrier_n.indx)[-1]-1)
#       barrier_nn.indx.lu <- c(cumsum(c(0,0,sapply(barrier_n.indx, length)[-1]))[1:n_tr],
#                               c(0,sapply(barrier_n.indx, length)[-1])) %>% as.integer()
#       barrier_nninfo <- list(type = "barrier",
#                              n.indx = barrier_n.indx,
#                              n.neighbors = m, nn.indx = barrier_nn.indx,
#                              nn.indx.lu = barrier_nn.indx.lu, ord = ord)
#       
#       barrier_m.s <- spNNGP(y_tr ~ 1, coords = coords_tr, starting = starting,
#                             method = "response", n.neighbors = m,
#                             tuning = tuning, priors = priors, 
#                             cov.model = "matern",
#                             n.samples = n.samples, n.omp.threads = 10,
#                             neighbor.info = barrier_nninfo, verbose = F)
#       barrier_p.s <- predict(barrier_m.s, X.0 = matrix(1, nrow = n_tt, ncol = 1),
#                              coords.0 = as.matrix(coords_tt),
#                              sub.sample = 
#                                list(start = burn+1, end = n.samples, thin = 1),
#                              nn.indx.0 = barrier_nn.indx.0, 
#                              n.omp.threads = 10, verbose = F)
#       
#       # time 
#       timenb_tmp[s,] <- as.numeric(time_nb)[1:3]
#       timeres_tmp[s,] <- as.numeric(barrier_m.s$run.time)[1:3] + 
#         as.numeric(barrier_p.s$run.time)[1:3]
#       
#       # summary results
#       ystarBRGP <- rowMeans(barrier_p.s$p.y.0)
#       yquantBRGP <- apply(barrier_p.s$p.y.0, 1, 
#                           function(x) quantile(x, probs = c(0.025, 0.975)))
#       predres_tmp[s,] <- data.frame(y_tt = y_tt, yhat = ystarBRGP, 
#                                     lower = yquantBRGP[1,], 
#                                     upper = yquantBRGP[2,]) %>% 
#         mutate(error = y_tt - yhat,
#                width = upper-lower,
#                cover = ifelse((y_tt > lower & y_tt < upper), 1, 0)) %>% 
#         summarise(rmspe = sqrt(mean(error^2)),
#                   mape = mean(abs(error)),
#                   coverage = mean(cover),
#                   meanwidth = mean(width)) %>% 
#         as.numeric()
#       
#       ## beta
#       paramres_tmp[s,1] <- mean(barrier_m.s$p.beta.samples)
#       
#       ## theta
#       paramres_tmp[s,2:4] <- colMeans(barrier_m.s$p.theta.samples[,1:3])
#       
#       cat(paste0(listname, ": ", s, "th simulation completed.\n"))
#     }
#     
#     timenb[[listname]] <- timenb_tmp
#     timeres[[listname]] <- timeres_tmp
#     predres[[listname]] <- predres_tmp
#     paramres[[listname]] <- paramres_tmp
#     
#     saveRDS(list(timenb = timenb_tmp, 
#                  timeres = timeres_tmp, 
#                  predres = predres_tmp, 
#                  paramres = paramres_tmp),
#             paste0(path, "sim/horseshoe/n", n_tr,"m",m,".RDS"))
#   }
# }
# 
# saveRDS(list(timenb = timenb, timeres = timeres, predres = predres, 
#              paramres = paramres, seedmat = seedmat),
#         paste0(path, "sim/horseshoe_BORAGP.RDS"))

# call results
resBRGP <- readRDS(paste0(path, "sim/horseshoe_BORAGP.RDS"))
resNNGP <- readRDS(paste0(path, "sim/horseshoe_NNGP.RDS"))
resinla <- readRDS(paste0(path, "sim/horseshoe_barrierSGF.RDS"))

# results
resnames <- c("RMSPE", "MAPE", "95% CI coverage", "Mean 95% CI width")
parnames <- c("beta", "sigmasq", "tausq", "phi")
dfBRGPL = dfNNGPL = dfinlaL <- list()
for (i in 1:4) {
  dfBRGPL[[i]] <- data.frame(n = rep(ntrlist, each = rplc*length(mlist)),
                             m = rep(rep(mlist, each = rplc), 
                                     times = length(ntrlist)),
                             model = "BORA-GP",
                             vals = c(sapply(resBRGP$predres, 
                                             function(x) x[,i])),
                             what = resnames[i],
                             vals2 = c(sapply(resBRGP$paramres, 
                                              function(x) x[,i])),
                             what2 = parnames[i])
  dfNNGPL[[i]] <- data.frame(n = rep(ntrlist, each = rplc*length(mlist)),
                             m = rep(rep(mlist, each = rplc), 
                                     times = length(ntrlist)),
                             model = "NNGP",
                             vals = c(sapply(resNNGP$predres, 
                                             function(x) x[,i])),
                             what = resnames[i],
                             vals2 = c(sapply(resNNGP$paramres, 
                                              function(x) x[,i])),
                             what2 = parnames[i])
  dfinlaL[[i]] <- data.frame(n = rep(rep(ntrlist, each = rplc), 
                                     times = length(mlist)),
                             m = rep(mlist, each = rplc*length(ntrlist)),
                             model = "Barrier SGF",
                             vals = rep(c(sapply(resinla$predres, 
                                                 function(x) x[,i])), times = length(mlist)),
                             what = resnames[i],
                             vals2 = rep(c(sapply(resinla$paramres, 
                                                  function(x) x[,i])), times = length(mlist)),
                             what2 = parnames[i])
}

dfBRGP <- do.call(rbind, dfBRGPL)
dfNNGP <- do.call(rbind, dfNNGPL)
dfinla <- do.call(rbind, dfinlaL)
dfres <- rbind(dfBRGP, dfNNGP, dfinla)

## prediction performance
sapply(resBRGP$predres, colMeans) %>% round(3)
sapply(resBRGP$predres, function(x) apply(x, 2, sd)) %>% round(3)
sapply(resNNGP$predres, colMeans) %>% round(3)
sapply(resNNGP$predres, function(x) apply(x, 2, sd)) %>% round(3)
sapply(resinla$predres, colMeans) %>% round(3)
sapply(resinla$predres, function(x) apply(x, 2, sd)) %>% round(3)

## prediction performance in plots
dfres_sum <- dfres %>%
  group_by(n,m,model,what) %>%
  summarise(vals_mean = mean(vals),
            vals_low = quantile(vals, prob = 0.025),
            vals_hi = quantile(vals, prob = 0.975))
dfres_sum$m[dfres_sum$model == "NNGP"] <- 
  dfres_sum$m[dfres_sum$model == "NNGP"] + 0.5
dfres_sum$m[dfres_sum$model == "BORA-GP"] <- 
  dfres_sum$m[dfres_sum$model == "BORA-GP"] - 0.5

sumplot1 <- dfres_sum %>%
  mutate(model = factor(model, levels = c("BORA-GP", "NNGP", "Barrier SGF"))) %>%
  filter(what == "RMSPE") %>%
  ggplot() +
  geom_line(aes(m, vals_mean, col = model, linetype = model)) +
  geom_errorbar(aes(m, ymin = vals_low, ymax = vals_hi, col = model, linetype = model), width = 0.5) +
  geom_point(aes(m, vals_mean, col = model), size = 2) +
  facet_grid(what ~
               factor(n, labels = c("n = 300", "n = 600", "n = 1000")), scales = "free") +
  labs(x="", y="", color="", linetype="") +
  scale_color_scico_d(palette = "batlow", begin = 0.2, end = 0.9, direction = -1) +
  theme(legend.position = "right",
        plot.margin=margin(t=0,l=-10, r=0, b=-10),
        legend.margin=margin(b=0,r=0,t=0,l=-3)) +
  scale_x_continuous(breaks = c(10,15,20), labels = c("m = 10", "m = 15", "m = 20"))

sumplot2 <- dfres %>%
  mutate(model = factor(model, levels = c("BORA-GP", "NNGP", "Barrier SGF")),
         m = factor(m, levels = c("10", "15", "20"), labels = c("m = 10", "m = 15", "m = 20"))) %>%
  filter(what == "95% CI coverage") %>%
  ggplot() +
  geom_hline(data = data.frame(y=0.95, what = "95% CI coverage"), aes(yintercept = y),
             col = "#601200", linetype = "dashed") +
  geom_boxplot(aes(m, vals, fill = model), outlier.shape = NA) +
  facet_grid(what ~
               factor(n, labels = c("n = 300", "n = 600", "n = 1000")), scales = "free") +
  labs(x="", y="", fill="") +
  scale_fill_scico_d(palette = "batlow", begin = 0.2, end = 0.9, direction = -1) +
  theme(legend.position = "right",
        plot.margin=margin(t=0,l=-10, r=0, b=-10),
        legend.margin=margin(b=0,r=0,t=0,l=-3))

gg <- ggpubr::ggarrange(sumplot1, sumplot2, nrow = 2)
# for (ext in extension) {
#   ggsave(plot = gg,
#          paste0(path, "plots/horseshoe_pred", ext),
#          width = 9, height = 4)
# }

sumplot3 <- dfres_sum %>%
  mutate(model = factor(model, levels = c("BORA-GP", "NNGP", "Barrier SGF"))) %>%
  filter(what %in% c("MAPE", "Mean 95% CI width")) %>%
  ggplot() +
  geom_line(aes(m, vals_mean, col = model, linetype = model)) +
  geom_errorbar(aes(m, ymin = vals_low, ymax = vals_hi, col = model, linetype = model), width = 0.5) +
  geom_point(aes(m, vals_mean, col = model), size = 2) +
  facet_grid(factor(what, levels = c("MAPE", "Mean 95% CI width"),
                    labels = c("MAPE", "95% CI width")) ~
               factor(n, labels = c("n = 300", "n = 600", "n = 1000")), scales = "free") +
  labs(x="", y="", color="", linetype="") +
  scale_color_scico_d(palette = "batlow", begin = 0.2, end = 0.9, direction = -1) +
  theme(legend.position = "right",
        plot.margin=margin(t=0,l=-10, r=0, b=-10),
        legend.margin=margin(b=0,r=0,t=0,l=-3)) +
  scale_x_continuous(breaks = c(10,15,20), labels = c("m = 10", "m = 15", "m = 20"))
# for (ext in extension) {
#   ggsave(plot = sumplot3,
#          paste0(path, "plots/horseshoe_pred_supp", ext),
#          width = 9, height = 3.5)
# }

## parameter estimation
sapply(resBRGP$paramres, colMeans) %>% round(3)
sapply(resBRGP$paramres, function(x) apply(x, 2, sd)) %>% round(3)
sapply(resNNGP$paramres, colMeans) %>% round(3)
sapply(resNNGP$paramres, function(x) apply(x, 2, sd)) %>% round(3)
sapply(resinla$paramres, colMeans) %>% round(3)
sapply(resinla$paramres, function(x) apply(x, 2, sd)) %>% round(3)

## time
((sapply(resBRGP$timenb, colMeans) + sapply(resBRGP$timeres, colMeans))/60) %>% 
  round(3)
(sapply(resNNGP$timeres, colMeans)/60) %>% round(3)
(sapply(resinla$timeres, colMeans)/60) %>% round(3)

#--------------#
# latent model #
#--------------#
# simulate data
n_tr <- ntrlist[1]
n_tt <- n - n_tr
m <- mlist[3]
seedlist <- seedmat %>% filter(M == m, N_tr == n_tr) %>% 
  select(seed) %>% unlist() %>% as.numeric()

set.seed(seedlist[27]) 
idx_tr <- sample.int(n, n_tr)

# training data
coords_tr <- data_hs[idx_tr, c("easting", "northing")]
rownames(coords_tr) <- NULL
w_tr <- data_hs[idx_tr,"w"]
y_tr <- data_hs[idx_tr,"y"]

# ordering (default)
ord <- order(coords_tr[,1])

# test data
coords_tt <- data_hs[-idx_tr, c("easting", "northing")]
rownames(coords_tt) <- NULL
w_tt <- data_hs[-idx_tr,"w"]
y_tt <- data_hs[-idx_tr,"y"]

coords_all <- rbind(coords_tr, coords_tt)

# # run models
# set.seed(seedlist[1])
# time_nb <- system.time(
#   barrier_nninfo_all <- barrier_neighbor(coords = coords_tr,
#                                          coords.0 = coords_tt,
#                                          ord = ord,
#                                          n.neighbors = m,
#                                          barrier = rest_sf,
#                                          cores = 20,
#                                          verbose = T,
#                                          debug = list(barrier_n.indx = NULL,
#                                                       barrier_dist = NULL,
#                                                       barrier_nn.indx.0_list = NULL,
#                                                       barrier_dist0 = NULL,
#                                                       ref_window = NULL,
#                                                       nonref_window = NULL,
#                                                       ref_fill = TRUE,
#                                                       nonref_fill = TRUE))
# )
# barrier_nn.indx.0_list <- barrier_nninfo_all$barrier_nn.indx.0_list
# barrier_nn.indx.0 <- do.call("rbind", barrier_nn.indx.0_list)
# barrier_n.indx <- barrier_nninfo_all$barrier_n.indx
# barrier_nn.indx <- as.integer(unlist(barrier_n.indx)[-1]-1)
# barrier_nn.indx.lu <- c(cumsum(c(0,0,sapply(barrier_n.indx, length)[-1]))[1:n_tr],
#                         c(0,sapply(barrier_n.indx, length)[-1])) %>% as.integer()
# barrier_nninfo <- list(type = "barrier",
#                        n.indx = barrier_n.indx,
#                        n.neighbors = m, nn.indx = barrier_nn.indx,
#                        nn.indx.lu = barrier_nn.indx.lu, ord = ord)
# 
# barrier_m.s <- spNNGP(y_tr ~ 1, coords = coords_tr, starting = starting,
#                       method = "latent", n.neighbors = m,
#                       tuning = tuning, priors = priors,
#                       cov.model = "matern",
#                       n.samples = n.samples, n.omp.threads = 10,
#                       neighbor.info = barrier_nninfo, verbose = T)
# barrier_p.s <- predict(barrier_m.s, X.0 = matrix(1, nrow = n_tt, ncol = 1),
#                        coords.0 = as.matrix(coords_tt),
#                        sub.sample =
#                          list(start = burn+1, end = n.samples, thin = 1),
#                        nn.indx.0 = barrier_nn.indx.0,
#                        n.omp.threads = 10, verbose = T)
# 
# m.s <- spNNGP(y_tr ~ 1, coords = coords_tr, starting = starting,
#               method = "latent", n.neighbors = m,
#               tuning = tuning, priors = priors, cov.model = "matern",
#               n.samples = n.samples, n.omp.threads = 10, ord = ord, verbose = T)
# p.s <- predict(m.s, X.0 = matrix(1, nrow = n_tt, ncol = 1),
#                sub.sample = list(start = burn+1, end = n.samples, thin = 1),
#                coords.0 = as.matrix(coords_tt), n.omp.threads = 10, verbose = T)
# 
# mesh
max.edge <- 0.2
mesh <- inla.mesh.2d(boundary = poly,
                     loc = coords_tr,
                     max.edge = c(1,5)*max.edge,
                     cutoff = 0.04,
                     offset = c(max.edge, 1.5))
# 
# # barrier model
# tl <- length(mesh$graph$tv[,1])
# # - the number of triangles in the mesh != mesh$n
# posTri <- matrix(0, tl, 2)
# for (t in 1:tl){
#   temp <- mesh$loc[mesh$graph$tv[t, ], ]
#   posTri[t,] <- colMeans(temp)[c(1,2)]
# }
# posTri <- SpatialPoints(posTri)
# # - compute the triangle positions
# 
# normal <- over(poly, SpatialPoints(posTri), returnList = T)
# barrier.triangles <- setdiff(1:tl, unlist(normal))
# poly.barrier <- inla.barrier.polygon(mesh, barrier.triangles)
# 
# # connect observations to mesh nodes
# A.obs <- inla.spde.make.A(mesh, loc = as.matrix(coords_tr))
# stk.obs <- inla.stack(data = list(y = y_tr),
#                       effects = list(s = 1:mesh$n, # spatial random effects
#                                      data.frame(int = rep(1,n_tr))),
#                       A = list(A.obs, 1),
#                       remove.unused = FALSE, tag = "obs")
# 
# # same for prediction
# proj.pred <- inla.mesh.projector(mesh, loc = as.matrix(coords_tt))
# A.pred <- inla.spde.make.A(mesh, loc = proj.pred$loc)
# stk.pred <- inla.stack(data = list(y = NA),
#                        A = list(A.pred, 1),
#                        effects = list(s = 1:mesh$n,
#                                       data.frame(int = rep(1,n_tt))),
#                        tag = "pred")
# stk <- inla.stack(stk.obs, stk.pred)
# 
# barrier.model <- inla.barrier.pcmatern(mesh, barrier.triangles = barrier.triangles,
#                                        prior.range = c(1, .5), # P(range < 1) = 0.5
#                                        prior.sigma = c(3, 0.01)) # P(sigma > 3) = 0.01
# 
# formula <- y~ -1 + int + f(s, model=barrier.model)
# time_inla <- system.time({
#   res <- inla(formula,
#               data=inla.stack.data(stk),
#               control.predictor=list(A=inla.stack.A(stk), compute = TRUE),
#               control.compute=list(return.marginals.predictor=TRUE),
#               family="gaussian",
#               control.inla= list(int.strategy = "eb"), num.threads = 10)
# })
# 
# saveRDS(list(barrier_m.s = barrier_m.s, barrier_p.s = barrier_p.s,
#              m.s = m.s, p.s = p.s, res = res),
#         paste0(path, "sim/horseshoe_latent.RDS"))
# 
# # summary 
# whatBRGP <- rowMeans(barrier_m.s$p.w.samples[,-c(1:burn)]) +
#   mean(barrier_m.s$p.beta.samples[-c(1:burn)])
# wstarBRGP <- rowMeans(barrier_p.s$p.w.0) +
#   mean(barrier_m.s$p.beta.samples[-c(1:burn)])
# whatNNGP <- rowMeans(m.s$p.w.samples[,-c(1:burn)]) +
#   mean(m.s$p.beta.samples[-c(1:burn)])
# wstarNNGP <- rowMeans(p.s$p.w.0) +
#   mean(m.s$p.beta.samples[-c(1:burn)])
# wproj <- inla.mesh.projector(mesh, loc = as.matrix(coords_all))
# winla <- inla.mesh.project(wproj, res$summary.random$s)[,"mean"] +
#   res$summary.fixed$mean
# saveRDS(list(whatBRGP = whatBRGP, 
#              wstarBRGP = wstarBRGP, 
#              whatNNGP = whatNNGP, 
#              wstarNNGP = wstarNNGP, 
#              wproj = wproj, 
#              winla = winla), 
#         paste0(path, "sim/horseshoe_latent_res.RDS"))

# horseshoe plots
latent_res <- readRDS(paste0(path, "sim/horseshoe_latent_res.RDS"))
whatBRGP <- latent_res$whatBRGP
wstarBRGP <- latent_res$wstarBRGP
whatNNGP <- latent_res$whatNNGP
wstarNNGP <- latent_res$wstarNNGP
winla <- latent_res$winla

gg_w <- data.frame(easting = rep(coords_all$easting, times = 4),
           northing = rep(coords_all$northing, times = 4),
           what = c(w_tr, w_tt, whatNNGP, wstarNNGP, whatBRGP, wstarBRGP, winla),
           model = rep(c("Truth", "NNGP", "BORA-GP", "Barrier SGF"), each = n)) %>%
  ggplot() +
  facet_wrap(~factor(model, 
                     levels = c("Truth", "BORA-GP", "NNGP", "Barrier SGF")), 
             nrow = 2) +
  geom_raster(aes(easting, northing, fill = what)) +
  geom_sf(data = hs_sf, fill = NA) +
  geom_contour(aes(x = easting, y = northing, z = what), 
               color = "black", breaks = seq(-4, 4, by=.5)) +
  scale_fill_scico(palette = "roma", direction =-1) +
  labs(x="", y="", fill = "w") +
  theme(plot.margin = margin(t = -10, l = -10, r = 0, b = -20), aspect.ratio = 1/2,
        legend.margin = margin(b = 0, r = 0, t = 0, l = -5))
# for (ext in extension) {
#   ggsave(plot = gg_w,
#          paste0(path, "plots/horseshoe_w", ext),
#          width = 6, height = 3.6)
# }

