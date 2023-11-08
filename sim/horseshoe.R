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
library(scico)
library(fdaPDE)
library(gdistance)
library(raster)
library(geoR)
library(sfsmisc)

# path 
path <- "~/boraGP-sss/"

# call closePD 
source(paste0(path, "sim/ClosePD_042418.R"))

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
coords_hs_sf <- st_as_sf(data_hs %>% dplyr::select(easting, northing), 
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
# timemesh = predres = paramres = timeres <- list()
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
#   timemesh_tmp <- matrix(0, nrow = rplc, ncol = 3)
#   colnames(timemesh_tmp) <- c("user", "system", "elapsed")
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
#     time_mesh <- system.time({
#       max.edge <- 0.2
#       mesh <- inla.mesh.2d(boundary = poly,
#                            loc = coords_tr,
#                            max.edge = c(1,5)*max.edge,
#                            cutoff = 0.04,
#                            offset = c(max.edge, 1.5))
#       
#       # barrier model
#       tl <- length(mesh$graph$tv[,1])
#       posTri <- matrix(0, tl, 2)
#       for (t in 1:tl){
#         temp <- mesh$loc[mesh$graph$tv[t, ], ]
#         posTri[t,] <- colMeans(temp)[c(1,2)]
#       }
#       posTri <- SpatialPoints(posTri)
#       
#       normal <- over(poly, SpatialPoints(posTri), returnList = T)
#       barrier.triangles <- setdiff(1:tl, unlist(normal))
#       poly.barrier <- inla.barrier.polygon(mesh, barrier.triangles)
#       
#       # connect observations to mesh nodes
#       A.obs <- inla.spde.make.A(mesh, loc = as.matrix(coords_tr))
#       stk.obs <- inla.stack(data = list(y = y_tr),
#                             effects = list(s = 1:mesh$n,                          # spatial random effects
#                                            data.frame(int = rep(1,n_tr))),
#                             A = list(A.obs, 1),
#                             remove.unused = FALSE, tag = "obs")
#       
#       # same for prediction
#       proj.pred <- inla.mesh.projector(mesh, loc = as.matrix(coords_tt))
#       A.pred <- inla.spde.make.A(mesh, loc = proj.pred$loc)
#       stk.pred <- inla.stack(data = list(y = NA),
#                              A = list(A.pred, 1),
#                              effects = list(s = 1:mesh$n,
#                                             data.frame(int = rep(1,n_tt))),
#                              tag = "pred")
#       stk <- inla.stack(stk.obs, stk.pred)
#       
#       barrier.model <- inla.barrier.pcmatern(mesh,
#                                              barrier.triangles = barrier.triangles,
#                                              prior.range = c(1, .5),              # P(range < 1) = 0.5
#                                              prior.sigma = c(3, 0.01))            # P(sigma > 3) = 0.01
#       
#       formula <- y~ -1 + int + f(s, model=barrier.model)
#     })
#     
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
#     timemesh_tmp[s,] <- as.numeric(time_mesh)[1:3]
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
#   timemesh[[listname]] <- timemesh_tmp
#   timeres[[listname]] <- timeres_tmp
#   predres[[listname]] <- predres_tmp
#   paramres[[listname]] <- paramres_tmp
# }
# 
# saveRDS(list(timemesh = timemesh, timeres = timeres, predres = predres,
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

# #--------------------#
# # soap film smoother #
# #--------------------#
# # link: https://fromthebottomoftheheap.net/2016/03/27/soap-film-smoothers/
#
# boundary
# hs_soap <- list(list(easting = fsb$x, northing = fsb$y))                        # estimate boundary values
# 
# # save results
# predres = paramres = timeres <- list()
# 
# for (i in 1:length(ntrlist)) {
#   n_tr <- ntrlist[i]
#   n_tt <- n - n_tr
# 
#   # for soap film smoother, varying m is meaningless
#   # use seed for each n at m = 10
#   m <- 10
#   listname <- paste0("n=",n_tr)
#   
#   # save results
#   timeres_tmp <- matrix(0, nrow = rplc, ncol = 3)
#   colnames(timeres_tmp) <- c("user", "system", "elapsed")
#   predres_tmp = paramres_tmp <- matrix(NA, nrow = rplc, ncol = 4)
#   colnames(predres_tmp) <- c("rmspe", "mape", "coverage", "meanwidth")
#   colnames(paramres_tmp) <- c("beta", "sigsq", "tausq", "phi")
#   seedlist <- seedmat %>% filter(M == m, N_tr == n_tr) %>%
#     dplyr::select(seed) %>% unlist() %>% as.numeric()
#   
#   for (s in 1:rplc) {
#     seed <- seedlist[s]
#     set.seed(seed)
#     idx_tr <- sample.int(n, n_tr)
#     
#     # training data
#     coords_tr <- data_hs[idx_tr, c("easting", "northing")]
#     rownames(coords_tr) <- NULL
#     y_tr <- data_hs[idx_tr,"y"]
#     data_tr <- tibble(coords_tr, y = y_tr)
#     
#     # test data
#     coords_tt <- data_hs[-idx_tr, c("easting", "northing")]
#     rownames(coords_tt) <- NULL
#     y_tt <- data_hs[-idx_tr,"y"]
#     data_tt <- tibble(coords_tt, y = y_tt)
#     
#     # knots 
#     test <- mgcv::fs.test(x = coords_tr$easting, y = coords_tr$northing, 
#                           r0 = 0.15) 
#     knots <- coords_tr[-which(is.na(test)),]
# 
#     # model fitting
#     time_soap <- system.time({
#       soap <- gam(y ~ s(easting, northing, bs = "so", xt = list(bnd = hs_soap)),
#                   data = data_tr, method = "REML", knots = knots)
#       soap_pred <- predict(soap, newdata = data_tt, se.fit = TRUE, type = "response")
#     })
#     
#     # time
#     timeres_tmp[s,] <- as.numeric(time_soap)[1:3]
#     
#     # summary results
#     y_soap <- as.numeric(soap_pred$fit)
#     soap_sd <- as.numeric(soap_pred$se.fit)
#     y_soap_low <- y_soap + qnorm(.025)*soap_sd
#     y_soap_high <- y_soap + qnorm(.975)*soap_sd
#     predres_tmp[s,] <- data.frame(y_tt = y_tt, yhat = y_soap,
#                                   lower = y_soap_low, upper = y_soap_high) %>%
#       mutate(error = y_tt - yhat,
#              width = upper-lower,
#              cover = ifelse((y_tt > lower & y_tt < upper), 1, 0)) %>%
#       summarise(rmspe = sqrt(mean(error^2)),
#                 mape = mean(abs(error)),
#                 coverage = mean(cover),
#                 meanwidth = mean(width)) %>%
#       as.numeric()
#     
#     ## beta
#     paramres_tmp[s,1] <- soap$coefficients[1]
#     
#     ## tausq
#     paramres_tmp[s,3] <- soap$reml.scale
# 
#     cat(paste0(listname, ": ", s, "th simulation completed.\n"))
#     
#     saveRDS(list(timeres = timeres_tmp,
#                  predres = predres_tmp,
#                  paramres = paramres_tmp),
#             paste0(path, "sim/horseshoe_other_soap/n", n_tr, ".RDS"))
#   }
#   
#   timeres[[listname]] <- timeres_tmp
#   predres[[listname]] <- predres_tmp
#   paramres[[listname]] <- paramres_tmp
#   
# }
# 
# saveRDS(list(timeres = timeres, predres = predres,
#              paramres = paramres, seedmat = seedmat),
#         paste0(path, "sim/horseshoe_soap.RDS"))

# #--------#
# # SR-PDE #
# #--------#
# # link
# # https://onlinelibrary.wiley.com/doi/epdf/10.1111/insr.12444
# # https://academic.oup.com/jrsssb/article/75/4/681/7075920
# 
# # boundary nodes
# hs_nodes <- cbind(fsb$x, fsb$y)
# colnames(hs_nodes) <- c("easting", "northing")
# hs_segments <- cbind(1:length(fsb$x), c(2:length(fsb$y),1))
# 
# # smoothing parameters
# lambda <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)
# 
# # save results
# timepre = predres = paramres = timeres = lambda_pde <- list()
# 
# for (i in 1:length(ntrlist)) {
#   n_tr <- ntrlist[i]
#   n_tt <- n - n_tr
# 
#   # for SR-PDE, varying m is meaningless
#   # use seed for each n at m = 10
#   m <- 10
#   listname <- paste0("n=",n_tr)
# 
#   # save results
#   lambda_tmp <- rep(0, rplc)
#   timeres_tmp = timepre_tmp <- matrix(0, nrow = rplc, ncol = 3)
#   colnames(timeres_tmp) = colnames(timepre_tmp) <-
#     c("user", "system", "elapsed")
#   predres_tmp = paramres_tmp <- matrix(NA, nrow = rplc, ncol = 4)
#   colnames(predres_tmp) <- c("rmspe", "mape", "coverage", "meanwidth")
#   colnames(paramres_tmp) <- c("beta", "sigsq", "tausq", "phi")
#   seedlist <- seedmat %>% filter(M == m, N_tr == n_tr) %>%
#     dplyr::select(seed) %>% unlist() %>% as.numeric()
# 
#   for (s in 1:rplc) {
#     seed <- seedlist[s]
#     set.seed(seed)
#     idx_tr <- sample.int(n, n_tr)
# 
#     # training data
#     coords_tr <- data_hs[idx_tr, c("easting", "northing")]
#     rownames(coords_tr) <- NULL
#     y_tr <- data_hs[idx_tr,"y"]
#     ybar <- mean(y_tr)
# 
#     # test data
#     coords_tt <- data_hs[-idx_tr, c("easting", "northing")]
#     rownames(coords_tt) <- NULL
#     y_tt <- data_hs[-idx_tr,"y"]
# 
#     # preprocessing
#     time_pre_pde <- system.time({
#       # create mesh
#       mesh_pde <- create.mesh.2D(nodes = rbind(hs_nodes, coords_tr, coords_tt),
#                                  segments = hs_segments)
# 
#       # create the FEM basis
#       fem_basis <- create.FEM.basis(mesh_pde)
#     })
# 
#     time_pde <- system.time({
#       # fit SR-PDE
#       pde <- smooth.FEM(locations = coords_tr, 
#                         observations = y_tr - ybar,                             # mean-center y
#                         FEMbasis = fem_basis, lambda = lambda, GCV = TRUE)
#     })
# 
#     # best lambda
#     best <- which.min(pde$GCV)
#     lambda_tmp[s] <- lambda[best]
# 
#     # time
#     timepre_tmp[s,] <- as.numeric(time_pre_pde)[1:3]
#     timeres_tmp[s,] <- as.numeric(time_pde)[1:3]
# 
#     # summary results
#     tau <- pde$stderr[best]
#     y_pde <- ybar + 
#       as.numeric(pde$fit.FEM$coeff[(nrow(hs_nodes) + n_tr) + 1:n_tt, best])
#     y_pde_low <- y_pde + qnorm(.025)*tau
#     y_pde_high <- y_pde + qnorm(.975)*tau
#     predres_tmp[s,] <- data.frame(y_tt = y_tt, yhat = y_pde,
#                                   lower = y_pde_low, upper = y_pde_high) %>%
#       mutate(error = y_tt - yhat,
#              width = upper-lower,
#              cover = ifelse((y_tt > lower & y_tt < upper), 1, 0)) %>%
#       summarise(rmspe = sqrt(mean(error^2)),
#                 mape = mean(abs(error)),
#                 coverage = mean(cover),
#                 meanwidth = mean(width)) %>%
#       as.numeric()
# 
#     ## tausq
#     paramres_tmp[s,3] <- tau^2
# 
#     cat(paste0(listname, ": ", s, "th simulation completed.\n"))
#   }
# 
#   lambda_pde[[listname]] <- lambda_tmp
#   timepre[[listname]] <- timepre_tmp
#   timeres[[listname]] <- timeres_tmp
#   predres[[listname]] <- predres_tmp
#   paramres[[listname]] <- paramres_tmp
# }
# 
# saveRDS(list(lambda_pde = lambda_pde, timepre = timepre,
#              timeres = timeres, predres = predres,
#              paramres = paramres, seedmat = seedmat),
#         paste0(path, "sim/horseshoe_pde.RDS"))

# # w plot: similar to Figure 6(b) in https://academic.oup.com/jrsssb/article/75/4/681/7075920
# what_pde <- as.numeric(pde$fit.FEM$coeff[nrow(hs_nodes) + 1:n_tr, best])
# wstar_pde <- as.numeric(pde$fit.FEM$coeff[(nrow(hs_nodes) + n_tr) + 1:n_tt, best])
# 
# tibble(rbind(coords_tr, coords_tt), what = c(what_pde, wstar_pde)) %>%
#   ggplot() +
#   geom_raster(aes(easting, northing, fill = what)) +
#   geom_sf(data = hs_sf, fill = NA) +
#   geom_contour(aes(x = easting, y = northing, z = what),
#                color = "black", breaks = seq(-4, 4, by=.5)) +
#   scale_fill_scico(palette = "roma", direction =-1) +
#   labs(x="", y="", fill = "w") +
#   theme(plot.margin = margin(t = -10, l = -10, r = 0, b = -20), aspect.ratio = 1/2,
#         legend.margin = margin(b = 0, r = 0, t = 0, l = -5))

#---------------------# 
# ClosePD and MDSdist #
#---------------------# 
# link 
# https://link.springer.com/article/10.1007/s11004-019-09791-y
# https://onlinelibrary.wiley.com/doi/10.1002/env.588

# # Non-Euclidean water distances
# # 1. create a domain as a TransitionLayer (converted from a RasterLayer)
# nrows <- 150
# ncols <- 300
# hs_rl <- raster(xmn = -0.9, xmx = 3.4, ymn = -0.9, ymx = 0.9,
#                 nrows = nrows, ncols = ncols)
# xgrid <- seq(-0.9, 3.4, length = ncols)
# ygrid <- seq(-0.9, 0.9, length = nrows)
# xxgrid <- rep(xgrid, nrows)
# yygrid <- rep(ygrid, rep(ncols, nrows))
# zgrid <- mgcv::fs.test(xxgrid, yygrid)
# zgrid[!is.na(zgrid)] <- 1
# values(hs_rl) <- as.vector(zgrid)
# plot(hs_rl)
# hs_tl <- transition(as.factor(hs_rl),
#                     transitionFunction = "areas", directions = 8)
# hs_tl <- geoCorrection(hs_tl, type="c")
# 
# # 2. water distance: costDistance(TransitionLayer, obs coords)
# time_pre_all <- system.time({
#   waterdist_all <- costDistance(hs_tl,
#                                 as.matrix(data_hs[, c("easting", "northing")]))
#   waterdist_max <- max(waterdist_all[waterdist_all < Inf])
# 
#   ## some are Inf with which eigenvalues are not computable for MDS-transformed distance
#   ## Replace Inf with maximum value
#   waterdist_all2 <- waterdist_all
#   waterdist_all2[waterdist_all2 == Inf] <- waterdist_max
# 
#   # MDS-transformed distance
#   mds_all <- cmdscale(waterdist_all2, k = n-1, eig = T, x.ret = TRUE)
#   # caluclating distance between MDS points (automatically limited to positive eigenvalues)
#   distMDS_all <- dist(mds_all$points)
# })
# 
# # save results
# timepre = predres = paramres = timeres <- list()
# timepre_mds = predres_mds = paramres_mds = timeres_mds <- list()
# 
# for (i in 1:length(ntrlist)) {
#   n_tr <- ntrlist[i]
#   n_tt <- n - n_tr
# 
#   # varying m is meaningless
#   # use seed for each n at m = 10
#   m <- 10
#   listname <- paste0("n=",n_tr)
# 
#   # save results
#   timeres_tmp = timepre_tmp = 
#     timeres_mds_tmp = timepre_mds_tmp <- matrix(0, nrow = rplc, ncol = 3)
#   colnames(timeres_tmp) = colnames(timepre_tmp) = 
#     colnames(timeres_mds_tmp) = colnames(timepre_mds_tmp) <-
#     c("user", "system", "elapsed")
#   predres_tmp = paramres_tmp = 
#     predres_mds_tmp = paramres_mds_tmp <- matrix(NA, nrow = rplc, ncol = 4)
#   colnames(predres_tmp) = colnames(predres_mds_tmp) <- 
#     c("rmspe", "mape", "coverage", "meanwidth")
#   colnames(paramres_tmp) = colnames(paramres_mds_tmp) <- 
#     c("beta", "sigsq", "tausq", "phi")
#   seedlist <- seedmat %>% filter(M == m, N_tr == n_tr) %>%
#     dplyr::select(seed) %>% unlist() %>% as.numeric()
# 
#   for (s in 1:rplc) {
#     seed <- seedlist[s]
#     set.seed(seed)
#     idx_tr <- sample.int(n, n_tr)
# 
#     # training data
#     coords_tr <- data_hs[idx_tr, c("easting", "northing")]
#     rownames(coords_tr) <- NULL
#     y_tr <- data_hs[idx_tr,"y"]
#     data_tr <- tibble(coords_tr, y = y_tr)
# 
#     # test data
#     coords_tt <- data_hs[-idx_tr, c("easting", "northing")]
#     rownames(coords_tt) <- NULL
#     y_tt <- data_hs[-idx_tr,"y"]
#     data_tt <- tibble(coords_tt, y = y_tt)
# 
#     # preprocessing
#     time_pre_tr <- system.time({
#       # estimate waterdistance
#       waterdist <- costDistance(hs_tl, as.matrix(coords_tr))
#       
#       ## some are Inf with which eigenvalues are not computable for MDS-transformed distance
#       ## Replace Inf with maximum value 
#       waterdist2 <- waterdist
#       waterdist2[waterdist2 == Inf] <- waterdist_max
#       
#       # MDS-transformed distance
#       mds <- cmdscale(waterdist2, k = n_tr-1, eig = T, x.ret = TRUE)
#       # caluclating distance between MDS points (automatically limited to positive eigenvalues)
#       distMDS <- dist(mds$points)                      
#     })
#     
#     data_tr_geo <- as.geodata(data_tr, coords.col = 1:2, data.col = 3)
#     
#     #---------# 
#     # ClosePD #
#     #---------# 
#     time_closePD_vario <- system.time({
#       # calculate variogram and semivariogram
#       variog_NEuc <- variogNE(geodata = data_tr_geo, 
#                               max.dist = max(waterdist2),
#                               option = "bin", NonEuc = TRUE, messages = FALSE,
#                               NonEucDist = waterdist2)
#       
#       ## using classic WLS to fit semivariogram function
#       vfit_NEuc <- variofit(variog_NEuc, cov.model = "matern", 
#                             fix.kappa = TRUE, kappa = 1, messages = FALSE)
#       
#       semi_NEuc <- semivarioNE(vario = vfit_NEuc, dist = waterdist_all2, 
#                                model = "matern", MDSgamma = FALSE)
#       
#       # transform semivariance matrix of Water distances to covariance matrix 
#       cov_NEuc <- max(semi_NEuc[[1]]) - semi_NEuc[[1]]
#       # eigen(cov_NEuc)$values # identifying negative eigen values
#       
#       # re-estimate covariance matrix using new eigen values based on the tolerance 
#       cov_closePD <- posdefify(cov_NEuc, eps.ev = 1e-07)                              # eps.ev = 1/tau in ClosePD notation
#       # eigen(cov_closePD)$values
#     })
#     
#     time_closePD <- system.time({
#       pred_closePD <- krigeNE(geodata = data_tr_geo, locations = coords_tt, 
#                               covariance = cov_closePD, 
#                               krige = krige.control(obj.model = vfit_NEuc), 
#                               cv = -idx_tr)
#     })
# 
#     # time
#     timepre_tmp[s,] <- as.numeric(time_pre_all + time_pre_tr + 
#                                     time_closePD_vario)[1:3]
#     timeres_tmp[s,] <- as.numeric(time_closePD)[1:3]
# 
#     # summary results
#     tau_closePD <- sqrt(pred_closePD$krige.var)
#     y_closePD <- pred_closePD$predict
#     y_closePD_low <- y_closePD + qnorm(.025)*tau_closePD
#     y_closePD_high <- y_closePD + qnorm(.975)*tau_closePD
#     predres_tmp[s,] <- data.frame(y_tt = y_tt, yhat = y_closePD,
#                                   lower = y_closePD_low, 
#                                   upper = y_closePD_high) %>%
#       mutate(error = y_tt - yhat,
#              width = upper-lower,
#              cover = ifelse((y_tt > lower & y_tt < upper), 1, 0)) %>%
#       summarise(rmspe = sqrt(mean(error^2)),
#                 mape = mean(abs(error)),
#                 coverage = mean(cover),
#                 meanwidth = mean(width)) %>%
#       as.numeric()
# 
#     ## covariance parameters, may not be positive definite
#     paramres_tmp[s,2] <- vfit_NEuc$cov.pars[1]
#     paramres_tmp[s,3] <- vfit_NEuc$nugget
#     paramres_tmp[s,4] <- 1/vfit_NEuc$cov.pars[2]
#     
#     #---------# 
#     # MDSdist #
#     #---------# 
#     time_MDSdist_vario <- system.time({
#       # calculate variogram and semivariogram
#       variog_MDSdist <- variogNE(geodata = data_tr_geo, 
#                                  max.dist = max(distMDS),
#                                  option = "bin", NonEuc = TRUE, messages = FALSE,
#                                  NonEucDist = distMDS)
#       
#       ## using classic WLS to fit semivariogram function
#       vfit_MDSdist <- variofit(variog_MDSdist, cov.model = "matern", 
#                                fix.kappa = TRUE, kappa = 1, messages = FALSE)
#       
#       semi_MDSdist <- semivarioNE(vario = vfit_MDSdist, dist = distMDS_all, 
#                                   model = "matern", MDSgamma = FALSE)
#       
#       # transform semivariance matrix of water distances to covariance matrix 
#       cov_MDSdist <- max(semi_MDSdist[[1]]) - semi_MDSdist[[1]]
#     })
#     
#     time_MDSdist <- system.time({
#       pred_MDSdist <- krigeNE(geodata = data_tr_geo, locations = coords_tt, 
#                               covariance = cov_MDSdist, 
#                               krige = krige.control(obj.model = vfit_MDSdist), 
#                               cv = -idx_tr)
#     })
#     
#     # time
#     timepre_mds_tmp[s,] <- as.numeric(time_pre_all + time_pre_tr + 
#                                         time_MDSdist_vario)[1:3]
#     timeres_mds_tmp[s,] <- as.numeric(time_MDSdist)[1:3]
#     
#     # summary results
#     tau_mds <- sqrt(pred_MDSdist$krige.var)
#     y_mds <- pred_MDSdist$predict
#     y_mds_low <- y_mds + qnorm(.025)*tau_mds
#     y_mds_high <- y_mds + qnorm(.975)*tau_mds
#     predres_mds_tmp[s,] <- data.frame(y_tt = y_tt, yhat = y_mds,
#                                       lower = y_mds_low, upper = y_mds_high) %>%
#       mutate(error = y_tt - yhat,
#              width = upper-lower,
#              cover = ifelse((y_tt > lower & y_tt < upper), 1, 0)) %>%
#       summarise(rmspe = sqrt(mean(error^2)),
#                 mape = mean(abs(error)),
#                 coverage = mean(cover),
#                 meanwidth = mean(width)) %>%
#       as.numeric()
#     
#     ## covariance parameters
#     paramres_mds_tmp[s,2] <- vfit_MDSdist$cov.pars[1]
#     paramres_mds_tmp[s,3] <- vfit_MDSdist$nugget
#     paramres_mds_tmp[s,4] <- 1/vfit_MDSdist$cov.pars[2]
#     
#     cat(paste0(listname, ": ", s, "th simulation completed.\n"))
#   }
# 
#   timepre[[listname]] <- timepre_tmp
#   timeres[[listname]] <- timeres_tmp
#   predres[[listname]] <- predres_tmp
#   paramres[[listname]] <- paramres_tmp
#   timepre_mds[[listname]] <- timepre_mds_tmp
#   timeres_mds[[listname]] <- timeres_mds_tmp
#   predres_mds[[listname]] <- predres_mds_tmp
#   paramres_mds[[listname]] <- paramres_mds_tmp
# }
# 
# saveRDS(list(timepre = timepre,
#              timeres = timeres, predres = predres,
#              paramres = paramres, seedmat = seedmat),
#         paste0(path, "sim/horseshoe_closePD.RDS"))
# 
# saveRDS(list(timepre = timepre_mds,
#              timeres = timeres_mds, predres = predres_mds,
#              paramres = paramres_mds, seedmat = seedmat),
#         paste0(path, "sim/horseshoe_MDSdist.RDS"))

# call results
resBRGP <- readRDS(paste0(path, "sim/horseshoe_BORAGP.RDS"))
resNNGP <- readRDS(paste0(path, "sim/horseshoe_NNGP.RDS"))
resinla <- readRDS(paste0(path, "sim/horseshoe_barrierSGF.RDS"))
resSOAP <- readRDS(paste0(path, "sim/horseshoe_soap.RDS"))
resPDE <- readRDS(paste0(path, "sim/horseshoe_pde.RDS"))
resMDS <- readRDS(paste0(path, "sim/horseshoe_MDSdist.RDS"))
resPD <- readRDS(paste0(path, "sim/horseshoe_closePD.RDS"))

# results
resnames <- c("RMSPE", "MAPE", "95% CI coverage", "Mean 95% CI width")
parnames <- c("beta", "sigmasq", "tausq", "phi")
dfBRGPL = dfNNGPL = dfinlaL = dfSOAPL = dfPDEL = dfMDSL = dfPDL <- list()
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
                                                 function(x) x[,i])),
                                        times = length(mlist)),
                             what = resnames[i],
                             vals2 = rep(c(sapply(resinla$paramres,
                                                  function(x) x[,i])),
                                         times = length(mlist)),
                             what2 = parnames[i])
  dfSOAPL[[i]] <- data.frame(n = rep(rep(ntrlist, each = rplc),
                                     times = length(mlist)),
                             m = rep(mlist, each = rplc*length(ntrlist)),
                             model = "Soap film smoother",
                             vals = rep(c(sapply(resSOAP$predres,
                                                 function(x) x[,i])),
                                        times = length(mlist)),
                             what = resnames[i],
                             vals2 = rep(c(sapply(resSOAP$paramres,
                                                  function(x) x[,i])),
                                         times = length(mlist)),
                             what2 = parnames[i])
  dfPDEL[[i]] <- data.frame(n = rep(rep(ntrlist, each = rplc),
                                    times = length(mlist)),
                            m = rep(mlist, each = rplc*length(ntrlist)),
                            model = "SR-PDE",
                            vals = rep(c(sapply(resPDE$predres,
                                                function(x) x[,i])),
                                       times = length(mlist)),
                            what = resnames[i],
                            vals2 = rep(c(sapply(resPDE$paramres,
                                                 function(x) x[,i])),
                                        times = length(mlist)),
                            what2 = parnames[i])
  dfMDSL[[i]] <- data.frame(n = rep(rep(ntrlist, each = rplc),
                                    times = length(mlist)),
                            m = rep(mlist, each = rplc*length(ntrlist)),
                            model = "MDSdist",
                            vals = rep(c(sapply(resMDS$predres,
                                                function(x) x[,i])),
                                       times = length(mlist)),
                            what = resnames[i],
                            vals2 = rep(c(sapply(resMDS$paramres,
                                                 function(x) x[,i])),
                                        times = length(mlist)),
                            what2 = parnames[i])
  dfPDL[[i]] <- data.frame(n = rep(rep(ntrlist, each = rplc),
                                   times = length(mlist)),
                           m = rep(mlist, each = rplc*length(ntrlist)),
                           model = "ClosePD",
                           vals = rep(c(sapply(resPD$predres,
                                               function(x) x[,i])),
                                      times = length(mlist)),
                           what = resnames[i],
                           vals2 = rep(c(sapply(resPD$paramres,
                                                function(x) x[,i])),
                                       times = length(mlist)),
                           what2 = parnames[i])
}

dfBRGP <- do.call(rbind, dfBRGPL)
dfNNGP <- do.call(rbind, dfNNGPL)
dfinla <- do.call(rbind, dfinlaL)
dfSOAP <- do.call(rbind, dfSOAPL)
dfPDE <- do.call(rbind, dfPDEL)
dfMDS <- do.call(rbind, dfMDSL)
dfPD <- do.call(rbind, dfPDL)
dfres <- rbind(dfBRGP, dfNNGP, dfinla, dfSOAP, dfPDE, dfMDS, dfPD)

## prediction performance
sapply(resBRGP$predres, colMeans) %>% round(3)
sapply(resBRGP$predres, function(x) apply(x, 2, sd)) %>% round(3)
sapply(resNNGP$predres, colMeans) %>% round(3)
sapply(resNNGP$predres, function(x) apply(x, 2, sd)) %>% round(3)
sapply(resinla$predres, colMeans) %>% round(3)
sapply(resinla$predres, function(x) apply(x, 2, sd)) %>% round(3)
sapply(resSOAP$predres, colMeans) %>% round(3)
sapply(resSOAP$predres, function(x) apply(x, 2, sd)) %>% round(3)
sapply(resPDE$predres, colMeans) %>% round(3)
sapply(resPDE$predres, function(x) apply(x, 2, sd)) %>% round(3)
sapply(resMDS$predres, colMeans) %>% round(3)
sapply(resMDS$predres, function(x) apply(x, 2, sd)) %>% round(3)
sapply(resPD$predres, colMeans) %>% round(3)
sapply(resPD$predres, function(x) apply(x, 2, sd)) %>% round(3)

## prediction performance in plots
dfres_sum <- dfres %>%
  group_by(n,m,model,what) %>%
  summarise(vals_mean = mean(vals),
            vals_low = quantile(vals, prob = 0.025),
            vals_hi = quantile(vals, prob = 0.975))
dfres_sum$m[dfres_sum$model == "BORA-GP"] <-
  dfres_sum$m[dfres_sum$model == "BORA-GP"] - 0.5
dfres_sum$m[dfres_sum$model == "SR-PDE"] <-
  dfres_sum$m[dfres_sum$model == "SR-PDE"] + 0.5

sumplot1 <- dfres_sum %>%
  filter(model %in% c("BORA-GP", "NNGP", "Barrier SGF", "SR-PDE")) %>%
  mutate(model = factor(model, levels =
                          c("BORA-GP", "NNGP", "Barrier SGF", "SR-PDE"))) %>%
  filter(what == "RMSPE") %>%
  ggplot() +
  geom_line(aes(m, vals_mean, col = model, linetype = model)) +
  geom_errorbar(aes(m, ymin = vals_low, ymax = vals_hi,
                    col = model, linetype = model), width = 0.5) +
  geom_point(aes(m, vals_mean, col = model), size = 2) +
  facet_grid(what ~
               factor(n, labels = c("n = 300", "n = 600", "n = 1000")), scales = "free") +
  labs(x="", y="", color="", linetype="", shape="") +
  scale_color_scico_d(palette = "batlow", begin = 0.2, end = 0.9, direction = -1) +
  theme(legend.position = "right",
        plot.margin=margin(t=0,l=-10, r=0, b=-10),
        legend.margin=margin(b=0,r=0,t=0,l=-3)) +
  scale_x_continuous(breaks = c(10,15,20), labels = c("m = 10", "m = 15", "m = 20"))

sumplot2 <- dfres %>%
  filter(model %in% c("BORA-GP", "NNGP", "Barrier SGF", "SR-PDE")) %>%
  mutate(model = factor(model, levels =
                          c("BORA-GP", "NNGP", "Barrier SGF", "SR-PDE")),
         m = factor(m, levels = c("10", "15", "20"),
                    labels = c("m = 10", "m = 15", "m = 20"))) %>%
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
#          paste0(path, "plots/horseshoe_pred_wSRPDE", ext),
#          width = 9, height = 4)
# }

sumplot3 <- dfres_sum %>%
  filter(model %in% c("BORA-GP", "NNGP", "Barrier SGF", "SR-PDE"),
         what %in% c("MAPE", "Mean 95% CI width")) %>%
  mutate(model = factor(model, levels = c("BORA-GP", "NNGP", "Barrier SGF", "SR-PDE"))) %>%
  ggplot() +
  geom_line(aes(m, vals_mean, col = model, linetype = model)) +
  geom_errorbar(aes(m, ymin = vals_low, ymax = vals_hi,
                    col = model, linetype = model), width = 0.5) +
  geom_point(aes(m, vals_mean, col = model), size = 2) +
  facet_grid(factor(what, levels = c("MAPE", "Mean 95% CI width"),
                    labels = c("MAPE", "95% CI width")) ~
               factor(n, labels = c("n = 300", "n = 600", "n = 1000")), scales = "free") +
  labs(x="", y="", color="", linetype="", shape = "") +
  scale_shape_manual(values = c(19, 17, 15, 8, 12)) +
  scale_color_scico_d(palette = "batlow", begin = 0.2, end = 0.9, direction = -1) +
  theme(legend.position = "right",
        plot.margin=margin(t=0,l=-10, r=0, b=-10),
        legend.margin=margin(b=0,r=0,t=0,l=-3)) +
  scale_x_continuous(breaks = c(10,15,20), labels = c("m = 10", "m = 15", "m = 20"))
# for (ext in extension) {
#   ggsave(plot = sumplot3,
#          paste0(path, "plots/horseshoe_pred_supp_wSRPDE", ext),
#          width = 9, height = 3.5)
# }

## parameter estimation
sapply(resBRGP$paramres, colMeans) %>% round(3)
sapply(resBRGP$paramres, function(x) apply(x, 2, sd)) %>% round(3)
sapply(resNNGP$paramres, colMeans) %>% round(3)
sapply(resNNGP$paramres, function(x) apply(x, 2, sd)) %>% round(3)
sapply(resinla$paramres, colMeans) %>% round(3)
sapply(resinla$paramres, function(x) apply(x, 2, sd)) %>% round(3)
sapply(resSOAP$paramres, colMeans) %>% round(3)
sapply(resSOAP$paramres, function(x) apply(x, 2, sd)) %>% round(3)
sapply(resPDE$paramres, colMeans) %>% round(3)
sapply(resPDE$paramres, function(x) apply(x, 2, sd)) %>% round(3)
sapply(resMDS$paramres, colMeans) %>% round(3)
sapply(resMDS$paramres, function(x) apply(x, 2, sd)) %>% round(3)
sapply(resPD$paramres, colMeans) %>% round(3)
sapply(resPD$paramres, function(x) apply(x, 2, sd)) %>% round(3)

## average elapsed time for BORA-GP neighbor search
sapply(resBRGP$timenb, colMeans) %>% round(3)

## average elapsed time for INLA mesh
sapply(resinla$timemesh, colMeans) %>% round(3)

## average elapsed time for SR-PDE mesh and FEM basis
sapply(resPDE$timepre, colMeans) %>% round(3)

## average elapsed time for MDSdist to estimate distance and covariance
sapply(resMDS$timepre, colMeans) %>% round(3)

## average elapsed time for ClosePD to estimate distance and covariance
sapply(resPD$timepre, colMeans) %>% round(3)

## average elapsed time for model fitting (sec/iter)
(sapply(resBRGP$timeres, colMeans)/n.samples) %>% round(3)
(sapply(resNNGP$timeres, colMeans)/n.samples) %>% round(3)

## average elapsed time for model fitting
sapply(resBRGP$timeres, colMeans) %>% round(3)
sapply(resNNGP$timeres, colMeans) %>% round(3)
sapply(resinla$timeres, colMeans) %>% round(3)
sapply(resSOAP$timeres, colMeans) %>% round(3)
sapply(resPDE$timeres, colMeans) %>% round(3)
sapply(resMDS$timeres, colMeans) %>% round(3)
sapply(resPD$timeres, colMeans) %>% round(3)

#--------------#
# latent model #
#--------------#
# simulate data
n_tr <- ntrlist[1]
n_tt <- n - n_tr
m <- mlist[3]
seedlist <- seedmat %>% filter(M == m, N_tr == n_tr) %>% pull(seed)

set.seed(seedlist[27])
idx_tr <- sample.int(n, n_tr)

# training data
coords_tr <- data_hs[idx_tr, c("easting", "northing")]
rownames(coords_tr) <- NULL
w_tr <- data_hs[idx_tr,"w"]
y_tr <- data_hs[idx_tr,"y"]
ybar <- mean(y_tr)

# ordering (default)
ord <- order(coords_tr[,1])

# test data
coords_tt <- data_hs[-idx_tr, c("easting", "northing")]
rownames(coords_tt) <- NULL
w_tt <- data_hs[-idx_tr,"w"]
y_tt <- data_hs[-idx_tr,"y"]

coords_all <- rbind(coords_tr, coords_tt)

# # run models
# set.seed(seedlist[27])
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
# # mesh
# max.edge <- 0.2
# mesh <- inla.mesh.2d(boundary = poly,
#                      loc = coords_tr,
#                      max.edge = c(1,5)*max.edge,
#                      cutoff = 0.04,
#                      offset = c(max.edge, 1.5))
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
# # soap film smoother
# data_tr <- tibble(coords_tr, y = y_tr)
# data_tt <- tibble(coords_tt, y = y_tt)
#
# # knots
# test <- mgcv::fs.test(x = coords_tr$easting, y = coords_tr$northing,
#                       r0 = 0.15)
# knots <- coords_tr[-which(is.na(test)),]
#
# # model fitting
# time_soap <- system.time({
#   soap <- gam(y ~ s(easting, northing, bs = "so", xt = list(bnd = hs_soap)),
#               data = data_tr, method = "REML", knots = knots)
# })
#
# # SR-PDE
# # preprocessing
# time_pre_pde <- system.time({
#   # create mesh
#   mesh_pde <- create.mesh.2D(nodes = rbind(hs_nodes, coords_tr, coords_tt),
#                              segments = hs_segments)
#
#   # create the FEM basis
#   fem_basis <- create.FEM.basis(mesh_pde)
# })
#
# time_pde <- system.time({
#   # fit SR-PDE
#   pde <- smooth.FEM(locations = coords_tr, observations = y_tr - ybar,
#                     FEMbasis = fem_basis, lambda = lambda, GCV = TRUE)
# })
# best <- which.min(pde$GCV)
#
# # MDSdist (since tau^2 is estimated to be 0)
# # estimate waterdistance
# waterdist <- costDistance(hs_tl, as.matrix(coords_tr))
#
# ## some are Inf with which eigenvalues are not computable for MDS-transformed distance
# ## Replace Inf with maximum value
# waterdist2 <- waterdist
# waterdist2[waterdist2 == Inf] <- waterdist_max
#
# # MDS-transformed distance
# mds <- cmdscale(waterdist2, k = n_tr-1, eig = T, x.ret = TRUE)
# # caluclating distance between MDS points (automatically limited to positive eigenvalues)
# distMDS <- dist(mds$points)
#
# data_tr_geo <- as.geodata(tibble(coords_tr, y = y_tr),
#                           coords.col = 1:2, data.col = 3)
# variog_MDSdist <- variogNE(geodata = data_tr_geo,
#                            max.dist = max(distMDS),
#                            option = "bin", NonEuc = TRUE, messages = FALSE,
#                            NonEucDist = distMDS)
# vfit_MDSdist <- variofit(variog_MDSdist, cov.model = "matern",
#                          fix.kappa = TRUE, kappa = 1)
# semi_MDSdist <- semivarioNE(vario = vfit_MDSdist,
#                             dist = distMDS_all,
#                             model = "matern", MDSgamma = FALSE)
# cov_MDSdist <- max(semi_MDSdist[[1]]) - semi_MDSdist[[1]]
#
# MDSdist <- krigeNE(geodata = data_tr_geo, locations = coords_tt,
#                    covariance = cov_MDSdist,
#                    krige = krige.control(obj.model = vfit_MDSdist),
#                    cv = -idx_tr)
#
# # ClosePD (since tau^2 is estimated to be 0)
# variog_NEuc <- variogNE(geodata = data_tr_geo,
#                         max.dist = max(waterdist2),
#                         option = "bin", NonEuc = TRUE, messages = FALSE,
#                         NonEucDist = waterdist2)
# vfit_NEuc <- variofit(variog_NEuc, cov.model = "matern",
#                       fix.kappa = TRUE, kappa = 1)
# semi_NEuc <- semivarioNE(vario = vfit_NEuc,
#                          dist = waterdist_all2,
#                          model = "matern", MDSgamma = FALSE)
# cov_NEuc <- max(semi_NEuc[[1]]) - semi_NEuc[[1]]
# cov_closePD <- posdefify(cov_NEuc, eps.ev = 1e-07)
#
# closePD <- krigeNE(geodata = data_tr_geo, locations = coords_tt,
#                    covariance = cov_closePD,
#                    krige = krige.control(obj.model = vfit_NEuc),
#                    cv = -idx_tr)
#
# saveRDS(list(barrier_m.s = barrier_m.s, barrier_p.s = barrier_p.s,
#              m.s = m.s, p.s = p.s, res = res,
#              pde = pde, soap = soap,
#              MDSdist = MDSdist,
#              closePD = closePD),
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
# whatSOAP <- soap$fitted.values - soap$coefficients[1]
# wstarSOAP <- as.numeric(predict(soap, newdata = data_tt, type = "link")) -
#   soap$coefficients[1]
# whatPDE <- as.numeric(pde$fit.FEM$coeff[nrow(hs_nodes) + 1:n_tr, best])
# wstarPDE <- as.numeric(pde$fit.FEM$coeff[(nrow(hs_nodes) + n_tr) + 1:n_tt, best])
# wstarMDS <- MDSdist$predict
# wstarPD <- closePD$predict
#
# saveRDS(list(whatBRGP = whatBRGP,
#              wstarBRGP = wstarBRGP,
#              whatNNGP = whatNNGP,
#              wstarNNGP = wstarNNGP,
#              winla = winla,
#              whatSOAP = whatSOAP,
#              wstarSOAP = wstarSOAP,
#              whatPDE = whatPDE,
#              wstarPDE = wstarPDE,
#              wstarMDS = wstarMDS,
#              wstarPD = wstarPD),
#         paste0(path, "sim/horseshoe_latent_res.RDS"))

# horseshoe plots
latent_res <- readRDS(paste0(path, "sim/horseshoe_latent_res.RDS"))
whatBRGP <- latent_res$whatBRGP
wstarBRGP <- latent_res$wstarBRGP
whatNNGP <- latent_res$whatNNGP
wstarNNGP <- latent_res$wstarNNGP
winla <- latent_res$winla
whatSOAP <- latent_res$whatSOAP
wstarSOAP <- latent_res$wstarSOAP
whatPDE <- latent_res$whatPDE
wstarPDE <- latent_res$wstarPDE
wstarMDS <- latent_res$wstarMDS
wstarPD <- latent_res$wstarPD

gg_w <- data.frame(easting = rep(coords_all$easting, times = 8),
                   northing = rep(coords_all$northing, times = 8),
                   what = c(w_tr, w_tt,
                            whatNNGP, wstarNNGP,
                            whatBRGP, wstarBRGP,
                            winla,
                            whatSOAP, wstarSOAP,
                            whatPDE, wstarPDE,
                            w_tr, wstarMDS,
                            w_tr, wstarPD),
                   model = rep(c("Truth", "NNGP", "BORA-GP", "Barrier SGF",
                                 "Soap film smoother", "SR-PDE",
                                 "MDSdist", "ClosePD"), each = n)) %>%
  ggplot() +
  facet_wrap(~factor(model,
                     levels = c("Truth", "NNGP", "BORA-GP", "Barrier SGF",
                                "Soap film smoother", "SR-PDE",
                                "MDSdist", "ClosePD")),
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
#          paste0(path, "plots/horseshoe_w_more", ext),
#          width = 12, height = 3.6)
# }