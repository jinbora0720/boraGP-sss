# Faults
rm(list = ls())

# dependencies
library(tidyverse)
theme_set(theme_bw())
library(sf)
library(fields)
library(rgeos)
library(rgdal)
library(boraGP)
# install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(scico)
library(INLA)

# path 
path <- "~/boraGP-sss/"

# plot extensions
extension <- c(".pdf", ".png")

# create faults polygons
# output is a square
local.square.polygon = function(xlim, ylim){
  xlim = range(xlim); ylim = range(ylim)
  corner1 = c(xlim[1], ylim[2])
  corner2 = c(xlim[2], ylim[1])
  poly = Polygon(rbind(corner1, c(corner1[1], corner2[2]), corner2, 
                       c(corner2[1], corner1[2]), corner1), hole = FALSE)
  return(SpatialPolygons(list(Polygons(list(poly), ID = runif(1)))))
}

# the width/thickness of the barrier
width <- 0.002 # 20 times thicker, width/2 = 0.001
poly1 <- local.square.polygon(xlim = c(0.3, 2), 
                              ylim = 0.6 + width*c(-1, 1))
poly2 <- local.square.polygon(xlim = c(0.7, 2), 
                              ylim = 1.5 + width*c(-1, 1))
poly.original <- SpatialPolygons(c(poly1@polygons, poly2@polygons))
fault <- st_as_sf(poly.original)
faults <- st_make_valid(st_combine(fault$geometry))
bbox <- st_bbox(c(xmin = 0, xmax = 2, ymin = 0, ymax = 2))

# create test points
ngrid <- 67
locs <- expand.grid(seq(0, 2, length = ngrid), seq(0, 2, length = ngrid))
colnames(locs) <- c("easting", "northing")

# create training points
locs_tr <- expand.grid(seq(0, 2, length = 12), seq(0, 2, length = 34))
colnames(locs_tr) <- c("easting", "northing")

# put together 
coords_all <- unique(rbind(locs_tr, locs))
n <- nrow(coords_all)

# ordering for all locations
ord_all <- order(coords_all[,2])

# neighbor 
m <- 15

# training vs. testing 
n_tr <- nrow(locs_tr)
n_tt <- n-n_tr
coords_tr <- coords_all[1:n_tr, ]
coords_tt <- coords_all[-(1:n_tr), ]

# ordering for training locations
ord <- order(coords_tr[,2])

# #---------------#
# # true function #
# #---------------#
# # barrier neighbor info for all as training
# barrier_neighbor_time <- system.time(
#   barrier_nninfo_truth <- barrier_neighbor(coords = coords_all,
#                                            coords.0 = NULL,
#                                            ord = ord_all,
#                                            n.neighbors = m,
#                                            barrier = faults,
#                                            cores = 20,
#                                            verbose = TRUE,
#                                            debug = list(barrier_n.indx = NULL,
#                                                         barrier_dist = NULL,
#                                                         barrier_nn.indx.0_list = NULL,
#                                                         barrier_dist0 = NULL,
#                                                         ref_window = NULL,
#                                                         nonref_window = NULL,
#                                                         ref_fill = TRUE,
#                                                         nonref_fill = TRUE))
# )
# 
# barrier_n.indx_truth <- barrier_nninfo_truth$barrier_n.indx
# barrier_nn.indx_truth <- as.integer(unlist(barrier_n.indx_truth)[-1]-1)
# barrier_nn.indx.lu_truth <- c(cumsum(c(0,0,sapply(barrier_n.indx_truth, length)[-1]))[1:n],
#                               c(0,sapply(barrier_n.indx_truth, length)[-1])) %>%
#   as.integer()
# barrier_nninfo_truth <- list(type = "barrier",
#                              n.indx = barrier_n.indx_truth,
#                              n.neighbors = m,
#                              nn.indx = barrier_nn.indx_truth,
#                              nn.indx.lu = barrier_nn.indx.lu_truth,
#                              ord = ord_all)
# 
# # barrier neighbor info for training vs. test w/ m
# time_nb <- system.time(
#   barrier_nninfo_all <- barrier_neighbor(coords = coords_tr,
#                                          coords.0 = coords_tt,
#                                          ord = ord,
#                                          n.neighbors = m,
#                                          barrier = faults,
#                                          cores = 20,
#                                          verbose = TRUE,
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
# # true parameters
# true_sig_sq <- 1
# true_nu <- 1.5
# true_phi <- 4
# true_tau_sq <- 0.1
# true_beta <- c(1, .5)
# 
# #------------------#
# # MCMC preparation #
# #------------------#
rplc <- 30
# set.seed(48)
# seedsave <- c(1, sample.int(1e3, size = rplc-1)) # 113, 115, 120
# 
# # mcmc
# n.samples <- 15000
# burn <- 10000
# 
# # priors
# sigma.sq <- true_sig_sq
# tau.sq <- true_tau_sq
# nu <- true_nu
# phi <- true_phi
# 
# ## matern
# maxdist <- max(dist(as.matrix(coords_tr[,c("easting", "northing")])))
# d.low <- 0.25*maxdist
# d.high <- 0.75*maxdist
# dist_mat <- matrix(c(0, d.low, d.high, d.low, 0, 0, d.high, 0, 0),
#                    nrow = 3, ncol = 3, byrow = T)
# phi.cand <- seq(0.1, 10, by=0.01)
# cor.low = cor.high <- rep(0, length(phi.cand))
# for (i in 1:length(phi.cand)) {
#   cor <- Cov_matern(dist = dist_mat, sigmasq = 1, phi = phi.cand[i], nu = nu)
#   cor.low[i] <- cor[1,2]
#   cor.high[i] <- cor[1,3]
# }
# phi.high <- phi.cand[which.min(abs(cor.low - 0.05))]
# phi.low <- phi.cand[which.min(abs(cor.high - 0.05))]
# 
# starting <- list("phi" = phi, "sigma.sq" = sigma.sq,
#                  "tau.sq" = tau.sq, "nu" = nu)
# tuning <- list("phi" = 0.5, "sigma.sq" = 0.1, "tau.sq" = 0.1, "nu" = 0)         # nu is fixed
# priors <- list("phi.Unif" = c(phi.low, phi.high),
#                "sigma.sq.IG" = c(2, sigma.sq), "tau.sq.IG" = c(2, tau.sq),
#                "nu.Unif" = c(nu-0.5,nu+0.5))
# 
# # mesh
# max.edge <- 0.2
# mesh <- inla.mesh.2d(loc = coords_tr,
#                      interior = inla.sp2segment(sf::as_Spatial(faults)),
#                      max.edge = max.edge,
#                      offset = 0.1)
# 
# ## barrier model
# tl <- length(mesh$graph$tv[,1])
# posTri <- matrix(0, tl, 2)
# for (t in 1:tl){
#   temp <- mesh$loc[mesh$graph$tv[t, ], ]
#   posTri[t,] <- colMeans(temp)[c(1,2)]
# }
# posTri <- SpatialPoints(posTri)
# normal <- over(sf::as_Spatial(faults), SpatialPoints(posTri), returnList = T)
# barrier.triangles <- unlist(normal)
# poly.barrier <- inla.barrier.polygon(mesh, barrier.triangles)
# # plot(poly.barrier, col = "gray")
# # plot(mesh, add=T)
# 
# ## connect observations to mesh nodes
# A.obs <- inla.spde.make.A(mesh, loc = as.matrix(coords_tr))
# ## same for prediction
# proj.pred <- inla.mesh.projector(mesh, loc = as.matrix(coords_tt))
# A.pred <- inla.spde.make.A(mesh, loc = proj.pred$loc)
# 
# barrier.model <- inla.barrier.pcmatern(mesh, barrier.triangles = barrier.triangles,
#                                        prior.range = c(0.9, .9), # P(range < 0.9) = 0.9
#                                        prior.sigma = c(1.5, 0.1)) # P(sigma > 1.5) = 0.1
# formula <- y~ -1 + int + x + f(s, model=barrier.model)
# 
# #------------------#
# # rplc simulations #
# #------------------#
# params_tab_list = w_ressum_tab_list = ressum_tab_list =
#   longpreddata_w_list = longpreddata_list <- list()
# 
# for (ii in 1:rplc) {
# 
#   set.seed(seedsave[ii])
#   x <- rnorm(n)
#   data0 <- rboraGP(coords = coords_all,
#                    neighbor.info = barrier_nninfo_truth,
#                    X = cbind(1, x),
#                    beta = true_beta,
#                    sig_sq = true_sig_sq,
#                    phi = true_phi,
#                    nu = true_nu,
#                    tau_sq = true_tau_sq,
#                    base_cov = "matern")
#   y_tr <- data0$y[1:n_tr]
#   y_tt <- data0$y[-(1:n_tr)]
#   w_tr <- data0$w[1:n_tr]
#   w_tt <- data0$w[-(1:n_tr)]
#   x_tr <- x[1:n_tr]
#   x_tt <- x[-(1:n_tr)]
# 
#   #---------#
#   # BORA-GP #
#   #---------#
#   ## latent model
#   barrier_m.s <- spNNGP(y_tr ~ x_tr, coords = coords_tr, starting = starting,
#                         method = "latent", n.neighbors = m,
#                         tuning = tuning, priors = priors, cov.model = "matern",
#                         n.samples = n.samples, n.omp.threads = 10,
#                         neighbor.info = barrier_nninfo, verbose = F)
#   barrier_p.s <- predict(barrier_m.s, X.0 = cbind(1, x_tt),
#                          coords.0 = as.matrix(coords_tt),
#                          sub.sample = list(start = burn+1, end = n.samples, thin = 1),
#                          nn.indx.0 = barrier_nn.indx.0, n.omp.threads = 10, verbose = F)
# 
#   #------#
#   # NNGP #
#   #------#
#   ## latent model
#   m.s <- spNNGP(y_tr ~ x_tr, coords = coords_tr, starting = starting,
#                 method = "latent", n.neighbors = m,
#                 tuning = tuning, priors = priors, cov.model = "matern",
#                 n.samples = n.samples, n.omp.threads = 10, verbose = F)
#   p.s <- predict(m.s, X.0 = cbind(1, x_tt),
#                  coords.0 = as.matrix(coords_tt),
#                  sub.sample = list(start = burn+1, end = n.samples, thin = 1),
#                  n.omp.threads = 10, verbose = F)
# 
#   #-------------#
#   # Barrier SGF #
#   #-------------#
#   stk.obs <- inla.stack(data = list(y = y_tr),
#                         effects = list(s = 1:mesh$n, # spatial random effects
#                                        data.frame(int = rep(1,n_tr), x = x_tr)),
#                         A = list(A.obs, 1),
#                         remove.unused = FALSE, tag = "obs")
#   stk.pred <- inla.stack(data = list(y = NA),
#                          A = list(A.pred, 1),
#                          effects = list(s = 1:mesh$n,
#                                         data.frame(int = rep(1,n_tt), x = x_tt)),
#                          tag = "pred")
#   stk <- inla.stack(stk.obs, stk.pred)
#   timeinla <- system.time({
#   inlares <- inla(formula,
#                   data=inla.stack.data(stk),
#                   control.predictor=list(A=inla.stack.A(stk), compute = TRUE),
#                   control.compute=list(return.marginals.predictor=TRUE),
#                   family="gaussian",
#                   control.inla= list(int.strategy = "eb"), num.threads = 10, verbose = F)
#   })
# 
#   #---------#
#   # Results #
#   #---------#
#   # time
#   ## boraGP
#   timeBRGP <- as.numeric(barrier_m.s$run.time)[1:3] + as.numeric(barrier_p.s$run.time)[1:3]
# 
#   ## NNGP
#   timeNNGP <- as.numeric(m.s$run.time)[1:3] + as.numeric(p.s$run.time)[1:3]
# 
#   # parameter estimation
#   ## boraGP
#   betaBRGP <- colMeans(barrier_m.s$p.beta.samples)
#   thetaBRGP <- colMeans(barrier_m.s$p.theta.samples[,1:3])
# 
#   ## NNGP
#   betaNNGP <- colMeans(m.s$p.beta.samples)
#   thetaNNGP <- colMeans(m.s$p.theta.samples[,1:3])
# 
#   ## Barrier SGF
#   betainla <- inlares$summary.fixed[,"mean"]
#   thetainla <- rep(0, 3) # ref: https://haakonbakkagit.github.io/btopic103.html
#   ### sigma = exp(theta1)
#   thetainla[1] <- inlares$internal.marginals.hyperpar[2] %>%
#     lapply(function(m) {
#       inla.tmarginal(function(x) exp(x)^2, m)
#     }) %>%
#     sapply(function(m)
#       unlist(inla.zmarginal(m, silent = TRUE))[1])
#   ### tausq
#   thetainla[2] <- inlares$internal.marginals.hyperpar[1] %>%
#     lapply(function(m) {
#       inla.tmarginal(function(x) 1/exp(x), m)
#     }) %>%
#     sapply(function(m)
#       unlist(inla.zmarginal(m, silent = TRUE))[1])
#   ### phi = sqrt(8)/r where spatial range r = exp(theta2)
#   thetainla[3] <- inlares$internal.marginals.hyperpar[3] %>%
#     lapply(function(m) {
#       inla.tmarginal(function(x) sqrt(8)/exp(x), m)
#     }) %>%
#     sapply(function(m)
#       unlist(inla.zmarginal(m, silent = TRUE))[1])
# 
#   params_tab <- matrix(c(betaBRGP, thetaBRGP, timeBRGP,
#                          betaNNGP, thetaNNGP, timeNNGP,
#                          betainla, thetainla, as.numeric(timeinla)[1:3]), nrow = 8)
#   colnames(params_tab) <- c("BORA-GP: latent, m=15",
#                             "NNGP: latent, m=15",
#                             "Barrier SGF")
#   rownames(params_tab) <- c(paste0("beta0=",true_beta[1]),
#                             paste0("beta1=",true_beta[2]),
#                             paste0("sigsq=",true_sig_sq),
#                             paste0("tausq=",true_tau_sq),
#                             paste0("phi=",true_phi),
#                             "user", "system", "elapsed")
#   params_tab_list[[ii]] <- params_tab %>% t()
# 
#   # prediction
#   ## boraGP
#   whatBRGP <- rowMeans(barrier_m.s$p.w.samples[,-c(1:burn)])
#   wstarBRGP <- rowMeans(barrier_p.s$p.w.0)
#   wquantBRGP <- apply(barrier_p.s$p.w.0, 1, function(x) quantile(x, probs = c(0.025, 0.975)))
#   ystarBRGP <- rowMeans(barrier_p.s$p.y.0)
#   yquantBRGP <- apply(barrier_p.s$p.y.0, 1, function(x) quantile(x, probs = c(0.025, 0.975)))
# 
#   ## NNGP
#   whatNNGP <- rowMeans(m.s$p.w.samples[,-c(1:burn)])
#   wstarNNGP <- rowMeans(p.s$p.w.0)
#   wquantNNGP <- apply(p.s$p.w.0, 1, function(x) quantile(x, probs = c(0.025, 0.975)))
#   ystarNNGP <- rowMeans(p.s$p.y.0)
#   yquantNNGP <- apply(p.s$p.y.0, 1, function(x) quantile(x, probs = c(0.025, 0.975)))
# 
#   ## inla
#   index.obs <- inla.stack.index(stk, "obs")$data
#   index.pred <- inla.stack.index(stk, "pred")$data
# 
#   y_inla <- inlares$summary.fitted.values[index.pred, "mean"]
#   inla_sd <- sqrt(inlares$summary.fitted.values[index.pred,"sd"]^2 # variability from the coefficient
#                   + 1/inlares$summary.hyperpar[1,"mean"]) # tausq
#   y_inla_low <- y_inla + qnorm(.025)*inla_sd
#   y_inla_high <- y_inla + qnorm(.975)*inla_sd
# 
#   whatinla <- inlares$summary.fitted.values[index.obs, "mean"] -
#     cbind(1, x_tr)%*%inlares$summary.fixed$mean
#   w_inla <- y_inla - cbind(1, x_tt)%*%inlares$summary.fixed$mean
#   w_inla_low <- w_inla + qnorm(.025)*inlares$summary.random$s$sd[index.pred]
#   w_inla_high <- w_inla + qnorm(.975)*inlares$summary.random$s$sd[index.pred]
# 
#   preddata <- data.frame(coords_tt,
#                          w_tt = w_tt,
#                          whatN = wstarNNGP,
#                          whatB = wstarBRGP,
#                          whatI = w_inla,
#                          wlowerN = wquantNNGP[1,], wupperN = wquantNNGP[2,],
#                          wlowerB = wquantBRGP[1,], wupperB = wquantBRGP[2,],
#                          wlowerI = w_inla_low, wupperI = w_inla_high,
# 
#                          y_tt = y_tt,
#                          yhatN = ystarNNGP,
#                          yhatB = ystarBRGP,
#                          yhatI = y_inla,
#                          lowerN = yquantNNGP[1,], upperN = yquantNNGP[2,],
#                          lowerB = yquantBRGP[1,], upperB = yquantBRGP[2,],
#                          lowerI = y_inla_low, upperI = y_inla_high)
# 
#   predres <- preddata %>%
#     mutate(errorN_w = w_tt - whatN,
#            errorB_w = w_tt - whatB,
#            errorI_w = w_tt - whatI,
#            widthN_w = wupperN-wlowerN,
#            widthB_w = wupperB-wlowerB,
#            widthI_w = wupperI-wlowerI,
#            coverN_w = ifelse((w_tt > wlowerN & w_tt < wupperN), 1, 0),
#            coverB_w = ifelse((w_tt > wlowerB & w_tt < wupperB), 1, 0),
#            coverI_w = ifelse((w_tt > wlowerI & w_tt < wupperI), 1, 0),
# 
#            errorN = y_tt - yhatN,
#            errorB = y_tt - yhatB,
#            errorI = y_tt - yhatI,
#            widthN = upperN-lowerN,
#            widthB = upperB-lowerB,
#            widthI = upperI-lowerI,
#            coverN = ifelse((y_tt > lowerN & y_tt < upperN), 1, 0),
#            coverB = ifelse((y_tt > lowerB & y_tt < upperB), 1, 0),
#            coverI = ifelse((y_tt > lowerI & y_tt < upperI), 1, 0))
# 
#   ressum <- predres %>%
#     summarise(rmspeB_w = sqrt(mean(errorB_w^2)),
#               rmspeN_w = sqrt(mean(errorN_w^2)),
#               rmspeI_w = sqrt(mean(errorI_w^2)),
#               mapeB_w = mean(abs(errorB_w)),
#               mapeN_w = mean(abs(errorN_w)),
#               mapeI_w = mean(abs(errorI_w)),
#               meanwidthB_w = mean(widthB_w),
#               meanwidthN_w = mean(widthN_w),
#               meanwidthI_w = mean(widthI_w),
#               coverageB_w = mean(coverB_w),
#               coverageN_w = mean(coverN_w),
#               coverageI_w = mean(coverI_w),
# 
#               rmspeB = sqrt(mean(errorB^2)),
#               rmspeN = sqrt(mean(errorN^2)),
#               rmspeI = sqrt(mean(errorI^2)),
#               mapeB = mean(abs(errorB)),
#               mapeN = mean(abs(errorN)),
#               mapeI = mean(abs(errorI)),
#               meanwidthB = mean(widthB),
#               meanwidthN = mean(widthN),
#               meanwidthI = mean(widthI),
#               coverageB = mean(coverB),
#               coverageN = mean(coverN),
#               coverageI = mean(coverI))
# 
#   ## w
#   w_ressum_tab <- ressum[1:12] %>% as.numeric() %>% matrix(., nrow = 4, byrow = T)
#   colnames(w_ressum_tab) <- c("BORA-GP: latent, m=15", 
#                               "NNGP: latent, m=15", 
#                               "Barrier SGF")
#   rownames(w_ressum_tab) <- c("RMSPE", "MAPE", "Mean 95% CI width", "95% CI coverage")
#   w_ressum_tab_list[[ii]] <- w_ressum_tab %>% t()
# 
#   ## y
#   ressum_tab <- ressum[-c(1:12)] %>% as.numeric() %>% matrix(., nrow = 4, byrow = T)
#   colnames(ressum_tab) <- c("BORA-GP: latent, m=15",
#                             "NNGP: latent, m=15",
#                             "Barrier SGF")
#   rownames(ressum_tab) <- c("RMSPE", "MAPE", "Mean 95% CI width", "95% CI coverage")
#   ressum_tab_list[[ii]] <- ressum_tab %>% t()
# 
#   # plot data
#   longpreddata_w_list[[ii]] <- data.frame(easting = rep(coords_all$easting, times = 4),
#                                           northing = rep(coords_all$northing, times = 4),
#                                           what = c(w_tr, w_tt,
#                                                    whatBRGP, wstarBRGP,
#                                                    whatNNGP, wstarNNGP,
#                                                    whatinla, w_inla),
#                                           model = rep(c("Truth", 
#                                                         "BORA-GP: latent, m=15",
#                                                         "NNGP: latent, m=15",
#                                                         "Barrier SGF"), each = n)) %>%
#     mutate(wres = what - rep(c(w_tr, w_tt), times = 4))
# 
#   longpreddata_list[[ii]] <- data.frame(easting = rep(coords_all$easting, times = 4),
#                                         northing = rep(coords_all$northing, times = 4),
#                                         yhat = c(y_tr, y_tt,
#                                                  y_tr, ystarBRGP,
#                                                  y_tr, ystarNNGP,
#                                                  y_tr, y_inla),
#                                         model = rep(c("Truth", 
#                                                       "BORA-GP: latent, m=15",
#                                                       "NNGP: latent, m=15",
#                                                       "Barrier SGF"), each = n))
#   cat(paste0(ii, "th simulation completed.\n"))
# }
# 
# saveRDS(list(time_nb = time_nb, 
#              params_tab_list = params_tab_list,
#              w_ressum_tab_list = w_ressum_tab_list,
#              ressum_tab_list = ressum_tab_list),
#         paste0(path, "sim/faults_thicker_tabs.RDS"))
# saveRDS(list(longpreddata_w_list = longpreddata_w_list,
#              longpreddata_list = longpreddata_list,
#              barrier_nninfo_all = barrier_nninfo_all),
#         paste0(path, "sim/faults_thicker_preddata.RDS"))

#---------#
# results #
#---------#
tabres <- readRDS(paste0(path, "sim/faults_thicker_tabs.RDS"))
otherres <- readRDS(paste0(path, "sim/faults_thicker_preddata.RDS"))

# elapsed time for BORA-GP neighbor search
time_tab <- rbind(as.numeric(tabres$time_nb[1:3]), 
                  as.numeric(tabres$time_nb1[1:3]), 
                  as.numeric(tabres$time_nb2[1:3]))
rownames(time_tab) <- c("BORA-GP: latent, m=15",
                        "BORA-GP: latent, m=10",
                        "BORA-GP: latent, m=20")
colnames(time_tab) <- c("user", "system", "elapsed")

# elapsed time for INLA mesh
tabres$time_mesh 

# average elapsed time for model fitting (sec/iter)
((Reduce("+", tabres$params_tab_list)/rplc)[,6:8]/n.samples) %>% 
  round(3)

# parameters
(Reduce("+", tabres$params_tab_list)/rplc)[,1:5] %>% round(3)
apply(simplify2array(tabres$params_tab_list), 1:2, sd)[,1:5] %>% round(3)

# y prediction
(Reduce("+", tabres$ressum_tab_list)/rplc) %>% round(3)
apply(simplify2array(tabres$ressum_tab_list), 1:2, sd) %>% round(3)

# plot w surface
longpreddata_w_list <- otherres$longpreddata_w_list

ii <- 19
datatmp <- longpreddata_w_list[[ii]]
datatmp[datatmp$model == "BORA-GP: latent, m=15", "model"] <- "BORA-GP"
datatmp[datatmp$model == "NNGP: latent, m=15", "model"] <- "NNGP"

contourplt <- datatmp %>%
  filter(model %in% c("Truth", "BORA-GP", "Barrier SGF")) %>%
  ggplot() +
  facet_wrap(~factor(model,
                     levels = c("Truth", "BORA-GP", "Barrier SGF")),
             nrow = 1) +
  geom_contour_filled(aes(x=easting, y=northing, z=what), color = "gray88",
                      binwidth = 0.5) +
  geom_sf(data = faults, col = "black", fill = "black", linewidth = 1) +
  scale_fill_scico_d(palette = "roma", direction = -1) +
  labs(x="", y="", fill = "w") +
  theme(plot.margin=margin(t=0,l=-10, r=0, b=-15),
        legend.margin=margin(b=0,r=0,t=0,l=-3))

residplt <- datatmp %>%
  filter(model %in% c("Truth", "BORA-GP", "Barrier SGF")) %>%
  ggplot() +
  facet_wrap(~factor(model,
                     levels = c("Truth", "BORA-GP", "Barrier SGF")),
             nrow = 1) +
  geom_contour_filled(aes(x=easting, y=northing, z=-wres), color = "gray10",
                      breaks = c(-Inf, seq(-2, -0.5, by = 0.5),
                                 seq(0.5, 2, by = 0.5))) +
  geom_sf(data = faults, col = "black", fill = "black", linewidth = 1) +
  scale_fill_scico_d(palette = "vik") +
  labs(x="", y="", fill = "Residual") +                                         # expression(paste("w-", hat(w), sep=""))
  theme(plot.margin=margin(t=0,l=-10, r=0, b=-15),
        legend.margin=margin(b=0,r=0,t=0,l=-3))

gg <- ggpubr::ggarrange(contourplt, residplt, nrow = 2)
# for (ext in extension) {
#   ggsave(plot = gg, paste0(path, "plots/faults_thicker_w", ext),
#          width = 6.75, height = 5.5)
# }
