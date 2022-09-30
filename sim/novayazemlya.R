rm(list = ls())

# dependencies
library(tidyverse)
theme_set(theme_bw())
library(sf) 
library(rnaturalearth)
library(rnaturalearthdata)
library(boraGP)
library(INLA)
library(scico)

# path 
path <- "~/boraGP-sss/"

# plot extensions
extension <- c(".pdf", ".png")

#################
# data cleaning #
#################
# boundary
crs <- st_crs("+proj=utm +zone=40N +datum=WGS84 +units=km")
world <- ne_countries(scale = "medium", returnclass = "sf")
rest <- st_crop(world, xmin = 29, xmax = 88, 
                ymin = 67.8, ymax = 81) 
rest_utm <- rest %>% 
  smoothr::densify(max_distance = 1) %>%
  st_transform(crs) %>% 
  st_combine()

# call data (Jul 20, 2017 - Jul 31, 2017)
arctic <- readRDS(paste0(path, "data/NovayaZemlya_SSS_Jul2017.RDS")) %>% 
  filter(lat > 67.8 & long > 29 & long < 88) 

arctic_avg <- arctic %>% 
  group_by(long, lat) %>% 
  summarise(sss = mean(sss)) %>%
  ungroup()

## conversion to utm 
arctic_utm <- st_transform(st_as_sf(arctic_avg, 
                                    coords = c("long", "lat"),
                                    crs = 4326), 
                           crs = crs) 

## make coords 
### delete locations on land
coords <- arctic_utm %>% st_coordinates() %>% data.frame() %>% 
  rename(easting = X, northing = Y)
dlt_idx <- unlist(st_intersects(rest_utm, 
                                arctic_utm$geometry)) %>% unique()              # 5 locations
coords <- coords[-dlt_idx,]
arctic_utm <- arctic_utm[-dlt_idx,]
rm(dlt_idx)

### delete locations too close to land
dlt_idx <- unlist(st_is_within_distance(rest_utm,
                                        arctic_utm$geometry,
                                        dist = 1))                              # 4 locations
coords <- coords[-dlt_idx,]
arctic_utm <- arctic_utm[-dlt_idx,]
rm(dlt_idx)

### remove locations outside the bbox
bbox <- st_bbox(rest_utm)
dlt_idx <- which(coords$northing >= bbox$ymax)                                  # 17 locations
coords <- coords[-dlt_idx,]
arctic_utm <- arctic_utm[-dlt_idx,]
rm(dlt_idx)

# more prediction along Novaya Zemlya
bdry_utm <- readRDS(paste0(path, "data/NovayaZemlya.RDS"))
grid0 <- expand.grid(easting = seq(bbox$xmin, bbox$xmax, length = 200), 
                     northing = seq(bbox$ymin, bbox$ymax, length = 200))
grid0_sf <- st_transform(st_as_sf(grid0*1000, 
                                  coords = c("easting", "northing"),
                                  crs = 32640), 
                         crs = crs) 
nb_idx2 <- unlist(st_is_within_distance(bdry_utm$geometry,
                                        grid0_sf$geometry,
                                        dist = 100)) 
grid01 <- grid0[nb_idx2, ]
grid01_sf <- grid0_sf[nb_idx2, ]

### delete grid on land
dlt_idx <- unique(c(unlist(st_intersects(bdry_utm, grid01_sf$geometry)), 
                    unlist(st_intersects(rest_utm, grid01_sf$geometry))))
grid02 <- grid01[-dlt_idx, ]
grid02_sf <- grid01_sf[-dlt_idx, ]

### not stuck between islands 
idx_cand <- grid02_sf$geometry %>% st_coordinates() %>% as.data.frame() %>% 
  filter(Y >= 7863, Y <= 7948, X <= 400) %>% rownames() %>% 
  as.numeric()
dlt_idx <- idx_cand[c(10, 11, 22, 33, 34, 46:47, 58:59, 69, 79, 90:92, 102:104, 114)]
grid03 <- grid02[-dlt_idx, ]
grid03_sf <- grid02_sf[-dlt_idx, ]

### not too close to a boundary 
dlt_idx <- unlist(st_is_within_distance(rest_utm,
                                        grid03_sf$geometry,
                                        dist = 2))
grid <- grid03[-dlt_idx, ]
grid_sf <- grid03_sf[-dlt_idx, ]
rm(dlt_idx)

# train vs test 
n <- nrow(arctic_utm) # 8418
n_tr <- floor(0.07*n) # 589
n_tt <- n - n_tr
coords_pred <- grid
n_pred <- nrow(coords_pred)

# number of neighbors
m <- 15

# ########
# # MCMC #
# ########
rplc <- 30
# n.samples <- 15000
# burn <- 10000
set.seed(720)
seedlist <- c(6, sample.int(100, rplc-1))
seedlist[20] <- 101

# priors
## variogram
arctic_utm$ssss <- scale(arctic_utm$sss)
m_sss <- attr(arctic_utm$ssss, "scaled:center")
sd_sss <- attr(arctic_utm$ssss, "scaled:scale")
sample_vario <- gstat::variogram(ssss ~ 1, data = arctic_utm)
vario_matern <- gstat::vgm(psill = 1.2,
                           model = "Mat",
                           range = 400,
                           nugget = 0.01,
                           kappa = 1)
fit_matern <- gstat::fit.variogram(sample_vario, vario_matern)

## estimated starting value
sigma.sq <- fit_matern[2,2]
phi <- 1/fit_matern[2,3]
nu <- fit_matern[2,4]
tau.sq <- 0.01
phi.low <- phi - phi/2
phi.high <- phi + phi/2

starting <- list("phi" = phi, "sigma.sq" = sigma.sq,
                 "tau.sq" = tau.sq, "nu" = nu)
tuning <- list("phi" = 0.001, "sigma.sq" = 0.5, "tau.sq" = 0.005, "nu" = 0)     # nu is fixed
priors <- list("phi.Unif" = c(phi.low, phi.high),
               "sigma.sq.IG" = c(2, sigma.sq), "tau.sq.IG" = c(2, tau.sq),
               "nu.Unif" = c(nu-0.5,nu+0.5))
# 
# # save results
# params_tab_list = ressum_tab_list = ressum_nb_tab_list =
#   longdata_nb_list = longdata_list <- list()
# time_barrier_tab = neg_tab <- matrix(0, rplc, ncol = 3)
# 
# for (ii in 1:rplc) {
#   seed <- seedlist[ii]
#   set.seed(seed)
#   idx_tr <- sample(1:n, n_tr, replace = F)
#   nb_idx <- unlist(st_is_within_distance(bdry_utm$geometry,
#                                          arctic_utm$geometry[-idx_tr],
#                                          dist = 100))
# 
#   # training data
#   y_tr <- arctic_utm$ssss[idx_tr]
#   orgy_tr <- arctic_utm$sss[idx_tr]
#   coords_tr <- coords[idx_tr,]
#   rownames(coords_tr) <- NULL
# 
#   # ordering
#   ord <- order(coords_tr[,1])
# 
#   # test data
#   y_tt <- arctic_utm$ssss[-idx_tr]
#   orgy_tt <- arctic_utm$sss[-idx_tr]
#   coords_tt <- coords[-idx_tr,]
#   rownames(coords_tt) <- NULL
# 
#   coords_all <- rbind(coords_tr, coords_tt, coords_pred)
# 
#   #------#
#   # INLA #
#   #------#
#   # mesh
#   mesh <- inla.mesh.2d(loc = coords_tr,
#                        interior = inla.sp2segment(sf::as_Spatial(rest_utm)),
#                        max.edge = 100,
#                        offset = 200)
# 
#   # barrier model
#   tl <- length(mesh$graph$tv[,1])
#   posTri <- matrix(0, tl, 2)
#   for (t in 1:tl){
#     temp <- mesh$loc[mesh$graph$tv[t, ], ]
#     posTri[t,] <- colMeans(temp)[c(1,2)]
#   }
#   posTri <- SpatialPoints(posTri)
#   posTri_sf <- st_as_sf(posTri)
#   st_crs(posTri_sf) <- st_crs(rest_utm)
#   barrier.triangles <- unlist(st_intersects(rest_utm, posTri_sf))
#   poly.barrier <- inla.barrier.polygon(mesh, barrier.triangles)
#   # plot(poly.barrier, col = "gray")
# 
#   # model
#   # connect observations to mesh nodes
#   A.obs <- inla.spde.make.A(mesh, loc = as.matrix(coords_tr))
#   stk.obs <- inla.stack(data = list(y = y_tr),
#                         effects = list(s = 1:mesh$n, # spatial random effects
#                                        data.frame(int = rep(1,n_tr))), # intercept
#                         A = list(A.obs, 1),
#                         remove.unused = FALSE, tag = "obs")
# 
#   proj.pred <- inla.mesh.projector(mesh,
#                                    loc = as.matrix(rbind(coords_tt, coords_pred)))
#   A.pred <- inla.spde.make.A(mesh, loc = proj.pred$loc)
#   stk.pred <- inla.stack(data = list(y = NA),
#                          A = list(A.pred, 1),
#                          effects = list(s = 1:mesh$n,
#                                         data.frame(int = rep(1,n_tt+n_pred))),
#                          tag = "pred")
#   stk <- inla.stack(stk.obs, stk.pred)
# 
#   barrier.model <- inla.barrier.pcmatern(mesh,
#                                          barrier.triangles = barrier.triangles,
#                                          prior.range = c(sqrt(8)/phi.low, 0.99),# P(range < sqrt(8)/phi.low) = 0.99
#                                          prior.sigma = c(1, 0.05))              # P(sigma > 1) = 0.05
#   formula <- y ~ -1 + int + f(s, model = barrier.model)
#   timeinla <- system.time(
#     inlares <- inla(formula,
#                     data = inla.stack.data(stk),
#                     control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
#                     control.compute = list(return.marginals.predictor = TRUE),
#                     family = "gaussian",
#                     control.inla = list(int.strategy = "eb"),
#                     num.threads = 20, verbose = FALSE)
#   )
# 
#   ## Barrier SGF
#   betainla <- inlares$summary.fixed[,"mean"]
#   thetainla <- rep(0, 5)                                                        # ref: https://haakonbakkagit.github.io/btopic103.html
#   thetainla[4] <- 1
#   ### sigma = exp(theta1)
#   sigmasqmarginals <- inlares$internal.marginals.hyperpar[2] %>%
#     lapply(function(m) {
#       inla.tmarginal(function(x) exp(x)^2, m)
#     })
#   thetainla[1] <- sigmasqmarginals %>%
#     sapply(function(m)
#       unlist(inla.zmarginal(m, silent = TRUE))[1])
#   ### tausq
#   tausq_marginals <- inlares$internal.marginals.hyperpar[1] %>%
#     lapply(function(m) {
#       inla.tmarginal(function(x) 1/exp(x), m)
#     })
#   if (sum(tausq_marginals[[1]][,1] == Inf) > 0) {
#     thetainla[2] <- 0
#   } else {
#     thetainla[2] <- tausq_marginals %>%
#       sapply(function(m)
#         unlist(inla.zmarginal(m, silent = TRUE))[1])
#   }
# 
#   ### phi = sqrt(8)/r where spatial range r = exp(theta2)
#   phimarginals <- inlares$internal.marginals.hyperpar[3] %>%
#     lapply(function(m) {
#       inla.tmarginal(function(x) sqrt(8)/exp(x), m)
#     })
#   thetainla[3] <- phimarginals %>%
#     sapply(function(m)
#       unlist(inla.zmarginal(m, silent = TRUE))[1])
#   ### identifiable sigma^2*phi^(2*nu)
#   thetainla[5] <- unlist(inla.zmarginal(sigmasqmarginals[[1]]*
#                                           (phimarginals[[1]]^2),
#                                         silent = TRUE))[1]
#   #---------#
#   # BORA-GP #
#   #---------#
#   # barrier neighbor info
#   time_barrier <- system.time(
#     barrier_nninfo_all <- barrier_neighbor(coords = coords_tr,
#                                            coords.0 = rbind(coords_tt, coords_pred),
#                                            ord = ord,
#                                            n.neighbors = m,
#                                            barrier = rest_utm,
#                                            cores = 20,
#                                            verbose = FALSE,
#                                            debug = list(barrier_n.indx = NULL,
#                                                         barrier_dist = NULL,
#                                                         barrier_nn.indx.0_list = NULL,
#                                                         barrier_dist0 = NULL,
#                                                         ref_window = NULL,
#                                                         nonref_window = NULL,
#                                                         ref_fill = TRUE,
#                                                         nonref_fill = TRUE))
#   )
#   time_barrier_tab[ii,] <- as.numeric(time_barrier)[1:3]
#   barrier_nn.indx.0_list <- barrier_nninfo_all$barrier_nn.indx.0_list
#   barrier_nn.indx.0 <- do.call("rbind", barrier_nn.indx.0_list)
#   barrier_n.indx <- barrier_nninfo_all$barrier_n.indx
# 
#   barrier_nn.indx <- as.integer(unlist(barrier_n.indx)[-1]-1)
#   barrier_nn.indx.lu <- c(cumsum(c(0,0,sapply(barrier_n.indx, length)[-1]))[1:n_tr],
#                           c(0,sapply(barrier_n.indx, length)[-1])) %>% as.integer()
#   barrier_nninfo <- list(type = "barrier",
#                          n.indx = barrier_n.indx,
#                          n.neighbors = m, nn.indx = barrier_nn.indx,
#                          nn.indx.lu = barrier_nn.indx.lu, ord = ord)
# 
#   # fitting
#   barrier_m.s <- spNNGP(y_tr ~ 1, coords = coords_tr, starting = starting,
#                         method = "response", n.neighbors = m,
#                         tuning = tuning, priors = priors, cov.model = "matern",
#                         n.samples = n.samples, n.omp.threads = 20,
#                         neighbor.info = barrier_nninfo,
#                         return.neighbor.info = TRUE,
#                         verbose = FALSE)
#   barrier_p.s <- predict(barrier_m.s, X.0 = matrix(1, nrow = n_tt+n_pred),
#                          coords.0 = as.matrix(rbind(coords_tt, coords_pred)),
#                          sub.sample =
#                            list(start = burn+1, end = n.samples, thin = 1),
#                          nn.indx.0 = barrier_nn.indx.0, n.omp.threads = 20,
#                          verbose = FALSE)
#   #------#
#   # NNGP #
#   #------#
#   m.s <- spNNGP(y_tr ~ 1, coords = coords_tr, starting = starting,
#                 method = "response", n.neighbors = m,
#                 tuning = tuning, priors = priors, cov.model = "matern",
#                 n.samples = n.samples, n.omp.threads = 20, ord = ord,
#                 return.neighbor.info = TRUE, verbose = FALSE)
#   nninfo <- m.s$neighbor.info
#   p.s <- predict(m.s, X.0 = matrix(1, nrow = n_tt + n_pred),
#                  coords.0 = as.matrix(rbind(coords_tt, coords_pred)),
#                  sub.sample = list(start = burn+1, end = n.samples, thin = 1),
#                  n.omp.threads = 20, verbose = FALSE)
# 
#   #---------#
#   # Results #
#   #---------#
#   # time
#   ## boraGP
#   timeBRGP <- as.numeric(barrier_m.s$run.time)[1:3] +
#     as.numeric(barrier_p.s$run.time)[1:3]
# 
#   ## NNGP
#   timeNNGP <- as.numeric(m.s$run.time)[1:3] + as.numeric(p.s$run.time)[1:3]
# 
#   # parameter estimation
#   ## boraGP
#   betaBRGP <- mean(barrier_m.s$p.beta.samples[-c(1:burn)])
#   thetaBRGP <- c(colMeans(barrier_m.s$p.theta.samples[-c(1:burn),]),
#                  mean(barrier_m.s$p.theta.samples[-c(1:burn),1]*
#                         (barrier_m.s$p.theta.samples[-c(1:burn),3]^2)))
# 
#   ## NNGP
#   betaNNGP <- mean(m.s$p.beta.samples[-c(1:burn)])
#   thetaNNGP <- c(colMeans(m.s$p.theta.samples[-c(1:burn),]),
#                  mean(m.s$p.theta.samples[-c(1:burn),1]*
#                         (m.s$p.theta.samples[-c(1:burn),3]^2)))
# 
#   params_tab <- matrix(c(betaBRGP, thetaBRGP, timeBRGP,
#                          betaNNGP, thetaNNGP, timeNNGP,
#                          betainla, thetainla, as.numeric(timeinla)[1:3]),
#                        ncol = 3)
#   colnames(params_tab) <- c("BORA-GP",
#                             "NNGP",
#                             "Barrier SGF")
#   rownames(params_tab) <- c("beta",
#                             "sigsq",
#                             "tausq",
#                             "phi",
#                             "nu",
#                             "theta",                                            # sigsq*phi^(2*nu)
#                             "user", "system", "elapsed")                        # elapsed time to report
#   params_tab_list[[ii]] <- params_tab
# 
#   # prediction
#   ## boraGP
#   ystarBRGP <- rowMeans(barrier_p.s$p.y.0)
#   yquantBRGP <- apply(barrier_p.s$p.y.0, 1,
#                       function(x) quantile(x, probs = c(0.025, 0.975)))
#   orgystarBRGP0 <- rowMeans(barrier_p.s$p.y.0*sd_sss + m_sss)
#   neg_tab[ii,1] <- sum(orgystarBRGP0 < 0)
#   orgystarBRGP <- ifelse(orgystarBRGP0 < 0, 0, orgystarBRGP0)
#   
#   ## NNGP
#   ystarNNGP <- rowMeans(p.s$p.y.0)
#   yquantNNGP <- apply(p.s$p.y.0, 1,
#                       function(x) quantile(x, probs = c(0.025, 0.975)))
#   orgystarNNGP0 <- rowMeans(p.s$p.y.0*sd_sss + m_sss)
#   neg_tab[ii,2] <- sum(orgystarNNGP0 < 0)
#   orgystarNNGP <- ifelse(orgystarNNGP0 < 0, 0, orgystarNNGP0)
#   
#   ## inla
#   index.obs <- inla.stack.index(stk, "obs")$data
#   index.pred <- inla.stack.index(stk, "pred")$data
#   y_inla <- inlares$summary.fitted.values[index.pred, "mean"]
#   inla_sd <- sqrt(inlares$summary.fitted.values[index.pred,"sd"]^2              # variability from the coefficient
#                   + 1/inlares$summary.hyperpar[1,"mean"])                       # tausq
#   y_inla_low <- y_inla + qnorm(.025)*inla_sd
#   y_inla_high <- y_inla + qnorm(.975)*inla_sd
#   orgy_inla0 <- sapply(inlares$marginals.fitted.values[index.pred],
#                       function(s) inla.emarginal(function(x) x*sd_sss + m_sss, s))
#   neg_tab[ii,3] <- sum(orgy_inla0 < 0)
#   orgy_inla <- ifelse(orgy_inla0 < 0, 0, orgy_inla0)
# 
#   # prediction summary
#   preddata <- data.frame(coords_tt,
#                          y_tt = y_tt,
#                          yhatN = ystarNNGP[1:n_tt],
#                          yhatB = ystarBRGP[1:n_tt],
#                          yhatI = y_inla[1:n_tt],
#                          lowerN = yquantNNGP[1,1:n_tt],
#                          upperN = yquantNNGP[2,1:n_tt],
#                          lowerB = yquantBRGP[1,1:n_tt],
#                          upperB = yquantBRGP[2,1:n_tt],
#                          lowerI = y_inla_low[1:n_tt],
#                          upperI = y_inla_high[1:n_tt])
# 
#   predres <- preddata %>%
#     mutate(errorN = y_tt - yhatN,
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
#     summarise(rmspeB = sqrt(mean(errorB^2)),
#               rmspeN = sqrt(mean(errorN^2)),
#               rmspeI = sqrt(mean(errorI^2)),
#               mapeB = mean(abs(errorB)),
#               mapeN = mean(abs(errorN)),
#               mapeI = mean(abs(errorI)),
#               coverageB = mean(coverB),
#               coverageN = mean(coverN),
#               coverageI = mean(coverI),
#               meanwidthB = mean(widthB),
#               meanwidthN = mean(widthN),
#               meanwidthI = mean(widthI))
# 
#   ressum_nb <- predres[nb_idx,] %>%
#     summarise(rmspeB = sqrt(mean(errorB^2)),
#               rmspeN = sqrt(mean(errorN^2)),
#               rmspeI = sqrt(mean(errorI^2)),
#               mapeB = mean(abs(errorB)),
#               mapeN = mean(abs(errorN)),
#               mapeI = mean(abs(errorI)),
#               coverageB = mean(coverB),
#               coverageN = mean(coverN),
#               coverageI = mean(coverI),
#               meanwidthB = mean(widthB),
#               meanwidthN = mean(widthN),
#               meanwidthI = mean(widthI))
# 
#   ## whole
#   ressum_tab <- ressum %>% as.numeric() %>% matrix(., nrow = 4, byrow = T)
#   colnames(ressum_tab) <- c("BORA-GP",
#                             "NNGP",
#                             "Barrier SGF")
#   rownames(ressum_tab) <- c("RMSPE", "MAPE",
#                             "95% CI coverage",
#                             "Mean 95% CI width")
#   ressum_tab_list[[ii]] <- ressum_tab
# 
#   ## near Novaya Zemlya
#   ressum_nb_tab <- ressum_nb %>% as.numeric() %>% matrix(., nrow = 4, byrow = T)
#   colnames(ressum_nb_tab) <- c("BORA-GP",
#                                "NNGP",
#                                "Barrier SGF")
#   rownames(ressum_nb_tab) <- c("RMSPE", "MAPE",
#                                "95% CI coverage",
#                                "Mean 95% CI width")
#   ressum_nb_tab_list[[ii]] <- ressum_nb_tab
# 
#   # plot data
#   longdata_nb <- tibble(easting = coords_all[c(n_tr + nb_idx, n+(1:n_pred)),1],
#                         northing = coords_all[c(n_tr + nb_idx, n+(1:n_pred)),2],
#                         "SMAP measurements" = c(orgy_tt[nb_idx], rep(NA, n_pred)),
#                         "BORA-GP" = orgystarBRGP[c(nb_idx, n_tt+(1:n_pred))],
#                         NNGP = orgystarNNGP[c(nb_idx, n_tt+(1:n_pred))],
#                         "Barrier SGF" = orgy_inla[c(nb_idx, n_tt+(1:n_pred))]) %>%
#     pivot_longer(-c("easting", "northing"),
#                  names_to = "model", values_to = "vals")
#   longdata_nb_list[[ii]] <- longdata_nb
# 
#   longdata <- data.frame(easting = rep(coords_all$easting, times = 4),
#                          northing = rep(coords_all$northing, times = 4),
#                          yhat = c(y_tr, y_tt, rep(NA, n_pred),
#                                   y_tr, ystarBRGP,
#                                   y_tr, ystarNNGP,
#                                   y_tr, y_inla),
#                          orgyhat = c(orgy_tr, orgy_tt, rep(NA, n_pred),
#                                      orgy_tr, orgystarBRGP,
#                                      orgy_tr, orgystarNNGP,
#                                      orgy_tr, orgy_inla),
#                          model = rep(c("Truth",
#                                        "BORA-GP",
#                                        "NNGP",
#                                        "Barrier SGF"), each = n + n_pred))
#   longdata_list[[ii]] <- longdata
# 
#   saveRDS(list(predres = predres,
#                neg = neg_tab[ii,],
#                nninfo = nninfo,
#                barrier_nninfo_all = barrier_nninfo_all,
#                time_barrier = time_barrier,
#                params_tab = params_tab,
#                ressum_tab = ressum_tab,
#                ressum_nb_tab = ressum_nb_tab,
#                longdata_nb = longdata_nb,
#                longdata = longdata),
#           paste0(path, "sim/novayazemlya/", ii, ".RDS"))
# 
#   cat(paste0(ii, "th simulation completed.\n"))
# }
# 
# saveRDS(list(params_tab_list = params_tab_list,
#              ressum_tab_list = ressum_tab_list,
#              ressum_nb_tab_list = ressum_nb_tab_list,
#              longdata_list = longdata_list,
#              longdata_nb_list = longdata_nb_list,
#              time_barrier_tab = time_barrier_tab, 
#              neg_tab = neg_tab),
#         paste0(path, "sim/novayazemlya.RDS"))

#------------------------------------------------------------------------------#
res <- readRDS(paste0(path, "sim/novayazemlya.RDS"))

# parameter estimation
(Reduce("+", res$params_tab_list)/rplc) %>% round(3)
apply(simplify2array(res$params_tab_list), 1:2, sd) %>% round(3)

# prediction
(Reduce("+", res$ressum_tab_list)/rplc) %>% round(3)
apply(simplify2array(res$ressum_tab_list), 1:2, sd) %>% round(3)

# plot
gg_4 <- res$longdata_list[[1]][c(1:n_tr, (n + n_pred) + 1:n,
                                 2*(n + n_pred) + 1:n, 3*(n + n_pred) + 1:n),] %>%
  ggplot() +
  geom_point(aes(easting, northing, color = orgyhat), size = 1) +
  geom_sf(data = rest_utm) +
  facet_wrap(~ factor(model, levels = c("Truth", "BORA-GP",
                                        "NNGP", "Barrier SGF"),
                      labels = c("Observations", "BORA-GP",
                                 "NNGP", "Barrier SGF")),
             nrow = 1) +
  scale_color_distiller(palette = "RdYlBu", na.value = "transparent") +
  labs(col = "SSS\n(psu)", x = "", y = "") +
  theme(plot.margin = margin(t = -5, l = -10, r = 2, b = 0),
        legend.position = "right")

## near Novaya Zemlya
gg_novaya <- res$longdata_nb_list[[1]] %>%
  filter(model != "SMAP measurements") %>%
  ggplot() +
  geom_point(aes(easting, northing, color = vals)) +
  geom_sf(data = bdry_utm) +
  facet_wrap(~ factor(model, levels = c("BORA-GP", "NNGP", "Barrier SGF")),
             nrow = 1) +
  scale_color_distiller(palette = "RdYlBu", na.value = "transparent") +
  labs(col = "SSS\n(psu)", x = "", y = "") +
  theme(plot.margin = margin(t = -15, l = -30, r = -50, b = -10),
        legend.margin = margin(b = 0, r = -10, t = 0, l = -1))

## prediction comparison near Novaya Zemlya
nngp_rmspe <- sapply(res$ressum_nb_tab_list, function(x) x[1,"NNGP"])
boragp_rmspe <- sapply(res$ressum_nb_tab_list, function(x) x[1,"BORA-GP"])
barrier_rmspe <- sapply(res$ressum_nb_tab_list, function(x) x[1,"Barrier SGF"])

nngp_mape <- sapply(res$ressum_nb_tab_list, function(x) x[2,"NNGP"])
boragp_mape <- sapply(res$ressum_nb_tab_list, function(x) x[2,"BORA-GP"])
barrier_mape <- sapply(res$ressum_nb_tab_list, function(x) x[2,"Barrier SGF"])

## from BORA-GP to NNGP
((nngp_rmspe - boragp_rmspe)/boragp_rmspe*100) %>% summary() %>% round(3)
((nngp_mape - boragp_mape)/boragp_mape*100) %>% summary() %>% round(3)
## from BORA-GP to Barrier SGF
((barrier_rmspe - boragp_rmspe)/boragp_rmspe*100) %>% summary() %>% round(3)
((barrier_mape - boragp_mape)/boragp_mape*100) %>% summary() %>% round(3)

gg_mape <- data.frame(alternative = rep(rep(c("NNGP", "Barrier SGF"), each = rplc),
                             times = 2),
           vals = c((nngp_rmspe - boragp_rmspe)/boragp_rmspe,
                    (barrier_rmspe - boragp_rmspe)/boragp_rmspe,
                    (nngp_mape - boragp_mape)/boragp_mape,
                    (barrier_mape - boragp_mape)/boragp_mape),
           idx = rep(1:rplc, times = 4),
           what = factor(rep(c("RMSPE", "MAPE"), each = 2*rplc),
                         levels = c("RMSPE", "MAPE"))) %>%
  filter(what == "MAPE") %>%
  ggplot() +
  geom_hline(yintercept = 0, color = "gray") +
  geom_boxplot(aes(factor(alternative, levels = c("NNGP", "Barrier SGF")),
                   vals, fill = alternative)) +
  labs(x = "", y = "Increases in MAPE from BORA-GP") +
  scale_fill_manual(values = c("NNGP" = "#97872E", "Barrier SGF" = "#184E60")) +
  scale_y_continuous(labels = scales::label_percent()) +
  theme(plot.margin = margin(t = -15, l = 2, r = -10, b = -10),
        legend.position = "none")

gg_nb <- ggpubr::ggarrange(gg_mape, gg_novaya, nrow = 1, widths = c(1, 2.5))
gg_4_zoom <- ggpubr::ggarrange(gg_4, gg_nb, nrow = 2)
# for (ext in extension) {
#   ggsave(plot = gg_4_zoom,
#          paste0(path, "plots/novayazemlya_pred", ext),
#          width = 11, height = 5.5)
# }

# different neighbors
res1 <- readRDS(paste0(path, "/sim/novayazemlya_seed6.RDS"))

## generate data for the seed
seed <- seedlist[1]
set.seed(seed)
idx_tr <- sample(1:n, n_tr, replace = F)
nb_idx <- unlist(st_is_within_distance(bdry_utm$geometry,
                                       arctic_utm$geometry[-idx_tr],
                                       dist = 100))

# training data
y_tr <- arctic_utm$sss[idx_tr]
coords_tr <- coords[idx_tr,]
rownames(coords_tr) <- NULL

# ordering
ord <- order(coords_tr[,1])

# test data
y_tt <- arctic_utm$sss[-idx_tr]
coords_tt <- coords[-idx_tr,]
rownames(coords_tt) <- NULL

coords_all <- rbind(coords_tr, coords_tt, coords_pred)

## BORA-GP neighbor
barrier_nninfo_all <- res1$barrier_nninfo_all
barrier_nn.indx.0_list <- barrier_nninfo_all$barrier_nn.indx.0_list
barrier_nn.indx.0 <- do.call("rbind", barrier_nn.indx.0_list)
barrier_n.indx <- barrier_nninfo_all$barrier_n.indx

barrier_nn.indx <- as.integer(unlist(barrier_n.indx)[-1]-1)
barrier_nn.indx.lu <- c(cumsum(c(0,0,sapply(barrier_n.indx, length)[-1]))[1:n_tr],
                        c(0,sapply(barrier_n.indx, length)[-1])) %>% as.integer()
barrier_nninfo <- list(type = "barrier",
                       n.indx = barrier_n.indx,
                       n.neighbors = m, nn.indx = barrier_nn.indx,
                       nn.indx.lu = barrier_nn.indx.lu, ord = ord)

## correlation
Ctilde <- boraGP:::create_Ctilde(coords = coords_tr,
                                 neighbor.info = res1$nninfo,
                                 sig_sq = 1,
                                 phi = res1$params_tab[4,2],
                                 nu = res1$params_tab[5,2],
                                 base_cov = "matern", cores = 20)

barrier_Ctilde <- boraGP:::create_Ctilde(coords = coords_tr,
                                         neighbor.info = barrier_nninfo,
                                         sig_sq = 1,
                                         phi = res1$params_tab[4,1],
                                         nu = res1$params_tab[5,1],
                                         base_cov = "matern", cores = 20)
## preprocessing
nnidx_tt <- RANN::nn2(coords_tr, coords_tt, k = m)$nn.idx
barrier_nn.indx.0_ord <- list()
for (i in 1:n_tt) {
  barrier_nn.indx.0_ord[[i]] <- order(ord)[barrier_nn.indx.0_list[[i]]]
}
coords_tr_ord <- coords_tr[ord,]
nn.indx.0_ord <- RANN::nn2(coords_tr_ord, coords_tt, k = m)$nn.idx

pt <- c(4803, 4175)
gg_cov_list <- list()
for (i in 1:length(pt)) {
  ptnn <- nnidx_tt[pt[i],]
  barrier_ptnn <- barrier_nn.indx.0_list[[pt[i]]]

  NNGP_cov <- boraGP:::NGPcov_m(v1 = coords_all[n_tr + pt[i],],
                                v2_mat = coords_all[1:n,],
                                coords = coords_tr,
                                neighbor.info = res1$nninfo,
                                sig_sq = 1,
                                phi = res1$params_tab[4,2],
                                nu = res1$params_tab[5,2],
                                base_cov = "matern",
                                Ctilde = Ctilde,
                                coords.0 = coords_tt,
                                nn.indx.0_ord = nn.indx.0_ord)
  
  boraGP_cov <- boraGP:::NGPcov_m(v1 = coords_all[n_tr + pt[i],],
                                  v2_mat = coords_all[1:n,],
                                  coords = coords_tr, neighbor.info = barrier_nninfo,
                                  sig_sq = 1,
                                  phi = res1$params_tab[4,1],
                                  nu = res1$params_tab[5,1],
                                  base_cov = "matern",
                                  Ctilde = barrier_Ctilde,
                                  coords.0 = coords_tt,
                                  nn.indx.0_ord = barrier_nn.indx.0_ord)
  
  gg_cov_list[[i]] <- tibble(coords_all[1:n,],  
                             SMAP = rep(NA, n), 
                             NNGP = NNGP_cov, 
                             "BORA-GP" = boraGP_cov) %>%
    pivot_longer(-c(easting, northing), values_to = "cor", names_to = "what") %>% 
    filter(cor > 0.2) %>% 
    ggplot() +
    facet_grid(~ factor(what, levels = c("SMAP", "BORA-GP", "NNGP"), 
                        labels = c("15 neighbors", "BORA-GP", "NNGP"))) + 
    geom_point(aes(easting, northing, col = cor)) +
    geom_sf(data = rest_utm) + 
    geom_point(data = data.frame(coords_all[ptnn, ], what = "SMAP"),
               aes(easting, northing), 
               col = "#888888", size = 3, shape = 19) +
    geom_point(data = data.frame(coords_all[barrier_ptnn,], what = "SMAP"),
               aes(easting, northing),
               col = "#DC3220", size = 2, shape = 17) +
    geom_point(data = coords_all[n_tr+pt[i],], aes(easting, northing),
               col = "black", size = 2) +
    labs(x = "", y = "", color = "Cor") +
    scale_color_scico(palette = "batlow", limits = c(0.2, 1),
                      na.value = "transparent") +
    theme(plot.margin = margin(t = 0, l = -5, r = 0, b = -5),
          legend.margin = margin(t = 0, l = 0, r = 0, b = 0))
}
gg_covall <- ggpubr::ggarrange(gg_cov_list[[1]],
                               gg_cov_list[[2]],
                               common.legend = T, nrow = 2, 
                               legend = "right")
# for (ext in extension) {
#   ggsave(plot = gg_covall,
#          paste0(path, "plots/novayazemlya_cov", ext),
#          width = 10, height = 5)
# }

# # UQ
# res1$predres %>% 
#   select(easting, northing, coverB, coverN, coverI) %>% 
#   pivot_longer(-c(easting, northing), 
#                values_to = "cover", 
#                names_to = "model") %>% 
#   mutate(model = factor(model, 
#                         levels = c("coverB", "coverN", "coverI"),
#                         labels = c("BORA-GP", "NNGP", "Barrier SGF")), 
#          cover = factor(cover, 
#                         levels = c("1", "0"),
#                         labels = c("Yes", "No"))) %>% 
#   ggplot() +
#   geom_point(aes(easting, northing, color = cover)) +
#   geom_sf(data = rest_utm) +
#   facet_wrap(~ model) +
#   scale_color_manual(values = c("No" = "#ffae42", "Yes" = "gray")) +
#   labs(x = "", y = "", color = "In 95% CI?") +
#   theme(plot.margin=margin(t=0,l=0, r=0, b=0),
#         legend.margin=margin(b=-5,r=0,t=-5,l=0),
#         legend.position = "bottom",
#         legend.text = element_text(size = 12))
