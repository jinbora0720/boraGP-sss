rm(list = ls())

# dependencies
library(tidyverse)
theme_set(theme_bw())
library(sf) 
sf_use_s2(FALSE)
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

########
# MCMC #
########
rplc <- 30
n.samples <- 15000
burn <- 10000
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

# # save results
# time_barrier_tab <- matrix(0, rplc, ncol = 3)
# neg_tab <- matrix(0, rplc, ncol = 2)
# 
# for (ii in 1) {
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
#   # test data
#   y_tt <- arctic_utm$ssss[-idx_tr]
#   orgy_tt <- arctic_utm$sss[-idx_tr]
#   coords_tt <- coords[-idx_tr,]
#   rownames(coords_tt) <- NULL
# 
#   coords_all <- rbind(coords_tr, coords_tt, coords_pred)
#
#   ###################################
#   # ordering 1: ascending in y-axis #
#   ###################################
#   ord <- order(coords_tr[,2])
# 
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
#                          betaNNGP, thetaNNGP, timeNNGP),
#                        ncol = 2)
#   colnames(params_tab) <- c("BORA-GP",
#                             "NNGP")
#   rownames(params_tab) <- c("beta",
#                             "sigsq",
#                             "tausq",
#                             "phi",
#                             "nu",
#                             "theta",                                            # sigsq*phi^(2*nu)
#                             "user", "system", "elapsed")                        # elapsed time to report
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
#   # prediction summary
#   preddata <- data.frame(coords_tt,
#                          y_tt = y_tt,
#                          yhatN = ystarNNGP[1:n_tt],
#                          yhatB = ystarBRGP[1:n_tt],
#                          lowerN = yquantNNGP[1,1:n_tt],
#                          upperN = yquantNNGP[2,1:n_tt],
#                          lowerB = yquantBRGP[1,1:n_tt],
#                          upperB = yquantBRGP[2,1:n_tt])
# 
#   predres <- preddata %>%
#     mutate(errorN = y_tt - yhatN,
#            errorB = y_tt - yhatB,
#            widthN = upperN-lowerN,
#            widthB = upperB-lowerB,
#            coverN = ifelse((y_tt > lowerN & y_tt < upperN), 1, 0),
#            coverB = ifelse((y_tt > lowerB & y_tt < upperB), 1, 0))
# 
#   ressum <- predres %>%
#     summarise(rmspeB = sqrt(mean(errorB^2)),
#               rmspeN = sqrt(mean(errorN^2)),
#               mapeB = mean(abs(errorB)),
#               mapeN = mean(abs(errorN)),
#               coverageB = mean(coverB),
#               coverageN = mean(coverN),
#               meanwidthB = mean(widthB),
#               meanwidthN = mean(widthN))
# 
#   ressum_nb <- predres[nb_idx,] %>%
#     summarise(rmspeB = sqrt(mean(errorB^2)),
#               rmspeN = sqrt(mean(errorN^2)),
#               mapeB = mean(abs(errorB)),
#               mapeN = mean(abs(errorN)),
#               coverageB = mean(coverB),
#               coverageN = mean(coverN),
#               meanwidthB = mean(widthB),
#               meanwidthN = mean(widthN))
# 
#   ## whole
#   ressum_tab <- ressum %>% as.numeric() %>% matrix(., nrow = 4, byrow = T)
#   colnames(ressum_tab) <- c("BORA-GP",
#                             "NNGP")
#   rownames(ressum_tab) <- c("RMSPE", "MAPE",
#                             "95% CI coverage",
#                             "Mean 95% CI width")
# 
#   ## near Novaya Zemlya
#   ressum_nb_tab <- ressum_nb %>% as.numeric() %>% matrix(., nrow = 4, byrow = T)
#   colnames(ressum_nb_tab) <- c("BORA-GP",
#                                "NNGP")
#   rownames(ressum_nb_tab) <- c("RMSPE", "MAPE",
#                                "95% CI coverage",
#                                "Mean 95% CI width")
# 
#   # plot data
#   longdata_nb <- tibble(easting = coords_all[c(n_tr + nb_idx, n+(1:n_pred)),1],
#                         northing = coords_all[c(n_tr + nb_idx, n+(1:n_pred)),2],
#                         "SMAP measurements" = c(orgy_tt[nb_idx], rep(NA, n_pred)),
#                         "BORA-GP" = orgystarBRGP[c(nb_idx, n_tt+(1:n_pred))],
#                         NNGP = orgystarNNGP[c(nb_idx, n_tt+(1:n_pred))]) %>%
#     pivot_longer(-c("easting", "northing"),
#                  names_to = "model", values_to = "vals")
# 
#   longdata <- data.frame(easting = rep(coords_all$easting, times = 3),
#                          northing = rep(coords_all$northing, times = 3),
#                          yhat = c(y_tr, y_tt, rep(NA, n_pred),
#                                   y_tr, ystarBRGP,
#                                   y_tr, ystarNNGP),
#                          orgyhat = c(orgy_tr, orgy_tt, rep(NA, n_pred),
#                                      orgy_tr, orgystarBRGP,
#                                      orgy_tr, orgystarNNGP),
#                          model = rep(c("Truth",
#                                        "BORA-GP",
#                                        "NNGP"), each = n + n_pred))
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
#           paste0(path, "sim/novayazemlya_y_1st.RDS"))
# 
#   #####################################
#   # ordering 2: ascending in x+y-axis #
#   #####################################
#   ord <- order(coords_tr[,1] + coords_tr[,2]) 
#   # bdry_utm %>% 
#   #   ggplot() +
#   #   geom_sf() +
#   #   geom_point(data = coords_tr[ord[1:m],], 
#   #              aes(easting, northing))
#   
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
#                          betaNNGP, thetaNNGP, timeNNGP),
#                        ncol = 2)
#   colnames(params_tab) <- c("BORA-GP",
#                             "NNGP")
#   rownames(params_tab) <- c("beta",
#                             "sigsq",
#                             "tausq",
#                             "phi",
#                             "nu",
#                             "theta",                                            # sigsq*phi^(2*nu)
#                             "user", "system", "elapsed")                        # elapsed time to report
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
#   # prediction summary
#   preddata <- data.frame(coords_tt,
#                          y_tt = y_tt,
#                          yhatN = ystarNNGP[1:n_tt],
#                          yhatB = ystarBRGP[1:n_tt],
#                          lowerN = yquantNNGP[1,1:n_tt],
#                          upperN = yquantNNGP[2,1:n_tt],
#                          lowerB = yquantBRGP[1,1:n_tt],
#                          upperB = yquantBRGP[2,1:n_tt])
#   
#   predres <- preddata %>%
#     mutate(errorN = y_tt - yhatN,
#            errorB = y_tt - yhatB,
#            widthN = upperN-lowerN,
#            widthB = upperB-lowerB,
#            coverN = ifelse((y_tt > lowerN & y_tt < upperN), 1, 0),
#            coverB = ifelse((y_tt > lowerB & y_tt < upperB), 1, 0))
#   
#   ressum <- predres %>%
#     summarise(rmspeB = sqrt(mean(errorB^2)),
#               rmspeN = sqrt(mean(errorN^2)),
#               mapeB = mean(abs(errorB)),
#               mapeN = mean(abs(errorN)),
#               coverageB = mean(coverB),
#               coverageN = mean(coverN),
#               meanwidthB = mean(widthB),
#               meanwidthN = mean(widthN))
#   
#   ressum_nb <- predres[nb_idx,] %>%
#     summarise(rmspeB = sqrt(mean(errorB^2)),
#               rmspeN = sqrt(mean(errorN^2)),
#               mapeB = mean(abs(errorB)),
#               mapeN = mean(abs(errorN)),
#               coverageB = mean(coverB),
#               coverageN = mean(coverN),
#               meanwidthB = mean(widthB),
#               meanwidthN = mean(widthN))
#   
#   ## whole
#   ressum_tab <- ressum %>% as.numeric() %>% matrix(., nrow = 4, byrow = T)
#   colnames(ressum_tab) <- c("BORA-GP",
#                             "NNGP")
#   rownames(ressum_tab) <- c("RMSPE", "MAPE",
#                             "95% CI coverage",
#                             "Mean 95% CI width")
# 
#   ## near Novaya Zemlya
#   ressum_nb_tab <- ressum_nb %>% as.numeric() %>% matrix(., nrow = 4, byrow = T)
#   colnames(ressum_nb_tab) <- c("BORA-GP",
#                                "NNGP")
#   rownames(ressum_nb_tab) <- c("RMSPE", "MAPE",
#                                "95% CI coverage",
#                                "Mean 95% CI width")
# 
#   # plot data
#   longdata_nb <- tibble(easting = coords_all[c(n_tr + nb_idx, n+(1:n_pred)),1],
#                         northing = coords_all[c(n_tr + nb_idx, n+(1:n_pred)),2],
#                         "SMAP measurements" = c(orgy_tt[nb_idx], rep(NA, n_pred)),
#                         "BORA-GP" = orgystarBRGP[c(nb_idx, n_tt+(1:n_pred))],
#                         NNGP = orgystarNNGP[c(nb_idx, n_tt+(1:n_pred))]) %>%
#     pivot_longer(-c("easting", "northing"),
#                  names_to = "model", values_to = "vals")
# 
#   longdata <- data.frame(easting = rep(coords_all$easting, times = 3),
#                          northing = rep(coords_all$northing, times = 3),
#                          yhat = c(y_tr, y_tt, rep(NA, n_pred),
#                                   y_tr, ystarBRGP,
#                                   y_tr, ystarNNGP),
#                          orgyhat = c(orgy_tr, orgy_tt, rep(NA, n_pred),
#                                      orgy_tr, orgystarBRGP,
#                                      orgy_tr, orgystarNNGP),
#                          model = rep(c("Truth",
#                                        "BORA-GP",
#                                        "NNGP"), each = n + n_pred))
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
#           paste0(path, "sim/novayazemlya_xpy_1st.RDS"))
#   
#   cat(paste0(ii, "th simulation completed.\n"))
# }

#------------------------------------------------------------------------------#
res_x <- readRDS(paste0(path, "sim/novayazemlya_1st.RDS"))
res_y <- readRDS(paste0(path, "sim/novayazemlya_y_1st.RDS"))
res_xpy <- readRDS(paste0(path, "sim/novayazemlya_xpy_1st.RDS"))

# parameter estimation
res_x$params_tab %>% round(3)
res_y$params_tab %>% round(3)
res_xpy$params_tab %>% round(3)

# prediction
res_x$ressum_tab %>% round(3)
res_y$ressum_tab %>% round(3)
res_xpy$ressum_tab %>% round(3)

# prediction near boundary
res_x$ressum_nb_tab %>% round(3)
res_y$ressum_nb_tab %>% round(3)
res_xpy$ressum_nb_tab %>% round(3)

# plot
n_tr
n_tt
n_pred
n
res_x$longdata_nb %>% nrow()
res_y$longdata_nb %>% nrow()
res_xpy$longdata_nb %>% nrow()

rbind(res_x$longdata_nb %>% 
        filter(model == "BORA-GP") %>% 
        mutate(order = "x-coordinates") %>%
        slice_tail(n = n_pred), 
      res_y$longdata_nb %>% 
        filter(model == "BORA-GP") %>% 
        mutate(order = "y-coordinates") %>%
        slice_tail(n = n_pred), 
      res_xpy$longdata_nb %>% 
        filter(model == "BORA-GP") %>% 
        mutate(order = "x+y-coordinates") %>%
        slice_tail(n = n_pred)) %>%
  ggplot() +
  geom_point(aes(easting, northing, color = vals)) +
  geom_contour(aes(x = easting, y = northing, z = vals), 
               color = "white", 
               breaks = c(1, 35)) + 
  geom_sf(data = bdry_utm) +
  facet_wrap(~ factor(order, 
                      levels = c("x-coordinates", "y-coordinates", "x+y-coordinates")),
             nrow = 1) +
  scale_color_distiller(palette = "RdYlBu", na.value = "transparent", 
                        limits = c(0, max(arctic_utm$sss))) +
  labs(col = "SSS\n(psu)", x = "", y = "") +
  theme(plot.margin = margin(t = 0, l = -10, r = 0, b = -10),
        legend.margin = margin(b = 0, r = 0, t = 0, l = -1))

for (ext in extension) {
  ggsave(paste0(path, "plots/novayazemlya_flip", ext),
         width = 6.5, height = 3)
}