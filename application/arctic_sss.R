# barrier = land 
rm(list = ls())

# dependencies
library(tidyverse)
theme_set(theme_bw())
library(sf) 
library(rnaturalearth)
library(rnaturalearthdata)
library(boraGP)
# install.packages("INLA",
#                  repos=c(getOption("repos"),
#                          INLA="https://inla.r-inla-download.org/R/testing"),
#                  dep=TRUE)
library(INLA)
library(scico)

# path 
path <- "~/boraGP-sss/"

# plot extensions
extension <- c(".pdf", ".png")

# Arctic land
# sf::sf_use_s2(FALSE) 
world <- ne_countries(scale = "medium", returnclass = "sf")
arctic_utm <- st_crop(world, xmin = -180, xmax = 180, 
                      ymin = 57, ymax = 90) %>% 
  smoothr::densify(max_distance = 1) %>%
  st_transform("+proj=utm +zone=40N +datum=WGS84 +units=km") %>% 
  st_combine()

# call data (August, 2020)
# sss data source: https://thredds.jpl.nasa.gov/thredds/dodsC/ncml_aggregation/SalinityDensity/smap/aggregate__SMAP_JPL_L3_SSS_CAP_MONTHLY_V5.ncml
# sea ice extent data source: https://nsidc.org/data/G02135/versions/3 
sss <- readRDS(paste0(path, "data/SSS_2020Aug.RDS")) # 52067
ice_utm <- readRDS(paste0(path, "data/Iceextent_2020Aug.RDS")) %>% 
  st_combine()

# UTM projection
crs <- st_crs("+proj=utm +zone=40N +datum=WGS84 +units=km")
sss_utm <- st_transform(st_as_sf(sss, 
                                 coords = c("long", "lat"),
                                 crs = 4326), crs = crs) 

# combine arctic and ice to make one barrier
a <- st_sf(geom = arctic_utm, crs = crs)
b <- st_sf(geom = ice_utm, crs = crs)
barrier_utm <- st_combine(rbind(a, b))
rm(a, b)

# delete ones on or too close to the land in the data
dlt_barrier <- unlist(st_is_within_distance(arctic_utm,
                                            sss_utm$geometry, 
                                            dist = 1)) %>% unique()
sss_utm <- sss_utm[-dlt_barrier,] # 52005
rm(dlt_barrier)

# SSS transformation
sss_utm$ssss <- sqrt(sss_utm$sss)
hist(sss_utm$ssss)

# filling the mismatch 
long <- seq(-179.875, 179.875, by = 0.25)
lat <- sort(unique(sss$lat))
coords_longlat <- expand.grid(long, lat) %>% 
  as.data.frame() %>% 
  rename(long = Var1, lat = Var2)
coords_all_utm <- st_transform(st_as_sf(coords_longlat, 
                                        coords = c("long", "lat"),
                                        crs = 4326), crs = crs) 
coords <- coords_all_utm %>% 
  st_coordinates() %>% 
  as.data.frame() %>% 
  rename(easting = X, northing = Y)
rm(coords_longlat)

# delete ones on or too close to the land or the ice in grid
dlt_overlap <- unlist(st_intersects(sss_utm$geometry, 
                                    coords_all_utm$geometry)) %>% unique()
dlt_barrier2 <- unlist(st_is_within_distance(barrier_utm,
                                             coords_all_utm$geometry, 
                                             dist = 9)) %>% unique()
dlt_blocked <- coords %>% 
  mutate(idx = 1:nrow(coords)) %>% 
  filter((easting > -1950 & easting < -800 & northing < 7800) | 
           (northing > 12800 & easting > 800 & easting < 2500) | 
           (easting > 3500 & northing < 11000) | 
           (northing < 8050 & easting > 1150 & easting < 1600)) %>% 
  pull(idx)
dlt_all <- unique(c(dlt_barrier2, dlt_overlap, dlt_blocked))

coords_pred <- coords[-dlt_all,] # 3702
n_pred <- nrow(coords_pred)

# delete blocked areas 
dlt_blocked <- sss_utm$geometry %>% 
  st_coordinates() %>% 
  as.data.frame() %>% 
  rename(easting = X, northing = Y) %>% 
  mutate(idx = 1:nrow(sss_utm)) %>% 
  filter((easting > -1950 & easting < -800 & northing < 7800) | 
           (northing > 12800 & easting > 800 & easting < 2500) | 
           (easting > 3500 & northing < 11000)| 
           (northing < 8050 & easting > 1150 & easting < 1600)) %>% 
  pull(idx)

sss_utm2 <- sss_utm[-dlt_blocked,]
coords_tr <- sss_utm2$geometry %>% 
  st_coordinates() %>% 
  as.data.frame() %>% 
  rename(easting = X, northing = Y)
n_tr <- nrow(coords_tr) # 50225
orgy_tr <- sss_utm2$sss
y_tr <- sss_utm2$ssss
coords_all <- rbind(coords_tr, coords_pred)

rm(dlt_all, dlt_barrier2, coords, coords_all_utm)

# SSS in the Arctic Ocean and in Canadian archipelago 
gg_arctic <- sss_utm2 %>%
  ggplot() +
  geom_sf(aes(color = sss)) +
  geom_sf(data = ice_utm, fill = "skyblue") +
  geom_sf_label(data = ice_utm, label = "sea-ice") + 
  geom_sf(data = arctic_utm) + 
  scale_color_distiller(palette = "RdYlBu", 
                        na.value = "transparent") + 
  labs(x = "", y = "", color = "SSS\n(psu)") +
  theme(plot.margin = margin(t = 0, l = 0, r = 0, b = 0))

gg_canada <- sss_utm2 %>%
  ggplot() +
  geom_sf(aes(color = sss)) +
  geom_sf(data = ice_utm, fill = "skyblue") +
  geom_sf(data = arctic_utm) + 
  scale_color_distiller(palette = "RdYlBu", 
                        na.value = "transparent") + 
  labs(x = "", y = "", color = "SSS\n(psu)") + 
  coord_sf(xlim = c(-1800, 817.1763),
           ylim = c(10291.98, 12658.29), expand = FALSE) +                      # Canadian Archipelago
  scale_x_continuous(breaks = seq(-1800, 817.1763, by = 100)) +
  theme(plot.margin = margin(t = 0, l = 0, r = 0, b = 0))

gg <- ggpubr::ggarrange(gg_arctic, gg_canada, 
                        common.legend = TRUE, legend = "right")
# for (ext in extension) {
#   ggsave(plot = gg, paste0(path, "plots/Arctic_noobs", ext),
#          width = 9.3, height = 4)
# }

#------------------------------------------------------------------------------#
# number of neighbors
m <- 15

# ordering 
ord <- order(-coords_tr[,1]*coords_tr[,2]) # makes sense physically 
gg_order <- arctic_utm %>% 
  ggplot() +
  geom_point(data = data.frame(coords_tr[ord,], ord = 1:n_tr),
             aes(easting, northing, color = ord)) + 
  geom_sf(data = ice_utm, fill = "skyblue") + 
  geom_sf_label(data = ice_utm, label = "sea-ice") + 
  geom_sf() +
  labs(x = "", y = "", color = "Order") + 
  theme(plot.margin = margin(t = 0, l = -10, r = -2, b = -10), 
        legend.margin = margin(b = 0, r = 0, t = 0, l = -2))
# for (ext in extension) {
#   ggsave(plot = gg_order, paste0(path, "plots/Arctic_currents", ext),
#          width = 5, height = 3.7)
# }

# priors
## variogram
sample_vario <- gstat::variogram(ssss ~ 1, data = sss_utm)
vario_matern <- gstat::vgm(psill = 0.3, 
                           model = "Mat", 
                           range = 2000, 
                           nugget = 0.1, 
                           kappa = 1)
fit_matern <- gstat::fit.variogram(sample_vario, vario_matern)
# png(file = paste0(path, "plots/Arctic_variogram.png"),
#     width = 6.3, height = 3.7, units = "in", res = 720)
# pdf(file = paste0(path, "plots/Arctic_variogram.pdf"),
#     width = 6.3, height = 3.7)
# par(mar = rep(0, 4))
# plot(sample_vario, fit_matern,
#      xlab = "Euclidean distance", ylab = "Semivariance",
#      pch = 19, col = 1)
# dev.off()

## estimated starting value
sigma.sq <- fit_matern[2,2]
phi <- 1/fit_matern[2,3]
nu <- fit_matern[2,4]
tau.sq <- fit_matern[1,2]
phi.low <- 0.0001
phi.high <- 0.0031

starting <- list("phi" = phi, "sigma.sq" = sigma.sq, 
                 "tau.sq" = tau.sq, "nu" = nu)
tuning <- list("phi" = 0.0005, "sigma.sq" = 0.05, "tau.sq" = 0.05, "nu" = 0.05) 
priors <- list("phi.Unif" = c(phi.low, phi.high),
               "sigma.sq.IG" = c(2, sigma.sq), "tau.sq.IG" = c(2, tau.sq), 
               "nu.Unif" = c(nu-0.5,nu+0.5))

# MCMC
n.samples <- 25000
burn <- 5000
thin <- 20
save <- (n.samples - burn)/thin

# #------#
# # NNGP #
# #------#
# set.seed(123)
# m.s <- spNNGP(y_tr ~ 1, coords = coords_tr, starting = starting,
#               method = "response", n.neighbors = m,
#               tuning = tuning, priors = priors, cov.model = "matern",
#               n.samples = n.samples, n.omp.threads = 10, ord = ord,
#               return.neighbor.info = TRUE, verbose = TRUE)
# nninfo <- m.s$neighbor.info
# p.s <- predict(m.s, X.0 = matrix(1, nrow = n_pred),
#                coords.0 = as.matrix(coords_pred),
#                sub.sample = list(start = burn+1, end = n.samples, thin = thin),
#                n.omp.threads = 10, verbose = TRUE)
# 
# saveRDS(list(m.s = m.s, nninfo = nninfo, p.s = p.s),
#         paste0(path, "application/arctic_sss_NNGP.RDS"))

NNGPres <- readRDS(paste0(path, "application/arctic_sss_NNGP.RDS"))

# posterior summary
timeNNGP <- NNGPres$m.s$run.time + NNGPres$p.s$run.time
betaNNGP <- mean(NNGPres$m.s$p.beta.samples[burn + thin*1:save])
thetaNNGP <- colMeans(NNGPres$m.s$p.theta.samples[burn + thin*1:save,])

ystarNNGP <- rowMeans(NNGPres$p.s$p.y.0)
yquantNNGP <- apply(NNGPres$p.s$p.y.0, 1,
                    function(x) quantile(x, probs = c(0.025, 0.975)))
orgystarNNGP <- rowMeans(NNGPres$p.s$p.y.0^2)

# #---------#
# # BORA-GP #
# #---------#
# # barrier neighbor info
# # first order neighbors are found for R and U
# first_time <- system.time(
#   first <- barrier_neighbor(coords = coords_tr,
#                             coords.0 = coords_pred,
#                             ord = ord,
#                             n.neighbors = m,
#                             barrier = arctic_utm,
#                             cores = 20,
#                             verbose = TRUE,
#                             debug = list(barrier_n.indx = NULL,
#                                          barrier_dist = NULL,
#                                          barrier_nn.indx.0_list = NULL,
#                                          barrier_dist0 = NULL,
#                                          ref_window = NULL,
#                                          nonref_window = NULL,
#                                          ref_fill = FALSE,
#                                          nonref_fill = FALSE))
# )
# 
# second_time <- system.time(
#   second <- barrier_neighbor(coords = coords_tr,
#                              coords.0 = coords_pred,
#                              ord = ord,
#                              n.neighbors = m,
#                              barrier = barrier_utm,
#                              cores = 20,
#                              verbose = TRUE,
#                              debug = list(barrier_n.indx =
#                                             first$barrier_n.indx,
#                                           barrier_dist =
#                                             first$barrier_dist,
#                                           barrier_nn.indx.0_list =
#                                             first$barrier_nn.indx.0_list,
#                                           barrier_dist0 =
#                                             first$barrier_dist0,
#                                           ref_window = NULL,
#                                           nonref_window = NULL,
#                                           ref_fill = TRUE,
#                                           nonref_fill = TRUE))
# )
# 
# barrier_nb <- second
# barrier_nn.indx.0_list <- barrier_nb$barrier_nn.indx.0_list
# barrier_dist0 <- barrier_nb$barrier_dist0
# barrier_n.indx <- barrier_nb$barrier_n.indx
# barrier_dist <- barrier_nb$barrier_dist
# 
# # create a list for neighbor info
# barrier_nn.indx <- as.integer(unlist(barrier_n.indx)[-1]-1)
# barrier_nn.indx.lu <- c(cumsum(c(0,0,sapply(barrier_n.indx, length)[-1]))[1:n_tr],
#                         c(0,sapply(barrier_n.indx, length)[-1])) %>% as.integer()
# barrier_nninfo <- list(type = "barrier",
#                        n.indx = barrier_n.indx,
#                        n.neighbors = m,
#                        nn.indx = barrier_nn.indx,
#                        nn.indx.lu = barrier_nn.indx.lu,
#                        ord = ord)
# barrier_nn.indx.0 <- do.call("rbind", barrier_nn.indx.0_list)

# set.seed(123)
# barrier_m.s <- spNNGP(y_tr ~ 1, coords = coords_tr, starting = starting,
#                       method = "response", n.neighbors = m,
#                       tuning = tuning, priors = priors, cov.model = "matern",
#                       n.samples = n.samples, n.omp.threads = 10,
#                       neighbor.info = barrier_nninfo,
#                       return.neighbor.info = TRUE)
# barrier_p.s <- predict(barrier_m.s, X.0 = matrix(1, nrow = n_pred),
#                        coords.0 = as.matrix(coords_pred),
#                        sub.sample =
#                          list(start = burn+1, end = n.samples, thin = thin),
#                        nn.indx.0 = barrier_nn.indx.0, n.omp.threads = 10)
# 
# saveRDS(list(barrier_m.s = barrier_m.s,
#              barrier_nb = barrier_nb,
#              barrier_nn_time = first_time + second_time,
#              barrier_p.s = barrier_p.s),
#         paste0(path, "application/arctic_sss_BORAGP.RDS"))

BRGPres <- readRDS(paste0(path, "application/arctic_sss_BORAGP.RDS"))

# posterior summary
timeBRGP <- BRGPres$barrier_m.s$run.time + BRGPres$barrier_p.s$run.time
betaBRGP <- mean(BRGPres$barrier_m.s$p.beta.samples[burn + thin*1:save])
thetaBRGP <- colMeans(BRGPres$barrier_m.s$p.theta.samples[burn + thin*1:save,])

ystarBRGP <- rowMeans(BRGPres$barrier_p.s$p.y.0)
yquantBRGP <- apply(BRGPres$barrier_p.s$p.y.0, 1,
                    function(x) quantile(x, probs = c(0.025, 0.975)))
orgystarBRGP <- rowMeans(BRGPres$barrier_p.s$p.y.0^2)

# ESS: effective sample size
betaess <- posterior::ess_basic(
  BRGPres$barrier_m.s$p.beta.samples[burn + thin*1:save]
  )
tausqess <- posterior::ess_basic(
  BRGPres$barrier_m.s$p.theta.samples[burn + thin*1:save,2]
  )
yess <- apply(BRGPres$barrier_p.s$p.y.0, 1, posterior::ess_basic)

# convergence
# trace plots and running mean
convp <- data.frame(tau_sq = c(cumsum(BRGPres$barrier_m.s$p.theta.samples[burn + thin*1:save, 2])/1:save,
                               BRGPres$barrier_m.s$p.theta.samples[burn + thin*1:save, 2]),
                    beta = c(cumsum(BRGPres$barrier_m.s$p.beta.samples[burn + thin*1:save])/1:save,
                             BRGPres$barrier_m.s$p.beta.samples[burn + thin*1:save]),
                    what = rep(c("Running-mean", "Trace"), each = save),
                    iter = rep(1:save, 2)) %>%
  pivot_longer(-c(what, iter), values_to = "draw", names_to = "param") %>%
  ggplot() +
  geom_line(aes(iter, draw)) +
  facet_grid(factor(param, 
                    labels = c(bquote(beta[0]), bquote(tau^2))) ~ what,
             labeller = label_parsed,
             scales = "free_y") +
  labs(x = "Iterations", y = "") +
  theme(plot.margin = margin(t = 1, l = -5, r = 0, b = 0), 
        strip.text.x = element_blank())
# for (ext in extension) {
#   ggsave(plot = convp, paste0(path, "plots/arctic_convergence", ext),
#          width = 7, height = 2.8)
# }

# plot
alldata <- tibble(coords_all, Observed = c(orgy_tr, rep(NA, n_pred)), 
           NNGP = c(orgy_tr, orgystarNNGP), 
           "BORA-GP" = c(orgy_tr, orgystarBRGP)) %>% 
  pivot_longer(-c("easting", "northing"), names_to = "what", 
               values_to = "sss") 

## Focus: Canadian Archipelago
gg_range <- alldata %>%
  ggplot() +
  facet_grid(~ factor(what, levels = c("Observed", "BORA-GP", "NNGP"), 
                      labels = c("SMAP", "BORA-GP", "NNGP"))) +
  geom_point(aes(easting, northing, color = sss)) +
  geom_sf(data = ice_utm, fill = "skyblue") +
  geom_sf(data = arctic_utm) +
  scale_color_distiller(palette = "RdYlBu",
                        na.value = "transparent",
                        limits = c(min(orgy_tr), max(orgy_tr))) +
  coord_sf(xlim = c(-1450, -200), 
           ylim = c(11000, 12500), expand = FALSE) +
  labs(x = "", y = "", color = "SSS\n(psu)") + 
  scale_y_continuous(breaks = c(65.5, 67.5)) +
  theme(plot.margin = margin(t = -20, l = -5, r = 0, b = -10), 
        legend.margin = margin(t = -10, l = 0, r = 0, b = 0),
        legend.position = "bottom")

gg_nn <- alldata %>%
  filter(what == "NNGP") %>% 
  ggplot() +
  facet_grid(~ what) +
  geom_point(aes(easting, northing, color = sss)) +
  geom_sf(data = ice_utm, fill = "skyblue") +
  geom_sf(data = arctic_utm) +
  scale_color_distiller(palette = "RdYlBu",
                        na.value = "transparent") +                              
  coord_sf(xlim = c(-1450, -200), 
           ylim = c(11000, 12500), expand = FALSE) +
  labs(x = "", y = "", color = "SSS\n(psu)") + 
  scale_y_continuous(breaks = c(65.5, 67.5)) +
  theme(plot.margin = margin(t = -20, l = -5, r = 1, b = -10), 
        legend.margin = margin(t = -10, l = 0, r = 0, b = 0),
        legend.position = "bottom")

gg_together <- ggpubr::ggarrange(gg_range, gg_nn, widths = c(3, 1.15), 
                                 common.legend = FALSE)
# for (ext in extension) {
#   ggsave(plot = gg_together,
#          paste0(path, "plots/arctic_pred", ext),                                # v2: na.value = "black"
#          width = 10, height = 4.5)
# }

# inspect neighbors
large_sss <- alldata %>%
  filter(what == "NNGP") %>% 
  mutate(idx = 1:(n_tr + n_pred)) %>% 
  filter(sss > max(orgy_tr)) %>% 
  pull(idx)

pred_idx <- large_sss[1]
nb_idx <- BRGPres$barrier_nb$barrier_nn.indx.0_list[[pred_idx - n_tr]]
nnidx_tt <- RANN::nn2(coords_tr, coords_pred, k = m)$nn.idx
nb_idx2 <- nnidx_tt[pred_idx - n_tr,]

coords_all[pred_idx,] %>% 
  ggplot() + 
  geom_point(aes(easting, northing)) +
  geom_point(data = coords_tr[nb_idx2,], aes(easting, northing), color = "red") +
  geom_sf(data = ice_utm, fill = "skyblue") +
  geom_sf_label(data = ice_utm, label = "sea-ice") + 
  geom_sf(data = arctic_utm)
