# barrier = land 
rm(list = ls())

# dependencies
library(tidyverse)
theme_set(theme_bw())
library(sf) 
sf_use_s2(FALSE)
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

# call additional functions
source(paste0(path, "sim/spNNGP_freeS.R"))

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
sss_utm2$ssss <- scale(sss_utm2$sss)
m_sss <- attr(sss_utm2$ssss, "scaled:center")
sd_sss <- attr(sss_utm2$ssss, "scaled:scale")
y_tr <- sss_utm2$ssss
coords_all <- rbind(coords_tr, coords_pred)

# MCMC
n.samples <- 25000
burn <- 5000
thin <- 2
save <- (n.samples - burn)/thin

# # reference set
# ## cov1: a few neighbors overlap
# # long_S <- seq(-179.874, 179.874, by = 1)
# # lat_S <- seq(56, 84, by = 2)
# # coords_longlat_S <- expand.grid(long_S, lat_S) %>%
# #   as.data.frame() %>%
# #   rename(long = Var1, lat = Var2)
# # coords_S_utm <- st_transform(st_as_sf(coords_longlat_S,
# #                                       coords = c("long", "lat"),
# #                                       crs = 4326), crs = crs)
# 
# ## cov2: most neighbors overlap
# # long_S <- seq(-179.874, 179.874, by = 0.5)
# # lat_S <- seq(57.5, 83.5, by = 1)
# # coords_longlat_S <- expand.grid(long_S, lat_S) %>%
# #   as.data.frame() %>%
# #   rename(long = Var1, lat = Var2)
# # coords_S_utm <- st_transform(st_as_sf(coords_longlat_S,
# #                                       coords = c("long", "lat"),
# #                                       crs = 4326), crs = crs)
# 
# ## cov3: neighbors rarely overlap
# long_S <- seq(-181, 181, by = 5)
# lat_S <- seq(57.5, 83.5, by = 0.7) # 1
# coords_longlat_S <- expand.grid(long_S, lat_S) %>%
#   as.data.frame() %>%
#   rename(long = Var1, lat = Var2)
# coords_S_utm <- st_transform(st_as_sf(coords_longlat_S,
#                                       coords = c("long", "lat"),
#                                       crs = 4326), crs = crs)
# 
# coords_S <- coords_S_utm %>%
#   st_coordinates() %>%
#   as.data.frame() %>%
#   rename(easting = X, northing = Y)
# rm(coords_longlat_S)
# 
# ## delete ones on or too close to the land or the ice in grid
# dlt_barrier2 <- unlist(st_is_within_distance(barrier_utm,
#                                              coords_S_utm$geometry,
#                                              dist = 9)) %>% unique()
# coords_S <- coords_S[-dlt_barrier2,]
# k <- nrow(coords_S)
# 
# # rm(dlt_all, dlt_overlap, dlt_barrier2, coords, coords_all_utm, sss_utm)
# 
# # number of neighbors
# m <- 15
# 
# # ordering
# ord <- order(-coords_S[,1]*coords_S[,2]) # makes sense physically
# arctic_utm %>%
#   ggplot() +
#   geom_point(data = data.frame(coords_S[ord,], ord = 1:k),
#              aes(easting, northing, color = ord)) +
#   geom_sf(data = ice_utm, fill = "skyblue") +
#   geom_sf_label(data = ice_utm, label = "sea ice") +
#   geom_sf() +
#   labs(x = "", y = "", color = "Order") +
#   coord_sf(xlim = c(-1450, -200),
#            ylim = c(11200, 12500), expand = FALSE) +
#   theme(plot.margin = margin(t = 0, l = -10, r = -2, b = -10),
#         legend.margin = margin(b = 0, r = 0, t = 0, l = -2))
# 
# #------------------------------------------------------------------------------#
# ## Focus: Canadian Archipelago
# # coords.T
# coords_tr_mini <- coords_tr %>%
#   filter(easting > -1450, easting < -200,
#          northing > 11200, northing < 12500)
# idx_mini <- coords_tr %>%
#   mutate(idx = 1:n_tr) %>%
#   filter(easting > -1450, easting < -200,
#          northing > 11200, northing < 12500) %>%
#   pull(idx)
# y_tr_mini <- sss_utm2$ssss[idx_mini]
# 
# # coords.S
# coords_S_mini <- coords_S %>%
#   filter(easting > -1450, easting < -200,
#          northing > 11200, northing < 12500)
# k_mini <- nrow(coords_S_mini) # 179
# ord_mini <- order(-coords_S_mini[,1]*coords_S_mini[,2])
# 
# # coords.0
# coords_pred_mini <- coords_pred %>%
#   filter(easting > -1450, easting < -200,
#          northing > 11200, northing < 12500)
# coords_0_mini <- rbind(coords_tr_mini, coords_pred_mini)
# 
# # nothing overlaps
# unique(rbind(coords_0_mini, coords_S_mini)) %>% nrow() ==
#   nrow(coords_0_mini) + k_mini
# 
# # coords_tr_mini %>%
# #   ggplot() +
# #   geom_point(aes(easting, northing)) +
# #   geom_sf(data = arctic_utm)
# # data.frame(coords_S_mini[ord_mini,], ord = 1:k_mini) %>%
# #   ggplot() +
# #   geom_point(aes(easting, northing, color = ord)) +
# #   geom_sf(data = arctic_utm)
# 
# # covariance
# BRGPres <- readRDS(paste0(path, "application/arctic_sss_BORAGP_processed.RDS"))
# thetaBRGP <- colMeans(BRGPres$barrier_m.s$p.theta.samples[burn + thin*1:save,])
# 
# NNGPres <- readRDS(paste0(path, "application/arctic_sss_NNGP_processed.RDS"))
# thetaNNGP <- colMeans(NNGPres$m.s$p.theta.samples[burn + thin*1:save,])
# 
# ## NNGP neighbors
# # list containing ord and n.indx
# coords_S_ord_mini <- coords_S_mini[ord_mini,]
# nninfo_mini <- list()
# nninfo_mini$ord <- ord_mini
# nninfo_mini$n.indx <- list()
# nninfo_mini$n.indx[[1]] <- NA
# for (i in 2:k_mini) {
#   nninfo_mini$n.indx[[i]] <- RANN::nn2(data = coords_S_ord_mini[1:(i-1),],      # where
#                                        query = coords_S_ord_mini[i,],           # whose
#                                        k = min(i-1, m))$nn.idx %>%
#     as.numeric()
# }
# Ctilde <- boraGP:::create_Ctilde(coords = coords_S_mini,
#                                  neighbor.info = nninfo_mini,
#                                  sig_sq = 1,
#                                  phi = thetaNNGP[3],
#                                  nu = thetaNNGP[4],
#                                  base_cov = "matern", cores = 20)
# 
# ## BORA-GP neighbors
# barrier_nb_mini <- barrier_neighbor(coords = coords_S_mini,
#                                     coords.0 = coords_0_mini,
#                                     ord = ord_mini,
#                                     n.neighbors = m,
#                                     barrier = barrier_utm,
#                                     cores = 20,
#                                     verbose = TRUE,
#                                     debug = list(barrier_n.indx = NULL,
#                                                  barrier_dist = NULL,
#                                                  barrier_nn.indx.0_list = NULL,
#                                                  barrier_dist0 = NULL,
#                                                  ref_window = NULL,
#                                                  nonref_window = NULL,
#                                                  ref_fill = TRUE,
#                                                  nonref_fill = FALSE))
# barrier_nninfo_mini <- list(type = "barrier",
#                             n.indx = barrier_nb_mini$barrier_n.indx,
#                             ord = ord_mini)
# 
# ### point of interest
# j <- 215 + nrow(coords_tr_mini)
# 
# ### replace jth nn.indx.0 with first and second order neighbors
# barrier_nn.indx.0_list <- barrier_nb_mini$barrier_nn.indx.0_list                # cov calculation with barrier_nb_mini$barrier_nn.indx.0_list
# # barrier_0nb_mini <- barrier_neighbor(coords = coords_S_mini,
# #                                     coords.0 = coords_0_mini[j,],
# #                                     ord = ord_mini,
# #                                     n.neighbors = m,
# #                                     barrier = barrier_utm,
# #                                     cores = 20,
# #                                     verbose = TRUE,
# #                                     debug = list(barrier_n.indx =
# #                                                    barrier_nb_mini$barrier_n.indx,
# #                                                  barrier_dist =
# #                                                    barrier_nb_mini$barrier_dist,
# #                                                  barrier_nn.indx.0_list = NULL,
# #                                                  barrier_dist0 = NULL,
# #                                                  ref_window = NULL,
# #                                                  nonref_window = NULL,
# #                                                  ref_fill = TRUE,
# #                                                  nonref_fill = TRUE))
# # barrier_nn.indx.0_list[[j]] <-
# #   barrier_0nb_mini$barrier_nn.indx.0_list[[1]]                                # neighbor plot with barrier_nn.indx.0_list
# 
# barrier_Ctilde <- boraGP:::create_Ctilde(coords = coords_S_mini,
#                                          neighbor.info = barrier_nninfo_mini,
#                                          sig_sq = 1,
#                                          phi = thetaBRGP[3],
#                                          nu = thetaBRGP[4],
#                                          base_cov = "matern", cores = 20)
# 
# ## preprocessing
# nnidx_0_mini <- RANN::nn2(coords_S_mini, coords_0_mini, k = m)$nn.idx
# nn.indx.0_ord_mini <- RANN::nn2(coords_S_ord_mini, coords_0_mini, k = m)$nn.idx
# 
# barrier_nn.indx.0_ord_mini <- list()
# for (i in 1:nrow(coords_0_mini)) {
#   barrier_nn.indx.0_ord_mini[[i]] <-
#     order(ord_mini)[barrier_nn.indx.0_list[[i]]]
# }
# 
# ptnn <- nnidx_0_mini[j,]
# barrier_ptnn <- barrier_nn.indx.0_list[[j]]
# NNGP_cov <- boraGP:::NGPcov_m(v1 = coords_0_mini[j,],
#                               v2_mat = coords_0_mini,
#                               coords = coords_S_mini,
#                               neighbor.info = nninfo_mini,
#                               sig_sq = 1,
#                               phi = thetaNNGP[3],
#                               nu = thetaNNGP[4],
#                               base_cov = "matern",
#                               Ctilde = Ctilde,
#                               coords.0 = coords_0_mini,
#                               nn.indx.0_ord = nn.indx.0_ord_mini)
# boraGP_cov <- boraGP:::NGPcov_m(v1 = coords_0_mini[j,],
#                                 v2_mat = coords_0_mini,
#                                 coords = coords_S_mini,
#                                 neighbor.info = barrier_nninfo_mini,
#                                 sig_sq = 1,
#                                 phi = thetaBRGP[3],
#                                 nu = thetaBRGP[4],
#                                 base_cov = "matern",
#                                 Ctilde = barrier_Ctilde,
#                                 coords.0 = coords_0_mini,
#                                 nn.indx.0_ord = barrier_nn.indx.0_ord_mini)
# 
# gg_cov <- tibble(coords_0_mini,
#                  SMAP = rep(NA, nrow(coords_0_mini)),
#                  NNGP = NNGP_cov,
#                  "BORA-GP" = boraGP_cov) %>%
#   pivot_longer(-c(easting, northing), values_to = "cor", names_to = "what") %>%
#   filter(cor > 0.1) %>%
#   ggplot() +
#   facet_grid(~ factor(what, levels = c("SMAP", "BORA-GP", "NNGP"))) +
#   geom_point(aes(easting, northing, col = cor)) +
#   geom_point(data = data.frame(coords_S_mini[ptnn, ], what = "SMAP"),
#              aes(easting, northing), shape = 19, color = "#888888", size = 2) +
#   geom_point(data = data.frame(coords_S_mini[barrier_ptnn,], what = "SMAP"),
#              aes(easting, northing), shape = 17, color = "#DC3220") +
#   scale_color_scico(palette = "batlow", limits = c(0.1, 1),
#                     na.value = "transparent") +
#   geom_point(data = coords_0_mini[j,], aes(easting, northing),
#              col = "black", size = 2) +
#   geom_sf(data = ice_utm, fill = "skyblue") +
#   geom_sf(data = arctic_utm) +
#   labs(x = "", y = "", color = "Cor") +
#   coord_sf(xlim = c(-1450, -200),
#            ylim = c(11200, 12500), expand = FALSE) +
#   scale_y_continuous(breaks = c(65.5, 67.5)) +
#   theme(plot.margin = margin(t = 0, l = -5, r = 0, b = -5),
#         legend.margin = margin(t = 0, l = -5, r = 0, b = 0),
#         legend.position = "right")
# 
# saveRDS(list(barrier_nninfo_mini = barrier_nninfo_mini,
#              barrier_nn.indx.0_list = barrier_nn.indx.0_list,
#              gg_cov = gg_cov),
#         paste0(path, "application/arctic_sss_freeS_cov3.RDS"))
  
#------# 
# plot #
#------#
cov1 <- readRDS(paste0(path, "application/arctic_sss_freeS_cov1.RDS"))
cov2 <- readRDS(paste0(path, "application/arctic_sss_freeS_cov2.RDS"))
cov3 <- readRDS(paste0(path, "application/arctic_sss_freeS_cov3.RDS"))

sum(is.null(boraGP_cov))

gg_cov_final <- ggpubr::ggarrange(cov1$gg_cov, cov3$gg_cov, nrow = 2) 
for (ext in extension) {
  ggsave(plot = gg_cov_final, paste0(path, "plots/arctic_cov_freeS", ext),
         width = 8, height = 6)
}






