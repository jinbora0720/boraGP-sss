# Sliding doors
rm(list = ls())

# dependencies
library(tidyverse)
theme_set(theme_bw())
library(sf)
library(fields)
library(rgeos)
library(rgdal)
library(boraGP)

# path 
path <- "~/boraGP-sss/"

# plot extensions
extension <- c(".pdf", ".png")

# create dam door polygons
# output is a square
local.square.polygon = function(xlim, ylim){
  xlim = range(xlim); ylim = range(ylim)
  corner1 = c(xlim[1], ylim[2])
  corner2 = c(xlim[2], ylim[1])
  poly = Polygon(rbind(corner1, c(corner1[1], corner2[2]), corner2, 
                       c(corner2[1], corner1[2]), corner1), hole = FALSE)
  return(SpatialPolygons(list(Polygons(list(poly), ID = runif(1)))))
}

# the width of the opening in the barrier
smalldist <- 1

# the width/thickness of the barrier
width <- 0.8
poly1 <- local.square.polygon(xlim = c(2, 5-smalldist/2), 
                              ylim = 5 + width*c(-.5, .5))
poly2 <- local.square.polygon(xlim = c(5+smalldist/2, 8), 
                              ylim = 5 + width*c(-.5, .5))
poly.original <- SpatialPolygons(c(poly1@polygons, poly2@polygons))
slidingdoor <- st_as_sf(poly.original)
slidingdoors <- st_make_valid(st_combine(slidingdoor$geometry))
slidingdoors %>% ggplot() + geom_sf()
bbox <- st_bbox(c(xmin = 2, xmax = 8, ymin = 2, ymax = 8))

# create test points
ngrid <- 57
locs <- expand.grid(seq(2, 8, length = ngrid), seq(2, 8, length = ngrid))
colnames(locs) <- c("easting", "northing")
locs_sf <- st_as_sf(locs, coords = c("easting", "northing"))

## not on doors or too close to doors
dlt_idx <- c(unlist(st_intersects(slidingdoors, locs_sf$geometry)), 
             unlist(st_is_within_distance(slidingdoors, locs_sf$geometry, dist = 0.03)))
locs <- locs[-dlt_idx, ]
locs_sf <- locs_sf[-dlt_idx,]

# create training points
ngrid_tr <- 15
locs_tr <- expand.grid(seq(2, 8, length = ngrid_tr), seq(2, 8, length = ngrid_tr))
colnames(locs_tr) <- c("easting", "northing")
locs_tr_sf <- st_as_sf(locs_tr, coords = c("easting", "northing"))

## not on doors or too close to doors
dlt_idx <- c(unlist(st_intersects(slidingdoors, locs_tr_sf$geometry)), 
             unlist(st_is_within_distance(slidingdoors, locs_tr_sf$geometry, dist = 0.03)))
locs_tr <- locs_tr[-dlt_idx, ]
locs_tr_sf <- locs_tr_sf[-dlt_idx,]

# put together 
coords_all <- unique(rbind(locs_tr, locs))

n <- nrow(coords_all)
n_tr <- nrow(locs_tr)
n_tt <- n-n_tr

coords_tr <- coords_all[1:n_tr, ]
coords_tt <- coords_all[-(1:n_tr), ]

# ordering 
ord <- order(coords_tr[,2])
coords_sf <- st_as_sf(coords_all, coords = c("easting", "northing"))

# number of neighbors
## make sure your ordering has s1, ..., s_m+1 not blocked one another by a significant barrier
m <- 10 
slidingdoors %>% ggplot() + geom_sf() + 
  geom_sf(data = coords_sf[ord[1:(m+1)],])

#---------# 
# BORA-GP #
#---------# 
# barrier neighbor info
barrier_nninfo_all <- barrier_neighbor(coords = coords_tr,
                                       coords.0 = coords_tt,
                                       ord = ord,
                                       n.neighbors = m,
                                       barrier = slidingdoors,
                                       cores = 20,
                                       verbose = TRUE,
                                       debug = list(barrier_n.indx = NULL,
                                                    barrier_dist = NULL,
                                                    barrier_nn.indx.0_list = NULL,
                                                    barrier_dist0 = NULL,
                                                    ref_window = NULL,
                                                    nonref_window = NULL,
                                                    ref_fill = TRUE,
                                                    nonref_fill = TRUE))
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

#################
# Cov behaviour #
#################
## parameters estimated from full GP
sig_sq <- 1
phi <- 0.5 # effective range = 3/phi
dist_all <- as.matrix(dist(coords_all))

# stationary covariance
Cov <- sig_sq*exp(-phi*dist_all)

# nonstationary covariance
barrier_Ctilde <- boraGP:::create_Ctilde(coords = coords_tr, 
                                         neighbor.info = barrier_nninfo,
                                         sig_sq = sig_sq, phi = phi, 
                                         base_cov = "exponential")
barrier_nn.indx.0_ord <- list()
for (i in 1:length(barrier_nn.indx.0_list)) {
  barrier_nn.indx.0_ord[[i]] <- order(ord)[barrier_nn.indx.0_list[[i]]]
}

# paper plot 
middleline <- which(coords_tt[,1] == 5)
selectmiddleline <- middleline[c(14, 21, 42)]
boraGP_Ctilde <- list()
for (i in 1:length(selectmiddleline)) {
  pt <- selectmiddleline[i]
  boraGP_Ctilde[[i]] <- boraGP:::NGPcov_m(v1 = coords_all[n_tr + pt,], 
                                          v2_mat = coords_all,
                                          coords = coords_tr, 
                                          neighbor.info = barrier_nninfo,
                                          sig_sq = sig_sq, phi = phi,
                                          base_cov = "exponential", 
                                          Ctilde = barrier_Ctilde,
                                          coords.0 = coords_tt, 
                                          nn.indx.0_ord = barrier_nn.indx.0_ord)
}
gg <- data.frame(rbind(coords_all, coords_all, coords_all, 
                       coords_all, coords_all, coords_all), 
                 Cov = c(Cov[n_tr + selectmiddleline[1], ], 
                         Cov[n_tr + selectmiddleline[2], ], 
                         Cov[n_tr + selectmiddleline[3], ], 
                         boraGP_Ctilde[[1]], 
                         boraGP_Ctilde[[2]], 
                         boraGP_Ctilde[[3]]), 
                 Type = rep(c("Stationary covariance", "BORA-GP covariance"), 
                            each = 3*n), 
                 Pt = c(rep(1:3, each = n), rep(1:3, each = n))) %>% 
  ggplot() + 
  geom_sf(data = slidingdoors) +
  facet_grid(Type ~ Pt) + 
  geom_contour_filled(aes(x = easting, y = northing, z = Cov), 
                      breaks = c(seq(0, 0.6, by = 0.1), 0.8, 1)) +
  scale_fill_scico_d(palette = "batlow", begin = 0, end = 1) + 
  geom_point(data = data.frame(coords_tr, Type = "BORA-GP covariance"), 
             aes(easting, northing),
             size = 0.5, color = "grey15", shape = 4) +
  geom_point(data = data.frame(coords_all[n_tr + selectmiddleline,], Pt = 1:3), 
             aes(easting, northing), size = 2) + 
  labs(x = "", y = "", fill = "Cov") + 
  theme(aspect.ratio = 1,
        plot.margin=margin(t = 1, l = -3, r = 1, b = -10),
        legend.margin=margin(b = 0, r = 0, t = 0, l = 0), 
        strip.background.x = element_blank(),
        strip.text.x = element_blank())
for (ext in extension) {
  ggsave(plot = gg, paste0(path, "plots/cov_slidingdoors", ext),
         width = 7, height = 3.9)
}
