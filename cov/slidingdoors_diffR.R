# Sliding doors
# different Rs
rm(list = ls())

# dependencies
library(tidyverse)
theme_set(theme_bw())
library(sf)
library(fields)
library(rgeos)
library(rgdal)
library(boraGP)
library(scico)

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
selectmiddleline_y <- seq(2, 8, length = ngrid)[c(19, 28, 56)]

## not on doors or too close to doors
dlt_idx <- c(unlist(st_intersects(slidingdoors, locs_sf$geometry)), 
             unlist(st_is_within_distance(slidingdoors, locs_sf$geometry, dist = 0.03)))
locs <- locs[-dlt_idx, ]
locs_sf <- locs_sf[-dlt_idx,]

# create training points 
J <- 5
locs_tr = locs_tr_sf <- list()
#-------------------------------#
# Scenario 1: randomly selected #
#-------------------------------#
set.seed(421)
idx_tr <- sort(sample.int(100^2, 15^2))
locs_tr[[1]] <- expand.grid(seq(2, 8, length = 100), seq(2, 8, length = 100))[idx_tr,]
colnames(locs_tr[[1]]) <- c("easting", "northing")
locs_tr_sf[[1]] <- st_as_sf(locs_tr[[1]], coords = c("easting", "northing"))

#-----------------------------------------------#
# Scenario 2: grid before and after the barrier #
#-----------------------------------------------#
locs_tr[[2]] <- expand.grid(seq(2, 8, length = 15)[-c(7:9)], 
                            seq(2, 8, length = 15)[-c(7:9)])
colnames(locs_tr[[2]]) <- c("easting", "northing")
locs_tr_sf[[2]] <- st_as_sf(locs_tr[[2]], coords = c("easting", "northing"))

#--------------------------#
# Scenario 3: sparser grid #
#--------------------------#
locs_tr[[3]] <- expand.grid(seq(2, 8, length = 15), seq(2, 8, length = 8))
colnames(locs_tr[[3]]) <- c("easting", "northing")
locs_tr_sf[[3]] <- st_as_sf(locs_tr[[3]], coords = c("easting", "northing"))

#-------------------------------#
# Scenario 4: much sparser grid #
#-------------------------------#
locs_tr[[4]] <- expand.grid(seq(2, 8, length = 8), seq(2, 8, length = 8))
colnames(locs_tr[[4]]) <- c("easting", "northing")
locs_tr_sf[[4]] <- st_as_sf(locs_tr[[4]], coords = c("easting", "northing"))

#---------------------------#
# Scenario 5: original grid #
#---------------------------#
locs_tr[[5]] <- expand.grid(seq(2, 8, length = 15), seq(2, 8, length = 15))
colnames(locs_tr[[5]]) <- c("easting", "northing")
locs_tr_sf[[5]] <- st_as_sf(locs_tr[[5]], coords = c("easting", "northing"))

## not on doors or too close to doors
coords_all = coords_tr = coords_tt = ord = coords_sf <- list()
n = n_tr = n_tt <- rep(0, J)
for (j in 1:J) {
  dlt_idx <- c(unlist(st_intersects(slidingdoors, locs_tr_sf[[j]]$geometry)), 
               unlist(st_is_within_distance(
                 slidingdoors, locs_tr_sf[[j]]$geometry, dist = 0.03)))
  if (length(dlt_idx) > 0) {
    locs_tr[[j]] <- locs_tr[[j]][-dlt_idx, ]
    locs_tr_sf[[j]] <- locs_tr_sf[[j]][-dlt_idx,]
  }
  
  # put together 
  coords_all[[j]] <- unique(rbind(locs_tr[[j]], locs))
  
  n[j] <- nrow(coords_all[[j]])
  n_tr[j] <- nrow(locs_tr[[j]])
  
  coords_tr[[j]] <- coords_all[[j]][1:n_tr[j], ]
  coords_tt[[j]] <- coords_all[[j]][-(1:n_tr[j]), ]
  
  # ordering 
  ord[[j]] <- order(coords_tr[[j]][,2]) 
  coords_sf[[j]] <- st_as_sf(coords_all[[j]], coords = c("easting", "northing"))
}
n_tt <- n-n_tr

# number of neighbors
## make sure your ordering has s1, ..., s_m+1 not blocked one another by a significant barrier
m <- 10 
for (j in 1:J) {
  plotm <- slidingdoors %>% ggplot() + geom_sf() + 
    geom_sf(data = coords_sf[[j]][ord[[j]][1:(m+1)],])
  plotR <- slidingdoors %>% ggplot() + geom_sf() + 
    geom_sf(data = coords_sf[[j]][ord[[j]][1:n_tr[j]],])
  gridExtra::grid.arrange(plotm, plotR, nrow = 1)
}

#---------# 
# BORA-GP #
#---------# 
barrier_nninfo = barrier_nn.indx.0_list <- list()
for (j in 1:J) {
  # barrier neighbor info
  barrier_nninfo_all <- barrier_neighbor(coords = coords_tr[[j]],
                                         coords.0 = coords_tt[[j]],
                                         ord = ord[[j]],
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
  barrier_nn.indx.0_list[[j]] <- barrier_nninfo_all$barrier_nn.indx.0_list
  barrier_nn.indx.0 <- do.call("rbind", barrier_nninfo_all$barrier_nn.indx.0_list)
  barrier_n.indx <- barrier_nninfo_all$barrier_n.indx
  barrier_nn.indx <- as.integer(unlist(barrier_n.indx)[-1]-1)
  barrier_nn.indx.lu <- c(cumsum(c(0,0,sapply(barrier_n.indx, length)[-1]))[1:n_tr[j]],
                          c(0,sapply(barrier_n.indx, length)[-1])) %>% as.integer()
  barrier_nninfo[[j]] <- list(type = "barrier",
                              n.indx = barrier_n.indx,
                              n.neighbors = m, nn.indx = barrier_nn.indx,
                              nn.indx.lu = barrier_nn.indx.lu, ord = ord[[j]])
}

#################
# Cov behaviour #
#################
## parameters estimated from full GP
sig_sq <- 1
phi <- 0.5 # effective range = 3/phi
dist_all <- as.matrix(dist(coords_all[[5]]))

# stationary covariance
Cov <- sig_sq*exp(-phi*dist_all)

# nonstationary covariance
boraGP_Ctilde <- list()
selectmiddleline <- matrix(0, ncol = J, nrow = 3)
for (j in 1:J) {
  barrier_Ctilde <- boraGP:::create_Ctilde(coords = coords_tr[[j]],
                                           neighbor.info = barrier_nninfo[[j]],
                                           sig_sq = sig_sq, phi = phi,
                                           base_cov = "exponential")
  barrier_nn.indx.0_ord <- list()
  for (i in 1:length(barrier_nn.indx.0_list[[j]])) {
    barrier_nn.indx.0_ord[[i]] <- order(ord[[j]])[barrier_nn.indx.0_list[[j]][[i]]]
  }
  
  # paper plot 
  middleline <- which(coords_tt[[j]][,1] == 5)
  selectmiddleline[,j] <- middleline[
    which(coords_tt[[j]][middleline,2] %in% selectmiddleline_y)]
  boraGP_Ctilde[[j]] <- matrix(0, 3, n[j])
  for (i in 1:3) {
    pt <- selectmiddleline[i,j]
    boraGP_Ctilde[[j]][i,] <- boraGP:::NGPcov_m(v1 = coords_all[[j]][n_tr[j] + pt,],
                                                v2_mat = coords_all[[j]],
                                                coords = coords_tr[[j]],
                                                neighbor.info = barrier_nninfo[[j]],
                                                sig_sq = sig_sq, phi = phi,
                                                base_cov = "exponential",
                                                Ctilde = barrier_Ctilde,
                                                coords.0 = coords_tt[[j]],
                                                nn.indx.0_ord = barrier_nn.indx.0_ord)
  }
}

# plot
gg <- list()
gg[[1]] <- data.frame(rbind(coords_tt[[1]], coords_tt[[1]], coords_tt[[1]]), 
                      Cov = c(boraGP_Ctilde[[1]][1,-c(1:n_tr[1])], 
                              boraGP_Ctilde[[1]][2,-c(1:n_tr[1])], 
                              boraGP_Ctilde[[1]][3,-c(1:n_tr[1])]), 
                      Type = rep("BORA-GP covariance", each = 3*n_tt[1]), 
                      Pt = rep(1:3, each = n_tt[1])) %>% 
  ggplot() + 
  geom_sf(data = slidingdoors) +
  facet_grid(Type ~ Pt) + 
  geom_contour_filled(aes(x = easting, y = northing, z = Cov), 
                      breaks = c(seq(0, 0.6, by = 0.1), 0.8, 1)) +
  scale_fill_scico_d(palette = "batlow", begin = 0, end = 1) + 
  geom_point(data = data.frame(coords_tr[[1]]), aes(easting, northing),
             size = 0.5, color = "grey15", shape = 4) +
  geom_point(data = data.frame(coords_all[[1]][n_tr[1] + selectmiddleline[,1],], 
                               Pt = 1:3), 
             aes(easting, northing), size = 2) + 
  labs(x = "", y = "", fill = "Cov") + 
  theme(aspect.ratio = 1,
        plot.margin=margin(t = 1, l = -3, r = 1, b = -10),
        legend.margin=margin(b = 0, r = 0, t = 0, l = 0), 
        strip.background.x = element_blank(),
        strip.text.x = element_blank())

for (j in 2:J) {
  gg[[j]] <- data.frame(rbind(coords_all[[j]], coords_all[[j]], coords_all[[j]]), 
                        Cov = c(boraGP_Ctilde[[j]][1,], 
                                boraGP_Ctilde[[j]][2,], 
                                boraGP_Ctilde[[j]][3,]), 
                        Type = rep("BORA-GP covariance", each = 3*n[j]), 
                        Pt = rep(1:3, each = n[j])) %>% 
    ggplot() + 
    geom_sf(data = slidingdoors) +
    facet_grid(Type ~ Pt) + 
    geom_contour_filled(aes(x = easting, y = northing, z = Cov), 
                        breaks = c(seq(0, 0.6, by = 0.1), 0.8, 1)) +
    scale_fill_scico_d(palette = "batlow", begin = 0, end = 1) + 
    geom_point(data = data.frame(coords_tr[[j]]), aes(easting, northing),
               size = 0.5, color = "grey15", shape = 4) +
    geom_point(data = data.frame(coords_all[[j]][n_tr[j] + selectmiddleline[,j],], 
                                 Pt = 1:3), 
               aes(easting, northing), size = 2) + 
    labs(x = "", y = "", fill = "Cov") + 
    theme(aspect.ratio = 1,
          plot.margin=margin(t = 1, l = -3, r = 1, b = -10),
          legend.margin=margin(b = 0, r = 0, t = 0, l = 0), 
          strip.background.x = element_blank(),
          strip.text.x = element_blank())
}
gg_stcov <- data.frame(rbind(coords_all[[j]], coords_all[[j]], coords_all[[j]]), 
                       Cov = c(Cov[n_tr[j] + selectmiddleline[1,j], ], 
                               Cov[n_tr[j] + selectmiddleline[2,j], ], 
                               Cov[n_tr[j] + selectmiddleline[3,j], ]), 
                       Type = rep("Stationary covariance", each = 3*n[j]), 
                       Pt = rep(1:3, each = n[j])) %>% 
  ggplot() + 
  geom_sf(data = slidingdoors) +
  facet_grid(Type ~ Pt) + 
  geom_contour_filled(aes(x = easting, y = northing, z = Cov), 
                      breaks = c(seq(0, 0.6, by = 0.1), 0.8, 1)) +
  scale_fill_scico_d(palette = "batlow", begin = 0, end = 1) + 
  geom_point(data = data.frame(coords_all[[j]][n_tr[j] + selectmiddleline[,j],], 
                               Pt = 1:3), 
             aes(easting, northing), size = 2) + 
  labs(x = "", y = "", fill = "Cov") + 
  theme(aspect.ratio = 1,
        plot.margin=margin(t = 1, l = -3, r = 1, b = -10),
        legend.margin=margin(b = 0, r = 0, t = 0, l = 0), 
        strip.background.x = element_blank(),
        strip.text.x = element_blank())

gg_all <- ggpubr::ggarrange(gg_stcov, gg[[5]], 
                            gg[[2]], gg[[3]], 
                            gg[[4]], gg[[1]], common.legend = T, nrow = 3, ncol = 2, 
                            legend = "right")
for (ext in extension) {
  ggsave(plot = gg_all, paste0(path, "plots/cov_slidingdoors_diffRs", ext), 
         width = 10.5, height = 5)
}
