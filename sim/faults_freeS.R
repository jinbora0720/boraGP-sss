# check if spNNGP_freeS works as designed 
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

# call additional functions
source(paste0(path, "sim/spNNGP_freeS.R"))

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
width <- 0.0001
poly1 <- local.square.polygon(xlim = c(0.3, 2), 
                              ylim = 0.6 + width*c(-1, 1))
poly2 <- local.square.polygon(xlim = c(0.7, 2), 
                              ylim = 1.5 + width*c(-1, 1))
poly.original <- SpatialPolygons(c(poly1@polygons, poly2@polygons))
fault <- st_as_sf(poly.original)
faults <- st_make_valid(st_combine(fault$geometry))
bbox <- st_bbox(c(xmin = 0, xmax = 2, ymin = 0, ymax = 2))

# create test points
ngrid <- 34
locs <- expand.grid(seq(0, 2, length = ngrid), seq(0, 2, length = ngrid))
colnames(locs) <- c("easting", "northing")

# create training points
locs_tr <- expand.grid(seq(0, 2, length = 12), seq(0, 2, length = 12))
colnames(locs_tr) <- c("easting", "northing")

# put together 
coords_all <- unique(rbind(locs_tr, locs))
n <- nrow(coords_all)

# neighbor 
m <- 15

# training vs. test
n_tr <- nrow(locs_tr)
n_tt <- n-n_tr
coords_tr <- coords_all[1:n_tr, ]
coords_tt <- coords_all[-(1:n_tr), ]

# reference 
coords_S <- expand.grid(seq(0.01, 1.99, length = 10), seq(0.01, 1.99, length = 10))
colnames(coords_S) <- c("easting", "northing")
k <- nrow(coords_S)

## nothing overlaps
unique(rbind(coords_tr, coords_tt, coords_S)) %>% nrow() == n + k

# ordering
ord_all <- order(coords_all[,2])
ord <- order(coords_S[,2])
ord_tr <- order(coords_tr[,2])

#---------------#
# true function #
#---------------#
# barrier neighbor info for all as training
barrier_neighbor_time <- system.time(
  barrier_nninfo_truth <- barrier_neighbor(coords = coords_all,
                                           coords.0 = NULL,
                                           ord = ord_all,
                                           n.neighbors = m,
                                           barrier = faults,
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
)

barrier_n.indx_truth <- barrier_nninfo_truth$barrier_n.indx
barrier_nn.indx_truth <- as.integer(unlist(barrier_n.indx_truth)[-1]-1)
barrier_nn.indx.lu_truth <- c(cumsum(c(0,0,sapply(barrier_n.indx_truth, length)[-1]))[1:n],
                              c(0,sapply(barrier_n.indx_truth, length)[-1])) %>%
  as.integer()
barrier_nninfo_truth <- list(type = "barrier",
                             n.indx = barrier_n.indx_truth,
                             n.neighbors = m,
                             nn.indx = barrier_nn.indx_truth,
                             nn.indx.lu = barrier_nn.indx.lu_truth,
                             ord = ord_all)

# true parameters
true_sig_sq <- 1
true_nu <- 1.5
true_phi <- 4
true_tau_sq <- 0.1
true_beta <- c(1, .5)

#------------------#
# MCMC preparation #
#------------------#
# mcmc
## may have to increase # of iterations for stabler results
## this script is simply to check if the code runs well.
n.samples <- 2000
burn <- 1000

# priors
sigma.sq <- true_sig_sq
tau.sq <- true_tau_sq
nu <- true_nu
phi <- true_phi

## matern
maxdist <- max(dist(as.matrix(coords_tr[,c("easting", "northing")])))
d.low <- 0.25*maxdist
d.high <- 0.75*maxdist
dist_mat <- matrix(c(0, d.low, d.high, d.low, 0, 0, d.high, 0, 0),
                   nrow = 3, ncol = 3, byrow = T)
phi.cand <- seq(0.1, 10, by=0.01)
cor.low = cor.high <- rep(0, length(phi.cand))
for (i in 1:length(phi.cand)) {
  cor <- Cov_matern(dist = dist_mat, sigmasq = 1, phi = phi.cand[i], nu = nu)
  cor.low[i] <- cor[1,2]
  cor.high[i] <- cor[1,3]
}
phi.high <- phi.cand[which.min(abs(cor.low - 0.05))]
phi.low <- phi.cand[which.min(abs(cor.high - 0.05))]

starting <- list("phi" = phi, "sigma.sq" = sigma.sq,
                 "tau.sq" = tau.sq, "nu" = nu)
tuning <- list("phi" = 0.5, "sigma.sq" = 0.1, "tau.sq" = 0.1, "nu" = 0)
priors <- list("phi.Unif" = c(phi.low, phi.high),
               "sigma.sq.IG" = c(2, sigma.sq), "tau.sq.IG" = c(2, tau.sq),
               "nu.Unif" = c(nu-0.5,nu+0.5))

# run model
set.seed(397)
x <- rnorm(n)
x_S <- rnorm(k)
data0 <- rboraGP(coords = coords_all,
                 neighbor.info = barrier_nninfo_truth,
                 X = cbind(1, x),
                 beta = true_beta,
                 sig_sq = true_sig_sq,
                 phi = true_phi,
                 nu = true_nu,
                 tau_sq = true_tau_sq,
                 base_cov = "matern")
y_tr <- data0$y[1:n_tr]
y_tt <- data0$y[-(1:n_tr)]
x_tr <- x[1:n_tr]
x_tt <- x[-(1:n_tr)]

# #---------------------------------#
# # BORA-GP with Robust Adaptive MH #
# #---------------------------------#
# # barrier neighbor info for training vs. test w/ m
# time_nb <- system.time(
#   barrier_nninfo_all <- barrier_neighbor(coords = coords_S,
#                                          coords.0 = coords_all,
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
# barrier_nn.indx.lu <- c(cumsum(c(0,0,sapply(barrier_n.indx, length)[-1]))[1:k],
#                         c(0,sapply(barrier_n.indx, length)[-1])) %>% as.integer()
# barrier_nninfo <- list(type = "barrier",
#                        n.indx = barrier_n.indx,
#                        n.neighbors = m, nn.indx = barrier_nn.indx,
#                        nn.indx.lu = barrier_nn.indx.lu, ord = ord,
#                        nn.indx.0 = barrier_nn.indx.0,                           # added
#                        time_nb = time_nb)
# 
# set.seed(623)
# barrier_res <- spNNGP_freeS(y = y_tr, X = cbind(1, x_tr), coords.T = coords_tr, 
#                             X.S = cbind(1, x_S), coords.S = coords_S, 
#                             X.0 = cbind(1, x_tt), coords.0 = coords_tt,
#                             starting = starting, tuning = tuning, priors = priors,
#                             cov.model = "matern",
#                             sub.sample = list(start = burn+1, end = n.samples, thin = 1),
#                             cores = 20, neighbor.info = barrier_nninfo, 
#                             RAMH = TRUE, verbose = TRUE, predict = TRUE,
#                             debug = list(fix_cov_params = FALSE))
# 
# #------------------#
# # Original BORA-GP #
# #------------------#
# # barrier neighbor info for training vs. test w/ m
# org_time_nb <- system.time(
#   org_barrier_nninfo_all <- barrier_neighbor(coords = coords_tr,
#                                          coords.0 = coords_tt,
#                                          ord = ord_tr,
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
# org_barrier_nn.indx.0_list <- org_barrier_nninfo_all$barrier_nn.indx.0_list
# org_barrier_nn.indx.0 <- do.call("rbind", org_barrier_nn.indx.0_list)
# org_barrier_n.indx <- org_barrier_nninfo_all$barrier_n.indx
# org_barrier_nn.indx <- as.integer(unlist(org_barrier_n.indx)[-1]-1)
# org_barrier_nn.indx.lu <- c(cumsum(c(0,0,sapply(org_barrier_n.indx, length)[-1]))[1:n_tr],
#                         c(0,sapply(org_barrier_n.indx, length)[-1])) %>% as.integer()
# org_barrier_nninfo <- list(type = "barrier",
#                            n.indx = org_barrier_n.indx,
#                            n.neighbors = m, nn.indx = org_barrier_nn.indx,
#                            nn.indx.lu = org_barrier_nn.indx.lu, ord = ord_tr,
#                            nn.indx.0 = org_barrier_nn.indx.0,                           # added
#                            time_nb = org_time_nb)
# 
# set.seed(623)
# org_barrier_m.s <- spNNGP(y_tr ~ x_tr, coords = coords_tr, starting = starting,
#                           method = "response", n.neighbors = m,
#                           tuning = tuning, priors = priors, cov.model = "matern",
#                           n.samples = n.samples, n.omp.threads = 10,
#                           neighbor.info = org_barrier_nninfo, verbose = TRUE)
# org_barrier_p.s <- predict(org_barrier_m.s, X.0 = cbind(1, x_tt),
#                            coords.0 = as.matrix(coords_tt),
#                            sub.sample = list(start = burn+1,
#                                              end = n.samples, thin = 1),
#                            nn.indx.0 = org_barrier_nn.indx.0,
#                            n.omp.threads = 10, verbose = TRUE)
# 
# saveRDS(list(barrier_res = barrier_res, barrier_nninfo = barrier_nninfo,
#              org_barrier_m.s = org_barrier_m.s, org_barrier_p.s = org_barrier_p.s,
#              org_barrier_nninfo = org_barrier_nninfo),
#         paste0(path, "sim/faults_freeS.RDS"))

#---------#
# results #
#---------#
# time
res <- readRDS(paste0(path, "sim/faults_freeS.RDS"))
res$barrier_nninfo$time_nb
res$barrier_res$run.time
res$org_barrier_m.s$run.time + res$org_barrier_p.s$run.time

# convergence
true_beta
summary(res$barrier_res$p.beta.samples)
summary(res$org_barrier_m.s$p.beta.samples)

plot(res$barrier_res$p.beta.samples[,1])
plot(res$barrier_res$p.beta.samples[,2])

true_sig_sq; true_tau_sq; true_phi; true_nu
summary(res$barrier_res$p.theta.samples)
summary(res$org_barrier_m.s$p.theta.samples)

plot(res$barrier_res$p.theta.samples[,1])
plot(res$barrier_res$p.theta.samples[,2])
plot(res$barrier_res$p.theta.samples[,3])

# prediction
y0_mean <- rowMeans(res$barrier_res$p.y.0)
y0_quant <- apply(res$barrier_res$p.y.0, 1,
                  function(x) quantile(x, probs = c(0.025, 0.975)))
org_y0_mean <- rowMeans(res$org_barrier_p.s$p.y.0)
org_y0_quant <- apply(res$org_barrier_p.s$p.y.0, 1,
                      function(x) quantile(x, probs = c(0.025, 0.975)))
# plot
tibble(coords_all,
       Truth = c(y_tr, y_tt),
       "BORA-GP: S!=T" = c(y_tr, y0_mean),
       "BORA-GP: S=T" = c(y_tr, org_y0_mean)) %>%
  pivot_longer(-c("easting", "northing"), names_to = "model", values_to = "y") %>%
  ggplot() +
  facet_wrap(~factor(model,
                     levels = c("Truth", "BORA-GP: S=T", "BORA-GP: S!=T")),
             nrow = 1) +
  geom_raster(aes(x=easting, y=northing, fill=y)) +
  geom_sf(data = faults, col = "black", linewidth = 1) +
  scale_fill_scico(palette = "roma", direction = -1) +
  labs(x="", y="", fill = "y") +
  theme(plot.margin=margin(t=0,l=-10, r=0, b=-15),
        legend.margin=margin(b=0,r=0,t=0,l=-3))

# prediction errors
data.frame(y_tt = rep(y_tt, times = 2),
           yhat = c(y0_mean, org_y0_mean),
           lower = c(y0_quant[1,], org_y0_quant[1,]),
           upper = c(y0_quant[2,], org_y0_quant[2,]),
           class = rep(c("S!=T", "S=T"), each = n_tt)) %>%
  mutate(error = y_tt - yhat,
         width = upper-lower,
         cover = ifelse((y_tt > lower & y_tt < upper), 1, 0)) %>%
  group_by(class) %>%
  summarise(rmspe = sqrt(mean(error^2)),
            mape = mean(abs(error)),
            coverage = mean(cover),
            meanwidth = mean(width))
