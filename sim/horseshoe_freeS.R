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

# path 
path <- "~/boraGP-sss/"

# call additional functions
source(paste0(path, "sim/spNNGP_freeS.R"))

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
tau_sq <- 0
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
tau.sq <- 0.01^2
starting <- list("phi" = phi, "sigma.sq" = sigma.sq, 
                 "tau.sq" = tau.sq, "nu" = nu)
tuning <- list("phi" = 0.5, "sigma.sq" = 1, "tau.sq" = 5e-05, "nu" = 0)         # nu is fixed 
priors <- list("phi.Unif" = c(phi.low, phi.high),                               # approx. 3/(0.55*maxdist), 3/(0.19*maxdist)
               "sigma.sq.IG" = c(2, sigma.sq), "tau.sq.IG" = c(2, tau.sq), 
               "nu.Unif" = c(nu-0.5,nu+0.5))

# mcmc
n.samples <- 10000
burn <- 5000

# simulation variables 
m <- 20
n_tr <- 300
n_tt <- n - n_tr

# training data
set.seed(984)
idx_tr <- sample.int(n, n_tr)
coords_tr <- data_hs[idx_tr, c("easting", "northing")]
rownames(coords_tr) <- NULL
y_tr <- data_hs[idx_tr,"y"]

# test data
coords_tt <- data_hs[-idx_tr, c("easting", "northing")]
rownames(coords_tt) <- NULL
y_tt <- data_hs[-idx_tr,"y"]

#----------------------#
# reference set (grid) #
#----------------------#
xmS <- 32
ynS <- 15
xS <- seq(-0.99, 3.99, length = xmS)
yS <- seq(-0.99, 0.99, length = ynS)
xxS <- rep(xS, ynS)
yyS <- rep(yS, rep(xmS, ynS))
zS <- mgcv::fs.test(xxS, yyS)
zzS <- matrix(zS, xmS, ynS)

## nothing overlap?
c(x, xS) %>% unique() %>% length() == xm + xmS
c(y, yS) %>% unique() %>% length() == yn + ynS

## is.na(w) = T outside the horseshoe
dataS <- data.frame(easting = xxS, northing = yyS, w = as.numeric(zzS)) %>% 
  na.omit()
dataS <- dataS %>%  mutate(y = w + sqrt(tau_sq)*rnorm(nrow(dataS)))
coords_S <- dataS[, c("easting", "northing")]
k <- nrow(coords_S)

## ordering (default)
ord <- order(coords_S[,1])

#-------------------------------------#
# second reference set (coarser grid)
#-------------------------------------#
xmS2 <- 20
ynS2 <- 8
xS2 <- seq(-0.99, 3.99, length = xmS2)
yS2 <- seq(-0.99, 0.99, length = ynS2)
xxS2 <- rep(xS2, ynS2)
yyS2 <- rep(yS2, rep(xmS2, ynS2))
zS2 <- mgcv::fs.test(xxS2, yyS2)
zzS2 <- matrix(zS2, xmS2, ynS2)

## nothing overlap?
c(x, xS2) %>% unique() %>% length() == xm + xmS2
c(y, yS2) %>% unique() %>% length() == yn + ynS2

## is.na(w) = T outside the horseshoe
dataS2 <- data.frame(easting = xxS2, northing = yyS2, w = as.numeric(zzS2)) %>% 
  na.omit()
dataS2 <- dataS2 %>%  mutate(y = w + sqrt(tau_sq)*rnorm(nrow(dataS2)))
coords_S2 <- dataS2[, c("easting", "northing")]
k2 <- nrow(coords_S2)

## ordering (default)
ord2 <- order(coords_S2[,1])

# show locations 
data.frame(rbind(coords_tr, coords_S, coords_S2)) %>% 
  mutate(cat = rep(c("T", "S1", "S2"), c(n_tr, k, k2))) %>% 
  ggplot() + 
  geom_point(aes(easting, northing, shape = cat)) +
  scale_shape_manual(values = c(4, 3, 19), name = "") +
  geom_sf(data = hs_sf, fill = NA) + 
  labs(x = "", y = "")

# # run models
# #--------------------#
# # BORA-GP w/ S1 != T #
# #--------------------#
# time_nb <- system.time(
#   barrier_nninfo_all <- barrier_neighbor(coords = coords_S,
#                                          coords.0 = rbind(coords_tr, coords_tt),
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
# barrier_nn.indx.lu <- c(cumsum(c(0,0,sapply(barrier_n.indx, length)[-1]))[1:k],
#                         c(0,sapply(barrier_n.indx, length)[-1])) %>% as.integer()
# barrier_nninfo <- list(type = "barrier",
#                        n.indx = barrier_n.indx,
#                        n.neighbors = m, nn.indx = barrier_nn.indx,
#                        nn.indx.lu = barrier_nn.indx.lu, ord = ord,
#                        nn.indx.0 = barrier_nn.indx.0,                           # added
#                        time_nb = time_nb)
# set.seed(984)
# barrier_res <- spNNGP_freeS(y = y_tr, coords.T = coords_tr,
#                             coords.S = coords_S, coords.0 = coords_tt,
#                             starting = starting, priors = priors,
#                             cov.model = "matern",
#                             sub.sample = list(start = burn+1, end = n.samples, thin = 1),
#                             cores = 20, neighbor.info = barrier_nninfo,
#                             RAMH = TRUE, verbose = TRUE, predict = TRUE,
#                             debug = list(fix_cov_params = FALSE))
# 
# #--------------------#
# # BORA-GP w/ S2 != T #
# #--------------------#
# time_nb2 <- system.time(
#   barrier_nninfo_all2 <- barrier_neighbor(coords = coords_S2,
#                                          coords.0 = rbind(coords_tr, coords_tt),
#                                          ord = ord2,
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
# barrier_nn.indx.0_list2 <- barrier_nninfo_all2$barrier_nn.indx.0_list
# barrier_nn.indx.02 <- do.call("rbind", barrier_nn.indx.0_list2)
# barrier_n.indx2 <- barrier_nninfo_all2$barrier_n.indx
# barrier_nn.indx2 <- as.integer(unlist(barrier_n.indx2)[-1]-1)
# barrier_nn.indx.lu2 <- c(cumsum(c(0,0,sapply(barrier_n.indx2, length)[-1]))[1:k2],
#                         c(0,sapply(barrier_n.indx2, length)[-1])) %>% as.integer()
# barrier_nninfo2 <- list(type = "barrier",
#                         n.indx = barrier_n.indx2,
#                         n.neighbors = m, nn.indx = barrier_nn.indx2,
#                         nn.indx.lu = barrier_nn.indx.lu2, ord = ord2,
#                         nn.indx.0 = barrier_nn.indx.02,                         # added
#                         time_nb = time_nb2)
#                        
# set.seed(984)
# barrier_res2 <- spNNGP_freeS(y = y_tr, coords.T = coords_tr,
#                              coords.S = coords_S2, coords.0 = coords_tt,
#                              starting = starting, priors = priors,
#                              cov.model = "matern",
#                              sub.sample = list(start = burn+1, end = n.samples, thin = 1),
#                              cores = 20, neighbor.info = barrier_nninfo2,
#                              RAMH = TRUE, verbose = TRUE, predict = TRUE,
#                              debug = list(fix_cov_params = FALSE))
# 
# #------------------#
# # BORA-GP w/ S = T #
# #------------------#
# # ordering (default)
# ord_tr <- order(coords_tr[,1])
# 
# org_time_nb <- system.time(
#   org_barrier_nninfo_all <- barrier_neighbor(coords = coords_tr,
#                                          coords.0 = coords_tt,
#                                          ord = ord_tr,
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
#                            nn.indx.0 = org_barrier_nn.indx.0,                   # added
#                            time_nb = org_time_nb)
# set.seed(984)                  
# barrier_m.s <- spNNGP(y_tr ~ 1, coords = coords_tr, starting = starting,
#                       method = "response", n.neighbors = m,
#                       tuning = tuning, priors = priors,
#                       cov.model = "matern",
#                       n.samples = n.samples, n.omp.threads = 10,
#                       neighbor.info = org_barrier_nninfo, verbose = T)
# barrier_p.s <- predict(barrier_m.s, X.0 = matrix(1, nrow = n_tt, ncol = 1),
#                        coords.0 = as.matrix(coords_tt),
#                        sub.sample =
#                          list(start = burn+1, end = n.samples, thin = 1),
#                        nn.indx.0 = org_barrier_nn.indx.0,
#                        n.omp.threads = 10, verbose = T)
# 
# saveRDS(list(barrier_res = barrier_res, barrier_nninfo = barrier_nninfo, 
#              barrier_res2 = barrier_res2, barrier_nninfo2 = barrier_nninfo2,
#              org_barrier_m.s = barrier_m.s, org_barrier_p.s = barrier_p.s),
#         paste0(path, "sim/horseshoe_freeS.RDS"))
# 
# # y prediction
# y0_mean <- rowMeans(barrier_res$p.y.0)
# y0_quant <- apply(barrier_res$p.y.0, 1,
#                   function(x) quantile(x, probs = c(0.025, 0.975)))
# y0_mean2 <- rowMeans(barrier_res2$p.y.0)
# y0_quant2 <- apply(barrier_res2$p.y.0, 1,
#                    function(x) quantile(x, probs = c(0.025, 0.975)))
# org_y0_mean <- rowMeans(org_barrier_p.s$p.y.0)
# org_y0_quant <- apply(org_barrier_p.s$p.y.0, 1,
#                       function(x) quantile(x, probs = c(0.025, 0.975)))
# 
# saveRDS(list(y0_mean = y0_mean, y0_quant = y0_quant, 
#              y0_mean2 = y0_mean2, y0_quant2 = y0_quant2, 
#              org_y0_mean = org_y0_mean, org_y0_quant = org_y0_quant),
# paste0(path, "sim/horseshoe_freeS_res.RDS"))

# # time
# res <- readRDS(paste0(path, "sim/horseshoe_freeS.RDS"))
# 
# ## k = 303
# res$barrier_nninfo$time_nb
# res$barrier_res$run.time
# 
# ## k = 92
# res$barrier_nninfo2$time_nb
# res$barrier_res2$run.time
# 
# # convergence
# ## S = T
# summary(res$org_barrier_m.s$p.theta.samples)
# summary(res$org_barrier_m.s$p.beta.samples)
# 
# ## k = 303
# summary(res$barrier_res$p.theta.samples)
# summary(res$barrier_res$p.beta.samples)
# plot(res$barrier_res$p.theta.samples[,1])
# plot(res$barrier_res$p.theta.samples[,2])
# plot(res$barrier_res$p.theta.samples[,3])
# plot(res$barrier_res$p.beta.samples)
# 
# ## k = 92
# summary(res$barrier_res2$p.theta.samples)
# summary(res$barrier_res2$p.beta.samples)
# plot(res$barrier_res2$p.theta.samples[,1])
# plot(res$barrier_res2$p.theta.samples[,2])
# plot(res$barrier_res2$p.theta.samples[,3])
# plot(res$barrier_res2$p.beta.samples)

# read summary results
y0s <- readRDS(paste0(path, "sim/horseshoe_freeS_res.RDS"))

data.frame(y_tt = rep(y_tt, times = 3),
           yhat = c(y0s$y0_mean, y0s$y0_mean2, y0s$org_y0_mean),
           lower = c(y0s$y0_quant[1,], y0s$y0_quant2[1,], y0s$org_y0_quant[1,]),
           upper = c(y0s$y0_quant[2,], y0s$y0_quant2[2,], y0s$org_y0_quant[2,]),
           class = rep(c("R1", "R2", "R"), 
                       each = n_tt)) %>%
  mutate(error = y_tt - yhat,
         width = upper-lower,
         cover = ifelse((y_tt > lower & y_tt < upper), 1, 0)) %>%
  group_by(class) %>%
  summarise(rmspe = sqrt(mean(error^2)),
            mape = mean(abs(error)),
            coverage = mean(cover),
            meanwidth = mean(width)) 

# plot
gg_y <- tibble(rbind(coords_tr, coords_tt),
       Truth = c(y_tr, y_tt),
       "R[1]" = c(y_tr, y0s$y0_mean),
       "R[2]" = c(y_tr, y0s$y0_mean2),
       "R == T" = c(y_tr, y0s$org_y0_mean)) %>%
  pivot_longer(-c("easting", "northing"), names_to = "model", values_to = "y") %>%
  ggplot() +
  facet_wrap(~factor(model,
                     levels = c("Truth", "R == T", "R[1]", "R[2]")),
             nrow = 2, labeller = label_parsed) +
  geom_raster(aes(x=easting, y=northing, fill=y)) +
  geom_sf(data = hs_sf, fill = NA) +
  geom_point(data = data.frame(coords_tr, model = "R == T"), 
             aes(easting, northing), size = 0.5, color = "grey15", shape = 4) + 
  geom_point(data = data.frame(coords_S, model = "R[1]"), 
             aes(easting, northing), size = 0.5, color = "grey15", shape = 4) + 
  geom_point(data = data.frame(coords_S2, model = "R[2]"), 
             aes(easting, northing), size = 0.5, color = "grey15", shape = 4) + 
  geom_contour(aes(x = easting, y = northing, z = y),
               color = "black", breaks = seq(-4, 4, by=.5)) +
  scale_fill_scico(palette = "roma", direction = -1) +
  labs(x="", y="", fill = "y") +
  theme(plot.margin = margin(t = -7, l = -10, r = 0, b = -20),
        aspect.ratio = 1/2,
        legend.margin = margin(b = 0, r = 0, t = 0, l = -5))
# for (ext in extension) {
#   ggsave(plot = gg_y,
#          paste0(path, "plots/horseshoe_freeS_y", ext),
#          width = 6, height = 3.6)
# }


