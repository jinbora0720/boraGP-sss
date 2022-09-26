rm(list = ls())

# dependencies
library(sf)
library(tidyverse)
theme_set(theme_bw())
library(rnaturalearth)
library(rnaturalearthdata)

# path 
path <- "~/boraGP-sss/"

# plot extensions
extension <- c(".pdf", ".png")

# Kara sea with SSS 
# data source: https://podaac-tools.jpl.nasa.gov/las/UI.vm#panelHeaderHidden=false;differences=false;autoContour=false;globalMin=0;globalMax=41.04291;xCATID=SMAP_JPL_L3_SSS_CAP_8DAY-RUNNINGMEAN_V5;xDSID=SMAP_JPL_L3_SSS_CAP_8DAY-RUNNINGMEAN_V5;varid=smap_sss-SMAP_JPL_L3_SSS_CAP_8DAY-RUNNINGMEAN_V5;imageSize=auto;over=xy;compute=Nonetoken;tlo=25-Jul-2017%2012:00:00;thi=25-Jul-2017%2012:00:00;catid=SMAP_JPL_L3_SSS_CAP_8DAY-RUNNINGMEAN_V5;dsid=SMAP_JPL_L3_SSS_CAP_8DAY-RUNNINGMEAN_V5;varid=smap_sss-SMAP_JPL_L3_SSS_CAP_8DAY-RUNNINGMEAN_V5;avarcount=0;xlo=29.53375;xhi=84.3775;ylo=56.25;yhi=88.59375;operation_id=Plot_2D_XY_zoom;view=xy;ferret_land_type=filled
# call data (Jul 20, 2017 - Jul 31, 2017)
arctic <- readRDS(paste0(path, "data/NovayaZemlya_SSS_Jul2017.RDS")) %>% 
  filter(lat > 65 & long > 29 & long < 88) 

arctic_avg <- arctic %>% 
  group_by(long, lat) %>% 
  summarise(sss = mean(sss), sst = mean(sst)) %>%
  ungroup()

# conversion to utm 
arctic_utm <- st_transform(st_as_sf(arctic_avg, coords = c("long", "lat"),
                                    crs = 4326), 
                           "+proj=utm +zone=40N +datum=WGS84 +units=km") 

# Ob' and Yenisey sampling location 
river <- data.frame(river = c("Ob'", "Yenisey"), 
                    long = c(66.60, 86.48), 
                    lat = c(66.63, 67.43))
river_utm <- st_transform(st_as_sf(river, coords = c("long", "lat"),
                                   crs = 4326), 
                          "+proj=utm +zone=40N +datum=WGS84 +units=km") 

# boundary
world <- ne_countries(scale = "medium", returnclass = "sf")

## land
rest <- st_crop(world, xmin = 29, xmax = 88, 
                ymin = 65, ymax = 81) 
rest_utm <- rest %>% 
  smoothr::densify(max_distance = 1) %>%
  st_transform("+proj=utm +zone=40N +datum=WGS84 +units=km")

# observed data 
sss_min <- 3.369137e-20
sss_max <- 42.34648
g1 <- arctic_utm %>% ggplot() + 
  geom_sf(aes(col = sss)) + 
  geom_sf(data = rest_utm) + 
  geom_sf(data = river_utm) + 
  scale_color_distiller(palette = "RdYlBu", limits = c(sss_min, sss_max)) +
  labs(color = "SSS\n(psu)", x = "Longitude", y = "Latitude") + 
  geom_text(x = 923, y = 7372, label = "Ob'") + 
  geom_text(x = 1650, y = 7730, label = "Yenisey") + 
  geom_text(x = 650, y = 8600, label = "Novaya\nZemlya") + 
  theme(plot.margin = margin(t = 1, l = -15, r = -10, b = 1), 
        legend.margin = margin(b = 0, r = 0, t = 0, l = -2))
for (ext in extension) {
  ggsave(plot = g1, paste0(path, "plots/trueSSS", ext),
         width = 6.4, height = 4)
}

# river discharge 
# https://arcticgreatrivers.org/discharge/
river <- readRDS(paste0(path, "data/RiverDischarge.RDS"))
cols <- c("Ob'" = "#0072B2", "Yenisey" = "#D55E00")
types <- c("Ob'" = 1, "Yenisey" = 2)
g2 <- river %>% 
  mutate(year = format(date, "%Y"), month = format(date, "%m")) %>% 
  filter(year >= 2017 & year <= 2020) %>% 
  group_by(river, month) %>% 
  summarize(avg_discharge = mean(discharge)) %>% 
  ggplot() + 
  geom_line(aes(month, avg_discharge, group = river,
                col = river, linetype = river), size = 1) + 
  labs(x = "Month", y = expression(Discharge (m^3/sec)), 
       col = "River", linetype = "River") + 
  scale_color_manual(values = cols) +
  theme(legend.position = c(0.87,0.8), 
        legend.background = element_rect(fill = "transparent"), 
        plot.margin = margin(t = 1, l = 0, r = 0, b = 1), 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15))
for (ext in extension) {
  ggsave(plot = g2, paste0(path, "plots/riverdischarge", ext),
         width = 6, height = 4)
}

