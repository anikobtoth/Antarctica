### This code prepares some of the data items needed in other scripts
### Outputs are saved so this script does not need to be run each time.

library(rgdal)
library(raster)
library(tidyverse)

antarctica <- readOGR("../Data/Base", "Antarctic_mainland")

load("./data/DB_habPix.RData")

habPix <- habPix %>% mutate(x = round(as.integer(x), -3), y = round(as.integer(y), -3))

l <- list.dirs("../Data/Species/final_results", recursive = FALSE) %>% 
  file.path("trend_basedist0.tif")
n <- list.dirs("../Data/Species/final_results", recursive = FALSE, full.names = FALSE)


SDMs <- raster::stack(l)
crs(SDMs) <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# ALL small pixels (centroid points)
hp <- readOGR("../Data/Habitats", "habPix_ALL_joinOrig")
hp <- SpatialPointsDataFrame(hp, data = data.frame(hp))

# extract SDM data
hpext <- raster::extract(SDMs, coordinates(hp)[,1:2]) %>% data.frame() %>% setNames(n) %>% scale()
rm(SDMs)

## Format SDM Data #####
hpdat <- cbind(coordinates(hp)[,1:2], hpext) %>% data.frame() %>%
  mutate(x = round(as.integer(coords.x1), -3), y = round(as.integer(coords.x2), -3)) %>% na.omit()

# join habitat data
v1 <- full_join(habPix, hpdat, by = c("x", "y")) %>% select(-genus, -records, -spec, -IFA)

# Add missing lat lons 
habPix_spt <- SpatialPointsDataFrame(coords = cbind(v1$x, v1$y), 
                                     data = v1, proj4string = CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
habPix_wgs <- spTransform(habPix_spt, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs ")
v1 <- v1 %>% mutate(lon = coordinates(habPix_wgs)[,1], lat = coordinates(habPix_wgs)[,2])
rownames(v1) <- paste(v1$x, v1$y, sep = "_")

smallPix <- v1

## 1x1km pixels with any data.
### Most have both SDM and abiotic data but some have only one.

save(smallPix, file = "./data/DB_smallPix_all.Rdata")
