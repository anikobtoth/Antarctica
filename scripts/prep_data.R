
library(rgdal)
library(raster)
library(tidyverse)

## Raster centroids
pix <- readOGR("../Data/Abiotic_Layers/Points_100m", "centroids_100m_IFA")

# environmental data

abiotic <- c("cloud", "wind", "meanTemp", "melt", "modT", "aspect", "elevation", "rock", "rugosity", "slope", "solar", "DDm5", "totPrecip")
l <- list.files("../Data/Abiotic_Layers", ".tif$", full.names = T)
layers <- purrr::map(l, raster) %>% setNames(abiotic)

abiotic_data <- purrr::map(layers, ~raster::extract(.x, coordinates(pix)))
abiotic_data <- cbind(coordinates(pix), reduce(abiotic_data, data.frame) %>% 
                 set_names(abiotic)) %>% tibble() %>% dplyr::select(-rock)

saveRDS(abiotic_data, "./data/abiotic_100m_extract.rds", compress = T)

# SDM data

biotic <- list.dirs("../Data/Species/final_results", recursive = F, full.names = F)
l <- list.files("../Data/Species/final_results", "trend_basedist0.tif", recursive = T, full.names = T)
layers <- raster::stack(l) 
projection(layers) <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m
+no_defs"

biotic_data <- raster::extract(layers, coordinates(pix))
biotic_data <- cbind(coordinates(pix), data.frame(biotic_data) %>% set_names(biotic)) %>% 
  tibble()
saveRDS(biotic_data, "./data/biotic_100m_extract.rds", compress = T)


## combine 
lls <- spTransform(pix, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>% 
  coordinates() %>% data.frame() %>% setNames(c("lon", "lat"))

v1 <- tibble(pixID = data.frame(pix) %>% pull(pointid), 
             lls, 
             full_join(abiotic_data, biotic_data, by = c("coords.x1", "coords.x2")))

saveRDS(v1, "./data/combined_100m_extract.rds", compress =T)
