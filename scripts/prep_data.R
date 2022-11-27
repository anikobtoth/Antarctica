
library(rgdal)
library(raster)
library(tidyverse)

### in GIS ###
# union of ice-free polygon layers.
# conversion to raster
# raster pixel centroids (raster to point)

## Read in raster centroids ####
pix <- readOGR("../Data/Abiotic_Layers/Points_100m", "centroids_100m_union")

# environmental data ####
abiotic_data <- extract_env_dat(coordinates(pix))
saveRDS(abiotic_data, "./data/abiotic_100m_extract_union.rds", compress = T)

# SDM data ####
biotic_data <- extract_biotic_dat(coordinates(pix), good_models)
saveRDS(biotic_data, "./data/biotic_100m_extract_union.rds", compress = T)

## combine ####
lls <- spTransform(pix, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>% 
  coordinates() %>% data.frame() %>% setNames(c("lon", "lat"))

v1 <- tibble(pixID = data.frame(pix) %>% pull(pointid), 
             lls, 
             full_join(abiotic_data, biotic_data, by = c("coords.x1", "coords.x2")))

saveRDS(v1, "./data/combined_100m_extract_union.rds", compress =T)
