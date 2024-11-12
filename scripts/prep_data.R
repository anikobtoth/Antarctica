
library(sp)
library(sf)
library(raster)
library(tidyverse)

source("./scripts/Helper_Functions.R", echo=TRUE)

### in GIS ###
# download raster union of ice-free polygon layers.https://data.aad.gov.au/metadata/AAS_4568_ice-free_rock_outcrop_union

r <- raster("../Data/Base/rock_union1.tif")
pix <- as(r, "SpatialPoints")

## Read in raster centroids ####
#pix1 <- readOGR("../Data/Abiotic_Layers/Points_100m", "centroids_100m_union")

# environmental data ####
abiotic_data <- extract_env_dat(coordinates(pix), datapath = "../Data/Abiotic_Layers")
saveRDS(abiotic_data, "./data/abiotic_100m_extract_union.rds", compress = T)

# Biotic data ####
biotic_data <- extract_biotic_dat(coordinates(pix), gm(), datapath = "../Data/Species/SDM_interpolated")
saveRDS(biotic_data, "./data/biotic_100m_extract_union.rds", compress = T)

## combine ####
lls <- spTransform(pix, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>% 
  coordinates() %>% data.frame() %>% setNames(c("lon", "lat"))

v1 <- tibble(pixID = data.frame(pix) %>% pull(pointid), 
             lls, 
             full_join(abiotic_data, biotic_data, by = c("coords.x1", "coords.x2")))

saveRDS(v1, "./data/combined_100m_extract_union.rds", compress =T)
