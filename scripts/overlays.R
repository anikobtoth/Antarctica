library(readxl)
library(raster)
library(fasterize)
library(sf)
library(readr)
library(purrr)
library(dplyr)
library(tidyr)

epsg3031 <- CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
antarctica <- st_read("./data/Base", "Coastline_high_res_polygon_v7.1") %>% st_simplify()
icefree <- st_read("../Data/Base", "rocks_Union_Land")
#typology <- raster("../Data/Typology/typv6_nibble") 

# VOLCANOES ####
volc <- st_read("./data/Volcanoes.kml") %>% 
 st_transform(epsg3031)
volc_land <- st_intersection(volc, antarctica %>% filter(surface == "land")) # clip to land
volc_ifland <- st_intersection(volc_land, icefree) # clip to ice-free areas
 # add back three active volcanoes with no ice free polygon
volc <- rbind(volc_ifland %>% select(Name, Description, gid, surface), volc_land %>% filter(Name %in% c("Mt Rittman", "Linden Island", "Paulet Island")))

# load classification
volc_categories <- read_csv("./data/Volcano_classification.csv")
volc <- left_join(volc, volc_categories %>% select(Name, Classification), by = "Name") %>% 
  filter(Classification != "Inactive") %>% arrange(Classification) %>% 
  mutate(class_num = as.numeric(as.factor(Classification)))

st_write(volc, "../Data/Volcanic", "volcano_class", driver = "ESRI Shapefile")

#rasterize in GIS (the below did not work)
#volclass100m <- rasterize(sf = volc, raster = typology, field = "class_num")

# ROOKERIES ####
# download rookeries data: https://www.penguinmap.com/mapppd/
penguins <- read_csv("./data/Rookeries_raw_fulldata.csv")

rookeries <- penguins %>% 
  pivot_wider(id_cols = c('site_name', 'site_id', "longitude_epsg_4326", 'latitude_epsg_4326', 'common_name', "year"), 
              names_from = "count_type", values_from = "penguin_count", values_fn = mean) %>% 
  group_by(site_name, site_id, longitude_epsg_4326, latitude_epsg_4326, common_name) %>%
  summarise(start = min(year), end = max(year), nests = mean(nests, na.rm = T), adults = mean(adults, na.rm = T), 
            chicks = mean(chicks, na.rm = T)) 

# compare numbers of adults, nests, and chicks
ratios <- rookeries %>%
  filter(common_name != "emperor penguin") %>%
  mutate(cn = chicks/nests, an = adults/nests, ac = adults/chicks) %>% 
  group_by(common_name) %>% 
  summarise(cn = mean(adjust_count(cn), na.rm = T), 
            an = mean(adjust_count(an), na.rm = T), 
            ac = mean(adjust_count(ac), na.rm = T)) %>%
  data.frame() %>% namerows() %>% t() 
ratios <- list(ratios[,1], ratios[,2], ratios[,3])

rookeries <- st_as_sf(rookeries, coords = c("longitude_epsg_4326", "latitude_epsg_4326"), crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) %>% 
  st_transform(epsg3031) %>% 
  filter(common_name != "emperor penguin") %>%
  split(.$common_name) %>%
  map2(ratios, ~.x %>% mutate(BP = BP_translate(nests, adults, chicks, ratios = .y))) %>%
  bind_rows() %>% filter(!is.na(BP)) %>% 
  mutate(est_HA =  ((1.0909*BP)+4641.9)/10000, 
            r   =  2*sqrt(est_HA*10000/pi))

st_write(rookeries, "../Data/Rookeries", "penguin_rookeries", driver = "ESRI Shapefile")
## ArcMap steps: 
## some points manually modified in arcMap
## snap to nearest land
## Buffer by r
## Manually modify colonies >80,000 breeding pairs to match guano stains visible in satellite imagery (Google Earth)
## Clip to land 
## rasterize 

### LAKES ####
# download lakes layer
lakes <- st_read("../Data/Lakes", "add_lakes_high_res_polygon_v7.3")
lakes_land <- st_intersection(lakes, antarctica %>% filter(surface == "land"))

st_write(lakes_land, "../Data/Lakes", "lakes_land", driver = "ESRI Shapefile")

##
# Final assembly was done in arcMap: classified layer was superimposed with geothermal, penguin, and finally lakes layer. 

#