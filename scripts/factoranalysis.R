library(rgdal)
library(tidyverse)

source('./scripts/Helper_Functions.R')

## Load data ####
#antarctica <- readOGR("../Data/Base", "Antarctic_mainland")

v1 <- readRDS("./data/combined_100m_extract_union.rds")

abiotic <- c("cloud", "wind", "meanTemp", "melt", 
             "elevation", "rugosity", "slope", "totPrecip", "solar", "DDm5")  #don't include ModT, aspect

good_models <- gm()
  
##### Factor analysis L1 (abiotic variables) ####

cldat <- v1 %>% dplyr::select(all_of(abiotic)) %>% na.omit()
pdat <- data.frame(na.omit(v1 %>% dplyr::select(pixID, contains("coords"), all_of(abiotic))) %>% 
                     dplyr::select(pixID, contains("coords")), 
                     factor_analysis(cldat, name = "L1", scale = 1, nfact = 5))
 
#combine results
v1 <- full_join(pdat, v1)
v1 <- v1 %>% mutate(x = coords.x1, y = coords.x2)
rownames(v1) <- paste(v1$x, v1$y, sep = "_")

###### Hierarchical Factor analysis L2 (SDM data) ############
cldat <- v1 %>% split(.$consensus) %>% purrr::map(~dplyr::select(., all_of(good_models)))
pdat <- map2(1:5, c(4,6,7,5,6), ~factor_analysis(cldat[[.x]], name = paste0("L2_E", .x), scale = TRUE, nfact = .y))
v1 <- lapply(pdat, cbind) %>% reduce(rbind) %>% data.frame() %>% rownames_to_column() %>% tibble() %>% 
  setNames(c("pix", "consensus2", "confidence2")) %>% full_join(v1 %>% rownames_to_column("pix"), by = "pix")

##### Create output rasters ###############
out <- v1 %>% dplyr::select(-all_of(good_models), -all_of(abiotic), -modT, -aspect) %>% 
  mutate(unit_h = paste0("E", consensus, "B", consensus2))

#rownames(out) <- paste(out$x, out$y, sep = "_")

## Make rasters: Full Tier 2 typology #
library(raster)

unitsV7 <- raster(xmn = -2661867, xmx = 2614735, ymn = -2490172, ymx = 2321930,
                  crs = CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"),
                  resolution = 100) %>% setValues(NA) %>%
  writeRaster(filename = "../Data/Typology/typ_V7.tif", overwrite = TRUE)

out$rastervalue <- out$unit_h %>% as.factor() %>% as.numeric()
out <- out %>% filter(rastervalue != 27)
cells <- cellFromXY(unitsV7, cbind(out$x, out$y))
values <- out$unit_h %>% as.factor() %>% as.numeric()
unitsV7[cells] <- values
writeRaster(unitsV7, filename = "../Data/Typology/typV6_raw_10v", 
            format = "GTiff", overwrite = TRUE)

#in arcMap: apply Nibble function

## make rasters for each major environmental unit (Tier1)
for(i in 1:5){
  print(i)
  outi <- out %>% filter(consensus == i)
  units <- raster(xmn = -2700100, xmx = 2600100, ymn = -2200100, ymx = 2100100,
                    crs = CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"),
                    nrows = 43002, ncols = 53002)
  units[cellFromXY(units, cbind(outi$x, outi$y))] <- outi$consensus2
  
  writeRaster(units, filename = paste0("../Data/Typology/typV5_fa_hier_Env", i), 
              format = "GTiff", overwrite = TRUE)
  
}

# Tier 1 confidence raster

confT1 <- r
cells <- cellFromXY(confT1, cbind(out$x, out$y))
confT1[cells] <- out$confidence
writeRaster(confT1, filename = "../Data/Typology/ConfT1", 
            format = "GTiff", overwrite = TRUE)

# Tier 2 confidence raster

confT2 <- r
cells <- cellFromXY(confT2, cbind(out$x, out$y))
confT2[cells] <- out$confidence2
writeRaster(confT2, filename = "../Data/Typology/ConfT2", 
            format = "GTiff", overwrite = TRUE)


## make raster key ###
typkey <- out %>% filter(!grepl("NA", unit_h)) %>% pull(unit_h) %>% unique() %>% sort() %>% 
  data.frame(VALUE = 1:28) %>% setNames(c("unit", "habitat")) %>%
  mutate(unit = str_replace(unit, "env", "E") %>% str_replace("_sdm", "B"))
typkey <- rbind(typkey, typkey %>% 
                        mutate(unit = paste0("G1", unit), 
                               habitat = 200 + habitat), 
                typkey %>% mutate(unit = paste0("G2", unit), 
                                  habitat = 100 + habitat), 
                data.frame(unit = c("G3E3B3", "E3B8", "L"), habitat = c(315, 400, 500)))

