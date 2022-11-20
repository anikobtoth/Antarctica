library(rgdal)
library(tidyverse)

source('./scripts/Helper_Functions.R')

## Load data ####
#antarctica <- readOGR("../Data/Base", "Antarctic_mainland")

v1 <- readRDS("./data/combined_100m_extract.rds")

abiotic <- c("cloud", "wind", "meanTemp", "melt", 
             "elevation", "rugosity", "slope", "totPrecip", "solar", "DDm5")  #don't include ModT, aspect

n <- list.files("../Data/Species/final_results", ".tif$", recursive = FALSE, full.names = FALSE)

# bad_models <- c("Chordata_Aves_Sphenisciformes_Spheniscidae_Pygoscelis_adeliae",
#                 "Chordata_Aves_Procellariiformes___",
#                 "Arthropoda_Entognatha_Poduromorpha___",
#                 "Bryophyta_Bryopsida_Grimmiales___",
#                 "Cyanobacteria_____",
#                 "Tardigrada_____",
#                 "Marchantiophyta_____",
#                 "Ascomycota_Lecanoromycetes_Lecanorales_Lecideaceae__",
#                 "Ascomycota_Lecanoromycetes_Umbilicariales_Umbilicariaceae__")
bad_models <- c("adeliae.tif",
                "Procellariiformes.tif",
                "Poduromorpha.tif",
                "Grimmiales.tif",
                "Cyanobacteria.tif",
                "Tardigrada.tif",
                "Marchantiophyta.tif",
                "Lecideaceae.tif",
                "Umbilicariaceae.tif")

good_models <- n[!n%in% bad_models] 

##### Factor analysis L1 (abiotic variables) ####

cldat <- v1 %>% dplyr::select(all_of(abiotic)) %>% na.omit()
pdat <- data.frame(na.omit(v1 %>% select(pixID, contains("coords"), all_of(abiotic))) %>% 
                     select(pixID, contains("coords")), 
                   consensus = factor_analysis(cldat, mincomp = 0.45, name = "L1"))
 
#combine results
v1 <- full_join(pdat, v1)
v1 <- v1 %>% mutate(x = coords.x1, y = coords.x2)
rownames(v1) <- paste(v1$x, v1$y, sep = "_")

# Classify unclassified L1 pixels 
v1 <- classify_by_neighbours(v1, consensus, res = 100)

###### Hierarchical Factor analysis L2 (SDM data) ############
cldat <- v1 %>% split(.$consensus) %>% purrr::map(~dplyr::select(., all_of(good_models)))
pdat <- map2(cldat, names(cldat), ~factor_analysis(.x, name = paste0("L2_E", .y)))
v1 <- lapply(pdat, cbind) %>% reduce(rbind) %>% data.frame() %>% setNames(c("consensus2")) %>% merge(v1, by = 0, all = TRUE) %>% namerows()
# Classify unclassified L2 pixels 
v1 <- v1 %>% split(.$consensus) %>% purrr::map(classify_by_neighbours, consensus2, maxdist = 3, res = 100) %>% bind_rows()

##### Create output rasters ###############
out <- v1 %>% select(-all_of(good_models), -all_of(abiotic)) %>% 
  mutate(unit_h = paste0("E", consensus, "B", consensus2)
         #unit_d = paste0("env", consensus, "_sdm", V1), 
         #unit_r = paste0("sdm", V1, "_env", consensusA2)
         )

rownames(out) <- paste(out$x, out$y, sep = "_")

## Make rasters #
library(raster)

# unitsV2 <- raster(xmn = -2653500, xmx =2592500, ymn = -2121500, ymx = 2073500, 
#                   crs = CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"), 
#                   nrows = 4195, ncols = 5246)

#unitsV5 <- raster(xmn = -2700100, xmx = 2600100, ymn = -2200100, ymx = 2100100,
#                  crs = CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"),
#                  nrows = 43002, ncols = 53002)

unitsV6 <- raster(xmn = -2661766, xmx = 2614634, ymn = -2490071, ymx = 2321829,
                  crs = CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"),
                  resolution = 100) %>% setValues(NA) %>%
  writeRaster(filename = "../Data/Typology/typ_V6_faHier_10v.tif")

cells <- cellFromXY(unitsV6, cbind(out$x, out$y))
values <- out$unit_h %>% as.factor() %>% as.numeric()

unitsV6 <- update(unitsV6, cell = cells, v = values)

#unitsV6[cellFromXY(unitsV6, cbind(out$x, out$y))] <- out$unit_h %>% as.factor() %>% as.numeric()

#writeRaster(unitsV6, filename = "../Data/Typology/typV6_faHier_10v", 
#            format = "GTiff", overwrite = TRUE)

detach("package:raster", unload = TRUE)

## make rasters for each supergroup
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

### make unit rasters (hierarchical fa) ####
unitname <- unique(out$unit_h)
for(i in seq_along(unitname)){
  unitsV2 <- raster(xmn = -2653500, xmx =2592500, ymn = -2121500, ymx = 2073500, 
                    crs = CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"), 
                    nrows = 4195, ncols = 5246)
  unitsV2[cellFromXY(unitsV2, cbind(out$x, out$y))] <- (out$unit_h == unitname[i])
  writeRaster(unitsV2, filename = paste0("../Data/Typology/unit_rasters_hierFA/", unitname[i]), 
              format = "GTiff", overwrite = TRUE)
}

## make raster key ###
typkey <- out %>% pull(unit_h) %>% unique() %>% sort() %>% 
  data.frame(VALUE = 1:33) %>% setNames(c("unit", "biotic_assemblage")) %>%
  mutate(unit = str_replace(unit, "env", "E") %>% str_replace("_sdm", "B"))
typkey <- rbind(typkey, typkey %>% 
                        mutate(unit = paste0("G1", unit), 
                               biotic_assemblage = 200 + biotic_assemblage), 
                typkey %>% mutate(unit = paste0("G2", unit), 
                                  biotic_assemblage = 100 + biotic_assemblage), 
                data.frame(unit = c("G3E3B3", "E3B8", "L"), biotic_assemblage = c(315, 400, 500)))

#####
## add pixels from new rock outcrop layer #####
pix <- readOGR("../Data/Abiotic_Layers/Points_100m", "centroids_100m_IFA7") %>% st_as_sf()
ovr <- raster::extract(raster("../Data/Typology/typV5_fa_hier_9v.tif"), pix)
 # Level 1
v2 <- readRDS("./data/combined_100m_extract_v7.rds") %>% ## data was created by running prep_data script on the pix above.
  lef_join(data.frame(pix, ovr), by = c("pixID" = "pointid")) %>% filter(is.na(ovr))
fa1 <- readRDS("./results/factor_analyses/fa_L1.rds")
old <- readRDS("./results/factor_analyses/data_L1.rds")

Env <- fapred(v2, abiotic, fa1, old)
v2 <- right_join(Env, v2, by = "pixID")

 # Level 2
fa2 <- list.files("./results/factor_analyses/", pattern = "fa_L2", full.names = T) %>% purrr::map(readRDS)
old2 <- list.files("./results/factor_analyses/", pattern = "data_L2", full.names = T) %>% purrr::map(readRDS)
BioAs <- v2 %>% split(.$typ) %>% 
  pmap_dfr(.l = list(., fa2, old2), 
           .f = ~fapred(..1, good_models, ..2, ..3)) 

v2 <- full_join(BioAs, v2, by = "pixID") %>% mutate(unit = paste0("E", typ.y, "B", typ.x))
l <- unique(v2$unit) %>% sort()
l <- c(l[1:12], "E3B1", "E3B2", l[13], "E3B4", l[14:30])

## make raster 
typ_rock7 <- raster(raster("../Data/Base/rockoutcrop") %>% extent(),
                  crs = CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"),
                  nrow = 48143, ncol = 52775)
typ_rock7[cellFromXY(typ_rock7, cbind(v2$coords.x1, v2$coords.x2))] <- v2$unit %>% factor(levels = l) %>% as.numeric()

writeRaster(typ_rock7, filename = "../Data/Typology/typV5_rock7ADD", 
            format = "GTiff", overwrite = TRUE)


######

