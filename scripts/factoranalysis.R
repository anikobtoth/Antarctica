library(rgdal)
library(tidyverse)

source('./scripts/Helper_Functions.R')

## Load data ####
#antarctica <- readOGR("../Data/Base", "Antarctic_mainland")

v1 <- readRDS("./data/combined_100m_extract.rds")

abiotic <- c("cloud", "wind", "meanTemp", "melt", 
             "elevation", "rugosity", "slope", "totPrecip", "solar")  #don't include ModT, aspect, DDm5

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
                   consensus = factor_analysis(cldat, mincomp = 0.3))
 
 #plot(antarctica)
 #points(y~x, data = pdat, col = hsv(as.numeric(pdat$consensus)/nfact), pch = 16, cex = 0.1)
 #legend("bottomleft", fill = hsv((1:nfact)/nfact), legend = 1:nfact)

v1 <- full_join(pdat, v1)
v1 <- v1 %>% mutate(x = coords.x1, y = coords.x2)
rownames(v1) <- paste(v1$x, v1$y, sep = "_")

# Classify unclassified L1 pixels 
v1 <- classify_by_neighbours(v1, consensus, res = 100)

# pixels with no classified neighbors fall in separate category
 #v1$consensus <- factor(v1$consensus, levels = 1:(nfact + 1))
 #v1$consensus[is.na(v1$consensus)] <- (nfact + 1)

##### Dual Factor analysis L1 (SDM data) ###############
cldat <- v1 %>% dplyr::select(all_of(good_models)) %>% na.omit()
v1 <- merge(v1, cbind(factor_analysis(cldat, mincomp = 0.3)), by = 0, all = TRUE) %>% namerows()
v1$V1 <- factor(v1$V1)
v1 <- classify_by_neighbours(v1, V1)

###### Hierarchical Factor analysis L2 (SDM data) ############
cldat <- v1 %>% split(.$consensus) %>% purrr::map(~dplyr::select(., all_of(good_models)))
pdat <- purrr::map(cldat, factor_analysis)
v1 <- lapply(pdat, cbind) %>% reduce(rbind) %>% data.frame() %>% setNames(c("consensus2")) %>% merge(v1, by = 0, all = TRUE) %>% namerows()
# Classify unclassified L2 pixels 
v1 <- v1 %>% split(.$consensus) %>% purrr::map(classify_by_neighbours, consensus2) %>% bind_rows()

###### Reverse Hierarchical Factor analysis L2 (abiotic variables) ############
cldat <- v1 %>% split(.$V1) %>% purrr::map(~dplyr::select(., all_of(abiotic)))
pdat <- purrr::map(cldat, factor_analysis)
v1 <- lapply(pdat, cbind) %>% reduce(rbind) %>% data.frame() %>% setNames(c("consensusA2")) %>% merge(v1, by = 0, all = TRUE) %>% namerows()
# Classify unclassified L2 pixels 
v1 <- v1 %>% split(.$V1) %>% purrr::map(classify_by_neighbours, consensusA2) %>% bind_rows()


##### Create output rasters ###############
out <- v1 %>% select(-all_of(n), -all_of(abiotic), -ModT, -aspect, -DDm5) %>% 
  mutate(unit_h = paste0("env", consensus, "_sdm", consensus2)
         #unit_d = paste0("env", consensus, "_sdm", V1), 
         #unit_r = paste0("sdm", V1, "_env", consensusA2)
         )

rownames(out) <- paste(out$x, out$y, sep = "_")

## Make rasters #
library(raster)

# unitsV2 <- raster(xmn = -2653500, xmx =2592500, ymn = -2121500, ymx = 2073500, 
#                   crs = CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"), 
#                   nrows = 4195, ncols = 5246)

unitsV4 <- raster(xmn = -2700100, xmx = 2600100, ymn = -2200100, ymx = 2100100,
                  crs = CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"),
                  nrows = 43002, ncols = 53002)
unitsV4[cellFromXY(unitsV4, cbind(out$x, out$y))] <- out$unit_h %>% as.factor() %>% as.numeric()

writeRaster(unitsV4, filename = "../Data/Typology/typV5_fa_hier_9v", 
            format = "GTiff", overwrite = TRUE)

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

####  Exploratory - canditate units
rownames(units) <- paste(units$x, units$y, sep = "_")
habdat <- merge(units, v1, by = 0, all = TRUE)
hd <- habdat %>% dplyr::select(unit_h, contains("Aves")) %>% 
  reshape2::melt(id.vars = "unit_h") %>% mutate(unit_h = factor(unit_h))

ggplot(hd, aes(x = unit_h, y = value)) + geom_boxplot() + facet_wrap(~variable, scales = "free")

hd <- hd %>% separate(unit_h, sep = "_", into= c("abiotic", "biotic"))
ggplot(hd, aes(x = abiotic, col = abiotic, y = value)) + geom_boxplot(outlier.size = .2) + facet_wrap(~variable, scales = "free")


#####

######

