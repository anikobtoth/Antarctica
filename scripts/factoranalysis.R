library(rgdal)
library(tidyverse)

source('./scripts/Helper_Functions.R')


## Load data ####
antarctica <- readOGR("../Data/Base", "Antarctic_mainland")
load("./data/DB_smallPix_all.RData")

abiotic <- c("cloud", "sumtemp", "wind", "temp", "melt", 
             "modT_0315", "elev", "rad", "rugos", "slope", "DDminus5", "precip") ## excl. "coast", "geoT"

l <- list.dirs("../Data/Species/final_results", recursive = FALSE) %>% 
  file.path("trend_basedist0.tif")
n <- list.dirs("../Data/Species/final_results", recursive = FALSE, full.names = FALSE)

bad_models <- c("Chordata_Aves_Sphenisciformes_Spheniscidae_Pygoscelis_adeliae",
                "Chordata_Aves_Procellariiformes___",
                "Arthropoda_Entognatha_Poduromorpha___",
                "Bryophyta_Bryopsida_Grimmiales___",
                "Cyanobacteria_____",
                "Tardigrada_____",
                "Marchantiophyta_____",
                "Ascomycota_Lecanoromycetes_Lecanorales_Lecideaceae__",
                "Ascomycota_Lecanoromycetes_Umbilicariales_Umbilicariaceae__")

good_models <- n[!n%in% bad_models] 


v1 <- smallPix

##### Factor analysis L1 (abiotic variables) ####

cldat <- v1 %>% dplyr::select(all_of(abiotic)) %>% na.omit()
pdat <- data.frame(na.omit(v1[,1:6]), consensus = factor_analysis(cldat, mincomp = 0.3))
 
 #plot(antarctica)
 #points(y~x, data = pdat, col = hsv(as.numeric(pdat$consensus)/nfact), pch = 16, cex = 0.1)
 #legend("bottomleft", fill = hsv((1:nfact)/nfact), legend = 1:nfact)

v1 <- full_join(pdat, v1)
rownames(v1) <- paste(v1$x, v1$y, sep = "_")

# Classify unclassified L1 pixels 
v1 <- classify_by_neighbours(v1, consensus)

# pixels with no classified neighbors fall in separate category
v1$consensus <- factor(v1$consensus, levels = 1:(nfact + 1))
v1$consensus[is.na(v1$consensus)] <- (nfact + 1)

##### Dual Factor analysis L1 (SDM data) ###############
cldat <- v1 %>% dplyr::select(all_of(good_models)) %>% na.omit()
v1 <- merge(v1, cbind(factor_analysis(cldat, mincomp = 0.3)), by = 0, all = TRUE) %>% namerows()
v1$V1 <- factor(v1$V1)
v1 <- classify_by_neighbours(v1, V1)

###### Hierarchical Factor analysis L2 (SDM data) ############
cldat <- v1 %>% split(.$consensus) %>% map(~dplyr::select(., all_of(good_models)))
pdat <- map(cldat, factor_analysis)
v1 <- lapply(pdat, cbind) %>% reduce(rbind) %>% data.frame() %>% setNames(c("consensus2")) %>% merge(v1, by = 0, all = TRUE) %>% namerows()
# Classify unclassified L2 pixels 
v1 <- v1 %>% split(.$consensus) %>% map(classify_by_neighbours, consensus2) %>% bind_rows()

###### Reverse Hierarchical Factor analysis L2 (abiotic variables) ############
cldat <- v1 %>% split(.$V1) %>% map(~dplyr::select(., all_of(abiotic)))
pdat <- map(cldat, factor_analysis)
v1 <- lapply(pdat, cbind) %>% reduce(rbind) %>% data.frame() %>% setNames(c("consensusA2")) %>% merge(v1, by = 0, all = TRUE) %>% namerows()
# Classify unclassified L2 pixels 
v1 <- v1 %>% split(.$V1) %>% map(classify_by_neighbours, consensusA2) %>% bind_rows()


out <- v1 %>% select(-all_of(n), -all_of(abiotic), -coast, -geoT) %>% 
  mutate(unit_h = paste0("env", consensus, "_sdm", consensus2),
         unit_d = paste0("env", consensus, "_sdm", V1), 
         unit_r = paste0("sdm", V1, "_env", consensusA2))

rownames(out) <- paste(out$x, out$y, sep = "_")
## Make rasters ####
library(raster)
hP <- raster("../Data/Habitats/habpix_ifaext")

unitsV2 <- raster(xmn = -2653500, xmx =2592500, ymn = -2121500, ymx = 2073500, 
                  crs = CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"), 
                  nrows = 4195, ncols = 5246)
unitsV2[cellFromXY(unitsV2, cbind(out$x, out$y))] <- out$unit %>% as.factor() %>% as.numeric()
#unitsV1 <- setValues(unitsV1, values = units$consensus2, index = units$cell)

writeRaster(unitsV2, filename = "../Data/Typology/typV2_fa_hier_12v", 
            format = "GTiff", overwrite = TRUE)

detach("package:raster", unload = TRUE)



######### OLD STUFF #########
#### Exploratory abiotic factor analyses  ##############
# Varying factor levels

for(i in 1:9){
  nfact <- i+1
  cldat <- v1 %>% dplyr::select(all_of(abiotic)) %>% na.omit() %>% scale()
  
  test <- fa(cldat, nfact, rotate = "varimax" )
  
  sc <- data.frame(test$scores)
  
  pdat <- data.frame(sc, na.omit(v1[,1:6]))
  pdat$consensus <- apply(sc, 1, which.max)
  
  # Plots ##
  pdf(paste0("../Abiotic_Factor_Analyses/Figs_",nfact,"_factors.pdf"))
  plot(antarctica)
  points(y~x, data = pdat, col = hsv(pdat$consensus/nfact), 
         pch = 16, cex = 0.1)
  
  
  habdat <- data.frame(cldat)
  habdat$consensus <- apply(sc, 1, which.max)
  
  #habdat %>% group_by(consensus) %>% summarise_all(list(min, mean, median, max))
  
  hd <- reshape2::melt(habdat, id.vars = "consensus") %>% mutate(consensus = factor(consensus))
  print(ggplot(hd, aes(group = consensus, y = value, fill = consensus)) + 
          geom_boxplot() + facet_wrap(~variable, scales = "free") +
          scale_fill_manual(values =  hsv(1:nfact/nfact)))
  
  dev.off()
  
  # Loadings ##
  l <- data.frame(matrix(as.numeric(loadings(test)), attributes(loadings(test))$dim, dimnames=attributes(loadings(test))$dimnames))
  write.csv(l, paste0("../Abiotic_Factor_Analyses/Loadings_",nfact,"_factors.csv"))
  write.csv(test$Vaccounted %>% data.frame(), file = paste0("../Abiotic_Factor_Analyses/Variance_",nfact,"_factors.csv") )
  # Create raster ##
  rstr <- raster(SDMs[[1]])
  
  rstr[cellFromXY(rstr, cbind(pdat$x, pdat$y))] <- pdat$consensus
  
  writeRaster(rstr, filename = paste0("../Abiotic_Factor_Analyses/FA", nfact), 
              format = "GTiff", overwrite = TRUE)
  
}



##############
#############

