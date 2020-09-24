library(fields)
library(tidyverse)
library(nFactors)
library(psych)
library(rgdal)

antarctica <- readOGR("../Data/Base", "Antarctic_mainland")
load("C:/Users/Aniko/OneDrive - UNSW/Antarctica/Antarctica/data/DB_habPix.RData")

habPix <- habPix %>% mutate(x = round(as.integer(x), -3), y = round(as.integer(y), -3))

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

abiotic <- c("coast", "geoT", "cloud", "sumtemp", "wind", "temp", "melt", 
             "modT_0315", "elev", "rad", "rugos", "slope", "DDminus5", "precip")

library(raster)
SDMs <- raster::stack(l)
crs(SDMs) <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"


# ALL small pixels (centroid points)
hp <- readOGR("../Data/Habitats", "habPix_ALL_joinOrig")
hp <- SpatialPointsDataFrame(hp, data = data.frame(hp))

# extract SDM data
hpext <- raster::extract(SDMs, coordinates(hp)[,1:2]) %>% data.frame() %>% setNames(n) %>% scale()
rm(SDMs)
detach("package:raster", unload = TRUE)

hpdat <- cbind(coordinates(hp)[,1:2], hpext) %>% data.frame() %>%
  mutate(x = round(as.integer(coords.x1), -3), y = round(as.integer(coords.x2), -3)) %>% na.omit()

# join habitat data
v1 <- full_join(habPix, hpdat, by = c("x", "y")) %>% select(-genus, -records, -spec, -IFA)

# Add lat lons 
habPix_spt <- SpatialPointsDataFrame(coords = cbind(v1$x, v1$y), 
                                     data = v1, proj4string = CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
habPix_wgs <- spTransform(habPix_spt, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs ")
v1 <- v1 %>% mutate(lon = coordinates(habPix_wgs)[,1], lat = coordinates(habPix_wgs)[,2])
rownames(v1) <- paste(v1$x, v1$y, sep = "_")

## Factor analysis h1
cldat <- v1 %>% dplyr::select(all_of(abiotic)) %>% na.omit()
dat <- cldat %>% select(-coast, -geoT)
dat <- sapply(dat, function(x) normalize(x, method = "range", range = c(min(x[x > 0])/10, 1-(min(x[x > 0])/10)), margin = 2)) %>% data.frame()
cldat <- sapply(dat, qnorm) %>% scale()

#cv <- lapply(2:11, function(x) fa(cldat, x, rotate = "varimax")$Vaccounted["Cumulative Var",]) %>% sapply(last)

nfact <- 7
fa1 <- fa(cldat, nfact, rotate = "varimax" )
sc <- data.frame(fa1$scores)
pdat <- data.frame(na.omit(v1[,1:6]), sc)

pdat$consensus <- apply(sc, 1, which.max) %>% factor()

plot(antarctica)
points(y~x, data = pdat, col = hsv(as.numeric(pdat$consensus)/nfact), 
       pch = 16, cex = 0.1)

v1 <- full_join(pdat %>% select(-starts_with("MR")), v1)
rownames(v1) <- paste(v1$x, v1$y, sep = "_")

##### Classify unclassified h1 pixels ######

## points without abiotic data ##
## If (1) centroid inside mainland and (2) an immediate neighbour is classified then match to nearest neighbor
## Else -- coastal pixels and islands receive their own h1 category.
z <- length(which(is.na(v1$consensus)))
v0 <- v1 %>% filter(is.na(consensus)) #%>% split(factor(rep_len(1:round(z/3000), length.out = z))) # split to decrease comp load
v2 <- v1 %>% filter(!is.na(consensus))
v2 <- v2 %>% filter(x %in% unique(v0$x, v0$x+1000, v0$x - 1000) & y %in% unique(v0$y, v0$y+1000, v0$y - 1000))

# split into sections to decrease comp load
for(v0m in v0){
  temp <- v2[rdist.earth(v0m[,c("lon", "lat")], v2[,c("lon", "lat")], miles = F) %>% apply(1, which.min),] %>% 
    select(x, y, lon, lat, consensus) %>% 
    mutate(dist = rdist.earth(v0m[,c("lon", "lat")], v2[,c("lon", "lat")], miles = F) %>% apply(1, min)) 
  temp <- data.frame(x = v0m$x, y = v0m$y, consensus = temp$consensus, dist = temp$dist) %>% filter(dist < 1.5)
  
  v1[paste(temp$x, temp$y, sep = "_"),]$consensus <- temp$consensus
}

# pixels off mainland or with no classified neighbors fall in separate category
v1$consensus <- factor(v1$consensus, levels = 1:(nfact + 1))
v1$consensus[is.na(v1$consensus)] <- (nfact + 1)

detach("package:fields", unload = TRUE)
detach("package:maps", unload = TRUE)

###### Hierarchical Factor analysis ############
cldat <- v1 %>% dplyr::select(all_of(good_models), consensus)

# variance analysis
cvs <- cldat %>% na.omit() %>% split(.$consensus) %>% map(function(y) map(2:12, function(x) fa(y %>% select(-consensus), nfactors = x, rotate = "varimax")$loadings %>% as.matrix() %>% apply(2, max)) %>% sapply(min))
nfact <- sapply(cvs, function(x) which(x >= 0.35) %>% last())

hfa <- map2(cldat %>% split(.$consensus), nfact, function(x, y) fa(x %>% select(-consensus,), nfactors = y, rotate = "varimax"))
hfasc <- map(hfa, ~.$scores) %>% map(na.omit)
consensus2 <- map(hfasc, ~apply(., 1, which.max)) %>% map(data.frame) %>% map(setNames, c("consensus2")) %>% bind_rows(.id = "consensus") 
 
out <- merge(v1, consensus2, by = 0, all = TRUE) %>% select(x, y, lon, lat, Prop_in_IFA, consensus.x, consensus2)
rownames(out) <- paste(out$x, out$y, sep = "_")
### Classify unclassified h2 pixels ####

v0 <- out %>% filter(is.na(consensus2))
v2 <- out %>% filter(!is.na(consensus2))
# to decrease computational burden 
v2 <- v2 %>% filter(x %in% unique(v0$x, v0$x+1000, v0$x - 1000) & y %in% unique(v0$y, v0$y+1000, v0$y - 1000))
library(fields)
dists <- rdist.earth(v0[,c("lon", "lat")], v2[,c("lon", "lat")], miles = F)
temp <- v2[dists %>% apply(1, which.min),] %>% mutate(dist = dists %>% apply(1, min)) 
temp <- data.frame(x0 = v0$x, y0 = v0$y, consensus = v0$consensus.x, temp) %>% filter(dist < 1.5 & consensus == consensus.x)

out[paste(temp$x0, temp$y0, sep = "_"),]$consensus2 <- temp$consensus2
out <- out %>% mutate(unit = paste0("h", consensus.x, "_sdm", consensus2))
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

## code to merge pdat in with SDM_extract code. ####
pointoccs <- pdat %>% select(-starts_with("MR")) %>% mutate(unit = paste(ACBR_Name, consensus, sep = "_"))


sppSDMpix <- readRDS("../Antarctica/data/Species/Spp_SDMpix_occ.rds")
occ <- read_csv("./data/Species/Spp_iceFree_occ.csv")

occ <- occ %>% select(scientific, vernacular, Functional_group, kingdom, phylum, 
                      class, order_, family, genus, species) %>% unique()

pointoccs <- merge(pointoccs, sppSDMpix, by.x = "SdmPix", by.y = "pointid", all.y = TRUE)
t <- merge( occ %>% select(scientific, Functional_group),pointoccs, all.y = TRUE)


PA1 <- reshape2::dcast(t, Functional_group~unit, fun.aggregate = length, value.var = "SdmPix") %>%
  na.omit() %>% namerows()
PA <- tobinary.single(PA1)
tPA <- t(PA)
g <- network_analysis(tPA, threshold = 0.995)

consensus <- merge(v1 %>% select(SdmPix, cell, weight, x, y), 
                   pdat %>% select(-starts_with("MR")), all = TRUE) %>%
  mutate(unit = interaction(ACBR_Name, consensus, sep = "_"))
units <- data.frame(unit = V(g[[1]])$name, consensus2 = g[[2]]$membership) %>% 
  merge(consensus, by = "unit", all = TRUE)



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

