library(raster)
library(tidyverse)


l <- list.dirs("E:/Antarctica/Data/Species/final_results", recursive = FALSE) %>% 
  file.path("trend_basedist0.tif")
n <- list.dirs("E:/Antarctica/Data/Species/final_results", recursive = FALSE, full.names = FALSE)

SDMs <- stack(l)
crs(SDMs) <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
v <- getValues(SDMs) %>% data.frame() %>% na.omit()
names(v) <- n

bad_models <- c("Chordata_Aves_Sphenisciformes_Spheniscidae_Pygoscelis_adeliae",
                "Chordata_Aves_Procellariiformes___",
                "Arthropoda_Entognatha_Poduromorpha___",
                "Bryophyta_Bryopsida_Grimmiales___",
                "Cyanobacteria_____",
                "Tardigrada_____",
                "Marchantiophyta_____",
                "Ascomycota_Lecanoromycetes_Lecanorales_Lecideaceae__",
                "Ascomycota_Lecanoromycetes_Umbilicariales_Umbilicariaceae__")
                
                

# Calculate thresholds at top 20% of output suitability range.
thresh <- sapply(1:34, function(i) scale(v[,i]) %>% range(na.rm = TRUE) %>% quantile(0.85))

# number of suitable pixels by group
sapply(1:34, function(i) length(which(scale(v[,i]) >= thresh[i]))) %>% setNames(n)

# get coordinates of pixels
pix <- xyFromCell(SDMs, rownames(v) %>% as.integer())

# Produce pseudo-presence table
PA <- lapply(1:34, function(i) scale(v[,i]) >= thresh[i]) %>% reduce(data.frame) %>%
  apply(2, as.numeric) %>% t() %>% data.frame()
rownames(PA) <- n

pairs <- simpairs_lgnum(PA)
el <- dist2edgelist(pairs, PA)

eln <- el %>% filter(Z.Score < 0.5)
elp <- el %>% filter(Z.Score > 0.5)

gp <- graph_from_data_frame(elp, directed= FALSE)

### SDm and HAbitat pixels join

t <- read_csv("E:/Antarctica/Data/IFA_Hab_SDM_join.txt")
t <- t %>% select(ID, x, y, pointid, POINT_X, POINT_Y, ORIG_FID, Distance, ACBR_ID, ACBR_Name) %>% setNames(c("HabPix", "hP_X", "hP_Y", "SdmPix", "sP_X", "sP_Y", "IFA", "hP_IFA_Dist_m", "ACBR_ID", "ACBR_Name"))

sdmPix_env <- t %>% select(HabPix, SdmPix) %>% 
  merge(envpred_norm, by.x = "HabPix", by.y = "V1") %>% 
  group_by(SdmPix) %>% 
  summarise_all(mean) %>% select(-HabPix) 

sdmPix_weights <- merge(req_var %>% select(V1, rck01_prop), spatialUnits %>% select(HabPix, SdmPix), by.x = "V1", by.y = "HabPix") %>% group_by(SdmPix) %>% summarise(weight = sum(rck01_prop))
vars <- c("elev", "rugos", "precip", "DDminus5", "coast", "wind", "rad")

sdm_hcl <- cmeans(sdmPix_env %>% select(SdmPix, vars), 
                  8, iter.max = 10000, verbose = TRUE, 
                  dist = "euclidean", method = "ufcl", m = 2, 
                  rate.par = 0.3, weights = sdmPix_weights$weight)
  
sdmPix_weights$hcl8W.1 <- sdm_hcl$cluster %>% as.factor()

x <- princomp(sdmPix_env %>% select(vars)) 

autoplot(x, data =sdmPix_weights , col = "hcl8W.1", 
         loadings = TRUE, loadings.label = TRUE, 
         loadings.label.size = 5, alpha = 0.4)

plot3d(x$scores[,1:3], col=sdmPix_weights$hcl6.1)
text3d(x$loadings[,1:3], texts=rownames(x$loadings), col="red")


#### 
# PCA of SDM pixels based on suitability of functional groups
vars <- c(eln$Sp1, eln$Sp2) %>% unique()
v1 <- merge(sdmPix_weights, sdmPix, by.x = "SdmPix", by.y = "id", all = TRUE) %>% merge(v, by.x = "cell", by.y = 0, all =TRUE)

x <- princomp(scale(v[,vars]))
autoplot(x, data = v1, col = "hcl8W.1", loadings = TRUE, 
         loadings.label = TRUE, loadings.label.size = 3, alpha = 0.6)

cm <- cmeans(v[,vars], centers = v[c(18644, 10258, 36573, 63187, 48772, 54463),vars], iter.max = 1000,verbose = TRUE, 
             dist = "euclidean", method = "ufcl", m = 2, 
             rate.par = 0.3, weights = 1)
d <- data.frame(cl6 = cm$cluster)

autoplot(x, data = d, col = "Low_Chord" , loadings = TRUE, 
         loadings.label = TRUE, loadings.label.size = 3, alpha = 0.6)


##### Experimental groupings ######  
rownames(d) <- rownames(v)
d$Low_Ochr <- v$Ochrophyta_____ < mean(range(v$Ochrophyta_____)) & 
  v$Chordata_Aves_Sphenisciformes_Spheniscidae_Pygoscelis_papua > mean(range(v$Chordata_Aves_Sphenisciformes_Spheniscidae_Pygoscelis_papua))
d$High_Moss <- v$Bryophyta_Bryopsida_Bryales___ > -6 &
  v$Bryophyta_Bryopsida_Hypnales___ > -6 &
  v$Bryophyta_Bryopsida_Polytrichales___ > -6
d$High_Ochr_Tar_Cyan <- v$Ochrophyta_____ > -7 &
  v$Tardigrada_____ > -7 &
  v$Cyanobacteria_____ > -7

# Areas that birds can't survive in
d$Low_Chord  <- v$Chordata_Aves_Procellariiformes___ < mean(range(v$Chordata_Aves_Procellariiformes___)) &
                 v$Chordata_Aves_Sphenisciformes_Spheniscidae_Pygoscelis_adeliae < mean(range(v$Chordata_Aves_Sphenisciformes_Spheniscidae_Pygoscelis_adeliae)) &
                  v$Chordata_Aves_Sphenisciformes_Spheniscidae_Pygoscelis_antarctica < mean(range(v$Chordata_Aves_Sphenisciformes_Spheniscidae_Pygoscelis_antarctica)) &
                   v$Chordata_Aves_Sphenisciformes_Spheniscidae_Pygoscelis_papua < mean(range(v$Chordata_Aves_Sphenisciformes_Spheniscidae_Pygoscelis_papua))
sdmPix_sub <- merge(sdmPix_env, sdmPix, by.x = "SdmPix", by.y = "id")

a <- merge(sdmPix_sub, d, by.x = "cell", by.y = 0, all = TRUE)
a <- pivot_longer(a, cols= c("coast", "geoT", "cloud", "sumtemp", "wind", "temp", "melt", "modT_0315", "elev", "rad", "rugos", "slope", "DDminus5", "precip"))
ggplot(a, aes(x = Low_Ochr, y = value)) + geom_boxplot() + facet_wrap(~name, ncol = 7, nrow = 2)



## PCA Plots colored by taxon suitability 

for(i in 1:length(n)) autoplot(x, data = v, 
                               col = n[i], loadings = TRUE, 
                               loadings.label = TRUE, loadings.label.size = 3, 
                               alpha = 0.3) + labs(col = "Suitability") + ggtitle(n[i])
