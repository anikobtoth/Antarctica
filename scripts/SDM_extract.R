library(raster)
library(tidyverse)
library(ggfortify)
library(rgl)

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
pix <- xyFromCell(SDMs, rownames(v) %>% as.integer()) %>% data.frame() %>% mutate(cell = rownames(v))

# Produce pseudo-presence table
PA <- lapply(1:34, function(i) scale(v[,i]) >= thresh[i]) %>% reduce(data.frame) %>%
  apply(2, as.numeric) %>% t() %>% data.frame()
rownames(PA) <- n
#PA[-rownames(PA) %in% bad_models,]

pairs <- simpairs_lgnum(PA)
el <- dist2edgelist(pairs, PA)

eln <- el %>% filter(Z.Score < 0.5)
elp <- el %>% filter(Z.Score > 0.5)

gn <- graph_from_data_frame(eln, directed= FALSE)

### SDm and HAbitat pixels join
load("E:/Antarctica/Antarctica/data/spatialUnits.RData")
t <- read_csv("E:/Antarctica/Data/IFA_Hab_SDM_join.txt")
t <- t %>% select(ID, x, y, pointid, POINT_X, POINT_Y, ORIG_FID, Distance, ACBR_ID, ACBR_Name) %>% setNames(c("HabPix", "hP_X", "hP_Y", "SdmPix", "sP_X", "sP_Y", "IFA", "hP_IFA_Dist_m", "ACBR_ID", "ACBR_Name"))

sdmPix_env <- t %>% select(HabPix, SdmPix) %>% 
  merge(envpred_norm, by.x = "HabPix", by.y = "V1") %>% 
  group_by(SdmPix) %>% 
  summarise_all(mean) %>% select(-HabPix) 

sdmPix_weights <- merge(req_var %>% select(V1, rck01_prop), spatialUnits %>% select(HabPix, SdmPix), by.x = "V1", by.y = "HabPix") %>% group_by(SdmPix) %>% summarise(weight = sum(rck01_prop))
vars <- c("elev", "rugos", "precip", "DDminus5", "coast", "wind", "rad")

sdm_hcl <- cmeans(sdmPix_env %>% select(vars), 
                  7, iter.max = 10000, verbose = TRUE, 
                  dist = "euclidean", method = "ufcl", m = 2, 
                  rate.par = 0.3, weights = sdmPix_weights$weight)
  
sdmPix_weights$hcl7W.3 <- sdm_hcl$cluster %>% as.factor()

x <- princomp(sdmPix_env %>% select(vars)) 

autoplot(x, data =sdmPix_weights , col = "hcl7W.3", 
         loadings = TRUE, loadings.label = TRUE, 
         loadings.label.size = 5, alpha = 0.4)

plot3d(x$scores[,1:3], col=sdmPix_weights$hcl8W.1)
text3d(x$loadings[,1:3], texts=rownames(x$loadings), col="red")


#### 
# PCA of SDM pixels based on suitability of functional groups
vars <- c(el$Sp1, el$Sp2) %>% unique()
v1 <- merge(sdmPix_weights, sdmPix %>% filter(cell %in% pix$cell), by.x = "SdmPix", by.y = "id", all = TRUE) %>% merge(v, by.x = "cell", by.y = 0, all =TRUE)
v2 <- merge(sdmPix_env, sdmPix %>% filter(cell %in% pix$cell), by.x = "SdmPix", by.y = "id", all = TRUE) %>% merge(v, by.x = "cell", by.y = 0, all =TRUE)

x <- princomp(scale(v[,vars]))
autoplot(x, data = v1, col = "hcl8W.1", loadings = TRUE, 
         loadings.label = TRUE, loadings.label.size = 3, alpha = 0.4) #+
  #scale_color_manual(values = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "red", "blue", "forestgreen", NA))

cm6 <- cmeans(v[,vars], centers = 6, #v[c(18644, 10258, 36573, 63187, 48772, 54463),vars], 
             iter.max = 10000,verbose = TRUE, 
             dist = "euclidean", method = "ufcl", m = 2, 
             rate.par = 0.3, weights = 1)
d <- data.frame(cl6 = cm6$cluster)

autoplot(x, data = d, col = "cl6" , loadings = TRUE, 
         loadings.label = TRUE, loadings.label.size = 3, alpha = 0.6) #+
  scale_color_manual(values = c(NA, "red"))
  
  plot3d(x$scores[,1:3], col=d$cl6)
  text3d(x$loadings[,1:3], texts=rownames(x$loadings), col="red")
  

i <- 25
autoplot(x, data = v1, col = vars[i] , loadings = TRUE, 
         loadings.label = TRUE, loadings.label.size = 3, alpha = 0.6) + 
  labs(col = "Suitability")+
  ggtitle(vars[i])


#sdmPix_weights <- merge(sdmPix, sdmPix_weights, by.x = "id", by.y = "SdmPix", all = TRUE) 
dat <- merge(merge(sdmPix, sdmPix_weights, by.x = "id", by.y = "SdmPix", all = TRUE), d, by.x = "cell", by.y = 0, all = TRUE) 

test <- dat %>% group_by(hcl8W.1, cl6) %>% summarise(count = length(weight))


ggplot(test %>% filter(!is.na(cl6)), aes(y = count, axis1 = hcl8W.1, axis2 = cl6)) + 
  geom_alluvium(aes(fill = "red"), width = 1/12) + 
  geom_stratum(width = 1/12, fill = "gray", color = "white") + 
  geom_label(stat = "stratum", infer.label= TRUE) + 
  scale_x_discrete(limits = c("habitat", "assemblage"), expand = c(0.05, 0.05)) + 
  scale_fill_brewer(type = "qual", palette = "Set1")+
  ggtitle("Distribution of suitability clusters over habitat clusters")


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



## PCA Plots colored by taxon suitability #####

for(i in 1:length(n)) autoplot(x, data = v, 
                               col = n[i], loadings = TRUE, 
                               loadings.label = TRUE, loadings.label.size = 3, 
                               alpha = 0.3) + labs(col = "Suitability") + ggtitle(n[i])




## ACBR interaction SDMpix 7 habclusters ######################
sppSDMpix <- readRDS("E:/Antarctica/Antarctica/data/Species/Spp_SDMpix_occ.rds")

test <- merge(sdmPix_weights, sdmPix, by.x = "SdmPix", by.y = "id", all.x = TRUE)
test <- merge(test, sppSDMpix, by.x = "SdmPix", by.y = "pointid")
test$unit <- paste(test$ACBR_Name.x, test$hcl7W.3, sep = "_")

sppDat <- occ %>% select(scientific, vernacular, Functional_group, kingdom, phylum, 
                         class, order_, family, genus, species) %>% unique()
t <- merge( sppDat %>% select(scientific, Functional_group),test, all.y = TRUE)

PA1 <- reshape2::dcast(t, Functional_group~unit, fun.aggregate = length, value.var = "SdmPix") %>%
  na.omit() %>% namerows()
PA <- tobinary.single(PA1)
tPA <- t(PA)
pairs <- simpairs(tPA)
el <- dist2edgelist(pairs, tPA)
g <- el %>% dplyr::filter(Z.Score > quantile(Z.Score, 0.2, na.rm = T)) %>% 
  graph_from_data_frame(directed = F)
cl <- cluster_fast_greedy(g)

plot(g, vertex.label = NA, vertex.size = 6, vertex.color = cl$membership)
##################
##############
#############

