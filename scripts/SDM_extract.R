library(raster)
library(ggfortify)
library(rgl)
library(parallel)
library(e1071)
library(data.table)
library(igraph)
library(fields)
library(tidyverse)
source('./scripts/Helper_Functions.R')

l <- list.dirs("E:/Antarctica/Data/Species/final_results", recursive = FALSE) %>% 
  file.path("trend_basedist0.tif")
n <- list.dirs("E:/Antarctica/Data/Species/final_results", recursive = FALSE, full.names = FALSE)

SDMs <- raster::stack(l)
crs(SDMs) <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
v <- getValues(SDMs) %>% data.frame() %>% na.omit() %>% setNames(n)
v <- apply(v, 2, function(x) (x-min(x))/(max(x)-min(x)))


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

# get coordinates of pixels
pix <- xyFromCell(SDMs, rownames(v) %>% as.integer()) %>% data.frame(v)

rm(SDMs)
detach("package:raster", unload = TRUE)

### Prepare pixel data ########
req_var = fread("./data/Habitats/req_var.csv")

envpred_norm = fread("./data/Habitats/envpred_norm.csv")

load("E:/Antarctica/Antarctica/data/spatialUnits.RData")
load("E:/Antarctica/Antarctica/data/DB_sdmPix.Rdata")
t <- read_csv("E:/Antarctica/Data/IFA_Hab_SDM_join.txt")
t <- t %>% dplyr::select(ID, x, y, pointid, POINT_X, POINT_Y, ORIG_FID, Distance, ACBR_ID, ACBR_Name) %>% setNames(c("HabPix", "hP_X", "hP_Y", "SdmPix", "sP_X", "sP_Y", "IFA", "hP_IFA_Dist_m", "ACBR_ID", "ACBR_Name"))

sdmPix_env <- t %>% dplyr::select(HabPix, SdmPix) %>% 
  merge(envpred_norm, by.x = "HabPix", by.y = "V1") %>% 
  group_by(SdmPix) %>% 
  summarise_all(mean) %>% dplyr::select(-HabPix) 

sdmPix_weights <- merge(req_var %>% dplyr::select(V1, rck01_prop), spatialUnits %>% dplyr::select(HabPix, SdmPix), by.x = "V1", by.y = "HabPix") %>% group_by(SdmPix) %>% summarise(weight = sum(rck01_prop))
v1 <- merge(sdmPix_weights, sdmPix %>% filter(cell %in% pix$cell), by = "SdmPix", all = TRUE) %>% 
  merge(v, by.x = "cell", by.y = 0, all =TRUE) %>% merge(sdmPix_env, by = "SdmPix", all = TRUE)

abvars <- c("elev", "rugos", "precip", "DDminus5", "coast", "wind")
bivars <- c("Ochrophyta_____", "Ascomycota_Lecanoromycetes_Lecanorales_Bacidiaceae__", "Rotifera_____", "Ascomycota_Lecanoromycetes_Not assigned_Rhizocarpaceae__", "Arthropoda_Arachnida_Mesostigmata___")

vars <- c(abvars, bivars)

### Multiple cmeans consensus analysis #####
cldat <- v1 %>% select(all_of(vars)) %>% na.omit()
centers <- 8
itmx <- 10000

cl <- makeCluster(detectCores()-2)
clusterExport(cl, c("cldat", "sdmPix_weights", "centers", "itmx"))
clusterEvalQ(cl, 
             library(e1071))

sdm_hcl <- parLapply(cl, 1:100, function(x) cmeans(cldat, centers = centers[x], iter.max = itmx, 
                     verbose = FALSE, dist = "euclidean", method = "ufcl", 
                     m = 2, rate.par = 0.3, weights = sdmPix_weights$weight))

stopCluster(cl)

sdmPix_weights <-  data.frame(sdmPix_weights, purrr::map(sdm_hcl, ~as.factor(.$cluster)) %>% reduce(data.frame)) %>%
  setNames(c(names(sdmPix_weights), paste0("hcl", centers, "W", 1:length(sdm_hcl), ".itmx", itmx)))
            
out <- lapply(data.frame(combn(3:ncol(sdmPix_weights), 2)), 
      function(y) table(sdmPix_weights[,y[1]], 
                        sdmPix_weights[,y[2]]) %>% 
        apply(1, function(x) c(which.max(x), max(x)/sum(x))) %>% t()) %>%
  setNames(sapply(data.frame(combn(3:ncol(sdmPix_weights), 2)), 
                  function(y) paste(y[1], y[2], sep = "_")) ) %>% 
  purrr::map(~data.frame(first = 1:nrow(.), .)) %>% bind_rows(.id = "reps") %>%
  separate(reps, into = c("repA", "repB"), sep = "_") %>% 
  mutate(cl1 = paste(repA, first, sep = "_"),
         cl2 = paste(repB, X1, sep = "_")) %>%
  dplyr::select(cl1, cl2, X2)

g <- out %>% graph_from_data_frame(directed = F)
g <- delete.edges(g, E(g)[which(E(g)$X2 < 0.84)])
plot(g, vertex.size = 7)
cl <- cluster_fast_greedy(g, weights = E(g)$X2)
l <- layout_nicely(g)
plot(g, vertex.size = 4, vertex.color = cl$membership, vertex.label = NA, layout = l)

cons <- data.frame(row.names = V(g)$name, consensus = cl$membership)
cons_class <- lapply(25:ncol(sdmPix_weights), function(x) paste(x, sdmPix_weights[,x], sep = "_")) %>%
  lapply(function(x) cons[as.character(x),1]) %>% reduce(data.frame)

temp <- apply(cons_class, 1, table)
sdmPix_weights$consensus <- temp %>% sapply(which.max) %>% names()
sdmPix_weights$agreement <- temp %>% sapply(max)

x <- princomp(cldat %>% select(all_of(vars))) 

autoplot(x, data =sdmPix_weights , col = "consensus", 
         loadings = TRUE, loadings.label = TRUE, 
         loadings.label.size = 5, alpha = 0.4)

plot3d(x$scores[,1:3], col=sdmPix_weights$consensus)
text3d(x$loadings[,1:3], texts=rownames(x$loadings), col="red")


## ACBR interaction SDMpix 7 habclusters ######################
sppSDMpix <- readRDS("E:/Antarctica/Antarctica/data/Species/Spp_SDMpix_occ.rds")
occ <- read_csv("./data/Species/Spp_iceFree_occ.csv")

pointoccs <- merge(sdmPix_weights[,c(1,2, ncol(sdmPix_weights)-1)], sdmPix, by = "SdmPix", all.x = TRUE)
pointoccs <- merge(pointoccs, sppSDMpix, by.x = "SdmPix", by.y = "pointid")
pointoccs$unit <- paste(pointoccs$ACBR_Name.x, pointoccs$consensus, sep = "_")

sppDat <- occ %>% select(scientific, vernacular, Functional_group, kingdom, phylum, 
                         class, order_, family, genus, species) %>% unique()
t <- merge( sppDat %>% select(scientific, Functional_group),pointoccs, all.y = TRUE)

PA1 <- reshape2::dcast(t, Functional_group~unit, fun.aggregate = length, value.var = "SdmPix") %>%
  na.omit() %>% namerows()
PA <- tobinary.single(PA1)
tPA <- t(PA)
g <- network_analysis(tPA, threshold = 0.985)

consensus <- merge(v1 %>% select(SdmPix, cell, weight, x, y, ACBR_Name), 
                   sdmPix_weights %>% select(SdmPix, consensus, agreement), all = TRUE) %>%
  mutate(unit = interaction(ACBR_Name, consensus, sep = "_"))
units <- data.frame(unit = V(g[[1]])$name, consensus2 = g[[2]]$membership) %>% 
  merge(consensus, by = "unit", all = TRUE)

# units without point occurrences receive their own consensus2 category.
units$consensus2[is.na(units$consensus2)] <- units$unit[is.na(units$consensus2)] %>% as.character %>% 
  as.factor() %>% as.numeric() %>% `+`(max(units$consensus2, na.rm = TRUE))

# add lat longs to use rdist.earth properly later #
sdmPix_spt <- SpatialPointsDataFrame(coords = cbind(units$x, units$y), 
                                     data = units, proj4string = CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
sdmPix_wgs <- spTransform(sdmPix_spt, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs ")
units <- data.frame(units, coordinates(sdmPix_wgs))
names(units)[11:12] <- c("lon", "lat")

rownames(units) <- units$cell

### Identify islands and coastal pixels ####

islands <- read_csv("E:/Antarctica/Data/Species/SDMpixels_LandPoly_join.csv") %>% 
  select(pointid, grid_code, POINT_X, POINT_Y, FID_2, gid, surface, Area_sqm, Distance)

# get closest LAND polygon to points on an ice shelf
iceshelf <- read_csv("E:/Antarctica/Data/Species/SDMpixels_iceshelf_landpoly_join.csv")
iceshelf <- iceshelf %>% select(pointid, grid_code, POINT_X, POINT_Y, FID_3, gid_1, surface_1, Area_sqm_1, Distance_1) %>% 
  setNames(c("pointid", "grid_code", "POINT_X", "POINT_Y", "FID_2", "gid", "surface", "Area_sqm", "Distance"))

is <- rbind(islands %>% filter(surface == "land"), iceshelf)
is <- left_join(is, sdmPix, by = c("pointid"="id")) %>% 
  filter(cell %in% units$cell[is.na(units$consensus2)]) # keep cells that have not been classified


# Group cells into 3 different sized islands
smallislands <- is[which(is$Area_sqm < 1000000),]$cell
medislands   <- is[which(is$Area_sqm >= 1000000 & is$Area_sqm < 100000000),]$cell
largeislands <- is[which(is$Area_sqm >= 100000000 & is$Area_sqm < 1000000000000),]$cell

units$unit <- as.character(units$unit)
units[units$cell %in% smallislands,]$unit <- "OffshoreIsland_0to1sqkm"
units[units$cell %in% medislands,]$unit <- "OffshoreIsland_1to100sqkm"
units[units$cell %in% largeislands,]$unit <- "OffshoreIsland_100plussqkm"

units[units$cell %in% smallislands,]$consensus2 <- max(units$consensus2, na.rm = T) + 1
units[units$cell %in% medislands,]$consensus2 <- max(units$consensus2, na.rm = T) + 1
units[units$cell %in% largeislands,]$consensus2 <- max(units$consensus2, na.rm = T) + 1

# Coastal cells whose nearest land is the mainland grouped with habitat of nearest pixel neighbour ####
x1 <- units[units$cell %in% (is[which(is$Area_sqm > 10^12),] %>% pull(cell)),] %>% select(lon, lat, cell)
x2 <- units[!units$cell %in% (is[which(is$Area_sqm > 10^12),] %>% pull(cell))  & 
              !units$unit %in% c("OffshoreIsland_0to1sqkm", 
                                 "OffshoreIsland_1to100sqkm", 
                                 "OffshoreIsland_100plussqkm") &
              !is.na(units$unit)
            ,] %>% 
  select(lon, lat, cell, unit, consensus2)
temp <- x2[rdist.earth(x1, x2, miles = F) %>% 
             apply(1, which.min),] %>% select(unit, consensus2) %>% 
  mutate(dist = rdist.earth(x1, x2, miles = F) %>% apply(1, min), cell = x1$cell)

units[as.character(temp$cell),]$unit <- temp$unit
units[as.character(temp$cell),]$consensus2 <- temp$consensus2

# pixels without nearby consensus-classified units are not classified by neighbor. 
units[as.character(temp$cell),][which(temp$dist > 6),]$consensus2 <- NA
units[as.character(temp$cell),][which(temp$dist > 6),]$unit <- NA


## Consensus using SDM pseduopresences ####
# used to classify pixels with no environmental data and merge units which are too small #
thresh <- sapply(1:34, function(i) scale(v[,i]) %>% range(na.rm = TRUE) %>% quantile(0.85))
PA <- lapply(1:34, function(i) scale(v[,i]) >= thresh[i]) %>% reduce(data.frame) %>%
    apply(2, as.numeric) %>% data.frame() %>% setNames(n)
rownames(PA) <- rownames(v)
PA <- PA %>% select(n[which(!n %in% bad_models)])

PA <- merge(units %>% select(cell, consensus2), PA, by.x = "cell", by.y = 0)
PA1 <- PA[is.na(PA$consensus2),] %>% namerows() %>% select(-consensus2)
PA <- PA %>% na.omit() %>% namerows() %>% group_by(consensus2) %>% 
  summarise_all(sum) %>% select(-consensus2)
PA <- tobinary.single(PA)

## Match unclassified pixels to most similar unit (within 100km and 0.98 similarity) ####
contable <- apply(PA1, 1, function(x) apply(PA, 1, function(y) rbind(x, y) %>% cont_table)) %>% 
  purrr::map(bind_rows, .id = "unit") %>% bind_rows(.id = "cell") %>% select(-Sp1, -Sp2)

contable$score <- FETmP(contable)

neighbors <- contable %>% split(.$cell) %>% purrr::map(~filter(., score > 0.98)) %>% 
  purrr::map(~.$unit[order(.$score, decreasing = TRUE)] %>% as.numeric())
neighbors <- neighbors[sapply(neighbors, length) > 0]

# visual check
i <- 1
plot(y~x, data = sdmPix, pch = 16, cex = 0.1)
points(y~x, data = units %>% filter(consensus2 == 54), col = "green", pch = 3, cex = 1)
points(y~x, data = units %>% filter(cell == names(neighbors)[1]), col = "red", pch = 3)

## Spatially nearest units ###

u <- list()
for(i in 1:length(neighbors)){
  p <- units %>% filter(cell == as.numeric(names(neighbors))[i]) %>% select(lon, lat)
  u[[i]] <- units %>% filter(consensus2 %in% neighbors[[i]]) %>% 
  select(lon, lat, consensus2, unit) %>% split(.$consensus2) %>% 
    sapply(function(z) rdist.earth(p, z, miles = F) %>% min())
  print(u[[i]])
}
u <- setNames(u, names(neighbors))

temp <- data.frame(
  cell = names(u),
  unit = sapply(u, function(x) names(x)[which.min(x)]),
  dist = sapply(u, min)) %>% 
  filter(dist < 100)

units[as.character(temp$cell),]$consensus2 <- temp$unit

# Remaining mainland cells grouped with habitat of nearest pixel neighbour (if within 6km) #
x1 <- units[is.na(units$consensus2),] %>% select(lon, lat, cell)
x2 <- units[!is.na(units$consensus2),] %>% select(lon, lat, cell, unit, consensus2)

temp <- x2[rdist.earth(x1, x2, miles = F) %>% 
             apply(1, which.min),] %>% select(unit, consensus2) %>% 
  mutate(dist = rdist.earth(x1, x2, miles = F) %>% apply(1, min), cell = x1$cell)

units[as.character(temp$cell),]$unit <- temp$unit
units[as.character(temp$cell),]$consensus2 <- temp$consensus2

# pixels without nearby consensus-classified units are not classified by neighbor. 
units[as.character(temp$cell),][which(temp$dist > 6),]$consensus2 <- NA
units[as.character(temp$cell),][which(temp$dist > 6),]$unit <- NA


# #### lump tiny units (decided not to: experts instead?)#####
# 
# thresh <- sapply(1:34, function(i) scale(v[,i]) %>% range(na.rm = TRUE) %>% quantile(0.7))
# PA <- lapply(1:34, function(i) scale(v[,i]) >= thresh[i]) %>% reduce(data.frame) %>%
#   apply(2, as.numeric) %>% data.frame() %>% setNames(n)
# rownames(PA) <- rownames(v)
# PA <- PA %>% select(n[which(!n %in% bad_models)])
# 
# PA <- merge(units %>% select(cell, consensus2), PA, by.x = "cell", by.y = 0)
# PA <- PA %>% na.omit() %>% namerows() %>% group_by(consensus2) %>% 
#   summarise_all(sum) %>% select(-consensus2)
# PA <- tobinary.single(PA)
# 
# ## Match smallest units to most similar unit (within 100km and 0.98 similarity) ####
# contable <- cont_table(PA)
# 
# contable$score <- FETmP(contable)
# 
# g <- contable %>% graph_from_data_frame(directed = F) %>% delete.edges(E(.)[E(.)$score < 0.9])
# V(g)$frequency <- table(units$consensus2)[V(g)$name]
# plot(g, vertex.size = log(V(g)$frequency), vertex.label = NA)





### produce raster and shapefiles ####
unitsV1 <- raster(SDMs[[1]])

units <- units[sort(units$cell) %>% as.character(),] 

unitsV1[cellFromXY(unitsV1, cbind(units$x, units$y))] <- units$consensus2
#unitsV1 <- setValues(unitsV1, values = units$consensus2, index = units$cell)

writeRaster(unitsV1, filename = "../Data/Typology/typV1", 
            format = "ascii", overwrite = TRUE)

##################
######## Summary stats for units ##########

s <- merge(units, v1 %>% dplyr::select(SdmPix, cell, weight, x, y, ACBR_Name, all_of(good_models), all_of(abvars)))
s1 <- s %>% dplyr::select(-SdmPix, -cell, -weight, -x, -y, -ACBR_Name, -unit, -consensus, -agreement, -lon, -lat) %>%
  group_by(consensus2) %>% summarise_all(list(.median = median, .sd = sd), na.rm = T) 

s2 <- s1 %>% pivot_longer(-consensus2, names_to = "variable", values_to = "value") %>% mutate(variable = as.character(variable)) %>% 
  separate(variable, into = c("variable", "statistic"), sep = "\\.")

sx <- s1 %>% filter(!is.na(consensus2)) %>% na.omit() %>% 
    dplyr::select(consensus2, grep("median", names(s1), value = TRUE))
sx <- data.frame(unit = sx$consensus2, scale(sx %>% dplyr::select(-consensus2)))

x <- prcomp(sx %>% dplyr::select(-unit))

autoplot(x, sx, col = "wind_.median", 
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 2, label = TRUE)

lump <- c(15, 48)
sx %>% dplyr::filter(unit %in% lump) %>% dplyr::select(-unit) %>% 
  t() %>% plot(pch = 16, ylab = lump[2], xlab = lump[1], cex = 0.5)

sx %>% filter(unit %in% lump) %>% dplyr::select(-unit) %>% 
  t() %>% data.frame() %>% setNames(c("first", "second")) %>% 
  lm(second~first, data = .) %>% summary()


d <- dist(x$x, method = "euclidean") %>% dist2edgelist(data.frame(row.names = sx$unit, x$x)) %>%
  setNames(c("unit1", "unit2", "dist", "id"))

d$rstderr <- apply(d, 1, function(x){
  lump <- c(x[1], x[2]) %>% as.numeric()
  sx %>% filter(unit %in% lump) %>% dplyr::select(-unit) %>% 
    t() %>% data.frame() %>% setNames(c("first", "second")) %>% 
    lm(I(second-first)~1, data = .) %>% summary() %>% `$`(sigma) %>% return()
})


g <- d %>% filter(dist < 2 | rstderr < 0.25) %>% 
  graph_from_data_frame(directed = F) 


load("E:/Antarctica/Antarctica/documents/typology_key.Rdata")
npix <- unitskey %>% group_by(consensus2) %>% summarise(pixels = sum(pixels))
plot(g, edge.width = 8*(1/(E(g)$dist)), vertex.size = log(npix[V(g)$name,]$pixels)*2)
####### OLD STUFF #######
### Pseudo presences network analysis #####
# Calculate thresholds at top 20% of output suitability range.
thresh <- sapply(1:34, function(i) scale(v[,i]) %>% range(na.rm = TRUE) %>% quantile(0.85))

# number of suitable pixels by group
sapply(1:34, function(i) length(which(scale(v[,i]) >= thresh[i]))) %>% setNames(n)

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


#### # PCA of SDM pixels based on suitability of functional groups #####
vars <- n[which(!n %in% bad_models)]
v1 <- merge(sdmPix_weights, sdmPix %>% filter(cell %in% pix$cell), by.x = "SdmPix", by.y = "id", all = TRUE) %>% merge(v, by.x = "cell", by.y = 0, all =TRUE)
v2 <- merge(sdmPix_env, sdmPix %>% filter(cell %in% pix$cell), by.x = "SdmPix", by.y = "id", all = TRUE) %>% merge(v, by.x = "cell", by.y = 0, all =TRUE)

x <- princomp(scale(v[,vars]))
temp <- x$loadings
matrix(as.numeric(temp), attributes(temp)$dim, dimnames = attributes(temp)$dimnames) %>% 
  data.frame() %>% select(1:3) %>% abs() %>% rowSums() %>% sort()


autoplot(x, data = v1, col = "hcl8W.1", loadings = TRUE, 
         loadings.label = TRUE, loadings.label.size = 3, alpha = 0.4)

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


i <- 1
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


#############

