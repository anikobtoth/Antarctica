library(tidyverse)
library(igraph)
library(reshape2)
library(RColorBrewer)
library(vegan)
library(rgdal)
source('./scripts/Helper_Functions.R')

# Units 
out <- readRDS("../results/out_TYP_fa_hier_dual_reverse_units.rds")

n <- list.dirs("../../Data/Species/final_results", recursive = FALSE, full.names = FALSE)

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


# SDM and environmental data
load("./data/DB_smallPix_all.RData")
data <- merge(out %>% dplyr::select(ID, unit_h, x, y, Prop_in_IFA), smallPix)

# Antarctica shapefile
antarctica <- readOGR("../Data/Base", "Antarctic_mainland")

# Occurrence data
occ <- read_csv("./data/Species/Spp_iceFree_occ.csv")
occurrences <- read_csv("./data/Species/Ant_Terr_Bio_Data_FINAL.csv")

occurrences <- occurrences %>% dplyr::select(scientificName, vernacularName, 
                                             decimalLongitude, decimalLatitude, 
                                             ACBR_ID, ASPA_ID, Date, year, Publish_YEAR, 
                                             kingdom, phylum , class, 
                                             order, family, genus, species, Snap_IFA, 
                                             Dist_IFA, coordinateUncertaintyInMetres, 
                                             individualCount)

## Species data ####
sppDat <- occ %>% dplyr::select(scientific, vernacular, Functional_group, kingdom, phylum, 
                                class, order_, family, genus, species) %>% unique()


indepFG <- sppDat %>% filter(!Functional_group %in% good_models) %>% pull(Functional_group) %>% na.omit() %>% unique()
occurrences <- full_join(sppDat[,c("scientific", "Functional_group")], occurrences, by = c("scientific" = "scientificName")) %>% unique()

occ <- SpatialPointsDataFrame(coords = occurrences[,c("decimalLongitude", "decimalLatitude")],
                              data = occurrences, 
                              proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
occ <- spTransform(occ, "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

library(raster)

typ_fah <- raster("../Data/Typology/typV2_fa_hier_12v.tif")
#typ_fad <- raster("../Data/Typology/typV2_fa_dual_12v.tif")
#typ_far <- raster("../Data/Typology/typV2_fa_revhier_12v.tif")

occ$fah <- raster::extract(typ_fah, occ)

PA <- occ %>% data.frame() %>% reshape2::dcast(scientific~fah, fun.aggregate = length) %>% namerows()

PArel <- apply(PA, 2, function(x) x/sum(x)) 


## Choose point occurrences that were not used in SDM factor analysis 
# indepFG <- sppDat %>% filter(!Functional_group %in% good_models) %>% pull(Functional_group) %>% na.omit() %>% unique()
# 
# indepOCC <- occurrences %>% filter(Functional_group %in% indepFG)
# 
# indocc <- SpatialPointsDataFrame(coords = indepOCC[,c("decimalLongitude", "decimalLatitude")],
#                                  data = indepOCC, 
#                                  proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
# 
# indocc <- spTransform(indocc, "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

# indepOCC$fah <- raster::extract(typ_fah, indocc)
# indepOCC$fad <- raster::extract(typ_fad, indocc)
# indepOCC$far <- raster::extract(typ_far, indocc)
# 
# PA_fah <- indepOCC %>% reshape2::dcast(scientific~fah, fun.aggregate = length) %>% namerows()
# PA_fad <- indepOCC %>% reshape2::dcast(scientific~fad, fun.aggregate = length) %>% namerows()
# PA_far <- indepOCC %>% reshape2::dcast(scientific~far, fun.aggregate = length) %>% namerows()


## Correlation tests between units in each candidate typology
cor(PA_fah, PA_fah) %>% as.dist() %>% median()
cor(PA_fad, PA_fad) %>% as.dist() %>% median()
cor(PA_far, PA_far) %>% as.dist() %>% median()

## simpairs beta-diversity test 
simpairs(t(PA_fah)) %>% unlist() %>% na.omit() %>% hist()
simpairs(t(PA_fad)) %>% unlist() %>% na.omit() %>% hist()
simpairs(t(PA_far)) %>% unlist() %>% na.omit() %>% hist()

# zero inflated poisson model for a test case
M3 <- zeroinfl(`1` ~ `2` | ## Predictor for the Poisson process
                 `2`, ## Predictor for the Bernoulli process;
               dist = 'poisson',
               data = PA_fah)

E2 <- resid(M3, type = "pearson")
N  <- nrow(PA_fah)
p  <- length(coef(M3)) + 1  # '+1' is for variance parameter in NB
sum(E2^2) / (N - p)

# Compositional uniqueness of candidate units (using non-SDM species only)

singletons <- rownames(PA_fah)[which(rowSums(PA_fah) == 1)]
doubletons <- rownames(PA_fah)[which(rowSums(PA_fah) == 2)]

faunas <- lapply(PA_fah, function(x) rownames(PA_fah)[which(x > 0)])
sapply(faunas, function(x) length(which(x %in% singletons)))
sapply(faunas, function(x) length(which(x %in% doubletons)))

incid <- apply(PA_far, 1, function(x) x[which(x >0)] %>% sort())
incid <- incid[which(sapply(incid, length) > 9)]
# Evenness of species distribution among candiate units
purrr::map(incid, ~./sum(.)) %>% purrr::map(~. *log(.)) %>% purrr::map_dbl(~ -sum(.)/log(length(.))) %>% mean(na.rm =T)

