library(tidyverse)
library(reshape2)
library(rgdal)
library(sp)
library(cowplot)

# Helper functions ####
source('./scripts/Helper_Functions.R')

# Load data ####
out <- readRDS("./results/out_TYP_fa_hier_dual_reverse_units.rds")

n <- list.dirs("../Data/Species/final_results", recursive = FALSE, full.names = FALSE)

abiotic <- c("cloud", "sumtemp", "wind", "temp", "melt", 
             "modT_0315", "elev", "rad", "rugos", "slope", "DDminus5", "precip") ## excl. "coast", "geoT"

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

# Antarctica shapefile
antarctica <- readOGR("../Data/Base", "Antarctic_landpoly")

# Occurrence data
occ <- read_csv("./data/Species/Spp_iceFree_occ.csv")

occurrences <- read_csv("./data/Species/Ant_Terr_Bio_Data_FINAL.csv")

# Unit raster
library(raster)
typ_fah <- raster("../Data/Typology/typV2_fa_hier_12v.tif")


# Format data #####
data <- merge(out %>% dplyr::select(ID, unit_h, x, y, Prop_in_IFA), smallPix)
#data <- data %>% filter(!unit_h %in% c("env1_sdmNA", "env2_sdmNA", "env3_sdmNA", "env4_sdmNA", "env5_sdmNA"))
occ <- occ %>% dplyr::select(scientific, vernacular, Functional_group, kingdom, phylum, 
                             class, order_, family, genus, species) %>% unique()
occurrences <- occurrences %>% dplyr::select(scientificName, vernacularName, decimalLongitude, decimalLatitude, 
                                             ACBR_ID, ASPA_ID, Date, year, Publish_YEAR, kingdom, phylum, class, 
                                             order, family, genus, species, Snap_IFA, Dist_IFA, 
                                             coordinateUncertaintyInMetres, individualCount)

occurrences <- full_join(occ[,c("scientific", "Functional_group")], occurrences, by = c("scientific" = "scientificName")) %>% unique()
occ <- SpatialPointsDataFrame(coords = occurrences[,c("decimalLongitude", "decimalLatitude")],
                              data = occurrences, 
                              proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
occ <- spTransform(occ, "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
occ$fah <- raster::extract(typ_fah, occ)

PA <- occ %>% data.frame() %>% reshape2::dcast(scientific~fah, fun.aggregate = length) %>% namerows()
PArel <- apply(PA, 2, function(x) x/sum(x)) 

typ <- as(typ_fah, "SpatialPixelsDataFrame")
typ_df <- as.data.frame(typ) %>% filter(!typV2_fa_hier_12v %in% c(6, 11, 17, 23, 29))

rm(out, smallPix)
detach("package:raster", unload = TRUE)

# Plots ####

ecosummary <- data %>% group_by(unit_h) %>% summarise(across(all_of(abiotic) , ~ mean(.x, na.rm = TRUE)))
sdmsummary <- data %>% group_by(unit_h) %>% summarise(across(all_of(good_models), ~ mean(.x, na.rm = TRUE)))
count <- 0

singletons <- rownames(PA)[which(rowSums(PA) == 1)]
doubletons <- rownames(PA)[which(rowSums(PA) == 2)]

faunas <- lapply(PA, function(x) rownames(PA)[which(x > 0)])
#sapply(faunas, function(x) length(which(x %in% singletons)))
#sapply(faunas, function(x) length(which(x %in% doubletons)))

commonPA <-PA[rowSums(PA) > 20,]
apply(commonPA, 1, which.max)

## Endemicity
spp_list <- rownames(PA)
#skip <- download_GBIF_all_species(spp_list, path = "../Data/GBIF/")

GBIF.ALL <- list.files("../Data/GBIF/", pattern = ".RData") %>%   
  
  ## Restrict them to just the strings that match the species list for each run
  subset(. %in% paste(spp_list, "_GBIF_records.RData", sep = "")) %>%
  
  ## Pipe the list into lapply
  lapply(load.GBIF, GBIF_path = "../Data/GBIF/") %>% bind_rows 


for(i in sort(unique(typ_df$typV2_fa_hier_12v))[1:10]) {
  count <- count + 1
  unitname <- sort(unique(data$unit_h))[i]
  unit <- typ_df %>% dplyr::filter(typV2_fa_hier_12v == i) %>% mutate(ecosystem = as.factor(typV2_fa_hier_12v))
  #ecodat <- data %>% dplyr::filter(unit_h == unitname)
  ecodat <- melt(data, id.vars = c("ID", "unit_h", "x", "y", "Prop_in_IFA", "lon", 'lat', 'coords.x1', 'coords.x2'))
  
  rmarkdown::render('./documents/Ecosystem_Descriptions.Rmd',  
                    output_file =  paste("report_", sort(unique(data$unit_h))[i], '_', Sys.Date(), ".html", sep=''), 
                    output_dir = './documents/reports')
  
  
}



