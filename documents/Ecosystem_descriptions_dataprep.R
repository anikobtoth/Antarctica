library(tidyverse)
library(reshape2)
library(rgdal)
library(sp)
library(cowplot)
library(rts)
library(ggfortify)
library(rgl)

# Helper functions ####
source('./scripts/Helper_Functions.R')

# Load data ####
out <- readRDS("./results/out_TYP_fa_hier_dual_reverse_units.rds")

biotic <- list.dirs("../Data/Species/final_results", recursive = FALSE, full.names = FALSE)

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

good_models <- biotic[!biotic%in% bad_models] 

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

# GBIF occurrence data
GBIF_clean <- readRDS("./data/Species/GBIF_clean_data.rds")

# verbal descriptions
descr <- read_file("./documents/Unit_descriptions.txt")
descr <- descr %>% strsplit("\r\n\r\n") %>% unlist()

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

sppDat <- occurrences %>% dplyr::select(scientific, Functional_group, vernacularName, kingdom, phylum, class, order, family, genus, species) %>% unique()
# Manually remove errors 
sppDat <- sppDat %>% filter(!(scientific == "Chamaesiphon subglobosus" & (!phylum == "Cyanobacteria" | !Functional_group == "Cyanobacteria_____")))
sppDat <- sppDat %>% filter(!(scientific == "Aphanocapsa muscicola" & (!phylum == "Cyanobacteria" | !Functional_group == "Cyanobacteria_____")))
sppDat <- sppDat %>% filter(!(scientific == "Catocarpus gerlachei" & (!class == "Lecanoromycetes" | !Functional_group == "Ascomycota_Lecanoromycetes_Incertae sedis_Rhizocarpaceae__")))
sppDat <- sppDat %>% filter(!(scientific == "Barbilophozia cf. hatcheri" & !family == "Jungermanniaceae"))
sppDat <- sppDat %>% filter(!(scientific == "Geomonhystera antarcticola" & !class == "Adenophorea"))
sppDat <- sppDat %>% filter(!(scientific == "Halteria sp." & !order == "Oligotrichida"))
sppDat <- sppDat %>% filter(!(scientific == "Lepadella patella" & !family == "Lepadellidae"))

PA <- occ %>% data.frame() %>% reshape2::dcast(scientific~fah, fun.aggregate = length) %>% namerows()

# weed out so only one occurrence per species per cell
typology <- raster("../Data/Typology/typV2_fa_hier_12V.tif")
PA_distinct <- occ %>% data.frame() %>% mutate(cell = cellFromXY(typology, occ)) %>% 
  dplyr::select(scientific, cell, fah) %>% distinct() %>%
  reshape2::dcast(scientific~fah, fun.aggregate = length) %>% namerows()

PArel <- apply(PA, 2, function(x) x/sum(x)) 
commonspp <- PA[which(rowSums(PA) >10),] %>% apply(1, which.max)
dom_pct <- 100*(PA[which(rowSums(PA) >10),] %>% apply(1, max)/PA[which(rowSums(PA) >10),] %>% rowSums())

typ <- as(typ_fah, "SpatialPixelsDataFrame")
typ_df <- as.data.frame(typ) %>% filter(!typV2_fa_hier_12v %in% c(6, 11, 17, 23, 29))

rm(out, smallPix)
detach("package:raster", unload = TRUE)

sppLat <- GBIF_clean %>% group_by(species) %>% summarise(minLat = min(decimallatitude), maxLat = max(decimallatitude))
restricted_spp <- sppLat$species[which(sppLat$maxLat < -50)]
restricted_spp <- c(restricted_spp, sppDat$scientific[which(!sppDat$scientific %in% GBIF_clean$species)]) #include species with no GBIF records

# Plots ####

ecosummary <- data %>% group_by(unit_h) %>% summarise(across(all_of(abiotic) , ~ mean(.x, na.rm = TRUE)))
sdmsummary <- data %>% group_by(unit_h) %>% summarise(across(all_of(good_models), ~ mean(.x, na.rm = TRUE)))


singletons <- rownames(PA)[which(rowSums(PA) == 1)]
doubletons <- rownames(PA)[which(rowSums(PA) == 2)]

faunas <- lapply(PA, function(x) rownames(PA)[which(x > 0)])
#sapply(faunas, function(x) length(which(x %in% singletons)))
#sapply(faunas, function(x) length(which(x %in% doubletons)))

commonPA <-PA[rowSums(PA) > 20,]
apply(commonPA, 1, which.max)

## Endemicity
spp_list <- rownames(PA)

ecodat <- melt(data, id.vars = c("ID", "unit_h", "x", "y", "Prop_in_IFA", "lon", 'lat', 'coords.x1', 'coords.x2'))
bio_key <- data.frame(variable = biotic, taxon = c("mites Mesostigmata ", "mites Sarcoptiformes",
                                                               "mites Trombidiformes", "Springtails slim",
                                                               "springtails round", "lichens Acarosporacid",
                                                               "lichens Candelarid ", "lichens, Bacidiacid",
                                                               "lichens Cladonid", "lichens Lecanorid",
                                                               "lichens Lecideacid", "lichens Parmelid",
                                                               "lichens Stereocaulid", "lichens Rhizocarpid",
                                                               "lichens Physcid (shadow)", "lichens Teloschistid",
                                                               "lichens Umbilicarid", "mosses Bryales",
                                                               "mosses Dicranales", "mosses Grimmiales",
                                                               "mosses Hypnales (feather)", "mosses Polytrichales",
                                                               "mosses Pottiales", "algae Green",
                                                               "Petrels", "penguins Adelie", 
                                                               "penguins Chinstrap", "penguins Gentoo",
                                                               "Cyanobacteria", "Liverworts",
                                                               "Nematodes", "Algae", "Rotifers","Tardigrades"
                                                               ))

ecodat <- ecodat %>% left_join(bio_key)

count <- 0
for(i in sort(unique(typ_df$typV2_fa_hier_12v))[25:31]) {
  count <- count + 1
  unitname <- sort(unique(data$unit_h))[i]
  unit <- typ_df %>% dplyr::filter(typV2_fa_hier_12v == i) %>% mutate(ecosystem = as.factor(typV2_fa_hier_12v))
 
  rmarkdown::render('./documents/Ecosystem_Descriptions.Rmd',  
                    output_file =  paste("report_", sort(unique(data$unit_h))[i], '_', Sys.Date(), ".html", sep=''), 
                    output_dir = './documents/reports')
  
  
}

#Render typology overview
rmarkdown::render('./documents/AntarcticTypology_Summary.Rmd',  
                  output_file =  "AntarcticTypology_Overview.html", 
                  output_dir = './documents')


##### Supergroup descriptions 
x <- data %>% dplyr::select(unit_h, all_of(abiotic)) %>% na.omit() 
x %>% dplyr::select(-unit_h) %>% princomp() %>% autoplot(x, col = word(x$unit_h, 1, 1, sep = "_") %>% as.factor() %>% as.numeric())

y <- x %>% dplyr::select(-unit_h) %>% princomp() 
plot3d(y$scores[,1:3], col = word(x$unit_h, 1, 1, sep = "_") %>% as.factor() %>% as.numeric())
plotly::plot_ly(x = y$scores[,1], y = y$scores[,2], z = y$scores[,3], color = word(x$unit_h, 1, 1, sep = "_") %>% as.factor(), size = 0.4) 



ggplot(ecodat %>% dplyr::filter(variable %in% abiotic), 
       aes(x = word(unit_h, 1, 1, sep = "_"), y = value, fill = word(unit_h, 1, 1, sep = "_"))) +
  geom_boxplot() + facet_wrap(~variable) + labs(fill = "Environmental group")

ggplot(ecodat %>% dplyr::filter(variable %in% good_models), 
       aes(x = word(unit_h, 1, 1, sep = "_"), y = value, fill = word(unit_h, 1, 1, sep = "_"))) +
  geom_boxplot() + facet_wrap(~taxon) + labs(fill = "Environmental group")

ecodat %>% filter(grepl("env6", unit_h), !grepl("NA", unit_h), variable %in% good_models) %>% 
  ggplot(aes(x = unit_h, y = value, fill = unit_h)) + geom_boxplot() + facet_wrap(~taxon)
