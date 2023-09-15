library(plyr)
library(tidyverse)
library(cowplot)
library(rts)
library(sf)
library(terra)

# Helper functions ####
source('./scripts/Helper_Functions.R')

# Load data ####
# verbal descriptions
descr <- read_csv("./documents/Unit_descriptionsV6.csv")
descr <- descr %>% filter(!grepl(pattern = "E[1-5]$|G$", LCODE))

# Ecological data
data <- readRDS("./results/ecodatV6_2.rds") %>% filter(y < 2e06)
typ_df <- data %>% dplyr::select(LCODE, typv6_v1pl, x, y)

# Unit raster

typ_fah <- rast("../data/Typology/typv6_v1pl")

#variables
abiotic <- c("cloud", "wind", "meanTemp", "melt", 
             "elevation", "rugosity", "slope", "totPrecip", "solar", "DDm5")  #don't include ModT, aspect

# SDM and environmental data
#biotic <- list.files("../Data/Species/final_results", ".tif$", recursive = FALSE, full.names = FALSE)
good_models <- grep(names(data), pattern = ".tif", value = T)
bad_models <- c("adeliae.tif",
                "Procellariiformes.tif",
                "Poduromorpha.tif",
                "Grimmiales.tif",
                "Cyanobacteria.tif",
                "Tardigrada.tif",
                "Marchantiophyta.tif",
                "Lecideaceae.tif",
                "Umbilicariaceae.tif")

#good_models <- biotic[!biotic%in% bad_models] 

# Antarctica shapefile
antarctica <- st_read("../data/Base", "Antarctic_landpoly") %>% 
  st_simplify(preserveTopology = TRUE, dTolerance = 2000)

#ACBRs
acbrs <- st_read("../data/ACBRs", "ACBRs_v2_2016") %>% st_buffer(dist = 100) %>% filter(!ACBR_ID %in% c(1,2))
acbr_area <- acbrs %>% mutate(area = st_area(geometry)) %>% group_by(ACBR_Name) %>% summarise(area = sum(area)) %>% st_drop_geometry()
acbr_ext <- acbrs %>% split(.$ACBR_Name) %>% 
  purrr::map(~.x %>% st_bbox())


# Occurrence data
occ <- read_csv("../data/Species/Spp_iceFree_occ.csv")
occurrences <- read_csv("../data/Species/Ant_Terr_Bio_Data_FINAL.csv")

occ <- occ %>% dplyr::select(scientific, vernacular, Functional_group, kingdom, phylum, 
                             class, order_, family, genus, species) %>% unique()
occurrences <- occurrences %>% dplyr::select(scientificName, vernacularName, decimalLongitude, decimalLatitude, 
                                             ACBR_ID, ASPA_ID, Date, year, Publish_YEAR, kingdom, phylum, class, 
                                             order, family, genus, species, Snap_IFA, Dist_IFA, 
                                             coordinateUncertaintyInMetres, individualCount)

occurrences <- full_join(occ[,c("scientific", "Functional_group")], occurrences, by = c("scientific" = "scientificName")) %>% unique()
occ <- st_as_sf(occurrences, coords = c("decimalLongitude", "decimalLatitude"), crs = st_crs(4326))

occ <- st_transform(occ, st_crs(3031))
occ$fah <- extract(typ_fah, occ)$typV6_v1pl

sppDat <- occurrences %>% dplyr::select(scientific, Functional_group, vernacularName, kingdom, phylum, class, order, family, genus, species) %>% unique()
# Manually remove errors 
sppDat <- sppDat %>% filter(!(scientific == "Chamaesiphon subglobosus" & (!phylum == "Cyanobacteria" | !Functional_group == "Cyanobacteria_____")))
sppDat <- sppDat %>% filter(!(scientific == "Aphanocapsa muscicola" & (!phylum == "Cyanobacteria" | !Functional_group == "Cyanobacteria_____")))
sppDat <- sppDat %>% filter(!(scientific == "Catocarpus gerlachei" & (!class == "Lecanoromycetes" | !Functional_group == "Ascomycota_Lecanoromycetes_Incertae sedis_Rhizocarpaceae__")))
sppDat <- sppDat %>% filter(!(scientific == "Barbilophozia cf. hatcheri" & !family == "Jungermanniaceae"))
sppDat <- sppDat %>% filter(!(scientific == "Geomonhystera antarcticola" & !class == "Adenophorea"))
sppDat <- sppDat %>% filter(!(scientific == "Halteria sp." & !order == "Oligotrichida"))
sppDat <- sppDat %>% filter(!(scientific == "Lepadella patella" & !family == "Lepadellidae"))

# Raw PA table
PA <- occ %>% data.frame() %>% filter(!is.na(fah)) %>% reshape2::dcast(scientific~fah, fun.aggregate = length) %>% namerows()
PA <- PA[,as.character(sort(unique(data$typv6_v1pl)))] # remove unclassified units

# weed out so only one occurrence per species per cell
PA_distinct <- occ %>% data.frame() %>% mutate(cell = terra::cellFromXY(typ_fah, st_coordinates(occ))) %>% 
  dplyr::select(scientific, cell, fah) %>% distinct() %>%
  reshape2::dcast(scientific~fah, fun.aggregate = length) %>% namerows()
PA_distinct <- PA_distinct[,as.character(sort(unique(data$typv6_v1pl)))] # remove unclassified units

# relative abundance sampled
#PArel <- apply(PA, 2, function(x) x/sum(x))

# PA of common species 
# PA_common <- PA[which(rowSums(PA) > 10),]
# commonspp <- names(PA)[apply(PA_common, 1, which.max)]
# dom_pct <- 100*(apply(PA_common, 1, max)/rowSums(PA_common))

# GBIF occurrence data
GBIF_clean <- readRDS("./data/Species/GBIF_clean_data.rds")
# restricted species
sppLat <- GBIF_clean %>% group_by(species) %>% summarise(minLat = min(decimallatitude), maxLat = max(decimallatitude))
restricted_spp <- sppLat$species[which(sppLat$maxLat < -50)]
restricted_spp <- c(restricted_spp, sppDat$scientific[which(!sppDat$scientific %in% GBIF_clean$species)]) #include species with no GBIF records

# Plots ####

#ecosummary <- data %>% group_by(unit_h) %>% summarise(across(all_of(abiotic) , ~ mean(.x, na.rm = TRUE)))
#sdmsummary <- data %>% group_by(unit_h) %>% summarise(across(all_of(good_models), ~ mean(.x, na.rm = TRUE)))


#singletons <- rownames(PA)[which(rowSums(PA) == 1)]
#doubletons <- rownames(PA)[which(rowSums(PA) == 2)]

#faunas <- lapply(PA, function(x) rownames(PA)[which(x > 0)])
#sapply(faunas, function(x) length(which(x %in% singletons)))
#sapply(faunas, function(x) length(which(x %in% doubletons)))


#apply(PA_common, 1, which.max)

## Endemicity
#spp_list <- rownames(PA)

# Ecological data - long format ####

ecodatAb <- data %>% dplyr::select(-all_of(good_models), -x, -y, -xR, -yR, -aspect) %>% 
  pivot_longer(cols = all_of(abiotic)) %>% na.omit()
ecodatBi <- data %>% dplyr::select(-all_of(abiotic), -x, -y, -xR, -yR, -aspect) %>% 
  #mutate(across(contains(".tif"), ~normalize(.x, method = "range", range = c(1,100)))) %>% # scale as in the factor analysis
  pivot_longer(cols = all_of(good_models)) %>% na.omit()

bio_key <- data.frame(name = good_models, taxon = c("lichens Acarosporacid", "penguins Chinstrap",
                                                    "lichens Bacidiacid", "mosses Bryales", "lichens Candelarid", 
                                                    "algae Green", "lichens Cladonid", "mosses Dicranales",
                                                    "Springtails slim", "mosses Hypnales",
                                                    "lichens Lecanorid", 
                                                    "mites Mesostigmata ", "Nematodes", "Algae", "penguins Gentoo",
                                                    "lichens Parmelid", "lichens Physcid", 
                                                    "mosses Polytrichales", "mosses Pottiales", "lichens Rhizocarpid",
                                                    "Rotifers", "mites Sarcoptiformes", "lichens Stereocaulid", 
                                                    "lichens Teloschistid", "mites Trombidiformes"))

ecodatBi <- ecodatBi %>% left_join(bio_key)

## raw elevations ###
elev_table <- ecodatAb %>% filter(name == "elevation") %>%
  group_by(LCODE) %>% summarise(min = max(0, min(value, na.rm = T)), 
                                max = max(value, na.rm = T), 
                                mean = mean(value, na.rm = T), 
                                median = median(value, na.rm = T),
                                lower_90 = quantile(value, 0.05, na.rm = T), 
                                upper_90 = quantile(value, 0.95, na.rm = T))

#### Antarctica basemap ####
basemap <- ggplot() +  
  geom_sf(data=antarctica, fill="gray50", color= NA, size=0.25) + 
  coord_sf(datum = st_crs(3031)) +
  geom_tile(data = slice(typ_df, sample(1:nrow(typ_df)), 2000000), aes(x = x, y = y), fill = "gray35", col = "gray35", lwd = 0.2)

#### Generate reports ####
units <- data %>% dplyr::select(typv6_v1pl, LCODE) %>% distinct() %>% arrange(LCODE)
typv6 <- units %>% pull(typv6_v1pl)
unitnames <- units %>% pull(LCODE)
count <- 0
for(i in typv6) {
  count <- count + 1
  unitname <- unitnames[count]
  groupname <- word(unitname, 1,1, sep = "B")
  if(unitname %in% c("G1", "G2", "G3")) groupname <- "G"
  if(unitname == "L") groupname <- "E1"
  unit <- typ_df %>% dplyr::filter(typv6_v1pl == i) 
  rmarkdown::render('./documents/Ecosystem_Descriptions_doc.Rmd',  
                    output_file =  paste("report_", unitnames[count], '_', Sys.Date(), ".docx", sep=''), 
                    output_dir = './documents/reports/short_doc')
  
  
}

#Render typology overview
rmarkdown::render('./documents/AntarcticTypology_Overview.Rmd',  
                  output_file =  "AntarcticTypology_Overview.doc", 
                  output_dir = './documents')

