library(tidyverse)
#library(reshape2)
library(rgdal)
library(sp)
library(cowplot)
library(rts)
library(ggfortify)
library(rgl)
library(rgeos)
library(sf)

# Helper functions ####
source('./scripts/Helper_Functions.R')

# Load data ####
out <- readRDS("./results/out_TYP_hier_units_V5.rds")

abiotic <- c("cloud", "wind", "meanTemp", "melt", 
             "elevation", "rugosity", "slope", "totPrecip", "solar")  #don't include ModT, aspect, DDm5

# SDM and environmental data
smallPix <- readRDS("./data/combined_100m_extract.rds")
biotic <- smallPix %>% dplyr::select(ends_with(".tif")) %>% names()
bad_models <- c("adeliae.tif",
                "Procellariiformes.tif",
                "Poduromorpha.tif",
                "Grimmiales.tif",
                "Cyanobacteria.tif",
                "Tardigrada.tif",
                "Marchantiophyta.tif",
                "Lecideaceae.tif",
                "Umbilicariaceae.tif")

good_models <- biotic[!biotic%in% bad_models] 

# Antarctica shapefile
antarctica <- st_read("../Data/Base", "Antarctic_landpoly") %>% 
  st_simplify(preserveTopology = TRUE, dTolerance = 2)

# Occurrence data
occ <- read_csv("./data/Species/Spp_iceFree_occ.csv")

occurrences <- read_csv("./data/Species/Ant_Terr_Bio_Data_FINAL.csv")

# Unit raster
library(raster)
typ_fah <- raster("../Data/Typology/typV5_fa_hier_9v.tif")

# GBIF occurrence data
GBIF_clean <- readRDS("./data/Species/GBIF_clean_data.rds")

# verbal descriptions
descr <- read_file("./documents/Unit_descriptionsV5.txt")
descr <- descr %>% strsplit("\r\n\r\n") %>% unlist()

# Format data #####
typ_df <- out %>% mutate(typV5 = as.factor(unit_h) %>% as.numeric()) %>% 
  filter(!is.na(consensus2)) %>% 
  dplyr::select(x, y, unit_h, typV5) %>% tibble()

data <- full_join(out %>% dplyr::select(pixID, unit_h, x, y), smallPix, by = "pixID")

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

# Raw PA table
PA <- occ %>% data.frame() %>% filter(!is.na(fah)) %>% reshape2::dcast(scientific~fah, fun.aggregate = length) %>% namerows()
PA <- PA[,as.character(sort(unique(typ_df$typV5)))] # remove unclassified units

# weed out so only one occurrence per species per cell
PA_distinct <- occ %>% data.frame() %>% mutate(cell = cellFromXY(typ_fah, occ)) %>% 
  dplyr::select(scientific, cell, fah) %>% distinct() %>%
  reshape2::dcast(scientific~fah, fun.aggregate = length) %>% namerows()
PA_distinct <- PA_distinct[,as.character(sort(unique(typ_df$typV5)))] # remove unclassified units

# relative abundance sampled
PArel <- apply(PA, 2, function(x) x/sum(x))

# PA of common species 
PA_common <- PA[which(rowSums(PA) > 10),]
commonspp <- names(PA)[apply(PA_common, 1, which.max)]
dom_pct <- 100*(apply(PA_common, 1, max)/rowSums(PA_common))

rm(out, smallPix)

sppLat <- GBIF_clean %>% group_by(species) %>% summarise(minLat = min(decimallatitude), maxLat = max(decimallatitude))
restricted_spp <- sppLat$species[which(sppLat$maxLat < -50)]
restricted_spp <- c(restricted_spp, sppDat$scientific[which(!sppDat$scientific %in% GBIF_clean$species)]) #include species with no GBIF records

# Plots ####

ecosummary <- data %>% group_by(unit_h) %>% summarise(across(all_of(abiotic) , ~ mean(.x, na.rm = TRUE)))
sdmsummary <- data %>% group_by(unit_h) %>% summarise(across(all_of(good_models), ~ mean(.x, na.rm = TRUE)))


singletons <- rownames(PA)[which(rowSums(PA) == 1)]
doubletons <- rownames(PA)[which(rowSums(PA) == 2)]

faunas <- lapply(PA, function(x) rownames(PA)[which(x > 0)])
sapply(faunas, function(x) length(which(x %in% singletons)))
sapply(faunas, function(x) length(which(x %in% doubletons)))


apply(PA_common, 1, which.max)

## Endemicity
spp_list <- rownames(PA)

ecodat <- data %>% dplyr::select(-all_of(bad_models)) %>% 
  pivot_longer(cols = all_of(c(abiotic, good_models)))
bio_key <- data.frame(name = biotic, taxon = c("lichens Acarosporacid", "penguins Adelie", "penguins Chinstrap",
                                               "lichens Bacidiacid", "mosses Bryales", "lichens Candelarid ", 
                                               "algae Green", "lichens Cladonid", "Cyanobacteria","mosses Dicranales",
                                               "Springtails slim", "mosses Grimmiales", "mosses Hypnales (feather)",
                                               "lichens Lecanorid", "lichens Lecideacid",  "Liverworts",
                                               "mites Mesostigmata ", "Nematodes", "Algae", "penguins Gentoo",
                                               "lichens Parmelid", "lichens Physcid (shadow)", "springtails round",
                                               "mosses Polytrichales", "mosses Pottiales", "Petrels", "lichens Rhizocarpid",
                                               "Rotifers", "mites Sarcoptiformes", "lichens Stereocaulid", "Tardigrades",
                                               "lichens Teloschistid", "mites Trombidiformes", "lichens Umbilicarid"))

ecodat <- ecodat %>% left_join(bio_key)

## raw elevations ###
elev_table <- ecodat %>% filter(name == "elevation", !grepl(unit_h, pattern = "NA")) %>%
  group_by(unit_h) %>%  summarise(min = min(value), 
                                  max = max(value), 
                                  mean = mean(value), 
                                  lower_90 = quantile(value, 0.05), 
                                  upper_90 = quantile(value, 0.95))
#### Generate reports ####

count <- 0
for(i in sort(unique(typ_df$typV5))[3:5]) {
  count <- count + 1
  unitname <- sort(unique(data$unit_h))[i]
  unit <- typ_df %>% dplyr::filter(typV5 == i) %>% mutate(ecosystem = as.factor(typV5))
  rmarkdown::render('./documents/Ecosystem_Descriptions_doc.Rmd',  
                    output_file =  paste("report_", sort(unique(data$unit_h))[i], '_', Sys.Date(), ".docx", sep=''), 
                    output_dir = './documents/reports/short_doc')
  
  
}

#Render typology overview
rmarkdown::render('./documents/AntarcticTypology_Overview.Rmd',  
                  output_file =  "AntarcticTypology_Overview.doc", 
                  output_dir = './documents')


##### Supergroup descriptions #####
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



######### Regional units ###############
load("./data/spatialUnits.RData")
acbr <- raster("../Data/ACBRs/acbrs.tif")
pixels <- SpatialPointsDataFrame(coords = typ_df[,c("x", "y")], data = typ_df, proj4string = CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
acbrs_extract <- raster::extract(acbr, pixels, buffer =)
acbrs <- data.frame(pixels) %>% 
  dplyr::select(typV2_fa_hier_12v, x, y) %>% 
  mutate(ACBR_ID = acbrs_extract) %>% 
  merge(spatialUnits %>% dplyr::select(ACBR_ID, ACBR_Name) %>% 
          distinct(), by = "ACBR_ID")

#### for each RU, find ID of neighbouring pixels
acbrs$RU_raw <- paste(acbrs$ACBR_Name, acbrs$typV2_fa_hier_12v, sep = "_")
RUs <- paste(acbrs$ACBR_Name, acbrs$typV2_fa_hier_12v, sep = "_") %>% unique()
acbrs <- data$unit_h %>% unique() %>% sort() %>% data.frame(unit_ID = 1:36) %>% setNames(c("unitName", "unit_ID")) %>% 
  merge(acbrs, by.x = "unit_ID", by.y = "typV2_fa_hier_12v", all.y = TRUE) %>% 
  mutate(envir = word(unitName, 1,1, sep = "_"))

v2 <- list()
for(i in seq_along(RUs)){
  runit <- acbrs %>% filter(RU_raw == RUs[i])
  others <- acbrs %>% filter(!RU_raw == RUs[i]) 
  # units with > 100 pixels: only merge if same unit type in neighbouring ACBR has significant adjacent pixels
  if(nrow(runit) >= 100) others <- others %>% filter(!ACBR_Name == runit$ACBR_Name[1] & unit_ID == runit$unit_ID[1])
  # units with < 100 pixels: merge if same env group in same or neighbouring ACBR has significant adjacent pixels
  if(nrow(runit) < 100) others <- others %>% filter(envir == runit$envir[1])
  
  v2[[i]] <- apply(runit %>% dplyr::select(x, y), 1, function(z)  others %>% filter(x %in% c(z["x"], z["x"]+1000, z["x"]-1000), y %in% c(z["y"], z["y"]+1000, z["y"]-1000))) %>%
    bind_rows()
}


v3 <- v2 %>% purrr::map(~table(.$RU_raw) %>% cbind) %>% setNames(RUs)
x <- merge(table(acbrs$RU_raw) %>% cbind(), v2 %>% purrr::map(~table(.$RU_raw) %>% cbind) %>% setNames(RUs) %>% purrr::map_int(nrow) %>% cbind(), by = 0) %>% setNames(c("Runit", "Area", "Potential_merges"))

x$merger <- NA
x$adjpix <- NA
x$merger[which(x$Potential_merges == 0 & x$Area < 100)] <- "Discard"
x$merger[which(x$Potential_merges == 0 & x$Area >= 100)] <- "None"

for(i in seq_along(v3)){
  if(length(v3[[i]]) == 1){
    x$merger[which(x$Runit == names(v3[i]))] <- rownames(v3[[i]])
    x$adjpix[which(x$Runit == names(v3[i]))] <- v3[[i]][1]
  }else if(length(v3[[i]]) > 1){
    # check which matches are in the same unit
    region <- word(names(v3[i]), 1, 1, sep = "_")
    unit <- word(names(v3[i]), 2, 2, sep = "_")
    
    unitmatch <- grep(pattern = paste0("_", unit, "$"), rownames(v3[[i]]))
    regionmatch <- grep(pattern = region, rownames(v3[[i]]))
    
    best_unitmatch <- which.max(v3[[i]][unitmatch,])
    best_regionmatch <- which.max(v3[[i]][regionmatch,])
    
    if(length(unitmatch)==0){
      x$merger[which(x$Runit == names(v3[i]))] <-  rownames(v3[[i]])[regionmatch[best_regionmatch]]
      x$adjpix[which(x$Runit == names(v3[i]))] <- v3[[i]][regionmatch[best_regionmatch]]
    }else{
      if(length(regionmatch) == 0){
        x$merger[which(x$Runit == names(v3[i]))] <- rownames(v3[[i]])[unitmatch[best_unitmatch]]
        x$adjpix[which(x$Runit == names(v3[i]))] <- v3[[i]][unitmatch[best_unitmatch]]
      }else{
        if(v3[[i]][regionmatch[best_regionmatch]] > 2* v3[[i]][unitmatch[best_unitmatch]]) {
          x$merger[which(x$Runit == names(v3[i]))] <-  rownames(v3[[i]])[regionmatch[best_regionmatch]]
          x$adjpix[which(x$Runit == names(v3[i]))] <- v3[[i]][regionmatch[best_regionmatch]]
        }else{
          x$merger[which(x$Runit == names(v3[i]))] <- rownames(v3[[i]])[unitmatch[best_unitmatch]]
          x$adjpix[which(x$Runit == names(v3[i]))] <- v3[[i]][unitmatch[best_unitmatch]]
        }
      }
    }
  }
}


library(igraph)
x2 <- x %>% dplyr::select(Runit, merger, Area, weight = adjpix) %>% 
  filter(!merger %in% c("Discard"))

g <- graph_from_data_frame(x2, directed = TRUE)
g <- delete.vertices(g, "None")

plot(g, vertex.size = log(x2$Area), vertex.label.size = .5, vertex.label = NA)

cl <- cluster_edge_betweenness(g)

m <- numeric()
for(i in seq_along(cl)){m[i] <- length(cl[[i]])}


# i <- 14
# acbrs %>% filter(RU_raw %in% cl[[i]]) %>% ggplot(aes(x = x, y = y, col = RU_raw, fill = RU_raw)) + 
#   geom_tile(lwd = 1) + coord_equal()
# 
# 
# for(i in seq_along(cl)){
#  ggplot() + 
#     geom_tile(acbrs %>% filter(RU_raw %in% cl[[i]]), aes(x = x, y = y, col = RU_raw, fill = RU_raw), lwd = 1) + 
#     coord_equal() 
# }
# 



ext <- list()
bg <- list()
bg_tile <- list()
ru_merge <- list()
p <- list()

for(i in seq_along(cl)){
  
  ext[[i]] <- geo_bounds(acbrs %>% filter(RU_raw %in% cl[[i]]) %>% dplyr::select(y, x), buff = 0.15)
  bg[[i]] <- crop(antarctica, extent(ext[[i]]))
  bg_tile[[i]] <- df_crop(typ_df, ext[[i]])
  ru_merge[[i]] <- acbrs %>% filter(RU_raw %in% cl[[i]])
  
  p[[i]] <- ggplot() + 
    geom_polygon(data = bg[[i]], aes(x=long, y=lat, group=group), fill="gray50", color= NA, size=0.25) +
    geom_tile(data = bg_tile[[i]], aes(x = x, y = y), fill = "gray35", col = "gray35", lwd = 0.8) +
    geom_tile(data = ru_merge[[i]], aes(x = x, y = y, col = RU_raw, fill = RU_raw), lwd = 1.5) + 
    coord_equal() +
    theme(axis.title = element_blank(),
          plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
    scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) 
  
  save(ext, bg, bg_tile, ru_merge, p, file = "./results/RegionalUnitsPlots.RData")
}

# Generate report
load("./results/RegionalUnitsPlots.RData")
rmarkdown::render('./documents/Regional_units.Rmd',  
                  output_file =  "Regional_units.html", 
                  output_dir = './documents')
