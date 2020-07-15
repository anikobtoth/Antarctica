library(tidyverse)
library(igraph)
library(ggalluvial)
library(reshape2)
library(RColorBrewer)
source('./scripts/Helper_Functions.R')

# Load data ####
occurrences <- read_csv("data/Species/Ant_Terr_Bio_Data_FINAL.csv")

occurrences <- occurrences %>% select(scientificName, vernacularName, 
                                      decimalLongitude, decimalLatitude, 
                                      ACBR_ID, ASPA_ID, Date, year, Publish_YEAR, 
                                      kingdom, phylum , class, order, 
                                      family, genus, species, Snap_IFA, 
                                      Dist_IFA, coordinateUncertaintyInMetres, 
                                      individualCount)

## Species data ####
sppDat <- occurrences %>% select(scientificName, vernacularName, kingdom, phylum, 
                                 class, order, family, genus, species) %>% unique()

# Occurrence in ice-free areas (Ice-free patches as sites).
occ <- read_csv("data/Species/Spp_iceFree_occ.csv")


PA0 <- occ %>% select(scientific, OBJECTID, ACBR_Name, Functional_group, year) %>% 
  filter(#year > 1960 & 
    OBJECTID > 0, !Functional_group %in% c("Unknown", "Not assigned"))
PA0 <- PA0[-grep(" cf. " ,PA0$scientific),]
PA0 <- PA0[-grep(" sp." ,PA0$scientific),]

Abund <- PA0 %>% split(.$phylum) %>% 
  lapply(reshape2::dcast, scientific~OBJECTID, fun.aggregate = length, value.var = "phylum") %>%
  lapply(namerows)
Abund <- Abund[!map_lgl(Abund, ~dim(.) %>% is.null())] %>% purrr::map(clean.empty, mincol = 2, minrow = 2)
Abund <- Abund[map_lgl(Abund, ~nrow(.)>=2 && ncol(.) >= 5)]

PA <- tobinary(Abund) %>% lapply(clean.empty, mincol = 2, minrow = 2)

## pairs analysis on assemblages ####
pairs <- lapply(PA, simpairs)
el <- map2(pairs, PA, dist2edgelist) #%>% bind_rows(.id = "decade")

g <- el %>% #purrr::map(filter, Z.Score > quantile(Z.Score, 0.70, na.rm = T)) %>% 
  purrr::map(graph_from_data_frame, directed = F)
g <- lapply(g, function(x) delete_edges(x, which(E(x)$Z.Score< length(V(x))^(.16))))
g <- lapply(g, function(x) delete_vertices(x, which(degree(x) == 0)))

g <- g[map(g, ~length(E(.))) >4]
mod <- purrr::map(g, ~cluster_fast_greedy(.))#, weights = E(.)$Z.Score))
g <- map2(g, mod, function(x, y) {V(x)$cluster <- as.factor(y$membership) 
return(x)})

par(mfrow = c(2, 3), mar = c(0, 0, 1, 0), oma = c(1,1,1,1))
map2(g, names(g), plotnet_)

