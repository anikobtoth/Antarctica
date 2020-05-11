library(tidyverse)
library(igraph)
source('C:/Users/Aniko/Desktop/Antarctica/scripts/Helper_Functions.R')

# Load data ####
occurrences <- read_csv("data/Species/Ant_Terr_Bio_Data_FINAL.csv")

occurrences <- occurrences %>% select(scientificName, vernacularName, 
                                      decimalLongitude, decimalLatitude, 
                                      ACBR_ID, ASPA_ID, Date, year, Publish_YEAR, 
                                      kingdom, phylum , class, order, 
                                      family, genus, species, Snap_IFA, 
                                      Dist_IFA, coordinateUncertaintyInMetres, 
                                      individualCount)
# Occurrence in ice-free areas (Ice-free patches as sites).
occ <- read_csv("data/Species/Spp_iceFree_occ.csv")

PA0 <- occ %>% select(scientific, OBJECTID, ACBR_Name, year) %>% 
  filter(year > 1960 & OBJECTID > 0)

Abund <- PA0 %>% split(.$year %/% 10) %>% 
  lapply(reshape2::dcast, scientific~OBJECTID, fun.aggregate = length, value.var = "year") %>%
  lapply(namerows) %>%
  purrr::map(clean.empty, mincol = 5, minrow = 2)

PA <- tobinary(Abund) %>% lapply(clean.empty, mincol = 2, minrow = 2)

## pairs analysis on assemblages ####
test <- lapply(PA, simpairs)
el <- map2(test, PA, dist2edgelist) #%>% bind_rows(.id = "decade")

g <- el %>% purrr::map(filter, Z.Score > quantile(Z.Score, 0.975, na.rm = T)) %>% 
  #purrr::map(mutate, weight = Z.Score) %>%
  purrr::map(graph_from_data_frame, directed = F)

mod <- purrr::map(g, cluster_fast_greedy)

s <- map2(mod, PA, getSites)

assemb <- lapply(s, lapply, data.frame) %>% purrr::map(~reduce(., multimerge))
assemb <- purrr::map(assemb, function(x) x %>% setNames(paste("X", 1:ncol(x), sep = "")))

## PCoA on assemblages ####
dist <- map(PA, forbesMatrix) %>% map(as.dist, upper = F) %>% map(~return(1-.)) 
ord <- map(dist, cmdscale, k = 2) %>% map(data.frame)
ord <- bind_rows(ord, .id = "decade")
#ord <- map(ord, ~merge(., y = sitedat[,c(8:16)], all.x = T, all.y = F, by.x = 0, by.y = "siteid"))

#####
#####
#####