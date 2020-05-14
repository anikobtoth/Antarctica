library(tidyverse)
library(igraph)
library(reshape2)
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

## Species data ####
sppDat <- occurrences %>% select(scientificName, vernacularName, kingdom, phylum, 
                                 class, order, family, genus, species) %>% unique()

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
pairs <- lapply(PA, simpairs)
el <- map2(pairs, PA, dist2edgelist) #%>% bind_rows(.id = "decade")

g <- el %>% purrr::map(filter, Z.Score > quantile(Z.Score, 0.98, na.rm = T)) %>% 
  #purrr::map(mutate, weight = Z.Score) %>%
  purrr::map(graph_from_data_frame, directed = F)

mod <- purrr::map(g, cluster_fast_greedy)

## Assemblage clusters similarity analysis #####

clust <- purrr::map(mod, ~data.frame(names = .$names, cluster = .$membership)) %>% 
  bind_rows(.id = "decade") %>% mutate(clustID = paste(decade, cluster, sep = "_"))
clustm <- dcast(clust, names~clustID, fun.aggregate = length, fill = 0) %>% namerows()

# ordination
cdist <- forbesMatrix(clustm) %>% as.dist(upper = F)
cord <- cmdscale(1-cdist, k = 2) %>% data.frame(cluster = colnames(clustm), .)
cord <- cord %>% mutate(decade = word(cluster, 1, 1, sep="_"), cID = word(cluster, 2, 2, sep="_"))

ggplot(cord, aes(x = X1, y = X2, col = decade, label = cID)) + geom_point() + geom_text(hjust = 0, vjust = 0)

#pairs
cpairs <- clustm %>% t() %>% simpairs()
cel <- dist2edgelist(cpairs, t(clustm))

cg <- cel %>% dplyr::filter(Z.Score > quantile(Z.Score, 0.98, na.rm = T)) %>% graph_from_data_frame(directed = F)
V(cg)$decade <- V(cg)$name %>% word(1,1, sep = "_")
cmod <- cluster_fast_greedy(cg)
plot(cg, vertex.label = NA, vertex.size = 4, vertex.color = cmod$membership)

cgroups <- data.frame(cluster = cmod$names, group = cmod$membership)
cgroups <- merge(cord, cgroups) %>% mutate(cID = as.numeric(cID))
## Characterise sites by assemblage ####

s_count <- map2(mod, PA, getSites, type = "count")
s_prcnt <- map2(mod, PA, getSites, type = "percent")

assemb1 <- lapply(s_count, lapply, data.frame) %>% purrr::map(~reduce(., multimerge))
assemb1 <- purrr::map(assemb1, function(x) x %>% setNames(paste("X", 1:ncol(x), sep = "")))
assemb1 <- assemb1 %>% purrr::map(~apply(.,1, function(x) which(x == max(x, na.rm = TRUE)) %>% 
                                data.frame()) %>% 
                                  bind_rows(.id = "IFA")) %>% bind_rows(.id = "decade")
asmb1 <- assemb1 %>% group_by(decade, IFA) %>% summarise(count = paste0(., collapse = ";"))

assemb2 <- lapply(s_prcnt, lapply, data.frame) %>% purrr::map(~reduce(., multimerge))
assemb2 <- purrr::map(assemb2, function(x) x %>% setNames(paste("X", 1:ncol(x), sep = "")))
assemb2 <- assemb2 %>% purrr::map(~apply(.,1, function(x) which(x == max(x, na.rm = TRUE)) %>% 
                                           data.frame()) %>% 
                                    bind_rows(.id = "IFA")) %>% bind_rows(.id = "decade")
asmb2 <- assemb2 %>% group_by(decade, IFA) %>% summarise(percent = paste0(., collapse = ";"))
assemblages <- merge(asmb1, asmb2, all = TRUE)


## PCoA on assemblages ####
dist <- map(PA, forbesMatrix) %>% map(as.dist, upper = F) %>% map(~return(1-.)) 
ord <- map(dist, cmdscale, k = 2) %>% map2(PA, function(x, y) data.frame(IFA = colnames(y), x))
ord <- bind_rows(ord, .id = "decade")

IFA <- merge(ord, assemblages, all = TRUE) %>% mutate(IFA = as.numeric(IFA))


cID <- map2(IFA$decade, IFA$count %>% strsplit(";"), paste, sep = "_")
IFA$grp.count <- lapply(cID, function(x) cgroups$group[which(cgroups$cluster %in% x)]) %>% 
  purrr::map(unique) %>% sapply(paste0, collapse = ";")
cID <- map2(IFA$decade, IFA$percent %>% strsplit(";"), paste, sep = "_")
IFA$grp.percent <- lapply(cID, function(x) cgroups$group[which(cgroups$cluster %in% x)]) %>% 
  purrr::map(unique) %>% sapply(paste0, collapse = ";")

IFA <- merge(IFA, occ %>% select(OBJECTID, ACBR_Name) %>% unique(), by.x = "IFA", by.y = "OBJECTID", all.x = TRUE)

## Species allocation in cluster groups #####
species <- melt(clustm %>% mutate(species = rownames(clustm))) %>% 
  filter(value > 0) %>% 
  merge(cgroups %>% select(cluster, group), by.x = "variable", by.y = "cluster", all = TRUE)

sppDat <- merge(sppDat, species %>% select(species, group), by.x = "scientificName", by.y = "species", all = TRUE)
# Calculate IFA habitats ####

hab_IFA <- read_csv("data/Habitats/Habitat_IFA_join.csv")
hab_IFA$Inv_Distance <- 1/(hab_IFA$Distance+1)

habifa <- dcast(hab_IFA, OBJECTID_1~n8_W, value.var = "Inv_Distance", fun.aggregate = sum) %>% namerows()
habifa$dom_habitat <- apply(habifa, 1, which.max)

master <- merge(habifa %>% select(dom_habitat), IFA, by.x =0, by.y = "IFA")

masterg <- master %>% group_by(dom_habitat, grp.count, decade) %>% summarise(count = length(count))
ggplot(masterg, aes(y = count, axis1 = dom_habitat, axis2 = grp.count)) + geom_alluvium(aes(fill = decade), width = 1/12) + geom_stratum(width = 1/12, fill = "gray", color = "white") + geom_label(stat = "stratum", infer.label= TRUE) + scale_x_discrete(limits = c("habitat", "assemblage"), expand = c(0.05, 0.05)) + scale_fill_brewer(type = "qual", palette = "Set1")
#####
#####
#####
#####