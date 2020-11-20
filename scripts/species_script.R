library(tidyverse)
library(igraph)
library(ggalluvial)
library(reshape2)
library(RColorBrewer)
library(vegan)
source('./scripts/Helper_Functions.R')

# Load data ####
occurrences <- read_csv("data/Species/Ant_Terr_Bio_Data_FINAL.csv")

occurrences <- occurrences %>% dplyr::select(scientificName, vernacularName, 
                                      decimalLongitude, decimalLatitude, 
                                      ACBR_ID, ASPA_ID, Date, year, Publish_YEAR, 
                                      kingdom, phylum , class, 
                                      order, family, genus, species, Snap_IFA, 
                                      Dist_IFA, coordinateUncertaintyInMetres, 
                                      individualCount)


# Occurrence in ice-free areas (Ice-free patches as sites).
occ <- read_csv("./data/Species/Spp_iceFree_occ.csv")

## Species data ####
sppDat <- occ %>% dplyr::select(scientific, vernacular, Functional_group, kingdom, phylum, 
                         class, order_, family, genus, species) %>% unique()


indepFG <- sppDat %>% filter(!Functional_group %in% good_models) %>% pull(Functional_group) %>% na.omit() %>% unique()
occurrences <- full_join(sppDat[,c("scientific", "Functional_group")], occurrences, by = c("scientific" = "scientificName")) %>% unique()
## Choose point occurrences that were not used in SDM factor analysis 
indepOCC <- occurrences %>% filter(Functional_group %in% indepFG)

indocc <- SpatialPointsDataFrame(coords = indepOCC[,c("decimalLongitude", "decimalLatitude")],
                                 data = indepOCC, 
                                 proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

indocc <- spTransform(indocc, "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

library(raster)
typ_fah <- raster("../Data/Typology/typV2_fa_hier_12v.tif")
typ_fad <- raster("../Data/Typology/typV2_fa_dual_12v.tif")
typ_far <- raster("../Data/Typology/typV2_fa_revhier_12v.tif")

indepOCC$fah <- raster::extract(typ_fah, indocc)
indepOCC$fad <- raster::extract(typ_fad, indocc)
indepOCC$far <- raster::extract(typ_far, indocc)

PA_fah <- indepOCC %>% reshape2::dcast(scientific~fah, fun.aggregate = length) %>% namerows()
PA_fad <- indepOCC %>% reshape2::dcast(scientific~fad, fun.aggregate = length) %>% namerows()
PA_far <- indepOCC %>% reshape2::dcast(scientific~far, fun.aggregate = length) %>% namerows()

# zero inflated poisson model for a test case
M3 <- zeroinfl(`1` ~ `2` | ## Predictor for the Poisson process
                 `2`, ## Predictor for the Bernoulli process;
               dist = 'poisson',
               data = PA_fah)

E2 <- resid(M3, type = "pearson")
N  <- nrow(PA_fah)
p  <- length(coef(M3)) + 1  # '+1' is for variance parameter in NB
sum(E2^2) / (N - p)

#### Old stuff #####
#### Make PA tables with species against ice free polygons ####
PA0 <- occ %>% select(scientific, Functional_group, OBJECTID, ACBR_Name, phylum, year) %>% 
  filter(#year > 1960 & 
    OBJECTID > 0, !phylum %in% c("Unknown", "Not assigned"))
PA0 <- PA0[-grep(" cf. " ,PA0$scientific),]
PA0 <- PA0[-grep(" sp." ,PA0$scientific),]

Abund <- PA0 %>% split(.$year %/% 10) %>% 
  lapply(reshape2::dcast, scientific~OBJECTID, fun.aggregate = length, value.var = "year") %>%
  lapply(namerows) %>%
  purrr::map(clean.empty, mincol = 5, minrow = 2)

Abund <- Abund[!map_lgl(Abund, ~dim(.) %>% is.null())] %>% purrr::map(clean.empty, mincol = 2, minrow = 2)
Abund <- Abund[map_lgl(Abund, ~nrow(.)>=2 && ncol(.) >= 5)]
#Abund <- Abund[!map_lgl(Abund, is.null)]
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


## Assemblage clusters similarity analysis #####

clust <- purrr::map(mod, ~data.frame(names = .$names, cluster = .$membership)) %>% 
  bind_rows(.id = "decade") %>% mutate(clustID = paste(decade, cluster, sep = "_"))
clustm <- dcast(clust, names~clustID, fun.aggregate = length, fill = 0, value.var = "clustID") %>% namerows()

# ordination
cdist <- forbesMatrix(clustm) %>% as.dist(upper = F)
cord <- cmdscale(1-cdist, k = 2) %>% data.frame(cluster = colnames(clustm), .)
cord <- cord %>% mutate(decade = word(cluster, 1, 1, sep="_"), cID = word(cluster, 2, 2, sep="_"))

ggplot(cord, aes(x = X1, y = X2, col = decade, label = cID)) + geom_point() + geom_text(hjust = 0, vjust = 0) +
  ggtitle("PCoA of clusters by decade based on species composition")

#pairs
cpairs <- clustm %>% t() %>% simpairs()
cel <- dist2edgelist(cpairs, t(clustm))

cg <- cel %>% dplyr::filter(Z.Score > quantile(Z.Score, 0.84, na.rm = T)) %>% 
  graph_from_data_frame(directed = F)
V(cg)$decade <- V(cg)$name %>% word(1,1, sep = "_")
cmod <- cluster_fast_greedy(cg, weights = E(cg)$Z.Score)

par(mfrow = c(1,1))
plotnet(cg, cmod, "Cluster analysis of species clusters")

cgroups <- data.frame(cluster = cmod$names, group = cmod$membership)
cgroups <- merge(cord, cgroups) %>% mutate(cID = as.numeric(cID))

ggplot(cgroups, aes(x = X1, y = X2, col = as.factor(group), label = decade)) + geom_point() + geom_text(hjust = 0, vjust = 0)

## Characterise sites by assemblage ####

s_count <- map2(mod, PA, getSites, type = "count")
s_prcnt <- map2(mod, PA, getSites, type = "percent")

assemb1 <- lapply(s_count, lapply, data.frame) %>% purrr::map(~reduce(., multimerge))
assemb1 <- purrr::map(assemb1, function(x) x %>% setNames(paste("X", 1:ncol(x), sep = "")))
assemb1 <- assemb1 %>% purrr::map(~apply(.,1, function(x) which(x == max(x, na.rm = TRUE)) %>% 
                                data.frame()) %>% 
                                  bind_rows(.id = "IFA")) %>% bind_rows(.id = "decade")
asmb1 <- assemb1 %>% group_by(decade, IFA) %>% summarise(count = paste0(sort(.), collapse = ";"))


lapply(s_count, lapply, data.frame) %>% lapply(lapply, rows2name) %>% lapply(bind_rows, .id = "cluster") %>% bind_rows(.id = "decade") %>% unite("id", decade, cluster, sep = "_") %>% spread(key = id, value = "X..i..")

assemb2 <- lapply(s_prcnt, lapply, data.frame) %>% purrr::map(~reduce(., multimerge))
assemb2 <- purrr::map(assemb2, function(x) x %>% setNames(paste("X", 1:ncol(x), sep = "")))
assemb2 <- assemb2 %>% purrr::map(~apply(.,1, function(x) which(x == max(x, na.rm = TRUE)) %>% 
                                           data.frame()) %>% 
                                    bind_rows(.id = "IFA")) %>% bind_rows(.id = "decade")
asmb2 <- assemb2 %>% group_by(decade, IFA) %>% summarise(percent = paste0(sort(.), collapse = ";"))
assemblages <- merge(asmb1, asmb2, all = TRUE)


## PCoA on assemblages ####
dist <- map(PA, forbesMatrix) %>% map(as.dist, upper = F) %>% map(~return(1-.)) 
ord <- map(dist, cmdscale, k = 2) %>% map2(PA, function(x, y) data.frame(IFA = colnames(y), x))
ord <- bind_rows(ord, .id = "decade")

IFA <- merge(ord, assemblages, all = TRUE) %>% mutate(IFA = as.numeric(IFA))


cID <- map2(IFA$decade, IFA$count %>% strsplit(";"), paste, sep = "_")
IFA$grp.count <- lapply(cID, function(x) cgroups$group[which(cgroups$cluster %in% x)]) %>% 
  purrr::map(unique) %>% lapply(sort) %>% sapply(paste0, collapse = ";")
cID <- map2(IFA$decade, IFA$percent %>% strsplit(";"), paste, sep = "_")
IFA$grp.percent <- lapply(cID, function(x) cgroups$group[which(cgroups$cluster %in% x)]) %>% 
  purrr::map(unique) %>% lapply(sort) %>% sapply(paste0, collapse = ";")

IFA <- merge(IFA, occ %>% select(OBJECTID, ACBR_Name) %>% unique(), by.x = "IFA", by.y = "OBJECTID", all.x = TRUE)

## Species allocation in cluster groups #####
species <- melt(clustm %>% mutate(species = rownames(clustm))) %>% 
  filter(value > 0) %>% 
  merge(cgroups %>% select(cluster, group), by.x = "variable", by.y = "cluster", all = TRUE)


sppDat <- merge(sppDat, species %>% select(species, group), by.x = "scientificName", by.y = "species", all = TRUE) %>% unique()
table(sppDat$phylum, sppDat$group)[which(rowSums(table(sppDat$phylum, sppDat$group))>0),]

# Calculate IFA habitats ####

hab_IFA <- read_csv("./data/Habitats/Habitat_IFA_join.csv")
hab_IFA$Inv_Distance <- 1/(hab_IFA$Distance+1)

habifa <- dcast(hab_IFA, OBJECTID_1~n8_W, value.var = "Inv_Distance", fun.aggregate = sum) %>% namerows()
habifa$dom_habitat <- apply(habifa, 1, which.max)

master <- merge(habifa %>% select(dom_habitat), IFA, by.x =0, by.y = "IFA")
master$ecogroup <- paste0("hab", master$dom_habitat, "_grp", master$grp.count)

master$ecogroup <- paste0("hab", master$dom_habitat, "_asm", master$grp.count)

masterg <- master %>% group_by(dom_habitat, grp.count, decade) %>% summarise(count = length(count))


ggplot(masterg %>% filter(!grp.count == ""), aes(y = count, axis1 = dom_habitat, axis2 = grp.count)) + 
  geom_alluvium(aes(fill = decade), width = 1/12) + 
  geom_stratum(width = 1/12, fill = "gray", color = "white") + 
  geom_label(stat = "stratum", infer.label= TRUE) + 
  scale_x_discrete(limits = c("habitat", "assemblage"), expand = c(0.05, 0.05)) + 
  scale_fill_brewer(type = "qual", palette = "Set1")+
  ggtitle("Distribution of typical assemblages in typical habitats by decade")


## Classify IFAs that weren't used to generate groups.
groups <- split(sppDat$scientificName, f = sppDat$group)

# TODO: redo Abund without getting rid of undated occurrences; may help classify some IFAs, even if it's without any temporal resolution.
PA <- tobinary(Abund) %>% lapply(clean.empty)

PA0 <- occ %>% select(scientific, OBJECTID, ACBR_Name, year) 
PA0 <- PA0[-grep(" cf. " ,PA0$scientific),]
PA0 <- PA0[-grep(" sp." ,PA0$scientific),]

PA <- PA0 %>% reshape2::dcast(scientific~OBJECTID, fun.aggregate = length, value.var = "year") %>% namerows() %>% clean.empty()

test <- purrr::map_int(PA, ~getGroups(groups, rownames(PA[which(.>0),])) %>% 
             which.max) %>% data.frame(IFA = names(.))
names(test)[1] <- "group"

ecogroups <- merge(test, habifa %>% select(dom_habitat), by.x = "IFA", by.y = 0, all.x = TRUE)
ecogroups$ecogroup <- paste0("hab", ecogroups$dom_habitat, "_grp", ecogroups$group)

# All together, Not by decade ####

PA1 <- PA0 %>% reshape2::dcast(Functional_group~OBJECTID, fun.aggregate = length, value.var = "year") %>% namerows() %>% clean.empty(mincol = 5, minrow = 6)
pairs1 <- simpairs(PA1)
el1 <- dist2edgelist(pairs1, PA1)
g1 <- el1 %>% filter(Z.Score > quantile(Z.Score, 0.50, na.rm = T)) %>% 
  graph_from_data_frame(directed = F)
g1 <- delete_edges(g1, e = E(g1)[which(E(g1)$Z.Score < quantile(E(g1)$Z.Score, 0.90, na.rm = T))])
l <- layout_with_graphopt(g1)
plot(g1, vertex.label = NA, vertex.size = 4, layout = l)
cl1 <- cluster_fast_greedy(g1, weights = E(g1)$Z.Score)

plot(g1, vertex.label = NA, vertex.size = 4, layout = l, vertex.color = cl1$membership)

##### Sub-clusters within major assemblage groups ####
PAg <- lapply(groups, function(x) PA[x,]) %>% lapply(clean.empty)

pairsg <- lapply(PAg, simpairs)
elg <- map2(pairsg, PAg, dist2edgelist)
gg <- elg %>% purrr::map(filter, Z.Score > quantile(Z.Score, 0.9, na.rm = T)) %>% 
  purrr::map(graph_from_data_frame, directed = F)

modg <- purrr::map(gg, ~cluster_fast_greedy(., weights = E(.)$Z.Score))


#####