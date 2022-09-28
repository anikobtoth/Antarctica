## Generate regional units
library(dplyr)
library(purrr)
library(tidyr)
library(raster)
library(sf)

epsg3031 <- CRS("+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
out <- readRDS("./results/out_TYP_hier_units_V5.rds")
acbrs <- st_read("./data/ACBRs", "ACBRs_v2_2016")
ACBR_key <- readRDS("./data/ACBR_key.rds")

units <- out %>% SpatialPointsDataFrame(coords = .[,c('x', 'y')], data = ., proj4string = epsg3031) %>%
  st_as_sf(crs = epsg3031)

# raw intersection of typology and ACBRs
ru <- st_intersection(units, acbrs, tolerance = 50)
ru$ru <- paste(ru$ACBR_Name, ru$unit_h, sep = "_")

# identify ACBR neighbours
acbr_nbrs <- acbrs %>% split(.$ACBR_ID) %>% 
  purrr::map(~.x %>% st_bbox() %>% extent() %>% as("SpatialPolygons")) %>% 
  reduce(rbind, makeUniqueIDs = TRUE) %>% 
  st_as_sf() %>% 
  st_buffer(200) %>%
  st_intersection() %>%
  filter(n.overlaps > 1) %>% 
  st_set_crs(epsg3031)


#ru_meds <- ru %>% split(.$consensus)
merge.out <- list()
#ru1 <- ru_meds[[1]]
for(MED in 1:length(unique(ru$consensus))){
  ru1 <- ru %>% filter(consensus == MED)
  
  rus <- unique(ru1$ru)
  merge.names <- list()
  
  for(i in 1:length(rus)){
    curr.runit <- ru1 %>% filter(ru == rus[i])
    curr.ecosys <- curr.runit$unit_h[1]
    curr.acbr <- curr.runit$ACBR_ID[1]
    search.area <- acbr_nbrs %>% filter(purrr::map_lgl(acbr_nbrs$origins, ~curr.acbr %in% .x))
    adjacent <- 0    
    
    # are there adjacent bioregions? 
    if(nrow(search.area) > 0){
      adjacent <- numeric(nrow(search.area)) %>% setNames(search.area$origins %>% reduce(c) %>% `[`(which(.!= curr.acbr)))
      # if so, are there adjacent pixels of the same unit in adjacent bioregions?
      for(j in 1:nrow(search.area)){
        curr.runit.sa <- st_intersection(curr.runit, search.area[j,])
        ru1.sa <- ru1 %>% filter(unit_h == curr.ecosys, ACBR_ID != curr.acbr) %>% 
          st_intersection(search.area[j,])
        adj.pix <- st_distance(curr.runit.sa, ru1.sa)
        adjacent[j] <- length(which(as.numeric(adj.pix) < 150))
      }}
    
    # if smaller runit, check for adjacent units in the same MED and ACBR
    if(nrow(curr.runit) < 5000){
      med.units <- unique(ru1$unit_h)
      adjacentB <- rep(0, length(med.units)) %>% setNames(med.units)
      for(j in 1:length(med.units)){
        if(med.units[j] != curr.ecosys){
          ru1.acbr <- ru1 %>% filter(ACBR_ID == curr.acbr, unit_h == med.units[j])
          adj.pix <- st_distance(curr.runit, ru1.acbr)
          adjacentB[j] <- length(which(as.numeric(adj.pix) < 150))
        }
      }
    }
    
    if(nrow(curr.runit) >= 5000 && sum(adjacent) > 0){
      # merge with highest adjacency neighbour.
      merge.acbr <- names(which.max(adjacent))
      merge.names[[i]] <- ru1$ru[which(ru1$unit_h == curr.ecosys & ru1$ACBR_ID %in% as.numeric(c(merge.acbr, curr.acbr)))] %>% unique()
    } else if (nrow(curr.runit) >= 5000 && sum(adjacent) == 0){
      merge.names[[i]] <- unique(curr.runit$ru)
    } 
    
    if(nrow(curr.runit) < 5000 && sum(adjacent) == 0) {
      if(sum(adjacentB) > 0){
        # merge with highest adjacency unit in same acbr and med
        merge.ru <- names(which.max(adjacentB))
        merge.names[[i]] <- ru1$ru[which(ru1$unit_h %in% c(curr.ecosys, merge.ru) & ru1$ACBR_ID == curr.acbr)] %>% unique()
      }else{
        # Discard
        merge.names[[i]] <- NA
      }
    }
    if(nrow(curr.runit) < 5000 && sum(adjacent) > 0){
      if(sum(adjacentB) == 0){
        # merge with highest adjacency A
        merge.acbr <- names(which.max(adjacent))
        merge.names[[i]] <- ru1$ru[which(ru1$unit_h == curr.ecosys & ru1$ACBR_ID %in% as.numeric(c(merge.acbr, curr.acbr)))] %>% unique()
      }else{
        if(max(adjacentB) > 2*max(adjacent)){
          # merge with most adjacent B
          merge.ru <- names(which.max(adjacentB))
          merge.names[[i]] <- ru1$ru[which(ru1$unit_h %in% c(curr.ecosys, merge.ru) & ru1$ACBR_ID == curr.acbr)] %>% unique()
        }else{
          # merge with most adjacent A
          merge.acbr <- names(which.max(adjacent))
          merge.names[[i]] <- ru1$ru[which(ru1$unit_h == curr.ecosys & ru1$ACBR_ID %in% as.numeric(c(merge.acbr, curr.acbr)))] %>% unique()
        }
      }
    }
  }
  merge.out <- c(merge.out, merge.names)
}

final_runits <- merge.out %>% reduce(makeunits) %>% arrange(final)
ru <- left_join(ru, final_runits, by = c("ru" = "raw"), )
saveRDS(ru, "./results/RegionalUnitsV5.rds")

makeunits <- function(a, b){
  #if(is.na(a)) {return(data.frame("raw" = b, "final" = 1))}
  if(is_character(a)){a <- data.frame("raw" = a, "final" = 1)}
  
  if(any(b %in% a$raw)){
    f <- a$final[which(a$raw %in% b)] %>% unique()
    if(length(f) == 1) b <- data.frame("raw" = b, "final" = f)
    if(length(f) > 1) {
      a$final[which(a$final == max(f))] <- min(f)
      b <- data.frame("raw" = b, "final" = min(f))}
  }else{
    b <- data.frame("raw" = b, "final" = max(a$final)+1)
  }
  return(rbind(a, b) %>% distinct())
}


## plotting

antarctica <- st_read("./data/Base", "Coastline_high_res_polygon_v7.1") %>% st_simplify()

for(i in unique(ru$final) %>% na.omit %>% as.numeric){
  curr.ru <- ru %>% filter(final == i)
  ext <- geo_bounds(curr.ru)
  bg <- df_crop(ru, ext, 5)
  bgpoly <- st_crop(antarctica, ext)

  
  p[[i]] <- ggplot() +  
    geom_sf(data = bgpoly, fill="gray50", color= NA, size=0.25) +
    geom_tile(data = bg, aes(x = x, y = y), fill = "gray35", col = "gray35", lwd = 0.8) +
    geom_tile(data = curr.ru, aes(x=x, y=y, col= ru, fill = ru), lwd = .8) +
    theme(plot.margin = unit(c(0,0,0,0), "cm"), 
          axis.text = element_text(size = 7), 
          axis.title = element_blank()) +
    scale_y_continuous(labels = scales::scientific) + scale_x_continuous(labels = scales::scientific) 
  
}

saveRDS(p, "./results/RegionalUnitsV5_plots.rds")

m <- numeric()
for(i in unique(ru$final) %>% na.omit %>% as.numeric){
  curr.ru <- ru %>% filter(final == i)
  m[i] <- length(unique(curr.ru$ru))
}

###
##
#