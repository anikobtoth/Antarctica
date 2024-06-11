## Generate regional units
library(dplyr)
library(purrr)
library(tidyr)
library(terra)
library(sf)

#typ_df from ecosystem descriptions dataprep
acbrs <- st_read("./data/ACBRs", "ACBRs_v2_2016") %>% st_buffer(dist = 100) %>% filter(ACBR_ID != 2)
epsg3031 <- st_crs(acbrs)
ACBR_key <- readRDS("./data/ACBR_key.rds")

units <- typ_df %>% st_as_sf(coords = c("x", "y")) %>% st_set_crs(epsg3031) %>% 
  mutate(unit = word(LCODE, 1,1, sep = "B"))


# raw intersection of typology and ACBRs
ru <- st_intersection(units, acbrs, tolerance = 50)
ru$ru <- paste(ru$ACBR_Name, ru$LCODE, sep = "_")
excepted_rus <- ru %>% filter(grepl(.$unit, pattern = "G|L") | LCODE == "E1B8") ## exclude overlays
ru <- ru %>% filter(!(grepl(.$unit, pattern = "G[1-3]|L") | LCODE == "E1B8")) ## exclude overlays

rm(data, units, typ_df, descr)
# identify ACBR neighbours
acbr_nbrs <- acbrs %>% split(.$ACBR_ID) %>% 
  purrr::map(~.x %>% st_bbox() %>% st_as_sfc() %>% st_as_sf()) %>% 
  bind_rows(.id = "ACBR_ID") %>% 
  st_buffer(200) %>%
  st_intersection() %>%
  filter(n.overlaps > 1) %>% 
  st_set_crs(epsg3031)

merge.out <- ru_logicTree(ru, 5000)

final_runits <- merge.out[!map_lgl(merge.out, is.null)] %>% 
  reduce(makeunits) %>% arrange(final)
ru <- left_join(ru, final_runits, by = c("ru" = "raw"), )

# number regional units with no merges
flagged_small <- ru %>% filter(is.na(final))
flagged_small$final <- flagged_small$ru %>% as.factor() %>% as.numeric() %>% `+`(max(ru$final, na.rm = T))
ru <- ru %>% filter(!is.na(final)) %>% rbind(flagged_small)

## add excepted rus back in.
excepted_rus$final <- excepted_rus$ru %>% as.factor() %>% as.numeric() %>% `+`(max(ru$final, na.rm = T))
ru <- rbind(ru, excepted_rus)

saveRDS(ru, "./results/RegionalUnitsV6_3.rds")

ru_summary <- ru %>% st_drop_geometry() %>% group_by(typv6_v1pl, unit, ACBR_ID, ACBR_Name, ru, final) %>% summarise(n = n())
BETs <- ru_summary %>% group_by(final) %>% 
  summarise(ID_Tier1 = paste(unique(unit), collapse = ", "), 
            ACBR_Name = paste(unique(ACBR_Name), collapse = ", "), 
            ID_Tier2 = word(ru, 2, 2, sep = "_") %>% unique() %>% paste(collapse = ", "), 
            ID_Tier3 = ru[which(n == max(n))] %>% word(2,2, "_"), size = sum(n), n_merge = n()) %>%
  mutate(ACBR_Name = str_replace(ACBR_Name, "south", "South"), 
         ACBR_Name = str_replace(ACBR_Name, "Adelie", "Adélie"))
hc_names <- read_csv("documents/Unit_descriptionsV6.csv")
ACBR_abbrevs <- read_csv("documents/ACBR_abbrevs.csv")
BETs <- left_join(BETs, select(hc_names, LCODE, Name), by = c("ID_Tier3" = "LCODE")) %>% 
  left_join(select(hc_names, LCODE, Name), by = c("ID_Tier1" = "LCODE"), suffix = c("_Tier2", "_Tier1")) %>% 
  left_join(ACBR_abbrevs) %>%       
  mutate(amalg = ifelse(n_merge > 1, "a", ""))

# deal with amalgamated ACBRs
BETs[which(grepl(",", BETs$ACBR_Name)),]$Abbreviation <- BETs[which(grepl(",", BETs$ACBR_Name)),] %>% 
  pull(ACBR_Name) %>% strsplit(", ") %>% 
  map(~ACBR_abbrevs %>% filter(ACBR_Name %in% .x) %>% pull(Abbreviation)) %>% 
  map(paste0, collapse = "_") %>% unlist()

# full Tier 3 id 
BETs <- BETs %>% mutate(ID_Tier3 = paste0(ID_Tier3, "_", Abbreviation, amalg), 
                        Name_Tier3 = paste0(Name_Tier2, " in ", ACBR_Name)) %>% 
  select(BET_code = final, ID_Tier3, Name_Tier3, ID_Tier1, Name_Tier1, ID_Tier2, Name_Tier2, ACBR_Name) 

# BET summary table of names
saveRDS(BETs, "./Results/RegionalUnits_Naming_Table.rds")

# Regional Units helper functions
ru_logicTree <- function(ru, size_threshold = 2000){
  merge.out <- list()
  
  for(MEU in unique(ru$unit)){
    ru1 <- ru %>% filter(unit == MEU)
    
    rus <- unique(ru1$ru)
    merge.names <- list()
    
    for(i in 1:length(rus)){
      curr.runit <- ru1 %>% filter(ru == rus[i])
      curr.ecosys <- curr.runit$LCODE[1]
      curr.acbr <- curr.runit$ACBR_ID[1]
      search.area <- acbr_nbrs %>% filter(purrr::map_lgl(acbr_nbrs$origins, ~curr.acbr %in% .x))
      adjacent <- 0    
      
      # are there adjacent bioregions? 
      if(nrow(search.area) > 0){
        adjacent <- numeric(nrow(search.area)) %>% setNames(search.area$origins %>% reduce(c) %>% `[`(which(.!= curr.acbr)))
        # if so, are there adjacent pixels of the same unit in adjacent bioregions?
        for(j in 1:nrow(search.area)){
          curr.runit.sa <- st_intersection(curr.runit, search.area[j,])
          ru1.sa <- ru1 %>% filter(LCODE == curr.ecosys, ACBR_ID != curr.acbr) %>% 
            st_intersection(search.area[j,])
          adj.pix <- st_distance(curr.runit.sa, ru1.sa)
          adjacent[j] <- length(which(as.numeric(adj.pix) < 150))
        }}
      
      # if smaller runit, check for adjacent units in the same MEU and ACBR
      if(nrow(curr.runit) < size_threshold){
        MEU.units <- unique(ru1$LCODE)
        adjacentB <- rep(0, length(MEU.units)) %>% setNames(MEU.units)
        for(j in 1:length(MEU.units)){
          if(MEU.units[j] != curr.ecosys){
            ru1.acbr <- ru1 %>% filter(ACBR_ID == curr.acbr, LCODE == MEU.units[j])
            adj.pix <- st_distance(curr.runit, ru1.acbr)
            adjacentB[j] <- length(which(as.numeric(adj.pix) < 150))
          }
        }
      }
      
      if(nrow(curr.runit) >= size_threshold && sum(adjacent) > 0){
        # merge with highest adjacency neighbour.
        merge.acbr <- names(which.max(adjacent))
        merge.names[[i]] <- ru1$ru[which(ru1$LCODE == curr.ecosys & ru1$ACBR_ID %in% as.numeric(c(merge.acbr, curr.acbr)))] %>% unique()
      } else if (nrow(curr.runit) >= size_threshold && sum(adjacent) == 0){
        merge.names[[i]] <- unique(curr.runit$ru)
      } 
      
      
      if(nrow(curr.runit) < size_threshold && sum(adjacent) == 0) {
        if(max(adjacentB) > nrow(curr.runit)*0.05){ # at least 5% of unit pixels must be adjacent with merge candidate.
          # merge with highest adjacency unit in same acbr and MEU
          merge.ru <- names(which.max(adjacentB))
          merge.names[[i]] <- ru1$ru[which(ru1$LCODE %in% c(curr.ecosys, merge.ru) & ru1$ACBR_ID == curr.acbr)] %>% unique()
        }else{
          # Discard
          merge.names[[i]] <- NA
        }
      }
      if(nrow(curr.runit) < size_threshold && sum(adjacent) > 0){
        if(sum(adjacentB) == 0){
          # merge with highest adjacency A
          merge.acbr <- names(which.max(adjacent))
          merge.names[[i]] <- ru1$ru[which(ru1$LCODE == curr.ecosys & ru1$ACBR_ID %in% as.numeric(c(merge.acbr, curr.acbr)))] %>% unique()
        }else{
          if(max(adjacentB) > 2*max(adjacent)){
            # merge with most adjacent B
            merge.ru <- names(which.max(adjacentB))
            merge.names[[i]] <- ru1$ru[which(ru1$LCODE %in% c(curr.ecosys, merge.ru) & ru1$ACBR_ID == curr.acbr)] %>% unique()
          }else{
            # merge with most adjacent A
            merge.acbr <- names(which.max(adjacent))
            merge.names[[i]] <- ru1$ru[which(ru1$LCODE == curr.ecosys & ru1$ACBR_ID %in% as.numeric(c(merge.acbr, curr.acbr)))] %>% unique()
          }
        }
      }
    }
    merge.out <- c(merge.out, merge.names)
  }
  return(merge.out)
}

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

saveRDS(p, "./results/RegionalUnitsV6_plots.rds")

m <- numeric()
for(i in unique(ru$final) %>% na.omit %>% as.numeric){
  curr.ru <- ru %>% filter(final == i)
  m[i] <- length(unique(curr.ru$ru))
}

###
##
#