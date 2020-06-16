library(raster)
library(tidyverse)
library(rts)

l <- list.dirs("E:/Antarctica/Data/Species/final_results", recursive = FALSE) %>% 
  file.path("cif_basedist0.tif")
n <- list.dirs("E:/Antarctica/Data/Species/final_results", recursive = FALSE, full.names = FALSE)

SDMs <- stack(l)
v <- getValues(SDMs) %>% data.frame() %>% na.omit()

# Calculate thresholds at top 20% of output suitability range.
thresh <- sapply(1:34, function(i) scale(v[,i]) %>% range(na.rm = TRUE) %>% quantile(0.80))

# number of suitable pixels by group
sapply(1:34, function(i) length(which(scale(v[,i]) >= thresh[i]))) %>% setNames(n)

# get coordinates of pixels
pix <- xyFromCell(SDMs, rownames(v) %>% as.integer())

# Produce pseudo-presence table
PA <- lapply(1:34, function(i) scale(v[,i]) >= thresh[i]) %>% reduce(data.frame) %>%
  apply(2, as.numeric) %>% t() %>% data.frame()

pairs <- simpairs(PA)
el <- dist2edgelist(pairs, PA)
