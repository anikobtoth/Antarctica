library(tidyverse)
occurrences <- read_csv("data/Species/Ant_Terr_Bio_Data_FINAL.csv")

occurrences <- occurrences %>% select(scientificName, vernacularName, 
                                      decimalLongitude, decimalLatitude, 
                                      ACBR_ID, ASPA_ID, Date, year, Publish_YEAR, 
                                      kingdom, phylum , class, order, 
                                      family, genus, species, Snap_IFA, 
                                      Dist_IFA, coordinateUncertaintyInMetres, 
                                      individualCount)
