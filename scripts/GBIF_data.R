#skip <- download_GBIF_all_species(spp_list, path = "../Data/GBIF/")

GBIF.ALL <- list.files("../Data/GBIF/", pattern = ".RData") %>%   
  
  ## Restrict them to just the strings that match the species list for each run
  subset(. %in% paste(spp_list, "_GBIF_records.RData", sep = "")) %>%
  
  ## Pipe the list into lapply
  lapply(load.GBIF, GBIF_path = "../Data/GBIF/") %>% bind_rows 



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1.FILTER RECORDS BEFORE 1950 OR LACKING COORDINATES ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

## Now filter the records using conditions which are not too restrictive
GBIF.CLEAN <- GBIF.ALL %>% 
  
  ## These filters are very forgiving...
  filter(!is.na(lon) & !is.na(lat),
         
         ## CULTIVATED == 'CULTIVATED', 
         ## #establishmentMeans!='MANAGED' | is.na(establishmentMeans),
         ## New.Taxonomic.status == 'Accepted'
         
         year >= 1950 & !is.na(year),
         
         # ensure lat and lon are realistic
         lon >= -180 & lon <= 180, 
         lat >= -90  & lat <=90)

## How many records were removed by filtering?
message(round((nrow(GBIF.CLEAN))/nrow(GBIF.ALL)*100, 2), 
        " % records retained using spatially valid records")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## 2). GIVE EACH ROW UNIQUE ID ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

GBIF.CLEAN$CC.OBS <- paste0(1:nrow(GBIF.CLEAN), "_CC_", GBIF.CLEAN$scientificName)
GBIF.CLEAN$CC.OBS <- gsub(" ",     "_",  GBIF.CLEAN$CC.OBS, perl = TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## 3). FLAG SUSPICIOUS RECORDS ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(CoordinateCleaner)
library(timetk)

COORDCLEAN <- GBIF.CLEAN %>% 
  
  dplyr::rename(coord_spp        = scientificName,
                decimallongitude = lon, 
                decimallatitude  = lat) %>%
  
  ## The create a tibble for running the spatial outlier cleaning
  as_tibble() %>% 
  
  ## Flag suspicious coordinates
  clean_coordinates(.,
                    verbose         = TRUE,
                    tests = c("capitals", "centroids", "equal", "gbif", 
                              "institutions", "zeros"), ## duplicates flagged too many
                    
                    capitals_rad    = 5000,        ## remove records within 5km  of capitals
                    centroids_rad   = 5000,        ## remove records within 5km of country centroids
                    value           = "spatialvalid") %>% data.frame()


CLEAN.TRUE  = subset(COORDCLEAN, .summary == "TRUE")
CLEAN.FALSE = subset(COORDCLEAN, .summary == "FALSE")

## What percentage of records are retained?
message(round(nrow(CLEAN.TRUE)/nrow(COORDCLEAN)*100, 2), " % records retained")      

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## 4).FLAG OUTLIERS ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

GBIF.SPAT.ALL <- CLEAN.TRUE %>% as.data.frame() %>% 
  
  dplyr::select(searchTaxon, species,decimallatitude, decimallongitude) %>% 
  
  cc_outl(GBIF.SPAT.ALL, method  = "quantile",
                         mltpl   = 5,
                         tdi     = 500,
                         value   = "clean",
                         verbose = "TRUE")

saveRDS(GBIF.SPAT.ALL, file = "./data/Species/GBIF_clean_data.rds", compress = TRUE)


