
library(sp)
library(fields)
library(BBmisc)
library(psych)
library(paran)
library(rlang)
#library(gmp)


##### DATA MANIPULATION #####
### Change character vectors in df to factors
tofac <- function(df){
  df[,map_lgl(df, is.character)] <- map(df[,map_lgl(df, is.character)], factor)
  return(df)
}

### Abundance to presence-absence 
# matrix list
tobinary <- function(PA.LIST){
  binary <- lapply(PA.LIST, function(x) {
    x <- x/x
    x[is.na(x)] <- 0
    return(x)
  })
  return(binary)
}

# single matrix
tobinary.single <- function(x) {
  x <- x/x
  x[is.na(x)] <- 0
  return(x)}

# First column to rownames
namerows <- function(table){
  rownames(table) <- table[,1]
  table <- table[,2:ncol(table)]
  return(table)
}

# Rowames to first column
rows2name <- function(table){
  table <- data.frame(name = rownames(table), table)
  table$name <- as.character(table$name)
  
  return(table)
}

# remove empty rows and columns in 1 matrix
# or remove rows and columns with too few observations
clean.empty <- function(x, mincol = 1, minrow = 1){
  x <- x[which(rowSums(x) >= minrow),which(colSums(x) >= mincol)]
  return(x)
}

# triangular distance matrix to long format
dist2edgelist <- function(z, sppDat){  #edge list with link types attached
  k3 = as.matrix(z)
  dimnames(k3) <- list(rownames(sppDat), rownames(sppDat)) 
  
  xy <- t(combn(colnames(k3), 2))
  k3 <- data.frame(xy, dist=k3[xy], stringsAsFactors = F)
  
  k3 <- data.frame(k3, id = paste(k3$X1, k3$X2, sep = "-"), stringsAsFactors = F)
  colnames(k3) <- c("Sp1", "Sp2", "Z.Score", "id")
  
  return(k3)
}

## process tables into text format
df2txt <- function(df, threshold = 1.25, name = "name", type = c("a", "b")){
  if(type == 'a'){
    if(threshold > 1) vars <- df[df$envgr_cont > threshold, name] %>% unlist() %>% as.character()
    if(threshold < 1)vars <- df[df$envgr_cont < threshold, name] %>% unlist() %>% as.character()
  }
  if(type == 'b'){
    if(threshold > 1) vars <- df[df$unit_envgr > threshold, name] %>% unlist() %>% as.character()
    if(threshold < 1)vars <- df[df$unit_envgr < threshold, name] %>% unlist() %>% as.character()
  }
  
  #clean var names
  vars <- vars %>% str_replace_all("_", " ") %>% trimws() %>% str_replace_all(" ", "_")
  
  if(length(vars) == 1) return(vars)
  if(length(vars) == 0 ) return ("no variables")
  if(length(vars) == 2) return(paste(vars[1:length(vars)-1], "and", vars[length(vars)]))
  if(length(vars) > 2) return(paste(vars[1:(length(vars)-1)], collapse = ", ") %>% paste("and", vars[length(vars)]))
}


# Functions to assist with penguin rookery size calculations
rm_outliers <- function(data){
  outliers <- boxplot(data, plot = FALSE)$out
  data_no_outlier <- data[-which(data %in% outliers)]
  return(data_no_outlier)
}

adjust_count <- function(data){
  data %>% log() %>% rm_outliers() %>% exp()
}

BP_translate <- function(nests, adults, chicks, ratios){
  nests[which(nests == 0)] <- NA
  chicks[which(chicks == 0)] <- NA
  adults[which(adults == 0)] <- NA
  
  BP <- nests
  BP <- ifelse(is.na(BP), chicks/ratios["cn"], BP)
  BP <- ifelse(is.na(BP), adults/ratios["an"], BP)
  
  return(BP)
  
}

# Outliers function
far_outliers <- function(v){
  c(which(v < quantile(v, 0.25) - (3.5* IQR(v))),
    which(v > quantile(v, 0.75) + (3.5* IQR(v))))
}

##### DATA EXTRACTION ####
# list of good model names 
gm <- function(){
  n <- list.files("../Data/Species/final_results", ".tif$", recursive = FALSE, full.names = FALSE)
  bad_models <- c("adeliae.tif","Procellariiformes.tif","Poduromorpha.tif", "Grimmiales.tif",
                  "Cyanobacteria.tif", "Tardigrada.tif","Marchantiophyta.tif","Lecideaceae.tif",
                  "Umbilicariaceae.tif")
  return(n[!n%in% bad_models])
}

# Extract abiotic data
extract_env_dat <- function(xy, datapath){
  abiotic <- c("cloud", "wind", "meanTemp", "melt", "modT", "aspect", "elevation", "rock", "rugosity", "slope", "solar", "DDm5", "totPrecip")
  l <- list.files(datapath, ".tif$", full.names = T)
  layers <- purrr::map(l, raster) %>% setNames(abiotic)
  
  abiotic_data <- purrr::map(layers, ~raster::extract(.x, xy))
  abiotic_data <- cbind(xy, reduce(abiotic_data, data.frame) %>% 
                          set_names(abiotic)) %>% tibble()
  return(abiotic_data)
}
# Extract biotic suitability data
extract_biotic_dat <- function(xy, good_models, datapath){
  l <- list.files(datapath, ".tif$", recursive = F, full.names = T)
  l <- l[list.files(datapath, ".tif$", recursive = F, full.names = F) %in% good_models]
  layers <- raster::stack(l) 
  projection(layers) <- "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m
+no_defs"
  
  biotic_data <- raster::extract(layers, xy)
  biotic_data <- cbind(xy, data.frame(biotic_data) %>% set_names(good_models)) %>% 
    tibble()
  return(biotic_data)
}

##### ANALYSES ######
#wrapper for fa() that scales the data and saves the results to file
factor_analysis <- function(dat, name, nfact, scale = TRUE){
  message("Preparing data")
  dat <- dat %>% na.omit()
  rn <- rownames(dat)
  if(scale) dat <- sapply(dat, function(x) normalize(x, method = "range", range = c(1,100)) %>% log()) %>% data.frame()
  
  #paranal <- paran(dat, cfa = TRUE, graph = TRUE, color = TRUE, centile = 95, iterations = 5000)
  #nfact <- paranal$Retained
  
  message(paste("Factor analysis with", nfact, "factors"))  
  fa1 <- fa(dat, nfact, rotate = "varimax")
  saveRDS(fa1, paste0("./results/factor_analyses/fa_", name, ".rds"))
  saveRDS(dat, paste0("./results/factor_analyses/data_", name, ".rds"))
  #saveRDS(paranal, paste0("./results/factor_analyses/paran_", name, ".rds"))
  sc <- data.frame(fa1$scores)
  
  message("Calculating categories and confidence")
  #select best fit and calculate confidence as the ratio of best minus second-best fit to best minus worst fit. 
  consensus <- apply(sc, 1, which.max) %>% factor() %>% setNames(rn)
  confidence <- apply(sc, 1, function(x) diff(sort(exp(x)))) %>% apply(2, function(x) last(x)/sum(x)) 
  
  return(data.frame(consensus, confidence))
}

# wrapper for predict function that cleans and scales the 
# input data while saving pixel IDs and then uses input fa object to predict scores
fapred <- function(data, cols, faobj, olddat){
  message("formatting data")
  cldat <- data %>% dplyr::select(pixID, all_of(cols)) %>% na.omit() 
  pixID <- cldat$pixID
  cldat <- cldat %>% dplyr::select(-pixID) %>% 
    sapply(function(x) normalize(x, method = "range", range = c(1,100)) %>% log()) %>% 
    data.frame()
  message("calculating predictions") 
  typ <- predict(faobj, data = cldat, old.data = olddat) %>% 
    apply(1, which.max) %>% data.frame(pixID, typ = .)
  
  return(typ)
}

# classifies unclassified pixels in column "var" with nearest neighbour
classify_by_neighbours <- function(dat, var, maxdist = 1.5, res = 100){
  library(fields)
  
  v0 <- dat %>% filter(is.na({{ var }})) 
  v2 <- dat %>% filter(!is.na({{ var }}))
  
  if(nrow(v0) > 0){
    v2 <- v2 %>% filter(x %in% unique(v0$x, v0$x+res, v0$x - res) & y %in% unique(v0$y, v0$y+res, v0$y - res))
    temp <- v2[rdist.earth(v0[,c("lon", "lat")], v2[,c("lon", "lat")], miles = F) %>% apply(1, which.min),] %>% 
      dplyr::select(x, y, lon, lat, {{var}}) %>% 
      mutate(dist = rdist.earth(v0[,c("lon", "lat")], v2[,c("lon", "lat")], miles = F) %>% apply(1, min)) %>% dplyr::select({{var}}, dist)
    temp <- data.frame(x = v0$x, y = v0$y, temp) %>% filter(dist < maxdist)
    message(paste("Classifying", nrow(temp), "unclassified pixels."))
    dat[paste(temp$x, temp$y, sep = "_"),][[as.character(enquo(var))[2]]] <- temp[[as.character(enquo(var))[2]]]
  }
  
  return(dat)
}

#### Extract sites of species #

getSites2 <- function(species, PA, type = "count"){
  cS <- PA %>% filter(rownames(.) %in% species)
  cS <- cS %>% select(which(colSums(cS) > 0))
  if(type == "count"){
    cS <- colSums(cS)
    return(cS)
  }else if(type == "percent"){
    cS <- colSums(cS)/nrow(cS)
    return(cS)
  }else{message("Type must be 'count' or 'percent'")}
}

getSites <- function(clusters, PA, type = "count"){
  out <- list()
  for(i in seq_along(clusters)){
    out[[i]] <- getSites2(clusters[[i]], PA, type)
  }
  return(out)
}


## special merge ##
multimerge <- function(x, y){
  merge(x, y, by = 0, all = TRUE) %>% namerows()
}

getGroups <- function(groups, species){
  purrr::map_int(groups, ~length(which(species %in% .)))
}


####### PLOTTING ############

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

map_mosaic <- function(df.grids){
  
  width <- df.grids %>% sapply(function(x) return(x["xmax"]-x["xmin"]))
  height <- df.grids %>% sapply(function(x) return(x["ymax"]-x["ymin"]))
  aspect <- width/height
  
  l <- length(height)
  o <- order(height)
  aspect <- aspect[o]
  height <- height[o]
  width  <-  width[o]
  
  nrows <<- (sum(aspect)/2.6) %>% round() 
  if(nrows == 0) nrows <<- 1
  if(nrows > 1) {
    newlines <- combn(x = 1:(l-1), m = nrows-1)
    
    f <- list()
    aspdiff <- numeric()
    widdiff <- numeric()
    for(i in 1:ncol(newlines)) {
      f[[i]] <- rep(letters[1:length(newlines[,i])], c(newlines[1,i], diff(newlines[,i])))
      f[[i]] <- c(f[[i]], rep("z", length(aspect)- length(f[[i]]))) %>% as.factor()
      aspdiff[i] <- split(aspect, f[[i]]) %>% sapply(sum) %>% diff() %>% abs() %>% sum()
      widdiff[i] <- split(width, f[[i]]) %>% sapply(sum) %>% diff() %>% abs() %>% sum()
    }
    fac <- f[[#((widdiff %>% abs() %>% scale()) + 
      (aspdiff %>% abs() %>% scale())#) 
      %>% which.min()]]
    #fac <- f[[(aspdiff %>% scale() %>% abs()) %>% which.min()]]
    
    # cmd <- map2(split(paste0("plot", o), fac) %>% lapply(paste, collapse = ","), 
    #     split(width, fac) %>% lapply(paste, collapse = ","), 
    #     function(x, y) paste0("plot_grid(", x, ", rel_widths = c(", y, "), nrow = 1)")) %>% 
    #  paste(collapse = ", ")
    cmd <- split(paste0("plot", o), fac) %>% lapply(paste, collapse = "+") %>% paste(" + plot_annotation(theme = theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'))) + plot_layout(nrow = 1)", collapse = ") / (")
  }else{
    #cmd <- paste0(paste0("plot", o, collapse = ",") %>% paste0(", rel_widths = c(", paste0(width, collapse = ","), ")"))
    cmd <- paste0(paste0("plot", o, collapse = "+"))
  }
  #plotcmd <- paste0("plot_grid(", cmd, ", nrow = ", nrows, ")")
  plotcmd <- paste0("(", cmd, ")")
  return(plotcmd)
}

map_mosaic2 <- function(df.grids){
  
  width <-  df.grids %>% purrr::map_dbl(~return(.x["xmax"] - .x['xmin']))
  height <-  df.grids %>% purrr::map_dbl(~return(.x["ymax"] - .x['ymin']))
  aspect <- width/height
  
  centrx <- df.grids %>% purrr::map_dbl(~return(mean(c(.x["xmax"], .x['xmin']))))
  centry <- df.grids %>% purrr::map_dbl(~return(mean(c(.x["ymax"], .x['ymin']))))
  
  wide <- which(aspect >= 1.1)
  tall <- which(aspect <= 0.9)
  square <- which(aspect > 0.9 & aspect < 1.1)
  
  if(length(wide) < length(tall)) wide <- c(wide, square) else tall <- c(tall, square)
  
  tall <- sort(centrx[tall])
  wide <- sort(centry[wide], decreasing = T)
  
  if(length(tall) < 3){tallside <- centrx[names(tall)] < 0} else { tallside <-tall < median(tall)}  ## TRUE - left; FALSE - right
  if(length(wide) < 3){wideside <- centry[names(wide)] < 0} else {wideside <- !wide > median(wide)} ## TRUE - bottom; FALSE - top
  
  top <- names(which(!wideside))[order(centrx[names(which(!wideside))])]
  bottom <- names(which(wideside))[order(centrx[ names(which(wideside))])]
  left <- names(which(tallside))[order(centry[names(which(tallside))], decreasing = T)]
  right <- names(which(!tallside))[order(centry[names(which(!tallside))], decreasing = T)]
  
  topHeight <- 0
  botHeight <- 0
  leftWidth <- 0
  rightWidth <- 0 
  
  if(length(top) > 0){
    top <-    paste("plot_grid(", paste0("insets$`", top,   "` + theme(legend.position = 'none') + labs(subtitle = '", top, "')", collapse = ","), ", nrow = 1)")
    topHeight <- 0.4
  } else top <- NA    
  if(length(bottom) > 0){
    bottom <- paste("plot_grid(", paste0("insets$`", bottom,"` + theme(legend.position = 'none') + labs(subtitle = '", bottom, "')", collapse = ","), ", nrow = 1)") 
    botHeight <- 0.4
  } else bottom <- NA  
  if(length(left) > 0){
    left <-   paste("plot_grid(", paste0("insets$`", left,  "` + theme(legend.position = 'none') + labs(subtitle = '", left, "')", collapse = ","), ", ncol = 1)") 
    leftWidth <- 0.35
  } else left <- NA    
  if(length(right) > 0){
    right <-  paste("plot_grid(", paste0("insets$`", right, "` + theme(legend.position = 'none') + labs(subtitle = '", right, "')", collapse = ","), ", ncol = 1)") 
    rightWidth <- 0.35
  }  else right <- NA
  
  top <- eval(parse(text = top))
  bottom <- eval(parse(text = bottom))
  left <- eval(parse(text = left))
  right <- eval(parse(text = right))
  
  if(length(tall)> length(wide)){
    plot_grid(left, plot_grid(top, T1, bottom, ncol = 1, rel_heights = c(topHeight, 1, botHeight)), right, nrow = 1, rel_widths = c(leftWidth, 1, rightWidth), align = "none")
  } else {
    plot_grid(top, plot_grid(left, T1, right, nrow = 1, rel_widths = c(leftWidth, 1, rightWidth)), bottom, ncol = 1, rel_heights = c(topHeight, 1, botHeight), align = "none")
  }
}

#### Mapping ######

geo_bounds <- function(df, buff = 0.1){
  bufferx <- diff(range(df$x))*buff
  buffery <- diff(range(df$y))*buff
  ext <- c(xmin = min(df$x) - bufferx, 
           xmax = max(df$x) + bufferx,
           ymin = min(df$y) - buffery,
           ymax = max(df$y) + buffery)
  return(ext)
}

df_crop <- function(df, ext, buff = 5){
  df %>% filter(x >= ext["xmin"]-buff, 
                x <= ext["xmax"] + buff, 
                y >= ext["ymin"] - buff, 
                y <= ext["ymax"]+ buff) %>% 
    return()
}

aspect_ratio <- function(df){
  rangex <- diff(range(df$x))
  rangey <- diff(range(df$y))
  
  return(rangex/rangey)
}


#### GBIF Download data 
gbif_parse <- function(x) {
  # x is a vector of species names
  library(RJSONIO)
  library(RCurl)
  library(plyr)
  u <- "http://api.gbif.org/v1/parser/name"
  res <- fromJSON(
    postForm(u,
             .opts = list(postfields = RJSONIO::toJSON(x),
                          httpheader = c('Content-Type' = 'application/json')))
  )
  do.call(rbind.fill, lapply(res, as.data.frame))  
}


download_GBIF_all_species = function (species_list, path) {
  
  ## create variables
  skip.spp.list       = list()
  GBIF.download.limit = 10000
  
  ## for every species in the list
  for(sp.n in species_list){
    
    ## 1). First, check if the file exists
    ## data\base\HIA_LIST\GBIF\SPECIES
    file_name = paste0(path,"/", sp.n, "_GBIF_records.RData")
    
    ## If it's already downloaded, skip
    if (file.exists (file_name)) {
      
      print(paste ("file exists for species", sp.n, "skipping"))
      next
      
    }
    #  create a dummy file
    dummy = data.frame()
    save (dummy, file = file_name)
    sampDat <- occ_search(scientificName = sp.n, limit = 1)
    
    ## 2). Then check the spelling...incorrect nomenclature will return NULL result
    if (is.null(sampDat$meta$count)) {
      
      ## now append the species which had incorrect nomenclature to the skipped list
      ## this is slow, but it works for now
      print (paste ("Possible incorrect nomenclature", sp.n, "skipping"))
      nomenclature = paste ("Possible incorrect nomenclature |", sp.n)
      skip.spp.list <- c(skip.spp.list, nomenclature)
      next
      
    }
    
    ## 3). Skip species with no records
    if (sampDat$meta$count < 1) {
      
      ## now append the species which had no records to the skipped list
      print (paste ("No GBIF records for", sp.n, "skipping"))
      records = paste ("No GBIF records |", sp.n)
      skip.spp.list <- c(skip.spp.list, records)
      next
      
    }
    
    ## 4). Check how many records there are, and skip if there are over 200k
    if (sampDat$meta$count > GBIF.download.limit) #{
      
      ## now append the species which had > 200k records to the skipped list
      print (paste ("Number of records > max for GBIF download via R (100,000)", sp.n))
      #max =  paste ("Number of records > 200,000 |", sp.n)
      
      
      ## and send a request to GBIF for download
      # message("Sending request to GBIF to download ", sp.n, " using rgbif :: occ_download")
      # key  <- name_backbone(name = sp.n, rank = 'species')$usageKey
      # GBIF = occ_download(paste('taxonKey = ', key),  user = "popple_1500", pwd = "Popple1500", email = "hugh.burley@mq.edu.au")
      # save(GBIF, file = paste(path, sp.n, "_GBIF_request.RData", sep = ""))
      # skip.spp.list <- c(skip.spp.list, max)
      
 #   } else {
      
      ## 5). Download ALL records from GBIF
      message("Downloading GBIF records for ", sp.n, " using rgbif :: occ_data")
      key <- name_backbone(name = sp.n, rank = 'species')$usageKey
      # x = name_lookup(sp.n)
      # keys = x$data$key
      # 
      # GBIF <- keys %>%         
      #   
      #   ## pipe the list into lapply
      #   lapply(function(x) {
      #     
      #     ## Create the character string
      #     f <- occ_data(taxonKey = x, limit = GBIF.download.limit)
      #     f =  as.data.frame(f$data)
      #     ## Load each .RData file
      #     
      #   }) %>%
      #   
      #   ## Finally, bind all the rows together
      #   dplyr::bind_rows
      
      GBIF <- occ_data(taxonKey = key, limit = GBIF.download.limit)
      GBIF <- as.data.frame(GBIF$data)
      
      cat("Synonyms returned for :: ",  sp.n, unique(GBIF$scientificName), sep="\n")
      cat("Names returned for :: ", sp.n, unique(GBIF$name),               sep="\n")
      cat("Takonkeys returned for :: ", sp.n, unique(GBIF$taxonKey),       sep="\n")
      
      ## Could also only use the key searched, but that could knock out a lot of species
      #GBIF = GBIF[GBIF$taxonKey %in% key, ]
      #View(GBIF[c("name", "scientificName", "taxonKey")])
      
      message(dim(GBIF)[1], " Records returned for ", sp.n)
      
      ## 6). save records to .Rdata file, note that using .csv files seemed to cause problems...
      #save(GBIF, file = paste(path, sp.n, "_GBIF_records.RData", sep = ""))
      save(GBIF, file = file_name)
      
 #   }
    
  }
  
  return(skip.spp.list)
  
}

load.GBIF <- function(x, GBIF_path) {
  
  ## Create a character string of each .RData file
  f <- sprintf(paste0(GBIF_path, "/%s"), x)
  ## Load each file
  d <- get(load(f))
  ## Now drop the columns which we don't need
  message ('Reading GBIF data for ', x)
  ## Check if the dataframes have data
  if (nrow(d) <= 2) {
    
    ## If the species has < 2 records, escape the loop
    print (paste ("No GBIF records for ", x, " skipping "))
    return (d)
    
  }
  dat <- data.frame(searchTaxon = x, d[, colnames(d) %in% gbif.keep],
                    stringsAsFactors = FALSE)
  
  if(!is.character(dat$gbifID)) {
    
    dat$gbifID <- as.character(dat$gbifID)
    
  }
  
  ## Need to print the object within the loop
  names(dat)[names(dat) == 'decimalLatitude']  <- 'lat'
  names(dat)[names(dat) == 'decimalLongitude'] <- 'lon'
  dat$searchTaxon = gsub("_GBIF_records.RData", "", dat$searchTaxon)
  return(dat)
  
}

get_gbif_taxonomy <- function (x, subspecies = TRUE, higherrank = FALSE, verbose = FALSE, 
                               fuzzy = TRUE, conf_threshold = 90, resolve_synonyms = TRUE) 
{
  matchtype = status = confidence = NULL
  temp <- taxize::get_gbifid_(x, messages = verbose)
  for (i in 1:length(temp)) {
    warning_i = ""
    synonym_i = FALSE
    if (nrow(temp[[i]]) == 0) {
      warning_i <- paste("No matching species concept!")
      temp[[i]] <- data.frame(scientificName = x[i], matchtype = "NONE", 
                              status = "NA", rank = "species")
    }
    if (!fuzzy & nrow(temp[[i]]) > 0) {
      temp[[i]] <- subset(temp[[i]], matchtype != "FUZZY")
      if (nrow(temp[[i]]) == 0) {
        warning_i <- paste(warning_i, "Fuzzy matching might yield results.")
      }
    }
    if (!is.null(conf_threshold) & nrow(temp[[i]]) > 0) {
      temp[[i]] <- subset(temp[[i]], confidence >= conf_threshold)
      if (nrow(temp[[i]]) == 0) {
        temp[[i]] <- data.frame(scientificName = x[i], 
                                matchtype = "NONE", status = "NA", rank = "species")
        warning_i <- paste(warning_i, "No match! Check spelling or lower confidence threshold!")
      }
    }
    if (any(temp[[i]]$status == "ACCEPTED")) {
      temp[[i]] <- subset(temp[[i]], status == "ACCEPTED")
      temp[[i]] <- subset(temp[[i]], temp[[i]]$confidence == 
                            max(temp[[i]]$confidence))
      if (nrow(temp[[i]]) > 1) {
        temp[[i]] <- temp[[i]][1, ]
        warning_i <- paste(warning_i, "Selected first of multiple equally ranked concepts!")
      }
    }
    if (!any(temp[[i]]$status == "ACCEPTED") & any(temp[[i]]$status == 
                                                   "SYNONYM")) {
      if (resolve_synonyms) {
        keep <- temp[i]
        temp[i] <- taxize::get_gbifid_(temp[[i]]$species[which.max(temp[[i]]$confidence)], 
                                       messages = verbose)
        if (temp[[i]][1, ]$status == "ACCEPTED") {
          temp[[i]] <- subset(temp[[i]], matchtype == 
                                "EXACT" & status == "ACCEPTED")
          temp[[i]] <- subset(temp[[i]], temp[[i]]$confidence == 
                                max(temp[[i]]$confidence))
          if (nrow(temp[[i]]) > 1) {
            temp[[i]] <- temp[[i]][1, ]
            warning_i <- paste(warning_i, "Selected first of multiple equally ranked concepts!")
          }
          warning_i <- paste(warning_i, "A synonym was mapped to the accepted species concept!", 
                             sep = " ")
          synonym_i = TRUE
        }
        else {
          status <- temp[[i]][1, ]$status
          temp[i] <- keep
          if (nrow(temp[[i]]) > 1) {
            temp[[i]] <- temp[[i]][1, ]
            warning_i <- paste(warning_i, "Selected first of multiple equally ranked concepts!")
          }
          warning_i <- paste0(warning_i, " Resolved synonym '", 
                              temp[[i]]$species, "' is labelled '", status, 
                              "'. Clarification required!")
        }
      }
      else {
        temp[[i]] <- subset(temp[[i]], status == "SYNONYM")
        temp[[i]] <- subset(temp[[i]], temp[[i]]$confidence == 
                              max(temp[[i]]$confidence))
        warning_i <- paste(warning_i, "The provided taxon seems to be a synonym of '", 
                           temp[[i]]$species, "'!", sep = "")
      }
    }
    if (all(temp[[i]]$status == "DOUBTFUL")) {
      temp[[i]] <- subset(temp[[i]], status == "DOUBTFUL")
      warning_i <- paste(warning_i, "Mapped concept is labelled 'DOUBTFUL'!")
      temp[[i]] <- subset(temp[[i]], temp[[i]]$confidence == 
                            max(temp[[i]]$confidence))
      if (nrow(temp[[i]]) > 1) {
        temp[[i]] <- temp[[i]][1, ]
        warning_i <- paste(warning_i, "Selected first of multiple equally ranked concepts!")
      }
    }
    rankorder <- c("kingdom", "phylum", "class", "order", 
                   "family", "genus", "species", "subspecies")
    if (match(temp[[i]]$rank, rankorder) > 7 & !subspecies) {
      if (length(strsplit(as.character(temp[[i]]$canonicalname), 
                          " ")[[1]]) > 2) {
        temp[i] <- taxize::get_gbifid_(paste(strsplit(names(temp[i]), 
                                                      " ")[[1]][1:2], collapse = " "), messages = verbose)
        temp[[i]] <- subset(temp[[i]], temp[[i]]$confidence == 
                              max(temp[[i]]$confidence))
        warning_i <- paste(warning_i, "Subspecies has been remapped to species concept!", 
                           sep = " ")
      }
      else {
        temp[[i]] <- data.frame(scientificName = x[i], 
                                matchtype = "NONE", rank = "subspecies")
        warning_i <- paste(warning_i, "No mapping of subspecies name to species was possible!", 
                           sep = " ")
      }
    }
    if (temp[[i]]$matchtype == "HIGHERRANK") {
      if (higherrank) {
        temp[[i]] <- subset(temp[[i]], temp[[i]]$confidence == 
                              max(temp[[i]]$confidence))
        warning_i <- paste(warning_i, "No matching species concept! Entry has been mapped to higher taxonomic level.")
      }
      else {
        temp[[i]] <- data.frame(scientificName = x[i], 
                                matchtype = "NONE", rank = "highertaxon")
        warning_i <- paste("No matching species concept!", 
                           warning_i)
      }
    }
    if (temp[[i]]$matchtype != "NONE") {
      temp[[i]] <- data.frame(scientificName = x[i], synonym = synonym_i, 
                              scientificNameStd = temp[[i]]$canonicalname, 
                              matchtype = temp[[i]]$matchtype,
                              author = sub(paste0(temp[[i]]$canonicalname, 
                                                  " "), "", temp[[i]]$scientificname), taxonRank = temp[[i]]$rank, 
                              confidence = temp[[i]]$confidence, kingdom = if (is.null(temp[[i]]$kingdom)) 
                                NA
                              else temp[[i]]$kingdom, phylum = if (is.null(temp[[i]]$phylum)) 
                                NA
                              else temp[[i]]$phylum, class = if (is.null(temp[[i]]$class)) 
                                NA
                              else temp[[i]]$class, order = if (is.null(temp[[i]]$order)) 
                                NA
                              else temp[[i]]$order, family = if (is.null(temp[[i]]$family)) 
                                NA
                              else temp[[i]]$family, genus = if (is.null(temp[[i]]$genus)) 
                                NA
                              else temp[[i]]$genus, taxonomy = "GBIF Backbone Taxonomy", 
                              taxonID = paste0("http://www.gbif.org/species/", 
                                               temp[[i]]$usagekey, ""), warnings = NA)
    }
    else {
      temp[[i]] <- data.frame(scientificName = x[i], warnings = NA)
    }
    temp[[i]]$warnings <- warning_i
    if (verbose & nchar(warning_i) >= 1) 
      warning(warning_i)
  }
  out <- data.table::rbindlist(temp, fill = TRUE)
  class(out) <- c("data.frame", "taxonomy")
  return(out)
}


gbif.keep <- c(## TAXONOMY
  "searchTaxon",
  "species",
  "scientificName",
  "taxonRank",
  "taxonKey",
  "genus",
  "family",
  
  ## CULTIVATION
  "cloc",
  "basisOfRecord",
  "locality",
  "establishmentMeans",
  "institutionCode",
  "datasetName",
  "habitat",
  "eventRemarks",
  
  ## RECORD ID
  "recordedBy",
  "identifiedBy",
  "gbifID",
  "catalogNumber",
  
  ## PLACE/TIME
  "lat",
  "lon",
  "decimalLatitude",
  "decimalLongitude",
  "country",
  "coordinateUncertaintyInMeters",
  "geodeticDatum",
  "year",
  "month",
  "day",
  "eventDate",
  "eventID")


#####
#####
####