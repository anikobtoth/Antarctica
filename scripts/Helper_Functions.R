
library(sp)
library(rlang)
library(gmp)
library(BBmisc)
library(psych)
library(fields)

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

## Contingency table 
cont_table <- function(x){ 
  samples = ncol(x)  #S
  a = matrix(nrow=nrow(x),ncol=nrow(x),data=0)
  occs = array()
  
  #Calculate overlap
  for (i in 2:nrow(x))  {
    for (j in 1:(i-1))
    {
      a[i,j] = length(which(x[i,] > 0 & x[j,] > 0)) # B
    }
  }
  
  a <- as.dist(a, diag = F, upper = F)
  
  l <- dist2edgelist(a, x)
  s <- rowSums(x) %>% data.frame()
  
  t <- merge(l, s, by.x = "Sp1", by.y = 0)
  t <- merge(t, s, by.x = "Sp2", by.y = 0)
  t <- t %>% select(Sp1, Sp2, ..x, ..y, Z.Score)
  t$samples <- ncol(x)
  names(t) <- c("Sp1", "Sp2", "presSp1", "presSp2", "presBoth", "samples")
  t$absSp1 <- t$samples - t$presSp1
  t$absSp2 <- t$samples - t$presSp2
  
  return(t)
}

##### ANALYSES ######
#wrapper for fa() that chooses factor number based on importance of factor loadings
factor_analysis <- function(dat, mincomp = 0.35){
  message("Choosing number of factors")
  dat <- dat %>% na.omit()
  rn <- rownames(dat)
  dat <- sapply(dat, function(x) normalize(x, method = "range", range = c(min(x[x > 0])/10, 1-(min(x[x > 0])/10)), margin = 2)) %>% data.frame()
  dat <- sapply(dat, qnorm) %>% scale()
  
  cvs <- purrr::map(2:ncol(dat), function(x) fa(dat, nfactors = x, rotate = "varimax")$loadings %>% as.matrix() %>% abs() %>% apply(2, max)) %>% sapply(min)
  nfact <<- 1+ length(which(cvs >=mincomp))
  
  message(paste("Running factor analysis with", nfact, "factors"))  
  fa1 <- fa(dat, nfact, rotate = "varimax" )
  sc <- data.frame(fa1$scores)
  
  consensus <- apply(sc, 1, which.max) %>% factor() %>% setNames(rn)
  
  return(consensus)
}

# classifies unclassified pixels in column "var" with nearest neighbour
classify_by_neighbours <- function(dat, var, maxdist = 1.5){
  library(fields)
 
  v0 <- dat %>% filter(is.na({{ var }})) 
  v2 <- dat %>% filter(!is.na({{ var }}))
  
  if(nrow(v0) > 0){
    v2 <- v2 %>% filter(x %in% unique(v0$x, v0$x+1000, v0$x - 1000) & y %in% unique(v0$y, v0$y+1000, v0$y - 1000))
    temp <- v2[rdist.earth(v0[,c("lon", "lat")], v2[,c("lon", "lat")], miles = F) %>% apply(1, which.min),] %>% 
      select(x, y, lon, lat, {{var}}) %>% 
      mutate(dist = rdist.earth(v0[,c("lon", "lat")], v2[,c("lon", "lat")], miles = F) %>% apply(1, min)) %>% select({{var}}, dist)
    temp <- data.frame(x = v0$x, y = v0$y, temp) %>% filter(dist < maxdist)
    message(paste("Classifying", nrow(temp), "unclassified pixels."))
    dat[paste(temp$x, temp$y, sep = "_"),][[as.character(enquo(var))[2]]] <- temp[[as.character(enquo(var))[2]]]
    }
  
  return(dat)
}

# FETmP
simpairs <- function(x){ #simpairs function, simpairs only out
  samples = ncol(x)  #S
  z = matrix(nrow=nrow(x),ncol=nrow(x),data=0)
  occs = array()
  
  #convert to P/A. Occs = rowsums of PA matrix.
  x <- x/x
  x[is.na(x)] <- 0
  occs <- rowSums(x)
  
  #SimPairs Algorithm
  for (i in 2:nrow(x))  {
    for (j in 1:(i-1))
    {
      a = length(which(x[i,] > 0 & x[j,] > 0)) # B
      
      #simpairs
      for (k in 0:a)
        z[i,j] = z[i,j] + choose(occs[j] , k) * choose(samples - occs[j] , occs[i] - k) / choose(samples , occs[i])
      z[i,j] = z[i,j] - choose(occs[j] , a) * choose(samples - occs[j] , occs[i] - a) / choose(samples , occs[i]) / 2
      if(z[i,j] == 1) z[i,j] <- 0.99999999999999994 # solve rounding issue.
      z[i,j] = qnorm(z[i,j])
      z[j,i] = z[i,j]
    }
  }
  print("check")
  return(as.dist(z, diag = F, upper = F))
}


simpairs_lgnum <- function(x){ #simpairs function, simpairs only out
  samples = ncol(x) %>% as.bigz()  #S
  z = matrix(nrow=nrow(x),ncol=nrow(x),data=0)
  occs = array()
  
  #convert to P/A. Occs = rowsums of PA matrix.
  x <- x/x
  x[is.na(x)] <- 0
  occs <- rowSums(x)
  
  #SimPairs Algorithm
  for (i in 2:nrow(x))  {
    for (j in 1:(i-1))
    {
      a = length(which(x[i,] > 0 & x[j,] > 0)) # B
      
      occsi <- as.bigz(occs[i])
      occsj <- as.bigz(occs[j])
      ovl <- occs[i] + occs[j] - ncol(x)
      
      #simpairs
      for (k in max(0, ovl):a){
       
        k <- as.bigz(k)
        
        A <- factorial(occsj)/(factorial(k)*factorial(occsj-k))
        B <- factorial(samples - occsj)/(factorial(occsi - k)*factorial(samples - occsj - occsi + k))
        C <- factorial(samples)/(factorial(occsi)*factorial(samples-occsi))
        p <- exp(log.bigz(A)+log.bigz(B)-log.bigz(C))
        z[i,j] <- z[i,j]+ p
      }
      z[i,j] <- z[i,j]-p/2
      
      z[j,i] = z[i,j]
    }
  }
  print("check")
  return(as.dist(z, diag = F, upper = F))
}

FETmP <- function(contable){
  #cl <- makeCluster(detectCores()-2)
  #clusterExport(cl, c("contable", "FETmP_"))
  #clusterEvalQ(cl, library(tidyverse))
  #out <- parRapply(cl, contable, function(x) {x <- as.numeric(x)
  #return(FETmP_(x[3], x[4], x[5], x[6]))})
  #stopCluster(cl)
 
  out <- apply(contable, 1, function(x) {x <- as.numeric(x)
 return(FETmP_(x[3], x[4], x[5], x[6]))})
  return(out)

}

FETmP_ <- function(presSp1, presSp2, presBoth, samples){
  absSp2 <- samples-presSp2
  minovl <- max(presSp1+presSp2-samples,0)
  p <- choose(presSp2, minovl:presBoth) * choose(absSp2, presSp1-minovl:presBoth)/choose(samples, presSp1)
  
  return(sum(p)-0.5*last(p))
}


cmeans_multi <- function(dat, weights=1, centres = 8, iter = 1000, reps = 100){
  cl <- makeCluster(detectCores()-2)
  clusterExport(cl, c("dat", "weights", "centres", "iter"), envir=environment())
  clusterEvalQ(cl, library(e1071))
  
  sdm_hcl <- parLapply(cl, 1:reps, function(x) cmeans(dat, centers = centres, iter.max = iter, 
                                                     verbose = FALSE, dist = "euclidean", method = "ufcl", 
                                                     m = 2, rate.par = 0.3, weights = weights$weight))
  stopCluster(cl)
  return(sdm_hcl)
}

network_analysis <- function(x, threshold = 0.8){
  pairs <- simpairs(x)
  el <- dist2edgelist(pairs, x)
  g <- el %>% #dplyr::filter(Z.Score > quantile(Z.Score, threshold, na.rm = T)) %>% 
    graph_from_data_frame(directed = F)
  g <- delete_edges(g, E(g)[which(E(g)$Z.Score < quantile(E(g)$Z.Score, threshold, na.rm = T))])
  clust <- cluster_fast_greedy(g) #weights = E(g)$Z.Score)
  plot(g, vertex.label = NA, vertex.size = 6, vertex.color = clust$membership)
  return(list(g, clust))
}

cont_table <- function(x){ #simpairs function, simpairs only out
  samples = ncol(x)  #S
  a = matrix(nrow=nrow(x),ncol=nrow(x),data=0)
  occs = array()
  
  #Calculate overlap
  for (i in 2:nrow(x))  {
    for (j in 1:(i-1))
    {
      a[i,j] = length(which(x[i,] > 0 & x[j,] > 0)) # B
    }
  }
  
  a <- as.dist(a, diag = F, upper = F)
  
  l <- dist2edgelist(a, x)
  s <- rowSums(x) %>% data.frame()
  
  t <- merge(l, s, by.x = "Sp1", by.y = 0)
  t <- merge(t, s, by.x = "Sp2", by.y = 0)
  t <- t %>% select(Sp1, Sp2, ..x, ..y, Z.Score)
  t$samples <- ncol(x)
  names(t) <- c("Sp1", "Sp2", "presSp1", "presSp2", "presBoth", "samples")
  t$absSp1 <- t$samples - t$presSp1
  t$absSp2 <- t$samples - t$presSp2
  
  #t$presSp1absSp2 <- t$presSp1- t$presBoth
  #t$presSp2absSp1 <- t$presSp2- t$presBoth
  
  #t$absBoth <- t$samples - t$presBoth - t$presSp1absSp2 - t$presSp2absSp1
  
  return(t)
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

plotnet_ <- function(g, title) {
  pal <- brewer.pal(length(unique(V(g)$cluster)),"Accent")
  plot.igraph(g, vertex.label = NA, vertex.size = 5, vertex.color = pal[V(g)$cluster], main = title)
  legend("topleft", legend = levels(as.factor(V(g)$cluster)),pt.cex = 2, fill = pal)
  box()
}

plotnet <- function(g, mod, title){
  V(g)$cluster <- as.factor(mod$membership)
  plotnet_(g, title)
}


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#### Mapping ######

geo_bounds <- function(df, buff = 0.1){
  bufferx <- diff(range(df$x))*buff
  buffery <- diff(range(df$y))*buff
  ext <- c(left = min(df$x) - bufferx, 
             right = max(df$x) + bufferx,
             bottom = min(df$y) - buffery,
             top = max(df$y) + buffery)
  return(ext)
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
  GBIF.download.limit = 100000
  
  ## for every species in the list
  for(sp.n in species_list){
    
    ## 1). First, check if the f*&%$*# file exists
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
    if (sampDat$meta$count > GBIF.download.limit) {
      
      ## now append the species which had > 200k records to the skipped list
      print (paste ("Number of records > max for GBIF download via R (100,000)", sp.n))
      max =  paste ("Number of records > 200,000 |", sp.n)
      
      
      ## and send a request to GBIF for download
      # message("Sending request to GBIF to download ", sp.n, " using rgbif :: occ_download")
      # key  <- name_backbone(name = sp.n, rank = 'species')$usageKey
      # GBIF = occ_download(paste('taxonKey = ', key),  user = "popple_1500", pwd = "Popple1500", email = "hugh.burley@mq.edu.au")
      # save(GBIF, file = paste(path, sp.n, "_GBIF_request.RData", sep = ""))
      # skip.spp.list <- c(skip.spp.list, max)
      
    } else {
      
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
      
    }
    
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