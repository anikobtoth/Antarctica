
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
#####
#####
####