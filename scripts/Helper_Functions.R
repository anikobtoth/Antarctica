
library(sp)
library(rlang)

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

# remove empty rows and columns in 1 matrix
# or remove rows and columns with too few observations
clean.empty <- function(x, mincol = 1, minrow = 1){
  x <- x[which(rowSums(x) > minrow-1),]
  x <- x[,which(colSums(x) > mincol-1)]
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


##### ANALYSES ######
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

### Forbes index functions by J. Alroy ###
##### Available from http://bio.mq.edu.au/~jalroy/Forbes.R ###
##### help page http://bio.mq.edu.au/~jalroy/Forbes.html ###

forbes<-function(x,y,corrected)	{
  if (missing(corrected))	{
    corrected <- T
  }
  if (is.numeric(x) && is.numeric(y) && min(x) == 0 && min(y) == 0 && length(x) == length(y))	{
    a <- length(which((x * y) > 0))
    b <- length(which(x > 0)) - a
    c <- length(which(y > 0)) - a
  } else	{
    a <- length(na.omit(match(x,y)))
    b <- length(x) - a
    c <- length(y) - a
  }
  n <- a + b + c
  if (corrected == T)	{
    return(a * (n + sqrt(n))/((a + b) * (a + c) + a * sqrt(n) + (b * c)/2))
  } else	{
    return(a * n/((a + b) * (a + c)))
  }
}

forbesMatrix<-function(x,corrected)	{
  x[is.na(x)] <- 0
  m = matrix(nrow=ncol(x),ncol=ncol(x))
  for (i in 1:ncol(x))	{
    for (j in 1:ncol(x))	{
      m[i,j] = forbes(x[,i],x[,j],corrected)
    }
  }
  m <- as.data.frame(m)
  rownames(m) <- colnames(x)
  colnames(m) <- colnames(x)
  return(m)
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
