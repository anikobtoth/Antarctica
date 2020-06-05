## This script classifies Antarctic ice-free area into habitats using fuzzy 
#c-means clustering.
#The input layers are on a 1km resolution (bedmap grid) and were compiled into
# a matrix for use here (see: "Habitat_matrix_updated_Dec18.R")

library(e1071)
library(data.table)
library(tidyverse)
library(ggfortify)

## Load data table required for appending to output at end (i.e. lat/lon,
# # of species, # of genera, # of records)
req_var = fread("./data/Habitats/req_var.csv")

## Load data table of environmental predictors with the IF cells only, 
#with all na values removed and normalised on a scale of 0 to 1
envpred_norm = fread("./data/Habitats/envpred_norm.csv")

# read in habitat data (ACBR join)
habDat <- read_csv("data/Habitats/Hab8_ACBR_join.csv")

########### PCA analysis ############################

x <- princomp(envpred_norm[,2:15] %>% select(-sumtemp, -rugos, -rad, -melt, -modT_0315, -DDminus5))

autoplot(x, data = habDat, col = "n8U", 
         loadings = TRUE, loadings.label = TRUE, 
         loadings.label.size = 5, alpha = 0.4)


library(rgl)

habDat$n8U <- cmeans(envpred_norm %>% select(-sumtemp, -rugos, -rad, -melt, -modT_0315, -DDminus5), 
       8, iter.max = 1000,verbose = FALSE, 
       dist = "euclidean", method = "ufcl", m = 2, 
       rate.par = 0.3, weights = 1)$cluster %>% as.factor()

plot3d(x$scores[,1:3], col=habDat$n8U)
text3d(x$loadings[,1:3], texts=rownames(x$loadings), col="red")

########### UNWEIGHTED UNSUPERVISED FUZZY c-MEANS ########################

##NOTE: All parameters included up to this point, though may want to
#consider removing some (such as either temp or modT_01315) if they are 
#highly correlated (otherwise they may end up with more emphasis in clustering)

### 6 groups Unweighted 1000 iterations
UFCM_n6_U = cmeans(envpred_norm %>% select(-sumtemp, -rugos, -rad, -melt, -modT_0315, -DDminus5), 
                   6, iter.max = 1000,verbose = FALSE, 
                   dist = "euclidean", method = "ufcl", m = 2, 
                   rate.par = 0.3, weights = 1)

UFCM_n6_U_DT =req_var
UFCM_n6_U_DT <- cbind(UFCM_n6_U_DT, UFCM_n6_U$membership)
UFCM_n6_U_DT$n6_U <-UFCM_n6_U$cluster
write.csv(UFCM_n6_U_DT, file="UFCM_n6_U.csv") 

### 6 groups Weigted 1000 iterations
UFCM_n6_W = cmeans(envpred_norm, 6, iter.max = 1000,verbose = FALSE, 
                   dist = "euclidean", method = "ufcl", m = 2, 
                   rate.par = 0.3, weights = req_var$rck01_prop)

UFCM_n6_W_DT =req_var
UFCM_n6_W_DT <-cbind(UFCM_n6_W_DT, UFCM_n6_W$membership)
UFCM_n6_W_DT$n6_W <-UFCM_n6_W$cluster
write.csv(UFCM_n6_W_DT, file="UFCM_n6_W.csv")

### 8 groups Unweighted 1000 iterations
UFCM_n8_U = cmeans(envpred_norm, 8, iter.max = 1000,verbose = FALSE, 
                   dist = "euclidean", method = "ufcl", m = 2, 
                   rate.par = 0.3, weights = 1)

UFCM_n8_U_DT =req_var
UFCM_n8_U_DT<- cbind(UFCM_n8_U_DT, UFCM_n8_U$membership)
UFCM_n8_U_DT$n8_U <-UFCM_n8_U$cluster
write.csv(UFCM_n8_U_DT, file="UFCM_n8_U.csv")

### 8 groups Weighted 1000 iterations
UFCM_n8_W = cmeans(epred_norm_narm, 8, iter.max = 1000,verbose = FALSE, 
                   dist = "euclidean", method = "ufcl", m = 2, 
                   rate.par = 0.3, weights = req_var$rck01_prop)

UFCM_n8_W_DT =req_var
UFCM_n8_W_DT$member <-UFCM_n8_W$membership
UFCM_n8_W_DT$n8_W <-UFCM_n8_W$cluster
write.csv(UFCM_n8_W_DT, file="UFCM_n8_W.csv")

### 10 groups Unweighted 1000 iterations
UFCM_n10_U = cmeans(epred_norm_narm, 10, iter.max = 1000,verbose = FALSE, 
                   dist = "euclidean", method = "ufcl", m = 2, 
                   rate.par = 0.3, weights = 1)

UFCM_n10_U_DT =req_var
UFCM_n10_U_DT$member <-UFCM_n10_U$membership
UFCM_n10_U_DT$n10_U <-UFCM_n10_U$cluster
write.csv(UFCM_n10_U_DT, file="UFCM_n10_U.csv")

### 10 groups Weighted 1000 iterations
UFCM_n10_W = cmeans(epred_norm_narm, 10, iter.max = 1000,verbose = FALSE, 
                   dist = "euclidean", method = "ufcl", m = 2, 
                   rate.par = 0.3, weights = req_var$rck01_prop)

UFCM_n10_W_DT =req_var
UFCM_n10_W_DT$member <-UFCM_n10_W$membership
UFCM_n10_W_DT$n10_W <-UFCM_n10_W$cluster
write.csv(UFCM_n10_W_DT, file="UFCM_n10_W.csv")

#NOTE: still hasent converged before 1000 iterations, so may wish to
#consider using more
