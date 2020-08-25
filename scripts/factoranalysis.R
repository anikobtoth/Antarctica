
library(nFactors)
library(psych)

cldat <- v1[,7:ncol(v1)] %>% na.omit() %>% dplyr::select(-geoT, -all_of(bad_models)) %>% scale()

test <- fa(scale(cldat), 7, rotate = "varimax" )

sc <- data.frame(test$scores)

pdat <- data.frame(sc, na.omit(v1[,1:6]))

plot(y~x, data = pdat, col = hsv(h = pnorm(scale(pdat$MR6))), 
     pch = 16, cex = pnorm(scale(pdat$MR6))^6)


pdat$cluster <- apply(sc, 1, which.max)

antarctica <- readOGR("../Data/Base", "Antarctic_mainland")
points(y~x, data = pdat, col = cluster, 
     pch = 16, cex = 0.1)

cv <- lapply(2:24, function(x) fa(scale(cldat), x, rotate = "varimax")$Vaccounted["Cumulative Var",]) %>%
  sapply(last)
