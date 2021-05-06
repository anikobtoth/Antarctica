
## Summary of unit areas distribution in ACBRs
ovl_area <- read_csv("../Data/Base/IFA_Unit_areas.txt")
area_summary <- ovl_area %>% group_by(gridcode, ACBR_Name) %>% summarise(area = sum(Shape_Area)/1000000)
area_summary <- data.frame(gridcode = unique(area_summary$gridcode), env = rep(1:6, times = c(5,4,5,5,5,7))) %>% 
  merge(area_summary, by = "gridcode")
ggplot(area_summary, aes(x = gridcode, fill = as.factor(env), y = area)) + geom_col() + 
  facet_wrap(~ACBR_Name, scales = "free") + scale_fill_viridis(discrete = TRUE) + labs(fill = "Env group") + 
  theme(legend.position = c(0.85, 0.1))

## summary of unit areas overall
area_summary <- area_summary %>% group_by(env, gridcode) %>% summarise(area = sum(area))
ggplot(area_summary, aes(x = gridcode, fill = as.factor(env), y = area)) + geom_col() + scale_fill_viridis(discrete = TRUE) + labs(fill = "Env group") 


## endemics 
PA01 <- tobinary.single(PA)
endemics <- merge(area_summary, PA01[which(rowSums(PA01) ==1),] %>% colSums(), by.x = "gridcode", by.y = 0) %>% rename(endemics = y)
ggplot(endemics, aes(x = gridcode, fill = as.factor(env), y = endemics)) + geom_col() + scale_fill_viridis(discrete = TRUE) + labs(fill = "Env group") 
