library(ggplot2)
library(viridis)
library(dplyr)

## Summary of unit areas distribution in ACBRs
  # ovl_area <- read_csv("../Data/Base/IFA_Unit_areas.txt")
  # area_summary <- ovl_area %>% group_by(gridcode, ACBR_Name) %>% summarise(area = sum(Shape_Area)/1000000)
  # area_summary <- data.frame(gridcode = unique(area_summary$gridcode), env = rep(1:6, times = c(5,4,5,5,5,7))) %>% 
  #   merge(area_summary, by = "gridcode")
  # ggplot(area_summary, aes(x = gridcode, fill = as.factor(env), y = area)) + geom_col() + 
  #   facet_wrap(~ACBR_Name, scales = "free") + scale_fill_viridis(discrete = TRUE) + labs(fill = "Env group") + 
  #   theme(legend.position = c(0.85, 0.1))

## Summary of unit areas distribution in ACBRs V5
ovl_area <- read_csv("../Data/Base/IFA_Unit_areas_V5.txt") # generated using TabulateAreas() in arcGIS
ACBR_key <- readRDS("./data/ACBR_key.rds")
area_summary <- ovl_area %>% pivot_longer(contains("VALUE_"), names_to = "ACBR", values_to = "area") %>% 
  select(-Rowid_) %>% mutate(ACBR = as.numeric(word(ACBR, 2,2, "_")), 
                             area = area/1000000) %>%
  left_join(ACBR_key, by = c("ACBR" = "ACBR_ID")) %>% 
  left_join(out %>% pull(unit_h) %>% unique() %>% sort() %>% 
              data.frame(VALUE = 1:33) %>% setNames(c("unit", "VALUE")), by = "VALUE") %>%
  mutate(env = word(unit, 1, 1, sep = "_"))

 ggplot(area_summary, aes(x = VALUE, fill = as.factor(env), y = area)) + geom_col() + 
   facet_wrap(~ACBR_Name, scales = "free") + scale_fill_viridis(discrete = TRUE) + 
   labs(fill = "Environmental group", y = "Area (square km)", x = "Ecosystem") + 
   theme(legend.position = c(0.85, 0.1))



## summary of unit areas overall
area_summary %>% filter(!grepl("sdmNA", unit, fixed = TRUE)) %>% 
  group_by(env, VALUE) %>% summarise(area = sum(area)) %>%
  ggplot(aes(x = VALUE, fill = as.factor(env), y = area)) + geom_col() + 
    scale_fill_viridis(discrete = TRUE) + 
    labs(fill = "Group", y = "Area (square km)", x = "Ecosystem")


## endemics 
PA01 <- tobinary.single(PA)
endemics <- merge(area_summary, PA01[which(rowSums(PA01) ==1),] %>% colSums(), by.x = "gridcode", by.y = 0) %>% rename(endemics = y)
ggplot(endemics, aes(x = gridcode, fill = as.factor(env), y = endemics)) + geom_col() + scale_fill_viridis(discrete = TRUE) + labs(fill = "Env group") 
