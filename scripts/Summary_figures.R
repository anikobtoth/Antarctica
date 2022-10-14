library(ggplot2)
library(viridis)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

## Summary of unit areas distribution in ACBRs
  # ovl_area <- read_csv("../Data/Base/IFA_Unit_areas.txt")
  # area_summary <- ovl_area %>% group_by(gridcode, ACBR_Name) %>% summarise(area = sum(Shape_Area)/1000000)
  # area_summary <- data.frame(gridcode = unique(area_summary$gridcode), env = rep(1:6, times = c(5,4,5,5,5,7))) %>% 
  #   merge(area_summary, by = "gridcode")
  # ggplot(area_summary, aes(x = gridcode, fill = as.factor(env), y = area)) + geom_col() + 
  #   facet_wrap(~ACBR_Name, scales = "free") + scale_fill_viridis(discrete = TRUE) + labs(fill = "Env group") + 
  #   theme(legend.position = c(0.85, 0.1))

## Summary of unit areas distribution in ACBRs V5
ovl_area <- read_csv("../Data/Base/typV6_acbr_unit_areas.txt") # generated using TabulateAreas() in arcGIS
ACBR_key <- readRDS("./data/ACBR_key.rds")
area_summary <- ovl_area %>% pivot_longer(contains("VALUE_"), names_to = "biotic_assemblage", values_to = "area") %>% 
  dplyr::select(-Rowid_) %>% mutate(biotic_assemblage = as.numeric(word(biotic_assemblage, 2,2, "_")), 
                             area = area/1000000) %>%
  left_join(ACBR_key, by = c("VALUE" = "ACBR_ID")) %>% 
  left_join(typkey, by = c("biotic_assemblage")) %>%
  mutate(plot_unit = ifelse(biotic_assemblage < 100, unit, 
                            ifelse(biotic_assemblage < 200, "G2", 
                                   ifelse(biotic_assemblage < 300, "G1", unit))), 
         env = str_sub(unit, 1,2)) %>%
  na.omit()
  
  

 ggplot(area_summary, aes(x = plot_unit, fill = as.factor(env), y = area)) + 
   geom_col() + facet_wrap(~ACBR_Name, scales = "free", ncol = 3) + 
   scale_fill_manual(values = c(viridis(5), "darkred", "tomato3", "pink", "dodgerblue" )) + 
   labs(fill = "Unit", y = "Area (square km)", x = "Biotic assemblage") + 
   theme(axis.text.x = element_blank(), 
         panel.grid.minor = element_blank(), 
         panel.grid.major.x = element_blank())

## summary of unit areas overall
area_summary %>% filter(!grepl("BNA", unit, fixed = TRUE)) %>% 
  group_by(env, plot_unit) %>% summarise(area = sum(area)) %>%
  ggplot(aes(x = plot_unit, fill = as.factor(env), y = area)) + geom_col() + 
  scale_fill_manual(values = c(viridis(5), "darkred", "tomato3", "pink", "dodgerblue" )) + 
  labs(fill = "Unit", y = "Area (square km)", x = "Biotic assemblage") +
  theme(axis.text.x  = element_text(angle = 45))


## endemics 
PA01 <- tobinary.single(PA)
endemics <- merge(area_summary, PA01[which(rowSums(PA01) ==1),] %>% colSums(), by.x = "gridcode", by.y = 0) %>% rename(endemics = y)
ggplot(endemics, aes(x = gridcode, fill = as.factor(env), y = endemics)) + geom_col() + scale_fill_viridis(discrete = TRUE) + labs(fill = "Env group") 
