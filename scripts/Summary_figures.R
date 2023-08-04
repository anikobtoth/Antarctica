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
ovl_area <- read_csv("../Data/Base/typV6_v1pl_areas.txt") # generated using TabulateAreas() in arcGIS
ACBR_key <- readRDS("./data/ACBR_key.rds")
habitats <- read_csv("./documents/Unit_descriptionsV6.csv") %>% 
  dplyr::select(-Distribution, -Description) %>% 
  filter(!is.na(Gridcode))

area_summary <- ovl_area %>% pivot_longer(contains("VALUE_"), names_to = "habitat", values_to = "area") %>% 
  dplyr::select(-Rowid_) %>% mutate(habitat = as.numeric(word(habitat, 2,2, "_")), 
                                    area = area/1000000) %>%
  left_join(ACBR_key, by = c("ACBR_NAME" = "ACBR_Name")) %>% 
  mutate(super_unit = (habitat %/% 100)*100, 
         T2 = habitat%%100, 
         habitat = ifelse(super_unit>0, super_unit, habitat)) %>%
  left_join(habitats, by = c("habitat" = "Gridcode")) %>%
  mutate(plot_unit = LCODE) %>%
  mutate(env = str_sub(LCODE, 1,2)) %>%
  mutate(env = ifelse(env %in% c("G1", "G2", "G3"), "G", env)) %>%
  filter(area > 0)
  
  

 ggplot(area_summary %>% filter(ACBR_NAME != "South Orkney Islands") %>% na.omit(), aes(x = as.factor(LCODE), fill = as.factor(env), y = area)) + 
   geom_col() + facet_wrap(~ACBR_NAME, scales = "free_y", ncol = 3) + 
   scale_fill_manual(values = c(viridis(5), "gray20", "dodgerblue"), 
                     labels = c("Mild lowlands", "Humid midlands", "Sunny/dry midlands", "High mountains", "High flatlands", "Geothermal", "Lakes")) + 
   labs(fill = "Unit", y = "Area (square km)", x = "Habitat complex") + 
   theme(axis.text.x = element_blank(), 
         panel.grid.minor = element_blank(), 
         panel.grid.major.x = element_blank())

## summary of unit areas overall
area_summary %>% #filter(ACBR_NAME != "South Orkney Islands") %>% 
  group_by(env, LCODE, Name) %>% dplyr::summarise(area = sum(area)) %>% na.omit() %>%
  ggplot(aes(x = LCODE, fill = as.factor(env), y = area)) + geom_col() + 
  facet_grid(.~as.factor(env), scales = "free_x", space = "free") +
  scale_fill_manual(values = c(viridis(5), "gray20", "dodgerblue"), 
                    labels = c("Mild lowlands", "Humid midlands", "Sunny/dry midlands", "High mountains", "High flatlands", "Geothermal", "Lakes")) + 
  labs(fill = "Unit", y = "Area (square km)", x = "Habitat complex") +
  theme(axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5))

# map of tier 1

typv6 <- read_csv("./results/typv6.txt")
typv6 <- typv6 %>% filter(RASTERVALU != -9999, RASTERVALU != 0, y < 2e06) %>% 
  mutate(unit = ifelse(RASTERVALU %in% 1:4, 2, ifelse(RASTERVALU %in% 6:11, 4, ifelse(RASTERVALU %in% c(12:18, 400), 1, ifelse(RASTERVALU %in% 20:24, 3, ifelse(RASTERVALU %in% 26:31, 5, ifelse(RASTERVALU %in% c(100, 200, 300), 6, 7)))))))
ggplot() +  
  geom_sf(data=antarctica, fill="gray50", color= NA, size=0.25) + 
  coord_sf(datum = st_crs(3031)) +
  geom_tile(data = typv6, aes(x = x, y = y, fill = unit, col = unit), lwd = 0.5)


# summary of geothermal units

ovl_area <- read_csv("../Data/Base/typV6_v2pl_areas.txt") # generated using TabulateAreas() in arcGIS

area_summary %>% filter(super_unit %in% c(100, 200, 300)) %>% 
  mutate(superunit = fct_recode(as.factor(super_unit), "volcanic" = "100", "dormant" = "200","geothermal" ="300")) %>%
  group_by(superunit, unit = LCODE, env) %>% summarise(area = sum(area)) %>%
  ggplot(aes(x = unit, fill = as.factor(str_sub(unit, 1,2)), y = area)) + geom_col() + 
  scale_fill_manual(values = c(viridis(5))) +
  facet_grid(~superunit, scales = "free", space = "free")+
  labs(fill = "Unit", y = "Area (square km)", x = "Habitat complex") +
  theme(axis.text.x  = element_text(angle = 90), 
        panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank()) +
  scale_y_log10()


## endemics 
PA01 <- tobinary.single(PA)
endemics <- merge(area_summary, PA01[which(rowSums(PA01) ==1),] %>% colSums(), by.x = "gridcode", by.y = 0) %>% rename(endemics = y)
ggplot(endemics, aes(x = gridcode, fill = as.factor(env), y = endemics)) + geom_col() + scale_fill_viridis(discrete = TRUE) + labs(fill = "Env group") 
