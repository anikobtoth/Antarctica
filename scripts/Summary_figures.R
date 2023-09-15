library(ggplot2)
library(viridis)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)


## Summary of unit areas distribution in ACBRs V6
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
 
 # Bubble plots ####
 area_summary %>% filter(ACBR_NAME != "South Orkney Islands", area >= 50) %>% 
   na.omit() %>%
   ggplot(aes(x = LCODE, y = ACBR_NAME, col = area, size = area)) + geom_point() +
   scale_color_viridis() +
   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# by proportion
 area_summary %>% group_by(ACBR_NAME) %>% summarise(ACBR_area = sum(area)) %>%
   full_join(area_summary) %>% filter(ACBR_NAME != "South Orkney Islands") %>% 
   na.omit() %>% mutate(proportion = area/ACBR_area) %>% 
   filter(proportion > 0.05) %>%
   ggplot(aes(x = LCODE, y = ACBR_NAME, col = proportion, size = area)) + geom_point() +
   scale_color_viridis() + scale_size_area(breaks = c(50, 250, 500, 1000, 1500, 2000)) +
   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
 
 
 area_summary %>% group_by(ACBR_NAME) %>% summarise(ACBR_area = sum(area)) %>%
   full_join(area_summary) %>% filter(ACBR_NAME != "South Orkney Islands") %>% 
   na.omit() %>% mutate(proportion = area/ACBR_area) %>% 
   filter(proportion > 0.05) %>%
   ggplot(aes(x = LCODE, y = ACBR_NAME, col = area, size = proportion)) + geom_point() +
   scale_color_viridis() +
   labs(x = "Habitat Complex") +
   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
         axis.title.y = element_blank())
 
## Name abbreviations####
 
area_summary <- area_summary %>% 
  mutate(Name_abbr = Name %>% 
            str_replace_all("land", "l") %>% 
            str_replace_all("Sunny", "Sun.") %>% 
            str_replace_all("exposed", "exp.") %>%
            str_replace_all("Exposed", "Exp.") %>%
            str_replace_all("sheltered", "shel.") %>%
            str_replace_all("Sheltered", "shel.") %>%
            str_replace_all("coastal", "cstl.") %>%
            str_replace_all("outcrops", "outc.") %>%
            str_replace_all("meltwater-eroded", "melt.") %>%
            str_replace_all("meltwater", "melt.") %>%
            str_replace_all("moraine fields", "morf.") %>%
            str_replace_all("Meltwater-eroded", "Melt.") %>%
            str_replace_all("escarpments", "escrp.") %>%
            str_replace_all("mountain slopes", "msl.") %>%
            str_replace_all("slopes", "sl.") %>%
            str_replace_all("and", "&") %>%
            str_replace_all("inclines", "incl.") %>%
            str_replace_all("plateaus", "plat.") %>%
            str_replace_all("nunatak faces &", "") %>% 
            str_replace_all("nunatak", "ntk.") %>% 
            str_replace_all("fringes", "") %>% 
            str_replace_all("screes", "scr.") %>% 
            str_replace_all("ecosystems", "") %>% 
            str_replace_all("geothermal", "geoth.") %>% 
            str_replace_all("volcanic", "volc") %>% 
            str_replace_all("transitional", "trans.") %>% 
            str_replace_all("  ", " "))
 
## summary of unit areas overall ####
area_summary %>% filter(ACBR_NAME != "South Orkney Islands") %>% 
  mutate(LCODE = ifelse(grepl("L", LCODE), "L1", LCODE)) %>% 
  group_by(env, LCODE, Name, Name_abbr) %>% dplyr::summarise(area = sum(area)) %>% na.omit() %>%
  ggplot(aes(x = paste(LCODE, Name_abbr, sep = "-"), fill = as.factor(env), y = area)) + geom_col() + 
  facet_grid(.~as.factor(env), scales = "free_x", space = "free") +
  scale_fill_manual(values = c(viridis(5), "gray20", "dodgerblue"), 
                    labels = c("Mild lowlands", "Humid midlands", "Sunny/dry midlands", "High mountains", "High flatlands", "Geothermal", "Lakes")) + 
  labs(fill = "Unit", y = "Area (square km)", x = "Habitat complex") +
  theme(axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0))

# map of tier 1 ####

typv6 <- read_csv("./results/typv6.txt")
typv6 <- typv6 %>% filter(RASTERVALU != -9999, RASTERVALU != 0, y < 2e06) %>% 
  mutate(unit = ifelse(RASTERVALU %in% 1:4, 2, ifelse(RASTERVALU %in% 6:11, 4, ifelse(RASTERVALU %in% c(12:18, 400), 1, ifelse(RASTERVALU %in% 20:24, 3, ifelse(RASTERVALU %in% 26:31, 5, ifelse(RASTERVALU %in% c(100, 200, 300), 6, 7)))))))
ggplot() +  
  geom_sf(data=antarctica, fill="gray50", color= NA, size=0.25) + 
  coord_sf(datum = st_crs(3031)) +
  geom_tile(data = typv6, aes(x = x, y = y, fill = unit, col = unit), lwd = 0.5)

# Bioregional Ecosystem Types ####
library(sf)
ru <- readRDS("./results/RegionalUnitsV6_3.rds")
ru_summary <- ru %>% st_drop_geometry() %>% 
  group_by(final, ACBR_Name, ACBR_ID, unit) %>% 
  summarise(area = n()/100, comp = paste(unique(LCODE), collapse = "+")) %>%
  mutate(unit = ifelse(grepl("G", unit), "G", unit))
lev <- ru_summary %>% filter(area >=10) %>% group_by(ACBR_Name) %>% summarise(vals = n()) %>% arrange(vals) %>% pull(ACBR_Name)
p <- ru_summary %>% filter(area > 10) %>% 
  mutate(ACBR_Name = factor(ACBR_Name, levels = lev)) %>%
  ggplot(aes(x = as.factor(final), y = area, fill = unit)) + 
  geom_col() + facet_wrap(.~ACBR_Name, scales = "free", dir = "v", ncol = 3) + 
  scale_fill_manual(values = c(viridis(5), "gray20", "dodgerblue"), 
                    labels = c("Mild lowlands", "Humid midlands", "Sunny/dry midlands", "High mountains", "High flatlands", "Geothermal", "Lakes")) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, size = 8), strip.text = element_text(size = 7)) +
  labs(x = "Bioregional Ecosystem Type", y = "Area (square km)", fill = "Major Env. Unit")

# convert ggplot object to grob object
gp <- ggplotGrob(p)

# optional: take a look at the grob object's layout
#gtable::gtable_show_layout(gp)

# get gtable columns corresponding to the facets (5 & 9, in this case)
facet.columns <- gp$layout$l[grepl("panel", gp$layout$name)]

# get the number of unique x-axis values per facet (1 & 3, in this case)
x.var <- sapply(ggplot_build(p)$layout$panel_scales_x,
                function(l) length(l$range$range))

# change the relative widths of the facet columns based on
# how many unique x-axis values are in each facet
gp$widths[facet.columns] <- gp$widths[facet.columns] * x.var

# plot result
grid::grid.draw(gp)

# summary of geothermal units ####

ovl_area <- read_csv("../Data/Base/typV6_v2pl_areas.csv") # generated using TabulateAreas() in arcGIS

ovl_area %>% pivot_longer(contains("VALUE_"), names_to = "habitat", values_to = "area") %>% 
  dplyr::select(-Rowid) %>% mutate(habitat = as.numeric(word(habitat, 2,2, "_")), 
                                    area = area/1000000) %>% 
  filter(habitat > 99, habitat < 400, ACBR_NAME != "South Orkney Islands", area > 0.5) %>%
  left_join(ACBR_key, by = c("ACBR_NAME" = "ACBR_Name")) %>% 
  mutate(super_unit = (habitat %/% 100)*100, 
         Gridcode = habitat%%100) %>%
  left_join(habitats, by = "Gridcode") %>%
  mutate(plot_unit = LCODE) %>%
  mutate(env = str_sub(LCODE, 1,2), 
         superunit = fct_recode(as.factor(super_unit), "G1" = "100", "G2" = "200","G3" ="300")) %>% 
  group_by(superunit, unit = LCODE, env) %>% summarise(area = sum(area)) %>%
  ggplot(aes(x = unit, fill = as.factor(str_sub(unit, 1,2)), y = area)) + geom_col() + 
  scale_fill_manual(values = c(viridis(5))) +
  facet_grid(~superunit, scales = "free", space = "free")+
  labs(fill = "Unit", y = "Area (square km)", x = "Habitat complex") +
  theme(axis.text.x  = element_text(angle = 90, vjust = 0.3), 
        panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank()) 


