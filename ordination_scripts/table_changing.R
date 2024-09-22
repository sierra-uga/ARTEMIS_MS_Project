library(readr)
library("dplyr")
library("tidyr")

water_properties_csv <- read_csv("ARTEMIS_github/required_files/water_properties.csv")
water_properties <- filter(water_properties_csv, Sample.Control == "True.Sample")
water_properties <- water_properties %>% select(c(-Sample.Control))

water_properties <- water_properties %>% 
  mutate(Community = replace(Filter_pores, Filter_pores == 0.2, "Free-living")) %>%
  mutate(Community = replace(Filter_pores, Filter_pores >= 2, "Particle-associated")) %>%
  mutate(Community = replace(Community, Community =="0.2", "Free-living")) 

water_properties <- water_properties %>% relocate(Community, .after = Sample)
write.csv(water_properties, "~/Documents/Research/Ordination analysis R scripts/ARTEMIS_github/final_graphics/water_properties_table.csv")