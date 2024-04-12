## ggoceanmaps testing
library(tidyverse)
if("package:PlotSvalbard" %in% search()) { detach("package:PlotSvalbard", unload=TRUE) }
library(ggOceanMaps)
library(oce)

# MAP of all stations
dt <- data
dt <- select(data, c("Station", "Latitude", "Longitude")) %>% na.omit(.) %>% distinct(., Station, .keep_all= TRUE) %>%
  rename(., c(lat = Latitude, lon = Longitude))

all_points <- basemap(data = dt, limits = c(-120, -107, -75, -71.5), bathymetry = TRUE, rotate = TRUE, glaciers = TRUE) +
  ggspatial::geom_spatial_point(data = dt, aes(x = lon, y = lat), color = "red") + ggspatial::geom_spatial_text_repel(
  data = dt, aes(x = lon, y = lat, label = Station)) + scale_color_manual(values = c("red", "yellow"))

all_points
ggsave("updated_map.pdf", width = 9, height = 8)

### importing casts
files <- dir(path="~/Documents/Research/Temperature_Salinity/Temp_Salinity_Plots/ctd/true_ctd", pattern='*.cnv')
setwd("~/Documents/Research/Temperature_Salinity/Temp_Salinity_Plots/ctd/true_ctd")
casts <- list()
for (ifile in 1:length(files)) {
  casts[[ifile]] <- read.ctd.sbe(files[ifile])
}
str(casts, 1)


# section stuff
# https://www.rdocumentation.org/packages/oce/versions/1.8-2/topics/section-class

### references
# https://stackoverflow.com/questions/77172828/plot-a-temperature-profile-with-contours-over-distance-and-depth-in-r

# https://www.clarkrichards.org/2016/04/25/making-section-plots-with-oce-and-imagep/
# cross-section
