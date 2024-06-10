## ggoceanmaps testing
library(tidyverse)
if("package:PlotSvalbard" %in% search()) { detach("package:PlotSvalbard", unload=TRUE) }
install.packages("ggpsatial")
library(ggOceanMaps)
library(ggspatial)
library(oce)

# MAP of all stations

metadata <- filter(metatable, Sample.Control == "True.Sample")# filter by transect
metadata <- filter(metadata, Iron != "NA")# filter by transect
metadata <- filter(metadata, More_Depth_Threshold == "Bottom")

dt <- select(metadata, c("Station", "Latitude", "Longitude", "Iron")) %>% na.omit(.) %>% distinct(., Station, .keep_all= TRUE) %>%
  rename(., c(lat = Latitude, lon = Longitude))

dt$lat <- as.numeric(dt$lat)
dt$lon <- as.numeric(dt$lon)
dt$Iron <- as.numeric(dt$Iron)

all_points <- basemap(data = dt, limits = c(-120, -107, -75, -71.5), bathymetry = FALSE, rotate = TRUE, glaciers = TRUE, legends = FALSE) 

bottom <- all_points +
  ggspatial::geom_spatial_point(data = dt, aes(x = lon, y = lat, fill=Iron), color = "black", pch=23, cex=3.5) + ggspatial::geom_spatial_text_repel(
  data = dt, aes(x = lon, y = lat, label = Station), cex=2.5, vjust=-1.2) + scale_fill_gradient(low="green", high="red", limits=c(0,1.2)) + ggtitle("Bottom water") +
   theme(legend.position = "bottom",
         text = element_text(family = "Helvetica"),
         axis.title.x = element_blank(),
         axis.title.y = element_blank())

ggsave("map_with_iron_overlay_deep.pdf", width = 10, height = 8)

surface

## mid-surface
metadata <- filter(metadata, More_Depth_Threshold == "Mid-Surface")
dt <- select(metadata, c("Station", "Latitude", "Longitude", "Iron")) %>% na.omit(.) %>% distinct(., Station, .keep_all= TRUE) %>%
  rename(., c(lat = Latitude, lon = Longitude))

dt$lat <- as.numeric(dt$lat)
dt$lon <- as.numeric(dt$lon)
dt$Iron <- as.numeric(dt$Iron)

all_points <- basemap(data = dt, limits = c(-120, -107, -75, -71.5), bathymetry = FALSE, rotate = TRUE, glaciers = TRUE, legends = FALSE) 

mid_surface <- all_points +
  ggspatial::geom_spatial_point(data = dt, aes(x = lon, y = lat, fill=Iron), color = "black", pch=23, cex=3.5) + ggspatial::geom_spatial_text_repel(
    data = dt, aes(x = lon, y = lat, label = Station), cex=2.5, vjust=-1.2) + scale_fill_gradient(low="green", high="red", limits=c(0,1.10)) + ggtitle("Mixed-Layer (90-300m)") +
  theme(legend.position = "right",
        text = element_text(family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

ps_combined <- ggarrange(
  surface, mid_surface, labels = NULL, ncol=1,
  common.legend = TRUE, legend = "right"
)



ggsave("map_with_iron_overlay2.pdf", width = 11, height = 15)

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
