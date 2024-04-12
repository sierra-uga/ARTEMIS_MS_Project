# load libraries
library(tidyverse)
if("package:PlotSvalbard" %in% search()) { detach("package:PlotSvalbard", unload=TRUE) }
library(ggOceanMaps)
library(oce)
library(cowplot)


## meltwater TRANSECT (1)
meltwater_stations <- c(casts[[116]], casts[[128]], casts[[216]]) # also cast [[66]], [[33]] casts[[116]],  #station 22/56b
meltwater <- as.section(meltwater_stations, sectionId = "meltwater")
# variables
meltwater_s <- sectionGrid(meltwater, p=seq(0,700,1), method="lm") # best method for CTD data
meltwater_nstation <- length(meltwater_s[['station']])
meltwater_p <- unique(meltwater_s[['pressure']])
meltwater_np <- length(meltwater_p)

# for basemap, points
meltwater_lat <- meltwater_lon <- array(NA, dim=c(3, 1)) # CHANGE DIMENSIONS DEPENDING ON NUMBER OF STATION
for (i in 1:meltwater_nstation) {
  meltwater_lat[i, ] <- meltwater_s[['station']][[i]][['latitude']]
  meltwater_lon[i, ] <- meltwater_s[['station']][[i]][['longitude']]
}
meltwater_lat <- as.data.frame(meltwater_lat)
meltwater_lon <- as.data.frame(meltwater_lon)
meltwater_stationz <- c("56", "68", "146")
meltwater_pointz <- do.call(rbind, Map(data.frame, lat=meltwater_lat, lon=meltwater_lon)) # takes points from lat and lon data frame
meltwater_pointz$station <- meltwater_stationz # add station as a column

# actual map of meltwater_stations
meltwater_map <- basemap(data = meltwater_pointz, limits = c(-115, -112, -74.5, -73.5), bathymetry = TRUE, bathy.style = "rcb", rotate = TRUE, glaciers = TRUE) +
  ggspatial::geom_spatial_point(data = meltwater_pointz, aes(x = lon, y = lat), color = "red") + ggspatial::geom_spatial_text_repel(
    data = meltwater_pointz, aes(x = lon, y = lat, label = station))

# 
meltwater_T <- meltwater_O <- array(NA, dim=c(meltwater_nstation, meltwater_np))
for (i in 1:meltwater_nstation) {
  meltwater_T[i, ] <- meltwater_s[['station']][[i]][['potential temperature']]
  meltwater_O[i, ] <- meltwater_s[['station']][[i]][['oxygen']]
}
meltwater_distance <- unique(meltwater_s[['distance']])
# output: [1]  0.00000000  0.02119082  0.02403822 16.33216538 40.58487518
meltwater_labels <- c("56", "68", "146")
meltwater_depth <- c(190, 190, 200)
#convert depth to pressure
meltwater_depths <- swPressure(meltwater_depth, latitude = -73) #convert depth to pressure
meltwater_ox_min <- 4.5
meltwater_ox_max <- 9

# colormaps
meltwater_Tcm <- colormap(meltwater_T, breaks=seq(-2, 1, 0.2), col=oceColorsTemperature) # parameters for Temperature colormap
meltwater_Ocm <- colormap(meltwater_O, breaks=seq(meltwater_ox_min, meltwater_ox_max, 0.5), col=oceColorsViridis) # parameters for Oxygen colormap

######

pdf(file="~/Documents/Research/Temperature_Salinity/Temp_Salinity_Plots/graphics/meltwater_plume_oxygen_map_draft.pdf", width=5, height=7)
par(mfrow=c(2, 1))
# plot 1 - small map
plot(meltwater, which='map')

# plot 2 - section plot of temp
imagep(meltwater_distance, meltwater_p, meltwater_T, colormap=meltwater_Tcm, flipy=TRUE,
       ylab='p [dbar]', filledContour=TRUE,
       zlab='temperature [degC]', drawContours = TRUE, xlim=c(min(meltwater_distance), max(meltwater_distance)))
points(meltwater_distance, y=meltwater_depths, pch=25, col="black", bg="red")
# plot 3 - section plot of oxygen
imagep(meltwater_distance, meltwater_p, meltwater_O, colormap=meltwater_Ocm, flipy=TRUE,
       xlab='distance [km]', ylab='p [dbar]', drawContours = TRUE, filledContour=TRUE, axes = TRUE,
       zlab='oxygen [mL/L]', xlim=c(min(meltwater_distance), max(meltwater_distance))) 
points(meltwater_distance, y=meltwater_depths, pch=25, col="black", bg="red", cex=1.2)

meltwater_map
dev.off()

#######

# one long plot (map) + 2 section plots
pdf(file="~/Documents/Research/Temperature_Salinity/Temp_Salinity_Plots/graphics/meltwater_plume_oxygen_map.pdf", width=9, height=6)
# specify the matrix
matrix(c(1, 1, 2, 3), nrow = 2, byrow = FALSE)

# 3 plots to be combined in 2 row/ 2 columns and arranged by columns
layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = FALSE), widths= c(1,2))

# specify the 3 plots
# plot 1 - small map
plot(meltwater, which='map', map.xlim=c(-115, -112), map.ylim=c(-74.5, -73.5), pch=25, col="red")
# plot 2 - section plot of temp
imagep(meltwater_distance, meltwater_p, meltwater_T, colormap=meltwater_Tcm, flipy=TRUE,
       ylab='p [dbar]', filledContour=TRUE,
       zlab='temperature [degC]', drawContours = TRUE, xlim=c(min(meltwater_distance), max(meltwater_distance)))
points(meltwater_distance, y=meltwater_depths, pch=25, col="black", bg="red")
# plot 3 - section plot of oxygen
imagep(meltwater_distance, meltwater_p, meltwater_O, colormap=meltwater_Ocm, flipy=TRUE,
       xlab='distance [km]', ylab='p [dbar]', drawContours = TRUE, filledContour=TRUE,
       zlab='oxygen [mL/L]', xlim=c(min(meltwater_distance), max(meltwater_distance))) 
points(meltwater_distance, y=meltwater_depths, pch=25, col="black", bg="red", cex=1.2)

meltwater_map
dev.off()




