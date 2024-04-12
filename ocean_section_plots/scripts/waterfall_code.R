# load libraries
library(gridBase)
library(grid)
library(patchwork)

## WATERFALL TRANSECT (1)
waterfall_stations <- c(casts[[269]], casts[[3]], casts[[5]], casts[[13]], casts[[16]]) # chosen by hand, because stations are a little bit off from bottle file
waterfall <- as.section(waterfall_stations, sectionId = "waterfall")
# variables
waterfall_s <- sectionGrid(waterfall, p=seq(0,1750,5), method="lm") # best method for CTD data
waterfall_nstation <- length(waterfall_s[['station']])
waterfall_p <- unique(waterfall_s[['pressure']])
waterfall_np <- length(waterfall_p)

# for basemap, points
waterfall_lat <- waterfall_lon <- array(NA, dim=c(5, 1)) # CHANGE DIMENSIONS DEPENDING ON NUMBER OF STATION
for (i in 1:waterfall_nstation) {
  waterfall_lat[i, ] <- waterfall_s[['station']][[i]][['latitude']]
  waterfall_lon[i, ] <- waterfall_s[['station']][[i]][['longitude']]
}
waterfall_lat <- as.data.frame(waterfall_lat)
waterfall_lon <- as.data.frame(waterfall_lon)
waterfall_stationz <- c("198", "2", "4", "12", "14") # unique to specific transect
waterfall_pointz <- do.call(rbind, Map(data.frame, lat=waterfall_lat, lon=waterfall_lon)) # takes points from lat and lon data frame
waterfall_pointz$station <- waterfall_stationz # add station as a column

# actual map of waterfall_stations
waterfall_map <- basemap(data = waterfall_pointz, limits = c(-120, -111, -75, -71.5), bathymetry = TRUE, bathy.style = "rcb", rotate = TRUE, glaciers = TRUE) +
  ggspatial::geom_spatial_point(data = waterfall_pointz, aes(x = lon, y = lat), color = "red") + ggspatial::geom_spatial_text_repel(
    data = waterfall_pointz, aes(x = lon, y = lat, label = station))

# 
waterfall_T <- waterfall_O <- array(NA, dim=c(waterfall_nstation, waterfall_np))
for (i in 1:waterfall_nstation) {
  waterfall_T[i, ] <- waterfall_s[['station']][[i]][['temperature']]
  waterfall_O[i, ] <- waterfall_s[['station']][[i]][['oxygen2']]
}
waterfall_distance <- unique(waterfall_s[['distance']])
waterfall_labels <- c("STN198","STN002", "STN004", "STN012", "STN014")
waterfall_depth <- c(375, 505, 545, 715, 700)
#convert depth to pressure
waterfall_depths <- swPressure(waterfall_depth, latitude = -73)
waterfall_ox_min <- min(waterfall_stations[[1]]@data[["oxygen"]])
waterfall_ox_max <- max(waterfall_stations[[1]]@data[["oxygen"]])

# colormaps
waterfall_Tcm <- colormap(waterfall_T, breaks=seq(-2, 2, 0.5), col=oceColorsTemperature) # parameters for Temperature colormap
waterfall_Ocm <- colormap(waterfall_O, breaks=seq(waterfall_ox_min, waterfall_ox_max, 0.5), col=oceColorsViridis) # parameters for Oxygen colormap

#pdf(file="waterfall_section_alla.pdf", width=8, height=6)
par(mfrow=c(2, 2))
#plot 1 - small map
plot(waterfall, which='map')
imagep(waterfall_distance, waterfall_p, waterfall_T, colormap=waterfall_Tcm, flipy=TRUE,
       ylab='p [dbar]', filledContour=TRUE,
       zlab='temperature [degC]', drawContours = TRUE, xlim=c(min(waterfall_distance), max(waterfall_distance)))
points(waterfall_distance, y=waterfall_depths, pch=25, col="black", bg="red")

plot(waterfall, which='potential temperature', xtype='distance',
     ztype='image', ylim=c(1700, 0), stationTicks=TRUE, grid=FALSE, coastline="coastlineWorld") 
points(waterfall_distance, y=waterfall_depth, pch=25, col="black", bg="red")

imagep(waterfall_distance, waterfall_p, waterfall_O, colormap=waterfall_Ocm, flipy=TRUE,
       xlab='distance [km]', ylab='p [dbar]', drawContours = TRUE, filledContour=TRUE,
       zlab='oxygen [mL/L]', xlim=c(min(waterfall_distance), max(waterfall_distance))) 
points(waterfall_distance, y=waterfall_depths, pch=25, col="black", bg="red", cex=1.2)

waterfall_map
dev.off()


# one long plot (map) + 2 section plots
pdf(file="~/Documents/Research/Temperature_Salinity/Temp_Salinity_Plots/graphics/waterfall_oxygen_map.pdf", width=9, height=6)
# specify the matrix
matrix(c(1, 1, 2, 3), nrow = 2, byrow = FALSE)

# 3 plots to be combined in 2 row/ 2 columns and arranged by columns
layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = FALSE), widths= c(1,2))

# specify the 3 plots
# plot 1 - small map
plot(waterfall, which='map', map.xlim=c(-120, -111), map.ylim=c(-75, -71.5), pch=25, col="red")
# plot 2 - section plot of temp
imagep(waterfall_distance, waterfall_p, waterfall_T, colormap=waterfall_Tcm, flipy=TRUE,
       ylab='p [dbar]', filledContour=TRUE,
       zlab='temperature [degC]', drawContours = TRUE, xlim=c(min(waterfall_distance), max(waterfall_distance)))
points(waterfall_distance, y=waterfall_depths, pch=25, col="black", bg="red")
# plot 3 - section plot of oxygen
imagep(waterfall_distance, waterfall_p, waterfall_O, colormap=waterfall_Ocm, flipy=TRUE,
       xlab='distance [km]', ylab='p [dbar]', drawContours = TRUE, filledContour=TRUE,
       zlab='oxygen [mL/L]', xlim=c(min(waterfall_distance), max(waterfall_distance))) 
points(waterfall_distance, y=waterfall_depths, pch=25, col="black", bg="red", cex=1.2)

waterfall_map
dev.off()