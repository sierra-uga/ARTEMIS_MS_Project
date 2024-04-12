# load libraries
library(gridBase)
library(grid)
library(patchwork)

## coastal TRANSECT (1)
coastal_stations <- c(casts[[153]], casts[[202]], casts[[169]], casts[[26]], casts[[16]], casts[[142]],
                      casts[[66]], casts[[128]], casts[[216]]) # chosen by hand, because stations are a little bit off from bottle file
coastal <- as.section(coastal_stations, sectionId = "coastal")
# variables
coastal_s <- sectionGrid(coastal, p=seq(0,1300,5), method="boxcar") # best method for CTD data
coastal_nstation <- length(coastal_s[['station']])
coastal_p <- unique(coastal_s[['pressure']])
coastal_np <- length(coastal_p)

# for basemap, points
coastal_lat <- coastal_lon <- array(NA, dim=c(9, 1)) # CHANGE DIMENSIONS DEPENDING ON NUMBER OF STATION
for (i in 1:coastal_nstation) {
  coastal_lat[i, ] <- coastal_s[['station']][[i]][['latitude']]
  coastal_lon[i, ] <- coastal_s[['station']][[i]][['longitude']]
}
coastal_lat <- as.data.frame(coastal_lat)
coastal_lon <- as.data.frame(coastal_lon)
coastal_stationz <- c("89", "132", "106", "20", "14", "78", "56", "68", "146") # unique to specific transect
coastal_pointz <- do.call(rbind, Map(data.frame, lat=coastal_lat, lon=coastal_lon)) # takes points from lat and lon data frame
coastal_pointz$station <- coastal_stationz # add station as a column

# actual map of coastal_stations
coastal_map <- basemap(data = coastal_pointz, limits = c(-115, -109.5, -75, -73), bathymetry = TRUE, bathy.style = "rcb", rotate = TRUE, glaciers = TRUE) +
  ggspatial::geom_spatial_point(data = coastal_pointz, aes(x = lon, y = lat), color = "red") + ggspatial::geom_spatial_text_repel(
    data = coastal_pointz, aes(x = lon, y = lat, label = station))

# 
coastal_T <- coastal_O <- array(NA, dim=c(coastal_nstation, coastal_np))
for (i in 1:coastal_nstation) {
  coastal_T[i, ] <- coastal_s[['station']][[i]][['temperature']]
  coastal_O[i, ] <- coastal_s[['station']][[i]][['oxygen2']]
}
coastal_distance <- unique(coastal_s[['distance']])
# output: [1]   0.00000  32.67378  48.56614  61.17885  62.78023  80.58525 101.50049 106.99770 133.80764
coastal_labels <-  c("89", "132", "106", "20", "14", "78", "56b", "68", "146")
coastal_depth <- c(200, 200, 200, 175, 175, 200, 180, 190, 200)
#convert depth to pressure
coastal_depths <- swPressure(coastal_depth, latitude = -73)
coastal_ox_min <- 4.5
coastal_ox_max <- 9

# colormaps
coastal_Tcm <- colormap(coastal_T, breaks=seq(-2, 2, 0.5), col=oceColorsTemperature) # parameters for Temperature colormap
coastal_Ocm <- colormap(coastal_O, breaks=seq(coastal_ox_min, coastal_ox_max, 0.5), col=oceColorsViridis) # parameters for Oxygen colormap

#pdf(file="coastal_section_alla.pdf", width=8, height=6)
par(mfrow=c(1, 1))
#plot 1 - small map
plot(coastal, which='map')
imagep(coastal_distance, coastal_p, coastal_T, colormap=coastal_Tcm, flipy=TRUE,
       ylab='p [dbar]', filledContour=TRUE,
       zlab='temperature [degC]', drawContours = TRUE, xlim=c(min(coastal_distance), max(coastal_distance)))
points(coastal_distance, y=coastal_depths, pch=25, col="black", bg="red")

imagep(coastal_distance, coastal_p, coastal_O, colormap=coastal_Ocm, flipy=TRUE,
       xlab='distance [km]', ylab='p [dbar]', drawContours = TRUE, filledContour=TRUE,
       zlab='oxygen [mL/L]', xlim=c(min(coastal_distance), max(coastal_distance))) 
points(coastal_distance, y=coastal_depths, pch=25, col="black", bg="red", cex=1.2)

coastal_map
dev.off()

plot(coastal, which='potential temperature', xtype='distance',
     ztype='image', ylim=c(1300, 0), stationTicks=TRUE, grid=FALSE, coastline="coastlineWorld") 
points(coastal_distance, y=coastal_depth, pch=25, col="black", bg="red")


dotson_distance <- c(62.78023, 80.58525, 101.50049)
dotson_depths <- c(176.8954, 202.1785, 181.9518)
# one long plot (map) + 2 section plots
pdf(file="~/Documents/Research/Temperature_Salinity/Temp_Salinity_Plots/graphics/coastal_oxygen_map.pdf", width=9, height=6)
# specify the matrix
matrix(c(1, 1, 2, 3), nrow = 2, byrow = FALSE) # makes the grid for the figures

# 3 plots to be combined in 2 row/ 2 columns and arranged by columns
layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = FALSE), widths= c(1,2))

# specify the 3 plots
# plot 1 - small map
plot(coastal, which='map', map.xlim=c(-115, -109.5), map.ylim=c(-75, -73), pch=25, col="red")
# plot 2 - section plot of temp
imagep(coastal_distance, coastal_p, coastal_T, colormap=coastal_Tcm, flipy=TRUE,
       ylab='p [dbar]', filledContour=TRUE,
       zlab='temperature [degC]', drawContours = TRUE, xlim=c(max(coastal_distance), min(coastal_distance))) # flip to follow the transect
points(coastal_distance, y=coastal_depths, pch=25, col="black", bg="red", cex=1.2)
points(dotson_distance, y=dotson_depths, pch=25, col="black", bg="dodgerblue", cex=1.2)
# plot 3 - section plot of oxygen
imagep(coastal_distance, coastal_p, coastal_O, colormap=coastal_Ocm, flipy=TRUE,
       xlab='distance [km]', ylab='p [dbar]', drawContours = TRUE, filledContour=TRUE,
       zlab='oxygen [mL/L]', xlim=c(max(coastal_distance), min(coastal_distance)))  # flip to follow the transect
points(coastal_distance, y=coastal_depths, pch=25, col="black", bg="red", cex=1.2)
points(dotson_distance, y=dotson_depths, pch=25, col="black", bg="dodgerblue", cex=1.2)

coastal_map
dev.off()