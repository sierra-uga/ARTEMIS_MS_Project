## coastal TRANSECT (1)
stations <- c(casts[[152]], casts[[201]], casts[[168]], casts[[25]], casts[[15]], casts[[141]],
              casts[[65]], casts[[127]], casts[[215]])
coastal <- as.section(stations, sectionId = "coastal")

# variables
s <- sectionGrid(coastal, p='levitus')
nstation <- length(s[['station']])
p <- unique(s[['pressure']])
np <- length(p)

# for basemap, points
lat <- lon <- array(NA, dim=c(9, 1))
for (i in 1:nstation) {
  lat[i, ] <- s[['station']][[i]][['latitude']]
  lon[i, ] <- s[['station']][[i]][['longitude']]
}
lat <- as.data.frame(lat)
lon <- as.data.frame(lon)
stationz <- c("89", "132", "106", "20", "14", "78", "56", "68", "146")
pointz <- do.call(rbind, Map(data.frame, lat=lat, lon=lon))
pointz$station <- stationz

# actual map of stations
coastal_map <- basemap(data = pointz, limits = c(-115, -109.5, -75, -73), bathymetry = TRUE, bathy.style = "rcb", rotate = TRUE, glaciers = TRUE) +
  ggspatial::geom_spatial_point(data = pointz, aes(x = lon, y = lat), color = "red") + ggspatial::geom_spatial_text_repel(
    data = pointz, aes(x = lon, y = lat, label = station))

# 
T <- S <- array(NA, dim=c(nstation, np))
for (i in 1:nstation) {
  T[i, ] <- s[['station']][[i]][['temperature']]
  S[i, ] <- s[['station']][[i]][['salinity']]
}

distance <- unique(s[['distance']])
coastal_labels <- c("89", "132", "106", "20", "14", "78", "56", "68", "146")

pdf(file="coastal_section_alla.pdf", width=8, height=6)
par(mfrow=c(2, 2))
plot(coastal, which='map')
Tcm <- colormap(T, breaks=seq(-2, 2, 0.3), col=oceColorsTemperature)
Scm <- colormap(S, breaks=seq(33, 35, 0.2), col=oceColorsSalinity)
imagep(distance, p, T, colormap=Tcm, flipy=TRUE,
       ylab='p [dbar]', filledContour=TRUE,
       zlab='temperature [degC]', drawContours = FALSE, xlim=c(min(distance), max(distance)))
plot(coastal, which='Conservative Temperature', xtype='distance',
     ztype='image', ylim=c(1200, 0), stationTicks=TRUE, grid=TRUE)
imagep(distance, p, S, colormap=Scm, flipy=TRUE,
       xlab='distance [km]', ylab='p [dbar]', drawContours = TRUE, filledContour=TRUE,
       zlab='salinity', xlim=c(min(distance), max(distance)))

coastal_map

par(mfrow=c(2, 1))
imagep(distance, p, T, col=oceColorsTemperature, flipy=TRUE)
imagep(distance, p, S, col=oceColorsSalinity, flipy=TRUE)
dev.off()
