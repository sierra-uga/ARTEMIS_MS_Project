library(oce)
#install.packages("ocedata")
library(ocedata)

#importing ctd data from artemis files
base_path <- "~/Documents/Research/Ordination analysis R scripts/ARTEMIS_github" # needed to use full path - so change this or remove from "files" vector to get the actual path
subdirectory <- "ocean_section_plots/required_files/ctd" # if downloaded from github
files <- dir(paste(base_path, subdirectory, sep = "/"), full.names = TRUE, pattern = ".cnv") #works

# loop through list of the files to get the specific cast files
casts <- list()
for (ifile in 1:length(files)) {
  casts[[ifile]] <- read.ctd.sbe(files[ifile])
}


# Define the list of casts (assuming you have something similar)
casts <- list()  # This should be your actual list of cast objects

# Target station values
target_stations <- c("89", "132", "106", "20", "198", "2", "4", "12", "115", "12.3", "14", "78", "56a", "56b", "22", "68", "146", "181", "174", "151.2", "153")

# Initialize a list to hold the results
matched_casts <- list()

# Loop through each cast
for (i in seq_along(casts)) {
  # Extract station metadata, assuming casts are S4 objects
  station_value <- casts[[i]]@metadata[["station"]]
  
  # Check if the station value is in the target list
  if (station_value %in% target_stations) {
    matched_casts[[length(matched_casts) + 1]] <- casts[[i]]
  }
}


# Now matched_casts contains all the casts whose station metadata matches the target values
# You can print or return this list
print(matched_casts)

updat <- c(matched_casts, casts[[3]])

plotTS(casts)

coastal_stations <- c(casts[[153]], casts[[202]], casts[[169]], casts[[26]], casts[[16]], casts[[142]],
                      casts[[66]], casts[[128]], casts[[216]])
meltwater_stations <- c(casts[[116]], casts[[128]], casts[[216]])
waterfall_stations <- c(casts[[269]], casts[[3]], casts[[5]], casts[[13]], casts[[16]], casts[[66]], casts[[33]], casts[[116]])  #station 22/56b

all_stations <- unique(c(coastal_stations, meltwater_stations, waterfall_stations))
map_stations <- as.section(casts, sectionId = "map")

station_numbers <- sapply(matched_casts, function(cast) cast@metadata[["station"]])


data <- read.table("ocean_section_plots/required_files/CTD_meta.csv", sep=",", header=TRUE)
dt <- select(data, c("Station", "Latitude", "Longitude")) %>% na.omit(.) %>% distinct(., Station, .keep_all= TRUE) %>%
  rename(., c(Latitude = lat, Longitude = lon))

test_ctd <- as.ctd(data$Salinity, data$Temperature, data$PrDM, data$Longitude, data$Latitude, data$"Sbeox0ML/L", station=data$Station)
test_sec <- as.section(data$Salinity, data$Temperature, data$PrDM, data$Longitude, data$Latitude, station=data$Station)
test_s <- sectionGrid(test_sec, p=seq(0,1500,5), method="boxcar") # best method for CTD data

drawPalette(colormap=test_s, zlab="Oxygen")

plotTS(all_s, pch=19, col=all_Ocm$zcol, mar=par("mar"), type="p") # the mar adjusts for the palette


all_points <- basemap(data = dt, limits = c(-120, -109, -75, -71.5), bathymetry = TRUE, rotate = TRUE, glaciers = TRUE, bathy.style="rbb") +
  ggspatial::geom_spatial_point(data = dt, aes(x = Longitude, y = Latitude), color = "black", pch=23, cex=1.2, fill="red") + 
  ggspatial::geom_spatial_text_repel(
    data = dt, aes(x = Longitude, y = Latitude, label = Station), cex=3) +
theme(legend.position = "none",
      text = element_text(family = "Helvetica"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())

pdf(file="ocean_section_plots/graphics/POSTER_station_map_legend.pdf", width=7, height=6)
drawPalette(colormap=cm, zlab = "Depth [m]")
dev.off()
ggsave(file="ocean_section_plots/graphics/POSTER_station_map_legend.pdf", width=7, height=5)


all <- as.section(filtered_casts, sectionId = "all")
plot(all, which='map')

Eastern_CC <- c("STN20", "STN089", "STN106", "STN132")
Dotson <- c("STN014", "STN022", "STN056", "STN078")
Open_polynya <- c("STN002", "STN004", "STN012", "STN115", "STN12.3", "STN174", "STN181")
Western_CC <- c("STN068", "STN146")
Cont_Shelf <- c("STN198")
Getz <- c("STN153", "STN151.2")

east <- dt %>% filter(Station %in% Eastern_CC)
dot <- dt %>% filter(Station %in% Dotson)
west <- dt %>% filter(Station %in% Western_CC)
open <- dt %>% filter(Station %in% Open_polynya)
cont <- dt %>% filter(Station %in% Cont_Shelf)
getz <- dt %>% filter(Station %in% Getz)

all_points <- basemap(data = dt, limits = c(-120, -109, -75, -71.5), bathymetry = TRUE, rotate = TRUE, glaciers = TRUE, bathy.style="rbb") +
  ggspatial::geom_spatial_point(data = dot, aes(x = Longitude, y = Latitude), color = "black", pch=23, cex=3, fill="#5AD0FC") +
  ggspatial::geom_spatial_point(data = east, aes(x = Longitude, y = Latitude), color = "black", pch=23, cex=3, fill="darkred") +
  ggspatial::geom_spatial_point(data = west, aes(x = Longitude, y = Latitude), color = "black", pch=23, cex=3, fill="red") +
  ggspatial::geom_spatial_point(data = open, aes(x = Longitude, y = Latitude), color = "black", pch=23, cex=3, fill="#09A20D") +
  ggspatial::geom_spatial_point(data = cont, aes(x = Longitude, y = Latitude), color = "black", pch=23, cex=3, fill="#A3DCA5") +
  ggspatial::geom_spatial_point(data = getz, aes(x = Longitude, y = Latitude), color = "black", pch=23, cex=3, fill="#006B93") +
  ggspatial::geom_spatial_text_repel(
    data = dt, aes(x = Longitude, y = Latitude, label = Station), cex=3.5) +
  theme(legend.position = "none",
        text = element_text(family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
ggsave(file="ocean_section_plots/graphics/POSTER_alt_station_map_legend.pdf", width=9, height=7)


antarctica <- basemap(-60, bathymetry = TRUE, glaciers = TRUE, legends= FALSE)
ggsave(file="ocean_section_plots/graphics/antarctica_inlet.png", width=9, height=7)

all_points <- basemap(data = dt, limits = c(-120, -109, -75, -71.5), bathymetry = TRUE, rotate = TRUE, glaciers = TRUE, bathy.style="rbb") +
  theme(legend.position = "none",
        text = element_text(family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

# Define constants

#calcualte gade line
#0.5 °C and 34.55  (via paper)
#1.1 and 34.8
Tocean <- 1
Socean <- 34.8
Lf <- 334
Cp <- 3.97

# Define a function for Tp(Sp) based on the equation
Tp <- function(Sp) {
  Tocean + Lf/Cp * (1 - Socean/Sp)
}

# Generate a range of salinities
Sp <- seq(32, 38, by = 0.01)

# Calculate corresponding temperatures using the Tp function
Tp_values <- Tp(Sp)


#plot for the actual line
# Plot Tp vs. Sp
plot(Sp, Tp_values, type = "l", xlab = "Salinity (Sp)", ylim=c(-2,2), ylab = "Temperature (Tp)", main = "Temperature vs. Salinity")



plotTS(all, pch=19, col=all_Ocm$zcol, mar=par("mar"), type="p") # the mar adjusts for the palette


all_s <- sectionGrid(all, p=seq(0,1500,5), method="boxcar") # best method for CTD data
all_nstation <- length(all_s[['station']])
all_p <- unique(all_s[['pressure']])

all_np <- length(all_p)

all_T <- all_O <- all_S <- all_P <- array(NA, dim=c(all_nstation, all_np))
for (i in 1:all_nstation) {
  all_T[i, ] <- all_s[['station']][[i]][['temperature']]
  all_S[i, ] <- all_s[['station']][[i]][['salinity']]
  all_O[i, ] <- all_s[['station']][[i]][['oxygen2']]
  all_P[i, ] <- all_s[['station']][[i]][['pressure']]
}

all_ox_min <- 4
all_ox_max <- 10

all_station <- all_s[['station']]

# colormaps
all_Tcm <- colormap(all_T, breaks=seq(-2, 2, 0.5), col=oceColorsTemperature) # parameters for Temperature colormap
all_Ocm <- colormap(all_O, breaks=seq(all_ox_min, all_ox_max, 0.25), col=oceColorsViridis) # parameters for Oxygen colormap


cm <- colormap(all_s[["oxygen2"]], breaks=seq(all_ox_min, all_ox_max, 0.25))

all_p <- all_s[['pressure']]

pdf(file="ocean_section_plots/graphics/POSTER_temperature_salinity_plot.pdf", width=7, height=6)
#drawPalette(colormap=cm, zlab=expression(paste("Oxygen [mL L"^"-1","]")))
drawPalette(colormap=cm, zlab=expression(paste("Oxygen [mL/L]")))
plotTS(all_s, inSitu=TRUE, add=TRUE, pch=16, col=cm$zcol, type="n", Slim = c(33.5, 35), referencePressure = all_p, eos = "gsw", drawFreezing = TRUE, mar=par("mar")+c(0, 0, 0, .2))
plotTS(all_s, inSitu=TRUE, pch=16, col=cm$zcol, type="p", Slim = c(33.5, 35), referencePressure = all_p, eos = "gsw", drawFreezing = TRUE, mar=par("mar")+c(0, 0, 0, .2))
#lines(Sp, Tp_values, col="gray")
dev.off()


# Define the prefix
prefix <- "d:\\data\\raw\\nbp2202_"

# Define target identifiers based on the format of hexfilename
target_identifiers <- c("ctd005", "ctd008", "ctd018", "ctd022", "ctd033", "ctd296", "ctd041", 
                        "ctd075", "ctd126", "ctd139", "ctd155", "ctd167", "ctd184", "ctd194", 
                        "ctd220", "ctd235", "ctd245", "ctd258", "ctd267", "ctd271", "ctd279")

# Concatenate prefix with each target identifier
target_identifiers_with_prefix <- paste0(prefix, target_identifiers, ".hex")
target_identifiers <- target_identifiers_with_prefix

# Initialize a list to hold the filtered casts
filtered_casts <- list()

# Loop through each cast
for (i in seq_along(casts)) {
  # Extract hexfilename from the metadata, assuming it's an S4 object
  hexfilename <- casts[[i]]@metadata[["hexfilename"]]
  
  # Extract the significant part of the filename
  # Adjusted regex for Windows paths and correct extraction
  significant_part <- sub(".*\\\\(CTD\\d+)\\.hex$", "\\1", basename(hexfilename))
  
  # Convert the extracted part to lowercase
  significant_part <- tolower(significant_part)
  
  # Check if the extracted part is in the target identifiers list
  if (significant_part %in% target_identifiers) {
    filtered_casts[[length(filtered_casts) + 1]] <- casts[[i]]
  }
}

# Print or inspect the filtered casts
if (length(filtered_casts) > 0) {
  print(filtered_casts)
} else {
  print("No casts matched the criteria.")
}



library(ocedata) #for the coastlineWorldFine data
data(coastlineWorldFine)


c(-120, -109, -75, -71.5)

mp <- function() {
  mapPlot(coastlineWorldFine, projection="+proj=stere +lon_0=-90 +lat_0=90",
          longitudelim = c(-115, -111),
          latitudelim = c(-74.5, -71), col='grey')
}


library(marmap)
b <- as.topo(getNOAA.bathy(-180, 0, -55, -90, keep=TRUE))


cm <- colormap(all_s[["pressure"]], col1="#08316b", col0="#f7fbff", x0=)

cm <- colormap(x0=c(0, 50, 300, 500, 1000, 1500, 2000),
               x1=c(50, 300, 500, 1000, 1500, 2000, 3500),
               col0=c("#f7fbff","#e6f0f9","#c9dfee","#9ecae1", "#60a4cf", "#2f70a6", "#08316b"),
               col1=c("#e6f0f9","#c9dfee","#9ecae1","#60a4cf", "#2f70a6", "#08316b", "#08316b"), zlim=c(0,3500))


cm <- colormap(seq(-4000, 0, 500), col=oceColorsJet)
drawPalette(colormap=cm)

pdf(file="ocean_section_plots/graphics/station_alt.pdf", width=7, height=6)
drawPalette(colormap=cm, zlab = "Depth [m]")
mp()
mapImage(b, col=oceColorsGebco, breaks=seq(-4000, 0, 500))
mapPolygon(coastlineWorldFine, col='grey')
mapGrid() 
mapPoints(all, col="black", pch=25, cex=0.7, bg="red")
dev.off()


