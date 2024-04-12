install.packages("devtools")
install.packages("ggTS")

no
# load libaries
install.packages("ggplot2")
install.packages("devtools")

library(ggplot2)
library(devtools)

# read in metadata
CTD <- read.table("CTD_meta.csv", sep=",", header=TRUE)
metadata <- read.table("artemis-metadata.csv", sep=",", header=TRUE)

install.packages("dplyr")
library(dplyr)

dat2 <- CTD %>% select(Station, Salinity, Temperature) %>%
  group_by(Station) %>%
  do(as.data.frame(spline(x= .[["Salinity"]], y= .[["Temperature"]], n = nrow(.)*10)))

ggplot(CTD, aes(Salinity, Temperature)) +
  geom_point(aes(color = factor(Station))) +
  geom_smooth(data = dat2, aes(x = x, y = y, color = factor(Station)))

ggplot(CTD, aes(x=Salinity,y=Temperature, col=Station))+geom_point()

# Generalized T-S plot
# https://hafezahmad.medium.com/making-temperature-salinity-diagrams-called-the-t-s-diagram-with-python-and-r-programming-5deec6378a29
library(marelac)
library(plot3D)
data <- read.table("CTD_meta.csv", sep=",", header=TRUE)
ts<-data[,c(5,10)]
mint=min(ts['Temperature'])
maxt=max(ts['Temperature'])
mins=min(ts['Salinity'])
maxs=max(ts['Salinity'])
salC<-seq(from=mins,to=maxs,length.out = 156)
tempC<-seq(from=mint,to=maxt,length.out = 156)
sigma.c<-outer(salC,tempC,FUN = function(S,t)sw_dens(S = S, t = t)-1000)
sigma.c

png(file = 'ts_diagram.png',width = 15,res=500,pointsize = 12,bg='white')

jpeg(file="ts_diagram.jpeg")

par(mar=c(5,5,4,6))
contour2D(x=salC,y=tempC,z=sigma.c,lwd=2,main='ARTEMIS (DNA Stations) TS Diagram',
          col='black',xlab=expression('Salinity'),ylab=expression('Temperature("~degree*C")'))
temp<-unlist(ts['Temperature'],use.names = FALSE)
sal<-unlist(ts['Salinity'],use.names = FALSE)
sigma_theta<-sw_dens(S=sal,t=temp)-1000
scatter2D(sal,temp,colvar = sigma_theta,pch=16,cex=1.25,add=TRUE,
          clim = range(sigma.c),colkey = FALSE)
colkey(clim = range(sigma.c),dist = 0.005,side=4,add=TRUE,
       clab = expression('Density(kg m"^-3")'),col.clab = 'black',
       side.clab = 4,line.clab = 2.5,length = 1,width = 0.8,
       col.axis = 'black',col.ticks = 'black',cex.axis = 0.9)
dev.off()


library(oce)
TS <- as.ctd(CTD$Salinity, CTD$Temperature, CTD$PrDM, CTD$Oxygen)
plotTS(TS, pch=16, col=colfunc(20))

library("ggplot2")

colfunc <- colorRampPalette(c("red", "blue"))

plot(1:20, 1:20, pch = 19, cex=2, col = colfunc(20))

library(oce)
install.packages("ocedata")
library(ocedata)

#importing ctd data from artemis files
ctd2 <- list.files(path="/Users/sierra/Documents/Research/Temperature_Salinity/Temp_Salinity_Plots/ctd/sv", pattern="\\.cnv$")
files = dir("./ctd/sv", full.names = TRUE, pattern = ".cnv") #works


station1 <- read.ctd.sbe("/Users/sierra/Documents/Research/Temperature_Salinity/Temp_Salinity_Plots/ctd/sv/SV_nbp2202_002.cnv")
station1%>%plotTS()


read.ctd.sbe(files)

# creates a bunch of plots
for(file in files){
  import <- read.ctd.sbe(file)
  print(plotTS(import, main=file))
}

###THIS SHIT WORKDS, had to transfer all the ctd files into main directory
files <- dir(pattern='*.cnv')
casts <- list()
for (ifile in 1:length(files)) {
  casts[[ifile]] <- read.ctd.sbe(files[ifile])
}

#structure
str(casts, 1)

T_all <- NULL
for (i in 1:length(casts)) {
  T_all <- c(T_all, casts[[i]][['temperature']])
}

T_all <- unlist(lapply(casts, function(x) x[['temperature']]))
S_all <- unlist(lapply(casts, function(x) x[['salinity']]))

plot(S_all, T_all, Slim=31.5,)
#plot to add my stations over the ctd data
pdf("CTD_TS_with_stations.pdf")
plotTS(casts, Slim=c(32.75,34.75),pch=16,col="dodgerblue3")
plotTS(TS, pch=16, cex=.75, col="black", add=TRUE)
dev.off()
#etc
CTD <- read.table("CTD_meta.csv", sep=",", header=TRUE)
TS <- as.ctd(CTD$Salinity, CTD$Temperature, CTD$PrDM)
plotTS(TS, pch=16, col="black")


data(section)
str(section@data$station, 1)


