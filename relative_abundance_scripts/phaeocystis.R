
ps_phaeo <- ps_noncontam_prev05 %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family == "Phaeocystis antarctica"
  )

ps_phaeo <- subset_samples(ps_phaeo, Sample.Control == "True.Sample")

ps_phaeo_rank <- ps_phaeo %>%
  tax_glom(taxrank = "Family")# agglomerate at Order level, can change to different taxonomic level!
# Transform to rel. abundance (normalize data

ps_phaeo_abun <- ps_phaeo_rank %>%
  psmelt() %>% # transform a phyloseq object into a data frame, otherwise graphs wont work
  filter(Abundance > 0.02) %>% # Filter out low abundance taxa
  arrange(Family)

ps_phaeo_agg <- aggregate(Abundance ~ Station * Family * More_Depth_Threshold * Latitude * Longitude, data = ps_phaeo_abun, FUN = mean)

ps_phaeo_abun_final <- ps_phaeo_agg %>% group_by(Station) %>% mutate(Relative = Abundance / sum(Abundance))

ps_phaeo_abun_mid <- filter(ps_phaeo_abun_final, More_Depth_Threshold == "Mid")

ps_phaeo_abun_surf <- filter(ps_phaeo_abun_final, More_Depth_Threshold == "Surface")

dt_surf <- select(ps_phaeo_abun_surf, c("Station", "Latitude", "Longitude", "Relative")) %>% na.omit(.) %>% distinct(., Station, .keep_all= TRUE) %>%
  rename(., c(lat = Latitude, lon = Longitude))

dt_surf$lat <- as.numeric(dt_surf$lat)
dt_surf$lon <- as.numeric(dt_surf$lon)


all_points <- basemap(data = dt, limits = c(-120, -107, -75, -71.5), bathymetry = FALSE, rotate = TRUE, glaciers = TRUE, legends = FALSE) 

surface <- all_points +
  ggspatial::geom_spatial_point(data = dt_surf, aes(x = lon, y = lat, fill=Relative), color = "black", pch=23, cex=3.5) + ggspatial::geom_spatial_text_repel(
    data = dt_surf, aes(x = lon, y = lat, label = Station), cex=2.5, vjust=-1.2) + scale_fill_gradient(low="dodgerblue", high="green", limits=c(0,1)) + ggtitle("Surface (>40m)") +
  theme(legend.position = "bottom",
        text = element_text(family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

dt_mid <- select(ps_phaeo_abun_mid, c("Station", "Latitude", "Longitude", "Relative")) %>% na.omit(.) %>% distinct(., Station, .keep_all= TRUE) %>%
  rename(., c(lat = Latitude, lon = Longitude))

dt_mid$lat <- as.numeric(dt_mid$lat)
dt_mid$lon <- as.numeric(dt_mid$lon)

all_points <- basemap(data = dt, limits = c(-120, -107, -75, -71.5), bathymetry = FALSE, rotate = TRUE, glaciers = TRUE, legends = FALSE) 

mid_thing <- all_points +
  ggspatial::geom_spatial_point(data = dt_mid, aes(x = lon, y = lat, fill=Relative), color = "black", pch=23, cex=3.5) + ggspatial::geom_spatial_text_repel(
    data = dt_mid, aes(x = lon, y = lat, label = Station), cex=2.5, vjust=-1.2) + scale_fill_gradient(low="dodgerblue", high="green", limits=c(0,1)) + ggtitle("Mid (190m - 475m)") +
  theme(legend.position = "bottom",
        text = element_text(family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


ps_combined <- ggarrange(
  surface, mid_thing, labels = NULL,
  common.legend = TRUE, legend = "right"
)

ggsave("phaeo_proxy_map.pdf", width = 11, height = 15)

#date
metadata <- filter(metatable, Sample.Control == "True.Sample")# filter by transect

dt <- select(metadata, c("Station", "Latitude", "Longitude", "Date")) %>% na.omit(.) %>% distinct(., Station, .keep_all= TRUE) %>%
  rename(., c(lat = Latitude, lon = Longitude))

dt$lat <- as.numeric(dt$lat)
dt$lon <- as.numeric(dt$lon)
dt$Date = as.Date(dt$Date, format = "%m/%d/%y")

dt <- dt %>%
  arrange(Date)

date <- all_points +
  ggspatial::geom_spatial_point(data = dt, aes(x = lon, y = lat, fill=Date, group=Date), pch=23, cex=3.5) + ggspatial::geom_spatial_text_repel(
    data = dt, aes(x = lon, y = lat, label = Station), cex=2.5, vjust=-1.2) + scale_fill_viridis_c(trans="date") + ggtitle("Stations with Date") +
  #scale_fill_gradient2(trans="date", low="yellow", high="purple")
  guides(fill = guide_colorbar(reverse=TRUE)) +
  theme(legend.position = "right",
        text = element_text(family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
ggsave("date_map.pdf", width = 10, height = 8)