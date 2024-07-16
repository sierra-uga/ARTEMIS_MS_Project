ps_sub <- ps_noncontam_prev05 %>%
  subset_taxa(
    #Kingdom == "Bacteria" &
      Family  != "Mitochondria" &
      Family   != "Chloroplast" &
      Order   != "Chloroplast" &
      Family  != "Mitochondria" 
  )

ps_sub <- subset_samples(ps_sub, Sample.Control == "True.Sample") %>% subset_samples(Coastal_Current_Name == "transect3") %>% subset_samples(sample_name != "STN12.3.732.fil.poly.S.r1") %>%
  subset_samples(sample_name != "STN089.200.fil.dura.r2") %>% 
  phyloseq_validate() %>% tax_fix() %>% prune_taxa(taxa_sums(.) > 0, .) 

ps_free <- ps_sub %>% subset_samples(Filter_pores == "free-living") %>% 
  prune_taxa(taxa_sums(.) > 0, .) 

ps_part <- ps_sub %>% subset_samples(Filter_pores == "particle-associated") %>% 
  prune_taxa(taxa_sums(.) > 0, .) 


ps_free_nmds <- ordinate(
  physeq = ps_free, 
  method = "NMDS", 
  distance = "bray"
)

shape_breaks <- unique(sample_data(ps_free)$True_Flow)
color_breaks <- unique(sample_data(ps_free)$Location)
myColors <- c("dodgerblue", "red", "darkred", brewer.pal(9, "Paired"), "#A43D27", "#497687", "#5E4987", "darkblue", "lightblue2", "darkgoldenrod", "dodgerblue", "seagreen", "red", "blue", "purple", "yellow")
names(myColors) <- levels(sample_data(ps_free)$Location)
custom_colors <- scale_colour_manual(name = "Location", values = myColors)

ps_free_nmds_plot <- plot_ordination(
  physeq = ps_free,
  ordination = ps_free_nmds,
  color = "Location",
  shape = "True_Flow",
  title = "Free-living Coastal Current Stations"
) +
  scale_color_manual(values = c(myColors, "darkgreen", "dodgerblue", "red2", "blueviolet", "aquamarine3", "gray"),
                     name = "Location",
                     breaks = color_breaks,
                     labels = color_breaks) +
  scale_shape_manual(values = c(8, 8, 1, 2),
                     name = "Flow",
                     breaks = c("Outflow"),
                     labels = c("Outflow") )+
  geom_point(shape = 16, aes(color = Location), alpha = 0.9, size = 3) +
  geom_point(aes(shape = True_Flow), alpha = 0.9, size = 4) +
  annotate("text", x = -.25, y = -1, label ="2D Stress: 0.079")

#takes away grid from ggplot
ps_free_nmds_plot + theme(legend.position = "bottom") + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave("ordination_scripts/graphics/Coastal_current_NMS_Location_flow_free_living.pdf", width = 7, height = 5, dpi = 150)

#### PART


ps_part_nmds <- ordinate(
  physeq = ps_part, 
  method = "NMDS", 
  distance = "bray"
)

shape_breaks <- unique(sample_data(ps_part)$True_Flow)
color_breaks <- unique(sample_data(ps_part)$Location)
myColors <- c("dodgerblue", "red", "darkred", brewer.pal(9, "Paired"), "#A43D27", "#497687", "#5E4987", "darkblue", "lightblue2", "darkgoldenrod", "dodgerblue", "seagreen", "red", "blue", "purple", "yellow")
names(myColors) <- levels(sample_data(ps_part)$Location)
custom_colors <- scale_colour_manual(name = "Location", values = myColors)

ps_part_nmds_plot <- plot_ordination(
  physeq = ps_part,
  ordination = ps_part_nmds,
  color = "Location",
  shape = "True_Flow",
  title = "Particle-associated Coastal Current Stations"
) +
  scale_color_manual(values = c(myColors, "darkgreen", "dodgerblue", "red2", "blueviolet", "aquamarine3", "gray"),
                     name = "Location",
                     breaks = color_breaks,
                     labels = color_breaks) +
  scale_shape_manual(values = c(8, 8, 1, 2),
                     name = "Flow",
                     breaks = c("Outflow"),
                     labels = c("Outflow") )+
  geom_point(shape = 16, aes(color = Location), alpha = 0.9, size = 3) +
  geom_point(aes(shape = True_Flow), alpha = 0.9, size = 4) +
  annotate("text", x = -1, y = -1, label ="2D Stress: 0.079")

#takes away grid from ggplot
ps_part_nmds_plot + theme(legend.position = "bottom") + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave("ordination_scripts/graphics/Coastal_current_NMS_Location_flow_part_assoc.pdf", width = 7, height = 5, dpi = 150)

### BY STATION

ps_free_nmds <- ordinate(
  physeq = ps_free, 
  method = "NMDS", 
  distance = "bray"
)

shape_breaks <- unique(sample_data(ps_free)$True_Flow)
color_breaks <- unique(sample_data(ps_free)$Station)
myColors <- c(brewer.pal(9, "Paired"), "#A43D27", "#497687", "#5E4987", "darkblue", "lightblue2", "darkgoldenrod", "dodgerblue", "seagreen", "red", "blue", "purple", "yellow")
names(myColors) <- levels(sample_data(ps_free)$Station)

ps_free_nmds_plot <- plot_ordination(
  physeq = ps_free,
  ordination = ps_free_nmds,
  color = "Station",
  shape = "True_Flow",
  title = "Free-living Coastal Current Stations"
) +
  scale_color_manual(values = c(myColors, "darkgreen", "dodgerblue", "red2", "blueviolet", "aquamarine3", "gray"),
                     name = "Station",
                     breaks = color_breaks,
                     labels = color_breaks) +
  scale_shape_manual(values = c(8, 8, 1, 2),
                     name = "Flow",
                     breaks = c("Outflow"),
                     labels = c("Outflow") )+
  geom_point(shape = 16, aes(color = Station), alpha = 0.9, size = 3) +
  geom_point(aes(shape = True_Flow), alpha = 0.9, size = 4) +
  annotate("text", x = -.25, y = -1, label ="2D Stress: 0.079")

#takes away grid from ggplot
ps_free_nmds_plot + theme(legend.position = "bottom") + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave("ordination_scripts/graphics/Coastal_current_NMS_STATION_flow_free_living.pdf", width = 7, height = 5, dpi = 150)


## part 
ps_part_nmds <- ordinate(
  physeq = ps_part, 
  method = "NMDS", 
  distance = "bray"
)

shape_breaks <- unique(sample_data(ps_part)$True_Flow)
color_breaks <- unique(sample_data(ps_part)$Station)
myColors <- c(brewer.pal(9, "Paired"), "#A43D27", "#497687", "#5E4987", "darkblue", "lightblue2", "darkgoldenrod", "dodgerblue", "seagreen", "red", "blue", "purple", "yellow")
names(myColors) <- levels(sample_data(ps_part)$Station)

ps_part_nmds_plot <- plot_ordination(
  physeq = ps_part,
  ordination = ps_part_nmds,
  color = "Station",
  shape = "True_Flow",
  title = "Particle-associated Coastal Current Stations"
) +
  scale_color_manual(values = c(myColors, "darkgreen", "dodgerblue", "red2", "blueviolet", "aquamarine3", "gray"),
                     name = "Station",
                     breaks = color_breaks,
                     labels = color_breaks) +
  scale_shape_manual(values = c(8, 8, 1, 2),
                     name = "Flow",
                     breaks = c("Outflow"),
                     labels = c("Outflow") )+
  geom_point(shape = 16, aes(color = Station), alpha = 0.9, size = 3) +
  geom_point(aes(shape = True_Flow), alpha = 0.9, size = 4) +
  annotate("text", x = -.25, y = -1, label ="2D Stress: 0.079")

#takes away grid from ggplot
ps_part_nmds_plot + theme(legend.position = "bottom") + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave("ordination_scripts/graphics/Coastal_current_NMS_STATION_flow_part_assoc.pdf", width = 7, height = 5, dpi = 150)


