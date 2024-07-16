
ps_sub <- ps_noncontam_prev05 %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "Mitochondria" &
      Family   != "Chloroplast" &
      Order   != "Chloroplast" &
      Family  != "Mitochondria" 
  )

ps_sub <- subset_samples(ps_sub, Sample.Control == "True.Sample") %>% subset_samples(Transect_Name == "transect1") %>% subset_samples(sample_name != "STN12.3.732.fil.poly.S.r1") %>%
  phyloseq_validate() %>% tax_fix() %>% prune_taxa(taxa_sums(.) > 0, .) 

ps_free_nmds <- ordinate(
  physeq = ps_sub, 
  method = "NMDS", 
  distance = "bray"
)

shape_breaks <- unique(sample_data(ps_sub)$Filter_pores)
color_breaks <- unique(sample_data(ps_sub)$Station)
myColors <- c("darkgreen", "lightgreen", brewer.pal(9, "Paired"), "#A43D27", "#497687", "#5E4987", "darkblue", "lightblue2", "darkgoldenrod", "dodgerblue", "seagreen", "red", "blue", "purple", "yellow")
names(myColors) <- levels(sample_data(ps_free)$Station)
custom_colors <- scale_colour_manual(name = "Station", values = myColors)

ps_free_nmds_plot <- plot_ordination(
  physeq = ps_sub,
  ordination = ps_free_nmds,
  color = "Station",
  shape = "Filter_pores",
  title = "Both Communities of CDW Waterfall Stations"
) +
  scale_color_manual(values = c(myColors, "darkgreen", "dodgerblue", "red2", "blueviolet", "aquamarine3", "gray"),
                     name = "Station",
                     breaks = color_breaks,
                     labels = color_breaks) +
  scale_shape_manual(values = c(16, 8, 1, 2),
                     name = "Community",
                     breaks = shape_breaks,
                     labels = shape_breaks) +
  geom_point(aes(shape = Filter_pores), alpha = 0.9, size = 3) +
annotate("text", x = -.25, y = -1, label ="2D Stress: 0.086")

#takes away grid from ggplot
ps_free_nmds_plot + theme(legend.position = "bottom") + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave("ordination_scripts/graphics/CDW_waterfall_NMS.pdf", width = 7, height = 5, dpi = 150)