###############
library(cowplot)

ps_free <- subset_samples(ps_free, sample.illumina != "089_200_FIL_R2") 

# Identify the unique stations in your data
unique_stations <- unique(sample_data(ps_free)$Station)

# Create a list to store sub-phyloseq objects
sub_physeq_list <- list()

colors <- c("darkgreen", "dodgerblue", "red2", "blueviolet", "aquamarine3", "gray")
# Loop through each station and subset the phyloseq object
for (station in unique_stations) {
  sub_physeq <- subset_samples(ps_free, station == Station)
  sub_physeq_list[[station]] <- sub_physeq
}

# Create a list to store ordination results and plots
ordination_list <- list()
plot_list <- list()

# Loop through each sub-phyloseq object
for (station_name in unique_stations) {
  # Ordinate using NMDS
  ord_result <- ordinate(
    physeq = sub_physeq_list[[station_name]], 
    method = "PCoA", 
    distance = "bray"
  )
  
  # Plot ordination
  plot <- plot_ordination(
    physeq = sub_physeq_list[[station_name]],
    ordination = ord_result,
    color = "Depth",
    shape = "Depth") +
    ggtitle(paste(station_name)) +
    theme(legend.position = "bottom") + theme_bw() + 
    geom_point(size = 2) + # Adjust point size
    scale_color_manual(values = colors,
                       name = "Depth") +
    scale_shape_manual(values = c(21, 22, 23, 24, 25, 8),
                       name = "Depth") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  # Store results in lists
  ordination_list[[station_name]] <- ord_result
  plot_list[[station_name]] <- plot
}
plot_grid(plotlist = plot_list, ncol = 5)
ggsave("graphics/replicate_ordination_test.pdf", width = 20, height = 10, dpi = 150)