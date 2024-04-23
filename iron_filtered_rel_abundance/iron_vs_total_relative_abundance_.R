library("dplyr")
library("tidyr")
library("phyloseq")
library("qiime2R")
library("ggplot2")
library("vegan")
library("plyr")
library("fantaxtic")
library("ggpubr")
library(tidyverse)
library(RColorBrewer)


# Extract abundance data from phyloseq object
all_ps_free <- ps_noncontam_prev05 %>% subset_samples(Filter_pores == "0.2") %>% subset_samples(Sample.Control == "True.Sample")
all_ps_part <- ps_noncontam_prev05 %>% subset_samples(Filter_pores >= "2") %>% subset_samples(Sample.Control == "True.Sample")

iron_ps_free <- filtered_ps %>% subset_samples(Filter_pores == "0.2") %>% prune_taxa(taxa_sums(.) > 0, .)
iron_ps_part <- filtered_ps %>% subset_samples(Filter_pores >= "2") %>% prune_taxa(taxa_sums(.) > 0, .) 

level_order_prev <- c("STN198", "STN002", "STN004", "STN181", "STN012", "STN115", "STN12.3", "STN20", "STN014", "STN089",
                 "STN132", "STN106", "STN078", "STN056a", "STN056b", "STN22", "STN068", "STN146", "STN174",
                 "STN151.2", "STN153") # in order, ish.

level_order <- c("STN089", "STN132", "STN106", "STN20", "STN198", "STN002", "STN004", "STN012", "STN115", "STN12.3", "STN014", "STN078", "STN056a", "STN056b", "STN22", "STN068", "STN146", "STN181", "STN174", "STN151.2", "STN153")


###############
# free-living #
###############
abundance_selected_free <- phyloseq::otu_table(iron_ps_free)
df_abundance_selected_free <- as.data.frame(abundance_selected_free)
selected_taxa_names_free <- rownames(df_abundance_selected_free) # List of specific taxa names, that i want to exclude

# Extract abundance data from phyloseq objects
abundance_all_free <- t(as.data.frame(otu_table(all_ps_free)))
# Extract transect numbers from metadata (assuming it's available in both data frames)
station_all_free <- as.factor(sample_data(all_ps_free)$Station)
depth_threshold_all_free <- as.factor(sample_data(all_ps_free)$More_Depth_Threshold)

# Call the function to calculate relative abundance based on Station and Depth_Threshold
result_ps_free <- calculate_relative_abundance(abundance_all_free, station_all_free, depth_threshold_all_free, selected_taxa_names_free)

# Access the debug data frame (debug information)
debug_df <- result_ps_free$Debug

# Access the debug data frame with relative abundance values
relative_result_free <- result_ps_free$Debug

mean_per_station_free <- aggregate(Abundance ~ Station * More_Depth_Threshold * data_source, data = relative_result_free, FUN = mean)

fill_free <- ggplot(mean_per_station_free, aes(x = factor(Station, level = level_order), y = Abundance, fill = data_source)) + 
  facet_grid(~factor(More_Depth_Threshold, levels=c("Surface", "Intermediate", "Bottom"))~.) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.3, width=0.9) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Transect Number", y = "Relative Abundance (%)", title = "Free-living") +
  scale_fill_manual(values = c("other" = "gray", "selected" = "salmon"),
                    labels = c("All other taxa", "Iron-related taxa"), 
                    name = "Taxa") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=9, angle=90),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.spacing.x = unit(0, "points"), # Reducing space between facets
        strip.background = element_blank(), # Optionally hide the strip background for a cleaner look
        panel.border = element_blank()) + # Optionally remove panel borders
  geom_vline(xintercept = c(4.5,11.5), linetype = "dashed", linewidth=0.6, color = "black")
ggsave("iron_filtered_rel_abundance/graphics/filtered_free_total_selected.pdf", width = 8, height = 4, dpi = 150)


###############
# particle    #
###############
abundance_selected_part <- phyloseq::otu_table(iron_ps_part)
df_abundance_selected_part <- as.data.frame(abundance_selected_part)
selected_taxa_names_part <- rownames(df_abundance_selected_part) # List of specific taxa names, that i want to exclude

# Extract abundance data from phyloseq objects
abundance_all_part <- t(as.data.frame(otu_table(all_ps_part)))
# Extract transect numbers from metadata (assuming it's available in both data frames)
station_all_part <- as.factor(sample_data(all_ps_part)$Station)
depth_threshold_all_part <- as.factor(sample_data(all_ps_part)$More_Depth_Threshold)

# Call the function to calculate relative abundance based on Station and Depth_Threshold
result_ps_part <- calculate_relative_abundance(abundance_all_part, station_all_part, depth_threshold_all_part, selected_taxa_names_part)

# Access the debug data frame (debug information)
debug_df <- result_ps_part$Debug

# Access the debug data frame with relative abundance values
relative_result_part <- result_ps_part$Debug
mean_per_station_part <- aggregate(Abundance ~ Station * More_Depth_Threshold * data_source, data = relative_result_part, FUN = mean)

fill_part <- ggplot(mean_per_station_part, aes(x = factor(Station, level = level_order), y = Abundance, fill = data_source)) + 
  facet_grid(~factor(More_Depth_Threshold, levels=c("Surface", "Intermediate", "Bottom"))~.) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.3, width=0.9) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Transect Number", y = "Relative Abundance (%)", title = "Particle-associated") +
  scale_fill_manual(values = c("other" = "gray", "selected" = "salmon"),
                    labels = c("All other taxa", "Iron-related taxa"), 
                    name = "Taxa") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=9, angle=90),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=9 , color="white"),
        legend.position = "none",
        panel.spacing.x = unit(0, "points"), # Reducing space between facets
        strip.background = element_blank(), # Optionally hide the strip background for a cleaner look
        panel.border = element_blank()) + # Optionally remove panel borders
  geom_vline(xintercept = c(4.5,11.5), linetype = "dashed", linewidth=0.6, color = "black") +# Add vertical lines
  #theme(plot.title = element_text(hjust = 0.5, size=17)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ggtitle("Particle-associated (>3 µm)")# Rotate x-axis labels for better readability
ggsave("iron_filtered_rel_abundance/graphics/filtered_free_total_selected.pdf", width = 8, height = 4, dpi = 150)

###################### 
#  stacked barplots  #
######################
#  BOTH COMMUNITIES  # 
######################

# combines the two graphs together
ps_combined <- ggarrange(
  fill_free, fill_part, labels = NULL,
  common.legend = TRUE, legend = "right"
)
annotate_figure(ps_combined, top = text_grob("Relative Abundance of Iron/Non-Iron Taxa for All Stations", 
                                             color = "black", face = "bold", size = 18))

ggsave("iron_filtered_rel_abundance/graphics/combined_barplot_diff_thresholds.pdf", width = 13, height = 7, dpi = 150)
