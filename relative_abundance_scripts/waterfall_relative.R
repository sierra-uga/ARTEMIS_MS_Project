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

##################
#  data culling  #
##################

# select only bacteria, remove chloroplasts
ps_sub <- ps_noncontam_prev05 %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "Mitochondria" &
      Class   != "Chloroplast" &
      Order   != "Chloroplast" &
      Family  != "Mitochondria" 
  )

# free-living phyloseq
waterfall_ps_free <- ps_sub %>% subset_samples(Filter_pores == "0.2") %>% prune_taxa(taxa_sums(.) > 0, .) 

# particle-associated phyloseq
waterfall_ps_part <- ps_sub %>% subset_samples(Filter_pores >= "2") %>% prune_taxa(taxa_sums(.) > 0, .) 

###################### 
#  stacked barplots  #
######################
#    FREE-LIVING     # 
######################

# Create a data frame for freeliving
waterfall_data_free <- waterfall_ps_free %>% subset_samples(Transect_Name == "transect1") %>%
  tax_glom(taxrank = "Order") %>% # agglomerate at genus level
  transform_sample_counts(function(x) {x/sum(x)} )# Transform to rel. abundance

waterfall_top_free <- top_taxa(waterfall_data_free, 
                               n_taxa = 11,
                               include_na_taxa = T)

waterfall_data_free <- waterfall_top_free$ps_obj %>%
  psmelt() %>%  
  arrange(Transect_Number)           # Melt to long format

#fix out of order transect numbers, do it manually.
# works
waterfall_data_free <- waterfall_data_free %>% 
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "1", 0.0000)) %>%
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "2", 62.80497)) %>%
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "3", 221.81425)) %>%
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "4", 298.40045)) %>%
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "5", 343.60091))

# particle-associated
waterfall_data_part <- waterfall_ps_part %>% subset_samples(Transect_Name == "transect1") %>%
  tax_glom(taxrank = "Order") %>% # agglomerate at genus level
  transform_sample_counts(function(x) {x/sum(x)} )# Transform to rel. abundance

waterfall_top_part <- top_taxa(waterfall_data_part, 
                               n_taxa = 11,
                               include_na_taxa = T)

waterfall_data_part <- waterfall_top_part$ps_obj %>%
  psmelt() %>%  
  arrange(Transect_Number)     

waterfall_data_part <- waterfall_data_part %>% 
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "1", 0.0000)) %>%
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "2", 62.80497)) %>%
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "3", 221.81425)) %>%
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "4", 298.40045)) %>%
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "5", 343.60091))


waterfall_plot_labels <- c("STN198","STN002", "STN004", "STN012", "STN014")
waterfall_plot_breaks <- unique(sample_data(waterfall_data_free)$Transect_Number)
waterfall_sec_labels <- seq(0 , 350.60, by=50)
waterfall_sec_breaks <- seq(0 , 350.60, by=50)

myColors <- c(brewer.pal(9, "Paired"), "#A43D27", "#497687", "#5E4987", "darkblue", "lightblue2", "darkgoldenrod", "dodgerblue", "seagreen")
waterfall_data_free$Order <- as.factor(waterfall_data_free$Order)
waterfall_data_part$Order <- as.factor(waterfall_data_part$Order)
names(myColors) <- levels(c(waterfall_data_free$Order, waterfall_data_part$Order))
custom_colors <- scale_colour_manual(name = "Order", values = myColors)
# Melt to long format

# Plot 
waterfall_barplot_free <- ggplot(waterfall_data_free, aes(x = Transect_Number, y = Abundance, fill = Order)) +
  geom_bar(stat = "identity", position="fill", width=15) + theme_classic() + ggtitle("\n Free-living (<0.2 µm)") +
  #geom_col(position = "dodge") + # changes to multiple bars
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = myColors, drop = FALSE) +
  scale_x_continuous(
    name = "Distance (km)",
    breaks = waterfall_sec_breaks,
    labels = waterfall_sec_labels,
    expand = c(0,0),
    sec.axis = dup_axis(
      name = "",
      labels = waterfall_plot_labels,
      breaks = waterfall_plot_breaks)
  ) +
  theme(plot.title = element_text(hjust = 0.5, size=14)) +
  theme(axis.title.x = element_text()) + # remove x title
  theme(axis.text.y = element_text()) + # remove y text
  theme(axis.title.y = element_text()) + # remove y title
  theme(legend.position = "right") + # position legent
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  theme(panel.spacing.y = unit(1, "lines")) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance")
ggsave("graphics/waterfall_order_rel_abundance_free.pdf", width = 6.5, height = 4, dpi = 150)

###################### 
#  stacked barplots  #
######################
#      PARTICLE      # 
######################

# Plot 
waterfall_barplot_part <- ggplot(waterfall_data_part, aes(x = Transect_Number, y = Abundance, fill = Order)) +
  geom_bar(stat = "identity", position="fill", width=15) + theme_classic() + ggtitle("\n Particle-associated (>2 µm)") +
  # geom_col(position = "dodge") + # changes to multiple bars
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = myColors) + # set manual colors
  scale_x_continuous(
    name = "Distance (km)",
    breaks = waterfall_sec_breaks,
    labels = waterfall_sec_labels,
    expand = c(0,0),
    sec.axis = dup_axis(
      name = "",
      labels = waterfall_plot_labels,
      breaks = waterfall_plot_breaks)
  ) +
  theme(plot.title = element_text(hjust = 0.5, size=14)) +
  theme(axis.title.x = element_text()) + # remove x title
  theme(axis.text.y = element_text()) + # remove y text
  theme(axis.title.y = element_text()) + # remove y title
  theme(legend.position = "right") + # position legent
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  theme(panel.spacing.y = unit(1, "lines")) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("")
ggsave("graphics/waterfall_order_rel_abundance_part.pdf", width = 6.5, height = 4, dpi = 150)

# combined plot

total <- rbind(waterfall_data_part, waterfall_data_free)
# make combined FAKE plot to grab legend from and to put in the comine plot :^)
legend_plot <- ggplot(total, aes(x = Transect_Number, y = Abundance, fill = Order)) +
  geom_bar(stat = "identity", position="fill", width=2) + theme_classic() +
  # geom_col(position = "dodge") + # changes to multiple bars
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = myColors) 

legend_combined <- get_legend(legend_plot)

waterfall_combined <- ggarrange(
  waterfall_barplot_free, waterfall_barplot_part, labels = NULL,
  common.legend = FALSE, legend = "right", legend.grob = legend_combined
)

annotate_figure(waterfall_combined, top = text_grob("\n CDW Waterfall (Transect 1)", 
                                                    color = "dodgerblue3", face = "bold", size = 18))

ggsave("graphics/decontam_waterfall_order_combined_relative.pdf", width = 13, height = 7, dpi = 150)