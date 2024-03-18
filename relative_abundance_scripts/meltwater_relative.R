# libraries
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

# select only bacteria, remove chloroplasts
ps_sub <- ps %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "Mitochondria" &
      Class   != "Chloroplast" &
      Order   != "Chloroplast" &
      Family  != "Mitochondria" 
  )

# free-living phyloseq
meltwater_ps_free <- ps_sub %>% subset_samples(Filter_pores == "0.2") %>% prune_taxa(taxa_sums(.) > 0, .)

# particle-associated phyloseq
meltwater_ps_part <- ps_sub %>% subset_samples(Filter_pores >= "2") %>% prune_taxa(taxa_sums(.) > 0, .) 

###################### 
#  stacked barplots  #
######################
#    FREE-LIVING     # 
######################

# Create a data frame for freeliving
meltwater_data_free <- meltwater_ps_free %>% subset_samples(Transect_Name == "transect2") %>%
  subset_samples(., Station != "STN22") %>% # Filter out low abundance taxa
  tax_glom(taxrank = "Order") %>% # agglomerate at genus level
  transform_sample_counts(function(x) {x/sum(x)} )# Transform to rel. abundance

meltwater_top_free <- top_taxa(meltwater_data_free, 
                               n_taxa = 16,
                               include_na_taxa = T)

meltwater_data_free <- meltwater_top_free$ps_obj %>%
  psmelt() %>%  
  arrange(Transect_Number)           # Melt to long format

myColors <- c(brewer.pal(9, "Paired"), "#A43D27", "#497687", "#5E4987", "darkblue", "lightblue2", "darkgoldenrod", "dodgerblue", "seagreen")
meltwater_data_free$Order <- as.factor(meltwater_data_free$Order)
meltwater_data_part$Order <- as.factor(meltwater_data_part$Order) # HAVE TO RUN data_part
names(myColors) <- levels(c(meltwater_data_free$Order, meltwater_data_part$Order))
custom_colors <- scale_colour_manual(name = "Order", values = myColors)

#fix out of order transect numbers, do it manually.
# works
meltwater_data_free <- meltwater_data_free %>% 
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "1", 0.00000000)) %>%
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "2", 16.33216538)) %>%
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "3", 40.58487518))

meltwater_plot_labels <- c("56b", "68", "146")
meltwater_plot_breaks <- unique(sample_data(meltwater_data_free)$Transect_Number)
meltwater_sec_labels <- seq(0 , 42, by=10)
meltwater_sec_breaks <- seq(0 , 42, by=10)

# Plot 
meltwater_barplot_free <- ggplot(meltwater_data_free, aes(x = Transect_Number, y = Abundance, fill = Order)) +
  geom_bar(stat = "identity", position="fill", width=2) + theme_classic() + ggtitle("\n Free-living (<0.2 µm)") +
  #geom_col(position = "dodge") + # changes to multiple bars
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = myColors, drop = FALSE) +
  scale_x_continuous(
    name = "Distance (km)",
    breaks = meltwater_sec_breaks,
    labels = meltwater_sec_labels,
    expand = c(0,0),
    sec.axis = dup_axis(
      name = "",
      labels = meltwater_plot_labels,
      breaks = meltwater_plot_breaks)
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
ggsave("graphics/meltwater_order_rel_abundance_free.pdf", width = 6.5, height = 4, dpi = 150)

###################### 
#  stacked barplots  #
######################
#      PARTICLE      # 
######################

meltwater_data_part <- meltwater_ps_part %>% subset_samples(Transect_Name == "transect2") %>%
  tax_glom(taxrank = "Order") %>% # agglomerate at genus level
  transform_sample_counts(function(x) {x/sum(x)} )# Transform to rel. abundance

meltwater_top_part <- top_taxa(meltwater_data_part, 
                               n_taxa = 16,
                               include_na_taxa = T)

meltwater_data_part <- meltwater_top_part$ps_obj %>%
  psmelt() %>%  
  filter(., Station != "STN22") %>%# Filter out low abundance taxa
  arrange(Transect_Number)     


meltwater_data_part <- meltwater_data_part %>% 
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "1", 0.00000000)) %>%
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "2", 16.33216538)) %>%
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "3", 40.58487518))


# Plot 
meltwater_barplot_part <- ggplot(meltwater_data_part, aes(x = Transect_Number, y = Abundance, fill = Order)) +
  geom_bar(stat = "identity", position="fill", width=2) + theme_classic() + ggtitle("\n Particle-associated (>2 µm)") +
  # geom_col(position = "dodge") + # changes to multiple bars
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = myColors) + # set manual colors
  scale_x_continuous(
    name = "Distance (km)",
    breaks = meltwater_sec_breaks,
    labels = meltwater_sec_labels,
    expand = c(0,0),
    sec.axis = dup_axis(
      name = "",
      labels = meltwater_plot_labels,
      breaks = meltwater_plot_breaks)
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
ggsave("graphics/meltwater_order_rel_abundance_part.pdf", width = 6.5, height = 4, dpi = 150)


# combined plot

total <- rbind(meltwater_data_part, meltwater_data_free)
# make combined FAKE plot to grab legend from and to put in the comine plot :^)
legend_plot <- ggplot(total, aes(x = Transect_Number, y = Abundance, fill = Order)) +
  geom_bar(stat = "identity", position="fill", width=2) + theme_classic() +
  # geom_col(position = "dodge") + # changes to multiple bars
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = myColors) 

legend_combined <- get_legend(legend_plot)

meltwater_combined <- ggarrange(
  meltwater_barplot_free, meltwater_barplot_part, labels = NULL,
  common.legend = FALSE, legend = "right", legend.grob = legend_combined
)

annotate_figure(meltwater_combined, top = text_grob("\n Meltwater Plume (Transect 2)", 
                                                    color = "dodgerblue3", face = "bold", size = 18))

ggsave("graphics/meltwater_order_combined_relative.pdf", width = 13, height = 7, dpi = 150)

