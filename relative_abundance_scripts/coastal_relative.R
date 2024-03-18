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
ps_sub <- ps %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "Mitochondria" &
      Class   != "Chloroplast" &
      Order   != "Chloroplast" &
      Family  != "Mitochondria" 
  )

# free-living phyloseq
coastal_ps_free <- ps_sub %>% subset_samples(Filter_pores == "0.2") %>% prune_taxa(taxa_sums(.) > 0, .) 

# particle-associated phyloseq
coastal_ps_part <- ps_sub %>% subset_samples(Filter_pores >= "2") %>% prune_taxa(taxa_sums(.) > 0, .) 

###################### 
#  stacked barplots  #
######################
#    FREE-LIVING     # 
######################


# Create a data frame for freeliving
coastal_data_free <- coastal_ps_free %>% subset_samples(Coastal_Current_Name == "transect3") %>%
  tax_glom(taxrank = "Order") %>% # agglomerate at genus level
  transform_sample_counts(function(x) {x/sum(x)} )# Transform to rel. abundance

coastal_top_free <- top_taxa(coastal_data_free, 
                             n_taxa = 16,
                             include_na_taxa = T)

coastal_data_free <- coastal_top_free$ps_obj %>%
  psmelt() %>%  
  filter(., Coastal_Current_Number != "4") %>%  
  arrange(Coastal_Current_Number)           # Melt to long format

#fix out of order transect numbers, do it manually.
# works
coastal_data_free <- coastal_data_free %>% 
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "1", 0.0000)) %>%
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "2", 32.67378)) %>%
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "3", 48.56614)) %>%
  #mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "4", 61.17885)) %>% # REMOVED
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "5", 62.78023)) %>%
  #mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "6", 80.58525)) %>%  # REMOVED BECAUSE doesn't have free-living.
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "7", 101.50049)) %>%
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "8", 106.99770)) %>%
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "9", 133.80764))

# particle-associated
coastal_data_part <- coastal_ps_part %>% subset_samples(Coastal_Current_Name == "transect3") %>%
  tax_glom(taxrank = "Order") %>% # agglomerate at genus level
  transform_sample_counts(function(x) {x/sum(x)} )# Transform to rel. abundance

coastal_top_part <- top_taxa(coastal_data_part, 
                             n_taxa = 16,
                             include_na_taxa = T)

coastal_data_part <- coastal_top_part$ps_obj %>%
  psmelt() %>%  
  filter(., Coastal_Current_Number != "4") %>%  
  arrange(Coastal_Current_Number)     

coastal_data_part <- coastal_data_part %>% 
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "1", 0.0000)) %>%
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "2", 32.67378)) %>%
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "3", 48.56614)) %>%
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "4", 61.17885)) %>%
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "5", 62.78023)) %>%
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "6", 80.58525)) %>%
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "7", 101.50049)) %>%
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "8", 106.99770)) %>%
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "9", 133.80764))

coastal_plot_labels <- c("89", "132", "106", "14", "78", "56b", "68", "146") # DIFFERENT BC OF FREE-LIVING
coastal_plot_breaks <- unique(coastal_data_part$Coastal_Current_Number) # HAVE TO CHANGE
coastal_sec_labels <- seq(0 , 140, by=20)
coastal_sec_breaks <- seq(0 , 140, by=20)

myColors <- c(brewer.pal(9, "Paired"), "#A43D27", "#497687", "#5E4987", "darkblue", "lightblue2", "darkgoldenrod", "dodgerblue", "seagreen")
coastal_data_free$Order <- as.factor(coastal_data_free$Order)
coastal_data_part$Order <- as.factor(coastal_data_part$Order)
names(myColors) <- levels(c(coastal_data_free$Order, coastal_data_part$Order))
custom_colors <- scale_colour_manual(name = "Order", values = myColors)
# Melt to long format

# Plot 
coastal_barplot_free <- ggplot(coastal_data_free, aes(x = Coastal_Current_Number, y = Abundance, fill = Order)) +
  geom_bar(stat = "identity", position="fill", width=6) + theme_classic() + ggtitle("\n Free-living (<0.2 µm)") +
  #geom_col(position = "dodge") + # changes to multiple bars
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = myColors, drop = FALSE) +
  scale_x_continuous(
    name = "Distance (km)",
    breaks = coastal_sec_breaks,
    labels = coastal_sec_labels,
    expand = c(0,0),
    sec.axis = dup_axis(
      name = "",
      labels = coastal_plot_labels,
      breaks = coastal_plot_breaks)
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
ggsave("graphics/coastal_order_rel_abundance_free.pdf", width = 6.5, height = 4, dpi = 150)

###################### 
#  stacked barplots  #
######################
#      PARTICLE      # 
######################

# Plot 
coastal_barplot_part <- ggplot(coastal_data_part, aes(x = Coastal_Current_Number, y = Abundance, fill = Order)) +
  geom_bar(stat = "identity", position="fill", width=6) + theme_classic() + ggtitle("\n Particle-associated (>2 µm)") +
  # geom_col(position = "dodge") + # changes to multiple bars
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = myColors) + # set manual colors
  scale_x_continuous(
    name = "Distance (km)",
    breaks = coastal_sec_breaks,
    labels = coastal_sec_labels,
    expand = c(0,0),
    sec.axis = dup_axis(
      name = "",
      labels = coastal_plot_labels,
      breaks = coastal_plot_breaks)
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
ggsave("graphics/coastal_order_rel_abundance_part.pdf", width = 6.5, height = 4, dpi = 150)

# combined plot

total <- rbind(coastal_data_part, coastal_data_free)
# make combined FAKE plot to grab legend from and to put in the comine plot :^)
legend_plot <- ggplot(total, aes(x = Coastal_Current_Number, y = Abundance, fill = Order)) +
  geom_bar(stat = "identity", position="fill", width=2) + theme_classic() +
  # geom_col(position = "dodge") + # changes to multiple bars
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = myColors) 

legend_combined <- get_legend(legend_plot)

coastal_combined <- ggarrange(
  coastal_barplot_free, coastal_barplot_part, labels = NULL,
  common.legend = FALSE, legend = "right", legend.grob = legend_combined
)

annotate_figure(coastal_combined, top = text_grob("\n Coastal Current (Transect 3)", 
                                                  color = "dodgerblue3", face = "bold", size = 18))

ggsave("graphics/coastal_order_combined_relative.pdf", width = 13, height = 7, dpi = 150)