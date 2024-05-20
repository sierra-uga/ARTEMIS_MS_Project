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
library("speedyseq")

ps_sub <- ps_noncontam_prev05 %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "Mitochondria" &
      Class   != "Chloroplast" &
      Order   != "Chloroplast" &
      Family  != "Mitochondria" 
  )

ps_sub <- subset_samples(ps_sub, Sample.Control == "True.Sample")

temp <- sample_data(ps_sub)
temp$More_Depth_Threshold[temp$More_Depth_Threshold == "Mid-Bottom"] <- "T-min"
temp$More_Depth_Threshold[temp$More_Depth_Threshold == "Mid-Surface"] <- "Mixed Layer"
sample_data(ps_sub) <- temp

# free-living phyloseq
ps_free <- ps_sub %>% subset_samples(Filter_pores == "0.2") %>% prune_taxa(taxa_sums(.) > 0, .) 

# particle-associated phyloseq
ps_part <- ps_sub %>% subset_samples(Filter_pores >= "2") %>% prune_taxa(taxa_sums(.) > 0, .) 

#created a unique code based on Depth_Threshold, so took the MEAN of each depth threshold from each STATION.
# potentially frowned upon
ps0 <- merge_samples2(ps_free, "unique_depth",
                      fun_otu = mean
)

ps1 <- merge_samples2(ps_part, "unique_depth",
                      fun_otu = mean
)

##################
#  data culling  #
##################

# Create a data frame for freeliving, agglomerate by Order, transform to rel.abundance
data_free <- ps_free %>%
  tax_glom(taxrank = "Order") %>% # agglomerate at Order level, can change to different taxonomic level!
  transform_sample_counts(function(x) {x/sum(x)}) # Transform to rel. abundance (normalize data

data_top_free <- data_free %>%
  psmelt() %>% # transform a phyloseq object into a data frame, otherwise graphs wont work
  filter(Abundance > 0.02) %>% # Filter out low abundance taxa
  arrange(Order)

data_top_free <- aggregate(Abundance ~ Station * Order * More_Depth_Threshold, data = data_top_free, FUN = mean)


# particle-associated
data_part <- ps_part %>%
  tax_glom(taxrank = "Order") %>% # agglomerate at Order level
  transform_sample_counts(function(x) {x/sum(x)} )  # Transform to rel. abundance (normalize data)

data_top_part <- data_part %>%
  psmelt() %>% # transform a phyloseq object into a data frame, otherwise graphs wont work
  filter(Abundance > 0.02) %>%
  arrange(Order)# Filter out low abundance taxa

data_top_part <- aggregate(Abundance ~ Station * Order * More_Depth_Threshold, data = data_top_part, FUN = mean)

level_order_prev <- c("STN198", "STN002", "STN004", "STN181", "STN012", "STN115", "STN12.3", "STN20", "STN014", "STN089",
                 "STN132", "STN106", "STN078", "STN056a", "STN056b", "STN22", "STN068", "STN146", "STN174",
                 "STN151.2", "STN153") # in order, ish.

level_order_prev_inf <- c("STN198", "STN002", "STN004", "STN153", "STN151.2", "STN174", "STN181", "STN012", "STN115", "STN12.3", "STN20", "STN146","STN068", "STN056a", "STN056b", "STN22", "STN078", "STN014", "STN106", "STN132", "STN089")

level_order <- c("STN089", "STN132", "STN106", "STN20", "STN198", "STN002", "STN004", "STN012", "STN115", "STN12.3", "STN014", "STN078", "STN056a", "STN056b", "STN22", "STN068", "STN146", "STN181", "STN174", "STN151.2", "STN153")

###################
# Plot variables! #
###################

plot_labels <- unique(sample_data(ps_free)$Station) #
plot_breaks <- unique(sample_data(ps_free)$Station) # 


# Hex color codes (as vector)
hex_colors <- c(
  "Alteromonadales" = "#b33317",      # iron - red
  "Burkholderiales" = "lightblue1",      # organic/etc - blue
  "Cellvibrionales" = "#435a9e",      # organic/etc - blue
  "Chitinophagales" = "pink",      # organic - green
  "Ectothiorhodospirales" = "#5E4987",# purple sulfur - purple
  "Flavobacteriales" = "#6aa84f",     # organic - green
  "HOC36" = "#FFFF00",                # sulfur - yellow
  "Microtrichales" = "#d71414",       # iron - red
  "Nitrospinales" = "#f3be87",        # nitrate - orange
  "Oceanospirillales" = "#215909",    # organic - green
  "OM182 clade" = "darkblue",          # organic/etc - blue
  "Phycisphaerales" = "seagreen",      # organic - green
  "Pirellulales" = "#80cdc1",         # organic/etc - blue
  "Pseudomonadales" = "#4169E1",      # organic/etc - blue
  "Rhizobiales" = "#ec9438",          # nitrate - orange
  "Rhodobacterales" = "darkred",      # iron - red
  "Rhodospirillales" = "#a37fff",     # purple sulfur - purple
  "SAR11 clade" = "#6fa8dc",          # organic/etc - blue
  "SAR86 clade" = "#497687",          # organic/etc - blue
  "Sphingomonadales" = "#99c487",     # photosyn - light green
  "Thiomicrospirales" = "#edc31f",    # sulfur - yellow
  "Verrucomicrobiales" = "#01665e",   # organic - green
  "Vibrionales" = "#80b900",          # organic - green
  "Bacteriovoracales" = "#0b5394",    # organic/etc - blue
  "Planctomycetales" = "#cc6800",     # nitrate - orange
  "Sphingobacteriales" = "#008000"    # organic - green
)


myColors <- c(brewer.pal(9, "Paired"),'#e66101','darkgreen','#fdb863','#5e3c99', '#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','darkred','#c7eae5','#80cdc1','#35978f','#01665e','#003c30', "#A43D27", "#497687", "#5E4987", "darkgoldenrod", "lightblue2", "darkblue", "dodgerblue", "seagreen", "purple", "black")
# this must equal the levels of the Order
data_top_free$Order <- as.factor(data_top_free$Order) # setting the Order columns to factor
data_top_part$Order <- as.factor(data_top_part$Order) # setting the Order columns to factor
names(myColors) <- levels(c(data_top_free$Order, data_top_part$Order)) # setting the names of the colors to coordinate with the Order columns of each dataframe

###################### 
#  stacked barplots  #
######################
#    FREE-LIVING     # 
######################

barplot_free <- ggplot(data_top_free, aes(x = factor(Station, level = level_order), y = Abundance, fill = Order, group = Order)) + facet_grid(~factor(More_Depth_Threshold, levels=c("Surface", "Mixed Layer", "Mid", "T-min", "Bottom"))~.) + # facet grid seperates by different levels, horizontally
  geom_bar(stat = "identity", position="fill", color= "black", linewidth=0.3, width=0.9) + theme_classic() + # adds black outline to boxes
  scale_y_continuous(expand = c(0, 0)) + # extends the barplots to the axies
  scale_fill_manual(values = hex_colors, drop = FALSE) + # set the colors with custom colors (myColors)
  scale_x_discrete(
    breaks = plot_breaks, # setting breaks
    labels = plot_labels, # settting levels
    drop = FALSE
  ) +
  theme(text = element_text(family = "Helvetica"), 
        plot.title = element_text(hjust = 0.6),
    axis.title.x = element_blank(),
        axis.text=element_text(size=7),
        axis.text.x = element_text(size=9, angle=90, vjust=0.5),
        axis.title.y = element_blank(),
        legend.position = "none"
        # Reducing space between facets
        ) + # Optionally remove panel borders
  geom_vline(xintercept = c(4.5,11.5), linetype = "dashed", linewidth=0.6, color = "black") +# Add vertical lines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # for the legend, if you want one
  #ylab("Relative Abundance (Order > 2%) \n") + # remove # if you want y-axis title
  ggtitle("Free-living (<0.2 µm)")
ggsave("graphics/free_living_barplot_3.pdf", width = 8, height = 6, dpi = 150)

###################### 
#  stacked barplots  #
######################
#      PARTICLE      # 
######################

# the following plot is basically the same as above, look at annotation for free-living barplot if confused about what each line does!
barplot_part <- ggplot(data_top_part, aes(x = factor(Station, level = level_order), y = Abundance, fill = Order, group = Order)) + facet_grid(~factor(More_Depth_Threshold, levels=c("Surface", "Mixed Layer", "Mid", "T-min", "Bottom"))~., scales = "free_x", space = "free_x") +
  geom_bar(stat = "identity", position="fill", color= "black", linewidth=0.3, width=0.9) + theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = hex_colors, drop = FALSE) +
  scale_x_discrete(
    breaks = plot_breaks,
    labels = plot_labels,
    drop = FALSE
  ) +
  theme(text = element_text(family = "Helvetica"), 
        plot.title = element_text(hjust = 0.6),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=9, angle=90, vjust=0.5),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=5.5 , color="white"),
        legend.position = "none") +# Reducing space between facets
 # Optionally hide the strip background for a cleaner look
# Optionally remove panel borders
  geom_vline(xintercept = c(4.5,11.5), linetype = "dashed", linewidth=0.6, color = "black") +# Add vertical lines
  #theme(plot.title = element_text(hjust = 0.5, size=17)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ggtitle("Particle-associated (>3 µm)")

ggsave("graphics/part_associated_barplot_3.pdf", width = 8, height = 6, dpi = 150)

###################### 
#  stacked barplots  #
######################
#  BOTH COMMUNITIES  # 
######################

total <- rbind(data_top_part, data_top_free)
# make combined FAKE plot to grab legend from and to put in the combine plot :^)
legend_plot <- ggplot(total, aes(x = Station, y = Abundance, fill = Order)) +
  geom_bar(stat = "identity", position="fill", width=2) + theme_classic() +
  # geom_col(position = "dodge") + # changes to multiple bars
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = hex_colors) +
  guides(fill = guide_legend(override.aes = list(color = "black", size = 1))) # adds black outline around legend

# get legend from the fake combined plot
legend_combined <- get_legend(legend_plot)

# combines the two graphs together
ps_combined <- ggarrange(
  barplot_free, barplot_part, labels = NULL,
  common.legend = TRUE, legend = "right", legend.grob = legend_combined
)
# to remove the white space between the two community plots, you'd have to play with the plot.margins of each plot individually!
# something like: plot.margin=unit(c(1,1,-0.5,1), "cm")), where the margins follow the following structure:
# unit(c(top, right, bottom, left), units).

annotate_figure(ps_combined, top = text_grob("Total Relative Abundance for All ARTEMIS Stations", 
                                             color = "black", hjust=.7, face = "bold", size = 18, family = "Helvetica"))

ggsave("graphics/POSTER_combined_barplot_mid_thresholds_diff_order.pdf", width = 13, height = 7, dpi = 150)