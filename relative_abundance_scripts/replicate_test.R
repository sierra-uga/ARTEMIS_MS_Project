ps_sub <- subset_samples(ps_sub, Sample.Control == "True.Sample")
wanted_stations <- c("STN181", "STN174", "STN12.3", "STN151.2", "STN153", 
                     "STN198") # make a vector of desired stations
ps_sub <- ps_sub %>% subset_samples(!Station %in% wanted_stations) 

# free-living phyloseq
ps_free <- ps_sub %>% subset_samples(Filter_pores == "0.2") %>% prune_taxa(taxa_sums(.) > 0, .) 

# particle-associated phyloseq
ps_part <- ps_sub %>% subset_samples(Filter_pores >= "2") %>% prune_taxa(taxa_sums(.) > 0, .) 

# remove stations with only one replicate
data_free <- ps_free %>%
  tax_glom(taxrank = "Order") %>% # agglomerate at Order level, can change to different taxonomic level!
  transform_sample_counts(function(x) {x/sum(x)}) # Transform to rel. abundance (normalize data

data_top_free <- data_free %>%
  psmelt() %>% # transform a phyloseq object into a data frame, otherwise graphs wont work
  filter(Abundance > 0.02) %>%  # Filter out low abundance taxa
  arrange(watertype)# arrange by Order

# particle-associated
data_part <- ps_part %>%
  tax_glom(taxrank = "Order") %>% # agglomerate at Order level
  transform_sample_counts(function(x) {x/sum(x)} )  # Transform to rel. abundance (normalize data)

data_top_part <- data_part %>%
  psmelt() %>% # transform a phyloseq object into a data frame, otherwise graphs wont work
  filter(Abundance > 0.02) %>%   # Filter out low abundance taxa
  arrange(Order)  # arrange by Order

# plot breaks for graph, each Station (3 total)
myColors <- c(brewer.pal(9, "Paired"),'#e66101','#fdb863','#b2abd2','#5e3c99', '#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e','#003c30', "#A43D27", "#497687", "#5E4987", "darkgoldenrod", "lightblue2", "darkblue", "dodgerblue", "seagreen", "purple", "black") # this must equal the levels of the Order
data_top_free$Order <- as.factor(data_top_free$Order) # setting the Order columns to factor
data_top_part$Order <- as.factor(data_top_part$Order) # setting the Order columns to factor
names(myColors) <- levels(c(data_top_free$Order)) # setting the names of the colors to coordinate with the Order columns of each dataframe


barplot_free <- ggplot(data_top_free, aes(x = Depth_Threshold, y = Abundance, fill = Order, group = Order)) + facet_wrap(~factor(Station, levels=level_order)~.) + # facet grid seperates by different levels, horizontally
  geom_bar(stat = "identity", position="fill", color= "black", linewidth=0.3, width=0.9) + theme_classic() + # adds black outline to boxes
  scale_y_continuous(expand = c(0, 0)) + # extends the barplots to the axies
  scale_fill_manual(values = myColors, drop = FALSE) + # set the colors with custom colors (myColors)
  scale_x_discrete(
    breaks = plot_breaks, # setting breaks
    labels = plot_labels, # settting levels
    drop = FALSE
  ) +
  #theme(plot.title = element_text(hjust = 0.5, size=17)) + # remove # if you want title
  theme(axis.title.x = element_blank()) + # removing x-axis title
  theme(axis.text.x = element_text(size=9, angle=90)) + # setting x-axis title
  theme(axis.title.y = element_blank()) + # removing y-axis title
  theme(legend.position = "none") + # remove legend, delete line if you want a legend
  theme(panel.spacing.y = unit(1, "lines")) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # for the legend, if you want one
  #ylab("Relative Abundance (Order > 2%) \n") + # remove # if you want y-axis title
  ggtitle("Free-living (<0.2 µm)")


library(purrr)
library(cowplot)
# Define a function to create plots for each station
create_station_plot <- function(station_data) {
  ggplot(station_data, aes(x = sample.illumina, y = Abundance, fill = Order, group = Order)) +
    geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.3, width = 0.9) +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = myColors, drop = FALSE) +
    scale_x_discrete(
      drop = FALSE
    ) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 9, angle = 90, color = ifelse(endsWith(plot_labels, "_R2"), "red", "black")),
      axis.title.y = element_blank(),
      panel.spacing.y = unit(1, "lines"),
      legend.position = "none"
    ) +
    guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
    ggtitle(paste(unique(station_data$Station))) +
    ylab("Relative Abundance (Order > 2%) \n")
}

# Create a list of plots for each station
station_plots <- map(unique(data_top_free$Station), ~create_station_plot(filter(data_top_free, Station == .x)))

# Display the plots
plot_grid(plotlist = station_plots, ncol = 6)


# Define a function to create plots for each station
create_station_plot <- function(station_data) {
  ggplot(station_data, aes(x = sample.illumina, y = Abundance, fill = Order, group = Order)) +
    geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.3, width = 0.9) +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = myColors, drop = FALSE) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 9, angle = 90),
      axis.title.y = element_blank(),
      panel.spacing.y = unit(1, "lines"),
      legend.position = "none"
    ) +
    guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
    ggtitle(paste(unique(station_data$Station))) +
    ylab("Relative Abundance (Order > 2%) \n") +
    scale_x_discrete(labels = function(x) {
      # Define colors for sample names
      colors <- ifelse(endsWith(x, "_R2"), paste0(station_data$Depth[endsWith(x, "_R2")], "_ONE"), paste0(station_data$Depth[!endsWith(x, "_R2")], "_TWO"))
      # Create a named vector with colors
      names(colors) <- x
      colors
    })
}

# Create a list of plots for each station
station_plots <- map(unique(data_top_free$Station), ~create_station_plot(filter(data_top_free, Station == .x)))

# Arrange the plots in a grid
plot_grid(plotlist = station_plots, ncol = 1)



legend_plot <- ggplot(data_top_free, aes(x = Station, y = Abundance, fill = Order)) +
  geom_bar(stat = "identity", position="fill", width=2) + theme_classic() +
  # geom_col(position = "dodge") + # changes to multiple bars
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = myColors) +
  guides(fill = guide_legend(override.aes = list(color = "black", size = 1))) # adds black outline around legend


legend_combined <- get_legend(legend_plot)

# combines the two graphs together
ps_combined <- ggarrange(
  plotlist=station_plots, labels = NULL, ncol= 8, nrow=2,
  common.legend = TRUE, legend = "right", legend.grob = legend_combined
)
ggsave("graphics/replicate_test.pdf", width = 20, height = 10, dpi = 150)
