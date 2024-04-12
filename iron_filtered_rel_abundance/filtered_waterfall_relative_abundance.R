library(tidyverse)
library(RColorBrewer)

# free-living phyloseq
waterfall_ps_free <- filtered_ps %>% subset_samples(Filter_pores == "0.2") %>% prune_taxa(taxa_sums(.) > 0, .)

# particle-associated phyloseq
waterfall_ps_part <- filtered_ps %>% subset_samples(Filter_pores >= "2") %>% prune_taxa(taxa_sums(.) > 0, .) 

###################### 
#  stacked barplots  #
######################
#    FREE-LIVING     # 
######################

# Create a data frame for freeliving
waterfall_data_free <- waterfall_ps_free %>% subset_samples(Transect_Name == "transect1") %>%
  tax_glom(taxrank = "TAX")  %>% # agglomerate at tax level
  transform_sample_counts(function(x) {x/sum(x)} ) # Transform to rel. abundance

waterfall_data_free_temp <- waterfall_data_free

waterfall_data_free <- waterfall_data_free %>%
  psmelt() #%>%  # Melt to long format
#arrange(Transect_Number) %>%
#group_by(Transect_Number) %>% 
#mutate(Relative = Abundance/sum(Abundance))  #### USE THIS WHEN NOT USING THE METHOD BELOW. 
### FOR STACKED BARPLOTS

# Check unique values of Transect_Number
unique_transects <- unique(waterfall_data_free$Transect_Number)
print(unique_transects)

# use agg_by_tax function to calculate rel abundance !
waterfall_data_free_1 <- agg_by_tax(Transect_Number, waterfall_data_free)

# Define a vector of transect numbers
Transect_Numbers <- c("1", "2", "3", "4", "5") # transect numbers for this transect (1) :^)

# Apply the function to each transect number
waterfall_data_free_list <- lapply(Transect_Numbers, agg_by_tax, main_data_frame = waterfall_data_free) # lapply to all transect numbers

# Combine the results into a single dataframe
waterfall_data_free_combined <- do.call(rbind, waterfall_data_free_list)


#fix out of order transect numbers, do it manually.
# works
waterfall_data_free <- waterfall_data_free_combined %>% 
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "1", 0.0000)) %>%
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "2", 62.80497)) %>%
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "3", 221.81425)) %>%
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "4", 298.40045)) %>%
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "5", 343.60091))


###################### 
#  stacked barplots  #
######################
#      PARTICLE      # 
######################

waterfall_data_part <- waterfall_ps_part %>% subset_samples(Transect_Name == "transect1") %>%
  subset_samples(., Station != "STN22") %>% # Filter out low abundance taxa
  tax_glom(taxrank = "TAX") %>% # agglomerate at genus level
transform_sample_counts(function(x) {x/sum(x)} )# Transform to rel. abundance

waterfall_data_part <- waterfall_data_part %>%
  psmelt() #%>%
#arrange(Transect_Number)    

# Check unique values of Transect_Number
unique_transects <- unique(waterfall_data_part$Transect_Number)
print(unique_transects)

# use agg_by_tax function to calculate rel abundance !
waterfall_data_part_1 <- agg_by_tax(Transect_Number, waterfall_data_part)

# Define a vector of transect numbers
Transect_Numbers <- c("1", "2", "3", "4", "5") # transect numbers for this transect (1) :^)

# Apply the function to each transect number
waterfall_data_part_list <- lapply(Transect_Numbers, agg_by_tax, main_data_frame = waterfall_data_part) # lapply to all transect numbers

# Combine the results into a single dataframe
waterfall_data_part_combined <- do.call(rbind, waterfall_data_part_list)

waterfall_data_part <- waterfall_data_part_combined %>% 
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "1", 0.0000)) %>%
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "2", 62.80497)) %>%
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "3", 221.81425)) %>%
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "4", 298.40045)) %>%
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "5", 343.60091))



## for plots
myColors <- c(brewer.pal(9, "Paired"), "#A43D27", "#497687", "#5E4987", "darkblue", "lightblue2", "darkgoldenrod", "dodgerblue", "seagreen")
waterfall_data_free$TAX <- as.factor(waterfall_data_free$TAX)
waterfall_data_part$TAX <- as.factor(waterfall_data_part$TAX)
names(myColors) <- levels(c(waterfall_data_free$TAX, waterfall_data_part$TAX))
custom_colors <- scale_colour_manual(name = "Order", values = myColors)

waterfall_plot_labels <- c("STN198","STN002", "STN004", "STN012", "STN014")
waterfall_plot_breaks <- unique(sample_data(waterfall_data_free)$Transect_Number)
waterfall_sec_labels <- seq(0 , 350.60, by=50)
waterfall_sec_breaks <- seq(0 , 350.60, by=50)
waterfall_data_free$Transect_Number <- as.numeric(waterfall_data_free$Transect_Number)
waterfall_data_part$Transect_Number <- as.numeric(waterfall_data_part$Transect_Number)

# free-living plot
waterfall_barplot_free <- ggplot(waterfall_data_free, aes(x = Transect_Number, y = Relative, fill = TAX, group = TAX)) +
  geom_bar(stat = "identity", position="fill", color = "black", linewidth = 0.3) + theme_classic() + ggtitle("\n Free-living (<0.2 µm)") +
  #geom_col(position = "dodge") + # changes to multiple bars
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
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
ggsave("graphics/filtered_waterfall_order_rel_abundance_free.pdf", width = 6.5, height = 4, dpi = 150)

# particle-associate plot 
waterfall_barplot_part <- ggplot(waterfall_data_part, aes(x = Transect_Number, y = Relative, fill = TAX, group = TAX)) +
  geom_bar(stat = "identity", position="fill", color = "black", linewidth = 0.3) + theme_classic() + ggtitle("\n Particle-associated (>2 µm)") +
  # geom_col(position = "dodge") + # changes to multiple bars
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_fill_manual(values = myColors, drop=FALSE) + # set manual colors
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
ggsave("graphics/filtered_waterfall_order_rel_abundance_part.pdf", width = 6.5, height = 4, dpi = 150)

# combined plot

total <- rbind(waterfall_data_part, waterfall_data_free)
# make combined FAKE plot to grab legend from and to put in the comine plot :^)
legend_plot <- ggplot(total, aes(x = Transect_Number, y = Abundance, fill = TAX)) +
  geom_bar(stat = "identity", position="fill", width=2, color = "black", linewidth = 0.3) + theme_classic() +
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

ggsave("graphics/filtered_waterfall_order_combined_relative.pdf", width = 11, height = 6, dpi = 150)

############## barplot of total abundance vs iron-related taxa ##########

# Extract OTU tables from phyloseq objects
# Extract abundance data from phyloseq objects
waterfall_ps_free <- waterfall_ps_free %>% subset_samples(Transect_Name == "transect1")
abundance_all <- phyloseq::otu_table(waterfall_ps_free)
abundance_selected <- phyloseq::otu_table(waterfall_data_free_temp)

# Convert to data frames for easier manipulation
df_abundance_all <- as.data.frame(abundance_all)
df_abundance_selected <- as.data.frame(abundance_selected)

# List of specific taxa names you want to exclude
selected_taxa_names <- rownames(df_abundance_selected)
total_removed_other <- df_abundance_all %>%
  filter(!(rownames(df_abundance_all) %in% selected_taxa_names))

# Calculate total abundance of all taxa excluding selected taxa
total_abundance_selected <- rowSums(df_abundance_selected)
total_abundance_other <- rowSums(total_removed_other)

plot_data <- data.frame(
  Taxa = c(rep("Other", length(total_abundance_other)), rep("Selected", length(total_abundance_selected))),
  Abundance = c(total_abundance_other, total_abundance_selected)
)

plot_data <- plot_data %>%
  mutate(Relative = Abundance / sum(Abundance))

# Calculate total abundance of selected taxa and other taxa
ggplot(plot_data, aes(x = Taxa, y = Relative, fill = Taxa)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Taxa", y = "Total Abundance", title = "Abundance Comparison") +
  scale_fill_manual(values = c("Other" = "lightblue", "Selected" = "salmon")) +
  theme_minimal()


########

# Extract abundance data from phyloseq object
all_ps_free <- ps_noncontam_prev05 %>% subset_samples(Transect_Name == "transect1") %>% subset_samples(Filter_pores == "0.2")
###########

# Extract abundance data from phyloseq objects
abundance_all <- t(as.data.frame(otu_table(all_ps_free)))
abundance_selected <- t(as.data.frame(otu_table(waterfall_data_free_temp)))

# Extract transect numbers from metadata (assuming it's available in both data frames)
transect_numbers_all <- as.factor(sample_data(all_ps_free)$Transect_Number)
transect_numbers_selected <- as.factor(sample_data(waterfall_data_free_temp)$Transect_Number)

calculate_relative_abundance <- function(abundance_df, transect_numbers, selected_taxa_names) {
  # Initialize empty data frame to store debug information
  debug_df <- data.frame(
    Transect = character(),
    Abundance = numeric(),
    data_source = character(),
    stringsAsFactors = FALSE
  )
  
  # Initialize vectors to store results
  total_abundance_other <- numeric(length = length(unique(transect_numbers)))
  total_abundance_selected <- numeric(length = length(unique(transect_numbers)))
  
  # Loop over unique transect numbers
  for (i in 1:length(unique(transect_numbers))) {
    transect <- unique(transect_numbers)[i]
    indices <- which(transect_numbers == transect)
    selected_indices <- colnames(abundance_df) %in% selected_taxa_names
    
    # Calculate total abundance of other taxa
    total_abundance_other[i] <- sum(rowSums(abundance_df[indices, !selected_indices, drop = FALSE]))
    
    # Calculate total abundance of selected taxa
    total_abundance_selected[i] <- sum(rowSums(abundance_df[indices, selected_indices, drop = FALSE]))
    
    # Determine data source based on column name
    if (length(selected_taxa_names) > 0) {
      # Create a row of debug information
      debug_row_other <- data.frame(
        Transect = as.character(transect),
        Abundance = total_abundance_other[i],
        data_source = "other",
        stringsAsFactors = FALSE
      )
      
      debug_row_selected <- data.frame(
        Transect = as.character(transect),
        Abundance = total_abundance_selected[i],
        data_source = "selected",
        stringsAsFactors = FALSE
      )
      
      # Append the debug rows to the debug data frame
      debug_df <- rbind(debug_df, debug_row_other, debug_row_selected)
    }
  }
  
  # Calculate relative abundance (percentage) of selected taxa
  relative_abundance <- total_abundance_selected / (total_abundance_other + total_abundance_selected) * 100
  
  # Create data frame with results
  result_df <- data.frame(
    Transect = as.character(unique(transect_numbers)),
    Relative_Abundance = relative_abundance
  )
  
  # Return both the result data frame and the debug data frame
  list(Result = result_df, Debug = debug_df)
}

fill <- ggplot(relative_abundance_all[["Debug"]], aes(x = Transect, y = Abundance, fill = data_source)) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.3) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Transect Number", y = "Relative Abundance (%)", title = "") +
  scale_fill_manual(values = c("other" = "gray", "selected" = "salmon"),
                    labels = c("All other taxa", "Iron-related taxa")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

dodge <- ggplot(relative_abundance_all[["Debug"]], aes(x = Transect, y = Abundance, fill = data_source)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Transect Number", y = "Total Abundance", title = "") +
  scale_fill_manual(values = c("other" = "gray", "selected" = "salmon"),
                    labels = c("All other taxa", "Iron-related taxa")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

abun <- ggarrange(
  fill, dodge, labels = NULL,
  common.legend = TRUE, legend = "right"
)

ggsave("graphics/filtered_waterfall_total_selected.pdf", width = 8, height = 4, dpi = 150)

