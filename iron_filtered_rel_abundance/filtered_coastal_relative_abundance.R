# free-living phyloseq
coastal_ps_free <- filtered_ps %>% subset_samples(Filter_pores == "0.2") %>% prune_taxa(taxa_sums(.) > 0, .)

# particle-associated phyloseq
coastal_ps_part <- filtered_ps %>% subset_samples(Filter_pores >= "2") %>% prune_taxa(taxa_sums(.) > 0, .) 

###################### 
#  stacked barplots  #
######################
#    FREE-LIVING     # 
######################

# Create a data frame for freeliving
coastal_data_free <- coastal_ps_free %>% subset_samples(Coastal_Current_Name == "transect3") %>%
  tax_glom(taxrank = "TAX") # agglomerate at tax level
#transform_sample_counts(function(x) {x/sum(x)} ) # Transform to rel. abundance
coastal_data_free_temp <- coastal_data_free  # from above  


coastal_data_free <- coastal_data_free %>%
  psmelt() %>% filter(., Coastal_Current_Number != "4") %>%  
  filter(Abundance != 0) %>%
  arrange(Coastal_Current_Number)
#group_by(Coastal_Current_Number) %>% 
#mutate(Relative = Abundance/sum(Abundance))  #### USE THIS WHEN NOT USING THE METHOD BELOW. 
### FOR STACKED BARPLOTS

# Check unique values of Coastal_Current_Number
unique_transects <- unique(coastal_data_free$Coastal_Current_Number)
print(unique_transects)

# use agg_by_tax function to calculate rel abundance !
#coastal_data_free_1 <- agg_by_tax_coastal(transect_number, coastal_data_free) #just for one transect number

# Define a vector of transect numbers
transect_numbers <- c("1", "2", "3", "5", "7", "8", "9") # transect numbers for this transect (1) :^)

# Apply the function to each transect number
coastal_data_free_list <- lapply(transect_numbers, agg_by_tax_coastal, main_data_frame = coastal_data_free) # lapply to all transect numbers

# Combine the results into a single dataframe
coastal_data_free_combined <- do.call(rbind, coastal_data_free_list)


#fix out of order transect numbers, do it manually.
# works
coastal_data_free <- coastal_data_free_combined %>% 
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "1", 0.0000)) %>%
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "2", 32.67378)) %>%
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "3", 48.56614)) %>%
  #mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "4", 61.17885)) %>% # REMOVED
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "5", 62.78023)) %>%
  #mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "6", 80.58525)) %>%  # REMOVED BECAUSE doesn't have free-living.
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "7", 101.50049)) %>%
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "8", 106.99770)) %>%
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "9", 133.80764))


###################### 
#  stacked barplots  #
######################
#      PARTICLE      # 
######################

coastal_data_part <- coastal_ps_part %>% subset_samples(Coastal_Current_Name == "transect3") %>% # Filter out low abundance taxa
  tax_glom(taxrank = "TAX") #%>% # agglomerate at genus level
#transform_sample_counts(function(x) {x/sum(x)} )# Transform to rel. abundance
coastal_data_part_temp <- coastal_data_part

coastal_data_part <- coastal_data_part %>%
  psmelt() %>% filter(., Coastal_Current_Number != "4") %>%  
  filter(Abundance != 0) %>% # filter out iron taxa = 0, before calculating relative abundance
  arrange(Coastal_Current_Number)
#arrange(Coastal_Current_Number)    

# Check unique values of Coastal_Current_Number
unique_transects <- unique(coastal_data_part$Coastal_Current_Number)
print(unique_transects)

# use agg_by_tax function to calculate rel abundance !
#coastal_data_part_1 <- agg_by_tax(transect_number, coastal_data_part)

# Define a vector of transect numbers
transect_numbers <- c("1", "2", "3", "5", "6", "7", "8", "9") # transect numbers for this transect (1) :^)

# Apply the function to each transect number
coastal_data_part_list <- lapply(transect_numbers, agg_by_tax_coastal, main_data_frame = coastal_data_part) # lapply to all transect numbers

# Combine the results into a single dataframe
coastal_data_part_combined <- do.call(rbind, coastal_data_part_list)

coastal_data_part <- coastal_data_part_combined %>% 
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "1", 0.0000)) %>%
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "2", 32.67378)) %>%
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "3", 48.56614)) %>%
  #mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "4", 61.17885)) %>%
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "5", 62.78023)) %>%
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "6", 80.58525)) %>%
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "7", 101.50049)) %>%
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "8", 106.99770)) %>%
  mutate(Coastal_Current_Number = replace(Coastal_Current_Number, Coastal_Current_Number == "9", 133.80764))

## for plots
myColors <- c(brewer.pal(9, "Paired"), "#A43D27") #"#497687", "#5E4987", "darkblue", "lightblue2", "darkgoldenrod", "dodgerblue", "seagreen")
coastal_data_free$TAX <- as.factor(coastal_data_free$TAX)
coastal_data_part$TAX <- as.factor(coastal_data_part$TAX)
names(myColors) <- levels(c(coastal_data_free$TAX, coastal_data_part$TAX))
custom_colors <- scale_colour_manual(name = "Order", values = myColors)

coastal_plot_labels <- c("89", "132", "106", "14", "78", "56b", "68", "146") # DIFFERENT BC OF FREE-LIVING
coastal_plot_breaks <- unique(coastal_data_part$Coastal_Current_Number) # HAVE TO CHANGE
coastal_sec_labels <- seq(0 , 140, by=20)
coastal_sec_breaks <- seq(0 , 140, by=20)
coastal_data_free$Coastal_Current_Number <- as.numeric(coastal_data_free$Coastal_Current_Number)
coastal_data_part$Coastal_Current_Number <- as.numeric(coastal_data_part$Coastal_Current_Number)

# free-living plot
coastal_barplot_free <- ggplot(coastal_data_free, aes(x = Coastal_Current_Number, y = Relative, fill = TAX)) +
  geom_bar(stat = "identity", position="fill", width=5) + theme_classic() + ggtitle("\n Free-living (<0.2 µm)") +
  #geom_col(position = "dodge") + # changes to multiple bars
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  #geom_vline(xintercept = 105, color = "black", linetype = "dotted") +
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
ggsave("graphics/filtered_coastal_order_rel_abundance_free.pdf", width = 6.5, height = 4, dpi = 150)

# particle-associate plot 
coastal_barplot_part <- ggplot(coastal_data_part, aes(x = Coastal_Current_Number, y = Relative, fill = TAX)) +
  geom_bar(stat = "identity", position="fill", width=5) + theme_classic() + ggtitle("\n Particle-associated (>2 µm)") +
  # geom_col(position = "dodge") + # changes to multiple bars
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  #geom_vline(xintercept = 105, color = "black", linetype = "dotted") +
  scale_fill_manual(values = myColors, drop=FALSE) + # set manual colors
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
ggsave("graphics/filtered_coastal_order_rel_abundance_part.pdf", width = 6.5, height = 4, dpi = 150)

# combined plot

total <- rbind(coastal_data_part, coastal_data_free)
# make combined FAKE plot to grab legend from and to put in the comine plot :^)
legend_plot <- ggplot(total, aes(x = Coastal_Current_Number, y = Abundance, fill = TAX)) +
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

ggsave("graphics/filtered_coastal_order_combined_relative.pdf", width = 13, height = 7, dpi = 150)


#### total vs iron-only ####

# Extract abundance data from phyloseq object
all_ps_free <- ps_noncontam_prev05 %>% subset_samples(Coastal_Current_Name == "transect3") %>% subset_samples(Filter_pores == "0.2")
all_ps_part <- ps_noncontam_prev05 %>% subset_samples(Coastal_Current_Name == "transect3") %>% subset_samples(Filter_pores >= "2")
#coastal_data_free_temp <- coastal_data_free  # from above  
#coastal_data_part_temp <- coastal_data_part

###############
# free-living #
###############
abundance_selected_free <- phyloseq::otu_table(coastal_data_free_temp)
df_abundance_selected_free <- as.data.frame(abundance_selected_free)
selected_taxa_names_free <- rownames(df_abundance_selected_free) # List of specific taxa names, that i want to exclude

# Extract abundance data from phyloseq objects
abundance_all_free <- t(as.data.frame(otu_table(all_ps_free)))
# Extract transect numbers from metadata (assuming it's available in both data frames)
transect_numbers_all_free <- as.factor(sample_data(all_ps_free)$Coastal_Current_Number)
result_coastal_free <- calculate_relative_abundance(abundance_all_free, transect_numbers_all_free, selected_taxa_names_free)

# Access the debug data frame with relative abundance values
coastal_relative_result_free <- result_coastal_free$Debug

fill_free <- ggplot(coastal_relative_result_free, aes(x = Transect, y = Abundance, fill = data_source)) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.3) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Transect Number", y = "Relative Abundance (%)", title = "Free-living") +
  scale_fill_manual(values = c("other" = "gray", "selected" = "salmon"),
                    labels = c("All other taxa", "Iron-related taxa"), 
                    name = "Taxa") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

dodge_free <- ggplot(coastal_relative_result_free, aes(x = Transect, y = Abundance, fill = data_source)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Transect Number", y = "Total Abundance", title = "") +
  scale_fill_manual(values = c("other" = "gray", "selected" = "salmon"),
                    labels = c("All other taxa", "Iron-related taxa")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

abun <- ggarrange(
  fill_free, dodge_free, labels = NULL,
  common.legend = TRUE, legend = "right"
)

#######################
# particle-associated #
#######################
abundance_selected_part <- phyloseq::otu_table(coastal_data_part_temp)
df_abundance_selected_part <- as.data.frame(abundance_selected_part)
selected_taxa_names_part <- rownames(df_abundance_selected_part) # List of specific taxa names, that i want to exclude

# Extract abundance data from phyloseq objects
abundance_all_part <- t(as.data.frame(otu_table(all_ps_part)))
# Extract transect numbers from metadata (assuming it's available in both data frames)
transect_numbers_all_part <- as.factor(sample_data(all_ps_part)$Coastal_Current_Number)
result_coastal_part <- calculate_relative_abundance(abundance_all_part, transect_numbers_all_part, selected_taxa_names_part)

# Access the debug data frame with relative abundance values
coastal_relative_result_part <- result_coastal_part$Debug

fill_part <- ggplot(coastal_relative_result_part, aes(x = Transect, y = Abundance, fill = data_source)) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.3) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Transect Number", y = "Relative Abundance (%)", title = "Particle-associated") +
  scale_fill_manual(values = c("other" = "gray", "selected" = "salmon"),
                    labels = c("All other taxa", "Iron-related taxa"), 
                    name = "Taxa") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

dodge_part <- ggplot(coastal_relative_result_part, aes(x = Transect, y = Abundance, fill = data_source)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Transect Number", y = "Total Abundance", title = "") +
  scale_fill_manual(values = c("other" = "gray", "selected" = "salmon"),
                    labels = c("All other taxa", "Iron-related taxa"), 
                    name = "Taxa") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

coastal_iron_abundance_together <- ggarrange(
  fill_free, fill_part, labels = NULL,
  common.legend = TRUE, legend = "right"
)

ggsave("iron_filtered_rel_abundance/graphics/coastal_total_vs_iron.pdf", width = 8, height = 4, dpi = 150)



