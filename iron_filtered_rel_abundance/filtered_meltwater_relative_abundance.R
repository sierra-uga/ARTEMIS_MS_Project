# free-living phyloseq
meltwater_ps_free <- filtered_ps %>% subset_samples(Filter_pores == "0.2") %>% prune_taxa(taxa_sums(.) > 0, .)

# particle-associated phyloseq
meltwater_ps_part <- filtered_ps %>% subset_samples(Filter_pores >= "2") %>% prune_taxa(taxa_sums(.) > 0, .) 

###################### 
#  stacked barplots  #
######################
#    FREE-LIVING     # 
######################

# Create a data frame for freeliving
meltwater_data_free <- meltwater_ps_free %>% subset_samples(Transect_Name == "transect2") %>%
  subset_samples(., Station != "STN22") %>% # Filter out low abundance taxa
  tax_glom(taxrank = "TAX") # agglomerate at tax level
#transform_sample_counts(function(x) {x/sum(x)} ) # Transform to rel. abundance

meltwater_data_free_temp <- meltwater_data_free

meltwater_data_free <- meltwater_data_free %>%
  psmelt() #%>%  # Melt to long format
#arrange(Transect_Number) %>%
#group_by(Transect_Number) %>% 
#mutate(Relative = Abundance/sum(Abundance))  #### USE THIS WHEN NOT USING THE METHOD BELOW. 
                                              ### FOR STACKED BARPLOTS

### trying to separate each station to calculate relative abundance
meltwater_data_free_56 <- subset(meltwater_data_free, (Transect_Number == "1")) 
meltwater_data_free_56_temp <- as.data.frame(meltwater_data_free_56$Abundance, meltwater_data_free_56$TAX)
meltwater_data_free_56_temp$TAX <- rownames(meltwater_data_free_56_temp)
meltwater_data_free_56_temp <- aggregate(. ~ TAX, data = meltwater_data_free_56_temp, FUN=sum) 
colnames(meltwater_data_free_56_temp) <- c("TAX", "Abundance")
meltwater_data_free_56_temp$Transect_Number <- "1"
meltwater_data_free_56 <- meltwater_data_free_56_temp %>%
  mutate(Relative = Abundance/sum(Abundance))

meltwater_data_free_68 <- subset(meltwater_data_free, (Transect_Number == "2"))
meltwater_data_free_68_temp <- as.data.frame(meltwater_data_free_68$Abundance, meltwater_data_free_68$TAX)
meltwater_data_free_68_temp$TAX <- rownames(meltwater_data_free_68_temp)
meltwater_data_free_68_temp <- aggregate(. ~ TAX, data = meltwater_data_free_68_temp, FUN=sum) 
colnames(meltwater_data_free_68_temp) <- c("TAX", "Abundance")
meltwater_data_free_68_temp$Transect_Number <- "2"
meltwater_data_free_68 <- meltwater_data_free_68_temp %>%
  mutate(Relative = Abundance/sum(Abundance))

meltwater_data_free_146 <- subset(meltwater_data_free, (Transect_Number == "3"))
meltwater_data_free_146_temp <- as.data.frame(meltwater_data_free_146$Abundance, meltwater_data_free_146$TAX)
meltwater_data_free_146_temp$TAX <- rownames(meltwater_data_free_146_temp)
meltwater_data_free_146_temp <- aggregate(. ~ TAX, data = meltwater_data_free_146_temp, FUN=sum) 
colnames(meltwater_data_free_146_temp) <- c("TAX", "Abundance")
meltwater_data_free_146_temp$Transect_Number <- "3"
meltwater_data_free_146 <- meltwater_data_free_146_temp %>%
  mutate(Relative = Abundance/sum(Abundance))

meltwater_data_free_merged <- rbind(meltwater_data_free_56, meltwater_data_free_68)
meltwater_data_free <- rbind(meltwater_data_free_merged, meltwater_data_free_146)

#fix out of order transect numbers, do it manually.
# works
meltwater_data_free <- meltwater_data_free %>% 
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "1", 0.00000000)) %>%
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "2", 16.33216538)) %>%
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "3", 40.58487518))


###################### 
#  stacked barplots  #
######################
#      PARTICLE      # 
######################

meltwater_data_part <- meltwater_ps_part %>% subset_samples(Transect_Name == "transect2") %>%
  subset_samples(., Station != "STN22") %>% # Filter out low abundance taxa
  tax_glom(taxrank = "TAX") #%>% # agglomerate at genus level
#transform_sample_counts(function(x) {x/sum(x)} )# Transform to rel. abundance
meltwater_data_part_temp <- meltwater_data_part

meltwater_data_part <- meltwater_data_part %>%
  psmelt() #%>%
#arrange(Transect_Number)     

### trying to separate each station to calculate relative abundance
meltwater_data_part_56 <- subset(meltwater_data_part, (Transect_Number == "1")) 
meltwater_data_part_56_temp <- as.data.frame(meltwater_data_part_56$Abundance, meltwater_data_part_56$TAX)
meltwater_data_part_56_temp$TAX <- rownames(meltwater_data_part_56_temp)
meltwater_data_part_56_temp <- aggregate(. ~ TAX, data = meltwater_data_part_56_temp, FUN=sum) 
colnames(meltwater_data_part_56_temp) <- c("TAX", "Abundance")
meltwater_data_part_56_temp$Transect_Number <- "1"
meltwater_data_part_56 <- meltwater_data_part_56_temp %>%
  mutate(Relative = Abundance/sum(Abundance))

meltwater_data_part_68 <- subset(meltwater_data_part, (Transect_Number == "2"))
meltwater_data_part_68_temp <- as.data.frame(meltwater_data_part_68$Abundance, meltwater_data_part_68$TAX)
meltwater_data_part_68_temp$TAX <- rownames(meltwater_data_part_68_temp)
meltwater_data_part_68_temp <- aggregate(. ~ TAX, data = meltwater_data_part_68_temp, FUN=sum) 
colnames(meltwater_data_part_68_temp) <- c("TAX", "Abundance")
meltwater_data_part_68_temp$Transect_Number <- "2"
meltwater_data_part_68 <- meltwater_data_part_68_temp %>%
  mutate(Relative = Abundance/sum(Abundance))

meltwater_data_part_146 <- subset(meltwater_data_part, (Transect_Number == "3"))
meltwater_data_part_146_temp <- as.data.frame(meltwater_data_part_146$Abundance, meltwater_data_part_146$TAX)
meltwater_data_part_146_temp$TAX <- rownames(meltwater_data_part_146_temp)
meltwater_data_part_146_temp <- aggregate(. ~ TAX, data = meltwater_data_part_146_temp, FUN=sum) 
colnames(meltwater_data_part_146_temp) <- c("TAX", "Abundance")
meltwater_data_part_146_temp$Transect_Number <- "3"
meltwater_data_part_146 <- meltwater_data_part_146_temp %>%
  mutate(Relative = Abundance/sum(Abundance))

meltwater_data_part_merged <- rbind(meltwater_data_part_56, meltwater_data_part_68)
meltwater_data_part <- rbind(meltwater_data_part_merged, meltwater_data_part_146)

meltwater_data_part <- meltwater_data_part %>% 
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "1", 0.00000000)) %>%
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "2", 16.33216538)) %>%
  mutate(Transect_Number = replace(Transect_Number, Transect_Number == "3", 40.58487518))
meltwater_data_part$Transect_Number <- as.numeric(meltwater_data_part$Transect_Number)


## for plots
myColors <- c(brewer.pal(9, "Paired"), "#A43D27", "#497687", "#5E4987", "darkblue", "lightblue2", "darkgoldenrod", "dodgerblue", "seagreen")
meltwater_data_free$TAX <- as.factor(meltwater_data_free$TAX)
meltwater_data_part$TAX <- as.factor(meltwater_data_part$TAX)
names(myColors) <- levels(c(meltwater_data_free$TAX, meltwater_data_part$TAX))
custom_colors <- scale_colour_manual(name = "Order", values = myColors)

meltwater_plot_labels <- c("56", "68", "146")
meltwater_plot_breaks <- unique(sample_data(meltwater_data_free)$Transect_Number)
meltwater_sec_labels <- seq(0 , 42, by=10)
meltwater_sec_breaks <- seq(0 , 42, by=10)
meltwater_data_free$Transect_Number <- as.numeric(meltwater_data_free$Transect_Number)

# free-living plot
meltwater_barplot_free <- ggplot(meltwater_data_free, aes(x = Transect_Number, y = Relative, fill = TAX)) +
  geom_bar(stat = "identity", position="dodge", width=7) + theme_classic() + ggtitle("\n Free-living (<0.2 µm)") +
  #geom_col(position = "dodge") + # changes to multiple bars
  scale_y_continuous(limits = c(0, .6), expand = c(0, 0)) +
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
ggsave("graphics/filtered_meltwater_order_rel_abundance_free.pdf", width = 6.5, height = 4, dpi = 150)

# particle-associate plot 
meltwater_barplot_part <- ggplot(meltwater_data_part, aes(x = Transect_Number, y = Relative, fill = TAX)) +
  geom_bar(stat = "identity", position="dodge", width=7) + theme_classic() + ggtitle("\n Particle-associated (>2 µm)") +
  # geom_col(position = "dodge") + # changes to multiple bars
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_fill_manual(values = myColors, drop=FALSE) + # set manual colors
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
ggsave("graphics/filtered_meltwater_order_rel_abundance_part.pdf", width = 6.5, height = 4, dpi = 150)

# combined plot

total <- rbind(meltwater_data_part, meltwater_data_free)
# make combined FAKE plot to grab legend from and to put in the comine plot :^)
legend_plot <- ggplot(total, aes(x = Transect_Number, y = Abundance, fill = TAX)) +
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

ggsave("graphics/filtered_meltwater_order_combined_relative.pdf", width = 13, height = 7, dpi = 150)


#### iron vs total ###

# Extract abundance data from phyloseq object
all_ps_free <- ps_noncontam_prev05 %>% subset_samples(Transect_Name == "transect2") %>% subset_samples(Filter_pores == "0.2")
all_ps_part <- ps_noncontam_prev05 %>% subset_samples(Transect_Name == "transect2") %>% subset_samples(Filter_pores >= "2")
#meltwater_data_free_temp <- meltwater_data_free  # from above  
#meltwater_data_part_temp <- meltwater_data_part

###############
# free-living #
###############
abundance_selected_free <- phyloseq::otu_table(meltwater_data_free_temp)
df_abundance_selected_free <- as.data.frame(abundance_selected_free)
selected_taxa_names_free <- rownames(df_abundance_selected_free) # List of specific taxa names, that i want to exclude

# Extract abundance data from phyloseq objects
abundance_all_free <- t(as.data.frame(otu_table(all_ps_free)))
# Extract transect numbers from metadata (assuming it's available in both data frames)
transect_numbers_all_free <- as.factor(sample_data(all_ps_free)$Transect_Number)
result_meltwater_free <- calculate_relative_abundance(abundance_all_free, transect_numbers_all_free, selected_taxa_names_free)

# Access the debug data frame with relative abundance values
meltwater_relative_result_free <- result_meltwater_free$Debug

fill_free <- ggplot(meltwater_relative_result_free, aes(x = Transect, y = Abundance, fill = data_source)) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.3) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Transect Number", y = "Relative Abundance (%)", title = "Free-living") +
  scale_fill_manual(values = c("other" = "gray", "selected" = "salmon"),
                    labels = c("All other taxa", "Iron-related taxa"), 
                    name = "Taxa") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

dodge_free <- ggplot(meltwater_relative_result_free, aes(x = Transect, y = Abundance, fill = data_source)) +
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

ggsave("graphics/filtered_meltwater_total_selected.pdf", width = 8, height = 4, dpi = 150)

#######################
# particle-associated #
#######################
abundance_selected_part <- phyloseq::otu_table(meltwater_data_part_temp)
df_abundance_selected_part <- as.data.frame(abundance_selected_part)
selected_taxa_names_part <- rownames(df_abundance_selected_part) # List of specific taxa names, that i want to exclude

# Extract abundance data from phyloseq objects
abundance_all_part <- t(as.data.frame(otu_table(all_ps_part)))
# Extract transect numbers from metadata (assuming it's available in both data frames)
transect_numbers_all_part <- as.factor(sample_data(all_ps_part)$Transect_Number)
result_meltwater_part <- calculate_relative_abundance(abundance_all_part, transect_numbers_all_part, selected_taxa_names_part)

# Access the debug data frame with relative abundance values
meltwater_relative_result_part <- result_meltwater_part$Debug

fill_part <- ggplot(meltwater_relative_result_part, aes(x = Transect, y = Abundance, fill = data_source)) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.3) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Transect Number", y = "Relative Abundance (%)", title = "Particle-associated") +
  scale_fill_manual(values = c("other" = "gray", "selected" = "salmon"),
                    labels = c("All other taxa", "Iron-related taxa"), 
                    name = "Taxa") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

dodge_part <- ggplot(meltwater_relative_result_part, aes(x = Transect, y = Abundance, fill = data_source)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Transect Number", y = "Total Abundance", title = "") +
  scale_fill_manual(values = c("other" = "gray", "selected" = "salmon"),
                    labels = c("All other taxa", "Iron-related taxa"), 
                    name = "Taxa") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability

meltwater_iron_abundance_together <- ggarrange(
  fill_free, fill_part, labels = NULL,
  common.legend = TRUE, legend = "right"
)

ggsave("iron_filtered_rel_abundance/graphics/meltwater_total_vs_iron.pdf", width = 8, height = 4, dpi = 150)



