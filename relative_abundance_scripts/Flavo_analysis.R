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

ps_sub <- ps_noncontam_prev05 %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "Mitochondria" &
      Class   != "Chloroplast" &
      Order   != "Chloroplast" &
      Family  != "Mitochondria" &
      Genus   != "uncultured"
  )

ps_sub <- subset_samples(ps_sub, Sample.Control == "True.Sample")

order_2_keep=c("Flavobacteriales")

pseudo_tax <- taxonomy %>%
  filter(Order %in% order_2_keep) %>% 
  rownames_to_column("id")

list_filtered_tax <- pseudo_tax$id # create a list of tax id's for taxonomy
ASV_temp <- as.data.frame(ASV) %>% rownames_to_column("id") # create dataframe from ASV, then create rowname for id
new_ASV <- filter(ASV_temp, id %in% list_filtered_tax) # filter ASV based on list of filtered taxonomy
rownames(new_ASV) <- new_ASV$id # make taxa rownames
new_ASV$id <- NULL # remove column

rownames(pseudo_tax) <- pseudo_tax$id # make taxa rownames
pseudo_tax$id <- NULL # remove column

#### NEW PHYLOSEQ OBJECT
filtered_ASV <- otu_table(new_ASV, taxa_are_rows = TRUE)
filtered_taxmat <- as.matrix(pseudo_tax)
filtered_TAX <- tax_table(filtered_taxmat) #covert taxonomy table to phyloseq object

filtered_ps <- merge_phyloseq(filtered_ASV, filtered_TAX, META)

# free-living phyloseq
ps_free <- filtered_ps %>% subset_samples(Filter_pores == "0.2") %>% prune_taxa(taxa_sums(.) > 0, .) 

data_top_free <- ps_free %>%
  psmelt() %>% # transform a phyloseq object into a data frame, otherwise graphs wont work
  filter(Abundance > 0.02) %>% 
  filter(Sample.Control == "True.Sample") # Filter out low abundance tax

data_top_free <- aggregate(Abundance ~ Station * Genus * Depth_Threshold, data = data_top_free, FUN = mean)

ps_part <- filtered_ps %>% subset_samples(Filter_pores >= "2") %>% prune_taxa(taxa_sums(.) > 0, .) 

data_top_part <- ps_part %>%
  psmelt() %>% # transform a phyloseq object into a data frame, otherwise graphs wont work
  filter(Abundance > 0.02) %>% 
  filter(Sample.Control == "True.Sample")# Filter out low abundance tax

data_top_part <- aggregate(Abundance ~ Station * Genus * Depth_Threshold, data = data_top_part, FUN = mean)

myColors <- c(brewer.pal(9, "Paired"),'#e66101','darkgreen','#fdb863','#5e3c99', '#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','darkred','#c7eae5','#80cdc1','#35978f','#01665e','#003c30', "#A43D27", "darkgreen", "#497687", "#5E4987", "darkgoldenrod", "lightblue2", "darkblue", "dodgerblue", "seagreen", "purple", "black", '#f6e8c3','darkred','#c7eae5','#80cdc1', "white") # this must equal the levels of the Order
data_top_free$Genus <- as.factor(data_top_free$Genus) # setting the Order columns to factor
data_top_part$Genus <- as.factor(data_top_part$Genus) # setting the Order columns to factor
names(myColors) <- levels(c(data_top_free$Genus, data_top_part$Genus)) # setting the names of the colors to coordinate with the Order columns of each dataframe

barplot_free <- ggplot(data_top_free, aes(x = factor(Station, level = level_order), y = Abundance, fill = Genus, group = Genus)) + facet_grid(~factor(Depth_Threshold, levels=c("Surface", "Intermediate", "Bottom_water"))~., scales = "free_x", space = "free_x") +
  geom_bar(stat = "identity", position="fill", color= "black", linewidth=0.3, width=0.9) + theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = myColors, drop = FALSE) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=9, angle=90),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=9 , color="white"),
        legend.position = "right",
        panel.spacing.x = unit(0, "points"), # Reducing space between facets
        strip.background = element_blank(), # Optionally hide the strip background for a cleaner look
        panel.border = element_blank()) + # Optionally remove panel borders
  geom_vline(xintercept = c(4.5,11.5), linetype = "dashed", linewidth=0.6, color = "black") +# Add vertical lines
  #theme(plot.title = element_text(hjust = 0.5, size=17)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ggtitle("Free-living (<0.2 µm)")

# particle-associated phyloseq
barplot_part <- ggplot(data_top_part, aes(x = factor(Station, level = level_order), y = Abundance, fill = Genus, group = Genus)) + facet_grid(~factor(Depth_Threshold, levels=c("Surface", "Intermediate", "Bottom_water"))~., scales = "free_x", space = "free_x") +
  geom_bar(stat = "identity", position="fill", color= "black", linewidth=0.3, width=0.9) + theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = myColors, drop = FALSE) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=9, angle=90),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=9 , color="white"),
        legend.position = "right",
        panel.spacing.x = unit(0, "points"), # Reducing space between facets
        strip.background = element_blank(), # Optionally hide the strip background for a cleaner look
        panel.border = element_blank()) + # Optionally remove panel borders
  geom_vline(xintercept = c(4.5,11.5), linetype = "dashed", linewidth=0.6, color = "black") +# Add vertical lines
  #theme(plot.title = element_text(hjust = 0.5, size=17)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ggtitle("Particle-associated (>3 µm)")

ps_combined <- ggarrange(
  barplot_free, barplot_part, labels = NULL,
  common.legend = TRUE, legend = "right"
)
# to remove the white space between the two community plots, you'd have to play with the plot.margins of each plot individually!
# something like: plot.margin=unit(c(1,1,-0.5,1), "cm")), where the margins follow the following structure:
# unit(c(top, right, bottom, left), units).

annotate_figure(ps_combined, top = text_grob("Total Relative Abundance for All ARTEMIS Stations", 
                                             color = "black", face = "bold", size = 18))

ggsave("relative_abundance_scripts/graphics/flavo_analysis.pdf", width = 13, height = 7, dpi = 150)

