library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggpubr)
library(ggprism)
library(patchwork)
require(grid)
library(LinDA)
#install.packages("IgAScores")
library(IgAScores) # gives relative abundance of count table
library(ggplot2)
##### PICRUSt setup #####
metatable <- read_delim("required_files/artemis-eDNA-metadata-final.tsv", delim="\t") 

metadata <- metatable %>% filter(., Sample.Control == "True.Sample") %>% filter(., sample_name != "STN089.200.fil.dura.r2") %>% filter(., sample_name != "STN078.1040.pre.poly.3.LG")#%>% group_by(Station) #%>% filter(., sample_name != "STN115.35.fil.dura.r1") #%>% distinct(Filter_pores, .keep_all = TRUE) 
# remove sample that isn't in kegg abundance for some reason
#filter(., sample_name %in% wanted_samples)

abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, col_types=c("c", "n"), trim_ws = TRUE) 
ColumnstoKeep <- c("function", metadata$sample_name) # set vector of list of names to keep + KO Name column
KO_abundance_data <- subset(abundance_data, select = ColumnstoKeep) # subset (select columns) based on ColumnstoKeep

#for iron

iron_KO <- read.csv("required_files/KO_Numbers_all_metabolism.csv") # read in KO_number reference
#iron_KO <- iron_KO %>% filter(., Metabolism == "Iron uptake and metabolism") # filter by iron metabolism only

KO_iron_numbers <- iron_KO$KO_Num # set vector for numbers
KO_iron_abundance_data <- KO_abundance_data[KO_abundance_data$`function` %in% KO_iron_numbers, ] #filter by KO_number using KO_iron_number ref

KO_joined <- data.frame(iron_KO$Name, iron_KO$KO_Num)
colnames(KO_joined) <- c("Name", "function")

KO_iron_abundance_type <- left_join(KO_iron_abundance_data, KO_joined) %>% as.data.frame()
KO_iron_abundance_type <- KO_iron_abundance_type[!duplicated(KO_iron_abundance_type$`function`), ]
row.names(KO_iron_abundance_type) <- KO_iron_abundance_type$"function"

# create a phyloseq object for PICRUST data
metadata$Filter_pores <- ifelse(metadata$Filter_pores >= 0.2 & metadata$Filter_pores <= 2.0, "free-living", 
                                ifelse(metadata$Filter_pores == 3.0, "particle-associated", metadata$Filter_pores))
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$sample_name
META <- sample_data(metadata)


#for iron
#KO_iron_abundance_type$Name <- NULL
#KO_iron_abundance_type$`function` <- NULL
#otumaty = as(KO_iron_abundance_type, "matrix")
#rownames(otumaty) <- abundance_data$"#NAME"
#OTUy = otu_table(otumaty, taxa_are_rows=TRUE)

#for normal
KO_abundance_data <- as.data.frame(KO_abundance_data)
KO_abundance_data <- KO_abundance_data %>% distinct(`function`, .keep_all = TRUE)
rownames(KO_abundance_data) <- KO_abundance_data$`function`
KO_abundance_data$`function` <- NULL
otumaty = as(KO_abundance_data, "matrix")
#rownames(otumaty) <- abundance_data$"#NAME"
OTUy = otu_table(otumaty, taxa_are_rows=TRUE)
TAX <- tax_table(iron_KO)
rownames(TAX) <- rownames(iron_KO)

TAX <- subset_taxa(TAX, rownames(OTUy) %in% TAX)

pi_ps <- merge_phyloseq(OTUy, META)

pseq <- pi_ps %>%
  phyloseq_validate()


### set up for plots ###
level_order <- c("STN089", "STN132", "STN106", "STN20", "STN198", "STN002", "STN004", "STN012", "STN115", "STN12.3", "STN014", "STN078", "STN056a", "STN056b", "STN22", "STN068", "STN146", "STN181", "STN174", "STN151.2", "STN153")
depth_order <- c("Surface", "Intermediate", "Bottom_water")

#### plots #####
ps_free <- pseq %>% subset_samples(Filter_pores == "free-living") %>% prune_taxa(taxa_sums(.) > 0, .) 

# particle-associated phyloseq
ps_part <- pseq %>% subset_samples(Filter_pores == "particle-associated") %>% prune_taxa(taxa_sums(.) > 0, .) 


data_free <- ps_free %>% # agglomerate at Order level, can change to different taxonomic level!
  prune_taxa(taxa_sums(.) > 0, .) %>%
  transform_sample_counts(function(x) {x/sum(x)})  # Transform to rel. abundance (normalize data

data_top_free <- data_free %>%
  psmelt() %>% # transform a phyloseq object into a data frame, otherwise graphs wont work
  filter(Abundance > 0.05) %>% # Filter out low abundance taxa
  arrange(unique)

data_top_free <- aggregate(Abundance ~ Station * unique, data = data_top_free, FUN = mean)

# particle-associated
data_part <- ps_part %>%
  prune_taxa(taxa_sums(.) > 0, .) %>%# agglomerate at Order level
  transform_sample_counts(function(x) {x/sum(x)} )  # Transform to rel. abundance (normalize data)

data_top_part <- data_part %>%
  psmelt() %>% # transform a phyloseq object into a data frame, otherwise graphs wont work
  filter(Abundance > 0.05) %>% # Filter out low abundance taxa
  arrange(unique)

# create linegraph
line_free <- ggplot(free_abundance_agg, aes(x = factor(Station, level = level_order), y = relative_abundance, color = factor(Depth_Threshold, level = depth_order), group = factor(Depth_Threshold, level = depth_order))) +
  geom_point(size = 1.5) +  # Add points
  geom_line(linewidth=0.7) +  
  ylim(c(0,0.04)) +# Connect points with lines
  labs(x = "",y = "Relative Abundance", color = "Depth") +
  theme_classic2() +
  scale_color_manual(values = c("Surface" = "seagreen", "Intermediate" = "dodgerblue", "Bottom_water" = "red")) +  # Customize color
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=9),
        legend.position = "none") +
  geom_vline(xintercept = c(4.5,11.5), linetype = "dashed", linewidth=0.6, color = "black") + # Add vertical lines
  ggtitle("Free-living (<0.2 µm)")

#### particle - associated ####
part_abundance <- all_stations %>% filter(Filter_pores >= 2) # filter by particle-associated
part_abundance_agg <- aggregate(column_sums ~ Depth_Threshold * Station, data = part_abundance, FUN = mean) # take mean per TYPE, for each station
total_count <- sum(part_abundance_agg$column_sums)

# Calculate relative abundance
part_abundance_agg$relative_abundance <- part_abundance_agg$column_sums / total_count

# create linegraph
line_part <- ggplot(part_abundance_agg, aes(x = factor(Station, level = level_order), y = relative_abundance, color = factor(Depth_Threshold, level = depth_order), group = factor(Depth_Threshold, level = depth_order))) +
  geom_point(size = 1.5) +  # Add points
  geom_line(linewidth=0.7) +  # Connect points with lines
  ylim(c(0,0.04)) +
  labs(x = "Station", y = "Relative Abundance", color = "Depth") +
  theme_classic2() +
  scale_color_manual(values = c("Surface" = "seagreen", "Intermediate" = "dodgerblue", "Bottom_water" = "red")) +  # Customize color
  theme(axis.text.x = element_text(size=9, angle=90),
        axis.text.y = element_text(size=9),
        legend.position = "none") +
  geom_vline(xintercept = c(4.5,11.5), linetype = "dashed", linewidth=0.6, color = "black") + # Add vertical lines
ggtitle("Particle-associated (>3 µm)")

# combine the two line plots
lines_combined <- ggarrange(
  line_free, line_part, labels = NULL,
  common.legend = TRUE, legend = "right", align="h", ncol=1
)
ggsave("PICRUSt_analysis/graphics/rel_abund_depth_and_pores_in_depth.pdf", width = 8, height = 7, dpi = 150)

### plot includes both.
all_stations_agg <- aggregate(column_sums ~ Depth_Threshold * Station * Filter_pores, data = all_stations, FUN = mean) # take mean per TYPE, for each station
total_count <- sum(all_stations_agg$column_sums)
# Calculate relative abundance
all_stations_agg$relative_abundance <- all_stations_agg$column_sums / total_count
all_stations_agg$Filter_pores <- as.factor(all_stations_agg$Filter_pores)

# plot
ggplot(all_stations_agg, aes(x = factor(Station, level = level_order), y = relative_abundance, color = Filter_pores, shape = Filter_pores, group = Filter_pores)) +
  geom_point(size = 2) +  # Add points
  geom_line() +
  facet_grid(~factor(Depth_Threshold, levels=c("Surface", "Intermediate", "Bottom_water"))~., axes="all", axis.labels = "margins") +# Connect points with lines
  labs(x = "Station", y = "Relative Abundance", color = "Community") +
  theme_classic() +
  scale_x_discrete(breaks = level_order,
                   labels = level_order,
                   drop = FALSE) +
  scale_color_manual(values = c("free-living" = "dodgerblue2", "particle-associated" = "#A43D27")) +  # Customize color
  scale_shape_manual(values = c("free-living" = 18, "particle-associated" = 20)) +  # Customize shape
  guides(color = guide_legend(override.aes = list(shape = c(18, 20)))) +  # Customize legend
  theme(axis.text.x = element_text(size=9, angle=90),
        axis.text.y = element_text(size=9),
        legend.position = "right") + 
  geom_vline(xintercept = c(4.5,11.5), linetype = "dashed", linewidth=0.6, color = "black") + # Add vertical lines
ggtitle("Rel. Abund. of All Iron-Related Pathways for All Stations")
ggsave("PICRUSt_analysis/graphics/relative_abund_station_and_depth.pdf", width = 8, height = 7, dpi = 150)






# random graphs
# graph that does it by station
ggplot(part_abundance_agg, aes(x = Depth_Threshold, y = relative_abundance, color = Depth_Threshold)) +
  geom_point(size = 1) +
  geom_line() +
  facet_wrap(~ Station) +
  labs(x = "Depth Threshold", y = "Relative Abundance", color = "Depth Threshold") +
  theme_minimal()

ggplot(part_abundance_agg, aes(x = factor(Station, level = level_order), y = relative_abundance, fill= factor(Depth_Threshold, level = depth_order))) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Stations", y = "Relative Abundance", fill = "Depth Threshold") +
  theme_minimal()


# BY WATERMASS

watermass_data <- metadata %>% select(sample_name, watertype, Station, Filter_pores) # make temp dataframe for combining
all_stations_mass <- left_join(updated_abundance, watermass_data) # join two dataframes by sample_name
all_stations_mass$Filter_pores <- ifelse(all_stations$Filter_pores >= 0.2 & all_stations$Filter_pores <= 2.0, "free-living", 
                                    ifelse(all_stations$Filter_pores == 3.0, "particle-associated", all_stations$Filter_pores))
all_stations_mass <- aggregate(column_sums ~ watertype * Station * Filter_pores, data = all_stations_mass, FUN = mean) # take mean per TYPE, for each station
total_count <- sum(all_stations_mass$column_sums)

# Calculate relative abundance
all_stations_mass$relative_abundance <- all_stations_mass$column_sums / total_count
ggplot(all_stations_mass, aes(x = factor(Station, level = level_order), y = relative_abundance, fill = Filter_pores)) + 
  facet_grid(~factor(watertype, levels=c("AASW", "AASW-WW", "WW", "WW-CDW", "CDW"))~.) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.3, width=0.9) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "", y = "Relative Abundance (%)", title = "") +
  scale_fill_manual(values = c("free-living" = "dodgerblue2", "particle-associated" = "#A43D27"),
                    labels = c("free-living", "particle-associated"), 
                    name = "Community") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=9, angle=90),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=9 , color="white"),
        panel.spacing.x = unit(0, "points"), # Reducing space between facets#
        panel.border = element_blank()) + # Optionally remove panel borders
  geom_vline(xintercept = c(4.5,11.5), linetype = "dashed", linewidth=0.6, color = "black") +# Add vertical lines
  #theme(plot.title = element_text(hjust = 0.5, size=17)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ggtitle("")# Rotate x-axis labels for better readability


### all stations -- barplot ##
test <- aggregate(. ~ Type, data = KO_iron_abundance_type, FUN = mean)
df_test <- as.data.frame(t(test)) 
colnames(df_test) <- df_test[1,]
df_test <- df_test[-1,]
df_test$sample_name <- rownames(df_test)
test_up <- left_join(df_test, dataframe_to_combine)
test_up_long <- test_up %>% # combines by type
  pivot_longer(cols = -c(Station, Depth_Threshold, Filter_pores, sample_name),
               names_to = "Type",
               values_to = "Abundance")
test_up_long$Abundance <- as.numeric(test_up_long$Abundance)
agg_test <- aggregate(Abundance ~ Depth_Threshold * Station * Filter_pores + Type, data = test_up_long, FUN = mean)

agg_test$Filter_pores <- ifelse(agg_test$Filter_pores >= 0.2 & agg_test$Filter_pores <= 2.0, "free-living", 
                                    ifelse(agg_test$Filter_pores == 3.0, "particle-associated", agg_test$Filter_pores))

pathway_colors <- c("darkred", "darkgoldenrod", "#e66101", "darkorange", "#35978f", "#c7eae5")
type_free$Type <- as.factor(type_free$Type) # setting the Order columns to factor
names(pathway_colors) <- levels(type_free$Type) 
wanted_types <- c("Fe Storage ", "Fe-S", "Fe2+ (ferrous)", "Fe3+ (ferric)", "Siderophore biosynthesis", "Siderophore uptake")
type_free <- agg_test %>% filter(Filter_pores == "free-living") %>% filter(Type %in% wanted_types)

barplot_free <- ggplot(type_free, aes(x = factor(Station, level = level_order), y = Abundance, fill = Type)) + facet_grid(~factor(Depth_Threshold, levels=c("Surface", "Intermediate", "Bottom_water"))~.) + # facet grid seperates by different levels, horizontally
  geom_bar(stat = "identity", position="fill", color= "black", linewidth=0.3, width=0.9) + theme_classic() + # adds black outline to boxes
  scale_y_continuous(expand = c(0, 0)) + # extends the barplots to the axies
 scale_fill_manual(values = pathway_colors, drop = FALSE) + # set the colors with custom colors (myColors)
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=9, angle=90),
        axis.title.y = element_blank(),
        panel.spacing.x = unit(0, "points"), # Reducing space between facets
        strip.background = element_blank(), # Optionally hide the strip background for a cleaner look
        panel.border = element_blank()) + # Optionally remove panel borders
  geom_vline(xintercept = c(4.5,11.5), linetype = "dashed", linewidth=0.6, color = "black") +# Add vertical lines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # for the legend, if you want one
  #ylab("Relative Abundance (Order > 2%) \n") + # remove # if you want y-axis title
  ggtitle("Free-living (<0.2 µm)")

type_part <- agg_test %>% filter(Filter_pores == "particle-associated") %>% filter(Type %in% wanted_types)
barplot_part <- ggplot(type_part, aes(x = factor(Station, level = level_order), y = Abundance, fill = Type)) + facet_grid(~factor(Depth_Threshold, levels=c("Surface", "Intermediate", "Bottom_water"))~.) + # facet grid seperates by different levels, horizontally
  geom_bar(stat = "identity", position="fill", color= "black", linewidth=0.3, width=0.9) + theme_classic() + # adds black outline to boxes
  scale_y_continuous(expand = c(0, 0)) + # extends the barplots to the axies
  scale_fill_manual(values = pathway_colors, drop = FALSE) + # set the colors with custom colors (myColors)
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=9, angle=90),
        axis.title.y = element_blank(),
        panel.spacing.x = unit(0, "points"), # Reducing space between facets
        strip.background = element_blank(), # Optionally hide the strip background for a cleaner look
        panel.border = element_blank()) + # Optionally remove panel borders
  geom_vline(xintercept = c(4.5,11.5), linetype = "dashed", linewidth=0.6, color = "black") +# Add vertical lines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # for the legend, if you want one
  #ylab("Relative Abundance (Order > 2%) \n") + # remove # if you want y-axis title
  ggtitle("Particle-associated (>3 µm)")

ps_combined <- ggarrange(
  barplot_free, barplot_part, labels = NULL,
  common.legend = TRUE, legend = "right"
)

ggsave("PICRUSt_analysis/graphics/rel_abun_by_type.pdf", width = 13, height = 7, dpi = 150)


#comparing abundance to iron vs all other related
KO_all_abundance <- KO_abundance_data 
KO_all_abundance$Type <- ifelse(KO_all_abundance$`#NAME` %in% KO_iron_numbers, "Iron-related", "Other")
KO_all_abundance$"#NAME" <- NULL
save <- KO_all_abundance$Type
KO_all_abundance$Type <- NULL
KO_all_abundance <- relabund(KO_all_abundance)
KO_all_abundance$Type <- save
threshold <- 0  # Adjust this value as needed

# Filter rows below the threshold
filtered_data <- KO_all_abundance %>%
  filter(if_all(where(is.numeric), ~. >= threshold))

all_test <- aggregate(. ~ Type, data = KO_all_abundance, FUN = mean)
all_df_test <- as.data.frame(t(all_test)) 
colnames(all_df_test) <- all_df_test[1,]
all_df_test <- all_df_test[-1,] # run twice
all_df_test$sample_name <- rownames(all_df_test)
all_test_up <- left_join(all_df_test, dataframe_to_combine)
all_test_up_long <- all_test_up %>% # combines by type
  pivot_longer(cols = -c(Station, Depth_Threshold, Filter_pores, sample_name),
               names_to = "Type",
               values_to = "Abundance")
all_test_up_long$Abundance <- as.numeric(all_test_up_long$Abundance)
agg_test <- aggregate(Abundance ~ Depth_Threshold * Station * Filter_pores + Type, data = all_test_up_long, FUN = mean)

agg_test$Filter_pores <- ifelse(agg_test$Filter_pores >= 0.2 & agg_test$Filter_pores <= 2.0, "free-living", 
                                ifelse(agg_test$Filter_pores == 3.0, "particle-associated", agg_test$Filter_pores))

type_part <- agg_test %>% filter(Filter_pores == "particle-associated")
barplot_part <- ggplot(type_part, aes(x = factor(Station, level = level_order), y = Abundance, fill = Type)) + facet_grid(~factor(Depth_Threshold, levels=c("Surface", "Intermediate", "Bottom_water"))~.) + # facet grid seperates by different levels, horizontally
  geom_bar(stat = "identity", position="fill", color= "black", linewidth=0.3, width=0.9) + theme_classic() + # adds black outline to boxes
  scale_y_continuous(expand = c(0, 0)) + # extends the barplots to the axies
  #scale_fill_manual(values = pathway_colors, drop = FALSE) + # set the colors with custom colors (myColors)
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=9, angle=90),
        axis.title.y = element_blank(),
        panel.spacing.x = unit(0, "points"), # Reducing space between facets
        strip.background = element_blank(), # Optionally hide the strip background for a cleaner look
        panel.border = element_blank()) + # Optionally remove panel borders
  geom_vline(xintercept = c(4.5,11.5), linetype = "dashed", linewidth=0.6, color = "black") +# Add vertical lines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # for the legend, if you want one
  #ylab("Relative Abundance (Order > 2%) \n") + # remove # if you want y-axis title
  ggtitle("Particle-associated (>3 µm)")

type_free <- agg_test %>% filter(Filter_pores == "free-living")
barplot_free <- ggplot(type_free, aes(x = factor(Station, level = level_order), y = Abundance, fill = Type)) + facet_grid(~factor(Depth_Threshold, levels=c("Surface", "Intermediate", "Bottom_water"))~.) + # facet grid seperates by different levels, horizontally
  geom_bar(stat = "identity", position="fill", color= "black", linewidth=0.3, width=0.9) + theme_classic() + # adds black outline to boxes
  scale_y_continuous(expand = c(0, 0)) + # extends the barplots to the axies
  #scale_fill_manual(values = pathway_colors, drop = FALSE) + # set the colors with custom colors (myColors)
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=9, angle=90),
        axis.title.y = element_blank(),
        panel.spacing.x = unit(0, "points"), # Reducing space between facets
        strip.background = element_blank(), # Optionally hide the strip background for a cleaner look
        panel.border = element_blank()) + # Optionally remove panel borders
  geom_vline(xintercept = c(4.5,11.5), linetype = "dashed", linewidth=0.6, color = "black") +# Add vertical lines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # for the legend, if you want one
  #ylab("Relative Abundance (Order > 2%) \n") + # remove # if you want y-axis title
  ggtitle("Free-living (<0.2 µm)")

coastal_combined <- ggarrange(
  barplot_free, barplot_part, labels = NULL,
  common.legend = TRUE, legend = "right"
)

ggsave("PICRUSt_analysis/graphics/all_vs_iron-related_barplot.pdf", width = 13, height = 7, dpi = 150)
