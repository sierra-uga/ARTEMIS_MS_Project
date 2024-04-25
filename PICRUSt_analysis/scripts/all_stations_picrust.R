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

##### PICRUSt setup #####
abundance_file <- "PICRUSt_analysis/required_files/pred_metagenome_unstrat.tsv" # read in abundance file (KO)
metadata <- read_delim(
  "required_files/artemis-eDNA-metadata-final.tsv",  # read in metadata file
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)
# test to make sure they sample names are the same
#list1 <- metadata$sample_name
#list2 <- colnames(abundance_data)
#list2[!(list1 %in% list2)] # output: [1] "STN115.35.fil.dura.r1"

metadata <- metadata %>% filter(., Sample.Control == "True.Sample") %>% filter(., sample_name != "STN115.35.fil.dura.r1") %>% 
  filter(., sample_name != "STN089.200.fil.dura.r2") %>% group_by(Station) #%>% distinct(Filter_pores, .keep_all = TRUE) 
# remove sample that isn't in kegg abundance for some reason

abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE) # read in KO dataset
ColumnstoKeep <- c("#NAME", metadata$sample_name) # set vector of list of names to keep + KO Name column
KO_abundance_data <- subset(abundance_data, select = ColumnstoKeep) # subset (select columns) based on ColumnstoKeep

#### Creating dataframe for all stations #######
# subset KO_abundance by KO_number in reference table
iron_KO <- read.csv("PICRUSt_analysis/required_files/KO_numbers_Sun_et_al.csv") # read in KO_number reference
iron_KO <- iron_KO %>% filter(., Metabolism == "Iron uptake and metabolism") # filter by iron metabolism only

KO_iron_numbers <- iron_KO$KO_Num # set vector for numbers
KO_iron_abundance_data <- KO_abundance_data[KO_abundance_data$`#NAME` %in% KO_iron_numbers, ] #filter by KO_number using KO_iron_number ref

KO_joined <- data.frame(iron_KO$Type, iron_KO$KO_Num)
colnames(KO_joined) <- c("Type", "#NAME")

KO_iron_abundance_type <- left_join(KO_iron_abundance_data, KO_joined) %>% column_to_rownames("#NAME")
test <- aggregate(. ~ Type, data = KO_iron_abundance_type, FUN = mean) # take mean per TYPE, for each station
# this one is for singular picrust (i.e. for Fe-S only)

######### doing relative abundance ########
KO_iron_abundance_type$Type <- NULL #remove for sum analysis
column_sums <- colSums(KO_iron_abundance_type) # take sum of each column
updated_abundance <- as.data.frame(column_sums) # create dataframe with sums
updated_abundance$sample_name <- rownames(updated_abundance) #create a column with rows (sample_names)
dataframe_to_combine <- metadata %>% select(sample_name, Depth_Threshold, Station, Filter_pores) # make temp dataframe for combining
all_stations <- left_join(updated_abundance, dataframe_to_combine) # join two dataframes by sample_name
all_stations$Filter_pores <- ifelse(all_stations$Filter_pores >= 0.2 & all_stations$Filter_pores <= 2.0, "free-living", 
                                        ifelse(all_stations$Filter_pores == 3.0, "particle-associated", all_stations$Filter_pores))

### set up for plots ###
level_order <- c("STN089", "STN132", "STN106", "STN20", "STN198", "STN002", "STN004", "STN012", "STN115", "STN12.3", "STN014", "STN078", "STN056a", "STN056b", "STN22", "STN068", "STN146", "STN181", "STN174", "STN151.2", "STN153")
depth_order <- c("Surface", "Intermediate", "Bottom_water")

#### plots #####
#### free - living ####
free_abundance <- all_stations %>% filter(Filter_pores == 0.2) # filter by freeicle-associated
free_abundance_agg <- aggregate(column_sums ~ Depth_Threshold * Station, data = free_abundance, FUN = mean) # take mean per TYPE, for each station
total_count <- sum(free_abundance_agg$column_sums)

# Calculate relative abundance
free_abundance_agg$relative_abundance <- free_abundance_agg$column_sums / total_count

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