# Install the devtools package if not already installed
# install.packages("devtools")

# Install ggpicrust2 from GitHub
# devtools::install_github("cafferychen777/ggpicrust2")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

# pkgs <- c("phyloseq", "ALDEx2", "SummarizedExperiment", "Biobase", "devtools", 
#          "ComplexHeatmap", "BiocGenerics", "BiocManager", "metagenomeSeq", 
#          "Maaslin2", "edgeR", "lefser", "limma", "KEGGREST", "DESeq2")

# for (pkg in pkgs) {
#  if (!requireNamespace(pkg, quietly = TRUE))
#   BiocManager::install(pkg)
# }

######################################## ggpicrust2 #################################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs <- c("phyloseq", "ALDEx2", "SummarizedExperiment", "Biobase", "devtools", 
          "ComplexHeatmap", "BiocGenerics", "BiocManager", "metagenomeSeq", 
          "Maaslin2", "edgeR", "lefser", "limma", "KEGGREST", "DESeq2")

for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggpubr)
library(ggprism)

library(KEGGREST)
library(patchwork)
require(grid)
#install.packages("IgAScores")
library(IgAScores) # gives relative abundance of count table
install.packages("ALDEx2")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ALDEx2")
library(ALDEx2)


abundance_file <- "PICRUSt_analysis/required_files/pred_metagenome_unstrat_descrip.tsv" # read in abundance file (KO)
metadata <- read_delim(
  "required_files/artemis-eDNA-metadata-final.tsv",  # read in outflow_metadata file
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)
# test to make sure they sample names are the same
#list1 <- outflow_metadata$sample_name
#list2 <- colnames(abundance_data)
#list2[!(list1 %in% list2)] # output: [1] "STN115.35.fil.dura.r1"

outflow_metadata <- outflow_metadata %>% filter(., Sample.Control == "True.Sample") %>% filter(., sample_name != "STN115.35.fil.dura.r1") %>% 
  filter(., sample_name != "STN089.200.fil.dura.r2") #%>% filter(., Transect_Name == "transect2") %>% filter(., Station != "STN22") #%>% distinct(Filter_pores, .keep_all = TRUE) 
# remove sample that isn't in kegg abundance for some reason

abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE) # read in KO dataset
ColumnstoKeep <- c("#NAME", outflow_metadata$sample_name) # set vector of list of names to keep + KO Name column
KO_abundance_data <- subset(abundance_data, select = ColumnstoKeep) # subset (select columns) based on ColumnstoKeep

#### Creating dataframe for all stations #######
# subset KO_abundance by KO_number in reference table
iron_KO <- read.csv("PICRUSt_analysis/required_files/KO_numbers_Sun_et_al.csv") # read in KO_number reference
iron_KO <- iron_KO %>% filter(., Metabolism == "Iron uptake and metabolism") # filter by iron metabolism only

KO_iron_numbers <- iron_KO$KO_Num # set vector for numbers
KO_iron_abundance_data <- KO_abundance_data[KO_abundance_data$`#NAME` %in% KO_iron_numbers, ] #filter by KO_number using KO_iron_number ref

######### for total rel abundance
# BASED ON KO NUMBER ONLY, not type
KO_iron_abundance_data_row <- KO_iron_abundance_data %>% column_to_rownames("#NAME")
KO_relative_abundance <- relabund(KO_iron_abundance_data_row, percentage = FALSE) # calculates relative abundance

KO_total_relative_abundance <- data.frame(KO_relative_abundance$sum, row.names=KO_iron_abundance_data$`#NAME`)
KO_total_relative_abundance <- relabund(KO_total_relative_abundance, percentage = FALSE) 

KO_total_relative_abundance <- arrange(KO_total_relative_abundance, KO_relative_abundance.sum)
KO_total_relative_abundance$pathway <- row.names(KO_total_relative_abundance)

your_order <- order(KO_total_relative_abundance$KO_relative_abundance.sum)
KO_total_relative_abundance$pathway <- factor(KO_total_relative_abundance$pathway, levels = KO_total_relative_abundance$pathway[your_order]) # setting order of pathways with the numerical value of the sum

p <- ggplot(data=KO_total_relative_abundance, aes(x=KO_relative_abundance.sum, y=pathway)) +
  geom_bar(stat="identity", width=0.5)

######################### ALL iron-related pathways #########################

# convert KO_iron_abundance_data to type
KO_joined <- data.frame(iron_KO$Type, iron_KO$KO_Num)
colnames(KO_joined) <- c("Type", "#NAME")

KO_iron_abundance_type <- left_join(KO_iron_abundance_data, KO_joined) %>% column_to_rownames("#NAME")
### all stations -- stacked barplot ##
test <- aggregate(. ~ Type, data = KO_iron_abundance_type, FUN = mean)
dataframe_to_combine <- outflow_metadata %>% select(sample_name, Transect_Number, Station, Filter_pores) 
df_test <- as.data.frame(t(test)) 
colnames(df_test) <- df_test[1,]
df_test <- df_test[-1,]
df_test$sample_name <- rownames(df_test)
test_up <- left_join(df_test, dataframe_to_combine)
test_up_long <- test_up %>% # combines by type
  pivot_longer(cols = -c(Station, Transect_Number, Filter_pores, sample_name),
               names_to = "Type",
               values_to = "Abundance")
test_up_long$Abundance <- as.numeric(test_up_long$Abundance)
agg_test <- aggregate(Abundance ~ Transect_Number * Filter_pores + Type, data = test_up_long, FUN = mean)
agg_test$Filter_pores <- ifelse(agg_test$Filter_pores >= 0.2 & agg_test$Filter_pores <= 2.0, "free-living", 
                                ifelse(agg_test$Filter_pores == 3.0, "particle-associated", agg_test$Filter_pores))

type_free <- agg_test %>% filter(Filter_pores == "free-living")

barplot_free <- ggplot(type_free, aes(x = Transect_Number, y = Abundance, fill = Type)) + # facet grid seperates by different levels, horizontally
  geom_bar(stat = "identity", position="fill", color= "black", linewidth=0.3, width=0.9) + theme_classic() + # adds black outline to boxes
  scale_y_continuous(expand = c(0, 0)) + # extends the barplots to the axies
  #scale_fill_manual(values = pathway_colors, drop = FALSE) + # set the colors with custom colors (myColors)
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=9, angle=90),
        axis.title.y = element_blank(),
        panel.spacing.x = unit(0, "points"), # Reducing space between facets
        strip.background = element_blank(), # Optionally hide the strip background for a cleaner look
        panel.border = element_blank()) + # Optionally remove panel borders
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # for the legend, if you want one
  #ylab("Relative Abundance (Order > 2%) \n") + # remove # if you want y-axis title
  ggtitle("Free-living (<0.2 µm)")

type_part <- agg_test %>% filter(Filter_pores == "particle-associated")
barplot_part <- ggplot(type_part, aes(x = Transect_Number, y = Abundance, fill = Type)) + # facet grid seperates by different levels, horizontally
  geom_bar(stat = "identity", position="fill", color= "black", linewidth=0.3, width=0.9) + theme_classic() + # adds black outline to boxes
  scale_y_continuous(expand = c(0, 0)) + # extends the barplots to the axies
  #scale_fill_manual(values = pathway_colors, drop = FALSE) + # set the colors with custom colors (myColors)
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=9, angle=90),
        axis.title.y = element_blank(),
        panel.spacing.x = unit(0, "points"), # Reducing space between facets
        strip.background = element_blank(), # Optionally hide the strip background for a cleaner look
        panel.border = element_blank()) + # Optionally remove panel borders
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # for the legend, if you want one
  #ylab("Relative Abundance (Order > 2%) \n") + # remove # if you want y-axis title
  ggtitle("Free-living (<0.2 µm)")


type_vector <- KO_iron_abundance_type$Type
KO_iron_abundance_type$Type <- NULL # temporarily remove column to figure out sums
KO_iron_abundance_type$sum <- rowSums(KO_iron_abundance_type)
KO_iron_abundance_type$Type <- type_vector

KO_total_data_frame <- data.frame(KO_iron_abundance_type$Type, KO_iron_abundance_type$sum)
colnames(KO_total_data_frame) <- c("Type", "sum")
KO_total_data_frame <- aggregate(. ~ Type, data = KO_total_data_frame, FUN = mean) # combines + summarizes

# calculate relative abundance for type
KO_total_data_frame_row <- KO_total_data_frame %>% column_to_rownames("Type")
KO_relative_abundance_type <- relabund(KO_total_data_frame_row, percentage = FALSE) # calculate rel. abundance

your_order <- order(KO_relative_abundance_type$sum)
KO_relative_abundance_type$Type <- rownames(KO_relative_abundance_type) #re-add type
KO_relative_abundance_type$Type <- factor(KO_relative_abundance_type$Type, levels = KO_relative_abundance_type$Type[your_order]) # setting order of pathways with the numerical value of the sum

p <- ggplot(data=KO_relative_abundance_type, aes(x=sum, y=Type)) +
  geom_bar(stat="identity", width=0.65) + ggtitle("All Communities") + ylab ("Iron-related capabilities \n") +
  xlab("Relative abundance")
ggsave("meltwater_total_iron_relative_abundance_plot.pdf", width=5, height=8)

#########################.      #########################

# Annotate pathway results without KO to KEGG conversion

metadata <- read_delim(
  "required_files/artemis-eDNA-metadata-final.tsv",  # read in outflow_metadata file
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

# Load KEGG pathway abundance
# data(kegg_abundance)
kegg_abundance <- ko2kegg_abundance(abundance_file) 
kegg_abundance$"#NAME" <- row.names(kegg_abundance)


outflow_metadata <- outflow_metadata %>% filter(., Sample.Control == "True.Sample") %>% filter(., sample_name != "STN115.35.fil.dura.r1") %>% 
  filter(., sample_name != "STN089.200.fil.dura.r2") #%>% filter(., Transect_Name == "transect2") %>% filter(., Station != "STN22") #%>% distinct(Filter_pores, .keep_all = TRUE) 
# remove sample that isn't in kegg abundance for some reason

abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE) # read in KO dataset
ColumnstoKeep <- c("#NAME", metadata$sample_name) # set vector of list of names to keep + KO Name column
KO_abundance_data <- subset(kegg_abundance, select = ColumnstoKeep) # subset (select columns) based on ColumnstoKeep
kegg_abundance$"#NAME" <- NULL


results_file_input <- ggpicrust2(data = abundance_data,
                                 metadata = metadata,
                                 group = "Location",
                                 pathway = "KO",
                                 daa_method = "LinDA",
                                 ko_to_kegg = TRUE,
                                 order = "pathway_class",
                                 p_values_bar = TRUE,
                                 x_lab = "pathway_name")

# Perform pathway differential abundance analysis (DAA) using ALDEx2 method
# Please change group to "your_group_column" if you are not using example dataset
metadata$Depth_Threshold <- as.character(metadata$Depth_Threshold)
metadata$Iron <- as.numeric(metadata$Iron)
rownames(kegg_abundance) <- NULL
daa_results_df <- pathway_daa(abundance = kegg_abundance %>% column_to_rownames("#NAME"), metadata = metadata, group = "Location", daa_method = "LinDA", select = NULL, p.adjust = "BH", reference = "Eastern_CC")

daa_results_df <- pathway_daa(abundance = kegg_abundance, metadata = metadata, group = "Location", daa_method = "ALDEx2", select = NULL, reference = NULL) 

data("metadata")

# Filter results for ALDEx2_Welch's t test method
# Please check the unique(daa_results_df$method) and choose one
daa_sub_method_results_df <- daa_results_df[daa_results_df$method == "LinDA", ]

# Annotate pathway results using KO to KEGG conversion
daa_annotated_sub_method_results_df <- pathway_annotation(pathway = "KO", daa_results_df = feature_with_p_0.05, ko_to_kegg = TRUE)


feature_with_p_0.05

feature_with_p_0.05 <- daa_results_df %>% 
  filter(p_adjust < 0.01)
data("metacyc_abundance")
kegg_abundance <- as_tibble(kegg_abundance)

# Create the heatmap
pathway_heatmap(abundance = kegg_abundance %>% filter("#NAME" %in% feature_with_p_0.05$feature) %>% column_to_rownames("#NAME"), metadata = metadata, group = "Location")







##### BOXPLOTS
station_type_abundance_mean <- aggregate(. ~ Type, data = KO_iron_abundance_type, FUN = mean) # organizes type by mean
station_type_abundance_mean$sum <- NULL
wanted_station <- data.frame(outflow_metadata$sample_name, outflow_metadata$Transect_Number, outflow_metadata$Filter_pores) 
colnames(wanted_station) <- c("sample_name", "Transect_Number", "Community")

station_type_abundance_intermediate <- t(station_type_abundance_mean)
colnames(station_type_abundance_intermediate) <- station_type_abundance_intermediate[1,] # make column names the Type
station_type_abundance_intermediate <- station_type_abundance_intermediate[-1, ] #remove duplicate first row
station_type_abundance_intermediate <- data.frame(station_type_abundance_intermediate)

#saved <- station_type_abundance_intermediate[39, ] 
#station_type_abundance_intermediate <- station_type_abundance_intermediate[-39, ]
#station_type_abundance_intermediate <- station_type_abundance_intermediate[-which(rownames(station_type_abundance_intermediate) == "sum"), ] #remove sum for analysis
station_type_abundance_intermediate$Transect_Number <- wanted_station$Transect_Number
station_type_abundance_intermediate$Community <- wanted_station$Community

station_type_abundance_intermediate <- mutate_all(station_type_abundance_intermediate, as.numeric)
station_type_abundance_intermediate$Transect_Number <- as.character(station_type_abundance_intermediate$Transect_Number)
station_type_abundance_intermediate$Community <- as.character(station_type_abundance_intermediate$Community)

station_type_abundance_intermediate_agg <- aggregate(. ~ Transect_Number * Community, data = station_type_abundance_intermediate, FUN = mean) 
#station_type_abundance_intermediate_agg$Transect_Number <- as.factor(station_type_abundance_intermediate_agg$Transect_Number)

transect_plot_df <- station_type_abundance_intermediate %>% 
  gather(variable, pathway_mean, -Transect_Number) %>%  
  mutate(Transect_Number = factor(Transect_Number, levels=unique(Transect_Number)))

transect_plot_df_community <- station_type_abundance_intermediate %>% 
  gather(variable, pathway_mean, -Community) %>% 
  mutate(Community = factor(Community, levels=unique(Community)))

transect_plot_df$Community <- transect_plot_df_community$Community 
transect_plot_df$pathway_mean <- as.numeric(transect_plot_df$pathway_mean)
transect_plot_df2 <- aggregate(pathway_mean ~ Transect_Number * variable * Community, data = transect_plot_df, FUN = mean) # take mean per TYPE, for each station
transect_plot_df2 <- transect_plot_df2 %>% filter(., variable != "Community")
total_count <- sum(transect_plot_df2$pathway_mean)
# Calculate relative abundance
transect_plot_df2$relative_abundance <- transect_plot_df2$pathway_mean / total_count

# variables for plots
remove_x_lab <- xlab("")
remove_y_lab <- ylab("")
bold_theme <- theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
                    axis.text.x = element_text(face="bold", size=10))
cols <- c("seagreen", "dodgerblue2")
add_color <- scale_fill_manual(values = cols, name = "Community", labels = c("Particle- \nassociated", "Free-living"))
add_outline <- scale_color_manual(values= c("seagreen", "dodgerblue2"))

# subset dataframe based on specific variables
Fe_S <- subset(transect_plot_df2, transect_plot_df2$variable == "Fe.S")
# bar
Fe_S_plot <- ggplot(Fe_S, aes(x = Transect_Number, y = relative_abundance, color = Community, fill = Community, group = factor(Community))) + bold_theme + add_color +
  geom_bar(stat="identity", position="dodge", width=0.65) + ggtitle("\n Fe-S") + remove_x_lab + remove_y_lab + add_outline + theme_classic2()
# line
Fe_S_plot <- ggplot(Fe_S, aes(x = Transect_Number, y = relative_abundance, color = Community, group = factor(Community))) + bold_theme +
  geom_point(size=1.5) + geom_line(linewidth=0.7) + ggtitle("Fe-S") + remove_x_lab + remove_y_lab + add_outline + theme_classic2()

# Siderophore uptake
Siderophore_uptake <- subset(transect_plot_df2, variable == "Siderophore.uptake")
Sid_uptake_plot <- ggplot(Siderophore_uptake, aes(x = Transect_Number, y = relative_abundance, color = Community, group = factor(Community))) +
  geom_point(size=1.5) + geom_line(linewidth=0.7) + 
  ggtitle("Siderophore uptake") + remove_x_lab + remove_y_lab + add_outline + theme_classic()

# Siderophore biosynthesis
sid_bio <- subset(transect_plot_df2, variable == "Siderophore.biosynthesis")
Sid_bio_plot <- ggplot(sid_bio, aes(x = Transect_Number, y = relative_abundance, color = Community, group = factor(Community))) +
  geom_point(size=1.5) + geom_line(linewidth=0.7) + 
  ggtitle("Siderophore biosynthesis") + remove_x_lab + remove_y_lab + add_outline + theme_classic()

# Fe3+ metabolism
ferric <- subset(transect_plot_df2, variable == "Fe3...ferric.")
ferric_plot <- ggplot(ferric, aes(x = Transect_Number, y = relative_abundance, color = Community, group = factor(Community))) +
  geom_point(size=1.5) + geom_line(linewidth=0.7) + 
  ggtitle("Fe3+ metabolism") + remove_x_lab + remove_y_lab + add_outline + theme_classic()

# Fur family (labeled incorrectly as Fe3+ uptake in the original question)
ferrous <- subset(transect_plot_df2, variable == "Fur.family")
ferrous_plot <- ggplot(ferrous, aes(x = Transect_Number, y = relative_abundance, color = Community, group = factor(Community))) +
  geom_point(size=1.5) + geom_line(linewidth=0.7) + 
  ggtitle("Fur family") + remove_x_lab + remove_y_lab + add_outline + theme_classic()

# Fe2+ metabolism
ferric_uptake <- subset(transect_plot_df2, variable == "Fe2...ferrous.")
ferric_uptake_plot <- ggplot(ferric_uptake, aes(x = Transect_Number, y = relative_abundance, color = Community, group = factor(Community))) +
  geom_point(size=1.5) + geom_line(linewidth=0.7) + 
  ggtitle("Fe2+ metabolism") + remove_x_lab + remove_y_lab + add_outline + theme_classic()

# Fe storage
storage <- subset(transect_plot_df2, variable == "Fe.Storage.")
storage_plot <- ggplot(storage, aes(x = Transect_Number, y = relative_abundance, color = Community, group = factor(Community))) +
  geom_point(size=1.5) + geom_line(linewidth=0.7) + 
  ggtitle("Fe storage") + remove_x_lab + remove_y_lab + add_outline + theme_classic()

# Heme
heme <- subset(transect_plot_df2, variable == "Heme")
heme_plot <- ggplot(heme, aes(x = Transect_Number, y = relative_abundance, color = Community, group = factor(Community))) +
  geom_point(size=1.5) + geom_line(linewidth=0.7) + 
  ggtitle("Heme") + remove_x_lab + remove_y_lab + add_outline + theme_classic()

# troR
tror <- subset(transect_plot_df2, variable == "troR")
tror_plot <- ggplot(tror, aes(x = Transect_Number, y = relative_abundance, color = Community, group = factor(Community))) +
  geom_point(size=1.5) + geom_line(linewidth=0.7) + 
  ggtitle("Iron homeostasis") + remove_x_lab + remove_y_lab + add_outline + theme_classic()

# actual boxplot
boxplot <- ggarrange(Sid_uptake_plot, ferric_plot, ferrous_plot, storage_plot, heme_plot, tror_plot, Fe_S_plot, Sid_bio_plot,
          ncol = 4, nrow = 2, common.legend = TRUE, legend = "right")
annotate_figure(boxplot, left = textGrob("Pathway abundance mean", rot = 90, vjust = 1, gp = gpar(cex = 1)), # adding x and y to actual annotations
                bottom = textGrob("Transect number", gp = gpar(cex = 1), vjust=.2))
ggsave("meltwater_iron_lineplots.pdf", width=10, height=8) # save as pdf

variable_list <- c("Fe.S", "Siderophore.uptake", "Siderophore.biosynthesis", "Fe3...ferric.", "Fur.family", "Fe2...ferrous.", "Fe.Storage.", "Heme")
subset_transect_plot_df <- transect_plot_df2[transect_plot_df2$variable %in% variable_list,]

CombinedPlotNOTCH=ggplot(subset_transect_plot_df, aes(x=Transect_Number, y=relative_abundance, color=Community, group = factor(Community))) + geom_point() +
  geom_line() + facet_grid(~variable)
CombinedPlotNOTCH
ggsave("ONE_PLOT_changed_transect_iron_line.pdf", width=10, height=8)


#### stacked barplot
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




### trying LINDA anaylsis
linda_station_type_abundance_mean <- station_type_abundance_mean
rownames(linda_station_type_abundance_mean) <- linda_station_type_abundance_mean$Type
linda_station_type_abundance_mean$Type <- NULL

linda.obj <- linda(linda_station_type_abundance_mean, outflow_metadata, formula = "~Transect_Number*Filter_pores",  feature.dat.type="proportion", p.adj.method = "BH",)
linda_plot <- linda.plot(linda.obj, c("Transect_Number"),
           titles = c('Transect # AND Community Diff Abund.'), alpha = 0.1, lfc.cut = 1,
           legend = TRUE, directory = NULL, width = 11, height = 8)

plot_print <- linda_plot[["plot.lfc"]]
ggsave("linda_outflow_Conf.pdf", width=10, height=6)

daa_results_df_transect <- pathway_daa(abundance = station_type_abundance_mean %>% column_to_rownames("Type"), outflow_metadata = outflow_metadata, group = "Transect_Number", daa_method = "LinDA", select = NULL, p.adjust = "BH", reference = "1")
daa_results_df_community <- pathway_daa(abundance = station_type_abundance_mean %>% column_to_rownames("Type"), outflow_metadata = outflow_metadata, group = "Filter_pores", daa_method = "LinDA", select = NULL, p.adjust = "BH", reference = "1")
