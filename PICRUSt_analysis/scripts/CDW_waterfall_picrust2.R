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

abundance_file <- "PICRUSt_analysis/required_files/pred_metagenome_unstrat.tsv" # read in abundance file (KO)
metadata <- read_delim(
  "required_files/artemis-eDNA-metadata-final.tsv",  # read in metadata file
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

list1 <- metadata$sample_name
list2 <- colnames(abundance_data)
list2[!(list1 %in% list2)] # output: [1] "STN115.35.fil.dura.r1"

metadata <- metadata %>% filter(., Sample.Control == "True.Sample") %>% filter(., sample_name != "STN115.35.fil.dura.r1") 
# remove sample that isn't in kegg abundance for some reason
inflow_metadata <- metadata %>% filter(., Transect_Name == "transect1") %>% filter(., Station != "STN22") #remove station 22
# for some reason it wont filter by a string or more than two samples, so had to do separately.

abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE) # read in KO dataset
ColumnstoKeep <- c("#NAME", inflow_metadata$sample_name) # set vector of list of names to keep + KO Name column
KO_abundance_data <- subset(abundance_data, select = ColumnstoKeep) # subset (select columns) based on ColumnstoKeep

#### following pipeline from picrust github ### inflow ONLY
#write.table(abundance_final, file='inflow_abundance.tsv', quote=TRUE, sep='\t', row.names=FALSE) #slay , works

# subset KO_abundance by KO_number in reference table
iron_KO <- read.csv("PICRUSt_analysis/required_files/KO_numbers_Sun_et_al.csv") # read in KO_number reference
iron_KO <- iron_KO %>% filter(., Metabolism == "Iron uptake and metabolism") # filter by iron metabolism only

KO_iron_numbers <- iron_KO$KO_Num # set vector
KO_iron_abundance_data <- KO_abundance_data[KO_abundance_data$`#NAME` %in% KO_iron_numbers, ] #filter by KO_number using KO_iron_number ref

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
ggsave("CDW_waterfall_total_iron_relative_abundance_plot.pdf", width=5, height=8)

#########################.      #########################

# Annotate pathway results without KO to KEGG conversion
daa_annotated_sub_method_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df, ko_to_kegg = TRUE)

##### BOXPLOTS
station_type_abundance_mean <- aggregate(. ~ Type, data = KO_iron_abundance_type, FUN = mean) # organizes type by mean
station_type_abundance_mean$sum <- NULL
wanted_station <- data.frame(inflow_metadata$sample_name, inflow_metadata$Transect_Number, inflow_metadata$Filter_pores) 
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

#station_type_abundance_intermediate_agg <- aggregate(. ~ Transect_Number, data = station_type_abundance_intermediate, FUN = mean) 
#station_type_abundance_intermediate_agg$Transect_Number <- as.factor(station_type_abundance_intermediate_agg$Transect_Number)

transect_plot_df <- station_type_abundance_intermediate %>% 
  gather(variable, pathway_mean, -Transect_Number) %>% 
  mutate(Transect_Number = factor(Transect_Number, levels=unique(Transect_Number)))

transect_plot_df_community <- station_type_abundance_intermediate %>% 
  gather(variable, pathway_mean, -Community) %>% 
  mutate(Community = factor(Community, levels=unique(Community)))

transect_plot_df$Community <- transect_plot_df_community$Community 
transect_plot_df$pathway_mean <- as.numeric(transect_plot_df$pathway_mean)
transect_plot_df$Transect_Number <- factor(transect_plot_df$Transect_Number, levels=c("1", "2", "3", "4", "5"), ordered = TRUE) # order transect numbers for the boxplots

# variables for plots
remove_x_lab <- xlab("")
remove_y_lab <- ylab("")
bold_theme <- theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
                    axis.text.x = element_text(face="bold", size=10))
cols <- c("darkseagreen", "dodgerblue2")
add_color <- scale_fill_manual(values = cols, name = "Community", labels = c("Particle- \nassociated", "Free-living"))
add_outline <- scale_color_manual(values= c("darkgreen", "darkblue"))

# subset dataframe based on specific variables
Fe_S <- subset(transect_plot_df, transect_plot_df$variable == "Fe.S")
Fe_S_plot <- ggplot(Fe_S, aes(x = Transect_Number, y = pathway_mean, fill = Community, color = Community)) + bold_theme +
  geom_boxplot(outlier.shape=16, outlier.size=1) + ggtitle("\n Fe-S") + remove_x_lab + remove_y_lab + add_outline + add_color 

Siderophore_uptake <- subset(transect_plot_df, transect_plot_df$variable == "Siderophore.uptake")
Sid_uptake_plot <- ggplot(Siderophore_uptake, aes(x = Transect_Number, y = pathway_mean, fill = Community, color = Community)) + bold_theme +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1) + ggtitle("Siderophore \n uptake") + remove_x_lab + remove_y_lab + add_outline + add_color

sid_bio <- subset(transect_plot_df, transect_plot_df$variable == "Siderophore.biosynthesis")
Sid_bio_plot <- ggplot(sid_bio, aes(x = Transect_Number, y = pathway_mean, fill = Community, color = Community)) + bold_theme +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1) + ggtitle("Siderophore \n biosynthesis") + remove_x_lab + remove_y_lab + add_outline + add_color

ferric <- subset(transect_plot_df, transect_plot_df$variable == "Fe3...ferric.")
ferric_plot <- ggplot(ferric, aes(x = Transect_Number, y = pathway_mean, fill = Community, color = Community)) + bold_theme +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1) + ggtitle("Fe3+ \n metabolism") + remove_x_lab + remove_y_lab + add_outline + add_color

ferrous <- subset(transect_plot_df, transect_plot_df$variable == "Fur.family")
ferrous_plot <- ggplot(ferrous, aes(x = Transect_Number, y = pathway_mean, fill = Community, color = Community)) + bold_theme +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1) + ggtitle("\n Fe3+ uptake") + remove_x_lab + remove_y_lab + add_outline + add_color

ferric_uptake <- subset(transect_plot_df, transect_plot_df$variable == "Fe2...ferrous.")
ferric_uptake_plot <- ggplot(ferric_uptake, aes(x = Transect_Number, y = pathway_mean, fill = Community, color = Community)) + bold_theme +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1) + ggtitle("Fe2+ metabolism") + remove_x_lab + remove_y_lab + add_outline + add_color

storage <- subset(transect_plot_df, transect_plot_df$variable == "Fe.Storage.")
storage_plot <- ggplot(storage, aes(x = Transect_Number, y = pathway_mean, fill = Community, color = Community)) + bold_theme +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1) + ggtitle("\n Fe storage") + remove_x_lab + remove_y_lab + add_outline + add_color

heme <- subset(transect_plot_df, transect_plot_df$variable == "Heme")
heme_plot <- ggplot(heme, aes(x = Transect_Number, y = pathway_mean, fill = Community, color = Community)) + bold_theme +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1) + ggtitle("\n Heme") + remove_x_lab + remove_y_lab + add_outline + add_color

tror <- subset(transect_plot_df, transect_plot_df$variable == "troR")
tror_plot <- ggplot(heme, aes(x = Transect_Number, y = pathway_mean, fill = Community, color = Community)) + bold_theme +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1) + ggtitle("Iron \n homeostasis") + remove_x_lab + remove_y_lab + add_outline + add_color

# actual boxplot
boxplot <- ggarrange(Sid_uptake_plot, ferric_plot, ferrous_plot, storage_plot, heme_plot, tror_plot, Fe_S_plot, Sid_bio_plot,
                     ncol = 4, nrow = 2, common.legend = TRUE, legend = "right")
annotate_figure(boxplot, left = textGrob("Pathway abundance mean", rot = 90, vjust = 1, gp = gpar(cex = 1)), # adding x and y to actual annotations
                bottom = textGrob("Transect number", gp = gpar(cex = 1), vjust=.2))
ggsave("CDW_waterfall_iron_boxplots.pdf", width=10, height=8) # save as pdf

variable_list <- c("Fe.S", "Siderophore.uptake", "Siderophore.biosynthesis", "Fe3...ferric.", "Fur.family", "Fe2...ferrous.", "Fe.Storage.", "Heme")
subset_transect_plot_df <- transect_plot_df[transect_plot_df$variable %in% variable_list,]

CombinedPlotNOTCH=ggplot(subset_transect_plot_df, aes(x=Transect_Number, y=pathway_mean, fill=Community)) + geom_boxplot() + facet_grid(~variable)
CombinedPlotNOTCH
ggsave("CDW_waterfall_all_boxplots.pdf", width=10, height=8)

### trying LINDA anaylsis
linda_station_type_abundance_mean <- station_type_abundance_mean
rownames(linda_station_type_abundance_mean) <- linda_station_type_abundance_mean$Type
linda_station_type_abundance_mean$Type <- NULL
library(MicrobiomeStat)
linda.obj <- linda(linda_station_type_abundance_mean, inflow_metadata, formula = "~Transect_Number+Filter_pores",  feature.dat.type="proportion", p.adj.method = "BH",)
linda_plot <- linda.plot(linda.obj, c("Transect_Number"),
                         titles = c('Transect # AND Community Diff Abund.'), alpha = 0.1, lfc.cut = 1,
                         legend = TRUE, directory = NULL, width = 11, height = 8)

plot_print <- linda_plot[["plot.lfc"]]
ggsave("linda_inflow_Conf.pdf", width=10, height=6)

daa_results_df_transect <- pathway_daa(abundance = station_type_abundance_mean %>% column_to_rownames("Type"), metadata = inflow_metadata, group = "Transect_Number", daa_method = "LinDA", select = NULL, p.adjust = "BH", reference = "1")
daa_results_df_community <- pathway_daa(abundance = station_type_abundance_mean %>% column_to_rownames("Type"), metadata = inflow_metadata, group = "Filter_pores", daa_method = "LinDA", select = NULL, p.adjust = "BH", reference = "1")
