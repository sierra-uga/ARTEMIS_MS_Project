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

abundance_file <- "pred_metagenome_unstrat.tsv" # read in abundance file (KO)
metadata <- read_delim(
  "artemis-eDNA-metadata-final.tsv",  # read in metadata file
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

list1 <- metadata$sample_name
list2 <- colnames(abundance_data)
list2[!(list1 %in% list2)] # output: [1] "STN115.35.fil.dura.r1"

metadata <- metadata %>% filter(., Sample.Control == "True.Sample") %>% filter(., sample_name != "STN115.35.fil.dura.r1") %>% filter(., sample_name != "STN089.200.fil.dura.r2")
# remove sample that isn't in kegg abundance for some reason
coastal_metadata <- metadata %>% filter(., Coastal_Current_Name == "transect3")
coastal_metadata_unique <- coastal_metadata %>%
  group_by(Station) %>%
  distinct(Filter_pores, .keep_all = TRUE)

abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE) # read in KO dataset
ColumnstoKeep <- c("#NAME", coastal_metadata_unique$sample_name) # set vector of list of names to keep + KO Name column
KO_abundance_data <- subset(abundance_data, select = ColumnstoKeep) # subset (select columns) based on ColumnstoKeep

#### following pipeline from picrust github ### coastal ONLY
#write.table(abundance_final, file='coastal_abundance.tsv', quote=TRUE, sep='\t', row.names=FALSE) #slay , works

# subset KO_abundance by KO_number in reference table
iron_KO <- read.csv("KO_numbers_Sun_et_al.csv") # read in KO_number reference
iron_KO <- iron_KO %>% filter(., Metabolism == "Iron uptake and metabolism") # filter by iron metabolism only

KO_iron_numbers <- iron_KO$KO_Num # set vector
KO_iron_abundance_data <- KO_abundance_data[KO_abundance_data$`#NAME` %in% KO_iron_numbers, ] #filter by KO_number using KO_iron_number ref

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
ggsave("coastal_current_total_iron_relative_abundance_plot.pdf", width=5, height=8)

#########################.      #########################

# Annotate pathway results without KO to KEGG conversion
daa_annotated_sub_method_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df, ko_to_kegg = TRUE)

##### BOXPLOTS
station_type_abundance_mean <- aggregate(. ~ Type, data = KO_iron_abundance_type, FUN = mean) # organizes type by mean
station_type_abundance_mean$sum <- NULL
wanted_station <- data.frame(coastal_metadata_unique$sample_name, coastal_metadata_unique$Coastal_Current_Number, coastal_metadata_unique$Filter_pores) 
colnames(wanted_station) <- c("sample_name", "Coastal_Current_Number", "Community")

station_type_abundance_intermediate <- t(station_type_abundance_mean)
colnames(station_type_abundance_intermediate) <- station_type_abundance_intermediate[1,] # make column names the Type
station_type_abundance_intermediate <- station_type_abundance_intermediate[-1, ] #remove duplicate first row
station_type_abundance_intermediate <- data.frame(station_type_abundance_intermediate)

#saved <- station_type_abundance_intermediate[39, ] 
#station_type_abundance_intermediate <- station_type_abundance_intermediate[-39, ]
#station_type_abundance_intermediate <- station_type_abundance_intermediate[-which(rownames(station_type_abundance_intermediate) == "sum"), ] #remove sum for analysis
station_type_abundance_intermediate$Coastal_Current_Number <- wanted_station$Coastal_Current_Number
station_type_abundance_intermediate$Community <- wanted_station$Community

station_type_abundance_intermediate <- mutate_all(station_type_abundance_intermediate, as.numeric)
station_type_abundance_intermediate$Coastal_Current_Number <- as.character(station_type_abundance_intermediate$Coastal_Current_Number)
station_type_abundance_intermediate$Community <- as.character(station_type_abundance_intermediate$Community)

#station_type_abundance_intermediate_agg <- aggregate(. ~ Coastal_Current_Number, data = station_type_abundance_intermediate, FUN = mean) 
#station_type_abundance_intermediate_agg$Coastal_Current_Number <- as.factor(station_type_abundance_intermediate_agg$Coastal_Current_Number)

transect_plot_df <- station_type_abundance_intermediate %>% 
  gather(variable, pathway_mean, -Coastal_Current_Number) %>% 
  mutate(Coastal_Current_Number = factor(Coastal_Current_Number, levels=unique(Coastal_Current_Number)))

transect_plot_df_community <- station_type_abundance_intermediate %>% 
  gather(variable, pathway_mean, -Community) %>% 
  mutate(Community = factor(Community, levels=unique(Community)))

transect_plot_df_community$pathway_mean <- as.numeric(transect_plot_df_community$pathway_mean)
ggplot(transect_plot_df_community, aes(x = pathway_mean, y = variable, fill = Community)) +
  geom_bar(stat = "identity") +
  labs(x = "Pathway Mean", y = "Variable", title = "Stacked Bar Plot") +
  theme_minimal()

transect_plot_df$Community <- transect_plot_df_community$Community 
transect_plot_df$pathway_mean <- as.numeric(transect_plot_df$pathway_mean)
transect_plot_df$Coastal_Current_Number <- factor(transect_plot_df$Coastal_Current_Number, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9"), ordered = TRUE) # order transect numbers for the boxplots

# variables for plots
remove_x_lab <- xlab("")
remove_y_lab <- ylab("")
bold_theme <- theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
                    axis.text.x = element_text(face="bold", size=10))
cols <- c("darkseagreen", "dodgerblue2")
add_color <- scale_fill_manual(values = cols, name = "Community", labels = c("Particle- \nassociated", "Free-living"))
add_outline <- scale_color_manual(values= c("darkgreen", "darkblue"))

variable_list <- c("Fe.S", "Siderophore.uptake", "Siderophore.biosynthesis", "Fe3...ferric.", "Fur.family", "troR", "Fe.Storage.", "Heme")
title_list <- c("Fe-S", "Siderophore uptake", "Siderophore biosynthesis", "Fe3+ metabolism", "Fur", "Fe homeostasis", "Fe storage", "Heme")  # Add your titles here
plots <- list()

# Loop through each unique combination and create a plot
for (i in seq_along(variable_list)) {
  variable <- variable_list[i]
  title <- title_list[i]
  
  # Subset transect_plot_df based on the current variable
  subset_data <- transect_plot_df[transect_plot_df$variable == variable, ]
  
  # Create a plot title
  plot_title <- paste("\n", title, sep = "")
  
  # Create the plot
  plot <- ggplot(subset_data, aes(x = Coastal_Current_Number, y = pathway_mean, fill = Community, color = Community)) +
    geom_boxplot(outlier.colour = "black", outlier.shape = 16, outlier.size = 1) +
    ggtitle(plot_title) +
    bold_theme +
    remove_x_lab +
    remove_y_lab +
    add_outline +
    add_color 
  
  # Add plot to the list
  plots[[title]] <- plot
}

boxplot_new <- ggarrange(plotlist = plots, ncol=4, nrow=2, common.legend = TRUE, legend = "right")
annotate_figure(boxplot, left = textGrob("Pathway abundance mean", rot = 90, vjust = 1, gp = gpar(cex = 1)), # adding x and y to actual annotations
                bottom = textGrob("Transect number", gp = gpar(cex = 1), vjust=.2))

##########################

# actual boxplot
boxplot <- ggarrange(Sid_uptake_plot, ferric_plot, ferrous_plot, storage_plot, heme_plot, tror_plot, Fe_S_plot, Sid_bio_plot,
                     ncol = 4, nrow = 2, common.legend = TRUE, legend = "right")
annotate_figure(boxplot, left = textGrob("Pathway abundance mean", rot = 90, vjust = 1, gp = gpar(cex = 1)), # adding x and y to actual annotations
                bottom = textGrob("Transect number", gp = gpar(cex = 1), vjust=.2))
ggsave("coastal_current_iron_boxplots.pdf", width=10, height=8) # save as pdf

variable_list <- c("Fe.S", "Siderophore.uptake", "Siderophore.biosynthesis", "Fe3...ferric.", "Fur.family", "Fe2...ferrous.", "Fe.Storage.", "Heme")
subset_transect_plot_df <- transect_plot_df[transect_plot_df$variable %in% variable_list,]

CombinedPlotNOTCH=ggplot(subset_transect_plot_df, aes(x=Coastal_Current_Number, y=pathway_mean, fill=Community)) + geom_boxplot() + facet_grid(~variable)
CombinedPlotNOTCH
ggsave("coastal_current_all_boxplots.pdf", width=10, height=8)

### trying LINDA anaylsis
linda_station_type_abundance_mean <- station_type_abundance_mean
rownames(linda_station_type_abundance_mean) <- linda_station_type_abundance_mean$Type
linda_station_type_abundance_mean$Type <- NULL
library(MicrobiomeStat)
linda.obj <- linda(linda_station_type_abundance_mean, coastal_metadata, formula = "~Coastal_Current_Number+Filter_pores",  feature.dat.type="proportion", p.adj.method = "BH",)
linda_plot <- linda.plot(linda.obj, c("Coastal_Current_Number"),
                         titles = c('Transect # AND Community Diff Abund.'), alpha = 0.1, lfc.cut = 1,
                         legend = TRUE, directory = NULL, width = 11, height = 8)

plot_print <- linda_plot[["plot.lfc"]]
ggsave("linda_coastal_Conf.pdf", width=10, height=6)

daa_results_df_transect <- pathway_daa(abundance = station_type_abundance_mean %>% column_to_rownames("Type"), metadata = coastal_metadata, group = "Coastal_Current_Number", daa_method = "LinDA", select = NULL, p.adjust = "BH", reference = "1")
daa_results_df_community <- pathway_daa(abundance = station_type_abundance_mean %>% column_to_rownames("Type"), metadata = coastal_metadata, group = "Filter_pores", daa_method = "LinDA", select = NULL, p.adjust = "BH", reference = "1")
