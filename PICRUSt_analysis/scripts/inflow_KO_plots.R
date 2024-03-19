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
library(ggprism)
library(patchwork)
#install.packages("IgAScores")
library(IgAScores) # gives relative abundance of count table

abundance_file <- "pred_metagenome_unstrat.tsv" # read in abundance file (KO)
metadata <- read_delim(
  "artemis-eDNA-metadata-final.tsv",  # read in metadata file
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)

abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE) # read in KO dataset

list1 <- metadata$sample_name
list2 <- colnames(abundance_data)
list2[!(list1 %in% list2)] # output: [1] "STN115.35.fil.dura.r1"

metadata <- metadata %>% filter(., Sample.Control == "True.Sample") %>% filter(., sample_name != "STN115.35.fil.dura.r1") 
# remove sample that isn't in kegg abundance for some reason
inflow_metadata <- metadata %>% filter(., Transect_Name == "transect1")
# for some reason it wont filter by a string or more than two samples, so had to do separately.
ColumnstoKeep <- c("#NAME", inflow_metadata$sample_name) # set vector of list of names to keep + KO Name column
KO_abundance_data <- subset(abundance_data, select = ColumnstoKeep) # subset (select columns) based on ColumnstoKeep

#### following pipeline from picrust github ### inflow ONLY
#write.table(abundance_final, file='inflow_abundance.tsv', quote=TRUE, sep='\t', row.names=FALSE) #slay , works

# subset KO_abundance by KO_number in reference table
iron_KO <- read.csv("KO_numbers_Sun_et_al.csv") # read in KO_number reference
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
  geom_bar(stat="identity", width=0.95) + ggtitle("All Communities")
ggsave("changed_total_iron_relative_abundance_plot.pdf", width=5, height=8)

#########################.      #########################

# Annotate pathway results without KO to KEGG conversion
daa_annotated_sub_method_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df, ko_to_kegg = TRUE)

##### BOXPLOTS
station_type_abundance_mean <- aggregate(. ~ Type, data = KO_iron_abundance_type, FUN = mean) # organizes type by mean
station_type_abundance_mean$sum <- NULL
wanted_station <- data.frame(inflow_metadata$sample_name, inflow_metadata$Transect_Number, inflow_metadata$Filter_pores) 
colnames(wanted_station) <- c("sample_name", "Transect_Number", "Community")

station_type_abundance_intermediate <- t(station_type_abundance_mean)
station_type_abundance_intermediate <- data.frame(station_type_abundance_intermediate)
colnames(station_type_abundance_intermediate) <- station_type_abundance_intermediate[1,] # make column names the Type
station_type_abundance_intermediate <- station_type_abundance_intermediate[-1, ] #remove duplicate first row

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

Fe_S <- subset(transect_plot_df, transect_plot_df$variable == "Fe-S")
Fe_S_plot <- ggplot(Fe_S, aes(x = Transect_Number, y = pathway_mean, fill = Community)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2) + ggtitle("Fe-S pathways") 

Siderophore_uptake <- subset(transect_plot_df, transect_plot_df$variable == "Siderophore uptake")
Sid_uptake_plot <- ggplot(Siderophore_uptake, aes(x = Transect_Number, y = pathway_mean, fill = Community)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2) + ggtitle("Siderophore Uptake")

sid_bio <- subset(transect_plot_df, transect_plot_df$variable == "Siderophore biosynthesis")
Sid_bio_plot <- ggplot(sid_bio, aes(x = Transect_Number, y = pathway_mean, fill = Community)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2) + ggtitle("Siderophore biosynthesis")

ferric <- subset(transect_plot_df, transect_plot_df$variable == "Fe3+ (ferric)")
ferric_plot <- ggplot(ferric, aes(x = Transect_Number, y = pathway_mean, fill = Community)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2) + ggtitle("Fe3+ (ferric) metabolism")

ferric_uptake <- subset(transect_plot_df, transect_plot_df$variable == "Fur family")
ferric_uptake_plot <- ggplot(ferric_uptake, aes(x = Transect_Number, y = pathway_mean, fill = Community)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2) + ggtitle("Fe3+ (ferric) uptake")

ferrous <- subset(transect_plot_df, transect_plot_df$variable == "Fe2+ (ferrous)")
ferrous_plot <- ggplot(ferrous, aes(x = Transect_Number, y = pathway_mean, fill = Community)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2) + ggtitle("Fe2+ (ferrous) metabolism")

ggarrange(Fe_S_plot, Sid_uptake_plot, Sid_bio_plot, ferric_plot, ferric_uptake_plot, ferrous_plot,
          ncol = 3, nrow = 2, common.legend = TRUE, legend = "right")
ggsave("changed_transect_iron_boxplots_INFLOW.pdf", width=10, height=8)

storage <- subset(transect_plot_df, transect_plot_df$variable == "Fe Storage ")
storage_plot <- ggplot(storage, aes(x = Transect_Number, y = pathway_mean, fill = Transect_Number)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2) + ggtitle("Fe storage")

heme <- subset(transect_plot_df, transect_plot_df$variable == "Heme")
heme_plot <- ggplot(heme, aes(x = Transect_Number, y = pathway_mean, fill = Community)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2) + ggtitle("Heme")

variable_list <- unique(transect_plot_df$variable)
subset_transect_plot_df <- transect_plot_df[transect_plot_df$variable %in% variable_list,]

CombinedPlotNOTCH=ggplot(subset_transect_plot_df, aes(x=Transect_Number, y=pathway_mean, fill=Community)) + geom_boxplot() + facet_grid(~variable)
CombinedPlotNOTCH
ggsave("ONE_PLOT_changed_transect_iron_boxplots_inflow.pdf", width=10, height=8)

### trying LINDA anaylsis
linda_station_type_abundance_mean <- station_type_abundance_mean
rownames(linda_station_type_abundance_mean) <- linda_station_type_abundance_mean$Type
linda_station_type_abundance_mean$Type <- NULL

linda.obj <- linda(linda_station_type_abundance_mean, inflow_metadata, formula = "~Transect_Number*Filter_pores",  feature.dat.type="proportion", p.adj.method = "BH",)
linda_plot <- linda.plot(linda.obj, c("Filter_pores3", "Transect_Number:Filter_pores3"),
                         titles = c('Transect Number AND Community log2fold Change'), alpha = 0.1, lfc.cut = 1,
                         legend = TRUE, directory = NULL, width = 11, height = 8)

plot_print <- linda_plot[["plot.lfc"]]

daa_results_df_transect <- pathway_daa(abundance = station_type_abundance_mean %>% column_to_rownames("Type"), metadata = inflow_metadata, group = "Transect_Number", daa_method = "LinDA", select = NULL, p.adjust = "BH", reference = "1")
daa_results_df_community <- pathway_daa(abundance = station_type_abundance_mean %>% column_to_rownames("Type"), metadata = inflow_metadata, group = "Filter_pores", daa_method = "LinDA", select = NULL, p.adjust = "BH", reference = "1")
