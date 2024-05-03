######################################## ggpicrust2 #################################################
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggpubr)
library(ggprism)
library(patchwork)
require(grid)
#install.packages("IgAScores")
library(IgAScores) # gives relative abundance of count table

abundance_file <- "PICRUSt_analysis/required_files/pred_metagenome_unstrat.tsv" # read in abundance file (KO)
coastal_metadata <- read_delim(
  "required_files/artemis-eDNA-metadata-final.tsv",  # read in coastal_metadata file
  delim = "\t",
  escape_double = FALSE,
  trim_ws = TRUE
)
# test to make sure they sample names are the same
#list1 <- coastal_metadata$sample_name
#list2 <- colnames(abundance_data)
#list2[!(list1 %in% list2)] # output: [1] "STN115.35.fil.dura.r1"

coastal_metadata <- coastal_metadata %>% filter(., Sample.Control == "True.Sample") %>% filter(., sample_name != "STN115.35.fil.dura.r1") %>% 
  filter(., sample_name != "STN089.200.fil.dura.r2") %>% filter(., Coastal_Current_Name == "transect3") #%>% distinct(Filter_pores, .keep_all = TRUE) 
# remove sample that isn't in kegg abundance for some reason

abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE) # read in KO dataset
ColumnstoKeep <- c("#NAME", coastal_metadata$sample_name) # set vector of list of names to keep + KO Name column
KO_abundance_data <- subset(abundance_data, select = ColumnstoKeep) # subset (select columns) based on ColumnstoKeep

#### Creating dataframe for all stations #######
# subset KO_abundance by KO_number in reference table
iron_KO <- read.csv("PICRUSt_analysis/required_files/KO_numbers_Sun_et_al.csv") # read in KO_number reference
iron_KO <- iron_KO %>% filter(., Metabolism == "Iron uptake and metabolism") # filter by iron metabolism only

KO_iron_numbers <- iron_KO$KO_Num # set vector for numbers
KO_iron_abundance_data <- KO_abundance_data[KO_abundance_data$`#NAME` %in% KO_iron_numbers, ] #filter by KO_number using KO_iron_number ref

# convert KO_iron_abundance_data to type
KO_joined <- data.frame(iron_KO$Type, iron_KO$KO_Num)
colnames(KO_joined) <- c("Type", "#NAME")

KO_iron_abundance_type <- left_join(KO_iron_abundance_data, KO_joined) %>% column_to_rownames("#NAME")
### all stations -- stacked barplot ##
test <- aggregate(. ~ Type, data = KO_iron_abundance_type, FUN = mean)
dataframe_to_combine <- coastal_metadata %>% select(sample_name, Coastal_Current_Number, Station, Filter_pores, Depth_Threshold) 
df_test <- as.data.frame(t(test)) 
colnames(df_test) <- df_test[1,]
df_test <- df_test[-1,]
df_test$sample_name <- rownames(df_test)
test_up <- left_join(df_test, dataframe_to_combine)
test_up_long <- test_up %>% # combines by type
  pivot_longer(cols = -c(Station, Coastal_Current_Number, Filter_pores, Depth_Threshold, sample_name),
               names_to = "Type",
               values_to = "Abundance")
test_up_long$Abundance <- as.numeric(test_up_long$Abundance)
agg_test <- aggregate(Abundance ~ Coastal_Current_Number * Filter_pores * Type * Depth_Threshold, data = test_up_long, FUN = mean)
agg_test$Filter_pores <- ifelse(agg_test$Filter_pores >= 0.2 & agg_test$Filter_pores <= 2.0, "free-living", 
                                ifelse(agg_test$Filter_pores == 3.0, "particle-associated", agg_test$Filter_pores))

wanted_types <- c("Fe Storage ", "Fe-S", "Fe2+ (ferrous)", "Fe3+ (ferric)", "Siderophore biosynthesis", "Siderophore uptake", "troR", "Fur family")
type_free <- agg_test %>% filter(Filter_pores == "free-living") %>% filter(Type %in% wanted_types)
pathway_colors <- c("darkred", "darkgoldenrod", "#e66101", "darkorange", "#35978f", "#c7eae5", "seagreen", "dodgerblue3")
type_free$Type <- as.factor(type_free$Type) # setting the Order columns to factor
names(pathway_colors) <- levels(type_free$Type)
type_free <- aggregate(Abundance ~ Coastal_Current_Number * Type * Depth_Threshold, data = type_free, FUN = mean)

barplot_free <- ggplot(type_free, aes(x = factor(Coastal_Current_Number), y = Abundance, fill = Type)) + 
  facet_grid(~factor(Depth_Threshold, levels=c("Surface", "Intermediate"))~.) + # facet grid seperates by different levels, horizontally
  geom_bar(stat = "identity", position="fill", color= "black", linewidth=0.3, width=0.9) + theme_classic() + 
  # adds black outline to boxes
  scale_y_continuous(expand = c(0, 0)) + # extends the barplots to the axies
  scale_fill_manual(values = pathway_colors, drop = FALSE) + # set the colors with custom colors (myColors)
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=9, angle=90),
        axis.title.y = element_blank(),
        panel.spacing.x = unit(0, "points"), # Reducing space between facets
        strip.background = element_blank(), # Optionally hide the strip background for a cleaner look
        panel.border = element_blank()) + # Optionally remove panel borders
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # for the legend, if you want one
  #ylab("Relative Abundance (Order > 2%) \n") + # remove # if you want y-axis title
  ggtitle("Free-living (<0.2 µm)")

type_part <- agg_test %>% filter(Filter_pores == "particle-associated") %>% filter(Type %in% wanted_types)
barplot_part <- ggplot(type_part, aes(x = factor(Coastal_Current_Number), y = Abundance, fill = Type)) + 
  facet_grid(~factor(Depth_Threshold, levels=c("Surface", "Intermediate"))~.) + # facet grid seperates by different levels, horizontally
  geom_bar(stat = "identity", position="fill", color= "black", linewidth=0.3, width=0.9) + theme_classic() + # adds black outline to boxes
  scale_y_continuous(expand = c(0, 0)) + # extends the barplots to the axies
  scale_fill_manual(values = pathway_colors, drop = FALSE) + # set the colors with custom colors (myColors)
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=9, angle=90),
        axis.title.y = element_blank(),
        panel.spacing.x = unit(0, "points"), # Reducing space between facets
        strip.background = element_blank(), # Optionally hide the strip background for a cleaner look
        panel.border = element_blank()) + # Optionally remove panel borders
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + # for the legend, if you want one
  #ylab("Relative Abundance (Order > 2%) \n") + # remove # if you want y-axis title
  ggtitle("Particle-associated (>3 µm)")

coastal_combined <- ggarrange(
  barplot_free, barplot_part, labels = NULL,
  common.legend = TRUE, legend = "right"
)

annotate_figure(coastal_combined, top = text_grob("\n Coastal Current (Transect 3)", 
                                                  color = "black", face = "bold", size = 18))

ggsave("PICRUSt_analysis/graphics/coastal_type_barplot.pdf", width = 13, height = 7, dpi = 150)
