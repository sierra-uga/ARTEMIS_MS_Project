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
outflow_metadata <- read_delim(
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
  filter(., sample_name != "STN089.200.fil.dura.r2") %>% filter(., Transect_Name == "transect2") %>% filter(., Station != "STN22") #%>% distinct(Filter_pores, .keep_all = TRUE) 
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

type_free <- agg_test %>% filter(Filter_pores == "free-living") %>% filter(Type %in% wanted_types)

barplot_free <- ggplot(type_free, aes(x = Transect_Number, y = Abundance, fill = Type)) + # facet grid seperates by different levels, horizontally
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
  ggtitle("Free-living (<0.2 µm)")

type_part <- agg_test %>% filter(Filter_pores == "particle-associated") %>% filter(Type %in% wanted_types)
total_count <- sum(type_part$Abundance)

# Calculate relative abundance
type_part$relative_abundance <- type_part$Abundance / total_count
barplot_part <- ggplot(type_part, aes(x = factor(Transect_Number), y = relative_abundance, color = Type, group=Type)) + # facet grid seperates by different levels, horizontally
  geom_line(linewidth=0.6) + theme_classic() + 
  geom_point() + # adds black outline to boxes
  scale_y_continuous(expand = c(0, 0)) + 
  scale_x_discrete(expand = c(0, .1)) +# extends the barplots to the axies
  scale_color_manual(values = pathway_colors, drop = TRUE) + # set the colors with custom colors (myColors)
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=9, angle=90),
        axis.title.y = element_blank(),
        panel.spacing.x = unit(0, "points"), # Reducing space between facets
        strip.background = element_blank(), # Optionally hide the strip background for a cleaner look
        panel.border = element_blank()) + # Optionally remove panel borders
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) + 
  geom_vline(xintercept = c(1,2,3), linetype = "dashed", linewidth=0.2, color = "lightgrey") +# for the legend, if you want one
  #ylab("Relative Abundance (Order > 2%) \n") + # remove # if you want y-axis title
  ggtitle("Particle-associated (>3 µm)")

ps_combined <- ggarrange(
  barplot_free, barplot_part, labels = NULL,
  common.legend = TRUE, legend = "right"
)

annotate_figure(ps_combined, top = text_grob("\n Outflow (Transect 2)", 
                                                  color = "black", face = "bold", size = 18))

ggsave("PICRUSt_analysis/graphics/outflow_type_bar_line_plot.pdf", width = 8, height = 7, dpi = 150)
