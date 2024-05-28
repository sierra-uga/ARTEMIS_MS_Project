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
library(microViz)

abundance_file <- "PICRUSt_analysis/required_files/pred_metagenome_unstrat_descrip.tsv" # read in abundance file (KO)

# test to make sure they sample names are the same
#list1 <- metadata$sample_name
#list2 <- colnames(abundance_data)
#list2[!(list1 %in% list2)] # output: [1] "STN115.35.fil.dura.r1"
metatable <- read.delim("required_files/artemis-eDNA-metadata-final.tsv", sep="\t", header=TRUE) 

metadata <- metatable %>% filter(., Sample.Control == "True.Sample") %>% filter(., Iron != "NA") %>%
  filter(., sample_name != "STN089.200.fil.dura.r2") #%>% group_by(Station) #%>% filter(., sample_name != "STN115.35.fil.dura.r1") #%>% distinct(Filter_pores, .keep_all = TRUE) 
# remove sample that isn't in kegg abundance for some reason

abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, col_types=c("c", "n"), trim_ws = TRUE) 
ColumnstoKeep <- c("description", metadata$sample_name) # set vector of list of names to keep + KO Name column
KO_abundance_data <- subset(abundance_data, select = ColumnstoKeep) # subset (select columns) based on ColumnstoKeep

#for iron

iron_KO <- read.csv("PICRUSt_analysis/required_files/KO_numbers_Sun_et_al.csv") # read in KO_number reference
iron_KO <- iron_KO %>% filter(., Metabolism == "Iron uptake and metabolism") # filter by iron metabolism only

KO_iron_numbers <- iron_KO$KO_Num # set vector for numbers
KO_iron_abundance_data <- KO_abundance_data[KO_abundance_data$`description` %in% KO_iron_numbers, ] #filter by KO_number using KO_iron_number ref

KO_joined <- data.frame(iron_KO$Name, iron_KO$KO_Num)
colnames(KO_joined) <- c("Name", "description")

KO_iron_abundance_type <- left_join(KO_iron_abundance_data, KO_joined) %>% column_to_rownames("description")

# create a phyloseq object for PICRUST data
metadata$Filter_pores <- ifelse(metadata$Filter_pores >= 0.2 & metadata$Filter_pores <= 2.0, "free-living", 
                          ifelse(metadata$Filter_pores == 3.0, "particle-associated", metadata$Filter_pores))
rownames(metadata) <- metadata$sample_name
META <- sample_data(metadata)

#for iron
rownames(KO_iron_abundance_type) <- KO_iron_abundance_type$Name
KO_iron_abundance_type$Name <- NULL
otumaty = as(KO_iron_abundance_type, "matrix")
#rownames(otumaty) <- abundance_data$"#NAME"
OTUy = otu_table(otumaty, taxa_are_rows=TRUE)

#for normal
KO_abundance_data <- as.data.frame(KO_abundance_data)
KO_abundance_data <- KO_abundance_data %>% distinct(description, .keep_all = TRUE)
rownames(KO_abundance_data) <- KO_abundance_data$description
KO_abundance_data$description <- NULL
otumaty = as(KO_abundance_data, "matrix")
#rownames(otumaty) <- abundance_data$"#NAME"
OTUy = otu_table(otumaty, taxa_are_rows=TRUE)

pi_ps <- merge_phyloseq(OTUy, META)

pseq <- pi_ps %>%
  phyloseq_validate()

pseq <- transform_sample_counts(pseq, function(OTU) OTU/sum(OTU) )
new_pseq = filter_taxa(pseq, function(x) sum(x) > 0.17, TRUE)


pseq %>%
  tax_transform(rank = "unique", trans = "rclr") %>%
  ord_calc(
    method = "auto"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Depth_Threshold", fill = "Depth_Threshold",
    shape = "circle", alpha = 0.5,
    size = 2
  ) +
  ggplot2::stat_ellipse(
    ggplot2::aes(colour = Depth_Threshold)
  )



ord_explore(pseq) 



htmp <- pseq %>%
  tax_transform(trans = "rclr", rank = "unique", zero_replace = "halfmin") %>%
  #tax_filter(min_prevalence = 0.3, use_counts = FALSE) %>%
  comp_heatmap(
    colors = heat_palette(sym = TRUE), grid_col = "black",
    tax_anno = taxAnnotation(
      Prev. = anno_tax_prev(bar_width = 0.3, size = grid::unit(1, "cm"))
    ),
    sample_side = "top", name = "Robust\nCLR",
    sample_anno = sampleAnnotation(
      "Filter_pores" = anno_sample("Filter_pores"),
      "Depth_Threshold" = anno_sample("Depth_Threshold"),
      col = list("Depth_Threshold" = c(
        "Bottom_water" = "black", "Intermediate" = "orange", "Surface" = "lightgrey"
      ))
    )
  )

dev.off()


htmp + htmp # to facet the thingie