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
library(phyloseq)
library(ComplexHeatmap)
library(DESeq2)
detach("package:microbiomeMarker", unload=TRUE)
abundance_file <- "PICRUSt_analysis/required_files/pred_metagenome_unstrat_descrip.tsv" # read in abundance file (KO)

# test to make sure they sample names are the same
#list1 <- metadata$sample_name
#list2 <- colnames(abundance_data)
#list2[!(list1 %in% list2)] # output: [1] "STN115.35.fil.dura.r1"
metatable <- read_delim("required_files/artemis-eDNA-metadata-final.tsv", delim="\t") 

metadata <- metatable %>% filter(., Sample.Control == "True.Sample") %>% #filter(., Iron != "NA") %>%
  filter(., sample_name != "STN089.200.fil.dura.r2") %>% filter(., sample_name != "STN078.1040.pre.poly.3.LG")#%>% group_by(Station) #%>% filter(., sample_name != "STN115.35.fil.dura.r1") #%>% distinct(Filter_pores, .keep_all = TRUE) 
# remove sample that isn't in kegg abundance for some reason

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
KO_iron_abundance_type$Name <- NULL
KO_iron_abundance_type$`function` <- NULL
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
TAX <- tax_table(iron_KO)
rownames(TAX) <- rownames(iron_KO)

TAX <- subset_taxa(TAX, rownames(OTUy) %in% TAX)

pi_ps <- merge_phyloseq(OTUy, META)

pseq <- pi_ps %>%
  phyloseq_validate() %>% subset_samples(., Location %in% c("Eastern_CC", "Western_CC")) %>% subset_samples(., More_Depth_Threshold %in% c("Surface", "Mid-Surface"))

ps_free <- subset_samples(pseq, Filter_pores == "free-living")
ps_part <- subset_samples(pseq, Filter_pores == "particle-associated")

#pseq <- transform_sample_counts(pseq, function(OTU) OTU/sum(OTU))
new_pseq = filter_taxa(pseq, function(x) sum(x) > 1, TRUE)
new_pseq <- transform_sample_counts(new_pseq, function(OTU) OTU * sum(OTU))

obj <- phyloseq_to_deseq2(ps_part, ~ More_Depth_Threshold + Location)
diagdds = DESeq(obj, test="Wald", fitType="mean", full= design(obj))

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(pseq)[rownames(sigtab), ], "matrix"))
head(sigtab)

KO_sig_numbers <- sigtab$KO_Num # set vector for numbers
KO_sig_tab <- iron_KO[iron_KO$KO_Num %in% KO_sig_numbers, ] 

iron_KO <- iron_KO[!duplicated(iron_KO$KO_Num), ]
row.names(iron_KO) <- iron_KO$KO_Num
sigtab$KO_Num <- sigtab$unique

sigtab_new <- dplyr::inner_join(sigtab, iron_KO, by="KO_Num")

ggplot(sigtab_new, aes(x=unique, y=log2FoldChange, color=Metabolism)) +
  geom_point(size=3, width = 0.2) + scale_x_discrete(labels = sigtab_new$Pathway.Final.Annotation.used.in.Fig.3.5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Eastern_CC (+) vs Western_CC (-)")
ggsave("relative_abundance_scripts/graphics/log_deseq_western_vs_eastern_particle.pdf", width=8, height=6)


KO_sig_numbers <- sigtab$KO_Num # set vector for numbers
KO_sig_tab <- iron_KO[iron_KO$KO_Num %in% KO_sig_numbers, ] 

ps_part_final <- subset_taxa(ps_part, unique %in% KO_sig_numbers)

colors = structure(1:31, names = unique(sigtab_new$Pathway.Final.Annotation.used.in.Fig.3.5))

colors2 = structure(1:31, names = unique(sigtab_new$Metabolism))
ha_right_ec2 = rowAnnotation(
  Pathway = sigtab_new$Pathway.Final.Annotation.used.in.Fig.3.5, border = TRUE, col  = list(Process = colors))

ha_right_ec1 = rowAnnotation(
  Nutrient = sigtab_new$Metabolism, border = TRUE, col  = list(Process = colors2))

row_labels_ec = sigtab_new$Pathway.Final.Annotation.used.in.Fig.3.5

htmp_part <- ps_part_final %>%
  tax_transform(trans = "rclr", rank = "unique", zero_replace = 1) %>%
  #tax_filter(min_prevalence = 0.1, use_counts = FALSE) %>%
  comp_heatmap(
    colors = heat_palette(sym = TRUE), grid_col = "black",
    tax_anno = taxAnnotation(
    Prev. = anno_tax_prev(bar_width = 0.3, size = grid::unit(1, "cm"))),
    sample_side = "top", name = "Robust\nCLR",
    sample_anno = sampleAnnotation(border = TRUE,
      "Location" = anno_sample("Location"),
      "Depth_Threshold" = anno_sample("Depth_Threshold"),
      col = list("Depth_Threshold" = c(
        "Bottom_water" = "black", "Intermediate" = "orange", "Surface" = "lightgrey"
      ))
    )
  ) + ha_right_ec2 + ha_right_ec1 #+ row_labels_ec

head(ps_part@otu_table)

ha <- rowAnnotation(
  Iron = sigtab_new$Metabolism, border = TRUE, show_legend=TRUE)

draw(ha)

lgd_sig = Legend(pch = "*", type = "points", labels = "< 0.01")
htmp_part %>% ComplexHeatmap::htmp_part + draw(ha)
)

htmp_part + htmp_free




dev.off()