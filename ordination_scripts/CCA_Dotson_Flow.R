# load in libraries
library("dplyr")
library("tidyr")
library("phyloseq")
library("qiime2R")
library("ggplot2")
# devtools::install_github("gmteunisse/fantaxtic")
library("fantaxtic")
library("vegan")
library("ggpubr")

##################
# phloseq setup  #
##################

# feature table
ASV <- qza_to_phyloseq(features="table.qza")

# read in metadata
metatable <- read.delim("artemis-eDNA-metadata-final.tsv", sep="\t", header=TRUE, row.names="Sample.illumina") 
# fixing metadata
metatable <- metatable[-159,] # removes row 89_200_FIL_R2 because taxasum was 0.
metatable <- filter(metatable, Sample.Control == "True.Sample") %>% 
  mutate_at(c(20:27, 29:38, 44), as.numeric)  # remove blanks, convert character columns for env. data to numeric
metatable <- metatable[(which(metatable$Station %in% c("STN056a", "STN056b", "STN22", "STN078", "STN014", "STN106"))),] # select stations for dotson review

# importing (continued)
row.names(metatable) <- metatable[["SampleID"]]
metatable <- metatable %>% select(SampleID, everything())
META <- sample_data(metatable)
list_true <- replace_na(META$True_Flow, "Other") #replace NA with "Other" for coloring
META$True_Flow <- list_true # changing actual column in dataframe

# importing taxonomy
taxonomy <- read.delim("taxonomy.tsv", sep="\t", header=TRUE) 
names(taxonomy) <- c("row", "tax", "Confidence") # rename columns
row.names(taxonomy) <- taxonomy[[1]] # making the row names the tax name
taxonomy <- taxonomy[,(-1)] # removing the duplicate column
# SILVA taxonomy is in one column, separate to be able to work with different taxonomic levels:
taxonomy <-  separate(taxonomy, tax, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", 
                                       "D7", "D8", "D9", "D10", "D11", 
                                       "D12", "D13", "D14"), sep = ";", fill = "right")
taxonomy <- taxonomy[,c(1:7)] # select only first 1-7 columns (to genus)
nm1 <- colnames(taxonomy)
taxonomy[nm1] <- lapply(taxonomy[nm1], gsub, pattern = "D_.__", replacement = "")
# turn the otu_table into a data.frame
taxonomy <- taxonomy[-1, ] # removes first row of taxonomy table
# change tax table
genus <- taxonomy$Genus 
taxonomy$TAX <- genus # make new column with 
taxonomy <- select(taxonomy, 8, 1:7) # moves new column to front
# OPTIONAL - makes actual row names the genus
.rowNamesDF(taxonomy, make.names=TRUE) <- genus
names <- row.names(taxonomy)
row.names(ASV) <- names
taxmat <- as.matrix(taxonomy)
TAX <- tax_table(taxmat) #covert taxonomy table to phyloseq object

# import rooted tree
TREE <- qza_to_phyloseq(tree="rooted-tree.qza")
TREE[["tip.label"]] <- names

# merge all imported objects into phyloseq
ps <- merge_phyloseq(ASV, TAX, META, TREE)
ps

##################
#  data culling  #
##################

# select only bacteria, remove chloroplasts
ps_sub <- ps %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "Mitochondria" &
      Class   != "Chloroplast" &
      Order   != "Chloroplast" &
      Family  != "Mitochondria" 
  )

ps_sub <- ps_sub %>% prune_taxa(taxa_sums(.) > 0, .) # remove 0 taxasums

# free-living phyloseq
ps_free <- ps_sub %>% subset_samples(Filter_pores == "0.2")

# particle-associated phyloseq
ps_part <- ps_sub %>% subset_samples(Filter_pores >= "2")

# selecting the top free-living taxa
top_free <- top_taxa(ps_free, 
                     tax_level = "Genus", 
                     n_taxa = 10)

# selecting the top particle-associated taxa
top_part <- top_taxa(ps_part, 
                     tax_level = "Genus", 
                     n_taxa = 10)

## aesthetics
# set ggplot shortcut to remove grid
remove_grid <- theme(legend.position = "bottom") + theme_bw() + # removes grid
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

# color for ggplot points
color_breaks <- unique(sample_data(ps_sub)$True_Flow)
color_point <- scale_color_manual(values = c("gray", "dodgerblue", "red2"),
                                  name = "True_Flow",
                                  breaks = color_breaks,
                                  labels = color_breaks)

##################
# CCA ordination #
##################

# ordinate actual_CCA_ref
cca1 <- ordinate(ps_free, method = "CCA")
cca2 <- ordinate(ps_part, method = "CCA")
# ordinate top_CCA_ref
top_cca1 <- ordinate(top_free$ps_obj, "CCA")
top_cca2 <- ordinate(top_part$ps_obj, "CCA")

# Free-living CCA 
free_CCA <- plot_ordination(ps_free, cca1,
                            type = "samples", color = "True_Flow")

# Particle-associated CCA
part_CCA <- plot_ordination(ps_part, cca2,
                            type = "samples", color = "True_Flow")

arrowmat_free <- scores(top_cca1, display = "species")
# Add labels, make a data.frame
arrowdf_1 <- data.frame(labels = rownames(arrowmat_free), arrowmat_free)
arrowmat_part <- scores(top_cca2, display = "species")
# Add labels, make a data.frame
arrowdf_2 <- data.frame(labels = rownames(arrowmat_part), arrowmat_part)

free_living_flow <- CCA_ord(top_cca1, free_CCA, "Top 10 (Genus)", free_living_flow)
part_living_flow <- CCA_ord(top_cca2, part_CCA, "Top 10 (Genus)", part_living_flow)

ggarrange(
  free_living_flow, part_living_flow, labels = c("F", "P"),
  common.legend = TRUE, legend = "right"
)
ggsave("graphics/split_Dotson_CCA.pdf", width = 11, height = 8, dpi = 150)



