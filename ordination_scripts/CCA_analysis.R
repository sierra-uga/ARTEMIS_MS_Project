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
# metatable <- metatable[-159,] # removes row 89_200_FIL_R2 because taxasum was 0.
metatable <- filter(metatable, Sample.Control == "True.Sample") %>% 
  mutate_at(c(20:27, 29:38, 44), as.numeric) # remove blanks, convert character columns for env. data to numeric
metatable <- metatable[-(which(metatable$Station %in% c("STN153", "STN174", "STN089", "STN151.2", "STN198", "STN002", "STN004"))),] 


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

Samples_toRemove <- c("STN089.200.fil.dura.r2", "STN078.1040.pre.poly.3.LG")  # remove from phyloseq
ps <- subset_samples(ps, !(SampleID %in% Samples_toRemove)) # remove from phyloseq

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
ps_free <- ps_sub %>% subset_samples(Filter_pores == "0.2") %>% prune_taxa(taxa_sums(.) > 0, .) 

# particle-associated phyloseq
ps_part <- ps_sub %>% subset_samples(Filter_pores >= "2") %>% prune_taxa(taxa_sums(.) > 0, .) 

# top taxa from all taxa
top_all <- top_taxa(ps_sub, 
                     tax_level = "Genus", 
                     n_taxa = 10)

# selecting the top free-living taxa
top_free <- top_taxa(ps_free, 
                tax_level = "Genus", 
                n_taxa = 10)

# selecting the top particle-associated taxa
top_part <- top_taxa(ps_part, 
                tax_level = "Genus", 
                n_taxa = 10)

top_taxa_free <- as.data.frame(top_free[["top_taxa"]])
top_taxa_part <- as.data.frame(top_part[["top_taxa"]])

##################
# CCA ordination #
##################

# Ordinate 2
cca <- ordinate(ps_sub, method = "CCA", formula = ~ watertype)
cca1 <- ordinate(ps_free, method = "CCA", formula = ~ watertype)
cca2 <- ordinate(ps_part, method = "CCA", formula = ~ watertype)
# Ordinate 3
top_cca <- ordinate(top_all$ps_obj, "CCA", formula= ~ watertype)
top_cca1 <- ordinate(top_free$ps_obj, "CCA", formula= ~ watertype)
top_cca2 <- ordinate(top_part$ps_obj, "CCA", formula= ~ watertype)

## aesthetics
# set ggplot shortcut to remove grid
remove_grid <- theme(legend.position = "bottom") + theme_bw() + # removes grid
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              axis.line = element_line(colour = "black"))

# color for ggplot points
color_breaks <- unique(sample_data(ps_sub)$watertype)
color_point <- scale_color_manual(values = c("darkgreen", "dodgerblue", "red2", "blueviolet", "aquamarine3", "gray"),
                   name = "Water Mass",
                   breaks = color_breaks,
                   labels = color_breaks)

# Both free-living + particle-associated CCA
all_CCA <- plot_ordination(ps_sub, cca,
                                   type = "samples", color = "watertype")
all_CCA <- all_CCA + geom_point(size = 2) + remove_grid + color_point
all_CCA
ggsave("graphics/all_CCA.pdf", width = 7, height = 6, dpi = 150)

# Free-living CCA 
free_CCA <- plot_ordination(ps_free, cca1,
                     type = "samples", color = "watertype")
f_CCA <- free_CCA + geom_point(size = 2) + remove_grid + color_point
f_CCA

# Particle-associated CCA
part_CCA <- plot_ordination(ps_part, cca2,
                                   type = "samples", color = "watertype")
p_CCA <- part_CCA + geom_point(size = 2) + remove_grid + color_point
p_CCA

# split plot for both communities
ggarrange(
  f_CCA, p_CCA, labels = c("F", "P"),
  common.legend = TRUE, legend = "right"
)
ggsave("graphics/split_CCA.pdf", width = 7, height = 6, dpi = 150)

# Mimic biplot top 40 taxa (split plot)
#p0 = plot_ordination(top_all$ps_obj, top_cca, type = "split", color = "watertype")
#p0 <- p0 + geom_point(size = 2)
#print(p0)

#######################
#  Using function to  #
#   make split graph  #
#######################

# use CCA_ord function (see files) to make graphs
final_free_plot <- CCA_ord(top_cca1, free_CCA, "CCA Top 10 Taxa (Genus)", final_free_plot)
final_part_plot <- CCA_ord(top_cca2, part_CCA, "CCA Top 10 Taxa (Genus)", final_part_plot)

# split plot for both communities
ggarrange(
  final_free_plot, final_part_plot, labels = c("F", "P"),
  common.legend = TRUE, legend = "right"
)
ggsave("graphics/split_CCA_genus.pdf", width = 11, height = 8, dpi = 150)


arrowmat1 <- scores(top_cca1, display = "species")
arrowmat2 <- scores(top_cca2, display = "species")
arrowdf1 <- data.frame(labels = rownames(arrowmat1), arrowmat1)

# Define the arrow aesthetic mapping
arrow_map1 = aes(xend = CCA1, yend = CCA2, x = 0, y = 0, shape = NULL, color = NULL, 
                label = labels)

# Add labels, make a data.frame
arrowdf_free <- data.frame(labels = rownames(arrowmat1), arrowmat1)

arrowdf_part <- data.frame(labels = rownames(arrowmat2), arrowmat2)


