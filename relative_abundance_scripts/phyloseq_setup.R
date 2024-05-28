# libraries
library("dplyr")
library("tidyr")
library("phyloseq")
install.packages("qiime2R")
library("qiime2R")
library("ggplot2")
library("vegan")
library("fantaxtic")
install.packages("devtools")
library("ggpubr")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages("BioConductor")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("decontam")
BiocManager::install("TreeSummarizedExperiment")

#BiocManager::install("decontam")
remotes::install_github("mikemc/speedyseq")
library("speedyseq")
library("decontam")
install.packages("oce")
install.packages("phyloseq")

if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")

if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")
## making relative abundance plots with distance as factor

##################
# phyloseq setup #
##################
setwd("~/Documents/Research/Ordination analysis R scripts/ARTEMIS_github")
# feature table
ASV <- qza_to_phyloseq(features="required_files/table.qza")
# read in metadata
metatable <- read.delim("required_files/artemis-eDNA-metadata-final.tsv", sep="\t", header=TRUE) 
#metatable <- filter(metatable, Sample.Control == "True.Sample")# filter by transect
metatable$is.neg <- metatable$Sample.Control == "Control.Sample"
metatable$Final_Qubit <- as.numeric(metatable$Final_Qubit) 
df <- metatable %>%
  group_by(Station, More_Depth_Threshold) %>%
  mutate(unique_code = paste0(Station, More_Depth_Threshold))# creates a column with a unique code for watertype + station
df <- df %>%
  group_by(Station, Depth_Threshold) %>%
mutate(unique_depth = paste0(Station, Depth_Threshold))

df <- df %>% mutate_at(c(21:26, 31:40), as.numeric) 
df <- as.data.frame(df) 

# Define a function to extract the replicate number
extract_replicate <- function(sample_name) {
  parts <- strsplit(sample_name, "\\.")[[1]]
  last_part <- tail(parts, 1)
  
  if (grepl("LG$", last_part)) {
    return("NA")  # Special case for "LG"
  } else {
    num <- gsub("[^0-9]", "", last_part)  # Extract numeric part
    return(num)
  }
}

# Apply the function to create the Replicate_Number column
df <- df %>%
  mutate(Replicate_Number = sapply(sample_name, extract_replicate))

rownames(df) <- df$sample_name
META <- sample_data(df)

# importing taxonomy
taxonomy <- read.delim("required_files/taxonomy.tsv", sep="\t", header=TRUE) 
names(taxonomy) <- c("row", "tax", "Confidence") # rename columns # CONFIDENCE LEVEL .7 in the QIIME2 script
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
#Order <- taxonomy$Order 
# family <- taxonomy$Family
#taxonomy$TAX <- Order # make new column with 
#taxonomy <- select(taxonomy, 8, 1:7) # moves new column to front
# OPTIONAL - makes actual row names the genus
#.rowNamesDF(taxonomy, make.names=TRUE) <- genus
names <- row.names(taxonomy)
row.names(ASV) <- names
# tax table cont.
taxmat <- as.matrix(taxonomy)
TAX <- tax_table(taxmat) #covert taxonomy table to phyloseq object

# import rooted tree
TREE <- qza_to_phyloseq(tree="required_files/rooted-tree.qza")
TREE[["tip.label"]] <- names

# merge all imported objects into phyloseq
ps <- merge_phyloseq(ASV, TAX, META, TREE)
ps <- subset_samples(ps, sample.illumina != "078_1040_PRE") # removed bc qubit value is 0 

# Combined method: here we increase the threshold to 0.5
ps_contamdf_comb05 <- isContaminant(ps, conc="Final_Qubit", neg="is.neg", threshold=0.5, detailed = TRUE, normalize = TRUE, method="combined")
table(ps_contamdf_comb05$contaminant) # it identified 564 potential contaminants
head(which(ps_contamdf_comb05$contaminant))
#ps_contamdf_comb05_list <- rownames_to_column(ps_contamdf_comb05, var = "ASV")
#write.csv(ps_contamdf_comb05_list, "ps_contamdf_comb05_list.csv")

# Manual inspection, comparing blanks/controls and list of potential contaminants given by decontam, suggests prev05 is the most appropriated method
# Remove contaminants from prune phyloseq object using prev05.
ps_noncontam_prev05 <- prune_taxa(!ps_contamdf_comb05$contaminant, ps)
ps_noncontam_prev05 #  4559 taxa and 279 samples (originally was 4692 taxa)

####
filtered_df <- as.data.frame(dist_bray_mtx)
threshold <- 0.2
# Use dplyr's mutate and ifelse to replace values
df_filtered <- filtered_df %>%
  mutate_all(~if_else(. >= threshold, NA_real_, .))
