# libraries
library("dplyr")
library("tidyr")
library("phyloseq")
library("qiime2R")
library("ggplot2")
library("vegan")
library("fantaxtic")
library("ggpubr")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("decontam")
library("decontam")

## making relative abundance plots with distance as factor

##################
# phyloseq setup #
##################

# feature table
ASV <- qza_to_phyloseq(features="required_files/table.qza")
# read in metadata
metatable <- read.delim("/Users/sierra/Documents/Research/PICRUST_KEGG_analysis/artemis-eDNA-metadata-final.tsv", sep="\t", header=TRUE, row.names="sample_name") 
#metatable <- filter(metatable, Sample.Control == "True.Sample")# filter by transect
metatable$is.neg <- metatable$Sample.Control == "Control.Sample"
metatable$Final_Qubit <- as.numeric(metatable$Final_Qubit) 
#metatable <- filter(metatable, sample.illumina != "078_1040_PRE")
# importing (continued)
#row.names(metatable) <- metatable[["SampleID"]]
#metatable <- metatable %>% select(SampleID, everything())
META <- sample_data(metatable)

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
ps_contamdf_comb05_list <- rownames_to_column(ps_contamdf_comb05, var = "ASV")
#write.csv(ps_contamdf_comb05_list, "ps_contamdf_comb05_list.csv")

# Manual inspection, comparing blanks/controls and list of potential contaminants given by decontam, suggests prev05 is the most appropriated method
# Remove contaminants from prune phyloseq object using prev05.
ps_noncontam_prev05 <- prune_taxa(!ps_contamdf_comb05$contaminant, ps)
ps_noncontam_prev05 #  4559 taxa and 279 samples (originally was 4692 taxa)


