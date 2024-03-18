# libraries
library("dplyr")
library("tidyr")
library("phyloseq")
library("qiime2R")
library("ggplot2")
library("vegan")
library("fantaxtic")
library("ggpubr")

## making relative abundance plots with distance as factor

##################
# phyloseq setup #
##################

# feature table
ASV <- qza_to_phyloseq(features="table.qza")

# read in metadata
metatable <- read.delim("/Users/sierra/Documents/Research/PICRUST_KEGG_analysis/artemis-eDNA-metadata-final.tsv", sep="\t", header=TRUE, row.names="sample_name") 
metatable <- filter(metatable, Sample.Control == "True.Sample")# filter by transect

# importing (continued)
#row.names(metatable) <- metatable[["SampleID"]]
#metatable <- metatable %>% select(SampleID, everything())
META <- sample_data(metatable)

# importing taxonomy
taxonomy <- read.delim("taxonomy.tsv", sep="\t", header=TRUE) 
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
#genus <- taxonomy$Genus 
# family <- taxonomy$Family
taxonomy$TAX <- genus # make new column with 
taxonomy <- select(taxonomy, 8, 1:7) # moves new column to front
# OPTIONAL - makes actual row names the genus
#.rowNamesDF(taxonomy, make.names=TRUE) <- genus
names <- row.names(taxonomy)
row.names(ASV) <- names
# tax table cont.
taxmat <- as.matrix(taxonomy)
TAX <- tax_table(taxmat) #covert taxonomy table to phyloseq object

# import rooted tree
TREE <- qza_to_phyloseq(tree="rooted-tree.qza")
TREE[["tip.label"]] <- names

# merge all imported objects into phyloseq
ps <- merge_phyloseq(ASV, TAX, META, TREE)
ps
