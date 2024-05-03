# libraries
library("dplyr")
library("tidyr")
library("phyloseq")
library("qiime2R")
library("ggplot2")
library("vegan")
library("plyr")
library("fantaxtic")
library("ggpubr")
library(tidyverse)

# list based on table !
phylum_to_keep=c("SAR324 clade(Marine group B)") #BACTEROIDOTA = Flavobacteriales order
class_to_keep=c("Zetaproteobacteria","Actinobacteria","Acidimicrobiia")
order_to_keep=c("Desulfuromonadales","Rhodobacterales","Alteromonadales","Nitrosomonadales",
                "Acidithiobacillales","SAR86 clade","Rhodospirillales","Acidiferrobacterales") # removed Flavobacteriales
family_to_keep=c("Acetobacteraceae","Alteromonadaceae", "Nitrosomonadaceae")
genus_to_keep=c("Caldithrix","Acidibater","SUP05 cluster","Acinetobacter") # replaced Pseudomonadales in order with acinetobacter
species_to_keep=c("Methylococcaceae bacterium SF-BR")

taxonomy <- as.data.frame(ps_noncontam_prev05@tax_table) #set taxonomy table from  as taxonomy (for next code)

tax_phylum <- taxonomy %>%
  filter(Phylum %in% phylum_to_keep) %>% 
  rownames_to_column("id")
Phylum <- tax_phylum$Phylum 
tax_phylum$TAX <- Phylum # make new column with tax ID
tax_phylum <- select(tax_phylum, 9, 1:8)

tax_class <- taxonomy %>% 
  filter(Class %in% class_to_keep) %>%
  rownames_to_column("id")
class <- tax_class$Class 
tax_class$TAX <- class # make new column with tax ID
tax_class <- select(tax_class, 9, 1:8)

tax_order <- taxonomy %>%
  filter(Order %in% order_to_keep) %>% 
  rownames_to_column("id")
order <- tax_order$Order 
tax_order$TAX <- order # make new column with tax ID
tax_order <- select(tax_order, 9, 1:8)

tax_family <- taxonomy %>% 
  filter(Family %in% family_to_keep) %>% 
  rownames_to_column("id")
family <- tax_family$Family
tax_family$TAX <- family # make new column with tax ID
tax_family <- select(tax_family, 9, 1:8)

tax_genus <- taxonomy %>% 
  filter(Genus %in% genus_to_keep) %>% 
  rownames_to_column("id")
genus <- tax_genus$Genus
tax_genus$TAX <- genus # make new column with tax ID
tax_genus <- select(tax_genus, 9, 1:8)

tax_species <- taxonomy %>% 
  filter(Species == "Methylococcaceae bacterium SF-BR") %>% 
  rownames_to_column("id")
species <- tax_species$Species
tax_species$TAX <- species # make new column with tax ID
tax_species <- select(tax_species, 9, 1:8)

merged_ps <- list(tax_phylum, tax_class, tax_order, tax_family, tax_genus, tax_species) # make list of the filtered taxonomy tables
merged_tax <- join_all(merged_ps, by = 'id', type = 'full') # merge taxonomy, deleting duplicates

list_filtered_tax <- merged_tax$id # create a list of tax id's for taxonomy
ASV_temp <- as.data.frame(ASV) %>% rownames_to_column("id") # create dataframe from ASV, then create rowname for id
new_ASV <- filter(ASV_temp, id %in% list_filtered_tax) # filter ASV based on list of filtered taxonomy
rownames(new_ASV) <- new_ASV$id # make taxa rownames
new_ASV$id <- NULL # remove column

rownames(merged_tax) <- merged_tax$id # make taxa rownames
merged_tax$id <- NULL # remove column

#### NEW PHYLOSEQ OBJECT
filtered_ASV <- otu_table(new_ASV, taxa_are_rows = TRUE)
filtered_taxmat <- as.matrix(merged_tax)
filtered_TAX <- tax_table(filtered_taxmat) #covert taxonomy table to phyloseq object

filtered_ps <- merge_phyloseq(filtered_ASV, filtered_TAX, META)