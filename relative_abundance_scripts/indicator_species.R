install.packages("indicspecies")
library(indicspecies)
library(phyloseq)
library(fantaxtic)
library(dplyr)
library(data.table)
library(microViz)

ps_sub <- ps_noncontam_prev05 %>%
  subset_taxa(
    #Kingdom == "Bacteria" &
      Family  != "Mitochondria" &
      Family   != "Chloroplast" &
      Order   != "Chloroplast" &
      Family  != "Mitochondria" 
  )

ps_sub <- subset_samples(ps_sub, Sample.Control == "True.Sample") %>% subset_samples(., sample.illumina != "089_200_FIL_R2") %>%
  phyloseq_validate() %>% tax_fix() %>% prune_taxa(taxa_sums(.) > 0, .) 

ps_sub <- tax_glom(ps_sub, "Genus", NArm = TRUE)

ps_sub <-label_duplicate_taxa(ps_sub,
                              tax_level = "Genus")

taxa_names <- as.data.frame(otu_table(ps_sub))
taxa_names <- taxa_names[rowSums(taxa_names) !=0,
                         c(which(colSums(taxa_names) !=0))] # REMOVES 0s for CCA
taxa_names <- na.omit(taxa_names) # remove NAs just in case
otu_table(ps_sub) <- taxa_names # re-inserts the OTU table for the phyloseq object


ps_free <- ps_sub %>% subset_samples(Filter_pores == "free-living") %>% prune_taxa(taxa_sums(.) > 0, .) 

ps_free <- ps_free %>%
  tax_glom(taxrank = "Genus") %>% # agglomerate at Order level, can change to different taxonomic level!# %>%
  transform_sample_counts(function(x) {x/sum(x)})  # Transform to rel. abundance (normalize data

ps_free@sam_data[["CTD_Depth"]] <- as.numeric(ps_free@sam_data[["CTD_Depth"]])

## upper 200m
ps_free_above <- ps_free %>% subset_samples(CTD_Depth <= 200)

## bellow 200m
ps_free_below <- ps_free %>% subset_samples(CTD_Depth >= 200)

#### FREE-LIVING #####

write.csv(ps_free_above@otu_table,'M3.SP.otus_free_above.csv')
taxa_free_above <- as.data.frame(ps_free_above@tax_table)
taxa_free_above$taxon <- rownames(taxa_free_above)
#Import phyloseq OTU table as an OTU table/dataframe
SpOTU_free_above <-read.csv('M3.SP.otus_free_above.csv')

#do some shuffling of the OTU table
SpOTUFlip_free_above <- as.data.frame(t(SpOTU_free_above)) #makes it a dataframe and puts x into y and y into x (flips it)
names(SpOTUFlip_free_above) <- as.matrix(SpOTUFlip_free_above[1, ]) # renames columns
SpOTUFlip_free_above<- SpOTUFlip_free_above[-1, ] #removes first row
SpOTUFlip_num_free_above<-as.data.frame(lapply(SpOTUFlip_free_above, as.numeric)) #convert from character to number
SpOTUFlip_num_free_above$sample_name<-row.names(SpOTUFlip_free_above) #puts row names as sample ID column
#OK now we have the OTU table that's somewhat in the way they like

#read in metadata
metadata_free_above <-ps_free_above@sam_data #read in metadata
head(metadata_free_above) # check

## Join based on SampleID
SpOTU_Final_free_above <-left_join(SpOTUFlip_num_free_above, metadata_free_above, by = c("sample_name" = "sample_name")) # join based on sample IDs, assuming they're the same for both OTU table and metadata
SPotus_free_above = SpOTU_Final_free_above[,1:325] #select just the ASV/OTU table part of the file (you may have to scroll to the back of the OTU file to find it...)


SPwat_free_above = SpOTU_Final_free_above$Location #the metadata column group you care about
SPotus_free_above <- SPotus_free_above[, colSums(SPotus_free_above) != 0]
indisp_free_above=multipatt(x=SPotus_free_above, cluster=SPwat_free_above, func = "r.g", print.perm = TRUE, control = how(nperm=9999))

#extract table of stats
indisp.sign_free_above <-as.data.table(indisp_free_above$sign, keep.rownames=TRUE)
#add adjusted p-value
indisp.sign_free_above[ ,p.value.bh:=p.adjust(p.value, method="BH")]
#now can select only the indicators with adjusted significant p-values
SPind_free_above <- indisp.sign_free_above[p.value.bh<=0.05, ]

SPind_free_above$taxon <- SPind_free_above$rn

merged_df_free_above <- merge(SPind_free_above, taxa_free_above, by = "taxon")

write.csv(merged_df_free_above,'Location_with_cont_open_free_above.csv')

######### BELOW #####

write.csv(ps_free_below@otu_table,'M3.SP.otus_free_below.csv')
taxa_free_below <- as.data.frame(ps_free_below@tax_table)
taxa_free_below$taxon <- rownames(taxa_free_below)
#Import phyloseq OTU table as an OTU table/dataframe
SpOTU_free_below <-read.csv('M3.SP.otus_free_below.csv')

#do some shuffling of the OTU table
SpOTUFlip_free_below <- as.data.frame(t(SpOTU_free_below)) #makes it a dataframe and puts x into y and y into x (flips it)
names(SpOTUFlip_free_below) <- as.matrix(SpOTUFlip_free_below[1, ]) # renames columns
SpOTUFlip_free_below<- SpOTUFlip_free_below[-1, ] #removes first row
SpOTUFlip_num_free_below<-as.data.frame(lapply(SpOTUFlip_free_below, as.numeric)) #convert from character to number
SpOTUFlip_num_free_below$sample_name<-row.names(SpOTUFlip_free_below) #puts row names as sample ID column
#OK now we have the OTU table that's somewhat in the way they like

#read in metadata
metadata_free_below <-ps_free_below@sam_data #read in metadata
head(metadata_free_below) # check

## Join based on SampleID
SpOTU_Final_free_below <-left_join(SpOTUFlip_num_free_below, metadata_free_below, by = c("sample_name" = "sample_name")) # join based on sample IDs, assuming they're the same for both OTU table and metadata

SPotus_free_below = SpOTU_Final_free_below[,1:325] #select just the ASV/OTU table part of the file (you may have to scroll to the back of the OTU file to find it...)
SPwat_free_below = SpOTU_Final_free_below$Location #the metadata column group you care about
SPotus_free_below <- SPotus_free_below[, colSums(SPotus_free_below) != 0]

indisp_free_below=multipatt(x=SPotus_free_below, cluster=SPwat_free_below, func = "r.g", print.perm = TRUE, control = how(nperm=9999))

#extract table of stats
indisp.sign_free_below <-as.data.table(indisp_free_below$sign, keep.rownames=TRUE)
#add adjusted p-value
indisp.sign_free_below[ ,p.value.bh:=p.adjust(p.value, method="BH")]
#now can select only the indicators with adjusted significant p-values
SPind_free_below <- indisp.sign_free_below[p.value.bh<=0.05, ]

SPind_free_below$taxon <- SPind_free_below$rn

merged_df_free_below <- merge(SPind_free_below, taxa_free_below, by = "taxon")

write.csv(merged_df_free_below,'indicator_location_free_below.csv')

#### PARTICLE-ASSOCIATED
# particle-associated phyloseq
ps_part <- ps_sub %>% subset_samples(Filter_pores == "particle-associated") %>% prune_taxa(taxa_sums(.) > 0, .)  #subset_samples(Location %in% c("Dotson", "Open_polynya")) %>%

ps_part <- ps_part %>%
  tax_glom(taxrank = "Genus") %>% # agglomerate at Order level, can change to different taxonomic level!# %>%
  transform_sample_counts(function(x) {x/sum(x)})

ps_part_above <- ps_part %>% subset_samples(CTD_Depth <= 200)

## bellow 200m
ps_part_below <- ps_part %>% subset_samples(CTD_Depth >= 200)

write.csv(ps_part_above@otu_table,'M3.SP.otus_part_above.csv')
taxa_part_above <- as.data.frame(ps_part_above@tax_table)
taxa_part_above$taxon <- rownames(taxa_part_above)
#Import phyloseq OTU table as an OTU table/dataframe
SpOTU_part_above <-read.csv('M3.SP.otus_part_above.csv')

#do some shuffling of the OTU table
SpOTUFlip_part_above <- as.data.frame(t(SpOTU_part_above)) #makes it a dataframe and puts x into y and y into x (flips it)
names(SpOTUFlip_part_above) <- as.matrix(SpOTUFlip_part_above[1, ]) # renames columns
SpOTUFlip_part_above<- SpOTUFlip_part_above[-1, ] #removes first row
SpOTUFlip_num_part_above<-as.data.frame(lapply(SpOTUFlip_part_above, as.numeric)) #convert from character to number
SpOTUFlip_num_part_above$sample_name<-row.names(SpOTUFlip_part_above) #puts row names as sample ID column
#OK now we have the OTU table that's somewhat in the way they like

#read in metadata
metadata_part_above <-ps_part_above@sam_data #read in metadata
head(metadata_part_above) # check

## Join based on SampleID
SpOTU_Final_part_above <-left_join(SpOTUFlip_num_part_above, metadata_part_above, by = c("sample_name" = "sample_name")) # join based on sample IDs, assuming they're the same for both OTU table and metadata

SPotus_part_above = SpOTU_Final_part_above[,1:449] #select just the ASV/OTU table part of the file (you may have to scroll to the back of the OTU file to find it...)
SPwat_part_above = SpOTU_Final_part_above$Location #the metadata column group you care about
SPotus_part_above <- SPotus_part_above[, colSums(SPotus_part_above) != 0]

indisp_part_above=multipatt(x=SPotus_part_above, cluster=SPwat_part_above, func = "r.g", print.perm = TRUE, control = how(nperm=9999))

#extract table of stats
indisp.sign_part_above <-as.data.table(indisp_part_above$sign, keep.rownames=TRUE)
#add adjusted p-value
indisp.sign_part_above[ ,p.value.bh:=p.adjust(p.value, method="BH")]
#now can select only the indicators with adjusted significant p-values
SPind_part_above <- indisp.sign_part_above[p.value.bh<=0.05, ]

SPind_part_above$taxon <- SPind_part_above$rn

merged_df_part_above <- merge(SPind_part_above, taxa_part_above, by = "taxon")

write.csv(merged_df_part_above,'Location_no_cont_open_part_above.csv')

######### BELOW #####

write.csv(ps_part_below@otu_table,'M3.SP.otus_part_below.csv')
taxa_part_below <- as.data.frame(ps_part_below@tax_table)
taxa_part_below$taxon <- rownames(taxa_part_below)
#Import phyloseq OTU table as an OTU table/dataframe
SpOTU_part_below <-read.csv('M3.SP.otus_part_below.csv')

#do some shuffling of the OTU table
SpOTUFlip_part_below <- as.data.frame(t(SpOTU_part_below)) #makes it a dataframe and puts x into y and y into x (flips it)
names(SpOTUFlip_part_below) <- as.matrix(SpOTUFlip_part_below[1, ]) # renames columns
SpOTUFlip_part_below<- SpOTUFlip_part_below[-1, ] #removes first row
SpOTUFlip_num_part_below<-as.data.frame(lapply(SpOTUFlip_part_below, as.numeric)) #convert from character to number
SpOTUFlip_num_part_below$sample_name<-row.names(SpOTUFlip_part_below) #puts row names as sample ID column
#OK now we have the OTU table that's somewhat in the way they like

#read in metadata
metadata_part_below <-ps_part_below@sam_data #read in metadata
head(metadata_part_below) # check

## Join based on SampleID
SpOTU_Final_part_below <-left_join(SpOTUFlip_num_part_below, metadata_part_below, by = c("sample_name" = "sample_name")) # join based on sample IDs, assuming they're the same for both OTU table and metadata

SPotus_part_below = SpOTU_Final_part_below[,1:449] #select just the ASV/OTU table part of the file (you may have to scroll to the back of the OTU file to find it...)
SPwat_part_below = SpOTU_Final_part_below$Location #the metadata column group you care about
SPotus_part_below <- SPotus_part_below[, colSums(SPotus_part_below) != 0]

indisp_part_below=multipatt(x=SPotus_part_below, cluster=SPwat_part_below, func = "r.g", print.perm = TRUE, control = how(nperm=9999))

#extract table of stats
indisp.sign_part_below <-as.data.table(indisp_part_below$sign, keep.rownames=TRUE)
#add adjusted p-value
indisp.sign_part_below[ ,p.value.bh:=p.adjust(p.value, method="BH")]
#now can select only the indicators with adjusted significant p-values
SPind_part_below <- indisp.sign_part_below[p.value.bh<=0.05, ]

SPind_part_below$taxon <- SPind_part_below$rn

merged_df_part_below <- merge(SPind_part_below, taxa_part_below, by = "taxon")

write.csv(merged_df_part_below,'indicator_location_part_below.csv')

### ADONIS
# FREE-LIVING
### adonis 
library(vegan)
bray_free_above <- phyloseq::distance(ps_free_above, method = "bray") # setting distance
sampledf_free_above <- data.frame(sample_data(ps_free_above))# make a data frame from the sample_data

#select from main data frame
adonis_frame_free_above <- dplyr::select(sampledf_free_above, Station, Latitude:CTD_Depth, Lab_NO3:Iron, Iron_Level:Pos_in_polynya)
adonis_frame_free_above$Location <- as.factor(adonis_frame_free_above$Location)

# Adonis test
adonis_free_above <- adonis2(bray_free_above ~ Location, data = adonis_frame_free_above, permutations = 999)
pairwise.adonis(bray_free_above ,factors=adonis_frame_free_above$Location)

# Post hoc for Location in polynya
beta_location_free_above <- betadisper(bray_free_above, adonis_frame_free_above$Location)

layout_matrix <- matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE)
layout(layout_matrix)

# Adjust margins
par(mar = c(5, 6, 4, 2) + 0.1)
permutest(beta_location_free_above)
boxplot(beta_location_free_above)
mod.HSD_free_above <- TukeyHSD(beta_location_free_above)
plot(mod.HSD_free_above, las=1)
plot(beta_location_free_above)

### BELOW
bray_free_below <- phyloseq::distance(ps_free_below, method = "bray") # setting distance
sampledf_free_below <- data.frame(sample_data(ps_free_below))# make a data frame from the sample_data

#select from main data frame
adonis_frame_free_below <- dplyr::select(sampledf_free_below, Station, Latitude:CTD_Depth, Lab_NO3:Iron, Iron_Level:Pos_in_polynya)
adonis_frame_free_below$Location <- as.factor(adonis_frame_free_below$Location)

# Adonis test
adonis_free_below <- adonis2(bray_free_below ~ Location, data = adonis_frame_free_below, permutations = 999)
pairwise.adonis(bray_free_below ,factors=adonis_frame_free_below$Location)

# Post hoc for location in polynya
beta_location_free_below <- betadisper(bray_free_below, adonis_frame_free_below$Location)

layout_matrix <- matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE)
layout(layout_matrix)

# Adjust margins
par(mar = c(5, 6, 4, 2) + 0.1)
permutest(beta_location_free_below)
boxplot(beta_location_free_below)
mod.HSD_free_below <- TukeyHSD(beta_location_free_below)
plot(mod.HSD_free_below, las=1)
plot(beta_location_free_below)

----
  
bray_part_above <- phyloseq::distance(ps_part_above, method = "bray") # setting distance
sampledf_part_above <- data.frame(sample_data(ps_part_above))# make a data frame from the sample_data

#select from main data frame
adonis_frame_part_above <- dplyr::select(sampledf_part_above, Station, Latitude:CTD_Depth, Lab_NO3:Iron, Iron_Level:Pos_in_polynya)
adonis_frame_part_above$Location <- as.factor(adonis_frame_part_above$Location)

# Adonis test
adonis_part_above <- adonis2(bray_part_above ~ Location, data = adonis_frame_part_above, permutations = 9999)
pairwise.adonis(bray_part_above ,factors=adonis_frame_part_above$Location)
# Post hoc for location in polynya
beta_location_part_above <- betadisper(bray_part_above, adonis_frame_part_above$Location)

layout_matrix <- matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE)
layout(layout_matrix)

# Adjust margins
par(mar = c(5, 6, 4, 2) + 0.1)
permutest(beta_location_part_above)
boxplot(beta_location_part_above)
mod.HSD_part_above <- TukeyHSD(beta_location_part_above)
plot(mod.HSD_part_above, las=1)
plot(beta_location_part_above)

### BELOW

adonis_part_below <- adonis2(bray_part_below ~ Location, data = adonis_frame_part_below, permutations = 999)
pairwise.adonis(bray_part_below ,factors=adonis_frame_part_below$Location)

TukeyHSD(adonis_part_below, "city")

#select from main data frame
adonis_frame_part_below <- dplyr::select(sampledf_part_below, Station, Latitude:CTD_Depth, Lab_NO3:Iron, Iron_Level:Pos_in_polynya)
adonis_frame_part_below$Location <- as.factor(adonis_frame_part_below$Location)

# Adonis test
adonis_part_below <- adonis2(bray_part_below ~ Location, data = adonis_frame_part_below)

# Post hoc for location in polynya
beta_location_part_below <- betadisper(bray_part_below, adonis_frame_part_below$Location)


# Adjust margins

pdf(file = "final_graphics/adonis_location.pdf", width = 12, height = 11) 
layout_matrix <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)
layout(layout_matrix)
mod.HSD_free_above <- TukeyHSD(beta_location_free_above)
mod.HSD_free_below <- TukeyHSD(beta_location_free_below)
mod.HSD_part_above <- TukeyHSD(beta_location_part_above)
mod.HSD_part_below <- TukeyHSD(beta_location_part_below)
par(mar = c(2, 8, 2, 1) + 0.5)
plot(mod.HSD_free_above, las=1, cex.axis = 0.8, tcl = -0.2)
par(mar = c(2, 2, 2, 5) + 0.5)
plot(mod.HSD_part_above, las=1, yaxt='n', cex.axis = 0.8, tcl = -0.2)
par(mar = c(2, 8, 2, 1) + 0.5)
plot(mod.HSD_free_below, las=1, cex.axis = 0.8, tcl = -0.2)
par(mar = c(2, 2, 2, 5) + 0.5)
plot(mod.HSD_part_below, las=1, yaxt='n', cex.axis = 0.8, tcl = -0.2)
dev.off()




