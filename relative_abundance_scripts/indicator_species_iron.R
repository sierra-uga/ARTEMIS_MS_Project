ps_sub <- ps_noncontam_prev05 %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "Mitochondria" &
      Family   != "Chloroplast" &
      Order   != "Chloroplast" &
      Family  != "Mitochondria" 
  )

ps_sub <- subset_samples(ps_sub, Sample.Control == "True.Sample") %>% subset_samples(Transect_Name == "transect1") %>% subset_samples(sample_name != "STN12.3.732.fil.poly.S.r1") %>%
  phyloseq_validate() %>% tax_fix() %>% prune_taxa(taxa_sums(.) > 0, .) 

ps_sub <- tax_glom(ps_sub, "Genus", NArm = TRUE)

ps_sub <-label_duplicate_taxa(ps_sub,
                              tax_level = "Genus")

taxa_names <- otu_table(ps_sub)
taxa_names <- taxa_names[rowSums(taxa_names) !=0,
                         c(which(colSums(taxa_names) !=0))] # REMOVES 0s for CCA
taxa_names <- na.omit(taxa_names) # remove NAs just in case
otu_table(ps_sub) <- taxa_names # re-inserts the OTU table for the phyloseq object

## FREE-LIVING ##
ps_free <- ps_sub %>% subset_samples(Filter_pores == "free-living") %>% 
  prune_taxa(taxa_sums(.) > 0, .) 

ps_free <- ps_free %>%
  tax_glom(taxrank = "Genus") %>% # agglomerate at Order level, can change to different taxonomic level!# %>%
  transform_sample_counts(function(x) {x/sum(x)})  # Transform to rel. abundance (normalize data

## upper 200m

#### FREE-LIVING #####

write.csv(ps_free@otu_table,'M3.SP.otus_free.csv')
taxa_free <- as.data.frame(ps_free@tax_table)
taxa_free$taxon <- rownames(taxa_free)
#Import phyloseq OTU table as an OTU table/dataframe
SpOTU_free <-read.csv('M3.SP.otus_free.csv')

#do some shuffling of the OTU table
SpOTUFlip_free <- as.data.frame(t(SpOTU_free)) #makes it a dataframe and puts x into y and y into x (flips it)
names(SpOTUFlip_free) <- as.matrix(SpOTUFlip_free[1, ]) # renames columns
SpOTUFlip_free<- SpOTUFlip_free[-1, ] #removes first row
SpOTUFlip_num_free<-as.data.frame(lapply(SpOTUFlip_free, as.numeric)) #convert from character to number
SpOTUFlip_num_free$sample_name<-row.names(SpOTUFlip_free) #puts row names as sample ID column
#OK now we have the OTU table that's somewhat in the way they like

#read in metadata
metadata_free <-ps_free@sam_data #read in metadata
head(metadata_free) # check

## Join based on SampleID
SpOTU_Final_free <-left_join(SpOTUFlip_num_free, metadata_free, by = c("sample_name" = "sample_name")) # join based on sample IDs, assuming they're the same for both OTU table and metadata

SPotus_free = SpOTU_Final_free[,1:223] #select just the ASV/OTU table part of the file (you may have to scroll to the back of the OTU file to find it...)
SPwat_free = SpOTU_Final_free$Event #the metadata column group you care about
SPotus_free <- SPotus_free[, colSums(SPotus_free) != 0]

indisp_free=multipatt(x=SPotus_free, cluster=SPwat_free, func = "r.g", print.perm = TRUE, control = how(nperm=9999))

#extract table of stats
indisp.sign_free <-as.data.table(indisp_free$sign, keep.rownames=TRUE)
#add adjusted p-value
indisp.sign_free[ ,p.value.bh:=p.adjust(p.value, method="BH")]
#now can select only the indicators with adjusted significant p-values
SPind_free <- indisp.sign_free[p.value.bh<=0.05, ]

SPind_free$taxon <- SPind_free$rn

merged_df_free <- merge(SPind_free, taxa_free, by = "taxon")

write.csv(merged_df_free,'indicator_Station_CDW_free_no_cont.csv')

### PARTICLE - ASSOCIATED

ps_part <- ps_sub %>% subset_samples(Filter_pores == "particle-associated") %>% 
 prune_taxa(taxa_sums(.) > 0, .) %>% subset_samples(watertype != "Other")

ps_part <- ps_part %>%
  tax_glom(taxrank = "Genus") %>% # agglomerate at Order level, can change to different taxonomic level!# %>%
  transform_sample_counts(function(x) {x/sum(x)})  # Transform to rel. abundance (normalize data

## upper 200m

#### FREE-LIVING #####
write.csv(ps_part@otu_table,'M3.SP.otus_part.csv')
taxa_part <- as.data.frame(ps_part@tax_table)
taxa_part$taxon <- rownames(taxa_part)
#Import phyloseq OTU table as an OTU table/dataframe
SpOTU_part <-read.csv('M3.SP.otus_part.csv')

#do some shuffling of the OTU table
SpOTUFlip_part <- as.data.frame(t(SpOTU_part)) #makes it a dataframe and puts x into y and y into x (flips it)
names(SpOTUFlip_part) <- as.matrix(SpOTUFlip_part[1, ]) # renames columns
SpOTUFlip_part<- SpOTUFlip_part[-1, ] #removes first row
SpOTUFlip_num_part<-as.data.frame(lapply(SpOTUFlip_part, as.numeric)) #convert from character to number
SpOTUFlip_num_part$sample_name<-row.names(SpOTUFlip_part) #puts row names as sample ID column
#OK now we have the OTU table that's somewhat in the way they like

#read in metadata
metadata_part <-ps_part@sam_data #read in metadata
head(metadata_part) # check

## Join based on SampleID
SpOTU_Final_part <-left_join(SpOTUFlip_num_part, metadata_part, by = c("sample_name" = "sample_name")) # join based on sample IDs, assuming they're the same for both OTU table and metadata

SPotus_part = SpOTU_Final_part[,1:253] #select just the ASV/OTU table part of the file (you may have to scroll to the back of the OTU file to find it...)
SPwat_part = SpOTU_Final_part$Event #the metadata column group you care about
SPotus_part <- SPotus_part[, colSums(SPotus_part) != 0]

indisp_part=multipatt(x=SPotus_part, cluster=SPwat_part, func = "r.g", print.perm = TRUE, control = how(nperm=9999))

#extract table of stats
indisp.sign_part <-as.data.table(indisp_part$sign, keep.rownames=TRUE)
#add adjusted p-value
indisp.sign_part[ ,p.value.bh:=p.adjust(p.value, method="BH")]
#now can select only the indicators with adjusted significant p-values
SPind_part <- indisp.sign_part[p.value.bh<=0.05, ]

SPind_part$taxon <- SPind_part$rn

merged_df_part <- merge(SPind_part, taxa_part, by = "taxon")

write.csv(merged_df_part,'indicator_CDW_all_stations.csv')


### adonis 
bray_free <- phyloseq::distance(ps_sub, method = "bray") # setting distance
sampledf_free <- data.frame(sample_data(ps_sub))# make a data frame from the sample_data

#select from main data frame
adonis_frame_free <- dplyr::select(sampledf_free, Station, Location, Filter_pores)
adonis_frame_free$Filter_pores <- as.factor(adonis_frame_free$Filter_pores)

# Adonis test
adonis_free <- adonis2(bray_free ~ Location * Filter_pores, data = adonis_frame_free, permutations = 9999)
pairwise.adonis(bray_free, factors=adonis_frame_free$Filter_pores)


# Post hoc for location in polynya
beta_location_free <- betadisper(bray_free, adonis_frame_free$Station)

layout_matrix <- matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE)
layout(layout_matrix)

ps_free_above
# Adjust margins
par(mar = c(5, 6, 4, 2) + 0.1)
permutest(beta_location_free)
boxplot(beta_location_free)
mod.HSD_free <- TukeyHSD(beta_location_free)
plot(mod.HSD_free, las=1)
plot(beta_location_free)

#### particle-associated

bray_part <- phyloseq::distance(ps_part, method = "bray") # setting distance
sampledf_part <- data.frame(sample_data(ps_part))# make a data frame from the sample_data

#select from main data frame
adonis_frame_part <- dplyr::select(sampledf_part, Station, Salinity:CTD_Depth, Lab_NO3:DOC, Location:watertype)
adonis_frame_part$watertype <- as.factor(adonis_frame_part$watertype)

# Adonis test
adonis_part <- adonis2(bray_part ~ watertype, data = adonis_frame_part)

# Post hoc for location in polynya
beta_location_part <- betadisper(bray_part, adonis_frame_part$watertype)

layout_matrix <- matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE)
layout(layout_matrix)

# Adjust margins
pdf(file = "final_graphics/adonis_watermass.pdf", width = 12, height = 8) 
par(mfrow=c(1, 2), mar = c(5, 6, 4, 2) + 0.1)
plot(mod.HSD_free, las=1)
plot(mod.HSD_part, las=1)
dev.off()

permutest(beta_location_part)
boxplot(beta_location_part)
mod.HSD_part <- TukeyHSD(beta_location_part)
plot(mod.HSD_part, las=1)
plot(beta_location_part)