# load in libraries
library("dplyr")
library("tidyr")
library("phyloseq")
library("qiime2R")
library("ggplot2")

### Phlyoseq object set up
# feature table
ASV <- qza_to_phyloseq(features="table.qza")

# read in metadata
metatable <- read.delim("artemis-eDNA-metadata-final.tsv", sep="\t", header=TRUE, row.names="Sample.illumina") 
# fixing metadata
metatable <- metatable[-159,] # removes row 89_200_FIL_R2 because taxasum was 0.
metatable <- filter(metatable, Sample.Control == "True.Sample") %>% 
  mutate_at(c(20:27, 29:38, 44), as.numeric) # remove blanks, convert character columns for env. data to numeric

# importing (continued)
row.names(metatable) <- metatable[["SampleID"]]
metatable <- metatable %>% select(SampleID, everything())
new_metatable <- metatable %>% 
  mutate(Community = replace(Filter_pores, Filter_pores == 0.2, "Free-living")) %>%
  mutate(Community = replace(Filter_pores, Filter_pores >= 2, "Particle-associated")) %>%
  mutate(Community = replace(Community, Community =="0.2", "Free-living")) 
# metatable <- metatable[-(which(metatable$Station %in% c("STN153", "STN174", "STN089", "STN151.2", "STN198"))),] 
META <- sample_data(new_metatable)

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
taxonomy[nm1] <- lapply(taxonomy[nm1], gsub, pattern = "D_.__", replacement = "") #remove beginning free of each column
taxmat <- as.matrix(taxonomy)
TAX <- tax_table(taxmat) #covert taxonomy table to phyloseq object

# import rooted tree
TREE <- qza_to_phyloseq(tree="rooted-tree.qza")

# merge all imported objects into phyloseq
ps <- merge_phyloseq(ASV, TAX, META, TREE)
ps

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
ps_free <- ps_sub %>% subset_samples(Filter_pores == "0.2")

# particle-associated phyloseq
ps_part <- ps_sub %>% subset_samples(Filter_pores >= "2")

######################################
# free-LIVING NMDS ordination + plot #
######################################

set.seed(1)

#unique station names = scale color
color_breaks <- unique(sample_data(ps_free)$watertype)

# Ordinate
ps_free_nmds <- ordinate(
  physeq = ps_free, 
  method = "NMDS", 
  distance = "bray"
)

ps_free_nmds_plot <- plot_ordination(
  physeq = ps_free,
  ordination = ps_free_nmds,
  color = "watertype",
  title = "NMS of free-living ASP Bacterial Communities based on Water Mass"
) +
  scale_color_manual(values = c("darkgreen", "dodgerblue", "red2", "blueviolet", "aquamarine3", "gray"),
                     name = "Water Mass",
                     breaks = color_breaks,
                     labels = color_breaks) +
  geom_point(shape = 16, aes(color = watertype), alpha = 0.9, size = 3) +
  annotate("text", x = -1.6, y = -1.7, label ="2D Stress: 0.113")

#takes away grid from ggplot
ps_free_nmds_plot + theme(legend.position = "bottom") + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave("graphics/free_living_NMS.pdf", width = 7, height = 6, dpi = 150)

##############################################
# PARTICLE-ASSOCIATED NMDS ordination + plot #
##############################################

set.seed(2)

#unique station names = scale color
color_breaks <- unique(sample_data(ps_part)$watertype)
shape_breaks

# Ordinate
ps_part_nmds <- ordinate(
  physeq = ps_part, 
  method = "NMDS", 
  distance = "bray"
)

ps_part_nmds_plot <- plot_ordination(
  physeq = ps_part,
  ordination = ps_part_nmds,
  color = "watertype",
  title = "NMS of particle-associated ASP Bacterial Communities based on Water Mass"
) +
  scale_color_manual(values = c("darkgreen", "dodgerblue", "red2", "blueviolet", "aquamarine3", "gray"),
                     name = "Water Mass",
                     breaks = color_breaks,
                     labels = color_breaks) +
  geom_point(shape = 16, aes(color = watertype), alpha = 0.9, size = 3) +
  annotate("text", x = -1.6, y = -1.7, label ="2D Stress: 0.118") + 
  theme(legend.position = "bottom") + theme_bw() + # removes grid
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave("graphics/particle_NMS.pdf", width = 7, height = 6, dpi = 150)


points_free <- ps_free_nmds[["points"]] 
ugga <- as.data.frame(points_free)


#### NMDS for all stations, transects as shapes
### NMDS for each transect, transect number as shape/color?
# coastal_current_name
# transect_name

# all_transect_ps_free <- ps_free %>% subset_samples(Transect_Name == "all_transect")

all_transect_ps_free <- ps_free
# Ordinate
all_transect_ps_free_nmds <- ordinate(
  physeq = all_transect_ps_free, 
  method = "NMDS", 
  distance = "bray"
)

all_transect_color_breaks <- c("transect1", "transect2")
#all_transect_shape_breaks <- unique(sample_data(all_transect_ps_free)$Transect_Number)

## ALL TRANSECT
#free-living
all_transect_ps_free_nmds_plot <- plot_ordination(
  physeq = all_transect_ps_free,
  ordination = all_transect_ps_free_nmds,
  color = "Transect_Name",
  title = "NMS of Free-living ASP Bacterial Communities based on Transect",
  shape = "Coastal_Current_Name"
) +
  scale_color_manual(values = c("red2", "dodgerblue2", "dodgerblue2", "blueviolet", "aquamarine3", "gray", "yellow", "green", "brown"),
                     name = "Transect 1 and 2",
                     breaks = all_transect_color_breaks,
                     labels = c("CDW Waterfall", "Meltwater Plume")) +
  scale_shape_manual(values = c(1, 12),
                     name = "Transect 3",
                     breaks = "transect3",
                     labels = "Coastal Current") +
  #guides(fill = guide_legend(title = "Sample Site", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  geom_point(shape = 16, aes(color = Transect_Name), alpha = 0.9, size = 3) +
  geom_point(color = "green3", aes(shape = Coastal_Current_Name), alpha = 1, size=3, stroke=1) + # CHANGE COLOR HERE
  ggrepel::geom_text_repel(aes(label = Transect_Number), size = 2, max.overlaps = 30) +
  ggrepel::geom_text_repel(nudge_x = -.04, nudge_y = .03, color = "green3", aes(label = Coastal_Current_Number), size = 2, max.overlaps = 30) +# adds transect numbers!!
  #geom_point(colour = "grey90", size = 1.5) +
  annotate("text", x = -1.6, y = -1.7, label ="2D Stress: 0.113")

#takes away grid from ggplot
all_transect_ps_free_nmds_plot + theme(legend.position = "bottom") + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave("graphics/all_transect_free_living_NMS.pdf", width = 9, height = 7, dpi = 150)

######################

#particle-associated
all_transect_ps_free <- ps_free
# Ordinate
all_transect_ps_free_nmds <- ordinate(
  physeq = all_transect_ps_free, 
  method = "NMDS", 
  distance = "bray"
)

all_transect_color_breaks <- c("transect1", "transect2")
#all_transect_shape_breaks <- unique(sample_data(all_transect_ps_free)$Transect_Number)

## ALL TRANSECT
all_transect_ps_part_nmds_plot <- plot_ordination(
  physeq = all_transect_ps_part,
  ordination = all_transect_ps_part_nmds,
  color = "Transect_Name",
  title = "NMS of Particle-Associated ASP Bacterial Communities based on Transect",
  shape = "Coastal_Current_Name"
) +
  scale_color_manual(values = c("red2", "dodgerblue2", "dodgerblue2", "blueviolet", "aquamarine3", "gray", "yellow", "green", "brown"),
                     name = "Transect 1 and 2",
                     breaks = all_transect_color_breaks,
                     labels = c("CDW Waterfall", "Meltwater Plume")) +
  scale_shape_manual(values = c(1, 12),
                     name = "Transect 3",
                     breaks = "transect3",
                     labels = "Coastal Current") +
  #guides(fill = guide_legend(title = "Sample Site", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  geom_point(shape = 16, aes(color = Transect_Name), alpha = 0.9, size = 3) +
  ggrepel::geom_text_repel(aes(label = Transect_Number), size = 2, max.overlaps = 30) +
  ggrepel::geom_text_repel(nudge_x = -.04, nudge_y = .03, color = "green3", aes(label = Coastal_Current_Number), size = 2, max.overlaps = 30) +# adds transect numbers!!
  #geom_point(colour = "grey90", size = 1.5) +
  geom_point(color = "green3", aes(shape = Coastal_Current_Name), alpha = 1, size=3, stroke=1) + # CHANGE COLOR HERE
  annotate("text", x = -1.6, y = -1.7, label ="2D Stress: 0.113")

#takes away grid from ggplot
all_transect_ps_part_nmds_plot + theme(legend.position = "bottom") + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave("graphics/all_transect_part_associated_NMS.pdf", width = 9, height = 7, dpi = 150)


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

#### CDW WATERFALL - TRANSECT 1 #### 
waterfall_ps_all <- ps_sub %>% subset_samples(Transect_Name == "transect1")

# Ordinate
waterfall_ps_all_nmds <- ordinate(
  physeq = waterfall_ps_all, 
  method = "NMDS", 
  distance = "bray"
)

waterfall_color_breaks <- c("transect1", "transect2")
#waterfall_shape_breaks <- unique(sample_data(waterfall_ps_all)$Transect_Number)

## ALL TRANSECT
waterfall_ps_all_nmds_plot <- plot_ordination(
  physeq = waterfall_ps_all,
  ordination = waterfall_ps_all_nmds,
  color = "Transect_Number",
  title = "NMS of All Communites for CDW Waterfall (Transect 1)",
  shape = "Community"
) +
  scale_color_gradient(low="pink", high="red",
                     name = "Transect 1")+
  scale_shape_manual(values = c(16, 17),
                     name = "Community",
                     breaks = c("Free-living", "Particle-associated"),
                     labels = c("Free-living", "Particle-associated")) +
  #guides(fill = guide_legend(title = "Sample Site", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  theme(legend.title = element_text(face = "bold", size =12), legend.title.align = 0.5) + # # adjusts the title of the legend
  #geom_point(shape = 16, aes(color = Transect_Name), alpha = 0.9, size = 3) +
  ggrepel::geom_text_repel(aes(label = Transect_Number), size = 2, max.overlaps = 30) +
  #ggrepel::geom_text_repel(nudge_x = -.04, nudge_y = .03, color = "green3", aes(label = Coastal_Current_Number), size = 2, max.overlaps = 30) +# adds transect numbers!!
  #geom_point(colour = "grey90", size = 1.5) +
  geom_point(aes(shape = Community), alpha = 0.9, size = 3) +
  #geom_point(color = "green3", aes(shape = Coastal_Current_Name), alpha = 1, size=3, stroke=1) + # CHANGE COLOR HERE
  annotate("text", x = -1.6, y = -1.7, label ="2D Stress: 0.113")

#takes away grid from ggplot
waterfall_ps_all_nmds_plot + theme(legend.position = "bottom") + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave("graphics/waterfall_all_NMS.pdf", width = 9, height = 7, dpi = 150)


#### TRANSECT 1 and 2 INFLOW/OUTFLOW #### 
# transect 1 + 2
transect_num <- c("transect1", "transect2")  # remove from phyloseq
inflow_outflow_ps <- subset_samples(ps_sub, (Transect_Name %in% transect_num)) # remove from phyloseq

# free-living phyloseq
inflow_outflow_ps_free <- inflow_outflow_ps %>% subset_samples(Filter_pores == "0.2")

# particle-associated phyloseq
inflow_outflow_ps_part <- inflow_outflow_ps %>% subset_samples(Filter_pores >= "2")

inflow_outflow_ps_free_transect1 <- inflow_outflow_ps_free %>% subset_samples(Transect_Name == "transect1")
inflow_outflow_ps_free_transect2 <- inflow_outflow_ps_free %>% subset_samples(Transect_Name == "transect2")
### FREE-LIVING
# Ordinate
inflow_outflow_ps_free_nmds <- ordinate(
  physeq = inflow_outflow_ps_free, 
  method = "NMDS", 
  distance = "bray"
)

waterfall_color_breaks <- c("transect1", "transect2")
#waterfall_shape_breaks <- unique(sample_data(waterfall_ps_all)$Transect_Number)
library(ggnewscale)


## ALL TRANSECT
inflow_outflow_ps_free_nmds_plot <- plot_ordination(
  physeq = inflow_outflow_ps_free,
  ordination = inflow_outflow_ps_free_nmds,
  color = "Transect_Name",
  title = "NMS of Free-living Communities at Inflow/Outflow Stations"
) +
  geom_point(aes(color=Transect_Name), size = 2) + 
  scale_fill_gradient(name = "Transect1",low="pink", high="red",
                      breaks=c("1","2","3","4","5")) +
  # guides(fill = guide_legend(title = "Sample Site", reverse = FALSE, keywidth = 1, keyheight = 1)) + # Plot the legend
  # ggrepel::geom_text_repel(aes(label = Transect_Number), size = 2, max.overlaps = 30) +
  # ggrepel::geom_text_repel(nudge_x = -.04, nudge_y = .03, color = "green3", aes(label = Coastal_Current_Number), size = 2, max.overlaps = 30) +# adds transect numbers!!
  # geom_point(colour = "grey90", size = 1.5) +
  ggrepel::geom_text_repel(aes(label = Transect_Number), size = 2, max.overlaps = 30) +
  # geom_point(aes(color = Transect_Number), alpha = 0.9, size = 3) +
  # geom_point(color = "green3", aes(shape = Coastal_Current_Name), alpha = 1, size=3, stroke=1) + # CHANGE COLOR HERE
  annotate("text", x = -1.6, y = -1.7, label ="2D Stress: 0.113")

#takes away grid from ggplot
waterfall_ps_all_nmds_plot + theme(legend.position = "bottom") + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggsave("graphics/waterfall_all_NMS.pdf", width = 9, height = 7, dpi = 150)

# NMS FOR free-living
# NMS for particle associated

