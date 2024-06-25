# load in libraries
library("dplyr")
library("tidyr")
library("phyloseq")
library("qiime2R")
library("ggplot2")
devtools::install_github("gmteunisse/fantaxtic")
library("fantaxtic")
library("vegan")
library("ggpubr")
install.packages("ggrepel")
library("ggrepel")
library("microbiome")
BiocManager::install("microbiome")
library("microbiome")

##################
# phloseq setup  #
##################

##################
#  data culling  #
##################

# select only bacteria, remove chloroplasts
ps_sub <- ps_noncontam_prev05 %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "Mitochondria" &
      Class   != "Chloroplast" &
      Order   != "Chloroplast" &
      Family  != "Mitochondria" 
  )

## filter by iron!

# Assuming ps_sub is your phyloseq object
ps_sub_depth <- subset_samples(ps_sub, More_Depth_Threshold %in% c("Bottom", "T-min"))
ps_sub_depth_no_polynya <- subset_samples(ps_sub_depth, Transect_Name %in% c("transect1"))
ps_sub_outflow <- subset_samples(ps_sub, True_Flow %in% c("Outflow"))

combined_ps <- merge_phyloseq(ps_sub_depth_no_polynya, ps_sub_outflow)
ps_sub_depth_no_polynya <- subset_samples(ps_sub_depth_no_polynya, !Location %in% c("Cont_Shelf"))
ps_sub_depth_no_polynya <- subset_samples(ps_sub_depth_no_polynya, !Location %in% c("Getz"))
temp <- sample_data(ps_sub)

temp$Local <- "Other"

# Update "Local" column where "Location" is "Open_polynya"
temp$Local[temp$Location == "Open_polynya"] <- "Open_Polynya"

temp$More_Depth_Threshold[temp$More_Depth_Threshold == "Mid-Bottom"] <- "T-min"
temp$More_Depth_Threshold[temp$More_Depth_Threshold == "Mid-Surface"] <- "Mixed Layer"

sample_data(ps_sub) <- temp

ps_sub <- ps_iron
ps_sub <- subset_samples(ps_sub, !Station %in% c("STN198"))
ps_sub <- subset_samples(ps_sub, !Station %in% c("STN174"))
ps_sub <- subset_samples(ps_sub, !Station %in% c("STN153"))

Samples_toRemove <- c("STN089.200.fil.dura.r2", "STN078.1040.pre.poly.3.LG")  # remove from phyloseq
ps_sub <- subset_samples(ps_sub, !(sample.illumina %in% Samples_toRemove)) # remove from phyloseq

ps_sub <- ps_sub %>% prune_taxa(taxa_sums(.) > 0, .) %>% prune_samples(sample_sums(.)>0, .)# remove 0 taxasums
ps_sub <- subset_samples(ps_sub, Sample.Control == "True.Sample") %>% tax_fix() %>% phyloseq_validate()
uwu <- otu_table(ps_sub)

uwu <- uwu[rowSums(uwu) !=0,
           c(which(colSums(uwu) !=0))] # REMOVES 0s for CCA

uwu <- na.omit(uwu)
otu_table(ps_sub) <- uwu

sample_data(ps_sub)$Latitude <- as.numeric(sample_data(ps_sub)$Latitude)
sample_data(ps_sub)$Longitude <- as.numeric(sample_data(ps_sub)$Longitude)
sample_data(ps_sub)$Temperature <- as.numeric(sample_data(ps_sub)$Temperature)
sample_data(ps_sub)$Sb_Oxygen <- as.numeric(sample_data(ps_sub)$Sb_Oxygen)
sample_data(ps_sub)$Salinity <- as.numeric(sample_data(ps_sub)$Salinity)
sample_data(ps_sub)$CTD_Depth <- as.numeric(sample_data(ps_sub)$CTD_Depth)
sample_data(ps_sub)$DOC <- as.numeric(sample_data(ps_sub)$DOC)
sample_data(ps_sub)$Lab_NO3 <- as.numeric(sample_data(ps_sub)$Lab_NO3)
sample_data(ps_sub)$Lab_NO2 <- as.numeric(sample_data(ps_sub)$Lab_NO2)
sample_data(ps_sub)$Lab_NH4 <- as.numeric(sample_data(ps_sub)$Lab_NH4)

ps_sub_genus <- tax_glom(ps_sub, "Genus", NArm = TRUE)
ps_sub_genus <- ps_sub

ps_sub_genus <-label_duplicate_taxa(ps_sub_genus,
                     tax_level = "Genus")


taxa_names(ps_sub_genus) <- tax_table(ps_sub_genus)[,"Genus"]

# free-living phyloseq
ps_free <- ps_sub_genus %>% subset_samples(Filter_pores == "free-living") #%>% prune_taxa(taxa_sums(.) > 0, .) 
ps_free <- combined_ps %>% subset_samples(Filter_pores == "0.2") %>% prune_taxa(taxa_sums(.) > 0, .) 

# particle-associated phyloseq
ps_part <- ps_sub_genus %>% subset_samples(Filter_pores == "particle-associated") #%>% prune_taxa(taxa_sums(.) > 0, .) 
ps_part <- combined_ps %>% subset_samples(Filter_pores >= "2") %>% prune_taxa(taxa_sums(.) > 0, .) 

ps_free <- subset_samples(ps_free, DOC != "NA")
ps_part <- subset_samples(ps_part, DOC != "NA")

ps_free <- subset_samples(ps_free, Iron != "NA")
ps_part <- subset_samples(ps_part, Iron != "NA")


ps_free <- subset_samples(ps_free, Location != "Open_polynya")
ps_part <- subset_samples(ps_part, Location != "Open_polynya")

ps_free <- subset_samples(ps_free, Location != "Cont_Shelf")
ps_part <- subset_samples(ps_part, Location != "Cont_Shelf")

ps_free <- subset_samples(ps_free, Station != "STN002")
ps_free <- subset_samples(ps_free, Station != "STN004")
ps_part <- subset_samples(ps_part, Station != "STN002")
ps_part <- subset_samples(ps_part, Station != "STN004")



# top taxa from all taxa
top_all <- top_taxa(ps_sub, 
                    tax_level = "Order", 
                    n_taxa = 10)

# selecting the top free-living taxa
top_free <- top_taxa(ps_free, 
                     tax_level = "Genus",
                     n_taxa = 15)

# selecting the top particle-associated taxa
top_part <- top_taxa(ps_part, 
                     tax_level = "Genus", 
                     n_taxa = 15)

top_taxa_free <- as.data.frame(top_free[["top_taxa"]])
top_taxa_part <- as.data.frame(top_part[["top_taxa"]])



##################
# CCA ordination #
##################

p

# Ordinate 2
cca <- ordinate(ps_sub_depth_no_polynya, method = "CCA", formula = ~ Latitude + Longitude + Salinity + Temperature + Sb_Oxygen)
cca1 <- ordinate(ps_free, method = "CCA", formula = ~ CTD_Depth + Salinity + Temperature + Sb_Oxygen + DOC + Lab_NO3 + Lab_NH4 + Lab_NO2 + Iron)
cca2 <- ordinate(ps_part, method = "CCA", formula = ~ CTD_Depth + Salinity + Temperature + Sb_Oxygen + DOC + Lab_NO3 + Lab_NH4 + Lab_NO2 + Iron)
# Ordinate 3
top_cca <- ordinate(top_all$ps_obj, "CCA", formula= ~ CTD_Depth + Salinity + Temperature + Sb_Oxygen + Iron)
top_cca1 <- ordinate(top_free$ps_obj, "CCA", formula= ~ CTD_Depth + Salinity + Temperature + Sb_Oxygen + Iron + Lab_NO3 + Lab_NH4 + Lab_NO2)
top_cca2 <- ordinate(top_part$ps_obj, "CCA", formula= ~ CTD_Depth + Salinity + Temperature + Sb_Oxygen + Iron + Lab_NO3 + Lab_NH4 + Lab_NO2)


ord_explore(ps_free)
ord_explore(ps_sub)
## aesthetics
# set ggplot shortcut to remove grid
remove_grid <- theme(legend.position = "bottom") + theme_bw() + # removes grid
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

# color for ggplot points
color_breaks <- c("Dotson", "Eastern_CC", "Cont_Shelf", "Western_CC", "Getz", "Open_polynya")
color_point <- scale_color_manual(values = c("#5AD0FC", "darkred","#A3DCA5", "red", "#006B93", "#09A20D"),
                                  name = "Location",
                                  breaks = color_breaks,
                                  labels = color_breaks)

color_point <- scale_color_manual(values = c("darkred", "yellow"),
                                  name = "Iron Level",
                                  breaks = c("High", "Low"),
                                  labels = c("High", "Low"))

color_point <- scale_color_manual(values = c("seagreen", "dodgerblue","purple", "darkred", "red2"),
                                  name = "Depth Threshold",
                                  breaks = color_breaks,
                                  labels = color_breaks)

cols <- rev(rainbow(7)[-7])
iron_point <- scale_color_gradientn(colors = c(cols, "darkred"), name = "dFe (nmol/kg)")
iron_point <- scale_color_gradient(low = "pink", high = "darkred",
                                    name = "dFe (nmol/kg)")


shapes <- scale_shape_manual(values = c(16,4,13),
                                  name = "Area",
                                  breaks = c("Outflow", "Inflow", "NA"),
                                  labels = c("Outflow", "Inflow", "NA"))

# Both free-living + particle-associated CCA
all_CCA <- plot_ordination(ps_sub, cca,
                           type = "samples", color="Station", shape="Filter_pores")
all_CCA <- all_CCA + geom_point(size = 2) + remove_grid
all_CCA
ggsave("graphics/all_CCA.pdf", width = 7, height = 6, dpi = 150)

# Free-living CCA 
free_CCA <- plot_ordination(ps_free, cca1,
                            type = "sites", col="More_Depth_Threshold")
f_CCA <- free_CCA +  geom_point(pch=16, size=2) + remove_grid +# color_point +#geom_label(size= 2, aes(label = Genus)) +  # X for other locations
 theme(text = element_text(family = "Helvetica"))
f_CCA

# Particle-associated CCA
part_CCA <- plot_ordination(ps_part, cca2,
                            type = "sites", col="More_Depth_Threshold")
p_CCA <- part_CCA +  geom_point(pch=16, size=2) + remove_grid + #+ color_point + # X for other locations
  theme(text = element_text(family = "Helvetica"))
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
library(ggrepel)

arrowmat = vegan::scores(cca, display = "bp")
labels <- c("Dotson", "Eastern CC", "Getz", "Open Polynya", "Western CC", "Salinity", "Temperature", "Oxygen", "Latitude", "Longitude")
# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
# Define the arrow aesthetic mapping
arrow_map = aes(xend = 1.6 * CCA1, yend = 1.6 * CCA2, x = 0, y = 0, shape = NULL, color = NULL, 
                label = labels)
label_map = aes(x = 1.8 * CCA1, y = 1.8 * CCA2, shape = NULL, color = NULL, 
                label = labels)
# Make a new graphic
arrowhead = arrow(length = unit(0.02, "npc"))
p0 = all_CCA + geom_segment(arrow_map, size = 0.6, data = arrowdf, color = "black", 
                            arrow = arrowhead) + geom_text_repel(label_map, size = 3, data = arrowdf, max.overlaps =  30)
p0

arrowmat = vegan::scores(cca1, display = "bp")
labels <- c("Depth", "Salinity", "Temperature", "Oxygen", "DOC") #"Dotson", "Eastern CC", "Getz", "Open Polynya", "Western CC", 
labels <- rownames(arrowmat)
# Add labels, make a data.frame
arrowdf <- data.frame(labels = labels, arrowmat)
# Define the arrow aesthetic mapping
arrow_map = aes(xend = 1.6 * CCA1, yend = 1.6 * CCA2, x = 0, y = 0, shape = NULL, color = NULL, 
                label = labels)
label_map = aes(x = 1.8 * CCA1, y = 1.8 * CCA2, shape = NULL, color = NULL, 
                label = labels)
# Make a new graphic
arrowhead = arrow(length = unit(0.02, "npc"))
p1 = f_CCA + geom_segment(arrow_map, size = 0.6, data = arrowdf, color = "black", 
                          arrow = arrowhead) + geom_text_repel(label_map, size = 3, data = arrowdf, max.overlaps =  30)
p1


arrowmat = vegan::scores(cca2, display = "bp")
#labels <- rownames(arrowmat)
# Add labels, make a data.frame
arrowdf <- data.frame(labels = labels, arrowmat)
# Define the arrow aesthetic mapping
arrow_map = aes(xend = 1.6 * CCA1, yend = 1.6 * CCA2, x = 0, y = 0, shape = NULL, color = NULL, 
                label = labels)
label_map = aes(x = 1.8 * CCA1, y = 1.8 * CCA2, shape = NULL, color = NULL, 
                label = labels)
# Make a new graphic
arrowhead = arrow(length = unit(0.02, "npc"))
p2 = p_CCA + geom_segment(arrow_map, size = 0.6, data = arrowdf, color = "black", 
                          arrow = arrowhead, arrow.fill="black") + geom_text_repel(label_map, size = 3, data = arrowdf)
p2

# use CCA_ord function (see files) to make graphs
final_free_plot <- CCA_ord(top_cca1, free_CCA, "CCA Top 10 Taxa (Order)", final_free_plot)
final_part_plot <- CCA_ord(top_cca2, part_CCA, "CCA Top 10 Taxa (Genus)", final_part_plot)

# split plot for both communities
combined <- ggarrange(
  p1, p2, labels = c("F", "P"),
  common.legend = TRUE, legend = "right"
)
annotate_figure(combined, top = text_grob("CCA of Samples by Depth and Env. Factors", 
                                             color = "black", face = "bold", size = 16, family = "Helvetica"))
ggsave("ordination_scripts/graphics/CCA_Depth_with_DOC_no_open.pdf", width = 11, height = 8, dpi = 150)

#canonical correspondence analysis

arrowmat1 <- scores(cca1, display = "species")
arrowmat2 <- scores(cca2, display = "species")
arrowdf1 <- data.frame(labels = rownames(arrowmat1), arrowmat1)

# Define the arrow aesthetic mapping
arrow_map1 = aes(xend = CCA1, yend = CCA2, x = 0, y = 0, shape = NULL, color = NULL, 
                 label = labels)

# Add labels, make a data.frame
arrowdf_free <- data.frame(labels = rownames(arrowmat1), arrowmat1)

arrowdf_part <- data.frame(labels = rownames(arrowmat2), arrowmat2)


arrowmat = vegan::scores(top_cca1, display = "species")
# Add labels, make a data.frame
labels <- rownames(arrowmat)
arrowdf <- data.frame(labels = labels, arrowmat)
# Define the arrow aesthetic mapping
arrow_map = aes(xend = 1.6 * CCA1, yend = 1.6 * CCA2, x = 0, y = 0, shape = NULL, color = NULL, 
                label = labels)
label_map = aes(x = 1.8 * CCA1, y = 1.8 * CCA2, shape = NULL, color = NULL, 
                label = labels)
# Make a new graphic
arrowhead = arrow(length = unit(0.01, "npc"))
p1 = f_CCA + geom_segment(arrow_map, size = 0.6, data = arrowdf, color = "gray20", 
                            arrow = arrowhead) + geom_text_repel(label_map, size = 2.5, data = arrowdf, max.overlaps = 20)
p1

arrowmat = vegan::scores(top_cca2, display = "species")
labels <- rownames(arrowmat)
# Add labels, make a data.frame
arrowdf <- data.frame(labels = labels, arrowmat)
# Define the arrow aesthetic mapping
arrow_map = aes(xend = 1.6 * CCA1, yend = 1.6 * CCA2, x = 0, y = 0, shape = NULL, color = NULL, 
                label = labels)
label_map = aes(x = 1.8 * CCA1, y = 1.8 * CCA2, shape = NULL, color = NULL, 
                label = labels)
# Make a new graphic
arrowhead = arrow(length = unit(0.01, "npc"))
p2 = p_CCA + geom_segment(arrow_map, size = 0.6, data = arrowdf, color = "gray20", 
                          arrow = arrowhead) + geom_text_repel(label_map, size = 2.5, data = arrowdf, max.overlaps =  30)
p2