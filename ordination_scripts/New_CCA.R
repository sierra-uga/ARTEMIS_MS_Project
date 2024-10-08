install.packages("ggtext")
library(ggtext)
library(ggrepel)
library(microViz)
## phloyseq setup (from phyloseq_setup.R)
# filtering
ps_sub <- ps_noncontam_prev05 %>% # filtering unwanted chloroplast/mitochondria
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "Mitochondria" &
      Class   != "Chloroplast" &
      Order   != "Chloroplast" &
      Family  != "Mitochondria" 
  )

ps_sub <- subset_samples(ps_sub, Sample.Control == "True.Sample") %>% # filtering by actual samples
  tax_fix() %>% 
  phyloseq_validate()

# fixing 0's for CCA
taxa_names <- otu_table(ps_sub)
taxa_names <- taxa_names[rowSums(taxa_names) !=0,
           c(which(colSums(taxa_names) !=0))] # REMOVES 0s for CCA
taxa_names <- na.omit(taxa_names) # remove NAs just in case
otu_table(ps_sub) <- taxa_names # re-inserts the OTU table for the phyloseq object

# changing environmental factors to numeric (archaic :) )
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
sample_data(ps_sub)$Iron <- as.numeric(sample_data(ps_sub)$Iron)

ps_sub <- subset_samples(ps_sub, Iron != "NA")
ps_sub <- subset_samples(ps_sub, DOC != "NA")

#ps_sub <- subset_samples(ps_sub, watertype == "CDW")

ps_sub <- subset_samples(ps_sub, Location != "Cont_Shelf") # if want to remove cont_shelf
#ps_sub <- subset_samples(ps_sub, Location != "Open_polynya") # if want to remove open polynya

# free-living
ps_free <- ps_sub %>% subset_samples(Filter_pores == "free-living") 

# particle-associated
ps_part <- ps_sub %>% subset_samples(Filter_pores == "particle-associated")

## CCA analysis
# free-living
## upper 200m
ps_free_above <- ps_free %>% subset_samples(CTD_Depth <= 200)
cca_free_above <- ordinate(ps_free_above, method = "CCA", formula = ~ CTD_Depth + Salinity + Temperature + Sb_Oxygen + Lab_NO3 + Lab_NH4 + Lab_NO2 + DOC)

## bellow 200m
ps_free_below <- ps_free %>% subset_samples(CTD_Depth >= 200)
cca_free_below <- ordinate(ps_free_below, method = "CCA", formula = ~ CTD_Depth +  Salinity + Temperature + Sb_Oxygen + Lab_NO3 + Lab_NH4 + Lab_NO2 + DOC)
# particle-associated
## upper 200m
ps_part_above <- ps_part %>% subset_samples(CTD_Depth <= 200)
cca_part_above <- ordinate(ps_part_above, method = "CCA", formula = ~ CTD_Depth +  Salinity + Temperature + Sb_Oxygen + Lab_NO3 + Lab_NH4 + Lab_NO2 + DOC)

## bellow 200m
ps_part_below <- ps_part %>% subset_samples(CTD_Depth >= 200)
cca_part_below <- ordinate(ps_part_below, method = "CCA", formula = ~ CTD_Depth +  Salinity + Temperature + Sb_Oxygen + Lab_NO3 + Lab_NH4 + Lab_NO2 + DOC)

remove_grid <- theme(legend.position = "bottom") + theme_bw() + # removes grid
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

# color for ggplot points
color_breaks <- c("Dotson", "Eastern_CC", "Western_CC", "Getz", "Open_polynya", "West_OP")
color_point <- scale_fill_manual(values = c("#5AD0FC", "darkred", "red", "#006B93", "#09A20D", "purple", "black"),
                                  name = "Location in Polynya",
                                  breaks = color_breaks,
                                  labels = c("Dotson", "East CC", "West CC", "Getz", "Open Polynya", "West OP"))

color_breaks <- c("Dotson", "Eastern_CC", "Western_CC", "Getz", "Open_polynya", "West_OP")
color_point2 <- scale_color_manual(values = c("#5AD0FC", "darkred", "red", "#006B93", "#09A20D", "purple", "black"),
                                 name = "Location in Polynya",
                                 breaks = color_breaks,
                                 labels = c("Dotson", "East CC", "West CC", "Getz", "Open Polynya", "West OP"))

shapes <- scale_shape_manual(values = c(8 , 13, 23, 22, 21),
                             name = "Water Mass",
                             breaks = c("AASW", "AASW-WW", "WW", "WW-CDW", "CDW"),
                             labels = c("AASW", "AASW-WW", "WW", "WW-CDW", "CDW"))

## ordination
# free-living
## upper 200m
free_above_CCA <- plot_ordination(ps_free_above, cca_free_above,
                            type = "sites", color="Location")
f_CCA <- free_above_CCA + geom_point(aes(fill = Location), size = 3, alpha = 0.7) + ggtitle("Free-living") +remove_grid + color_point + color_point2 +#geom_label(size= 2, aes(label = Genus)) +  # X for other locations
  theme(text = element_text(family = "Helvetica"), plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
  annotate("richtext", x = Inf, y = Inf, label = "Above 200m",
           hjust = 1.1, vjust = 1.5, size = 4, color = "black",
           fontface = "bold", fill = "white") +
  theme(plot.margin = margin(10, 10, 10, 10, "pt"))
f_CCA

arrowmat = vegan::scores(cca_free_above, display = "bp")
labels <- c("Depth", "Salinity", "Temperature", "Oxygen", "Nitrate", "Ammonia", "Nitrite", "DOC") #"Dotson", "Eastern CC", "Getz", "Open Polynya", "Western CC", 
#labels <- rownames(arrowmat)
# Add labels, make a data.frame
arrowdf <- data.frame(labels = labels, arrowmat)
# Define the arrow aesthetic mapping
arrow_map = aes(xend = 1.6 * CCA1, yend = 1.6 * CCA2, x = 0, y = 0, shape = NULL, color = NULL, 
                label = labels)
label_map = aes(x = 2 * CCA1, y = 2 * CCA2, shape = NULL, color = NULL, 
                label = labels)
# Make a new graphic
arrowhead = arrow(length = unit(0.02, "npc"))
f1 = f_CCA + geom_segment(arrow_map, size = 0.6, data = arrowdf, color = "#b54d04", 
                          arrow = arrowhead) + geom_text_repel(label_map, size = 3.5, data = arrowdf, max.overlaps =  30, color = "black")
f1

## lower 200m
free_below_CCA <- plot_ordination(ps_free_below, cca_free_below,
                            type = "sites", color="Location")
f_CCA2 <- free_below_CCA + geom_point(aes(fill = Location), size=3, alpha=0.7) + remove_grid + color_point + color_point2+#geom_label(size= 2, aes(label = Genus)) +  # X for other locations
  theme(text = element_text(family = "Helvetica")) +
  annotate("richtext", x = Inf, y = Inf, label = "Below 200m",
           hjust = 1.1, vjust = 1.5, size = 4, color = "black",
           fontface = "bold", fill = "white") +
  theme(plot.margin = margin(10, 10, 10, 10, "pt"))
f_CCA2

arrowmat = vegan::scores(cca_free_below, display = "bp")
labels <- c("Depth", "Salinity", "Temperature", "Oxygen", "Nitrate", "Ammonia", "Nitrite", "DOC") #"Dotson", "Eastern CC", "Getz", "Open Polynya", "Western CC", 
#labels <- rownames(arrowmat)
# Add labels, make a data.frame
arrowdf <- data.frame(labels = labels, arrowmat)
# Define the arrow aesthetic mapping
arrow_map = aes(xend = 1.6 * CCA1, yend = 1.6 * CCA2, x = 0, y = 0, shape = NULL, color = NULL, 
                label = labels)
label_map = aes(x = 2 * CCA1, y = 2 * CCA2, shape = NULL, color = NULL, 
                label = labels)
# Make a new graphic
arrowhead = arrow(length = unit(0.02, "npc"))
f2 = f_CCA2 + geom_segment(arrow_map, size = 0.6, data = arrowdf, color = "#b54d04", 
                          arrow = arrowhead) + geom_text_repel(label_map, size = 3.5, data = arrowdf, max.overlaps =  30, color = "black")
f2


# Particle-associated CCA
part_above_CCA <- plot_ordination(ps_part_above, cca_part_above,
                                  type = "sites", color="Location")
p_CCA <- part_above_CCA + geom_point(aes(fill = Location), size = 3, alpha = 0.7) + ggtitle("Particle-associated") + remove_grid + color_point + color_point2 +#geom_label(size= 2, aes(label = Genus)) +  # X for other locations
  theme(text = element_text(family = "Helvetica"), plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
  annotate("richtext", x = Inf, y = Inf, label = "Above 200m",
           hjust = 1.1, vjust = 1.5, size = 4, color = "black",
           fontface = "bold", fill = "white") +
  theme(plot.margin = margin(10, 10, 10, 10, "pt"))
p_CCA

arrowmat = vegan::scores(cca_part_above, display = "bp")
labels <- c("Depth", "Salinity", "Temperature", "Oxygen", "Nitrate", "Ammonia", "Nitrite", "DOC") #"Dotson", "Eastern CC", "Getz", "Open Polynya", "Western CC", 
#labels <- rownames(arrowmat)
# Add labels, make a data.frame
arrowdf <- data.frame(labels = labels, arrowmat)
# Define the arrow aesthetic mapping
arrow_map = aes(xend = 1.6 * CCA1, yend = 1.6 * CCA2, x = 0, y = 0, shape = NULL, color = NULL, 
                label = labels)
label_map = aes(x = 2 * CCA1, y = 2 * CCA2, shape = NULL, color = NULL, 
                label = labels)

# Make a new graphic
arrowhead = arrow(length = unit(0.02, "npc"))
p1 = p_CCA + geom_segment(arrow_map, size = 0.6, data = arrowdf, color = "#b54d04", 
                          arrow = arrowhead) + geom_text_repel(label_map, size = 3.5, data = arrowdf, max.overlaps =  30, color = "black")
p1

## lower 200m
part_below_CCA <- plot_ordination(ps_part_below, cca_part_below,
                                  type = "sites", color="Location")
p_CCA2 <- part_below_CCA + geom_point(aes(fill = Location), size=3, alpha=0.7) + remove_grid + color_point + color_point2+#geom_label(size= 2, aes(label = Genus)) +  # X for other locations
  theme(text = element_text(family = "Helvetica"))+
  annotate("richtext", x = Inf, y = Inf, label = "Below 200m",
           hjust = 1.1, vjust = 1.5, size = 4, color = "black",
           fontface = "bold", fill = "white") +
  theme(plot.margin = margin(10, 10, 10, 10, "pt"))
p_CCA2

arrowmat = vegan::scores(cca_part_below, display = "bp")
labels <- c("Depth", "Salinity", "Temperature", "Oxygen", "DOC", "Nitrate", "Ammonia", "Nitrite") #"Dotson", "Eastern CC", "Getz", "Open Polynya", "Western CC", 
#labels <- rownames(arrowmat)
# Add labels, make a data.frame
arrowdf <- data.frame(labels = labels, arrowmat)
# Define the arrow aesthetic mapping
arrow_map = aes(xend = 1.6 * CCA1, yend = 1.6 * CCA2, x = 0, y = 0, shape = NULL, color = NULL, 
                label = labels)
label_map = aes(x = 2 * CCA1, y = 2 * CCA2, shape = NULL, color = NULL, 
                label = labels)
# Make a new graphic
arrowhead = arrow(length = unit(0.02, "npc"))
p2 = p_CCA2 + geom_segment(arrow_map, size = 0.6, data = arrowdf, color = "#b54d04", 
                           arrow = arrowhead) + geom_text_repel(label_map, size = 3.5, data = arrowdf, max.overlaps =  30, color = "black")
p2

### make a fake plot to grab the legend from
data <- data.frame(
  x = rnorm(12),
  y = rnorm(12),
  Location = rep(c("Dotson", "Eastern_CC", "Western_CC", "Getz", "Open_polynya", "West_OP"), each = 2)
)

fake_plot <- ggplot(data, aes(x = x, y = y)) +
  geom_point(aes(color = Location), size = 4) +
  color_point2 +
  theme_void() +  # Hides the axes and background
  theme(legend.position = "right")

legend_combined <- get_legend(fake_plot)

## combine plots

library(ggpubr)
ggarrange(f1, legend.grob=legend_combined, legend="right", labels = c("A", "B", "C", "D"))

ggsave("final_graphics/LOCATION_LEGEND.pdf", width = 10, height = 8, dpi = 1200)
