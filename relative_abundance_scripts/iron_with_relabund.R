# metatable 
metatable <- read.delim("required_files/artemis-eDNA-metadata-final.tsv", sep="\t", header=TRUE) 
#metatable <- filter(metatable, Sample.Control == "True.Sample")# filter by transect
metatable$is.neg <- metatable$Sample.Control == "Control.Sample"
metatable$Final_Qubit <- as.numeric(metatable$Final_Qubit) 
metatable <- metatable %>% filter(., Iron != "NA") # for iron
metatable$More_Depth_Threshold[metatable$More_Depth_Threshold == "Mid-Bottom"] <- "T-min"
metatable$More_Depth_Threshold[metatable$More_Depth_Threshold == "Mid-Surface"] <- "Mixed Layer"
metatable$Filter_pores <- ifelse(metatable$Filter_pores >= 0.2, "free-living", 
                                    ifelse(metatable$Filter_pores == 3.0 & metatable$Filter_pores <= 2.0, "particle-associated", metatable$Filter_pores))
metatable$Siderophore <- as.numeric(metatable$Siderophore)

agg_iron <- metatable %>%
  group_by(Station, More_Depth_Threshold) %>%
  summarize(mean_iron = mean(Iron, na.rm = TRUE))

agg_iron$More_Depth_Threshold[agg_iron$More_Depth_Threshold == "Mid-Bottom"] <- "T-min"
agg_iron$More_Depth_Threshold[agg_iron$More_Depth_Threshold == "Mid-Surface"] <- "Mixed Layer"

level_order <- c("STN198", "STN012", "STN089", "STN132", "STN106", "STN20", "STN014", "STN078", "STN056a", "STN056b", "STN22", "STN068", "STN174", "STN153")

# phyloseq
ps_sub1 <- subset_samples(ps_noncontam_prev05, Sample.Control == "True.Sample")

ps_sub <- ps_sub1 %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "Mitochondria" &
      Class   != "Chloroplast" &
      Order   != "Chloroplast" &
      Family  != "Mitochondria" 
  )

pseq <- ps_sub %>%
  tax_fix() %>%
  phyloseq_validate()

ord_explore(ps_free) 

# free-living dataframe
ps_free <- pseq %>% subset_samples(Filter_pores == "free-living") %>% prune_taxa(taxa_sums(.) > 0, .) 
ord_explore(ps_free) 
# particle-associated phyloseq
ps_part <- pseq %>% subset_samples(Filter_pores == "particle-associated") %>% prune_taxa(taxa_sums(.) > 0, .) 
ord_explore(ps_part) 
# iron line graph
iron_line <- ggplot(agg_iron, aes(x = factor(Station, level = level_order), y = mean_iron, group=as.character(More_Depth_Threshold))) + 
  facet_grid(~factor(More_Depth_Threshold, levels=c("Surface", "Mixed Layer", "Mid", "T-min", "Bottom"))~., axes="all", axis.labels = "margins") +# Connect points with lines
  geom_bar(stat="identity", fill="lightgray") +
  ylim(c(0,1.3)) +
  geom_point(pch=16, size=1.2, color="darkred")+
  geom_line(color="red2") +
  theme_classic2() +
  theme(text = element_text(family = "Helvetica"), 
        plot.title = element_text(hjust = 0.6),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=9, angle=90, vjust=0.5))+
  geom_vline(xintercept = c(2.5,6.5,12.5), linetype = "dashed", linewidth=0.6, color = "black")# Add vertical lines

iron_line <- ggplot(agg_iron, aes(x = factor(Station, level = level_order), y = Siderophore, group=as.character(More_Depth_Threshold))) + 
  facet_grid(~factor(More_Depth_Threshold, levels=c("Surface", "Mixed Layer", "Mid", "T-min", "Bottom"))~., axes="all", axis.labels = "margins") +# Connect points with lines
  geom_bar(stat="identity", fill="lightgray") +
  geom_point(pch=16, size=1.2, color="darkred")+
  geom_line(color="red2") +
  theme_classic2() +
  theme(text = element_text(family = "Helvetica"), 
        plot.title = element_text(hjust = 0.6),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=9, angle=90, vjust=0.5))+
  geom_vline(xintercept = c(2.5,6.5,12.5), linetype = "dashed", linewidth=0.6, color = "black")# Add vertical lines

iron_bar <- ggplot(agg_iron, aes(x = factor(Station, level = level_order), y = mean_iron, group=as.character(More_Depth_Threshold))) + 
  facet_grid(~factor(More_Depth_Threshold, levels=c("Surface", "Mixed Layer", "Mid", "T-min", "Bottom"))~., axes="all", axis.labels = "margins") +# Connect points with lines
  geom_bar(stat="identity") +
  ylim(c(0,1.3)) +
  theme_classic2() +
  theme(text = element_text(family = "Helvetica"), 
        plot.title = element_text(hjust = 0.6),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=9, angle=90, vjust=0.5))+
  geom_vline(xintercept = c(2.5,6.5,12.5), linetype = "dashed", linewidth=0.6, color = "black")# Add vertical lines


# free relabund

ps_agg <- ps_free %>% transform_sample_counts(function(x) {x/sum(x)}) %>% ps_melt() %>% filter(Abundance > 0.03)

ps_agg <- ps_agg %>% aggregate(Abundance ~ Station * OTU * More_Depth_Threshold, data = ., FUN = mean)

myPal <- tax_palette(
  data = ps_free, rank = "unique", n = 41, pal = "brewerPlus",
  add = c(Other = "white")
)

myPal["Seq215"] <- "purple"
myPal["Seq169"] <- "darkred"

tax_palette_plot(myPal)

ps_agg$More_Depth_Threshold[ps_agg$More_Depth_Threshold == "Mid-Bottom"] <- "T-min"
ps_agg$More_Depth_Threshold[ps_agg$More_Depth_Threshold == "Mid-Surface"] <- "Mixed Layer"

barplot_free <- ggplot(ps_agg, aes(x = factor(Station, level = level_order), y = Abundance, fill = OTU, group = OTU)) + facet_grid(~factor(More_Depth_Threshold, levels=c("Surface", "Mixed Layer", "Mid", "T-min", "Bottom"))~.) + # facet grid seperates by different levels, horizontally
  geom_bar(stat = "identity", position="fill", color= "black", linewidth=0.3, width=0.9) + theme_classic() + # adds black outline to boxes
  scale_y_continuous(expand = c(0, 0)) + # extends the barplots to the axies
  scale_fill_manual(values = myPal, drop = FALSE) +
  scale_x_discrete(
    breaks = level_order, # setting breaks
    labels = level_order, # settting levels
    drop = FALSE
  ) +# set the colors with custom colors (myColors)
  theme(text = element_text(family = "Helvetica"), 
        plot.title = element_text(hjust = 0.6),
        axis.title.x = element_blank(),
        axis.text=element_text(size=7),
        axis.text.x = element_text(size=9, angle=90, vjust=0.5),
        axis.title.y = element_blank()
        # Reducing space between facets
  ) + # Optionally remove panel borders
  geom_vline(xintercept = c(2.5,6.5,12.5), linetype = "dashed", linewidth=0.6, color = "black") + # Add vertical lines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1))

tax_free <- tax_table(ps_free)
tax_free <- as.data.frame(ps_part@sam_data)

# part relabund
ps_agg <- ps_part %>% transform_sample_counts(function(x) {x/sum(x)}) %>% ps_melt() %>% filter(Abundance > 0.005)

ps_agg <- ps_agg %>% aggregate(Abundance ~ Station * Order * More_Depth_Threshold, data = ., FUN = mean) 

myPal <- tax_palette(
  data = ps_part, rank = "Order", n = 41, pal = "brewerPlus",
  add = c(Other = "white")
)

myPal["Seq215"] <- "purple"
myPal["Seq169"] <- "darkred"

tax_palette_plot(myPal)

ps_agg$More_Depth_Threshold[ps_agg$More_Depth_Threshold == "Mid-Bottom"] <- "T-min"
ps_agg$More_Depth_Threshold[ps_agg$More_Depth_Threshold == "Mid-Surface"] <- "Mixed Layer"

barplot_part <- ggplot(ps_agg, aes(x = factor(Station, level = level_order), y = Abundance, fill = Order, group = Order)) + facet_grid(~factor(More_Depth_Threshold, levels=c("Surface", "Mixed Layer", "Mid", "T-min", "Bottom"))~.) + # facet grid seperates by different levels, horizontally
  geom_bar(stat = "identity", position="fill", color= "black", linewidth=0.3, width=0.9) + theme_classic() + # adds black outline to boxes
  scale_y_continuous(expand = c(0, 0)) + # extends the barplots to the axies
  scale_fill_manual(values = myPal, drop = FALSE) +
  scale_x_discrete(
    breaks = level_order, # setting breaks
    labels = level_order, # settting levels
    drop = FALSE
  ) +# set the colors with custom colors (myColors)
  theme(text = element_text(family = "Helvetica"), 
        plot.title = element_text(hjust = 0.6),
        axis.title.x = element_blank(),
        axis.text=element_text(size=7),
        axis.text.x = element_text(size=9, angle=90, vjust=0.5),
        axis.title.y = element_blank()
        # Reducing space between facets
  ) + # Optionally remove panel borders
  geom_vline(xintercept = c(2.5,6.5,12.5), linetype = "dashed", linewidth=0.6, color = "black") + # Add vertical lines
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1))

tax_part <- tax_table(ps_free)
tax_part <- as.data.frame(tax_part)


#### TRUE_FLOW INFLOW OUTFLOW RDA
ps_sub %>%
  tax_transform(rank = "unique", trans = "hellinger") %>%
  ord_calc(
    constraints = c("Latitude", "Longitude", "Salinity", "Temperature", "CTD_Depth", "Lab_NO3", "Lab_PO4", "Lab_NO2", "Lab_NH4", "Sb_Oxygen", "Iron"),
    method = "auto"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "True_Flow", fill = "True_Flow",
    shape = "More_Depth_Threshold", alpha = 0.5,
    size = 2
  ) + 
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse(
    ggplot2::aes(colour = True_Flow)
  )

### 
ps_free %>%
tax_transform(rank = "unique", trans = "identity") %>%
  ord_calc(
    constraints = c("Iron", "Latitude", "Longitude", "Salinity", "Temperature", "CTD_Depth", "Lab_NO3", "Lab_PO4", "Lab_NO2", "Lab_NH4", "Sb_Oxygen"),
    method = "auto"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    plot_taxa = 1:8,
    colour = "watertype", fill = "watertype",
    shape = "More_Depth_Threshold", alpha = 0.5,
    size = 2
  ) + 
  scale_shape_girafe_filled()

ps_part %>%
  tax_transform(rank = "unique", trans = "identity") %>%
  ord_calc(
    constraints = c("Iron", "Latitude", "Longitude", "Salinity", "Temperature", "CTD_Depth", "Lab_NO3", "Lab_PO4", "Lab_NO2", "Lab_NH4", "Sb_Oxygen"),
    method = "auto"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    plot_taxa = 1:8,
    colour = "watertype", fill = "watertype",
    shape = "More_Depth_Threshold", alpha = 0.5,
    size = 2
  ) + 
  scale_shape_girafe_filled()
