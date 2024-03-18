### code for phyloseq in other script ###

ps_sub <- ps %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "Mitochondria" &
      Class   != "Chloroplast" &
      Order   != "Chloroplast" &
      Family  != "Mitochondria" 
  )

# free-living phyloseq
outflow_ps_free <- ps_sub %>% subset_samples(Filter_pores == "0.2") %>% prune_taxa(taxa_sums(.) > 0, .) 

# particle-associated phyloseq
outflow_ps_part <- ps_sub %>% subset_samples(Filter_pores >= "2") %>% prune_taxa(taxa_sums(.) > 0, .) 

colors <- c("red3", "lightblue2", "seagreen", "yellow2", "darkblue", "purple", "violet", "dodgerblue", "darkgrey", "darkgoldenrod", "black", "grey", "darkred", "darkgrey", "green", "dodgerblue4", "darkgoldenrod1", "pink", "green3")

#### 
# Create a data frame for freeliving
outflow_data_free <- outflow_ps_free %>%  
  subset_samples(Station != "STN068") %>% #remove station 68 from outflow
  subset_samples(True_Flow == "Outflow") %>%
  tax_glom(taxrank = "Genus") %>% # agglomerate at genus level
  transform_sample_counts(function(x) {x/sum(x)} )# Transform to rel. abundance

outflow_top_free <- top_taxa(outflow_data_free, 
                               n_taxa = 10,
                               include_na_taxa = F)

outflow_data_free <- outflow_top_free$ps_obj %>%
  psmelt() %>%  
  #filter(., Station == "transect1") %>%
  arrange(Station)   

# Plot 
outflow_barplot_part <- ggplot(outflow_data_free, aes(x = Station, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position="fill", width=1) + theme_classic() +
  # geom_col(position = "dodge") + # changes to multiple bars
  scale_fill_manual(values = colors) + # set manual colors
  theme(plot.title = element_text(hjust = 0.5, size=14)) +
  theme(axis.title.x = element_text()) + # remove x title
  theme(axis.text.y = element_text()) + # remove y text
  theme(axis.title.y = element_text()) + # remove y title
  theme(legend.position = "right") + # position legent
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  theme(panel.spacing.y = unit(1, "lines")) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Top 10 Genus) \n")


######## INFLOW ########### 
inflow_data_free <- outflow_ps_free %>%  
  subset_samples(Station != "STN068") %>% #remove station 68 from inflow
  subset_samples(True_Flow == "Inflow") %>%
  tax_glom(taxrank = "Genus") %>% # agglomerate at genus level
  transform_sample_counts(function(x) {x/sum(x)} )# Transform to rel. abundance

inflow_top_free <- top_taxa(inflow_data_free, 
                             n_taxa = 10,
                             include_na_taxa = F)

inflow_data_free <- inflow_top_free$ps_obj %>%
  psmelt() %>%  
  #filter(., Station == "transect1") %>%
  arrange(Station)   

# Plot 
inflow_barplot_part <- ggplot(inflow_data_free, aes(x = Station, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position="fill", width=1) + theme_classic() +
  # geom_col(position = "dodge") + # changes to multiple bars
  scale_fill_manual(values = colors) + # set manual colors
  theme(plot.title = element_text(hjust = 0.5, size=14)) +
  theme(axis.title.x = element_text()) + # remove x title
  theme(axis.text.y = element_text()) + # remove y text
  theme(axis.title.y = element_text()) + # remove y title
  theme(legend.position = "right") + # position legent
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  theme(panel.spacing.y = unit(1, "lines")) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Top 10 Genus) \n")
