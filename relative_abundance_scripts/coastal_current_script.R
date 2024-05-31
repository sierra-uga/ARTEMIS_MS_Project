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

test <- subset_samples(pseq, Transect_Name == "transect1") 
coastal_list <- c(unique(test@sam_data[["Station"]]), "STN20")
depth_threshold <- 250

pseq1 <- pseq %>% subset_samples(Station %in% coastal_list)

numeric_list <- c("Lab_NO3", "Lab_NO2", "Lab_NH4", "Longitude", "Lab_PO4", "Salinity", "Temperature", "CTD_Depth", "Latitude")
for (col in numeric_list) {
  pseq1@sam_data[[col]] <- as.numeric(pseq1@sam_data[[col]])
}

### NOTES: maybe add outflow to filtering as well??, then plot?
pseq1 <- pseq1 %>% subset_samples(More_Depth_Threshold %in% c("Mid-Bottom", "Bottom"))

ps_free <- pseq1 %>% subset_samples(Filter_pores == "free-living") %>% prune_taxa(taxa_sums(.) > 0, .) 
ord_explore(ps_part)


# particle-associated phyloseq
ps_part <- pseq1 %>% subset_samples(Filter_pores == "particle-associated") %>% prune_taxa(taxa_sums(.) > 0, .) 

ps_agg_free <- ps_free %>% transform_sample_counts(function(x) {x/sum(x)}) %>% ps_melt() %>% filter(Abundance > 0.03)

ps_agg_free <- ps_agg_free %>% aggregate(Abundance ~ Station * OTU * More_Depth_Threshold, data = ., FUN = mean)



myPal <- tax_palette(
  data = ps_free, rank = "unique", n = 41, pal = "brewerPlus",
  add = c(Other = "white")
)

myPal["Seq215"] <- "purple"
myPal["Seq169"] <- "darkred"


## PLOTS

your_phyloseq %>%
  tax_transform(rank = "unique", trans = "identity") %>%
  dist_calc(dist = "bray") %>%
  ord_calc(
    method = "auto"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    colour = "Location", fill = "Location",
    shape = "Depth_Threshold", alpha = 0.5,
    size = 2
  ) + 
  scale_shape_girafe_filled() +
  ggplot2::stat_ellipse(
    ggplot2::aes(colour = Location)
  )

ps_free %>%
  tax_transform(rank = "unique", trans = "identity") %>%
  ord_calc(
    method = "auto"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    plot_taxa = 1:5,
    colour = "Location", fill = "Location",
    shape = "Location", alpha = 0.5,
    size = 2
  ) + 
  scale_shape_girafe_filled()