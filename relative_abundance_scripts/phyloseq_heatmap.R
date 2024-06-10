if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"), update = FALSE, force=TRUE)

install.packages(
  "microViz",
  repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
)
install.packages("ggtext") # for rotated labels on ord_plot() 
install.packages("ggraph") # for taxatree_plots()
install.packages("DT") # for tax_fix_interactive()
install.packages("corncob") # for beta binomial models in tax_model()
library(microViz)

ps_sub <- subset_samples(ps_noncontam_prev05, Sample.Control == "True.Sample")

ps_sub <- ps_sub %>%
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

pseq <- pseq %>% tax_fix(unknowns = c("hydrothermal vent metagenome", "hydrothermal vent metagenome Class", "marine metagenome", "marine metagenome Class", "uncultured", "uncultured deep-sea bacterium", "uncultured planctomycete"))
pseq <- subset_samples(pseq, sample_name != "STN115.35.fil.dura.r2")

ps_free <- subset_samples(pseq, Filter_pores == "free-living")

ps_part <- subset_samples(pseq, Filter_pores == "particle-associated")
ps_part <- subset_samples(ps_part, sample_name != "STN198.20.pre.poly.3.S")


ord_explore(ps_part) # gif generated with microViz version 0.7.4 (plays at 1.75x speed)


pseq %>%
  tax_transform(rank = "Order", trans = "clr") %>%
  ord_calc(
    method = "auto"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    plot_taxa = 1:23,
    colour = "More_Depth_Threshold", fill = "More_Depth_Threshold",
    shape = "Filter_pores", alpha = 0.5,
    size = 2
  ) + 
  scale_shape_girafe_filled()

pseq %>%
  tax_transform(rank = "unique", trans = "identity") %>%
  ord_calc(
    method = "auto"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    plot_taxa = 1:9,
    colour = "Iron", fill = "Iron",
    shape = "Depth_Threshold", alpha = 0.5,
    size = 2
  ) + 
  scale_shape_girafe_filled()

pseq <- pseq %>% tax_fix(unknowns = c("metagenome Family", "hydrothermal vent metagenome Family", "marine metagenome Family", "metagenome Family", "seawater metagenome", "uncultured Family", "uncultured Rickettsiales bacterium", "uncultured bacterium", "uncultured bacterium Family", "uncultured deep-sea bacterium Family", "uncultured marine bacterium", "uncultured marine bacterium Family", "uncultured marine microorganism", "uncultured marine microorganism Family", "uncultured organism", "uncultured organism Family", "uncultured proteobacterium", "unidentified marine bacterioplankton"))


pseq <- pseq %>%
  tax_fix(
    min_length = 7,
    unknowns = c("metagenome Family"),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "classified"
  )

htmp <- pseq %>%
  tax_transform(trans = "clr", rank = "unique", zero_replace = "halfmin") %>%
  tax_filter(min_prevalence = 0.5, use_counts = FALSE) %>%
  comp_heatmap(
    colors = heat_palette(sym = TRUE), grid_col = "darkgrey",
    tax_anno = taxAnnotation(
      Prev. = anno_tax_prev(bar_width = 0.3, size = grid::unit(1, "cm"))
    ),
    sample_side = "top", name = "Robust\nCLR",
    sample_anno = sampleAnnotation(
      "Location" = anno_sample("Location"),
      "Depth_Threshold" = anno_sample("Depth_Threshold"),
      col = list("Depth_Threshold" = c(
        "Bottom_water" = "black", "Intermediate" = "orange", "Surface" = "lightgrey"
      ))
    )
  )


