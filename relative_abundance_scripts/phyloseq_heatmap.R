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


htmp <- ps_noncontam_prev05 %>%
  tax_transform("clr", zero_replace = "halfmin") %>%
  comp_heatmap(
    taxa = taxa, samples = samples, colors = heat_palette(sym = TRUE),
    sample_anno = sampleAnnotation(
      Location = anno_sample_cat("Location", legend_title = "Location in Polynya")
    )
  )