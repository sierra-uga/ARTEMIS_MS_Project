install.packages("MiscMetabar")
BiocManager::install("dada2", version = "3.19")
library(dada2)
library(MiscMetabar)

pseq@tax_table <- phyloseq::tax_table(cbind(
  pseq@tax_table,
  "taxon" = taxa_names(pseq)
))

res_height <- ancombc_pq(pseq, fact="Depth_Threshold", levels_fact = c("Surface", "Intermediate", "Intermediate"), tax_level = "unique")


df <- res_height$res
ggplot(
  res_height$res,
  aes(
    y = reorder(taxon, lfc_Depth_ThresholdIntermediate),
    x = lfc_Depth_ThresholdIntermediate,
    color = diff_Depth_ThresholdIntermediate
  )
) +
  geom_vline(xintercept = 0) +
  geom_segment(aes(
    xend = 0, y = reorder(taxon, lfc_Depth_ThresholdIntermediate),
    yend = reorder(taxon, lfc_Depth_ThresholdIntermediate)
  ), color = "darkgrey") +
  geom_point()