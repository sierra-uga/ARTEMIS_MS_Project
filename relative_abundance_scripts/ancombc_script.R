install.packages("MiscMetabar")
BiocManager::install("dada2", version = "3.19")
library(dada2)
library(MiscMetabar)

new_meta <- as.data.frame(sample_data(ps_noncontam_prev05))
new_meta <- new_meta %>%
  mutate(Iron_Level = ifelse(Iron >= 0.5, "High_iron", "Low_iron"))

ps <- ps_noncontam_prev05
#ps <- subset_samples(ps_noncontam_prev05, Pos_in_polynya %in% c("Eastern_CC", "Western_CC", "Dotson_CC", "Getz_CC", "Inflow", "Outflow"))
ps <- ps %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "Mitochondria" &
      Class   != "Chloroplast" &
      Order   != "Chloroplast" &
      Family  != "Mitochondria" 
  )
ps_sub <- subset_samples(ps, Sample.Control == "True.Sample") %>% subset_samples(., Iron_Level != "NA")
ps_sub <- ps_sub %>% tax_fix() %>%
  phyloseq_validate()

tax_agg(ps = ps_sub, "Genus")

ps_ps <- ps_sub %>% 
  tax_fix(unknowns = c("hydrothermal vent metagenome", "hydrothermal vent metagenome Class", "marine metagenome", "marine metagenome Class", "uncultured", "uncultured deep-sea bacterium", "uncultured planctomycete"))


pseq <- ps_ps %>%
  phyloseq_validate() %>% subset_samples(., Pos_in_polynya %in% c("Outflow", "Western_CC")) %>% 

pseq@tax_table <- phyloseq::tax_table(cbind(
  pseq@tax_table,
  "taxon" = taxa_names(pseq)
))

physeqGenus = tax_glom(pseq, "Genus")

ps_free <- subset_samples(pseq, Filter_pores == "free-living")
ps_part <- subset_samples(pseq, Filter_pores == "particle-associated")


below_free_tax <- as.data.frame(tax_table(ps_free_below))
below_free_tax$taxon <- rownames(below_free_tax)


print(out$res)

res_height <- ancombc_pq(ps_free_below, fact="Location", tax_level="Genus", levels_fact = c("Open_polynya", "Dotson"))


tse = mia::makeTreeSummarizedExperimentFromPhyloseq(ps_free_below)
library(ANCOMBC)

sample_data(ps_free_below)$Location <- as.factor(sample_data(ps_free_below)$Location)

out = ancombc2(data = tse, assay_name = "counts",
               fix_formula = "Location", group = "Location")

tax_level = "taxon"

new <- signif_ancombc(out)
rownames(sigtab) <- sigtab$taxon

new = cbind(as(new, "data.frame"), as(tax_table(ps_free_below)[rownames(new), ], "matrix"))

sig_taxa <- new$taxon

rownames(ps_free@tax_table) <- ps_free@tax_table$Genus

ps_free2 <- tax_glom(ps_free, taxrank = "Genus")

ps_free_final <- subset_taxa(ps_free_below, rownames(ps_free_below@tax_table) %in% sig_taxa)



df <- res_height$res
ggplot(
  res_height$res,
  aes(
    y = reorder(taxon, lfc_Pos_in_polynyaWestern_CC),
    x = lfc_Pos_in_polynyaWestern_CC,
    color = lfc_Pos_in_polynyaWestern_CC
  )
) +
  geom_vline(xintercept = 0) +
  geom_segment(aes(
    xend = 0, y = reorder(taxon, lfc_Pos_in_polynyaWestern_CC),
    yend = reorder(taxon, lfc_Pos_in_polynyaWestern_CC)
  ), color = "darkgrey") +
  geom_point()