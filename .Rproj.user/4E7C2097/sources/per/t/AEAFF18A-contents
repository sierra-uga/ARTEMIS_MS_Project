library(readr)
library(ggpicrust2) # required for differential abundance + pathway annotation
library(tibble)
library(tidyverse)
library(ggpubr)
library(ggprism)
library(patchwork)
require(grid)
library(LinDA)
library(DESeq2)
#install.packages("IgAScores")
library(IgAScores) # gives relative abundance of count table

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

ASV.hel = as.data.frame(
  apply(ASV, 2, function(x) sqrt(x / sum(x)))) # hellringer transformation

ps_sub <- ps_iron

sample_data(ASV.hel) <- ps_sub

deseq <- phyloseq_to_deseq2(ps_sub, ~Iron)

diagdds = DESeq(deseq)  #, fitType='local')
res = results(diagdds)
res = res[order(res$padj, na.last = NA), ]
alpha = 0.01
keepOTUs = rownames(res[res$padj > alpha, ])[1:50]
kosticTrimvs = prune_taxa(keepOTUs, ps_sub)
kosticTrim0 = prune_taxa(keepOTUs, ps_sub)
plot_heatmap(kosticTrimvs, taxa.order = "Phylum", taxa.label = "Order", sample.label = "Iron", 
             sample.order = "Iron")

diagdds <- DESeq(diagdds, test="Wald", fitType="parametric")



deseq <- phyloseq_to_DESeq(ps_sub)

tse_genus <- tax_glom(ps_sub, taxrank = "Order")


ds2 <- DESeqDataSet(tse_genus, ~Location)


# Annotate pathway results without KO to KEGG conversion
daa_annotated_sub_method_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_sub_method_results_df, ko_to_kegg = TRUE)

# differential abundance  
linda_station_type_abundance_mean <- station_type_abundance_mean
rownames(linda_station_type_abundance_mean) <- linda_station_type_abundance_mean$Type
linda_station_type_abundance_mean$Type <- NULL
library(MicrobiomeStat)
linda.obj <- linda(linda_station_type_abundance_mean, inflow_metadata, formula = "~Transect_Number+Filter_pores",  feature.dat.type="proportion", p.adj.method = "BH",)
linda_plot <- linda.plot(linda.obj, c("Transect_Number"),
                         titles = c('Transect # AND Community Diff Abund.'), alpha = 0.1, lfc.cut = 1,
                         legend = TRUE, directory = NULL, width = 11, height = 8)

plot_print <- linda_plot[["plot.lfc"]]
ggsave("linda_inflow_Conf.pdf", width=10, height=6)

daa_results_df_transect <- pathway_daa(abundance = station_type_abundance_mean %>% column_to_rownames("Type"), metadata = inflow_metadata, group = "Transect_Number", daa_method = "LinDA", select = NULL, p.adjust = "BH", reference = "1")
daa_results_df_community <- pathway_daa(abundance = station_type_abundance_mean %>% column_to_rownames("Type"), metadata = inflow_metadata, group = "Filter_pores", daa_method = "LinDA", select = NULL, p.adjust = "BH", reference = "1")


# extend Deseq2 function?


# sources:
 # https://cran.r-project.org/web/packages/ggpicrust2/vignettes/using_ggpicrust2.html ## GGPICRUST2 deseq method?
 # https://www.r-bloggers.com/2016/09/deseq2-course-work/ -- deseq2 heatmap
 # https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html -- deseq2 explanation
 # https://microbiome.github.io/course_2021_radboud/differential-abundance-analysis.html - comparison of diff. abundance methods
