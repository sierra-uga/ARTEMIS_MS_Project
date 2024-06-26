library(ComplexHeatmap)
library(tibble)

# Import Aldex2 results and rename the X variable by OTU, all three experiments
aldex2_ec_result <- read_excel("data/EC_data.xlsx", sheet = "ec_matrix") # 47 observations (EC) and 7 variables (EC description)
aldex2_ec_result <- as.data.frame(aldex2_ec_result)
rownames(aldex2_ec_result) <- aldex2_ec_result[, 1]
aldex2_ec_result <- aldex2_ec_result[, -(1)]


install.packages("zCompositions")
install.packages("seriation")
library("seriation")
library(zCompositions)
# Adjusting zeros on the matrix, all three experiments
aldex2_ec_result_czm <- cmultRepl(t(aldex2_ec_result),  label=0, method="CZM")
aldex2_ec_result_czm_trans <- t(apply(aldex2_ec_result_czm, 1, function(x){log(x) - mean(log(x))}))

############ Checking first lines of object, transpose it, and then create a heatmap according to the tax rank
head(aldex2_ec_result_czm_trans)
aldex2_ec_result_czm_trans_trav <- t(aldex2_ec_result_czm_trans)
aldex2_ec_result_czm_trans[, order(colnames(aldex2_ec_result_czm_trans))]
heatmap(aldex2_ec_result_czm_trans_trav, scale = "none", col = bluered(100),
        Colv = NA)


#Clean up presentation, all three experiments
ec_info <- data.frame(tax_table(ps_subset))
ec_info <- ec_info %>% rownames_to_column(var = "KO_Num")
ec_info <- dplyr::inner_join(ec_info, iron_KO)

ec_count <- data.frame(otu_table(ps_subset))
ec_count <- ec_count %>% rownames_to_column(var = "KO_Num")

ec_sample_tab <- data.frame(sample_data(ps_subset))
#ec_sample_tab <- ec_sample_tab %>% rownames_to_column(var = "sample_name")

# Define palette color
col_matrix <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
col_matrix2 <- colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
col_matrix3 <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Combine the heatmap and the annotation
aldex2_ec_result_czm_trans_trav_group <- scale(aldex2_ec_result_czm_trans_trav)

# Define colors for each level of qualitative variables, i.e. soil type and habitat
# Create the heatmap annotation
ha_top_ec = HeatmapAnnotation(border = TRUE,
                              Depth = as.vector(ec_sample_tab$Depth_Threshold),
                              Sample = as.vector(ec_sample_tab$Pos_in_polynya),
                              col = list(
                                Depth = c("Surface" = "lightgreen", "Intermediate" = "purple", "Bottom_water" = "darkred"),
                                Sample = c(c("Eastern_CC" = "darkblue", "Western_CC" = "lightblue", "Open_polynya_surf" = "green", "Dotson_CC" = "blue", 
                                             "Inflow" = "pink", "Outflow"="red"))
                              ),
                              annotation_legend_param = list(
                                Depth = list(
                                  title = "Depth",
                                  at = c("Surface", "Intermediate", "Bottom_water"),
                                  labels = c ("Surface", "Intermediate", "Bottom")
                                ),
                                Sample = list(
                                  title = "Location",
                                  at = c("Eastern_CC", "Western_CC", "Open_polynya_surf", "Dotson_CC", "Inflow", "Outflow"),
                                  labels = c("East CC", "West CC", "Open polynya", "Dotson CC", "Inflow", "Outflow")
                                )
                              ))

# Organize total abundance and taxa name
aldex2_ec_result_name <- as.data.frame(aldex2_ec_result)
str(aldex2_ec_result_name)
aldex2_ec_result_name$KO_Num <- rownames(aldex2_ec_result_name)
aldex2_ec_result_name <- left_join(aldex2_ec_result_name, ec_info, by="KO_Num")
head(aldex2_ec_result_name)
aldex2_ec_result_name <- aldex2_ec_result_name[, -(1:137)]

ec_count_total <- as.data.frame(ec_count)
ec_count_total$Total <- rowSums(ec_count[, -1])
#ec_count_total$KO_Num <- rownames(ec_count_total)
head(ec_count_total)
ec_count_total <- ec_count_total[, -(2:138)]
head(ec_count_total)

filter <- row.names(aldex2_ec_result_czm_trans_trav_group)

rownames(aldex2_ec_result_name) <- aldex2_ec_result_name$KO_Num
z_score_aldex2_ec_result_name_count_total <- left_join(aldex2_ec_result_name, ec_count_total)

z_score_aldex2_ec_result_name_count_total <- filter(z_score_aldex2_ec_result_name_count_total, KO_Num %in% filter)
head(z_score_aldex2_ec_result_name_count_total)


# Defining color scheme for row annotations

install.packages('circlize')
library("circlize")
col_matrix3 <- colorRamp2(c(-2.5, 0, 1.5), c("blue", "white", "red"))
abundance_col_fun = colorRamp2(c(0, 300, 3000, 30000, 300000),
                               c("#f5f5f5",
                                 "#c7eae5",
                                 "#80cdc1",
                                 "#35978f",
                                 "#01665e"))
process_col_fun = colorRamp2(c("Carbon", "Nitrogen", "Sulfur"),
                             c("#af8dc3", "#f7f7f7", "#7fbf7b"))

process_col_discr = gpar(c("Carbon", "Nitrogen", "Sulfur"),
                         c("#af8dc3", "#f7f7f7", "#7fbf7b"))
colors = structure(1:9, names = unique(z_score_aldex2_ec_result_name_count_total$Pathway.Final.Annotation.used.in.Fig.3.5))

ha_right_ec1 = rowAnnotation(
  Abundance = z_score_aldex2_ec_result_name_count_total$Total, border = TRUE, col = list(Abundance = abundance_col_fun))

ha_right_ec2 = rowAnnotation(
  Process = z_score_aldex2_ec_result_name_count_total$Pathway.Final.Annotation.used.in.Fig.3.5, border = TRUE, col  = list(Process = colors))
row_labels_ec = z_score_aldex2_ec_result_name_count_total$Name
row_labels_ec_lab = z_score_aldex2_ec_result_name_count_total$Pathway.Final.Annotation.used.in.Fig.3.5
#rownames(z_score_aldex2_ec_result_name_count_total) <- z_score_aldex2_ec_result_name_count_total$Name

ec_sample_tab <- ec_sample_tab %>% arrange(match(Pos_in_polynya, col_order))

col_order <- c("Inflow", "Outflow", "Open_polynya_surf", "Western_CC", "Dotson_CC", "Eastern_CC")

hm_ec_sel <- Heatmap(aldex2_ec_result_czm_trans_trav_group, name = "Z-score, CLR", col = col_matrix3,
                     #column_title = "Experiments & Sample Types", 
                     #column_title_gp = gpar(fontface = "bold", fontsize = 14),
                     #cluster_columns = cluster_within_group(as.vector(ec_sample_tab$SampleType)),
                     column_split = as.vector(ec_sample_tab$Filter_pores),
                     border = TRUE,
                     top_annotation = ha_top_ec,
                     right_annotation = ha_right_ec1,
                     #row_title = "",
                     row_labels = row_labels_ec,
                     row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                     row_names_side="left",
                     row_names_gp = gpar(fontsize = 6),
                     column_names_gp = gpar(fontsize = 8),
                     row_order = order(row_labels_ec_lab),
                     #column_order = order(ec_sample_tab$Pos_in_polynya),
                     rect_gp = gpar(col = "white", lwd = 1),
                     show_column_names = FALSE,
                     show_row_names = TRUE,
                     show_heatmap_legend = TRUE) + ha_right_ec2

draw(hm_ec_sel, padding = unit(c(1, 30, 1, 1), "mm"))
hm_ec_sel




##TAXA

row_names <- rownames(sigtab)

# Use gsub to extract the "SeqXXX" part from each row name
wanted_samples <- gsub(".*\\.(Seq[0-9]+)$", "\\1", row_names)

# Filter out any entries that do not match the pattern
wanted_samples <- wanted_samples[grepl("Seq[0-9]+$", wanted_samples)]

wanted_samples <- unique(wanted_samples)

taxonomy_table <- tax_table(ps_ps)
subset_taxonomy_table <- taxonomy_table[wanted_samples, , drop = FALSE]

ps_subset <- phyloseq(otu_table(ps_ps), sample_data(ps_ps), subset_taxonomy_table, phy_tree(ps_ps))

ps_subset <- subset_samples(ps_subset, Pos_in_polynya != "Cont_Shelf") %>% subset_samples(., sample_name != "STN089.200.fil.dura.r1") %>% subset_samples(., sample_name != "STN089.200.fil.dura.r2")

ps_subset <- ps_subset %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "Mitochondria" &
      Class   != "Chloroplast" &
      Order   != "Chloroplast" &
      Family  != "Mitochondria" 
  )

otu_table <- data.frame(otu_table(ps_free_final))

aldex2_ec_result <- otu_table

install.packages("zCompositions")
install.packages("seriation")
library("seriation")
library(zCompositions)
# Adjusting zeros on the matrix, all three experiments
aldex2_ec_result_czm <- cmultRepl(t(aldex2_ec_result),  label=0, method="CZM")
aldex2_ec_result_czm <- t(aldex2_ec_result)
aldex2_ec_result_czm_trans <- t(apply(aldex2_ec_result_czm, 1, function(x){log(x) - mean(log(x))}))

############ Checking first lines of object, transpose it, and then create a heatmap according to the tax rank
head(aldex2_ec_result_czm_trans)
aldex2_ec_result_czm_trans_trav <- t(aldex2_ec_result_czm_trans)
aldex2_ec_result_czm_trans[, order(colnames(aldex2_ec_result_czm_trans))]
heatmap(aldex2_ec_result_czm_trans_trav, scale = "none", col = bluered(100),
        Colv = NA)



#Clean up presentation, all three experiments
ec_info <- data.frame(tax_table(ps_free_final))
ec_info <- ec_info %>% rownames_to_column(var = "Seq_Num")
#ec_info <- dplyr::inner_join(ec_info, iron_KO)

ec_count <- data.frame(otu_table(ps_free_final))
ec_count <- ec_count %>% rownames_to_column(var = "Seq_Num")

ec_sample_tab <- data.frame(sample_data(ps_free_final))
ec_sample_tab <- ec_sample_tab %>% rownames_to_column(var = "Seq_Num")

# Define palette color
col_matrix <- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
col_matrix2 <- colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
col_matrix3 <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Combine the heatmap and the annotation
aldex2_ec_result_czm_trans_trav_group <- scale(aldex2_ec_result_czm_trans_trav)

# Define colors for each level of qualitative variables, i.e. soil type and habitat
# Create the heatmap annotation
ha_top_ec = HeatmapAnnotation(border = TRUE,
                              Depth = as.vector(ec_sample_tab$More_Depth_Threshold),
                              Sample = as.vector(ec_sample_tab$Location),
                              Iron = as.vector(ec_sample_tab$Iron_Level),
                              col = list(
                                Depth = c("Surface" = "lightgreen", "Mid-Surface" = "darkgreen", "Mid" = "purple", "Mid-Bottom" = "darkorange", "Bottom" = "darkred"),
                                Sample = c(c("Eastern_CC" = "darkblue", "Western_CC" = "lightblue", "Open_polynya" = "green", "Dotson" = "blue", 
                                             "Inflow" = "pink", "Outflow"="red", "Getz" = "yellow", "Cont_Shelf" = "darkgray")),
                                Iron = c("High" = "brown", "Low" = "yellow")
                              ),
                              annotation_legend_param = list(
                                Depth = list(
                                  title = "Depth",
                                  at = c("Surface", "Mid-Surface", "Mid","Mid-Bottom", "Bottom"),
                                  labels = c ("Surface", "Mixed Layer", "Mid", "T-min", "Bottom")
                                ),
                                Sample = list(
                                  title = "Location",
                                  at = c("Eastern_CC", "Western_CC", "Open_polynya", "Dotson", "Inflow", "Outflow", "Getz", "Cont_Shelf"),
                                  labels = c("East CC", "West CC", "Open polynya", "Dotson", "Inflow", "Outflow", "Getz", "Cont. Shelf")
                                ),
                                Iron = list(
                                  title = "Iron",
                                  at = c("High", "Low"),
                                  labels = c("High dFe >= 0.5 nmol/kg", "Low dFe < 0.5 nmol/kg")
                                )
                              ))

# Organize total abundance and taxa name
aldex2_ec_result_name <- as.data.frame(aldex2_ec_result)
str(aldex2_ec_result_name)
aldex2_ec_result_name$Seq_Num <- rownames(aldex2_ec_result_name)
aldex2_ec_result_name <- left_join(aldex2_ec_result_name, ec_info, by="Seq_Num")
head(aldex2_ec_result_name)
aldex2_ec_result_name <- aldex2_ec_result_name[, -(1:62)]

ec_count_total <- as.data.frame(ec_count)
ec_count_total$Total <- rowSums(ec_count[, -1])
#ec_count_total$KO_Num <- rownames(ec_count_total)
head(ec_count_total)
ec_count_total <- ec_count_total[, -(2:63)]
head(ec_count_total)

filter <- row.names(aldex2_ec_result_czm_trans_trav_group)

#rownames(aldex2_ec_result_name) <- aldex2_ec_result_name$Seq_Num
z_score_aldex2_ec_result_name_count_total <- left_join(aldex2_ec_result_name, ec_count_total)

z_score_aldex2_ec_result_name_count_total <- filter(z_score_aldex2_ec_result_name_count_total, Seq_Num %in% filter)
head(z_score_aldex2_ec_result_name_count_total)


# Defining color scheme for row annotations

install.packages('circlize')
library("circlize")
col_matrix3 <- colorRamp2(c(-2, 0, 3), c("blue", "white", "red"))
abundance_col_fun = colorRamp2(c(0, 300, 3000, 30000, 300000),
                               c("#f5f5f5",
                                 "#c7eae5",
                                 "#80cdc1",
                                 "#35978f",
                                 "#01665e"))
process_col_fun = colorRamp2(c("Carbon", "Nitrogen", "Sulfur"),
                             c("#af8dc3", "#f7f7f7", "#7fbf7b"))

process_col_discr = gpar(c("Carbon", "Nitrogen", "Sulfur"),
                         c("#af8dc3", "#f7f7f7", "#7fbf7b"))

myColors <- c(brewer.pal(9, "Paired"),'#e66101','darkgreen','#fdb863','#5e3c99', '#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','darkred','#c7eae5','#80cdc1','#35978f','#01665e','#003c30', "#A43D27", "#497687", "#5E4987", "darkgoldenrod", "lightblue2", "darkblue", "dodgerblue", "seagreen", "purple", "black")

colors = structure(1:10, names = unique(z_score_aldex2_ec_result_name_count_total$Phylum))

top_ec1 = rowAnnotation(
  Abundance = z_score_aldex2_ec_result_name_count_total$Total, border = TRUE, col = list(Abundance = abundance_col_fun))

ha_right_ec1 = rowAnnotation(
  Abundance = z_score_aldex2_ec_result_name_count_total$Total, border = TRUE, col = list(Abundance = abundance_col_fun))

unique_phyla <- unique(z_score_aldex2_ec_result_name_count_total$Class)

if (length(unique_phyla) > length(myColors)) {
  stop("Not enough colors specified in myColors for the unique phylum names.")
}

# Create a named vector for phylum colors
phylum_colors <- setNames(myColors[1:length(unique_phyla)], unique_phyla)

# Create the row annotation
ha_right_ec2 <- rowAnnotation(
  Class = z_score_aldex2_ec_result_name_count_total$Class, 
  border = TRUE, 
  col = list(Class = phylum_colors)
)

ha_right_ec2 = rowAnnotation(
  Class = z_score_aldex2_ec_result_name_count_total$Class, border = TRUE, col = list(myColors = myColors))


row_labels_ec = z_score_aldex2_ec_result_name_count_total$Genus

#rownames(z_score_aldex2_ec_result_name_count_total) <- z_score_aldex2_ec_result_name_count_total$Name

ec_sample_tab <- ec_sample_tab %>% arrange(match(Pos_in_polynya, col_order))

col_order <- c("Inflow", "Outflow", "Open_polynya_surf", "Western_CC", "Dotson_CC", "Eastern_CC")

as.vector(ec_sample_tab$Station)

pdf("relative_abundance_scripts/graphics/iron_free_heatmap.pdf", width=15, height = 12)

hm_ec_sel <- Heatmap(aldex2_ec_result_czm_trans_trav_group, name = "Z-score, CLR", col = col_matrix3,
                     column_title = "Particle-associated Community for High/Low dFe", 
                     column_title_gp = gpar(fontface = "bold", fontsize = 14),
                     #cluster_columns = cluster_within_group(as.vector(ec_sample_tab$SampleType)),
                     column_split = as.vector(ec_sample_tab$Iron_Level),
                     border = TRUE,
                     top_annotation = ha_top_ec,
                     left_annotation = ha_right_ec2,
                     row_title = "Genus",
                     row_labels = row_labels_ec,
                     row_title_gp = gpar(fontsize = 8, fontface = "bold"),
                     row_names_side="left",
                     row_names_gp = gpar(fontsize = 7),
                     column_names_gp = gpar(fontsize = 7),
                     row_order = order(z_score_aldex2_ec_result_name_count_total$Class),
                     column_order = order(ec_sample_tab$Depth_Threshold),
                     rect_gp = gpar(col = "white", lwd = 1),
                     show_column_names = FALSE,
                     show_row_names = TRUE,
                     show_heatmap_legend = TRUE) + ha_right_ec1
dev.off()

draw(hm_ec_sel, padding = unit(c(1, 30, 1, 1), "mm"))
hm_ec_sel

Heatmap(aldex2_ec_result_czm_trans_trav_group, name = "Z-score, CLR")
