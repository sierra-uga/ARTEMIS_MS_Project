library(phyloseq)
library(DESeq2)
library(dplyr)
library(microViz)

# for TAXA
ps <- ps_noncontam_prev05
ps <- subset_samples(ps_noncontam_prev05, Pos_in_polynya %in% c("Eastern_CC", "Western_CC", "Dotson_CC", "Getz_CC", "Inflow", "Outflow"))
ps <- ps %>%
  subset_taxa(
    Kingdom == "Bacteria" &
      Family  != "Mitochondria" &
      Class   != "Chloroplast" &
      Order   != "Chloroplast" &
      Family  != "Mitochondria" 
  )
ps_sub <- subset_samples(ps, Sample.Control == "True.Sample")
ps_sub <- ps_sub %>% tax_fix() %>%
  phyloseq_validate()

tax_agg(ps = ps_sub, "Genus")

ps_ps <- ps_sub %>% 
  tax_fix(unknowns = c("hydrothermal vent metagenome", "hydrothermal vent metagenome Class", "marine metagenome", "marine metagenome Class", "uncultured", "uncultured deep-sea bacterium", "uncultured planctomycete"))

  
ps_ps@otu_table <- ps_ps@otu_table + 1 #adding pseudocount

wanted_samples <- ps_sub@sam_data[["sample_name"]]

POS_IN_POLYNYA <- c(
  "Eastern_CC, Western_CC",
  "Eastern_CC, Dotson_CC",
  "Eastern_CC, Getz_CC",
  "Western_CC, Dotson_CC",
  "Western_CC, Getz_CC",
  "Dotson_CC, Getz_CC"
)


POS_IN_POLYNYA <- c(
  "Eastern_CC, Western_CC",
  "Eastern_CC, Open_polynya_surf",
  "Eastern_CC, Dotson_CC",
  "Eastern_CC, Inflow",
  "Eastern_CC, Outflow",
  "Western_CC, Open_polynya_surf",
  "Western_CC, Dotson_CC",
  "Western_CC, Inflow",
  "Western_CC, Outflow",
  "Open_polynya_surf, Dotson_CC",
  "Open_polynya_surf, Inflow",
  "Open_polynya_surf, Outflow",
  "Dotson_CC, Inflow",
  "Dotson_CC, Outflow",
  "Inflow, Outflow")

# Define your function
combine_sigtabs <- function(pos_in_polynya_combinations, ps_ps) {
  
  # Initialize an empty list to store results
  combined_sigtabs <- list()
  
  # Loop through each combination of POS_IN_POLYNYA
  for (POS_IN_POLYNYA in pos_in_polynya_combinations) {
    
    cat("Processing combination:", POS_IN_POLYNYA, "\n")
    
    # Subset phyloseq object based on POS_IN_POLYNYA
    pseq2 <- ps_ps %>%
      phyloseq_validate() %>%
      subset_samples(Pos_in_polynya %in% c(strsplit(POS_IN_POLYNYA, ", ")[[1]]))
    
    # Subset samples for free-living and particle-associated
    ps_free <- subset_samples(pseq2, Filter_pores == "free-living")
    ps_part <- subset_samples(pseq2, Filter_pores == "particle-associated")
    
    # DESeq analysis for free-living
    obj_free <- phyloseq_to_deseq2(ps_free, ~ Depth + Pos_in_polynya)
    diagdds_free <- DESeq(obj_free, test = "Wald", fitType = "mean", full = design(obj_free))
    res_free <- results(diagdds_free, cooksCutoff = FALSE)
    sigtab_free <- res_free[which(res_free$padj < 0.01), ]
    sigtab_free <- cbind(as(sigtab_free, "data.frame"), as(tax_table(pseq2)[rownames(sigtab_free), ], "matrix"))
    
    # DESeq analysis for particle-associated
    obj_part <- phyloseq_to_deseq2(ps_part, ~ Pos_in_polynya)
    diagdds_part <- DESeq(obj_part, test = "Wald", fitType = "mean", full = design(obj_part))
    res_part <- results(diagdds_part, cooksCutoff = FALSE)
    sigtab_part <- res_part[which(res_part$padj < 0.01), ]
    sigtab_part <- cbind(as(sigtab_part, "data.frame"), as(tax_table(pseq2)[rownames(sigtab_part), ], "matrix"))
    
    # Combine sigtab_free and sigtab_part
    combined_sigtab <- bind_rows(sigtab_free, sigtab_part)
    
    # Store combined result in the list
    combined_sigtabs[[POS_IN_POLYNYA]] <- combined_sigtab
  }
  
  # Combine all results into a single dataframe
  final_combined <- do.call(rbind, combined_sigtabs)
  
  return(final_combined)
}

# Define your POS_IN_POLYNYA combinations
POS_IN_POLYNYA_combinations <- c(
  "Eastern_CC, Western_CC",
  "Eastern_CC, Dotson_CC",
  "Eastern_CC, Getz_CC",
  "Western_CC, Dotson_CC",
  "Western_CC, Getz_CC",
  "Dotson_CC, Getz_CC"
)


POS_IN_POLYNYA_combinations <- c(
  "Eastern_CC, Western_CC",
  "Eastern_CC, Open_polynya_surf",
  "Eastern_CC, Dotson_CC",
  "Eastern_CC, Inflow",
  "Eastern_CC, Outflow",
  "Western_CC, Open_polynya_surf",
  "Western_CC, Dotson_CC",
  "Western_CC, Inflow",
  "Western_CC, Outflow",
  "Open_polynya_surf, Dotson_CC",
  "Open_polynya_surf, Inflow",
  "Open_polynya_surf, Outflow",
  "Dotson_CC, Inflow",
  "Dotson_CC, Outflow",
  "Inflow, Outflow"
)

# Assuming pi_ps is your phyloseq object
# Call the function to get the combined results
combined_results <- combine_sigtabs(POS_IN_POLYNYA_combinations, ps_ps)

filtered_results <- combined_results %>%
  filter(Phylum != "Bacteria Kingdom")

KO_sig_numbers <- filtered_results$Genus # set vector for numbers
KO_sig_tab <- iron_KO[iron_KO$KO_Num %in% KO_sig_numbers, ] 

KO_sig_numbers <- KO_sig_tab$Genus

ps_part_final <- subset_samples(ps_ps, Pos_in_polynya != "Cont_Shelf") %>% subset_taxa(ps_ps, Genus %in% KO_sig_numbers) %>% tax_agg("Genus") 


library(phyloseq)
library(DESeq2)
library(dplyr)

# Define your function
combine_sigtabs <- function(pos_in_polynya_combinations, ps_ps) {
  
  # Initialize empty lists to store results
  sigtab_free_list <- list()
  sigtab_part_list <- list()
  
  # Loop through each combination of POS_IN_POLYNYA
  for (POS_IN_POLYNYA_comb in pos_in_polynya_combinations) {
    
    # Split the current combination into individual parts
    pos_in_polynya <- strsplit(POS_IN_POLYNYA_comb, ", ")[[1]]
    
    cat("Processing combination:", POS_IN_POLYNYA_comb, "\n")
    
    # Subset phyloseq object based on current POS_IN_POLYNYA combination parts
    sample_ids <- rownames(ps_ps@sam_data[apply(ps_ps@sam_data$Pos_in_polynya %in% pos_in_polynya, 1, all), ])
    pseq <- subset_samples(ps_ps, sample_name %in% sample_ids)
    
    # Check the number of samples in pseq
    cat("Number of samples in pseq:", sample_sums(pseq), "\n")
    
    # Subset samples for free-living and particle-associated
    ps_free <- subset_samples(pseq, Filter_pores == "free-living")
    ps_part <- subset_samples(pseq, Filter_pores == "particle-associated")
    
    # Check the number of samples in ps_free and ps_part
    cat("Number of samples in ps_free:", sample_sums(ps_free), "\n")
    cat("Number of samples in ps_part:", sample_sums(ps_part), "\n")
    
    # DESeq analysis for free-living
    obj_free <- phyloseq_to_deseq2(ps_free, ~ More_Depth_Threshold + Pos_in_polynya)
    diagdds_free <- DESeq(obj_free, test = "Wald", fitType = "mean", full = design(obj_free))
    res_free <- results(diagdds_free, cooksCutoff = FALSE)
    sigtab_free <- res_free[which(res_free$padj < 0.01), ]
    sigtab_free <- cbind(as(sigtab_free, "data.frame"), as(tax_table(pseq)[rownames(sigtab_free), ], "matrix"))
    
    # Store sigtab_free in the list based on POS_IN_POLYNYA_comb
    sigtab_free_list[[POS_IN_POLYNYA_comb]] <- sigtab_free
    
    # DESeq analysis for particle-associated
    obj_part <- phyloseq_to_deseq2(ps_part, ~ More_Depth_Threshold + Location)
    diagdds_part <- DESeq(obj_part, test = "Wald", fitType = "mean", full = design(obj_part))
    res_part <- results(diagdds_part, cooksCutoff = FALSE)
    sigtab_part <- res_part[which(res_part$padj < 0.01), ]
    sigtab_part <- cbind(as(sigtab_part, "data.frame"), as(tax_table(pseq)[rownames(sigtab_part), ], "matrix"))
    
    # Store sigtab_part in the list based on POS_IN_POLYNYA_comb
    sigtab_part_list[[POS_IN_POLYNYA_comb]] <- sigtab_part
    
    # Print message to verify each iteration
    cat("Completed processing:", POS_IN_POLYNYA_comb, "\n")
  }
  
  # Return lists of sigtab_free and sigtab_part
  return(list(sigtab_free = sigtab_free_list, sigtab_part = sigtab_part_list))
}

# Define your POS_IN_POLYNYA combinations
POS_IN_POLYNYA_combinations <- c(
  "Eastern_CC, Western_CC",
  "Eastern_CC, Open_polynya_surf",
  "Eastern_CC, Dotson_CC",
  "Eastern_CC, Inflow",
  "Eastern_CC, Outflow",
  "Western_CC, Open_polynya_surf",
  "Western_CC, Dotson_CC",
  "Western_CC, Inflow",
  "Western_CC, Outflow",
  "Open_polynya_surf, Dotson_CC",
  "Open_polynya_surf, Inflow",
  "Open_polynya_surf, Outflow",
  "Dotson_CC, Inflow",
  "Dotson_CC, Outflow",
  "Inflow, Outflow"
)

# Assuming pi_ps is your phyloseq object
# Call the function to get the separated results
separated_results <- combine_sigtabs(POS_IN_POLYNYA_combinations, ps_ps)

# Access sigtab_free and sigtab_part results separately
sigtab_free_results <- separated_results$sigtab_free
sigtab_part_results <- separated_results$sigtab_part