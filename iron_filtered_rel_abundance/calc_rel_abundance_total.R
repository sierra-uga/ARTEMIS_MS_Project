calculate_relative_abundance <- function(abundance_df, station, More_Depth_Threshold, selected_taxa_names) {
  # Initialize empty data frame to store debug information
  debug_df <- data.frame(
    Station = character(),
    More_Depth_Threshold = character(),
    Abundance = numeric(),
    data_source = character(),
    stringsAsFactors = FALSE
  )
  
  # Convert station and More_Depth_Threshold to factors and then create a unique identifier for each combination
  unique_combinations <- unique(paste(station, More_Depth_Threshold, sep = "_"))
  
  # Initialize vectors to store results
  total_abundance_other <- numeric(length = length(unique_combinations))
  total_abundance_selected <- numeric(length = length(unique_combinations))
  relative_abundance <- numeric(length = length(unique_combinations))
  
  # Loop over unique combinations of Station and More_Depth_Threshold
  for (i in seq_along(unique_combinations)) {
    comb <- unique_combinations[i]
    indices <- which(paste(station, More_Depth_Threshold, sep = "_") == comb)
    selected_indices <- colnames(abundance_df) %in% selected_taxa_names
    
    # Calculate total abundance of other taxa
    total_abundance_other[i] <- sum(rowSums(abundance_df[indices, !selected_indices, drop = FALSE]))
    
    # Calculate total abundance of selected taxa
    total_abundance_selected[i] <- sum(rowSums(abundance_df[indices, selected_indices, drop = FALSE]))
    
    # Calculate relative abundance (percentage) for the current combination
    relative_abundance[i] <- total_abundance_selected[i] / (total_abundance_other[i] + total_abundance_selected[i]) * 100
    
    # Split combination into station and More_Depth_Threshold
    split_comb <- strsplit(comb, "_")[[1]]
    
    # Debug information data rows
    debug_df <- rbind(debug_df, data.frame(
      Station = split_comb[1],
      More_Depth_Threshold = split_comb[2],
      Abundance = total_abundance_other[i],
      data_source = "other",
      stringsAsFactors = FALSE
    ))
    
    debug_df <- rbind(debug_df, data.frame(
      Station = split_comb[1],
      More_Depth_Threshold = split_comb[2],
      Abundance = total_abundance_selected[i],
      data_source = "selected",
      stringsAsFactors = FALSE
    ))
  }
  
  # Create data frame with results
  result_df <- data.frame(
    Combination = unique_combinations,
    Relative_Abundance = relative_abundance
  )
  
  # Return both the result data frame and the debug data frame
  list(Result = result_df, Debug = debug_df)
}
