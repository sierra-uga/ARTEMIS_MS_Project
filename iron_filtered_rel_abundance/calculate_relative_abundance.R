calculate_relative_abundance <- function(abundance_df, transect_numbers, selected_taxa_names) {
  # Initialize empty data frame to store debug information
  debug_df <- data.frame(
    Transect = character(),
    Abundance = numeric(),
    data_source = character(),
    stringsAsFactors = FALSE
  )
  
  # Initialize vectors to store results
  total_abundance_other <- numeric(length = length(unique(transect_numbers)))
  total_abundance_selected <- numeric(length = length(unique(transect_numbers)))
  
  # Loop over unique transect numbers
  for (i in 1:length(unique(transect_numbers))) {
    transect <- unique(transect_numbers)[i]
    indices <- which(transect_numbers == transect)
    selected_indices <- colnames(abundance_df) %in% selected_taxa_names
    
    # Calculate total abundance of other taxa
    total_abundance_other[i] <- sum(rowSums(abundance_df[indices, !selected_indices, drop = FALSE]))
    
    # Calculate total abundance of selected taxa
    total_abundance_selected[i] <- sum(rowSums(abundance_df[indices, selected_indices, drop = FALSE]))
    
    # Determine data source based on column name
    if (length(selected_taxa_names) > 0) {
      # Create a row of debug information
      debug_row_other <- data.frame(
        Transect = as.character(transect),
        Abundance = total_abundance_other[i],
        data_source = "other",
        stringsAsFactors = FALSE
      )
      
      debug_row_selected <- data.frame(
        Transect = as.character(transect),
        Abundance = total_abundance_selected[i],
        data_source = "selected",
        stringsAsFactors = FALSE
      )
      
      # Append the debug rows to the debug data frame
      debug_df <- rbind(debug_df, debug_row_other, debug_row_selected)
    }
  }
  
  # Calculate relative abundance (percentage) of selected taxa
  relative_abundance <- total_abundance_selected / (total_abundance_other + total_abundance_selected) * 100
  
  # Create data frame with results
  result_df <- data.frame(
    Transect = as.character(unique(transect_numbers)),
    Relative_Abundance = relative_abundance
  )
  
  # Return both the result data frame and the debug data frame
  list(Result = result_df, Debug = debug_df)
}