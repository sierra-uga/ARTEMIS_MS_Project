agg_by_tax <- function(transect_number, main_data_frame) {
  # Subset data for the specified transect number
  transect_data <- subset(main_data_frame, (Transect_Number == as.character(transect_number)))
  
  # Check if any rows exist for the specified transect number
  if (nrow(transect_data) == 0) {
    warning("No data found for Transect ", transect_number)
    return(NULL)
  }
  
  # Create a temporary dataframe with relevant columns
  transect_data_temp <- as.data.frame(transect_data[c("Abundance", "TAX")])
  #transect_data_temp$TAX <- rownames(transect_data_temp)
  
  # Aggregate abundance by TAX
  transect_data_temp <- aggregate(. ~ TAX, data = transect_data_temp, FUN = mean)
  colnames(transect_data_temp) <- c("TAX", "Abundance")
  
  # Add Transect_Number column
  transect_data_temp$Transect_Number <- as.character(transect_number)
  
  # Calculate relative abundance
  transect_data_temp <- transect_data_temp %>%
    mutate(Relative = Abundance / sum(Abundance))
  
  return(transect_data_temp)
}

agg_by_tax_coastal <- function(transect_number, main_data_frame) {
  # Subset data for the specified transect number
  transect_data <- subset(main_data_frame, (Coastal_Current_Number == as.character(transect_number)))
  
  # Check if any rows exist for the specified transect number
  if (nrow(transect_data) == 0) {
    warning("No data found for Transect ", transect_number)
    return(NULL)
  }
  
  # Create a temporary dataframe with relevant columns
  transect_data_temp <- as.data.frame(transect_data[c("Abundance", "TAX")])
  #transect_data_temp$TAX <- rownames(transect_data_temp)
  
  # Aggregate abundance by TAX
  transect_data_temp <- aggregate(. ~ TAX, data = transect_data_temp, FUN = mean)
  colnames(transect_data_temp) <- c("TAX", "Abundance")
  
  # Add Coastal_Current_Number column
  transect_data_temp$Coastal_Current_Number <- as.character(transect_number)
  
  # Calculate relative abundance
  transect_data_temp <- transect_data_temp %>%
    mutate(Relative = Abundance / sum(Abundance))
  
  return(transect_data_temp)
}
