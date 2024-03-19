filtered_df <- KO_iron_abundance_type %>%
  group_by(Type) %>% select(-sum)

gathered_df <- filtered_df %>%
  gather(key = "sample_name", value = "Abundance", -Type)

merged_df <- gathered_df %>%
  left_join(select(coastal_metadata_unique, sample_name, Filter_pores, Coastal_Current_Number), by = "sample_name")

# If you want to reorder columns
merged_df <- merged_df[, c("sample_name", "Type", "Abundance", "Filter_pores", "Coastal_Current_Number")]
merged_df$Coastal_Current_Number <- factor(merged_df$Coastal_Current_Number, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9"), ordered = TRUE) # order transect numbers for the boxplots

Fe_S <- subset(merged_df, merged_df$Type == "Fe-S")
Fe_S_plot <- ggplot(Fe_S, aes(x = Coastal_Current_Number, y = Abundance, fill = Filter_pores)) + 
  geom_bar(stat="identity", width=0.65) + ggtitle("Fe-S pathways") + add_color + facet_wrap(vars(Filter_pores)) +
  theme_classic2()

Siderophore_uptake <- subset(merged_df, merged_df$Type == "Siderophore uptake")
Sid_uptake_plot <- ggplot(Siderophore_uptake, aes(x = Coastal_Current_Number, y = Abundance, fill = Filter_pores)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2) + ggtitle("Siderophore Uptake") +
  ylim(0,100000)