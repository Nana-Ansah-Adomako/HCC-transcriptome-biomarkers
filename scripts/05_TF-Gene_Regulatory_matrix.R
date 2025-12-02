# ----------------------------------------------------------------------
# 1. Data Preparation
# ----------------------------------------------------------------------

# Input data: genes and their TF lists
new_input_data <- data.frame(
  Gene = c("BUB1B", "CCNB1", "CDC20", "MCM4", "RRM2", "TPX2", "UBE2C", "UBE2S", 
           "CCNB2", "ECT2", "KPNA2", "PTTG1", "SAE1", "TTK", 
           "CDKN3", "NEK2", "PBK", 
           "ANXA2", "PKM", 
           "PHGDH", 
           "G6PD"),
  Set_List_String = c(
    "FOXM1, CENPA, MYBL2, TFDP1, DNMT1, HMGA1, E2F7, PA2G4, ZNF367, E2F1",
    "FOXM1, CENPA, MYBL2, TFDP1, DNMT1, HMGA1, E2F7, PA2G4, ZNF367, E2F1",
    "FOXM1, CENPA, MYBL2, TFDP1, DNMT1, HMGA1, E2F7, PA2G4, ZNF367, E2F1",
    "FOXM1, CENPA, MYBL2, TFDP1, DNMT1, HMGA1, E2F7, PA2G4, ZNF367, E2F1",
    "FOXM1, CENPA, MYBL2, TFDP1, DNMT1, HMGA1, E2F7, PA2G4, ZNF367, E2F1",
    "FOXM1, CENPA, MYBL2, TFDP1, DNMT1, HMGA1, E2F7, PA2G4, ZNF367, E2F1",
    "FOXM1, CENPA, MYBL2, TFDP1, DNMT1, HMGA1, E2F7, PA2G4, ZNF367, E2F1",
    "FOXM1, CENPA, MYBL2, TFDP1, DNMT1, HMGA1, E2F7, PA2G4, ZNF367, E2F1",
    "FOXM1, CENPA, MYBL2, TFDP1, DNMT1, E2F7, ZNF367, E2F1",
    "FOXM1, CENPA, MYBL2, TFDP1, DNMT1, E2F7, ZNF367, E2F1",
    "FOXM1, CENPA, MYBL2, TFDP1, DNMT1, E2F7, ZNF367, E2F1",
    "FOXM1, CENPA, MYBL2, TFDP1, DNMT1, E2F7, ZNF367, E2F1",
    "FOXM1, CENPA, MYBL2, TFDP1, DNMT1, HMGA1, E2F7, PA2G4, E2F1",
    "FOXM1, CENPA, MYBL2, TFDP1, DNMT1, E2F7, ZNF367, E2F1",
    "FOXM1, CENPA, MYBL2, TFDP1, DNMT1, E2F7, ZNF367, E2F1",
    "FOXM1, CENPA, MYBL2, TFDP1, DNMT1, E2F7, ZNF367, E2F1",
    "FOXM1, MYBL2, TFDP1, HMGA1, E2F7, E2F1",
    "FOXM1, MYBL2, TFDP1, HMGA1, E2F7, E2F1",
    "FOXM1, MYBL2, TFDP1, HMGA1, E2F1",
    "E2F1"
  ),
  stringsAsFactors = FALSE
)

# Convert list of TFs into long format
df_long <- new_input_data %>%
  separate_rows(Set_List_String, sep = ",\\s*") %>%
  rename(Set = Set_List_String) %>%
  select(Set, Gene) %>%
  distinct(Set, Gene)

# ----------------------------------------------------------------------
# 2. Matrix Creation
# ----------------------------------------------------------------------

# Create binary TF–Gene matrix
intersection_matrix <- df_long %>%
  mutate(Presence = 1) %>%
  pivot_wider(
    names_from = Gene,
    values_from = Presence,
    values_fill = 0
  ) %>%
  column_to_rownames("Set")

# Add row sums
intersection_matrix$Gene_Count <- rowSums(intersection_matrix)

# Print matrix
print("--- Binary Intersection Matrix ---")
print(intersection_matrix)

# ----------------------------------------------------------------------
# 3. Visualization
# ----------------------------------------------------------------------

# Bar plot data
plot_data_bars <- intersection_matrix %>%
  rownames_to_column("Set") %>%
  select(Set, Gene_Count) %>%
  mutate(Set = factor(Set, levels = rownames(intersection_matrix)))

# Heatmap data
plot_data_heatmap <- intersection_matrix %>%
  rownames_to_column("Set") %>%
  select(-Gene_Count) %>%
  pivot_longer(
    cols = -Set,
    names_to = "Gene",
    values_to = "Presence"
  ) %>%
  mutate(Set = factor(Set, levels = rownames(intersection_matrix)))

# Heatmap
intersection_heatmap <- ggplot(plot_data_heatmap, aes(x = Gene, y = Set, fill = factor(Presence))) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("0" = "gray95", "1" = "#2a7886")) +
  scale_y_discrete(expand = expansion(mult = c(0, 0))) +
  coord_fixed(ratio = 1) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 9, face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    panel.grid.major = element_blank()
  ) +
  labs(
    x = "Key Genes", 
    y = "Transcription Factor (TF) / Regulator",
    title = "Transcription Factor — Gene Regulatory Matrix"
  )

# Align bar plot order with heatmap
plot_data_bars$Set <- factor(plot_data_bars$Set, levels = levels(intersection_heatmap$data$Set))

# Bar plot
side_bar_plot <- ggplot(plot_data_bars, aes(x = Gene_Count, y = Set, fill = Gene_Count)) +
  geom_col(width = 0.6, show.legend = FALSE, fill = "#fcae1e") +
  geom_text(aes(label = Gene_Count), hjust = -0.5, size = 3) +
  scale_y_discrete(expand = expansion(mult = c(0, 0))) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 9, face = "bold"),
    panel.grid.major = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  ) +
  labs(x = "Gene Count", y = "") +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, max(plot_data_bars$Gene_Count) * 1.40)
  )

# Combined plot
combined_plot <- intersection_heatmap + side_bar_plot +
  plot_layout(widths = c(7.5, 3), heights = c(12, 7))

print("--- Combined Visualization ---")
print(combined_plot)
