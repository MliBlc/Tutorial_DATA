library(ggplot2)
library(dplyr)
library(stringr)

# Load the data
Data_degisecek <- read.csv("degisecek.csv", sep = ",")

# Calculate the -log10(Enrichment FDR)
Data_degisecek$neg_log_FDR <- -log10(Data_degisecek$Enrichment.FDR)

# Filter the data for FDR < 0.05
Data_degisecek_filtered <- Data_degisecek[Data_degisecek$Enrichment.FDR < 0.05, ]

# Select the top 10 pathways according to Fold Enrichment
Data_degisecek_top10 <- Data_degisecek_filtered %>%
  arrange(desc(Fold.Enrichment)) %>%
  head(10)

# Wrap pathway names longer than 20 characters
Data_degisecek_top10$Pathway <- str_wrap(Data_degisecek_top10$Pathway, width = 20)

# Create the plot
Data_degisecek_pl <- ggplot(Data_degisecek_top10, aes(x = Fold.Enrichment, y = reorder(Pathway, Fold.Enrichment))) +
  geom_segment(aes(x = 0, xend = Fold.Enrichment, yend = Pathway, color = neg_log_FDR), size = 1.3 * 0.6) +
  geom_point(aes(size = nGenes, color = neg_log_FDR)) +
  scale_color_gradient2(low = "blue", mid = "purple", high = "red", midpoint = median(Data_degisecek_top10$neg_log_FDR)) +
  guides(
    color = guide_colorbar(order = 1),
    size = guide_legend(order = 2)
  ) +
  labs(x = "Fold Enrichment", y = "",
       title = "Top 10 degisecek Enrichment Plot", # Removed underscores
       color = "-log10(Enrichment FDR)",
       size = "nGenes") +
  theme_minimal() +
  theme(
    text = element_text(family = "Times New Roman"), # Set font to Times New Roman
    plot.title = element_text(face = "bold"), # Bold title
    axis.text.y = element_text(size = 12, face = "bold", color = "black"), # Bold y-axis labels
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid = element_blank()
  )

# Save the plot as a 300 DPI TIFF file
ggsave("degisecek_top10.tiff", plot = Data_degisecek_pl, device = "tiff", width = 10, height = 8, dpi = 300)
