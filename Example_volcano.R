ggplot(as.data.frame(res1), aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = abs(log2FoldChange) > log2fc_threshold & padj < padj_threshold), alpha = 0.5) +
  scale_color_manual(values = c("black", "red"), guide = FALSE) +  # Red points for significant changes
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "blue") +  # Horizontal threshold line
  geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "blue") +  # Vertical threshold lines
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  # Remove gridlines
  xlim(c(-6, 6)) + 
  ylim(c(0, 10)) +
  ggtitle("Volcano Plot: Control vs Treatment") +
  xlab("Log2 Fold Change") + 
  ylab("-Log10 Adjusted p-value")

ggsave("Control_vs_Treatment_volcano.png", width = 8, height = 6)
