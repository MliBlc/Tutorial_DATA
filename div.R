#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
})

## ------------------------------------------------------------------
## User settings
## ------------------------------------------------------------------
counts_file <- "kraken_species_counts.tsv"   # Species x Samples
meta_file   <- "sample_metadata.tsv"         # Optional Sample, Group
outdir      <- "diversity_results"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------------
## Read data
## ------------------------------------------------------------------
cat("Reading counts from:", counts_file, "\n")
counts_df <- read_tsv(counts_file, col_types = cols())

# First column = Species, remaining = samples
counts_mat <- counts_df %>%
  column_to_rownames("Species") %>%
  as.matrix()

# Ensure numeric
storage.mode(counts_mat) <- "numeric"

# vegan expects community matrix with rows = samples, cols = species
comm <- t(counts_mat)

if (any(comm < 0, na.rm = TRUE)) {
  stop("Negative values detected in count matrix. Check your input.")
}

cat("Matrix dimensions (samples x species):", dim(comm)[1], "x", dim(comm)[2], "\n")

## ------------------------------------------------------------------
## Alpha diversity
## ------------------------------------------------------------------

cat("Computing alpha diversity...\n")

richness    <- specnumber(comm)                        # Observed species
shannon     <- diversity(comm, index = "shannon")      # Shannon index
simpson     <- diversity(comm, index = "simpson")      # Simpson index
inv_simpson <- diversity(comm, index = "invsimpson")   # Inverse Simpson

evenness <- shannon / log(richness)                    # Pielou's evenness

alpha_df <- data.frame(
  Sample      = rownames(comm),
  Richness    = richness,
  Shannon     = shannon,
  Simpson     = simpson,
  InvSimpson  = inv_simpson,
  Evenness    = evenness,
  row.names   = NULL
)

# Add metadata if available
if (file.exists(meta_file)) {
  cat("Joining sample metadata from:", meta_file, "\n")
  meta <- read_tsv(meta_file, col_types = cols())
  alpha_df <- alpha_df %>%
    left_join(meta, by = c("Sample" = "Sample"))
}

alpha_out <- file.path(outdir, "alpha_diversity.tsv")
write_tsv(alpha_df, alpha_out)
cat("Alpha diversity written to:", alpha_out, "\n")

## ------------------------------------------------------------------
## Beta diversity (distance matrices)
## ------------------------------------------------------------------

cat("Computing beta diversity (Bray-Curtis, Jaccard)...\n")

# Bray-Curtis on abundance
bray <- vegdist(comm, method = "bray")

# Jaccard on presence/absence
jaccard <- vegdist(comm, method = "jaccard", binary = TRUE)

saveRDS(bray,    file.path(outdir, "beta_bray.rds"))
saveRDS(jaccard, file.path(outdir, "beta_jaccard.rds"))

# Also export as long table (sample1, sample2, distance)
bray_df <- as.data.frame(as.matrix(bray)) %>%
  rownames_to_column("Sample1") %>%
  pivot_longer(-Sample1, names_to = "Sample2", values_to = "BrayCurtis")

jaccard_df <- as.data.frame(as.matrix(jaccard)) %>%
  rownames_to_column("Sample1") %>%
  pivot_longer(-Sample1, names_to = "Sample2", values_to = "Jaccard")

write_tsv(bray_df,    file.path(outdir, "beta_bray_long.tsv"))
write_tsv(jaccard_df, file.path(outdir, "beta_jaccard_long.tsv"))

cat("Beta diversity tables written.\n")

## ------------------------------------------------------------------
## Ordination: PCoA on Bray-Curtis
## ------------------------------------------------------------------

cat("Running PCoA on Bray-Curtis distances...\n")

bray_pcoa <- cmdscale(bray, eig = TRUE, k = 2)

pcoa_df <- data.frame(
  Sample = rownames(bray_pcoa$points),
  Axis1  = bray_pcoa$points[, 1],
  Axis2  = bray_pcoa$points[, 2],
  row.names = NULL
)

# Percent variance explained
var_explained <- round(100 * bray_pcoa$eig[1:2] / sum(bray_pcoa$eig[bray_pcoa$eig > 0]), 1)

if (file.exists(meta_file)) {
  meta <- read_tsv(meta_file, col_types = cols())
  pcoa_df <- pcoa_df %>%
    left_join(meta, by = c("Sample" = "Sample"))
}

pcoa_out <- file.path(outdir, "bray_pcoa_scores.tsv")
write_tsv(pcoa_df, pcoa_out)
cat("PCoA scores written to:", pcoa_out, "\n")

## ------------------------------------------------------------------
## Simple PCoA plot
## ------------------------------------------------------------------

plot_file <- file.path(outdir, "bray_pcoa_plot.pdf")

p <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2, label = Sample)) +
  geom_point(aes(color = if ("Group" %in% colnames(pcoa_df)) Group else NULL),
             size = 3) +
  geom_text(vjust = -0.8, size = 3) +
  theme_bw() +
  labs(
    title = "PCoA (Bray-Curtis)",
    x = paste0("Axis 1 (", var_explained[1], "%)"),
    y = paste0("Axis 2 (", var_explained[2], "%)")
  )

ggsave(plot_file, p, width = 6, height = 4)
cat("PCoA plot saved to:", plot_file, "\n")

cat("Done.\n")
