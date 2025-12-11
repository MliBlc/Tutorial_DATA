#!/usr/bin/env Rscript

## ================================================================
## Kraken2 mock data + alpha/beta diversity demo
## - Generates mock species counts (50 species)
## - Generates mock sample metadata
## - Writes TSV files (so students see “real” inputs)
## - Reads them back in and runs diversity analyses
## ================================================================


library(tidyverse)
library(vegan)
library(ggplot2)


set.seed(42)

outdir <- "diversity_demo_output"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------
## 1. Generate MOCK Kraken2 species count matrix
## ------------------------------------------------
# Parameters
n_species <- 50
n_samples <- 12

species <- paste0("Species_", sprintf("%02d", 1:n_species))
samples <- paste0("Sample", sprintf("%02d", 1:n_samples))

# Groups for teaching: e.g. 4 controls, 4 TreatmentA, 4 TreatmentB
group_vec <- rep(c("Control", "TreatmentA", "TreatmentB"), each = 4)
names(group_vec) <- samples

# Simulate different underlying abundances by group
# (so beta-diversity / PCoA actually shows separation)
base_abund <- matrix(NA,
                     nrow = n_species,
                     ncol = length(unique(group_vec)),
                     dimnames = list(species, unique(group_vec)))

# Control: fairly even community
base_abund[, "Control"]     <- rgamma(n_species, shape = 2,  rate = 0.5)

# TreatmentA: some species enriched
base_abund[, "TreatmentA"]  <- rgamma(n_species, shape = 2,  rate = 0.5)
enriched_A <- sample(species, 10)
base_abund[enriched_A, "TreatmentA"] <- base_abund[enriched_A, "TreatmentA"] * 3

# TreatmentB: different subset enriched
base_abund[, "TreatmentB"]  <- rgamma(n_species, shape = 2,  rate = 0.5)
enriched_B <- sample(setdiff(species, enriched_A), 10)
base_abund[enriched_B, "TreatmentB"] <- base_abund[enriched_B, "TreatmentB"] * 3

# Simulate read counts per sample (e.g. 50k–80k reads)
total_reads <- sample(50000:80000, n_samples, replace = TRUE)
names(total_reads) <- samples

# For each sample, draw counts from a multinomial based on its group profile
counts_mat <- matrix(0,
                     nrow = n_species,
                     ncol = n_samples,
                     dimnames = list(species, samples))

for (s in samples) {
  g <- group_vec[s]
  probs <- base_abund[, g]
  probs <- probs / sum(probs)
  counts_mat[, s] <- as.vector(rmultinom(1, size = total_reads[s], prob = probs))
}

# Cast as a data frame like a Kraken2 species table
counts_df <- data.frame(
  Species = rownames(counts_mat),
  counts_mat,
  check.names = FALSE
)

counts_file <- file.path(outdir, "kraken_species_counts_mock.tsv")
write.table(counts_df, counts_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Mock species counts written to:\n  ", counts_file, "\n\n")

## ------------------------------------------------
## 2. Generate MOCK metadata
## ------------------------------------------------
meta_df <- data.frame(
  Sample = samples,
  Group  = group_vec,
  stringsAsFactors = FALSE
)

meta_file <- file.path(outdir, "sample_metadata_mock.tsv")
write.table(meta_df, meta_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Mock metadata written to:\n  ", meta_file, "\n\n")

## ------------------------------------------------
## 3. Read data back in (like real analysis)
## ------------------------------------------------
cat("Reading counts and metadata back in (as if real data)...\n")

counts_df_in <- read_tsv(counts_file, col_types = cols())
meta_in      <- read_tsv(meta_file,   col_types = cols())

counts_mat_in <- counts_df_in %>%
  column_to_rownames("Species") %>%
  as.matrix()

storage.mode(counts_mat_in) <- "numeric"

# community matrix: rows = samples, cols = species
comm <- t(counts_mat_in)

# Ensure sample order matches metadata
meta_in <- meta_in %>% filter(Sample %in% rownames(comm))
comm    <- comm[meta_in$Sample, , drop = FALSE]

cat("Community matrix dimensions (samples x species):",
    nrow(comm), "x", ncol(comm), "\n\n")

## ------------------------------------------------
## 4. Alpha diversity
## ------------------------------------------------
cat("Calculating alpha diversity...\n")

richness    <- specnumber(comm)                        # Observed species
shannon     <- diversity(comm, index = "shannon")      # Shannon index
simpson     <- diversity(comm, index = "simpson")      # Simpson index
inv_simpson <- diversity(comm, index = "invsimpson")   # Inverse Simpson
evenness    <- shannon / log(richness)                 # Pielou's evenness

alpha_df <- data.frame(
  Sample      = rownames(comm),
  Richness    = richness,
  Shannon     = shannon,
  Simpson     = simpson,
  InvSimpson  = inv_simpson,
  Evenness    = evenness,
  stringsAsFactors = FALSE
) %>%
  left_join(meta_in, by = "Sample")

alpha_file <- file.path(outdir, "alpha_diversity_mock.tsv")
write_tsv(alpha_df, alpha_file)

cat("Alpha diversity table written to:\n  ", alpha_file, "\n\n")

## Simple alpha diversity boxplot by group (Shannon)
alpha_plot_file <- file.path(outdir, "alpha_shannon_by_group.pdf")
p_alpha <- ggplot(alpha_df, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.1, size = 2) +
  theme_bw() +
  labs(title = "Alpha diversity (Shannon index)",
       y = "Shannon index", x = "Group")

ggsave(alpha_plot_file, p_alpha, width = 6, height = 4)
cat("Alpha diversity plot saved to:\n  ", alpha_plot_file, "\n\n")

## ------------------------------------------------
## 5. Beta diversity (Bray–Curtis & Jaccard)
## ------------------------------------------------
cat("Calculating beta diversity...\n")

bray    <- vegdist(comm, method = "bray")
jaccard <- vegdist(comm, method = "jaccard", binary = TRUE)

saveRDS(bray,    file.path(outdir, "beta_bray_mock.rds"))
saveRDS(jaccard, file.path(outdir, "beta_jaccard_mock.rds"))

# Long format distance tables (nice to show in Excel/R)
bray_df <- as.data.frame(as.matrix(bray)) %>%
  rownames_to_column("Sample1") %>%
  pivot_longer(-Sample1, names_to = "Sample2", values_to = "BrayCurtis")

jaccard_df <- as.data.frame(as.matrix(jaccard)) %>%
  rownames_to_column("Sample1") %>%
  pivot_longer(-Sample1, names_to = "Sample2", values_to = "Jaccard")

write_tsv(bray_df,    file.path(outdir, "beta_bray_long_mock.tsv"))
write_tsv(jaccard_df, file.path(outdir, "beta_jaccard_long_mock.tsv"))

cat("Beta diversity distance tables written.\n\n")

## ------------------------------------------------
## 6. PCoA on Bray–Curtis
## ------------------------------------------------
cat("Running PCoA on Bray–Curtis...\n")

bray_pcoa <- cmdscale(bray, eig = TRUE, k = 2)

pcoa_df <- data.frame(
  Sample = rownames(bray_pcoa$points),
  Axis1  = bray_pcoa$points[, 1],
  Axis2  = bray_pcoa$points[, 2],
  stringsAsFactors = FALSE
) %>%
  left_join(meta_in, by = "Sample")

# Variance explained
var_expl <- round(
  100 * bray_pcoa$eig[1:2] / sum(bray_pcoa$eig[bray_pcoa$eig > 0]),
  1
)

pcoa_file <- file.path(outdir, "bray_pcoa_scores_mock.tsv")
write_tsv(pcoa_df, pcoa_file)
cat("PCoA scores written to:\n  ", pcoa_file, "\n\n")

# PCoA plot
p_pcoa <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2, color = Group, label = Sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.7, size = 3) +
  theme_bw() +
  labs(
    title = "PCoA on Bray–Curtis (mock data)",
    x = paste0("Axis 1 (", var_expl[1], "%)"),
    y = paste0("Axis 2 (", var_expl[2], "%)")
  )

pcoa_plot_file <- file.path(outdir, "bray_pcoa_mock.pdf")
ggsave(pcoa_plot_file, p_pcoa, width = 6, height = 4)
cat("PCoA plot saved to:\n  ", pcoa_plot_file, "\n\n")

## ------------------------------------------------
## 7. Optional: PERMANOVA (adonis2) by group
## ------------------------------------------------
cat("Running PERMANOVA (adonis2) by Group...\n")
permanova_res <- adonis2(bray ~ Group, data = meta_in)
permanova_file <- file.path(outdir, "permanova_bray_by_group.txt")
capture.output(permanova_res, file = permanova_file)
cat("PERMANOVA result written to:\n  ", permanova_file, "\n\n")

## ------------------------------------------------
## 8. Optional: relative abundance barplot (top 10 species)
## ------------------------------------------------
cat("Creating relative abundance barplot (top 10 species)...\n")

rel_abund <- sweep(comm, 1, rowSums(comm), FUN = "/")
rel_abund_df <- as.data.frame(rel_abund) %>%
  rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "Species", values_to = "RelAbund") %>%
  left_join(meta_in, by = "Sample")

# Pick top 10 species overall
top_species <- rel_abund_df %>%
  group_by(Species) %>%
  summarize(MeanAbund = mean(RelAbund)) %>%
  arrange(desc(MeanAbund)) %>%
  slice_head(n = 10) %>%
  pull(Species)

rel_abund_top <- rel_abund_df %>%
  filter(Species %in% top_species)

barplot_file <- file.path(outdir, "top10_species_relative_abundance_mock.pdf")
p_bar <- ggplot(rel_abund_top,
                aes(x = Sample, y = RelAbund, fill = Species)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Top 10 species: relative abundance (mock data)",
    y = "Relative abundance",
    x = "Sample"
  )
ggsave(barplot_file, p_bar, width = 8, height = 4)
cat("Relative abundance barplot saved to:\n  ", barplot_file, "\n\n")

cat("=== DEMO COMPLETE ===\n")
