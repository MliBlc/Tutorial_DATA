############################
### 1. LOAD REQUIRED PACKAGES
############################

# ape      : phylogenetic trees, genetic distances, FST
# pegas    : haplotypes, haplotype networks, Hamming distance
# msa      : multiple sequence alignment
# ggtree   : phylogenetic tree visualization
# ggplot2  : plotting backend
# haplotypes: haplotype frequency calculations

install.packages("BiocManager")
BiocManager::install(c("msa", "ggtree"))
cran_packages <- c("tidyverse", "ape", "pegas", "ggplot2", "haplotypes", "adegenet", "hierfstat") 

install.packages(cran_packages) 


library(tidyverse)
library(ape)
library(pegas)
library(msa)
library(ggtree)
library(ggplot2)
library(haplotypes)
library(adegenet)
library(hierfstat)

# Input FASTA file (unaligned or aligned)
fasta_file <- "S2_Appendix.fas"


############################
### 3. READ FASTA AND ALIGN SEQUENCES
############################

# Read DNA sequences from FASTA file
dna <- readDNAStringSet(fasta_file)

# Perform multiple sequence alignment (MSA)
# This step may take time for large datasets
aln <- msa(dna)

# Convert alignment into DNAbin format (required by ape & pegas)
nbin <- as.DNAbin(aln)


############################
### 4. EXTRACT POPULATION NAMES FROM SAMPLE LABELS
############################

# Assumes sample names are formatted as:
# Population_IndividualNumber (e.g., Aksu_1)
# Everything before "_" is treated as population name
pop <- sub("_.*", "", rownames(nbin))

# Convert to factor for downstream analyses
pop_factor <- as.factor(pop)


############################
### 5. GENETIC DISTANCE AND NJ TREE (K80 MODEL)
############################

# Compute Kimura 2-Parameter (K80) genetic distances
dist_k80 <- dist.dna(nbin, model = "K80")

# Build Neighbor-Joining (NJ) phylogenetic tree
tree_nj <- nj(dist_k80)


############################
### 6. PLOT NJ TREE (LINEAR)
############################

pdf("NJ_tree.pdf", width = 8, height = 6)

ggtree(tree_nj) +
  geom_tiplab(size = 2) +        # display sample names
  geom_treescale()               # add scale bar

dev.off()


############################
### 7. PLOT NJ TREE (CIRCULAR, POPULATION-COLORED)
############################

pdf("NJ_tree_circular.pdf", width = 7, height = 7)

ggtree(tree_nj, layout = "circular") %<+%
  data.frame(label = rownames(nbin), population = pop) +
  aes(color = population) +      # color tips by population
  geom_tiplab(size = 2) +
  theme(legend.position = "right")

dev.off()

############################
### 8. HAPLOTYPE EXTRACTION
############################

haps <- pegas::haplotype(nbin)
# Number of individuals carrying each haplotype
hap_freq <- attr(haps, "freq")


############################
### 9. HAPLOTYPE–INDIVIDUAL AND POPULATION MATRICES
############################

# Get index linking individuals to haplotypes
hap_index <- attr(haps, "index")

# Haplotype composition per individual
ind_hap <- table(
  haplotype  = rep(seq_along(hap_index), lengths(hap_index)),
  individual = rownames(nbin)[unlist(hap_index)]
)

#Populasyon matrisi
pop_hap <- table(
  haplotype  = rep(seq_along(hap_index), lengths(hap_index)),
  population = pop[unlist(hap_index)]
)

write.table(
  as.data.frame(pop_hap),
  "haplotype_population_frequencies.txt",
  sep = "\t",
  quote = FALSE
)

############################
### 10. HAPLOTYPE FREQUENCY PER POPULATION
############################
# Calculate haplotype frequencies for each population
hf <- table(
  haplotype  = rep(seq_along(hap_index), lengths(hap_index)),
  population = pop[unlist(hap_index)]
)

hf_df <- as.data.frame.matrix(hf)

write.table(
  hf_df,
  "haplotype_population_frequencies.txt",
  sep = "\t",
  quote = FALSE
)

###########################
### 11. HAMMING DISTANCE AND HEATMAP
############################

# Convert haplotypes to matrix form
hap_mat <- as.matrix(haps)

# Calculate Hamming distance between haplotypes
hap_dist <- dist.hamming(hap_mat)

# Save Hamming distance matrix
write.table(
  as.matrix(hap_dist),
  "haplotype_hamming_matrix.txt",
  sep = "\t",
  quote = FALSE
)

# Plot heatmap of haplotype distances
pdf("haplotype_heatmap.pdf", width = 7, height = 7)

heatmap(
  as.matrix(hap_dist),
  scale = "none",
  col = heat.colors(100),
  symm = TRUE
)

dev.off()


############################
### 12. HAPLOTYPE NETWORK – INDIVIDUAL LEVEL
############################

# Build haplotype network
net <- haploNet(haps)

pdf("haplotype_network_individuals.pdf", width = 9, height = 5)

plot(
  net,
  size = attr(net, "freq"),   # node size proportional to frequency
  pie  = ind_hap,             # individual composition
  show.mutation = TRUE,
  scale.ratio = 2,
  cex = 0.7
)

dev.off()


############################
### 13. HAPLOTYPE NETWORK – POPULATION LEVEL
############################

pdf("haplotype_network_populations.pdf", width = 9, height = 5)

plot(
  net,
  size = attr(net, "freq"),
  pie  = pop_hap,             # population composition
  show.mutation = TRUE,
  scale.ratio = 2,
  cex = 0.7
)

dev.off()

############################
### 14. FST (WEIR & COCKERHAM)
############################

gen <- DNAbin2genind(nbin, pop = pop_factor)

#Global Weir & Cockerham FST
fst_wc <- wc(gen)
fst_global <- fst_wc$FST
fst_global

############################
### 15. PAIRWISE FST BETWEEN POPULATIONS
############################

#Weir & Cockerham (hierfstat) 
pairwise_fst <- pairwise.WCfst(gen)

#dosyaya yazma 
write.table(
  pairwise_fst,
  "pairwise_FST.txt",
  sep = "\t",
  quote = FALSE
)

############################
### 16. PAIRWISE FST HEATMAP
############################

pdf("pairwise_FST_heatmap.pdf", width = 7, height = 7)

heatmap(
  as.matrix(pairwise_fst),
  scale = "none",
  col = colorRampPalette(c("white", "orange", "red"))(100),
  symm = TRUE
)

dev.off()


############################
### 17. HAPLOTYPE NJ TREE WITH BOOTSTRAP SUPPORT
############################

#Haplotype dizileri (DNAbin)
hap_mat <- as.DNAbin(haps)

############################
### 17. HAPLOTYPE NJ TREE WITH BOOTSTRAP SUPPORT
############################

# Hamming distances among haplotypes
hap_dist <- dist.hamming(haps)

# Build NJ tree
hap_tree <- nj(hap_dist)

# Bootstrap support (100 replicates)
boot <- boot.phylo(
    hap_tree,
    haps,
    B = 100,
    function(x) nj(dist.hamming(x))
)

# Assign bootstrap to node labels (numeric)
hap_tree$node.label <- boot

# Convert tree for ggtree
tree_df <- as_tibble(hap_tree)

# PDF output
pdf("haplotype_tree_bootstrap.pdf", width = 8, height = 5)

ggtree(hap_tree) +
    geom_tiplab(size = 3) +
    geom_nodelab(aes(label = label), hjust = -0.3, size = 3) +
    theme_tree2()

dev.off()
