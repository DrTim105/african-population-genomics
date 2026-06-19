# pca-visualisation.R =============================================
# Purpose: visualise PLINK2 PCA results for 1000 Genomes chr22,
#          coloured by 1000 Genomes superpopulation label
# ==================================================================

# Setup ----
library(tidyverse)

BASE_DIR <- "/Users/tim/bioinformatics-portfolio/01-african-genomics"

# Load PCA results ----
eigenvec <- read.table(
  file.path(BASE_DIR, "results/chr22_pca.eigenvec"),
  header = FALSE,
  stringsAsFactors = FALSE
)

# This file has a single ID column (IID only, no FID) - confirmed by
# inspecting the raw header: "#IID  PC1  PC2 ..."
colnames(eigenvec) <- c("IID", paste0("PC", 1:(ncol(eigenvec) - 1)))

cat("IID sample:\n")
print(head(eigenvec$IID, 3))

eigenval <- scan(file.path(BASE_DIR, "results/chr22_pca.eigenval"))

# Load sample metadata ----

# Peek at raw structure before parsing
readLines(
  file.path(BASE_DIR, "data/raw/integrated_call_samples_v3.20130502.ALL.panel"),
  n = 4
)

# col_types = cols_only(...) parses ONLY these four named columns and
# silently ignores the two unnamed trailing columns the panel file has
# (stray trailing tabs on some lines) - avoids the readr "New names"
# warning entirely, instead of select()-ing them away after the fact
samples <- read_tsv(
  file.path(BASE_DIR, "data/raw/integrated_call_samples_v3.20130502.ALL.panel"),
  col_types = cols_only(
    sample    = col_character(),
    pop       = col_character(),
    super_pop = col_character(),
    gender    = col_character()
  )
)

head(samples)
table(samples$super_pop)

# Merge PCA results with population metadata ----
eigenvec$IID  <- as.character(eigenvec$IID)
samples$sample <- as.character(samples$sample)

pca_data <- eigenvec %>%
  left_join(samples, by = c("IID" = "sample"))

# Diagnose the join BEFORE dropping anything ----
cat("Total rows after join:", nrow(pca_data), "\n")
cat("Rows with no superpopulation match:", sum(is.na(pca_data$super_pop)), "\n")

# These unmatched rows are 1000 Genomes high-coverage trio-completion
# samples (3,202 total in this VCF release) absent from the original
# 2,504-sample Phase 3 panel. They're also close relatives (parents/
# children) of samples already in the cohort, which violates the
# "unrelated individuals" assumption PCA relies on - so dropping them
# is the statistically correct choice here, not just a labelling fix.
pca_data <- pca_data %>%
  filter(!is.na(super_pop))

cat("Rows after filtering unmatched/related samples:", nrow(pca_data), "\n")

# Variance explained ----
pve <- data.frame(
  PC = 1:length(eigenval),
  variance_explained = eigenval / sum(eigenval) * 100
)

# Build the PCA plot ----
p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = super_pop)) +
  geom_point(size = 2, alpha = 0.75) +
  scale_color_brewer(palette = "Set1", name = "Superpopulation") +
  labs(
    title = "Population Structure: 1000 Genomes Chromosome 22",
    subtitle = sprintf(
      "Principal Component Analysis of %s individuals",
      format(nrow(pca_data), big.mark = ",")
    ),
    x = sprintf("PC1 (%.1f%% variance)", pve$variance_explained[1]),
    y = sprintf("PC2 (%.1f%% variance)", pve$variance_explained[2]),
    caption = "Data: 1000 Genomes Project Phase 3"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

# Save and inspect ----
dir.create(file.path(BASE_DIR, "results/figures"), recursive = TRUE, showWarnings = FALSE)

ggsave(
  file.path(BASE_DIR, "results/figures/pca_plot.png"),
  p_pca, width = 10, height = 7, dpi = 300
)

print(p_pca)

dim(pca_data)
head(pca_data)
sum(is.na(pca_data$PC1))
sum(is.na(pca_data$super_pop))