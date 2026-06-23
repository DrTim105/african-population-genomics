# scripts/04_within_superpop_pca.R ===================================
# Purpose: PCA run separately within AFR and within EUR, using
#          PLINK2's --keep to subset samples - directly comparable to
#          1b's within-Africa / within-Europe sklearn PCA, but using
#          PLINK2's genetics-aware PCA instead of raw-dosage PCA
# =====================================================================

# Setup ----
library(tidyverse)

BASE_DIR <- "/Users/tim/bioinformatics-portfolio/01-african-genomics"

samples <- read_tsv(
  file.path(BASE_DIR, "data/raw/integrated_call_samples_v3.20130502.ALL.panel"),
  col_types = cols_only(
    sample    = col_character(),
    pop       = col_character(),
    super_pop = col_character(),
    gender    = col_character()
  )
) %>%
  filter(!is.na(super_pop))

# Write one --keep file per superpopulation of interest ----
write_keep_file <- function(pop_code, samples_df, out_path) {
  samples_df %>%
    filter(super_pop == pop_code) %>%
    transmute(`#IID` = sample) %>%
    write_tsv(out_path)
}

write_keep_file("AFR", samples, file.path(BASE_DIR, "results/keep_afr.txt"))
write_keep_file("EUR", samples, file.path(BASE_DIR, "results/keep_eur.txt"))

# Within-Africa PCA, continental populations only - removes ACB/ASW
# (admixed) so PC1/PC2 aren't dominated by the admixture axis ----
continental_afr_pops <- c("ESN", "GWD", "LWK", "MSL", "YRI")

samples %>%
  filter(pop %in% continental_afr_pops) %>%
  transmute(`#IID` = sample) %>%
  write_tsv(file.path(BASE_DIR, "results/keep_afr_continental.txt"))


cat("AFR samples:", sum(samples$super_pop == "AFR"), "\n")
cat("EUR samples:", sum(samples$super_pop == "EUR"), "\n")


# =============================================================================
# STOP — RUN THESE THREE COMMANDS IN TERMINAL BEFORE CONTINUING
# (run from project root: 01-african-genomics/)
#
# The R code above created three --keep files (keep_afr.txt, keep_eur.txt,
# keep_afr_continental.txt). PLINK2 now runs PCA restricted to each subset.
#
# Within-AFR PCA (all 7 populations including ACB/ASW)
# plink2 --pfile results/chr22_pruned \
#   --keep results/keep_afr.txt \
#   --pca 10 --out results/chr22_pca_afr
#
# Within-EUR PCA (5 European populations)
# plink2 --pfile results/chr22_pruned \
#   --keep results/keep_eur.txt \
#   --pca 10 --out results/chr22_pca_eur
#
# Within-AFR PCA, continental only (ACB and ASW excluded)
# plink2 --pfile results/chr22_pruned \
#   --keep results/keep_afr_continental.txt \
#   --pca 10 --out results/chr22_pca_afr_continental
#
# Expected outputs per run: .eigenvec and .eigenval files
# =============================================================================


# Within-AFR PCA ----
eigenvec_afr <- read.table(
  file.path(BASE_DIR, "results/chr22_pca_afr.eigenvec"),
  header = FALSE, stringsAsFactors = FALSE
)
colnames(eigenvec_afr) <- c("IID", paste0("PC", 1:(ncol(eigenvec_afr) - 1)))
eigenvec_afr$IID <- as.character(eigenvec_afr$IID)

eigenval_afr <- scan(file.path(BASE_DIR, "results/chr22_pca_afr.eigenval"))
pve_afr <- data.frame(PC = 1:length(eigenval_afr),
                      variance_explained = eigenval_afr / sum(eigenval_afr) * 100)

pca_afr <- eigenvec_afr %>%
  left_join(samples, by = c("IID" = "sample"))

p_pca_afr <- ggplot(pca_afr, aes(x = PC1, y = PC2, color = pop)) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_brewer(palette = "Set1", name = "Population") +
  labs(
    title = "Within-Africa PCA — 1000 Genomes AFR (chr22)",
    subtitle = sprintf("%s individuals across %d populations",
                       format(nrow(pca_afr), big.mark = ","),
                       n_distinct(pca_afr$pop)),
    x = sprintf("PC1 (%.1f%% variance)", pve_afr$variance_explained[1]),
    y = sprintf("PC2 (%.1f%% variance)", pve_afr$variance_explained[2]),
    caption = "Data: 1000 Genomes Project Phase 3, chr22"
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(BASE_DIR, "results/figures/pca_within_afr.png"),
       p_pca_afr, width = 9, height = 7, dpi = 300)
print(p_pca_afr)

# Continental only Within-AFR PCA
eigenvec_afr_c <- read.table(
  file.path(BASE_DIR, "results/chr22_pca_afr_continental.eigenvec"),
  header = FALSE, stringsAsFactors = FALSE
)
colnames(eigenvec_afr_c) <- c("IID", paste0("PC", 1:(ncol(eigenvec_afr_c) - 1)))
eigenvec_afr_c$IID <- as.character(eigenvec_afr_c$IID)

eigenval_afr_c <- scan(file.path(BASE_DIR, "results/chr22_pca_afr_continental.eigenval"))
pve_afr_c <- data.frame(PC = 1:length(eigenval_afr_c),
                        variance_explained = eigenval_afr_c / sum(eigenval_afr_c) * 100)

pca_afr_c <- eigenvec_afr_c %>% left_join(samples, by = c("IID" = "sample"))

p_pca_afr_c <- ggplot(pca_afr_c, aes(x = PC1, y = PC2, color = pop)) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_brewer(palette = "Set1", name = "Population") +
  labs(
    title = "Within-Africa PCA — continental populations only",
    subtitle = "ACB and ASW (admixed) excluded to isolate geographic structure",
    x = sprintf("PC1 (%.1f%% variance)", pve_afr_c$variance_explained[1]),
    y = sprintf("PC2 (%.1f%% variance)", pve_afr_c$variance_explained[2]),
    caption = "Data: 1000 Genomes Phase 3, chr22"
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(BASE_DIR, "results/figures/pca_within_afr_continental.png"),
       p_pca_afr_c, width = 9, height = 7, dpi = 300)

print(p_pca_afr_c)



# Within-EUR PCA — identical structure, different files ----
eigenvec_eur <- read.table(
  file.path(BASE_DIR, "results/chr22_pca_eur.eigenvec"),
  header = FALSE, stringsAsFactors = FALSE
)
colnames(eigenvec_eur) <- c("IID", paste0("PC", 1:(ncol(eigenvec_eur) - 1)))
eigenvec_eur$IID <- as.character(eigenvec_eur$IID)

eigenval_eur <- scan(file.path(BASE_DIR, "results/chr22_pca_eur.eigenval"))
pve_eur <- data.frame(PC = 1:length(eigenval_eur),
                      variance_explained = eigenval_eur / sum(eigenval_eur) * 100)

pca_eur <- eigenvec_eur %>%
  left_join(samples, by = c("IID" = "sample"))

p_pca_eur <- ggplot(pca_eur, aes(x = PC1, y = PC2, color = pop)) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_brewer(palette = "Set1", name = "Population") +
  labs(
    title = "Within-Europe PCA — 1000 Genomes EUR (chr22)",
    subtitle = sprintf("%s individuals across %d populations",
                       format(nrow(pca_eur), big.mark = ","),
                       n_distinct(pca_eur$pop)),
    x = sprintf("PC1 (%.1f%% variance)", pve_eur$variance_explained[1]),
    y = sprintf("PC2 (%.1f%% variance)", pve_eur$variance_explained[2]),
    caption = "Data: 1000 Genomes Project Phase 3, chr22"
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(BASE_DIR, "results/figures/pca_within_eur.png"),
       p_pca_eur, width = 9, height = 7, dpi = 300)
print(p_pca_eur)
