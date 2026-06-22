# fst-validation.R =================================================
# Purpose: compute per-variant Fst using PLINK2's Weir & Cockerham
#          estimator, and compare against project 1b's hand-rolled
#          Fst calculation (same biological question, rigorous tool)
# ===================================================================

# Setup ----
library(tidyverse)

BASE_DIR <- "/Users/tim/bioinformatics-portfolio/01-african-genomics"

# Load population labels, same filtering logic as the PCA script ----
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

# Reshape into the format PLINK2's --pheno flag expects: a column
# named "#IID" (the sample ID PLINK2 will match against) plus one
# column per category we want it to use
fst_covar <- samples %>%
  transmute(`#IID` = sample, super_pop = super_pop)

write_tsv(fst_covar, file.path(BASE_DIR, "results/fst_covariate.txt"))

cat("Wrote covariate file for", nrow(fst_covar), "samples\n")


# Read all 10 pairwise Fst files and average across pairs, as an
# approximation of "overall differentiation" comparable to 1b's
# single global multi-population Fst ----
fst_files <- list.files(
  file.path(BASE_DIR, "results"),
  pattern = "chr22_fst\\..*\\.fst\\.var$",
  full.names = TRUE
)

fst_files  # sanity check - should list all 10

fst_list <- map(fst_files, function(f) {
  read_tsv(f, show_col_types = FALSE) %>%
    mutate(comparison = basename(f))
})

fst_all <- bind_rows(fst_list)

# Average the 10 pairwise Fst values together, per variant
fst_avg <- fst_all %>%
  group_by(ID) %>%
  summarise(mean_pairwise_fst = mean(WC_FST, na.rm = TRUE))  # swap WC_FST if your file uses a different name

cat("=== PLINK2: mean of 10 pairwise WC Fst comparisons ===\n")
cat("Mean:  ", round(mean(fst_avg$mean_pairwise_fst, na.rm = TRUE), 4), "\n")
cat("Median:", round(median(fst_avg$mean_pairwise_fst, na.rm = TRUE), 4), "\n")
cat("Max:   ", round(max(fst_avg$mean_pairwise_fst, na.rm = TRUE), 4), "\n")

# Also report AFR vs EUR alone - the single most contrastive pair,
# and the most relevant one given this project's focus
fst_afr_eur <- read_tsv(
  file.path(BASE_DIR, "results/chr22_fst.AFR.EUR.fst.var"),
  show_col_types = FALSE
)

cat("\n=== PLINK2: AFR vs EUR pairwise Fst only ===\n")
cat("Mean:  ", round(mean(fst_afr_eur$WC_FST, na.rm = TRUE), 4), "\n")
cat("Median:", round(median(fst_afr_eur$WC_FST, na.rm = TRUE), 4), "\n")
cat("Max:   ", round(max(fst_afr_eur$WC_FST, na.rm = TRUE), 4), "\n")

# Visualise Fst results ----

# Plot 1: Distribution of AFR vs EUR Fst, one bar per SNP - mirrors
# 1b's Fst histogram, but using the rigorous Weir & Cockerham values ----
p_fst_hist <- ggplot(fst_afr_eur, aes(x = WC_FST)) +
  geom_histogram(bins = 60, fill = "#D4600A", alpha = 0.8) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "#1A3A5C") +
  geom_vline(xintercept = 0.3, linetype = "dashed", color = "#1A3A5C") +
  labs(
    title = "Distribution of AFR vs EUR Fst across chromosome 22 variants",
    subtitle = "Most variants are near-identical between groups (low Fst); a small tail are strong AIMs",
    x = "Fst (Weir & Cockerham, AFR vs EUR)",
    y = "Number of variants",
    caption = "Data: 1000 Genomes Phase 3, chr22, LD-pruned"
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(BASE_DIR, "results/figures/fst_afr_eur_distribution.png"),
       p_fst_hist, width = 9, height = 6, dpi = 300)
print(p_fst_hist)

# Plot 2: Mean Fst across all 10 superpopulation pairs - the "general"
# picture, showing which pairs are most/least differentiated ----
fst_pair_summary <- fst_all %>%
  group_by(comparison) %>%
  summarise(mean_fst = mean(WC_FST, na.rm = TRUE)) %>%
  mutate(
    pair_label = comparison %>%
      str_remove("chr22_fst\\.") %>%
      str_remove("\\.fst\\.var") %>%
      str_replace("\\.", " vs ")
  )

p_fst_pairs <- ggplot(fst_pair_summary, aes(x = reorder(pair_label, mean_fst), y = mean_fst)) +
  geom_col(fill = "#1A3A5C") +
  coord_flip() +
  labs(
    title = "Mean Fst by superpopulation pair",
    subtitle = "Higher = more genetically differentiated on chromosome 22",
    x = NULL, y = "Mean Fst (Weir & Cockerham)",
    caption = "Data: 1000 Genomes Phase 3, chr22, LD-pruned"
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(BASE_DIR, "results/figures/fst_pairwise_comparison.png"),
       p_fst_pairs, width = 9, height = 6, dpi = 300)
print(p_fst_pairs)

# print the top 10 AIMs
fst_afr_eur %>% arrange(desc(WC_FST)) %>% head(10)
