# scripts/07_within_africa_fst.R ======================================
# Purpose: pairwise Fst between the 7 African 1000 Genomes populations
#          specifically - completes the within-Africa trio
#          (PCA + ADMIXTURE + Fst, all restricted to Africa)
# =====================================================================

library(tidyverse)
BASE_DIR <- "/Users/tim/bioinformatics-portfolio/01-african-genomics"

samples <- read_tsv(
  file.path(BASE_DIR, "data/raw/integrated_call_samples_v3.20130502.ALL.panel"),
  col_types = cols_only(sample = col_character(), pop = col_character(),
                        super_pop = col_character(), gender = col_character())
) %>% filter(!is.na(super_pop))

# Covariate file using `pop` (not super_pop) as the grouping column,
# restricted to AFR samples - same shape as script 03's file, just a
# different column and a smaller sample set ----
fst_covar_pop <- samples %>%
  filter(super_pop == "AFR") %>%
  transmute(`#IID` = sample, pop = pop)

write_tsv(fst_covar_pop, file.path(BASE_DIR, "results/fst_covariate_afr_pop.txt"))


fst_files_afr <- list.files(
  file.path(BASE_DIR, "results"),
  pattern = "chr22_fst_within_afr\\..*\\.fst\\.var$",
  full.names = TRUE
)

length(fst_files_afr)  # should print 21

fst_afr_summary <- map(fst_files_afr, function(f) {
  read_tsv(f, show_col_types = FALSE) %>%
    mutate(comparison = basename(f))
}) %>%
  bind_rows() %>%
  group_by(comparison) %>%
  summarise(mean_fst = mean(WC_FST, na.rm = TRUE)) %>%   # confirm column name first, as before
  mutate(
    pair_label = comparison %>%
      str_remove("chr22_fst_within_afr\\.") %>%
      str_remove("\\.fst\\.var") %>%
      str_replace("\\.", " vs ")
  )

p_fst_within_afr <- ggplot(fst_afr_summary, aes(x = reorder(pair_label, mean_fst), y = mean_fst)) +
  geom_col(fill = "#1A3A5C") +
  coord_flip() +
  labs(
    title = "Mean Fst between African 1000 Genomes populations",
    subtitle = "Higher = more genetically differentiated on chromosome 22",
    x = NULL, y = "Mean Fst (Weir & Cockerham)",
    caption = "Data: 1000 Genomes Phase 3, chr22, LD-pruned"
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(BASE_DIR, "results/figures/fst_within_africa.png"),
       p_fst_within_afr, width = 9, height = 8, dpi = 300)
print(p_fst_within_afr)