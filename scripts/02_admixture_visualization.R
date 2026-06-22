# admixture-visualisation.R =======================================
# Purpose: visualise ADMIXTURE K=5 ancestry proportions for 1000
#          Genomes chr22, faceted by superpopulation
# ==================================================================

# Setup ----
library(tidyverse)

BASE_DIR <- "/Users/tim/bioinformatics-portfolio/01-african-genomics"
K <- 5  # change this if  CV-error analysis points to a different K

# Load ADMIXTURE Q matrix ----
q_data <- read.table(
  file.path(BASE_DIR, sprintf("results/admixture/chr22_pruned.%d.Q", K))
)
colnames(q_data) <- paste0("Ancestry_", 1:K)

# Attach sample IDs (positional match - ADMIXTURE preserves .fam row order,
# there is no ID column in the .Q file itself) ----
fam <- read.table(file.path(BASE_DIR, "results/chr22_pruned.fam"))
q_data$sample_id <- fam$V2

# Load sample metadata ----
samples <- read_tsv(
  file.path(BASE_DIR, "data/raw/integrated_call_samples_v3.20130502.ALL.panel"),
  col_types = cols_only(
    sample    = col_character(),
    pop       = col_character(),
    super_pop = col_character(),
    gender    = col_character()
  )
)

# Merge and drop unmatched/related samples ----
q_data <- q_data %>%
  left_join(samples, by = c("sample_id" = "sample"))

cat("Rows with no superpopulation match:", sum(is.na(q_data$super_pop)), "\n")

q_data <- q_data %>%
  filter(!is.na(super_pop))

# Reshape to long format for ggplot ----
q_long <- q_data %>%
  pivot_longer(
    cols = starts_with("Ancestry_"),
    names_to = "component",
    values_to = "proportion"
  ) %>%
  arrange(super_pop, pop, sample_id)

# Lock bar order to the superpopulation grouping above, not alphabetical
# sample ID order (ggplot defaults to alphabetical factor levels otherwise)
q_long$sample_id <- factor(q_long$sample_id, levels = unique(q_long$sample_id))

# Build the stacked bar plot ----
p_adm <- ggplot(q_long, aes(x = sample_id, y = proportion, fill = component)) +
  geom_col(width = 1) +
  facet_grid(~super_pop, scales = "free_x", space = "free_x") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = sprintf("Ancestry Proportions (ADMIXTURE, K=%d)", K),
    x = NULL, y = "Ancestry proportion",
    fill = "Component"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold")
  )

# Save and inspect ----
dir.create(file.path(BASE_DIR, "results/figures"), recursive = TRUE, showWarnings = FALSE)

ggsave(
  file.path(BASE_DIR, sprintf("results/figures/admixture_K%d.png", K)),
  p_adm, width = 14, height = 5, dpi = 300
)

print(p_adm)


cv_errors <- data.frame(
  K     = 2:6,
  error = c(0.49710, 0.47744, 0.47178, 0.46659, 0.46576)
)

p_cv <- ggplot(cv_errors, aes(x = K, y = error)) +
  geom_line(color = "#1A3A5C", linewidth = 1) +
  geom_point(size = 4, color = "#1A3A5C") +
  geom_point(data = filter(cv_errors, K == 5),
             size = 5, color = "#D4600A") +
  scale_x_continuous(breaks = 2:6) +
  labs(
    title    = "ADMIXTURE cross-validation error by K",
    subtitle = "Orange point = selected K=5 (elbow: K5→K6 drop is 6x smaller than K4→K5)",
    x        = "K (number of ancestral populations)",
    y        = "Cross-validation error",
    caption  = "CV errors from ADMIXTURE --cv flag; lower = better fit"
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold"))

ggsave(
  file.path(BASE_DIR, "results/figures/admixture_cv_error.png"),
  p_cv, width = 7, height = 5, dpi = 300
)
print(p_cv)