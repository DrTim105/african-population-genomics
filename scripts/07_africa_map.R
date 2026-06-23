# =============================================================================
# scripts/07_africa_map.R
# Purpose : Geographic map of all 7 AFR-superpopulation sampling locations
#           from 1000 Genomes Phase 3, coloured by mean Fst from YRI.
# =============================================================================
# TERMINAL PREREQUISITES — this script reads 21 within-Africa pairwise Fst
# files created by script 06's PLINK2 command. Ensure the following have
# both been completed before running this script:
#
# 1. Script 04 run completely (creates results/keep_afr.txt)
# 2. Script 06 run to the write_tsv step, then terminal command:
#
# plink2 \
#   --pfile results/chr22_pruned \
#   --keep results/keep_afr.txt \
#   --pheno results/fst_covariate_afr_pop.txt \
#   --fst pop method=wc report-variants \
#   --out results/chr22_fst_within_afr
#
# Output files this script reads:
#   results/chr22_fst_within_afr.*.fst.var   (21 pairwise files)
#
# Dependencies: install.packages("maps") if not already installed
# =============================================================================

library(tidyverse)
install.packages("maps")
library(maps)   # for map_data("world")

BASE_DIR <- "/Users/tim/bioinformatics-portfolio/01-african-genomics"


# Population sampling coordinates -----------------------------------------
# Source: igsr_populations.tsv downloaded from internationalgenome.org
# These are the OFFICIAL IGSR-published coordinates, verified against
# the igsr_populations.tsv file provided in this project's data/raw folder.
# All previous Google Maps approximations replaced with these values.
#
# Note on ESN: IGSR coordinate (9.07N, 7.48E) falls in the Abuja/Niger
# State area rather than the Esan homeland in Edo State. This is the
# official published centroid.
#
# Note on ASW: "African Ancestry in SW USA" has no single sampling city.
# IGSR publishes (35.483, -97.533) as its centroid - included here.

pop_locations <- tribble(
  ~pop,  ~location,                   ~lat,       ~lon,
  "YRI", "Ibadan, Nigeria",            7.400000,   3.920000,
  "ESN", "Esan in Nigeria",            9.066660,   7.483333,
  "GWD", "Western Division, Gambia",  13.454876, -16.579032,
  "MSL", "Mende in Sierra Leone",      8.480000, -13.230000,
  "LWK", "Webuye, Kenya",             -1.270000,  36.610000,
  "ACB", "Barbados",                  13.100000, -59.620000,
  "ASW", "Southwest USA",             35.483000, -97.533330
)


# Compute mean Fst from YRI for each population ---------------------------
# Reads the 21 within-Africa pairwise Fst files from script 06, extracts
# only the 6 pairs involving YRI, then identifies the non-YRI partner.
# YRI itself gets Fst = 0 (it is the reference population).

fst_files_afr <- list.files(
  file.path(BASE_DIR, "results"),
  pattern = "chr22_fst_within_afr\\..*\\.fst\\.var$",
  full.names = TRUE
)

cat("Fst files found:", length(fst_files_afr), "(expect 21)\n")

# Read all 21 files, stack into one table, average per comparison
fst_afr_summary <- purrr::map(fst_files_afr, function(f) {
  read_tsv(f, show_col_types = FALSE) %>%
    mutate(comparison = basename(f))
}) %>%
  bind_rows() %>%
  group_by(comparison) %>%
  summarise(mean_fst = mean(WC_FST, na.rm = TRUE)) %>%
  mutate(
    pair_label = comparison %>%
      str_remove("chr22_fst_within_afr\\.") %>%   # strip filename prefix
      str_remove("\\.fst\\.var") %>%               # strip filename suffix
      str_replace("\\.", " vs ")                   # "ACB.YRI" -> "ACB vs YRI"
  )

# From pairs containing YRI, extract the non-YRI partner population
fst_from_yri <- fst_afr_summary %>%
  filter(str_detect(pair_label, "YRI")) %>%
  mutate(
    pop_a = str_split_fixed(pair_label, " vs ", 2)[, 1],
    pop_b = str_split_fixed(pair_label, " vs ", 2)[, 2],
    # whichever side is NOT YRI is the population we want
    other_pop = if_else(pop_a == "YRI", pop_b, pop_a)
  ) %>%
  select(other_pop, mean_fst)

print(fst_from_yri)   # should show 6 rows: ACB, ASW, ESN, GWD, LWK, MSL

# Merge coordinates with Fst values; YRI gets Fst = 0 by definition
map_data_pops <- map_data_pops %>%
  mutate(
    label_lon = lon + case_when(
      pop == "YRI" ~ -4.0,   # push YRI label left
      pop == "ESN" ~  4.0,   # push ESN label right
      TRUE         ~  0.0
    ),
    label_lat = lat + case_when(
      pop == "YRI" ~ -3.0,   # push YRI label below the dot
      TRUE         ~  2.0    # all others go above the dot
    )
  )

cat("\nMap data (verify all 7 rows have non-NA Fst):\n")
print(map_data_pops)


# Build the map -----------------------------------------------------------
# Map extent is deliberately wide enough to show all 7 populations:
# Africa (lon 3-37), Caribbean/Barbados (lon -59), SW USA (lon -97).
# xlim = c(-110, 50) captures all; ylim = c(-12, 52) captures all latitudes.

world <- map_data("world")

p_map <- ggplot() +
  
  # Draw country polygons from the maps package.
  # group = group is essential - without it ggplot connects every point
  # into one tangled line instead of separate country shapes.
  # linewidth = 0.2 keeps country borders thin and unobtrusive.
  geom_polygon(
    data  = world,
    aes(x = long, y = lat, group = group),
    fill  = "grey88",
    color = "white",
    linewidth = 0.2
  ) +
  
  # Population location dots - size AND color both encode Fst from YRI,
  # making high-Fst populations visually prominent two ways at once.
  # This is intentional redundant encoding, not an error.
  geom_point(
    data  = map_data_pops,
    aes(x = lon, y = lat, color = mean_fst, size = mean_fst),
    alpha = 0.9
  ) +
  
  # Population code labels - vjust = -1.2 places them just above each dot.
  # fontface = "bold" makes them legible against the map background.
  geom_text(
    data  = map_data_pops,
    aes(x = label_lon, y = label_lat, label = pop),
    vjust    = -1.2,
    size     = 3.8,
    fontface = "bold",
    color    = "grey20"
  ) +
  
  # Color scale: low Fst = dark blue (genetically close to YRI);
  # high Fst = orange (more differentiated from YRI).
  # Uses the same two colors as the rest of this project.
  scale_color_gradient(
    low  = "#1A3A5C",
    high = "#D4600A",
    name = "Mean Fst\nfrom YRI"
  ) +
  
  # Size scale: range = c(3, 9) maps 0 Fst to dot size 3,
  # maximum Fst to dot size 9. guide = "none" hides the size legend
  # since color already tells the same story.
  scale_size_continuous(
    range = c(3, 9),
    guide = "none"
  ) +
  
  # Crop to the window containing all 7 populations.
  # coord_fixed() also locks the lat/lon aspect ratio so the
  # map doesn't look stretched - important for geographic accuracy.
  coord_fixed(
    xlim = c(-110, 50),
    ylim = c(-12,  52)
  ) +
  
  labs(
    title    = "1000 Genomes AFR superpopulation sampling locations",
    subtitle = paste0(
      "Coloured by mean Fst from YRI (reference = 0)\n",
      "· YRI: Yoruba (Nigeria)  ·  ESN: Esan (Nigeria)\n",
      "· GWD: Gambian Mandinka  ·  MSL: Mende (Sierra Leone)\n",
      "· LWK: Luhya (Kenya)  ·  ACB: African Caribbean (Barbados)\n",
      "· ASW: African-American (SW USA) \n"
    ),
    x = NULL,
    y = NULL,
    caption = paste0(
      "Coordinates: IGSR official population metadata (igsr_populations.tsv). ",
      "ASW coordinate is IGSR-published centroid.\n",
      "Fst: Weir & Cockerham, PLINK2, chr22 LD-pruned variants (script 06)."
    )
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.caption  = element_text(size = 8, color = "grey50"),
    panel.grid    = element_blank(),
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    legend.position = "right"
  )

ggsave(
  file.path(BASE_DIR, "results/figures/africa_population_map.png"),
  p_map,
  width  = 12,
  height = 7,
  dpi    = 300
)

print(p_map)