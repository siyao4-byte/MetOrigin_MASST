####################################################################################################
# Graph 2: domainMASST Sankey + UpSet plots for MetOrigin-unmatched matched features
#
# Purpose
# - Visualize all domainMASST matched features.
# - Show domain overlap with an UpSet-style plot.
# - Show feature -> domain -> source summary with a Graph-1-style Sankey plot.
# - For phylogenetic/taxonomic results, collapse taxonomy to phylum level only.
#
# Inputs expected
#   outputs/domainmasst_unmatched_ms2/testing_summary/domainmasst_testing_summary_feature_domain_presence.csv
#   outputs/domainmasst_unmatched_ms2/testing_summary/domainmasst_testing_summary_json_subclasses_long.csv
#
# If these files are missing, first run:
#   .\.venv310\Scripts\python.exe scripts\summarize_domainmasst_results.py
####################################################################################################


# ==================================================================================================
# 0. Libraries
# ==================================================================================================

required_packages <- c(
  "tidyverse",
  "ggplot2",
  "readr",
  "stringr",
  "forcats",
  "scales",
  "patchwork"
)

options(repos = c(CRAN = "https://cloud.r-project.org"))

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(tidyverse)
library(ggplot2)
library(readr)
library(stringr)
library(forcats)
library(scales)
library(patchwork)


# ==================================================================================================
# 1. User settings
# ==================================================================================================

summary_dir <- "outputs/domainmasst_unmatched_ms2/testing_summary"

feature_domain_file <- file.path(summary_dir, "domainmasst_testing_summary_feature_domain_presence.csv")
json_long_file <- file.path(summary_dir, "domainmasst_testing_summary_json_subclasses_long.csv")
feature_metadata_file <- "outputs/domainmasst_unmatched_ms2/MetOrigin_unmatched_features_with_MS2.csv"

out_dir <- "Graph2_domainMASST_R_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# `combined` is an aggregate domain produced by domainMASST, so it duplicates domain-specific evidence.
# Keep it out of overlap/Sankey plots unless you explicitly want the aggregate tree.
domains_to_show <- c("plant", "microbe", "microbiome", "tissue", "food")
domainmasst_intersection_colour <- "#7B61C9"

# Taxonomic/phylogenetic domains are collapsed to phylum.
taxonomy_domains <- c("plant", "microbe")
taxonomy_rank_to_show <- "phylum"

# For non-taxonomic domains, use domain-specific source labels from the JSON tree.
non_taxonomy_domains <- c("microbiome", "tissue", "food")

# Keep the Sankey readable. It uses one feature tick per matched feature, so multi-origin features
# visibly branch into multiple domain/source flows. Lower-frequency right-side labels can be collapsed.
top_n_source_labels_per_domain <- 6
collapse_other_source_labels <- TRUE
show_feature_labels_in_sankey <- TRUE

# Optional domain-specific overrides for the right-most Sankey labels.
# This keeps low-count domains from producing many tiny overlapping labels.
top_n_source_labels_override <- c(
  plant = 2,
  microbe = 4,
  microbiome = 6,
  tissue = 6,
  food = 2
)

# Feature labels can explode the left side. They are exported in the CSV even when hidden in the plot.
max_feature_label_chars <- 52

# If TRUE, a feature-domain edge is repeated once for each phylum/source observed in that domain.
# This makes within-domain multi-origin evidence visible, e.g. one feature -> Plant -> multiple phyla.
show_within_domain_multi_source_branching <- TRUE

# Plot sizes
sankey_width <- 15
sankey_height <- 12
upset_width <- 12
upset_height <- 8
plot_dpi <- 300


# ==================================================================================================
# 2. Helper functions
# ==================================================================================================

parse_bool <- function(x) {
  str_to_lower(str_trim(as.character(x))) %in% c("true", "yes", "y", "1")
}

short_label <- function(x, max_chars = 40) {
  x <- as.character(x)
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 3), "..."), x)
}

clean_source_label <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("_", " ") %>%
    str_replace_all("\\s+", " ") %>%
    str_trim()
}

is_generic_source_label <- function(x) {
  z <- str_to_lower(str_trim(as.character(x)))
  z %in% c(
    "", "na", "nan", "none", "unknown", "root",
    "origin", "host", "health phenotype",
    "plants", "microbes", "microbiome", "tissue",
    "no matched json subclass", "json parse failed"
  )
}

make_bezier_points <- function(x0, y0, x1, y1, n = 80) {
  cx0 <- x0 + 0.10
  cy0 <- y0
  cx1 <- x1 - 0.10
  cy1 <- y1

  t <- seq(0, 1, length.out = n)

  tibble(
    x = (1 - t)^3 * x0 + 3 * (1 - t)^2 * t * cx0 + 3 * (1 - t) * t^2 * cx1 + t^3 * x1,
    y = (1 - t)^3 * y0 + 3 * (1 - t)^2 * t * cy0 + 3 * (1 - t) * t^2 * cy1 + t^3 * y1
  )
}

scale_lw <- function(value, max_value, min_lw = 0.5, max_lw = 8.0) {
  if (is.na(max_value) || max_value <= 0) return(rep(min_lw, length(value)))
  min_lw + (max_lw - min_lw) * value / max_value
}

allocate_layer_positions <- function(node_df, count_col = "n_features", top = 0.97, bottom = 0.03, gap = 0.006) {
  node_df <- node_df %>% mutate(.count = .data[[count_col]])

  if (nrow(node_df) == 0) {
    return(node_df %>% mutate(y = numeric(0), h = numeric(0), top = numeric(0), bottom = numeric(0)))
  }

  if (nrow(node_df) == 1) {
    return(node_df %>% mutate(y = 0.50, h = 0.18, top = 0.59, bottom = 0.41) %>% select(-.count))
  }

  total_space <- top - bottom
  usable_space <- total_space - gap * (nrow(node_df) - 1)
  heights <- usable_space * node_df$.count / sum(node_df$.count)

  y_top <- top
  tops <- numeric(nrow(node_df))
  bottoms <- numeric(nrow(node_df))
  centers <- numeric(nrow(node_df))

  for (i in seq_len(nrow(node_df))) {
    tops[i] <- y_top
    bottoms[i] <- y_top - heights[i]
    centers[i] <- (tops[i] + bottoms[i]) / 2
    y_top <- bottoms[i] - gap
  }

  node_df %>%
    mutate(y = centers, h = heights, top = tops, bottom = bottoms) %>%
    select(-.count)
}

make_domain_colors <- function(domains) {
  base <- c(
    plant      = "#5DAE61",
    microbe    = "#3F88C5",
    microbiome = "#8A63D2",
    tissue     = "#C95E5E",
    food       = "#D99B2B"
  )
  base[domains]
}


# ==================================================================================================
# 3. Read input
# ==================================================================================================

if (!file.exists(feature_domain_file)) {
  stop("Missing feature-domain summary: ", feature_domain_file,
       "\nRun: .\\.venv310\\Scripts\\python.exe scripts\\summarize_domainmasst_results.py")
}

if (!file.exists(json_long_file)) {
  stop("Missing JSON subclass summary: ", json_long_file,
       "\nRun: .\\.venv310\\Scripts\\python.exe scripts\\summarize_domainmasst_results.py")
}

feature_domain <- read_csv(feature_domain_file, show_col_types = FALSE) %>%
  mutate(
    feature_id = as.character(feature_id),
    domain = str_to_lower(as.character(domain)),
    has_json_match = parse_bool(has_json_match),
    has_counts_rows = parse_bool(has_counts_rows),
    has_domain_evidence = parse_bool(has_domain_evidence),
    counts_rows = suppressWarnings(as.numeric(counts_rows)),
    json_total_matched_size = suppressWarnings(as.numeric(json_total_matched_size))
  )

json_long <- read_csv(json_long_file, show_col_types = FALSE) %>%
  mutate(
    feature_id = as.character(feature_id),
    domain = str_to_lower(as.character(domain)),
    subclass_name = clean_source_label(subclass_name),
    subclass_rank = str_to_lower(as.character(subclass_rank)),
    matched_size = suppressWarnings(as.numeric(matched_size)),
    group_size = suppressWarnings(as.numeric(group_size)),
    occurrence_fraction = suppressWarnings(as.numeric(occurrence_fraction))
  )

if (file.exists(feature_metadata_file)) {
  feature_metadata_raw <- read_csv(feature_metadata_file, show_col_types = FALSE) %>%
    mutate(feature_id = as.character(feature_id))

  compound_col <- case_when(
    "compound_name" %in% names(feature_metadata_raw) ~ "compound_name",
    "compound.name" %in% names(feature_metadata_raw) ~ "compound.name",
    TRUE ~ NA_character_
  )

  mz_col <- if ("mz" %in% names(feature_metadata_raw)) "mz" else NA_character_
  rt_col <- if ("rt" %in% names(feature_metadata_raw)) "rt" else NA_character_

  feature_labels <- feature_metadata_raw %>%
    transmute(
      feature_id,
      compound_name_for_label = if (!is.na(compound_col)) as.character(.data[[compound_col]]) else NA_character_,
      mz_for_label = if (!is.na(mz_col)) suppressWarnings(as.numeric(.data[[mz_col]])) else NA_real_,
      rt_for_label = if (!is.na(rt_col)) suppressWarnings(as.numeric(.data[[rt_col]])) else NA_real_
    ) %>%
    mutate(
      compound_name_for_label = str_trim(compound_name_for_label),
      compound_name_missing = is.na(compound_name_for_label) |
        compound_name_for_label == "" |
        str_to_lower(compound_name_for_label) %in% c("na", "nan", "none", "null"),
      mz_rt_feature_label = paste0(
        if_else(is.na(mz_for_label), "mzNA", format(round(mz_for_label, 4), nsmall = 4, trim = TRUE)),
        "_",
        if_else(is.na(rt_for_label), "rtNA", format(round(rt_for_label, 4), nsmall = 4, trim = TRUE)),
        "_",
        feature_id
      ),
      feature_label_raw = if_else(
        compound_name_missing,
        mz_rt_feature_label,
        paste0(compound_name_for_label, " | ", feature_id)
      ),
      feature_label = short_label(feature_label_raw, max_feature_label_chars)
    ) %>%
    distinct(feature_id, .keep_all = TRUE) %>%
    select(feature_id, compound_name_for_label, mz_for_label, rt_for_label, mz_rt_feature_label, feature_label_raw, feature_label)
} else {
  warning("Feature metadata file not found: ", feature_metadata_file,
          ". Feature labels will fall back to feature IDs.")
  feature_labels <- feature_domain %>%
    distinct(feature_id) %>%
    mutate(
      compound_name_for_label = NA_character_,
      mz_for_label = NA_real_,
      rt_for_label = NA_real_,
      mz_rt_feature_label = feature_id,
      feature_label_raw = paste0("Feature ", feature_id),
      feature_label = short_label(feature_label_raw, max_feature_label_chars)
    )
}

matched_feature_domain <- feature_domain %>%
  filter(
    domain %in% domains_to_show,
    has_domain_evidence
  ) %>%
  distinct(feature_id, domain, .keep_all = TRUE) %>%
  left_join(feature_labels, by = "feature_id")

if (nrow(matched_feature_domain) == 0) {
  stop("No matched feature-domain rows found after filtering. Check domains_to_show and input files.")
}

cat("Matched unique features:", n_distinct(matched_feature_domain$feature_id), "\n")
cat("Matched feature-domain rows:", nrow(matched_feature_domain), "\n")


# ==================================================================================================
# 4. UpSet-style domain overlap plot
# ==================================================================================================

upset_wide <- matched_feature_domain %>%
  distinct(feature_id, domain) %>%
  mutate(value = TRUE) %>%
  pivot_wider(
    names_from = domain,
    values_from = value,
    values_fill = FALSE
  )

for (d in domains_to_show) {
  if (!d %in% names(upset_wide)) {
    upset_wide[[d]] <- FALSE
  }
}

upset_wide <- upset_wide %>%
  select(feature_id, all_of(domains_to_show)) %>%
  mutate(
    domain_set = pmap_chr(
      across(all_of(domains_to_show)),
      function(...) {
        vals <- c(...)
        active <- domains_to_show[as.logical(vals)]
        paste(active, collapse = " + ")
      }
    )
  )

upset_intersections <- upset_wide %>%
  group_by(domain_set) %>%
  summarise(
    n_features = n_distinct(feature_id),
    .groups = "drop"
  ) %>%
  filter(domain_set != "") %>%
  arrange(desc(n_features), domain_set) %>%
  mutate(intersection_id = row_number())

upset_matrix <- upset_intersections %>%
  select(intersection_id, domain_set) %>%
  separate_rows(domain_set, sep = " \\+ ") %>%
  rename(domain = domain_set) %>%
  mutate(
    domain = factor(domain, levels = rev(domains_to_show)),
    present = TRUE
  )

all_matrix_points <- expand_grid(
  intersection_id = upset_intersections$intersection_id,
  domain = factor(domains_to_show, levels = rev(domains_to_show))
) %>%
  left_join(upset_matrix, by = c("intersection_id", "domain")) %>%
  mutate(present = replace_na(present, FALSE))

domain_sizes <- matched_feature_domain %>%
  group_by(domain) %>%
  summarise(n_features = n_distinct(feature_id), .groups = "drop") %>%
  mutate(domain = factor(domain, levels = rev(domains_to_show)))

domain_colors <- make_domain_colors(domains_to_show)

upset_bar <- ggplot(upset_intersections, aes(x = intersection_id, y = n_features)) +
  geom_col(fill = "#2F3437", width = 0.72) +
  geom_text(aes(label = n_features), vjust = -0.35, size = 3.2) +
  scale_x_continuous(breaks = upset_intersections$intersection_id, labels = upset_intersections$intersection_id) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Features") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank()
  )

upset_dot <- ggplot(all_matrix_points, aes(x = intersection_id, y = domain)) +
  geom_point(aes(fill = present), shape = 21, size = 4.0, color = domainmasst_intersection_colour, stroke = 0.5) +
  geom_line(
    data = all_matrix_points %>% filter(present) %>% group_by(intersection_id) %>% filter(n() > 1) %>% ungroup(),
    aes(group = intersection_id),
    color = domainmasst_intersection_colour,
    linewidth = 0.7
  ) +
  scale_fill_manual(values = c(`TRUE` = domainmasst_intersection_colour, `FALSE` = "#D6DADD"), guide = "none") +
  scale_x_continuous(breaks = upset_intersections$intersection_id, labels = upset_intersections$intersection_id) +
  labs(x = "Domain intersection", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(color = "#2F3437")
  )

domain_size_plot <- ggplot(domain_sizes, aes(x = n_features, y = domain, fill = domain)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = n_features), hjust = -0.15, size = 3.1) +
  scale_fill_manual(values = domain_colors, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Domain size", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

upset_plot <- (
  plot_spacer() + upset_bar +
    domain_size_plot + upset_dot
) +
  plot_layout(
    widths = c(0.23, 0.77),
    heights = c(0.45, 0.55)
  ) +
  plot_annotation(
    title = "domainMASST Feature Overlap",
    subtitle = "All matched MetOrigin-unmatched features; combined aggregate domain excluded"
  )

ggsave(file.path(out_dir, "Graph2a_domainMASST_upset.png"), upset_plot,
       width = upset_width, height = upset_height, dpi = plot_dpi)

write_csv(upset_intersections, file.path(out_dir, "Graph2a_domainMASST_upset_intersections.csv"))
write_csv(upset_wide, file.path(out_dir, "Graph2a_domainMASST_feature_domain_matrix.csv"))


# ==================================================================================================
# 5. Graph-1-style Sankey plot: feature -> domain -> phylum/source
# ==================================================================================================

taxonomy_sources <- json_long %>%
  filter(
    domain %in% taxonomy_domains,
    subclass_rank == taxonomy_rank_to_show,
    matched_size > 0,
    !is_generic_source_label(subclass_name)
  ) %>%
  group_by(feature_id, domain, subclass_name) %>%
  summarise(
    evidence_weight = sum(matched_size, na.rm = TRUE),
    max_occurrence_fraction = max(occurrence_fraction, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(source_label = paste0(str_to_title(domain), " phylum: ", subclass_name))

non_taxonomy_sources <- json_long %>%
  filter(
    domain %in% non_taxonomy_domains,
    matched_size > 0,
    !is_generic_source_label(subclass_name)
  ) %>%
  group_by(feature_id, domain, subclass_name) %>%
  summarise(
    evidence_weight = sum(matched_size, na.rm = TRUE),
    max_occurrence_fraction = max(occurrence_fraction, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(source_label = paste0(str_to_title(domain), ": ", subclass_name))

sankey_long <- bind_rows(taxonomy_sources, non_taxonomy_sources) %>%
  semi_join(matched_feature_domain %>% select(feature_id, domain), by = c("feature_id", "domain")) %>%
  left_join(feature_labels %>% select(feature_id, feature_label, feature_label_raw, mz_rt_feature_label), by = "feature_id") %>%
  mutate(
    feature_group_label = "All matched domainMASST features",
    domain_label = str_to_title(domain),
    source_label = clean_source_label(source_label)
  )

# If a matched feature-domain has no usable phylum/source node, keep it visible.
missing_sankey_rows <- matched_feature_domain %>%
  anti_join(sankey_long %>% distinct(feature_id, domain), by = c("feature_id", "domain")) %>%
  transmute(
    feature_id,
    domain,
    subclass_name = paste0(domain, " source not specified"),
    evidence_weight = pmax(json_total_matched_size, counts_rows, 1, na.rm = TRUE),
    max_occurrence_fraction = NA_real_,
    source_label = paste0(str_to_title(domain), ": source not specified"),
    feature_label,
    feature_label_raw,
    mz_rt_feature_label,
    feature_group_label = "All matched domainMASST features",
    domain_label = str_to_title(domain)
  )

sankey_long <- bind_rows(sankey_long, missing_sankey_rows)

if (collapse_other_source_labels) {
  top_sources <- sankey_long %>%
    group_by(domain, source_label) %>%
    summarise(n_features = n_distinct(feature_id), .groups = "drop") %>%
    arrange(domain, desc(n_features), source_label) %>%
    group_by(domain) %>%
    mutate(
      rank_source = row_number(),
      top_n_for_domain = if_else(
        domain %in% names(top_n_source_labels_override),
        as.integer(top_n_source_labels_override[domain]),
        as.integer(top_n_source_labels_per_domain)
      )
    ) %>%
    ungroup() %>%
    filter(rank_source <= top_n_for_domain) %>%
    select(domain, source_label)

  sankey_long <- sankey_long %>%
    left_join(top_sources %>% mutate(is_top_source = TRUE), by = c("domain", "source_label")) %>%
    mutate(
      source_label_plot = if_else(
        !is.na(is_top_source) & is_top_source,
        source_label,
        paste0(str_to_title(domain), ": other source/phylum")
      )
    ) %>%
    select(-is_top_source)
} else {
  sankey_long <- sankey_long %>%
    mutate(source_label_plot = source_label)
}

sankey_edge_domain_source <- sankey_long %>%
  distinct(feature_id, domain, domain_label, source_label_plot) %>%
  group_by(domain, domain_label, source_label_plot) %>%
  summarise(n_features = n_distinct(feature_id), .groups = "drop")

source_order <- sankey_edge_domain_source %>%
  group_by(domain, source_label_plot) %>%
  summarise(n_features = sum(n_features), .groups = "drop") %>%
  arrange(domain, desc(n_features), source_label_plot) %>%
  pull(source_label_plot) %>%
  unique()

sankey_edge_feature_domain <- if (show_within_domain_multi_source_branching) {
  sankey_long %>%
    distinct(feature_id, domain, domain_label, source_label_plot, feature_label, feature_label_raw, mz_rt_feature_label)
} else {
  matched_feature_domain %>%
    distinct(feature_id, domain, feature_label, feature_label_raw, mz_rt_feature_label) %>%
    mutate(domain_label = str_to_title(domain), source_label_plot = NA_character_)
}

sankey_edge_feature_domain <- sankey_edge_feature_domain %>%
  mutate(
    edge_value = 1L
  )

x_feature <- 0.06
x_domain <- 0.45
x_source <- 0.88

# Layer 1: one feature tick per matched feature. This is the key layer for showing
# multi-origin behavior: features with multiple domains have multiple outgoing edges.
feature_nodes <- matched_feature_domain %>%
  distinct(feature_id, feature_label, feature_label_raw, mz_rt_feature_label) %>%
  left_join(
    matched_feature_domain %>%
      group_by(feature_id) %>%
      summarise(n_domains = n_distinct(domain), .groups = "drop"),
    by = "feature_id"
  ) %>%
  arrange(desc(n_domains), suppressWarnings(as.numeric(feature_id)), feature_id) %>%
  mutate(
    node_id = paste0("feature__", feature_id),
    node_label = if (show_feature_labels_in_sankey) feature_label else "",
    n_features = 1L,
    x = x_feature,
    y = if (n() == 1) 0.50 else seq(0.955, 0.045, length.out = n()),
    h = 0.0032,
    top = y + h / 2,
    bottom = y - h / 2,
    fill_colour = "#2F3437",
    anchor_out = x - 0.006
  )

domain_nodes <- sankey_edge_feature_domain %>%
  group_by(domain, domain_label) %>%
  summarise(n_features = n_distinct(feature_id), .groups = "drop") %>%
  mutate(
    domain = factor(domain, levels = domains_to_show),
    node_id = paste0("domain__", as.character(domain)),
    node_label = paste0(domain_label, "\n", n_features),
    x = x_domain
  ) %>%
  arrange(domain) %>%
  allocate_layer_positions(count_col = "n_features", top = 0.92, bottom = 0.08, gap = 0.020) %>%
  mutate(
    fill_colour = domain_colors[as.character(domain)],
    anchor_in = x - 0.040,
    anchor_out = x - 0.030
  )

source_nodes <- sankey_edge_domain_source %>%
  group_by(domain, domain_label, source_label_plot) %>%
  summarise(n_features = max(n_features), .groups = "drop") %>%
  mutate(
    domain = factor(domain, levels = domains_to_show),
    source_label_short = short_label(source_label_plot, 42),
    node_id = paste0("source__", row_number()),
    node_label = paste0(source_label_short, "\n", n_features),
    x = x_source
  ) %>%
  arrange(domain, desc(n_features), source_label_plot) %>%
  allocate_layer_positions(count_col = "n_features", top = 0.955, bottom = 0.045, gap = 0.010) %>%
  mutate(
    fill_colour = domain_colors[as.character(domain)],
    anchor_in = x - 0.048
  )

edges_1_2 <- sankey_edge_feature_domain %>%
  left_join(feature_nodes %>% select(feature_id, from_id = node_id, x0 = anchor_out, y0 = y), by = "feature_id") %>%
  left_join(domain_nodes %>% select(domain, to_id = node_id, x1 = anchor_in, y1 = y), by = "domain") %>%
  mutate(
    edge_type = "feature_to_domain",
    edge_colour = domain_colors[domain],
    line_width = if_else(show_within_domain_multi_source_branching, 0.38, 0.50),
    alpha_value = if_else(show_within_domain_multi_source_branching, 0.12, 0.16)
  ) %>%
  select(edge_type, from_id, to_id, feature_id, domain, source_label_plot, edge_value, edge_colour, line_width, alpha_value, x0, y0, x1, y1)

edges_2_3 <- sankey_edge_domain_source %>%
  left_join(domain_nodes %>% select(domain, from_id = node_id, x0 = anchor_out, y0 = y), by = "domain") %>%
  left_join(
    source_nodes %>% select(domain, source_label_plot, to_id = node_id, x1 = anchor_in, y1 = y),
    by = c("domain", "source_label_plot")
  ) %>%
  mutate(
    edge_type = "domain_to_source",
    edge_value = n_features,
    edge_colour = domain_colors[domain],
    line_width = scale_lw(n_features, max(n_features, na.rm = TRUE), min_lw = 0.6, max_lw = 7.0),
    alpha_value = 0.42
  ) %>%
  mutate(feature_id = NA_character_) %>%
  mutate(source_label_plot = source_label_plot) %>%
  select(edge_type, from_id, to_id, feature_id, domain, source_label_plot, edge_value, edge_colour, line_width, alpha_value, x0, y0, x1, y1)

sankey_edges <- bind_rows(edges_1_2, edges_2_3) %>%
  mutate(edge_id = row_number())

sankey_edges_curve <- sankey_edges %>%
  group_by(edge_id) %>%
  group_modify(~ make_bezier_points(.x$x0, .x$y0, .x$x1, .x$y1) %>%
                 mutate(
                   edge_type = .x$edge_type,
                   domain = .x$domain,
                   edge_value = .x$edge_value,
                   edge_colour = .x$edge_colour,
                   line_width = .x$line_width,
                   alpha_value = .x$alpha_value
                 )) %>%
  ungroup()

feature_tick_nodes <- feature_nodes %>%
  mutate(
    tick_xmin = x - 0.010,
    tick_xmax = x - 0.006,
    tick_ymin = y - 0.0013,
    tick_ymax = y + 0.0013,
    label_x = x - 0.014
  )

domain_bar_nodes <- domain_nodes %>%
  mutate(
    bar_xmin = x - 0.040,
    bar_xmax = x - 0.030,
    label_x = x - 0.022
  )

source_bar_nodes <- source_nodes %>%
  mutate(
    bar_xmin = x - 0.048,
    bar_xmax = x - 0.040,
    label_x = x - 0.036
  )

sankey_plot <- ggplot() +
  geom_path(
    data = sankey_edges_curve,
    aes(x = x, y = y, group = edge_id, colour = edge_colour, linewidth = line_width, alpha = alpha_value),
    lineend = "butt"
  ) +
  scale_colour_identity() +
  scale_linewidth_identity() +
  scale_alpha_identity() +

  geom_rect(
    data = feature_tick_nodes,
    aes(xmin = tick_xmin, xmax = tick_xmax, ymin = tick_ymin, ymax = tick_ymax),
    fill = "#2F3437",
    colour = "#2F3437",
    linewidth = 0
  ) +
  geom_text(
    data = feature_tick_nodes %>% filter(show_feature_labels_in_sankey),
    aes(x = label_x, y = y, label = feature_label),
    hjust = 1,
    size = 1.7,
    lineheight = 0.85,
    colour = "#202426"
  ) +

  geom_rect(
    data = domain_bar_nodes,
    aes(xmin = bar_xmin, xmax = bar_xmax, ymin = bottom, ymax = top, fill = fill_colour),
    colour = NA
  ) +
  geom_text(
    data = domain_bar_nodes,
    aes(x = label_x, y = y, label = node_label),
    hjust = 1,
    size = 3.0,
    fontface = "bold",
    lineheight = 0.92,
    colour = "#202426"
  ) +

  geom_rect(
    data = source_bar_nodes,
    aes(xmin = bar_xmin, xmax = bar_xmax, ymin = bottom, ymax = top, fill = fill_colour),
    colour = NA
  ) +
  geom_text(
    data = source_bar_nodes,
    aes(x = label_x, y = y, label = node_label),
    hjust = 1,
    size = 2.45,
    lineheight = 0.88,
    colour = "#202426"
  ) +
  scale_fill_identity() +

  annotate("text", x = x_feature, y = 1.005, label = "Feature", fontface = "bold", size = 3.7, hjust = 0.5) +
  annotate("text", x = x_domain, y = 1.005, label = "Domain", fontface = "bold", size = 3.7, hjust = 0.5) +
  annotate("text", x = x_source, y = 1.005, label = paste0("Source / ", str_to_title(taxonomy_rank_to_show)), fontface = "bold", size = 3.7, hjust = 0.5) +

  labs(
    title = "domainMASST Feature to Source Evidence",
    subtitle = "Graph-1-style Sankey; plant and microbe phylogeny collapsed to phylum",
    x = NULL,
    y = NULL
  ) +
  coord_cartesian(xlim = c(0.015, 0.965), ylim = c(0.02, 1.02), clip = "off") +
  theme_void(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10.5),
    plot.margin = margin(18, 70, 18, 70)
  )

ggsave(file.path(out_dir, "Graph2b_domainMASST_sankey_phylum_level.png"), sankey_plot,
       width = sankey_width, height = sankey_height, dpi = plot_dpi, bg = "white")

write_csv(sankey_long, file.path(out_dir, "Graph2b_domainMASST_sankey_long_table.csv"))
write_csv(sankey_edge_feature_domain, file.path(out_dir, "Graph2b_domainMASST_edges_feature_to_domain.csv"))
write_csv(sankey_edge_domain_source, file.path(out_dir, "Graph2b_domainMASST_edges_domain_to_source.csv"))
write_csv(sankey_edges, file.path(out_dir, "Graph2b_domainMASST_sankey_edges_graph1_style.csv"))


# ==================================================================================================
# 6. Console summary
# ==================================================================================================

domain_summary <- matched_feature_domain %>%
  group_by(domain) %>%
  summarise(n_features = n_distinct(feature_id), .groups = "drop") %>%
  arrange(desc(n_features))

phylum_summary <- taxonomy_sources %>%
  group_by(domain, subclass_name) %>%
  summarise(
    n_features = n_distinct(feature_id),
    total_matched_size = sum(evidence_weight, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(domain, desc(n_features), desc(total_matched_size))

write_csv(domain_summary, file.path(out_dir, "Graph2_domainMASST_domain_summary.csv"))
write_csv(phylum_summary, file.path(out_dir, "Graph2_domainMASST_phylum_summary.csv"))

cat("\nDomain summary:\n")
print(domain_summary)

cat("\nTop phylum summary:\n")
print(phylum_summary %>% group_by(domain) %>% slice_head(n = 10) %>% ungroup())

cat("\nWrote outputs to:", out_dir, "\n")
