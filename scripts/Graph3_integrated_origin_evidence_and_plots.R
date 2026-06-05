####################################################################################################
# Graph 3: Integrated MetOrigin + domainMASST origin evidence
#
# Outputs
# - Integrated long evidence table
# - Integrated one-row-per-feature summary table
# - Merged source-domain UpSet-style plot
# - MetOrigin-vs-domainMASST normalized-domain agreement matrix
# - Integrated Graph-1-style Sankey: feature -> evidence source -> normalized domain -> detail
#
# Rollback note
# - Separate-source working versions were backed up under:
#   backups/separate_source_ver1/
####################################################################################################


# ==================================================================================================
# 0. Libraries
# ==================================================================================================

required_packages <- c("tidyverse", "ggplot2", "readr", "stringr", "forcats", "scales", "patchwork")
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

graph1_long_file <- "Graph1_MetOrigin_R_outputs/Graph1a_feature_general_subclass_detail_long_table.csv"
graph2_masst_long_file <- "Graph2_domainMASST_R_outputs/Graph2b_domainMASST_sankey_long_table.csv"
metorigin_full_file <- "outputs/metorigin_matches/MAPS_with_MetOriginDB_origin_FULL.csv"
masst_merged_file <- "outputs/domainmasst_unmatched_ms2/MAPS_with_MetOriginDB_plus_domainMASST_origin_evidence.csv"

out_dir <- "Graph3_integrated_origin_R_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

source_order <- c("MetOrigin", "domainMASST")

# Right-most Sankey detail labels are collapsed by source/domain to avoid unreadable plots.
top_n_detail_per_source_domain <- 5
max_upset_intersections_to_plot <- Inf

plot_dpi <- 600
sankey_width <- 17
sankey_height <- 12
upset_width <- 18
upset_height <- 8.5
matrix_width <- 10
matrix_height <- 8


# ==================================================================================================
# 2. Helper functions
# ==================================================================================================

parse_bool <- function(x) {
  str_to_lower(str_trim(as.character(x))) %in% c("true", "yes", "y", "1")
}

short_label <- function(x, max_chars = 48) {
  x <- as.character(x)
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 3), "..."), x)
}

safe_chr <- function(x) {
  out <- as.character(x)
  out[is.na(out)] <- ""
  out
}

clean_joined_values <- function(x, sep = ";") {
  vals <- safe_chr(x)
  vals <- vals[vals != "" & !str_to_lower(vals) %in% c("na", "nan", "none", "null")]
  if (length(vals) == 0) return("")
  paste(sort(unique(vals)), collapse = sep)
}

normalize_domain <- function(source, domain) {
  source <- as.character(source)
  domain <- str_to_lower(as.character(domain))
  case_when(
    source == "MetOrigin" & domain == "bacteria" ~ "Microbe",
    source == "MetOrigin" & domain == "human" ~ "Tissue/Human",
    source == "MetOrigin" & domain == "environment" ~ "Environment",
    source == "MetOrigin" & domain == "drug" ~ "Drug",
    source == "MetOrigin" & domain == "food" ~ "Food",
    source == "MetOrigin" & domain == "plant" ~ "Plant",
    source == "MetOrigin" & domain == "animal" ~ "Animal",
    source == "domainMASST" & domain == "microbe" ~ "Microbe",
    source == "domainMASST" & domain == "microbiome" ~ "Microbiome",
    source == "domainMASST" & domain == "tissue" ~ "Tissue/Human",
    source == "domainMASST" & domain == "plant" ~ "Plant",
    source == "domainMASST" & domain == "food" ~ "Food",
    TRUE ~ str_to_title(domain)
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

scale_lw <- function(value, max_value, min_lw = 0.6, max_lw = 8.0) {
  if (is.na(max_value) || max_value <= 0) return(rep(min_lw, length(value)))
  min_lw + (max_lw - min_lw) * value / max_value
}

allocate_layer_positions <- function(node_df, count_col = "n_features", top = 0.96, bottom = 0.04, gap = 0.008) {
  node_df <- node_df %>% mutate(.count = .data[[count_col]])
  if (nrow(node_df) == 0) {
    return(node_df %>% mutate(y = numeric(0), h = numeric(0), top = numeric(0), bottom = numeric(0)))
  }
  if (nrow(node_df) == 1) {
    return(node_df %>% mutate(y = 0.50, h = 0.16, top = 0.58, bottom = 0.42) %>% select(-.count))
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
  node_df %>% mutate(y = centers, h = heights, top = tops, bottom = bottoms) %>% select(-.count)
}

source_colors <- c(MetOrigin = "#4A4F52", domainMASST = "#7B61C9")
integrated_upset_intersection_colour <- "#2F3437"
domain_colors <- c(
  Food = "#D99B2B",
  Microbe = "#3F88C5",
  Microbiome = "#8A63D2",
  Plant = "#5DAE61",
  Animal = "#A77C5B",
  `Tissue/Human` = "#C95E5E",
  Environment = "#7F8C8D",
  Drug = "#B08AD8"
)

source_domain_colors <- c(
  `MetOrigin::Food` = "#B97712",
  `domainMASST::Food` = "#F0C46A",
  `MetOrigin::Microbe` = "#1F5F91",
  `domainMASST::Microbe` = "#8EC5EA",
  `MetOrigin::Microbiome` = "#6E48B8",
  `domainMASST::Microbiome` = "#C4A8EE",
  `MetOrigin::Plant` = "#2F7D32",
  `domainMASST::Plant` = "#A8D99B",
  `MetOrigin::Animal` = "#7A513A",
  `domainMASST::Animal` = "#D5B49B",
  `MetOrigin::Tissue/Human` = "#9F3030",
  `domainMASST::Tissue/Human` = "#F2A0A0",
  `MetOrigin::Environment` = "#546164",
  `domainMASST::Environment` = "#BAC4C7",
  `MetOrigin::Drug` = "#8054B6",
  `domainMASST::Drug` = "#D1B4EF"
)

get_source_domain_colour <- function(evidence_source, domain) {
  key <- paste0(as.character(evidence_source), "::", as.character(domain))
  out <- unname(source_domain_colors[key])
  fallback <- unname(domain_colors[as.character(domain)])
  ifelse(is.na(out), fallback, out)
}


# ==================================================================================================
# 3. Read inputs
# ==================================================================================================

for (f in c(graph1_long_file, graph2_masst_long_file)) {
  if (!file.exists(f)) stop("Missing required input: ", f)
}

graph1_long <- read_csv(graph1_long_file, show_col_types = FALSE) %>%
  mutate(feature_id = as.character(feature_id))

masst_long <- read_csv(graph2_masst_long_file, show_col_types = FALSE) %>%
  mutate(feature_id = as.character(feature_id))

metorigin_extra <- if (file.exists(metorigin_full_file)) {
  read_csv(metorigin_full_file, show_col_types = FALSE) %>%
    mutate(feature_id = as.character(feature_id)) %>%
    select(
      feature_id,
      any_of(c("confidence_score", "metorigin_match_method", "metorigin_n_database_records", "metorigin_matched"))
    ) %>%
    distinct(feature_id, .keep_all = TRUE)
} else {
  tibble(feature_id = character())
}

masst_extra <- if (file.exists(masst_merged_file)) {
  read_csv(masst_merged_file, show_col_types = FALSE) %>%
    mutate(feature_id = as.character(feature_id)) %>%
    select(
      feature_id,
      any_of(c("masst_searched", "masst_n_matches_rows", "masst_n_dataset_rows", "masst_n_library_rows", "masst_n_count_domain_rows"))
    ) %>%
    distinct(feature_id, .keep_all = TRUE)
} else {
  tibble(feature_id = character())
}


# ==================================================================================================
# 4. Integrated long evidence table
# ==================================================================================================

metorigin_evidence <- graph1_long %>%
  left_join(metorigin_extra, by = "feature_id") %>%
  transmute(
    feature_id,
    feature_label,
    compound_name = safe_chr(compound_name),
    evidence_source = "MetOrigin",
    evidence_type = "identity_database_origin",
    origin_domain_original = safe_chr(general_class),
    origin_domain_normalized = normalize_domain(evidence_source, origin_domain_original),
    origin_subdomain = safe_chr(layer3_plot),
    origin_detail = safe_chr(layer4_plot),
    confidence_level = safe_chr(confidence_level),
    confidence_score = if ("confidence_score" %in% names(.)) suppressWarnings(as.numeric(confidence_score)) else NA_real_,
    metorigin_match_method = if ("metorigin_match_method" %in% names(.)) safe_chr(metorigin_match_method) else "",
    metorigin_n_database_records = if ("metorigin_n_database_records" %in% names(.)) suppressWarnings(as.numeric(metorigin_n_database_records)) else NA_real_,
    masst_metric = NA_real_,
    masst_occurrence_fraction = NA_real_,
    source_rank = "",
    source_label = origin_detail,
    source_label_for_plot = paste0("MetOrigin: ", short_label(origin_subdomain, 36))
  ) %>%
  distinct()

masst_evidence <- masst_long %>%
  left_join(masst_extra, by = "feature_id") %>%
  transmute(
    feature_id,
    feature_label,
    compound_name = "",
    evidence_source = "domainMASST",
    evidence_type = "spectral_source_context",
    origin_domain_original = safe_chr(domain),
    origin_domain_normalized = normalize_domain(evidence_source, origin_domain_original),
    origin_subdomain = safe_chr(subclass_name),
    origin_detail = safe_chr(source_label),
    confidence_level = "",
    confidence_score = NA_real_,
    metorigin_match_method = "",
    metorigin_n_database_records = NA_real_,
    masst_metric = suppressWarnings(as.numeric(evidence_weight)),
    masst_occurrence_fraction = suppressWarnings(as.numeric(max_occurrence_fraction)),
    source_rank = if_else(str_detect(str_to_lower(origin_detail), "phylum"), "phylum/source", "source"),
    source_label = origin_detail,
    source_label_for_plot = safe_chr(source_label_plot)
  ) %>%
  distinct()

evidence_long <- bind_rows(metorigin_evidence, masst_evidence) %>%
  mutate(
    evidence_source = factor(evidence_source, levels = source_order),
    source_domain_token = paste0(as.character(evidence_source), ": ", origin_domain_normalized)
  )

write_csv(evidence_long, file.path(out_dir, "Graph3_integrated_origin_evidence_long.csv"))


# ==================================================================================================
# 5. One-row-per-feature summary and agreement status
# ==================================================================================================

domain_sets <- evidence_long %>%
  distinct(feature_id, evidence_source, origin_domain_normalized) %>%
  group_by(feature_id, evidence_source) %>%
  summarise(
    domains = clean_joined_values(origin_domain_normalized),
    n_domains = n_distinct(origin_domain_normalized),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = evidence_source,
    values_from = c(domains, n_domains),
    values_fill = list(domains = "", n_domains = 0),
    names_sep = "_"
  )

all_feature_labels <- evidence_long %>%
  group_by(feature_id) %>%
  summarise(
    feature_label = first(feature_label[feature_label != ""]),
    .groups = "drop"
  )

feature_summary <- all_feature_labels %>%
  left_join(domain_sets, by = "feature_id") %>%
  mutate(
    domains_MetOrigin = replace_na(domains_MetOrigin, ""),
    domains_domainMASST = replace_na(domains_domainMASST, ""),
    n_domains_MetOrigin = replace_na(n_domains_MetOrigin, 0),
    n_domains_domainMASST = replace_na(n_domains_domainMASST, 0),
    metorigin_domain_set = str_split(domains_MetOrigin, ";"),
    masst_domain_set = str_split(domains_domainMASST, ";"),
    n_overlap_domains = map2_int(metorigin_domain_set, masst_domain_set, ~ length(intersect(.x[.x != ""], .y[.y != ""]))),
    combined_domains = map2_chr(metorigin_domain_set, masst_domain_set, ~ clean_joined_values(c(.x, .y))),
    combined_n_domains = map_int(str_split(combined_domains, ";"), ~ length(.x[.x != ""])),
    origin_agreement_status = case_when(
      n_domains_MetOrigin > 0 & n_domains_domainMASST == 0 ~ "metorigin_only",
      n_domains_MetOrigin == 0 & n_domains_domainMASST > 0 ~ "domainmasst_only",
      n_domains_MetOrigin > 0 & n_domains_domainMASST > 0 & domains_MetOrigin == domains_domainMASST ~ "both_same_normalized_domains",
      n_domains_MetOrigin > 0 & n_domains_domainMASST > 0 & n_overlap_domains > 0 ~ "partial_overlap",
      n_domains_MetOrigin > 0 & n_domains_domainMASST > 0 & n_overlap_domains == 0 ~ "different_domains",
      TRUE ~ "no_origin_evidence"
    )
  ) %>%
  select(
    feature_id,
    feature_label,
    metorigin_domains = domains_MetOrigin,
    metorigin_n_domains = n_domains_MetOrigin,
    domainmasst_domains = domains_domainMASST,
    domainmasst_n_domains = n_domains_domainMASST,
    n_overlap_domains,
    combined_domains,
    combined_n_domains,
    origin_agreement_status
  )

write_csv(feature_summary, file.path(out_dir, "Graph3_integrated_feature_summary.csv"))


# ==================================================================================================
# 6. Merged UpSet-style source-domain plot
# ==================================================================================================

source_domain_order <- evidence_long %>%
  distinct(source_domain_token, evidence_source, origin_domain_normalized) %>%
  mutate(
    evidence_source = factor(evidence_source, levels = source_order),
    source_domain_token = factor(source_domain_token)
  ) %>%
  arrange(evidence_source, origin_domain_normalized) %>%
  pull(source_domain_token) %>%
  unique()

merged_upset_wide <- evidence_long %>%
  distinct(feature_id, source_domain_token) %>%
  mutate(value = TRUE) %>%
  pivot_wider(names_from = source_domain_token, values_from = value, values_fill = FALSE)

for (token in source_domain_order) {
  if (!token %in% names(merged_upset_wide)) merged_upset_wide[[token]] <- FALSE
}

merged_upset_wide <- merged_upset_wide %>%
  select(feature_id, all_of(source_domain_order)) %>%
  mutate(
    source_domain_set = pmap_chr(
      across(all_of(source_domain_order)),
      function(...) {
        vals <- c(...)
        paste(source_domain_order[as.logical(vals)], collapse = " + ")
      }
    )
  )

merged_intersections_complete <- merged_upset_wide %>%
  group_by(source_domain_set) %>%
  summarise(n_features = n_distinct(feature_id), .groups = "drop") %>%
  filter(source_domain_set != "") %>%
  mutate(n_source_domains = str_count(source_domain_set, fixed(" + ")) + 1L) %>%
  mutate(intersection_kind = "observed_feature_set")

missing_singleton_intersections <- tibble(
  source_domain_set = source_domain_order,
  n_features = 0L,
  n_source_domains = 1L,
  intersection_kind = "zero_count_singleton"
) %>%
  anti_join(
    merged_intersections_complete %>% filter(n_source_domains == 1),
    by = "source_domain_set"
  )

merged_intersections_complete <- bind_rows(
  merged_intersections_complete,
  missing_singleton_intersections
) %>%
  arrange(desc(n_features), n_source_domains, source_domain_set) %>%
  mutate(intersection_id = row_number())

if (is.finite(max_upset_intersections_to_plot)) {
  merged_intersections <- merged_intersections_complete %>%
    slice_head(n = max_upset_intersections_to_plot) %>%
    mutate(intersection_id = row_number())
} else {
  merged_intersections <- merged_intersections_complete
}

merged_matrix <- merged_intersections %>%
  select(intersection_id, source_domain_set) %>%
  separate_rows(source_domain_set, sep = " \\+ ") %>%
  rename(source_domain_token = source_domain_set) %>%
  mutate(source_domain_token = factor(source_domain_token, levels = rev(source_domain_order)), present = TRUE)

all_merged_matrix_points <- expand_grid(
  intersection_id = merged_intersections$intersection_id,
  source_domain_token = factor(source_domain_order, levels = rev(source_domain_order))
) %>%
  left_join(merged_matrix, by = c("intersection_id", "source_domain_token")) %>%
  mutate(present = replace_na(present, FALSE))

source_domain_sizes <- evidence_long %>%
  distinct(feature_id, source_domain_token, evidence_source) %>%
  group_by(source_domain_token, evidence_source) %>%
  summarise(n_features = n_distinct(feature_id), .groups = "drop") %>%
  mutate(source_domain_token = factor(source_domain_token, levels = rev(source_domain_order)))

merged_upset_bar <- ggplot(merged_intersections, aes(x = intersection_id, y = n_features)) +
  geom_col(fill = "#2F3437", width = 0.72) +
  geom_text(aes(label = n_features), vjust = -0.35, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Features") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank())

merged_upset_dot <- ggplot(all_merged_matrix_points, aes(x = intersection_id, y = source_domain_token)) +
  geom_point(aes(fill = present), shape = 21, size = 3.2, color = integrated_upset_intersection_colour, stroke = 0.45) +
  geom_line(
    data = all_merged_matrix_points %>% filter(present) %>% group_by(intersection_id) %>% filter(n() > 1) %>% ungroup(),
    aes(group = intersection_id),
    color = integrated_upset_intersection_colour,
    linewidth = 0.6
  ) +
  scale_fill_manual(values = c(`TRUE` = integrated_upset_intersection_colour, `FALSE` = "#D6DADD"), guide = "none") +
  labs(x = "Source-domain intersection", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size = 8))

source_domain_size_plot <- ggplot(source_domain_sizes, aes(x = n_features, y = source_domain_token, fill = evidence_source)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = n_features), hjust = -0.15, size = 2.8) +
  scale_fill_manual(values = source_colors, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Source-domain size", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 8))

merged_upset_plot <- (
  plot_spacer() + merged_upset_bar +
    source_domain_size_plot + merged_upset_dot
) +
  plot_layout(widths = c(0.32, 0.68), heights = c(0.42, 0.58)) +
  plot_annotation(
    title = "Integrated MetOrigin + domainMASST Source-Domain Overlap",
    subtitle = paste0(
      "Complete observed source-domain intersections; ",
      nrow(merged_intersections_complete),
      " feature-set combinations counted"
    )
  )

ggsave(file.path(out_dir, "Graph3a_integrated_source_domain_upset.png"), merged_upset_plot,
       width = upset_width, height = upset_height, dpi = plot_dpi, bg = "white")

write_csv(merged_upset_wide, file.path(out_dir, "Graph3a_integrated_source_domain_feature_matrix.csv"))
write_csv(merged_intersections_complete, file.path(out_dir, "Graph3a_integrated_source_domain_upset_intersections_complete.csv"))
write_csv(merged_intersections, file.path(out_dir, "Graph3a_integrated_source_domain_upset_intersections.csv"))


# ==================================================================================================
# 7. Agreement matrix
# ==================================================================================================

metorigin_domains_long <- evidence_long %>%
  filter(evidence_source == "MetOrigin") %>%
  distinct(feature_id, metorigin_domain = origin_domain_normalized)

masst_domains_long <- evidence_long %>%
  filter(evidence_source == "domainMASST") %>%
  distinct(feature_id, masst_domain = origin_domain_normalized)

agreement_pairs <- metorigin_domains_long %>%
  inner_join(masst_domains_long, by = "feature_id") %>%
  group_by(metorigin_domain, masst_domain) %>%
  summarise(n_features = n_distinct(feature_id), .groups = "drop")

if (nrow(agreement_pairs) == 0) {
  agreement_matrix_plot <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = paste(
        "No feature IDs currently have both MetOrigin and domainMASST evidence.",
        "This is expected for the present workflow because domainMASST was run on MetOrigin-unmatched MS2 features.",
        sep = "\n"
      ),
      size = 4.2,
      lineheight = 1.05,
      color = "#202426"
    ) +
    labs(title = "MetOrigin vs domainMASST Normalized-Domain Agreement") +
    xlim(-1, 1) +
    ylim(-1, 1) +
    theme_void(base_size = 11) +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.title = element_text(face = "bold", size = 14)
    )
} else {
  agreement_matrix_plot <- ggplot(agreement_pairs, aes(x = masst_domain, y = metorigin_domain, fill = n_features)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = n_features), size = 3.4, color = "#202426") +
    scale_fill_gradient(low = "#F3F5F5", high = "#4A4F52", name = "Features") +
    labs(
      title = "MetOrigin vs domainMASST Normalized-Domain Agreement",
      x = "domainMASST normalized domain",
      y = "MetOrigin normalized domain"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1))
}

ggsave(file.path(out_dir, "Graph3b_metorigin_vs_domainmasst_agreement_matrix.png"), agreement_matrix_plot,
       width = matrix_width, height = matrix_height, dpi = plot_dpi, bg = "white")

write_csv(agreement_pairs, file.path(out_dir, "Graph3b_metorigin_vs_domainmasst_agreement_pairs.csv"))


# ==================================================================================================
# 8. Integrated Graph-1-style Sankey
# ==================================================================================================

detail_counts <- evidence_long %>%
  distinct(feature_id, evidence_source, origin_domain_normalized, source_label_for_plot) %>%
  group_by(evidence_source, origin_domain_normalized, source_label_for_plot) %>%
  summarise(n_features = n_distinct(feature_id), .groups = "drop") %>%
  arrange(evidence_source, origin_domain_normalized, desc(n_features), source_label_for_plot) %>%
  group_by(evidence_source, origin_domain_normalized) %>%
  mutate(rank_detail = row_number()) %>%
  ungroup()

top_details <- detail_counts %>%
  filter(rank_detail <= top_n_detail_per_source_domain) %>%
  select(evidence_source, origin_domain_normalized, source_label_for_plot)

sankey_evidence <- evidence_long %>%
  distinct(feature_id, feature_label, evidence_source, origin_domain_normalized, source_label_for_plot) %>%
  left_join(top_details %>% mutate(is_top_detail = TRUE),
            by = c("evidence_source", "origin_domain_normalized", "source_label_for_plot")) %>%
  mutate(
    detail_label_plot = if_else(
      !is.na(is_top_detail) & is_top_detail,
      source_label_for_plot,
      paste0(as.character(evidence_source), ": ", origin_domain_normalized, " other detail")
    )
  ) %>%
  select(-is_top_detail)

x_source <- 0.055
x_feature <- 0.30
x_domain <- 0.57
x_detail <- 0.90

feature_nodes <- sankey_evidence %>%
  distinct(feature_id, feature_label) %>%
  left_join(
    sankey_evidence %>% group_by(feature_id) %>% summarise(n_source_domains = n_distinct(paste(evidence_source, origin_domain_normalized)), .groups = "drop"),
    by = "feature_id"
  ) %>%
  arrange(desc(n_source_domains), suppressWarnings(as.numeric(feature_id)), feature_id) %>%
  mutate(
    node_id = paste0("feature__", feature_id),
    node_label = short_label(feature_label, 42),
    n_features = 1L,
    x = x_feature,
    y = if (n() == 1) 0.50 else seq(0.955, 0.045, length.out = n()),
    h = 0.0028,
    top = y + h / 2,
    bottom = y - h / 2,
    anchor_in = x - 0.010,
    anchor_out = x - 0.006
  )

source_nodes <- sankey_evidence %>%
  group_by(evidence_source) %>%
  summarise(n_features = n_distinct(feature_id), .groups = "drop") %>%
  mutate(
    evidence_source = factor(evidence_source, levels = source_order),
    node_id = paste0("source__", evidence_source),
    node_label = paste0(evidence_source, "\n", n_features),
    x = x_source
  ) %>%
  arrange(evidence_source) %>%
  allocate_layer_positions(count_col = "n_features", top = 0.86, bottom = 0.14, gap = 0.035) %>%
  mutate(
    fill_colour = source_colors[as.character(evidence_source)],
    anchor_out = x - 0.026
  )

domain_nodes <- sankey_evidence %>%
  group_by(evidence_source, origin_domain_normalized) %>%
  summarise(n_features = n_distinct(feature_id), .groups = "drop") %>%
  mutate(
    evidence_source = factor(evidence_source, levels = source_order),
    node_id = paste0("domain__", row_number()),
    node_label = paste0(origin_domain_normalized, "\n", n_features),
    x = x_domain
  ) %>%
  arrange(evidence_source, origin_domain_normalized) %>%
  allocate_layer_positions(count_col = "n_features", top = 0.94, bottom = 0.06, gap = 0.012) %>%
  mutate(
    fill_colour = get_source_domain_colour(evidence_source, origin_domain_normalized),
    anchor_in = x - 0.038,
    anchor_out = x - 0.029
  )

detail_nodes <- sankey_evidence %>%
  group_by(evidence_source, origin_domain_normalized, detail_label_plot) %>%
  summarise(n_features = n_distinct(feature_id), .groups = "drop") %>%
  mutate(
    evidence_source = factor(evidence_source, levels = source_order),
    node_id = paste0("detail__", row_number()),
    node_label = paste0(short_label(detail_label_plot, 38), "\n", n_features),
    x = x_detail
  ) %>%
  arrange(evidence_source, origin_domain_normalized, desc(n_features), detail_label_plot) %>%
  allocate_layer_positions(count_col = "n_features", top = 0.965, bottom = 0.035, gap = 0.006) %>%
  mutate(
    fill_colour = get_source_domain_colour(evidence_source, origin_domain_normalized),
    anchor_in = x - 0.050
  )

edge_source_feature <- sankey_evidence %>%
  distinct(feature_id, evidence_source) %>%
  left_join(source_nodes %>% select(evidence_source, from_id = node_id, x0 = anchor_out, y0 = y), by = "evidence_source") %>%
  left_join(feature_nodes %>% select(feature_id, to_id = node_id, x1 = anchor_in, y1 = y), by = "feature_id") %>%
  mutate(
    edge_type = "source_to_feature",
    edge_value = 1L,
    edge_colour = source_colors[as.character(evidence_source)],
    line_width = 0.42,
    alpha_value = 0.14,
    origin_domain_normalized = ""
  )

edge_feature_domain <- sankey_evidence %>%
  distinct(feature_id, evidence_source, origin_domain_normalized) %>%
  mutate(edge_value = 1L) %>%
  left_join(feature_nodes %>% select(feature_id, from_id = node_id, x0 = anchor_out, y0 = y), by = "feature_id") %>%
  left_join(domain_nodes %>% select(evidence_source, origin_domain_normalized, to_id = node_id, x1 = anchor_in, y1 = y),
            by = c("evidence_source", "origin_domain_normalized")) %>%
  mutate(
    edge_type = "feature_to_domain",
    edge_colour = get_source_domain_colour(evidence_source, origin_domain_normalized),
    line_width = 0.46,
    alpha_value = 0.16
  )

edge_domain_detail <- sankey_evidence %>%
  distinct(feature_id, evidence_source, origin_domain_normalized, detail_label_plot) %>%
  group_by(evidence_source, origin_domain_normalized, detail_label_plot) %>%
  summarise(edge_value = n_distinct(feature_id), .groups = "drop") %>%
  left_join(domain_nodes %>% select(evidence_source, origin_domain_normalized, from_id = node_id, x0 = anchor_out, y0 = y),
            by = c("evidence_source", "origin_domain_normalized")) %>%
  left_join(detail_nodes %>% select(evidence_source, origin_domain_normalized, detail_label_plot, to_id = node_id, x1 = anchor_in, y1 = y),
            by = c("evidence_source", "origin_domain_normalized", "detail_label_plot")) %>%
  mutate(
    edge_type = "domain_to_detail",
    edge_colour = get_source_domain_colour(evidence_source, origin_domain_normalized),
    line_width = scale_lw(edge_value, max(edge_value, na.rm = TRUE), min_lw = 0.6, max_lw = 6.5),
    alpha_value = 0.44,
    feature_id = NA_character_
  )

integrated_edges <- bind_rows(
  edge_source_feature,
  edge_feature_domain,
  edge_domain_detail
) %>%
  select(edge_type, from_id, to_id, feature_id, evidence_source, origin_domain_normalized,
         edge_value, edge_colour, line_width, alpha_value, x0, y0, x1, y1) %>%
  mutate(edge_id = row_number())

integrated_edges_curve <- integrated_edges %>%
  group_by(edge_id) %>%
  group_modify(~ make_bezier_points(.x$x0, .x$y0, .x$x1, .x$y1) %>%
                 mutate(
                   edge_type = .x$edge_type,
                   edge_colour = .x$edge_colour,
                   line_width = .x$line_width,
                   alpha_value = .x$alpha_value
                 )) %>%
  ungroup()

feature_tick_nodes <- feature_nodes %>%
  mutate(tick_xmin = x - 0.010, tick_xmax = x - 0.006, tick_ymin = y - 0.0012, tick_ymax = y + 0.0012, label_x = x - 0.014)

source_bar_nodes <- source_nodes %>% mutate(bar_xmin = x - 0.036, bar_xmax = x - 0.026, label_x = x - 0.020)
domain_bar_nodes <- domain_nodes %>% mutate(bar_xmin = x - 0.038, bar_xmax = x - 0.029, label_x = x - 0.022)
detail_bar_nodes <- detail_nodes %>% mutate(bar_xmin = x - 0.050, bar_xmax = x - 0.042, label_x = x - 0.036)

integrated_sankey_plot <- ggplot() +
  geom_path(
    data = integrated_edges_curve,
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
    data = feature_tick_nodes,
    aes(x = label_x, y = y, label = node_label),
    hjust = 1,
    size = 1.35,
    lineheight = 0.84,
    colour = "#202426"
  ) +
  geom_rect(data = source_bar_nodes, aes(xmin = bar_xmin, xmax = bar_xmax, ymin = bottom, ymax = top, fill = fill_colour), colour = NA) +
  geom_text(data = source_bar_nodes, aes(x = label_x, y = y, label = node_label), hjust = 1, size = 3.1, fontface = "bold", lineheight = 0.92) +
  geom_rect(data = domain_bar_nodes, aes(xmin = bar_xmin, xmax = bar_xmax, ymin = bottom, ymax = top, fill = fill_colour), colour = NA) +
  geom_text(data = domain_bar_nodes, aes(x = label_x, y = y, label = node_label), hjust = 1, size = 2.7, fontface = "bold", lineheight = 0.90) +
  geom_rect(data = detail_bar_nodes, aes(xmin = bar_xmin, xmax = bar_xmax, ymin = bottom, ymax = top, fill = fill_colour), colour = NA) +
  geom_text(data = detail_bar_nodes, aes(x = label_x, y = y, label = node_label), hjust = 0, size = 2.15, lineheight = 0.86) +
  scale_fill_identity() +
  annotate("text", x = x_feature, y = 1.005, label = "Feature", fontface = "bold", size = 3.5, hjust = 0.5) +
  annotate("text", x = x_source, y = 1.005, label = "Evidence Source", fontface = "bold", size = 3.5, hjust = 0.5) +
  annotate("text", x = x_domain, y = 1.005, label = "Normalized Domain", fontface = "bold", size = 3.5, hjust = 0.5) +
  annotate("text", x = x_detail, y = 1.005, label = "Subclass / Detail", fontface = "bold", size = 3.5, hjust = 0.5) +
  labs(
    title = "Integrated Origin Evidence",
    subtitle = "MetOrigin identity/database evidence and domainMASST spectral/source-context evidence kept as separate evidence sources",
    x = NULL,
    y = NULL
  ) +
  coord_cartesian(xlim = c(0.015, 0.985), ylim = c(0.02, 1.02), clip = "off") +
  theme_void(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10.2),
    plot.margin = margin(18, 170, 18, 115)
  )

ggsave(file.path(out_dir, "Graph3c_integrated_origin_sankey.png"), integrated_sankey_plot,
       width = sankey_width, height = sankey_height, dpi = plot_dpi, bg = "white")

write_csv(sankey_evidence, file.path(out_dir, "Graph3c_integrated_sankey_long_table.csv"))
write_csv(integrated_edges, file.path(out_dir, "Graph3c_integrated_sankey_edges_graph1_style.csv"))


# ==================================================================================================
# 9. Console summary
# ==================================================================================================

source_domain_summary <- evidence_long %>%
  distinct(feature_id, evidence_source, origin_domain_normalized) %>%
  group_by(evidence_source, origin_domain_normalized) %>%
  summarise(n_features = n_distinct(feature_id), .groups = "drop") %>%
  arrange(evidence_source, desc(n_features))

agreement_status_summary <- feature_summary %>%
  group_by(origin_agreement_status) %>%
  summarise(n_features = n_distinct(feature_id), .groups = "drop") %>%
  arrange(desc(n_features))

write_csv(source_domain_summary, file.path(out_dir, "Graph3_integrated_source_domain_summary.csv"))
write_csv(agreement_status_summary, file.path(out_dir, "Graph3_integrated_agreement_status_summary.csv"))

cat("\nIntegrated source-domain summary:\n")
print(source_domain_summary)

cat("\nAgreement status summary:\n")
print(agreement_status_summary)

cat("\nWrote Graph 3 integrated outputs to:", out_dir, "\n")
