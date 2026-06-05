
####################################################################################################
# Graph 1a: Four-layer Sankey-style origin plot for MAPS + MetOriginDB
#
# Final design (multi-origin aware):
#   Layer 1: matched feature / compound
#   Layer 2: general source class
#   Layer 3: category-specific subclass
#   Layer 4: most specific source detail / next layer
#
# IMPORTANT INTERPRETATION
# - One feature can contribute to multiple origin classes (e.g., Food AND Drug).
# - Therefore, one feature can branch into multiple flows.
# - Counts in classes / subclasses are NOT mutually exclusive.
# - Node bar length = number of unique features in that node.
# - Edge width = number of unique features flowing between two adjacent nodes.
#
# INPUT EXPECTED
#   MAPS_with_MetOriginDB_origin_EXCEL_PREVIEW_truncated_details.csv
#
# CORE METORIGIN COLUMNS USED
#   feature_id
#   compound_name
#   metorigin_matched
#   confidence_level
#   from_food
#   from_bacteria
#   from_drug
#   from_environment
#   from_human
#   from_plant
#   from_animal
#   from_which_food
#   bacteria_genus / from_which_bacteria
#   bacteria_species (if available)
#   from_which_drug
#   from_which_environment
#   from_which_part
#   from_which_plant
#   from_which_animal
#
# OPTIONAL OUTPUTS
# - PNG
# - PDF
# - long table
# - edge tables
####################################################################################################


# ==================================================================================================
# 0. Libraries
# ==================================================================================================

required_packages <- c("tidyverse", "ggplot2", "stringr", "readr", "forcats", "scales")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(tidyverse)
library(ggplot2)
library(stringr)
library(readr)
library(forcats)
library(scales)


# ==================================================================================================
# 1. User settings
# ==================================================================================================

input_file <- "C:\\Users\\siyao4\\OneDrive - The University of Melbourne\\Documents\\MAPS_MetOrigin_Matching_Project\\MAPS_MetOrigin_Matching_Project\\outputs\\metorigin_matches\\MAPS_with_MetOriginDB_origin_EXCEL_PREVIEW_truncated_details.csv"

out_dir <- "Graph1_MetOrigin_R_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Number of Sankey layers to show
# -----------------------------
# 2 = Feature -> General source class
# 3 = Feature -> General source class -> Subclass
# 4 = Feature -> General source class -> Subclass -> Detail
n_layers_to_show <- 3
if (!n_layers_to_show %in% c(2, 3, 4)) {
  stop("n_layers_to_show must be 2, 3, or 4.")
}

# -----------------------------
# Filtering
# -----------------------------
use_metorigin_matched_filter <- TRUE
use_confidence_filter <- TRUE
max_confidence_level <- 3

# -----------------------------
# Feature display
# -----------------------------
max_features_to_show <- Inf
# max_features_to_show <- 180

feature_ranking_method <- "n_classes"  # "n_classes" or "as_is"

use_compound_name_in_feature_label <- TRUE
truncate_feature_labels <- TRUE
feature_label_max_chars <- 48

truncate_right_labels <- TRUE
right_label_max_chars <- 42

# -----------------------------
# Class inclusion
# -----------------------------
major_classes_to_show <- c(
  "Food",
  "Bacteria",
  "Drug",
  "Environment",
  "Human",
  "Plant",
  "Animal"
)

# For a cleaner main-figure version you can try:
# major_classes_to_show <- c("Food", "Bacteria", "Drug", "Environment")

# -----------------------------
# Node collapsing
# -----------------------------
top_n_layer3_per_class <- 8
top_n_layer4_per_layer3 <- 5

collapse_other_layer3 <- TRUE
collapse_other_layer4 <- TRUE

# -----------------------------
# Plot size
# -----------------------------
plot_width <- 16
plot_height <- 10
plot_dpi <- 300

# -----------------------------
# Layout
# -----------------------------
x_feature <- 0.055
x_general <- 0.34
x_layer3  <- 0.62
x_layer4  <- 0.90

# -----------------------------
# edge width
# -----------------------------
feature_to_general_alpha <- 0.18
general_to_layer3_alpha <- 0.50
layer3_to_layer4_alpha  <- 0.65

# Re-space columns depending on whether you choose 2, 3, or 4 layers.
# The underlying full long table is still generated; this only changes what is drawn.
if (n_layers_to_show == 2) {
  x_general <- 0.82
}
if (n_layers_to_show == 3) {
  x_general <- 0.45
  x_layer3  <- 0.86
}
if (n_layers_to_show == 4) {
  x_general <- 0.34
  x_layer3  <- 0.62
  x_layer4  <- 0.90
}

# Used for vertical allocation
top_y <- 0.97
bottom_y <- 0.03

# -----------------------------
# Visual style
# -----------------------------
feature_to_general_alpha <- 0.10
general_to_layer3_alpha <- 0.38
layer3_to_layer4_alpha  <- 0.55

class_colors <- c(
  Food        = "#E5C774",
  Bacteria    = "#6AA9D6",
  Drug        = "#B08AD8",
  Environment = "#97A3A8",
  Human       = "#D88A86",
  Plant       = "#84BD84",
  Animal      = "#B99A7B"
)

yes_values <- c("yes", "y", "true", "t", "1")

# -----------------------------
# Output paths
# -----------------------------
graph1a_long_csv <- file.path(out_dir, "Graph1a_feature_general_subclass_detail_long_table.csv")
graph1a_cleaned_input_csv <- file.path(out_dir, "Graph1a_cleaned_input.csv")

graph1a_feature_general_csv <- file.path(out_dir, "Graph1a_edges_feature_to_general.csv")
graph1a_general_layer3_csv <- file.path(out_dir, "Graph1a_edges_general_to_subclass.csv")
graph1a_layer3_layer4_csv  <- file.path(out_dir, "Graph1a_edges_subclass_to_detail.csv")


# ==================================================================================================
# 2. Helper functions
# ==================================================================================================

clean_names_simple <- function(x) {
  x %>%
    str_trim() %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_remove("^_") %>%
    str_remove("_$")
}

is_yes <- function(x) {
  str_to_lower(str_trim(as.character(x))) %in% yes_values
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

safe_write_csv <- function(x, path) {
  tryCatch(
    write_csv(x, path),
    error = function(e) warning("Could not write CSV '", path, "': ", conditionMessage(e), call. = FALSE)
  )
}

safe_ggsave <- function(filename, plot, ...) {
  had_warning <- FALSE
  err <- tryCatch(
    {
      withCallingHandlers(
        ggsave(filename, plot, ...),
        warning = function(w) {
          had_warning <<- TRUE
          warning(conditionMessage(w), call. = FALSE)
          invokeRestart("muffleWarning")
        }
      )
      NULL
    },
    error = function(e) e
  )

  if (!is.null(err) || had_warning) {
    stamp <- format(Sys.time(), "%H%M%S")
    short_base <- tools::file_path_sans_ext(basename(filename)) %>%
      str_replace("^Graph1a_", "G1a_") %>%
      str_replace("_sankey$", "")
    fallback <- file.path(
      dirname(filename),
      paste0(short_base, "_aligned_", stamp, ".", tools::file_ext(filename))
    )
    warning(
      "Could not write plot '", filename, "'. Writing fallback copy '", fallback, "'.",
      call. = FALSE
    )
    tryCatch(
      ggsave(fallback, plot, ...),
      error = function(e) warning("Could not write fallback plot '", fallback, "': ", conditionMessage(e), call. = FALSE)
    )
  }
}

short_label <- function(x, max_chars = 40) {
  x <- as.character(x)
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 3), "..."), x)
}

split_origin_detail <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x)]
  
  if (length(x) == 0) return(character(0))
  
  x <- str_replace_all(x, fixed("{"), "|")
  x <- str_replace_all(x, fixed("}"), "|")
  x <- str_replace_all(x, ";", "|")
  x <- str_replace_all(x, "\\.\\.\\.$", "")
  
  parts <- str_split(x, "\\|") %>%
    unlist() %>%
    str_trim()
  
  parts <- parts[
    !is.na(parts) &
      parts != "" &
      !str_to_lower(parts) %in% c(
        "unknown", "na", "nan", "none", "null",
        "other", "others", "not available", "detail not specified"
      )
  ]
  
  unique(parts)
}

make_feature_label <- function(compound_name, feature_id) {
  compound_name <- as.character(compound_name)
  feature_id <- as.character(feature_id)
  
  bad <- is.na(compound_name) | compound_name == "" | str_to_lower(compound_name) == "na"
  compound_name[bad] <- paste0("Feature_", feature_id[bad])
  
  if (use_compound_name_in_feature_label) {
    label <- paste0(compound_name, " | ", feature_id)
  } else {
    label <- feature_id
  }
  
  if (truncate_feature_labels) {
    label <- short_label(label, feature_label_max_chars)
  }
  
  label
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

scale_lw <- function(value, max_value, min_lw = 1.0, max_lw = 12.0) {
  min_lw + (max_lw - min_lw) * value / max_value
}

allocate_layer_positions <- function(node_df, count_col = "n_features", top = 0.97, bottom = 0.03, gap = 0.006) {
  node_df <- node_df %>% mutate(.count = .data[[count_col]])
  
  if (nrow(node_df) == 0) {
    return(node_df %>% mutate(y = numeric(0), h = numeric(0), top = numeric(0), bottom = numeric(0)))
  }
  
  if (nrow(node_df) == 1) {
    return(node_df %>% mutate(y = 0.50, h = 0.10, top = 0.55, bottom = 0.45) %>% select(-.count))
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
    mutate(
      y = centers,
      h = heights,
      top = tops,
      bottom = bottoms
    ) %>%
    select(-.count)
}


# ==================================================================================================
# 3. Category-specific next-layer logic
# ==================================================================================================
# Goal:
#   General class -> category-specific subclass -> most specific source detail
#
# Layer 3 and Layer 4 are defined differently for different general classes.

classify_food_group <- function(x) {
  z <- str_to_lower(x)
  
  case_when(
    str_detect(z, "coffee|tea|cocoa|chocolate|beverage|drink|espresso|latte") ~ "Beverage",
    str_detect(z, "milk|cheese|yogurt|yoghurt|dairy") ~ "Dairy",
    str_detect(z, "fish|marine|shellfish|seafood") ~ "Seafood",
    str_detect(z, "ferment|kimchi|sauerkraut|pickle|kombucha|miso|tempeh") ~ "Fermented food",
    str_detect(z, "fruit|berry|citrus|apple|banana|grape|pear|orange|lemon") ~ "Fruit",
    str_detect(z, "vegetable|leaf|spinach|brassica|cabbage|lettuce|broccoli|kale") ~ "Vegetable",
    str_detect(z, "grain|wheat|rice|barley|oat|cereal|bread") ~ "Grain / cereal",
    str_detect(z, "beef|veal|cattle|bison|buffalo|deer|chicken|duck|pigeon|pork|meat") ~ "Meat / animal-derived food",
    TRUE ~ "Other food sources"
  )
}

classify_drug_group <- function(x) {
  z <- str_to_lower(x)
  
  case_when(
    str_detect(z, "antibiotic|penicillin|cephalosporin|macrolide|tetracycline") ~ "Antibiotic",
    str_detect(z, "analgesic|paracetamol|acetaminophen|ibuprofen|aspirin|naproxen") ~ "Analgesic / anti-inflammatory",
    str_detect(z, "antifungal|azole|fluconazole|ketoconazole") ~ "Antifungal",
    str_detect(z, "metformin|statin|antihypertensive|lisinopril|amlodipine") ~ "Common medication",
    TRUE ~ "Other drug sources"
  )
}

classify_environment_group <- function(x) {
  z <- str_to_lower(x)
  
  case_when(
    str_detect(z, "pesticide|insecticide|herbicide|fungicide|imidacloprid|propamocarb") ~ "Pesticide / agrochemical",
    str_detect(z, "plastic|phthalate|bisphenol|polymer") ~ "Plastic-related",
    str_detect(z, "solvent|industrial|dye|detergent") ~ "Industrial chemical",
    str_detect(z, "pollutant|contaminant|environment") ~ "Environmental contaminant",
    TRUE ~ "Other environmental sources"
  )
}

classify_human_group <- function(x) {
  z <- str_to_lower(x)
  
  case_when(
    str_detect(z, "blood|plasma|serum") ~ "Blood / circulation",
    str_detect(z, "urine|kidney|renal") ~ "Urinary system",
    str_detect(z, "feces|faeces|colon|gut|intestin|stool") ~ "Gastrointestinal",
    str_detect(z, "saliva|oral") ~ "Oral / saliva",
    str_detect(z, "brain|csf|cerebrospinal|neural") ~ "Nervous system",
    str_detect(z, "placenta|pregnancy|uterus|ovar|testi") ~ "Reproductive / pregnancy",
    TRUE ~ "Other human sources"
  )
}

classify_plant_group <- function(x) {
  # Use first word as a genus-like grouping
  out <- str_extract(x, "^[A-Za-z]+")
  ifelse(is.na(out) | out == "", "Other plant sources", out)
}

classify_animal_group <- function(x) {
  z <- str_to_lower(x)
  
  case_when(
    str_detect(z, "fish|salmon|tuna|marine|shellfish") ~ "Fish / marine animal",
    str_detect(z, "cow|cattle|mammal|milk|dairy|deer|bison|buffalo") ~ "Mammal",
    str_detect(z, "chicken|duck|bird|avian|pigeon") ~ "Bird",
    TRUE ~ "Other animal sources"
  )
}

make_layer3_layer4_for_class <- function(data, general_class_i, detail_col_i) {
  # Returns a table with:
  # feature_id, feature_label, compound_name, general_class, layer3, layer4
  
  if (nrow(data) == 0) return(tibble())
  
  # Special handling for bacteria:
  #   layer3 = genus
  #   layer4 = species if available; otherwise reuse detail
  if (general_class_i == "Bacteria") {
    genus_col <- if ("bacteria_genus" %in% names(data)) "bacteria_genus" else detail_col_i
    species_col <- if ("bacteria_species" %in% names(data)) "bacteria_species" else NA_character_
    
    genus_long <- data %>%
      mutate(genus_list = map(.data[[genus_col]], split_origin_detail)) %>%
      select(feature_id, feature_label, compound_name, any_of(c("confidence_level", "annotation_type", "cid", "smiles")), genus_list) %>%
      unnest_longer(genus_list, values_to = "layer3", keep_empty = TRUE) %>%
      mutate(
        layer3 = if_else(is.na(layer3) | layer3 == "", "Bacteria: genus not specified", layer3)
      )
    
    if (!is.na(species_col) && species_col %in% names(data)) {
      species_long <- data %>%
        mutate(species_list = map(.data[[species_col]], split_origin_detail)) %>%
        select(feature_id, species_list) %>%
        unnest_longer(species_list, values_to = "species_value", keep_empty = TRUE) %>%
        mutate(
          species_value = if_else(is.na(species_value) | species_value == "", "Bacteria: species not specified", species_value)
        )
      
      out <- genus_long %>%
        left_join(species_long, by = "feature_id") %>%
        mutate(
          general_class = general_class_i,
          layer4 = species_value
        )
    } else {
      out <- genus_long %>%
        mutate(
          general_class = general_class_i,
          layer4 = layer3
        )
    }
    
    return(out %>% select(feature_id, feature_label, compound_name, any_of(c("confidence_level", "annotation_type", "cid", "smiles")), general_class, layer3, layer4))
  }
  
  if (!detail_col_i %in% names(data)) {
    return(
      data %>%
        transmute(
          feature_id,
          feature_label,
          compound_name,
          confidence_level = if ("confidence_level" %in% names(data)) confidence_level else NA,
          annotation_type  = if ("annotation_type" %in% names(data)) annotation_type else NA,
          cid              = if ("cid" %in% names(data)) cid else NA,
          smiles           = if ("smiles" %in% names(data)) smiles else NA,
          general_class    = general_class_i,
          layer3           = paste0(general_class_i, ": detail not specified"),
          layer4           = paste0(general_class_i, ": detail not specified")
        )
    )
  }
  
  detail_long <- data %>%
    mutate(detail_list = map(.data[[detail_col_i]], split_origin_detail)) %>%
    select(feature_id, feature_label, compound_name, any_of(c("confidence_level", "annotation_type", "cid", "smiles")), detail_list) %>%
    unnest_longer(detail_list, values_to = "detail_value", keep_empty = TRUE) %>%
    mutate(
      detail_value = if_else(is.na(detail_value) | detail_value == "", paste0(general_class_i, ": detail not specified"), detail_value),
      general_class = general_class_i,
      layer4 = detail_value,
      layer3 = case_when(
        general_class_i == "Food"        ~ classify_food_group(detail_value),
        general_class_i == "Drug"        ~ classify_drug_group(detail_value),
        general_class_i == "Environment" ~ classify_environment_group(detail_value),
        general_class_i == "Human"       ~ classify_human_group(detail_value),
        general_class_i == "Plant"       ~ classify_plant_group(detail_value),
        general_class_i == "Animal"      ~ classify_animal_group(detail_value),
        TRUE ~ paste0(general_class_i, ": subclass")
      )
    ) %>%
    select(feature_id, feature_label, compound_name, any_of(c("confidence_level", "annotation_type", "cid", "smiles")), general_class, layer3, layer4)
  
  detail_long
}


# ==================================================================================================
# 4. Read and clean input
# ==================================================================================================

df <- read_csv(input_file, show_col_types = FALSE)
names(df) <- clean_names_simple(names(df))

cat("Input rows:", nrow(df), "\n")
cat("Input columns:", ncol(df), "\n")

if (!"feature_id" %in% names(df)) {
  warning("feature_id column not found. Using row number as feature_id.")
  df <- df %>% mutate(feature_id = as.character(row_number()))
} else {
  df <- df %>% mutate(feature_id = as.character(feature_id))
}

if (!"compound_name" %in% names(df)) {
  df <- df %>% mutate(compound_name = paste0("Feature_", feature_id))
}

if (use_metorigin_matched_filter && "metorigin_matched" %in% names(df)) {
  df <- df %>%
    filter(str_to_lower(as.character(metorigin_matched)) %in% c("true", "yes", "1"))
}

if (use_confidence_filter && "confidence_level" %in% names(df)) {
  df <- df %>%
    mutate(confidence_level_numeric = safe_numeric(confidence_level)) %>%
    filter(is.na(confidence_level_numeric) | confidence_level_numeric <= max_confidence_level)
}

df <- df %>%
  mutate(feature_label = make_feature_label(compound_name, feature_id))

safe_write_csv(df, graph1a_cleaned_input_csv)

cat("Rows after filtering:", nrow(df), "\n")


# ==================================================================================================
# 5. Build long table: feature -> general class -> layer3 -> layer4
# ==================================================================================================

bacteria_detail_col <- if ("bacteria_genus" %in% names(df)) "bacteria_genus" else "from_which_bacteria"

source_column_map <- tribble(
  ~general_class, ~flag_col,           ~detail_col,
  "Food",         "from_food",         "from_which_food",
  "Bacteria",     "from_bacteria",     bacteria_detail_col,
  "Drug",         "from_drug",         "from_which_drug",
  "Environment",  "from_environment",  "from_which_environment",
  "Human",        "from_human",        "from_which_part",
  "Plant",        "from_plant",        "from_which_plant",
  "Animal",       "from_animal",       "from_which_animal"
) %>%
  filter(general_class %in% major_classes_to_show) %>%
  filter(flag_col %in% names(df)) %>%
  mutate(general_class = factor(general_class, levels = major_classes_to_show)) %>%
  arrange(general_class)

cat("General classes used:\n")
print(source_column_map)

plot_long_list <- list()

for (i in seq_len(nrow(source_column_map))) {
  general_class_i <- as.character(source_column_map$general_class[i])
  flag_col_i <- source_column_map$flag_col[i]
  detail_col_i <- source_column_map$detail_col[i]
  
  tmp <- df %>% filter(is_yes(.data[[flag_col_i]]))
  if (nrow(tmp) == 0) next
  
  plot_long_list[[general_class_i]] <- make_layer3_layer4_for_class(
    data = tmp,
    general_class_i = general_class_i,
    detail_col_i = detail_col_i
  )
}

plot_long <- bind_rows(plot_long_list) %>%
  distinct(feature_id, general_class, layer3, layer4, .keep_all = TRUE)

if (nrow(plot_long) == 0) {
  stop("No feature-general-layer3-layer4 records were created. Check the source flag columns.")
}

# Optional feature filtering
feature_scores <- plot_long %>%
  group_by(feature_id, feature_label) %>%
  summarise(
    n_classes = n_distinct(general_class),
    n_layer3 = n_distinct(layer3),
    n_layer4 = n_distinct(layer4),
    .groups = "drop"
  )

if (is.finite(max_features_to_show) && nrow(feature_scores) > max_features_to_show) {
  if (feature_ranking_method == "n_classes") {
    keep_features <- feature_scores %>%
      arrange(desc(n_classes), desc(n_layer3), desc(n_layer4), feature_label) %>%
      slice_head(n = max_features_to_show) %>%
      pull(feature_id)
  } else {
    keep_features <- feature_scores %>%
      slice_head(n = max_features_to_show) %>%
      pull(feature_id)
  }
  
  plot_long <- plot_long %>% filter(feature_id %in% keep_features)
}


# ==================================================================================================
# 6. Collapse minor Layer 3 and Layer 4 nodes
# ==================================================================================================

layer3_counts <- plot_long %>%
  group_by(general_class, layer3) %>%
  summarise(n_features = n_distinct(feature_id), .groups = "drop") %>%
  mutate(general_class = factor(general_class, levels = major_classes_to_show)) %>%
  arrange(general_class, desc(n_features)) %>%
  group_by(general_class) %>%
  mutate(rank_layer3 = row_number()) %>%
  ungroup()

if (collapse_other_layer3) {
  top_layer3 <- layer3_counts %>%
    filter(rank_layer3 <= top_n_layer3_per_class) %>%
    select(general_class, layer3)
  
  plot_long2 <- plot_long %>%
    left_join(top_layer3 %>% mutate(is_top_layer3 = TRUE), by = c("general_class", "layer3")) %>%
    mutate(
      layer3_plot = if_else(!is.na(is_top_layer3) & is_top_layer3, layer3,
        paste0(general_class, ": other subclasses")
      )
    ) %>%
    select(-is_top_layer3)
} else {
  plot_long2 <- plot_long %>% mutate(layer3_plot = layer3)
}

layer4_counts <- plot_long2 %>%
  group_by(general_class, layer3_plot, layer4) %>%
  summarise(n_features = n_distinct(feature_id), .groups = "drop") %>%
  arrange(general_class, layer3_plot, desc(n_features)) %>%
  group_by(general_class, layer3_plot) %>%
  mutate(rank_layer4 = row_number()) %>%
  ungroup()

if (collapse_other_layer4) {
  top_layer4 <- layer4_counts %>%
    filter(rank_layer4 <= top_n_layer4_per_layer3) %>%
    select(general_class, layer3_plot, layer4)
  
  plot_long2 <- plot_long2 %>%
    left_join(top_layer4 %>% mutate(is_top_layer4 = TRUE), by = c("general_class", "layer3_plot", "layer4")) %>%
    mutate(
      layer4_plot = if_else(!is.na(is_top_layer4) & is_top_layer4, layer4,
        paste0(layer3_plot, ": other details")
      )
    ) %>%
    select(-is_top_layer4)
} else {
  plot_long2 <- plot_long2 %>% mutate(layer4_plot = layer4)
}

plot_long2 <- plot_long2 %>%
  distinct(feature_id, feature_label, general_class, layer3_plot, layer4_plot, .keep_all = TRUE)

safe_write_csv(plot_long2, graph1a_long_csv)


# ==================================================================================================
# 7. Edge tables
# ==================================================================================================

edge_feature_general <- plot_long2 %>%
  distinct(feature_id, feature_label, general_class) %>%
  mutate(value = 1L)

edge_general_layer3 <- plot_long2 %>%
  group_by(general_class, layer3_plot) %>%
  summarise(n_features = n_distinct(feature_id), .groups = "drop") %>%
  mutate(general_class = factor(general_class, levels = major_classes_to_show)) %>%
  arrange(general_class, desc(n_features))

edge_layer3_layer4 <- plot_long2 %>%
  group_by(general_class, layer3_plot, layer4_plot) %>%
  summarise(n_features = n_distinct(feature_id), .groups = "drop") %>%
  mutate(general_class = factor(general_class, levels = major_classes_to_show)) %>%
  arrange(general_class, layer3_plot, desc(n_features))

safe_write_csv(edge_feature_general, graph1a_feature_general_csv)
safe_write_csv(edge_general_layer3, graph1a_general_layer3_csv)
safe_write_csv(edge_layer3_layer4, graph1a_layer3_layer4_csv)

cat("Unique features in plot:", n_distinct(plot_long2$feature_id), "\n")
cat("Layer 3 nodes:", n_distinct(plot_long2$layer3_plot), "\n")
cat("Layer 4 nodes:", n_distinct(plot_long2$layer4_plot), "\n")


# ==================================================================================================
# 8. Node tables
# ==================================================================================================

# Layer 1 feature nodes
feature_nodes <- edge_feature_general %>%
  distinct(feature_id, feature_label) %>%
  arrange(feature_label) %>%
  mutate(
    node_id = paste0("feature__", feature_id),
    node_label = feature_label,
    n_features = 1L,
    x = x_feature,
    y = if (n() == 1) 0.50 else seq(top_y, bottom_y, length.out = n()),
    h = 0.0025,
    top = y + h / 2,
    bottom = y - h / 2,
    anchor_out = x - 0.006
  )

# Layer 2 general nodes
general_nodes <- edge_feature_general %>%
  group_by(general_class) %>%
  summarise(n_features = n_distinct(feature_id), .groups = "drop") %>%
  mutate(
    general_class = factor(general_class, levels = major_classes_to_show)
  ) %>%
  arrange(general_class) %>%
  mutate(
    node_id = paste0("general__", as.character(general_class)),
    node_label = paste0(as.character(general_class), "\n", n_features, " features"),
    x = x_general
  ) %>%
  allocate_layer_positions(count_col = "n_features", top = 0.94, bottom = 0.06, gap = 0.020) %>%
  mutate(
    anchor_in = x - 0.040,
    anchor_out = x - 0.030
  )

# Layer 3 nodes
layer3_nodes <- edge_general_layer3 %>%
  group_by(general_class, layer3_plot) %>%
  summarise(n_features = max(n_features), .groups = "drop") %>%
  mutate(
    general_class = factor(general_class, levels = major_classes_to_show),
    layer3_label = if (truncate_right_labels) short_label(layer3_plot, right_label_max_chars) else layer3_plot
  ) %>%
  arrange(general_class, desc(n_features), layer3_plot) %>%
  mutate(
    node_id = paste0("layer3__", row_number()),
    node_label = paste0(layer3_label, "\n", n_features),
    x = x_layer3
  ) %>%
  allocate_layer_positions(count_col = "n_features", top = 0.965, bottom = 0.035, gap = 0.009) %>%
  mutate(
    anchor_in = x - 0.040,
    anchor_out = x - 0.031
  )

# Layer 4 nodes
layer4_nodes <- edge_layer3_layer4 %>%
  group_by(general_class, layer3_plot, layer4_plot) %>%
  summarise(n_features = max(n_features), .groups = "drop") %>%
  mutate(
    general_class = factor(general_class, levels = major_classes_to_show),
    layer4_label = if (truncate_right_labels) short_label(layer4_plot, right_label_max_chars) else layer4_plot
  ) %>%
  arrange(general_class, layer3_plot, desc(n_features), layer4_plot) %>%
  mutate(
    node_id = paste0("layer4__", row_number()),
    node_label = paste0(layer4_label, "\n", n_features),
    x = x_layer4
  ) %>%
  allocate_layer_positions(count_col = "n_features", top = 0.970, bottom = 0.030, gap = 0.007) %>%
  mutate(anchor_in = x - 0.050)



# ==================================================================================================
# 9. Edge curves, with optional 2/3/4-layer display
# ==================================================================================================

show_layer3 <- n_layers_to_show >= 3
show_layer4 <- n_layers_to_show >= 4

# Feature -> General
edges_1_2 <- edge_feature_general %>%
  left_join(feature_nodes %>% select(feature_id, from_id = node_id, x0 = anchor_out, y0 = y), by = "feature_id") %>%
  mutate(general_class_chr = as.character(general_class)) %>%
  left_join(
    general_nodes %>%
      mutate(general_class_chr = as.character(general_class)) %>%
      select(general_class_chr, to_id = node_id, x1 = anchor_in, y1 = y),
    by = "general_class_chr"
  ) %>%
  mutate(
    edge_type = "feature_to_general",
    edge_value = 1L,
    edge_colour = class_colors[general_class_chr],
    line_width = 0.55,
    alpha_value = feature_to_general_alpha
  ) %>%
  select(edge_type, from_id, to_id, general_class = general_class_chr, edge_value, edge_colour, line_width, alpha_value, x0, y0, x1, y1)

# General -> Layer 3
edges_2_3 <- edge_general_layer3 %>%
  mutate(general_class_chr = as.character(general_class)) %>%
  left_join(
    general_nodes %>%
      mutate(general_class_chr = as.character(general_class)) %>%
      select(general_class_chr, from_id = node_id, x0 = anchor_out, y0 = y),
    by = "general_class_chr"
  ) %>%
  left_join(
    layer3_nodes %>%
      mutate(general_class_chr = as.character(general_class)) %>%
      select(general_class_chr, layer3_plot, to_id = node_id, x1 = anchor_in, y1 = y),
    by = c("general_class_chr", "layer3_plot")
  ) %>%
  mutate(
    edge_type = "general_to_layer3",
    edge_value = n_features,
    edge_colour = class_colors[general_class_chr],
    line_width = scale_lw(n_features, max(n_features, na.rm = TRUE), min_lw = 0.50, max_lw = 6.2),
    alpha_value = general_to_layer3_alpha
  ) %>%
  select(edge_type, from_id, to_id, general_class = general_class_chr, edge_value, edge_colour, line_width, alpha_value, x0, y0, x1, y1)

# Layer 3 -> Layer 4
edges_3_4 <- edge_layer3_layer4 %>%
  mutate(general_class_chr = as.character(general_class)) %>%
  left_join(
    layer3_nodes %>%
      mutate(general_class_chr = as.character(general_class)) %>%
      select(general_class_chr, layer3_plot, from_id = node_id, x0 = anchor_out, y0 = y),
    by = c("general_class_chr", "layer3_plot")
  ) %>%
  left_join(
    layer4_nodes %>%
      mutate(general_class_chr = as.character(general_class)) %>%
      select(general_class_chr, layer3_plot, layer4_plot, to_id = node_id, x1 = anchor_in, y1 = y),
    by = c("general_class_chr", "layer3_plot", "layer4_plot")
  ) %>%
  mutate(
    edge_type = "layer3_to_layer4",
    edge_value = n_features,
    edge_colour = class_colors[general_class_chr],
    line_width = scale_lw(n_features, max(n_features, na.rm = TRUE), min_lw = 0.50, max_lw = 6.2),
    alpha_value = layer3_to_layer4_alpha
  ) %>%
  select(edge_type, from_id, to_id, general_class = general_class_chr, edge_value, edge_colour, line_width, alpha_value, x0, y0, x1, y1)

edge_list <- list(edges_1_2)
if (show_layer3) edge_list <- c(edge_list, list(edges_2_3))
if (show_layer4) edge_list <- c(edge_list, list(edges_3_4))

edges <- bind_rows(edge_list) %>%
  mutate(edge_id = row_number())

edges_curve <- edges %>%
  group_by(edge_id) %>%
  group_modify(~ make_bezier_points(.x$x0, .x$y0, .x$x1, .x$y1) %>%
                 mutate(
                   edge_type = .x$edge_type,
                   general_class = .x$general_class,
                   edge_value = .x$edge_value,
                   edge_colour = .x$edge_colour,
                   line_width = .x$line_width,
                   alpha_value = .x$alpha_value
                 )) %>%
  ungroup()


# ==================================================================================================
# 10. Plot, with optional 2/3/4-layer display
# ==================================================================================================

feature_tick_nodes <- feature_nodes %>%
  mutate(
    tick_xmin = x - 0.010,
    tick_xmax = x - 0.006,
    tick_ymin = y - 0.0010,
    tick_ymax = y + 0.0010,
    label_x = x - 0.014
  )

general_bar_nodes <- general_nodes %>%
  mutate(
    general_class_chr = as.character(general_class),
    fill_colour = class_colors[general_class_chr],
    bar_xmin = x - 0.040,
    bar_xmax = x - 0.030,
    label_x = x - 0.022
  )

layer3_bar_nodes <- layer3_nodes %>%
  mutate(
    general_class_chr = as.character(general_class),
    fill_colour = class_colors[general_class_chr],
    bar_xmin = x - 0.040,
    bar_xmax = x - 0.031,
    label_x = x - 0.024
  )

layer4_bar_nodes <- layer4_nodes %>%
  mutate(
    general_class_chr = as.character(general_class),
    fill_colour = class_colors[general_class_chr],
    bar_xmin = x - 0.050,
    bar_xmax = x - 0.042,
    label_x = x - 0.038
  )

subtitle_by_layer <- case_when(
  n_layers_to_show == 2 ~ "Matched feature / compound → general source class",
  n_layers_to_show == 3 ~ "Matched feature / compound → general source class → category-specific subclass",
  n_layers_to_show == 4 ~ "Matched feature / compound → general source class → category-specific subclass → most specific source detail"
)

p_graph1a <- ggplot() +
  geom_path(
    data = edges_curve,
    aes(x = x, y = y, group = edge_id, colour = edge_colour, linewidth = line_width, alpha = alpha_value),
    lineend = "butt"
  ) +
  scale_colour_identity() +
  scale_linewidth_identity() +
  scale_alpha_identity() +

  # Layer 1: feature ticks and labels
  geom_rect(
    data = feature_tick_nodes,
    aes(xmin = tick_xmin, xmax = tick_xmax, ymin = tick_ymin, ymax = tick_ymax),
    fill = "black",
    colour = "black",
    linewidth = 0
  ) +
  geom_text(
    data = feature_tick_nodes,
    aes(x = label_x, y = y, label = node_label),
    hjust = 1,
    vjust = 0.5,
    size = 1.35,
    colour = "#2F2F2F"
  ) +

  # Layer 2: general class bars
  geom_rect(
    data = general_bar_nodes,
    aes(xmin = bar_xmin, xmax = bar_xmax, ymin = bottom, ymax = top, fill = fill_colour),
    colour = "black",
    linewidth = 0.25
  ) +
  geom_text(
    data = general_bar_nodes,
    aes(x = label_x, y = y, label = node_label),
    hjust = 0,
    vjust = 0.5,
    size = 2.7,
    fontface = "bold",
    lineheight = 0.9,
    colour = "#1F1F1F"
  ) +
  scale_fill_identity() +

  annotate("text", x = x_feature + 0.02, y = 1.006,
           label = paste0("Layer 1: ", n_distinct(edge_feature_general$feature_id), " features"),
           fontface = "bold", size = 3.0, hjust = 0.5) +
  annotate("text", x = x_general, y = 1.006, label = "Layer 2: general class",
           fontface = "bold", size = 3.0, hjust = 0.5)

if (show_layer3) {
  p_graph1a <- p_graph1a +
    geom_rect(
      data = layer3_bar_nodes,
      aes(xmin = bar_xmin, xmax = bar_xmax, ymin = bottom, ymax = top, fill = fill_colour),
      colour = "black",
      linewidth = 0.22
    ) +
    geom_text(
      data = layer3_bar_nodes,
      aes(x = label_x, y = y, label = node_label),
      hjust = 0,
      vjust = 0.5,
      size = 2.05,
      lineheight = 0.85,
      colour = "#222222"
    ) +
    annotate("text", x = x_layer3, y = 1.006, label = "Layer 3: subclass",
             fontface = "bold", size = 3.0, hjust = 0.5)
}

if (show_layer4) {
  p_graph1a <- p_graph1a +
    geom_rect(
      data = layer4_bar_nodes,
      aes(xmin = bar_xmin, xmax = bar_xmax, ymin = bottom, ymax = top, fill = fill_colour),
      colour = "black",
      linewidth = 0.20
    ) +
    geom_text(
      data = layer4_bar_nodes,
      aes(x = label_x, y = y, label = node_label),
      hjust = 0,
      vjust = 0.5,
      size = 1.65,
      lineheight = 0.82,
      colour = "#222222"
    ) +
    annotate("text", x = x_layer4, y = 1.006, label = "Layer 4: most specific layer",
             fontface = "bold", size = 3.0, hjust = 0.5)
}

p_graph1a <- p_graph1a +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1.035), expand = FALSE, clip = "off") +
  theme_void(base_size = 12) +
  labs(
    title = paste0("Graph 1a. Origin classification of MAPS features (", n_layers_to_show, " layers shown)"),
    subtitle = subtitle_by_layer,
    caption = "Node bar length and edge width represent unique feature counts. Features with multiple MetOriginDB source annotations contribute to multiple flows, so counts across source classes are not mutually exclusive."
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 10.5, hjust = 0.5, colour = "#444444"),
    plot.caption = element_text(size = 8.5, hjust = 0, colour = "#444444"),
    plot.margin = margin(t = 24, r = 20, b = 20, l = 95)
  )

# Add layer number to file names so you can generate 2-, 3-, and 4-layer versions without overwriting.
graph1a_png_layers <- file.path(out_dir, paste0("Graph1a_", n_layers_to_show, "layer_sankey.png"))

safe_ggsave(graph1a_png_layers, p_graph1a, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")

cat("\nGraph 1a saved:\n")
cat(graph1a_png_layers, "\n")

cat("\nTables saved:\n")
cat(graph1a_long_csv, "\n")
cat(graph1a_feature_general_csv, "\n")
cat(graph1a_general_layer3_csv, "\n")
cat(graph1a_layer3_layer4_csv, "\n")
cat(graph1a_cleaned_input_csv, "\n")

####################################################################################################
# Notes on n_layers_to_show
####################################################################################################
# n_layers_to_show <- 2
#   Feature / compound -> General source class
#
# n_layers_to_show <- 3
#   Feature / compound -> General source class -> Category-specific subclass
#
# n_layers_to_show <- 4
#   Feature / compound -> General source class -> Category-specific subclass -> Most specific detail
#
# The full Graph1a_feature_general_subclass_detail_long_table.csv is still exported in all cases.
# The n_layers_to_show setting only controls what is drawn in the figure.
####################################################################################################
# ==================================================================================================
# 11. Circular 2-layer network with origin nodes in the middle + UpSet plot
# ==================================================================================================

required_packages_extra <- c("ggraph", "igraph", "tidygraph", "ComplexHeatmap", "circlize")

for (pkg in required_packages_extra) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(ggraph)
library(igraph)
library(tidygraph)
library(ComplexHeatmap)
library(circlize)
library(grid)

graph1a_circular_center_png <- file.path(out_dir, "Graph1a_2layer_circular_center_nodes.png")
graph1b_upset_png <- file.path(out_dir, "Graph1b_origin_upset.png")

# --------------------------------------------------------------------------------------------------
# A. Circular 2-layer network
# --------------------------------------------------------------------------------------------------

origin_order <- c("Food", "Plant", "Human", "Bacteria", "Animal", "Drug", "Environment")
origin_order <- origin_order[origin_order %in% unique(as.character(edge_feature_general$general_class))]

network_edges <- edge_feature_general %>%
  mutate(
    feature_node = paste0("Feature: ", feature_label),
    origin_node = paste0("Origin: ", as.character(general_class)),
    general_class_chr = as.character(general_class)
  ) %>%
  distinct(feature_node, origin_node, general_class_chr) %>%
  transmute(
    from = feature_node,
    to = origin_node,
    general_class = general_class_chr
  )

feature_nodes_network <- network_edges %>%
  distinct(name = from) %>%
  arrange(name) %>%
  mutate(
    node_type = "Feature",
    general_class = "metabolite",
    node_label = "",
    x = cos(seq(0, 2 * pi, length.out = n() + 1)[-1]),
    y = sin(seq(0, 2 * pi, length.out = n() + 1)[-1]),
    node_size = 0.55
  )

origin_angles <- seq(pi / 2, pi / 2 - 2 * pi, length.out = length(origin_order) + 1)[-length(origin_order) - 1]
names(origin_angles) <- origin_order

origin_nodes_network <- network_edges %>%
  distinct(general_class) %>%
  mutate(general_class = factor(general_class, levels = origin_order)) %>%
  arrange(general_class) %>%
  mutate(
    name = paste0("Origin: ", as.character(general_class)),
    node_type = "Origin",
    n_features = map_int(as.character(general_class), ~ n_distinct(network_edges$from[network_edges$general_class == .x])),
    node_label = as.character(general_class),
    angle = origin_angles[as.character(general_class)],
    x = 0.43 * cos(angle),
    y = 0.43 * sin(angle),
    node_size = scales::rescale(n_features, to = c(8, 18))
  ) %>%
  select(name, node_type, general_class, node_label, x, y, node_size, n_features)

nodes_network <- bind_rows(
  feature_nodes_network %>%
    select(name, node_type, general_class, node_label, x, y, node_size),
  origin_nodes_network %>%
    select(name, node_type, general_class, node_label, x, y, node_size)
) %>%
  mutate(
    general_class = as.character(general_class),
    label_x = if_else(node_type == "Origin", x * 1.08, x),
    label_y = if_else(node_type == "Origin", y * 1.08, y)
  )

node_lookup <- nodes_network %>%
  mutate(node_index = row_number()) %>%
  select(name, node_index)

edges_network_indexed <- network_edges %>%
  left_join(node_lookup, by = c("from" = "name")) %>%
  rename(from_idx = node_index) %>%
  left_join(node_lookup, by = c("to" = "name")) %>%
  rename(to_idx = node_index) %>%
  filter(!is.na(from_idx), !is.na(to_idx))

graph_network <- tbl_graph(
  nodes = nodes_network,
  edges = edges_network_indexed %>% transmute(from = from_idx, to = to_idx, general_class),
  directed = FALSE
)

p_graph1a_circular_center <- ggraph(graph_network, layout = "manual", x = x, y = y) +
  geom_edge_link(
    aes(colour = general_class),
    alpha = 0.16,
    linewidth = 0.18,
    show.legend = FALSE
  ) +
  geom_node_point(
    data = function(x) x %>% filter(node_type == "Feature"),
    aes(x = x, y = y),
    size = 0.35,
    colour = "#222222",
    alpha = 0.72
  ) +
  geom_node_point(
    data = function(x) x %>% filter(node_type == "Origin"),
    aes(x = x, y = y, fill = general_class, size = node_size),
    shape = 21,
    colour = "white",
    stroke = 0.9,
    show.legend = FALSE
  ) +
  geom_node_text(
    data = function(x) x %>% filter(node_type == "Origin"),
    aes(x = label_x, y = label_y, label = node_label),
    size = 5,
    fontface = "bold",
    colour = "#111111"
  ) +
  scale_edge_colour_manual(values = class_colors) +
  scale_fill_manual(values = class_colors) +
  scale_size_identity() +
  coord_equal(xlim = c(-1.12, 1.12), ylim = c(-1.12, 1.12), clip = "off") +
  theme_void(base_size = 12) +
  labs(
    title = "Graph 1a. Circular two-layer origin network",
    subtitle = paste0(
      n_distinct(edge_feature_general$feature_id),
      " MAPS features linked to ",
      n_distinct(edge_feature_general$general_class),
      " MetOriginDB source classes"
    ),
    caption = "Outer ring = MAPS features. Middle nodes = general MetOriginDB source classes. Features with multiple source annotations connect to multiple class nodes."
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 10.5, hjust = 0.5, colour = "#444444"),
    plot.caption = element_text(size = 8.5, hjust = 0, colour = "#444444"),
    plot.margin = margin(20, 20, 20, 20)
  )

ggsave(
  graph1a_circular_center_png,
  p_graph1a_circular_center,
  width = 12,
  height = 12,
  dpi = plot_dpi,
  bg = "white"
)

cat("\nCircular centre-node 2-layer PNG saved:\n")
cat(graph1a_circular_center_png, "\n")


# --------------------------------------------------------------------------------------------------
# B. Graph 1b UpSet plot of major origin-category overlaps
# --------------------------------------------------------------------------------------------------

upset_wide <- edge_feature_general %>%
  distinct(feature_id, general_class) %>%
  mutate(value = 1L) %>%
  pivot_wider(
    names_from = general_class,
    values_from = value,
    values_fill = 0
  )

upset_classes <- intersect(origin_order, names(upset_wide))

upset_mat <- upset_wide %>%
  select(all_of(upset_classes)) %>%
  as.data.frame()

upset_mat[] <- lapply(upset_mat, function(x) as.numeric(x) > 0)

comb_mat <- make_comb_mat(upset_mat)

png(
  graph1b_upset_png,
  width = 12,
  height = 7,
  units = "in",
  res = plot_dpi,
  bg = "white"
)

UpSet(
  comb_mat,
  top_annotation = upset_top_annotation(
    comb_mat,
    add_numbers = TRUE,
    gp = gpar(fill = "#5A5A5A")
  ),
  right_annotation = upset_right_annotation(
    comb_mat,
    add_numbers = TRUE,
    gp = gpar(fill = "#5A5A5A")
  ),
  row_names_gp = gpar(fontsize = 10),
  column_title = "Graph 1b. Overlap among MetOriginDB source classes",
  column_title_gp = gpar(fontsize = 15, fontface = "bold")
)

dev.off()

cat("\nGraph 1b UpSet PNG saved:\n")
cat(graph1b_upset_png, "\n")
