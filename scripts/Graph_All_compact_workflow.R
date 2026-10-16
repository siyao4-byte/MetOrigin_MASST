####################################################################################################
# Compact all-graphs workflow
#
# Runs the current working graph scripts in order:
#   1. Graph 1: MetOrigin origin Sankey/summary plots
#   2. Graph 2: domainMASST Sankey/UpSet plots
#   3. Graph 3: integrated MetOrigin + domainMASST evidence plots
#
####################################################################################################


# ==================================================================================================
# 1. Settings
# ==================================================================================================

run_graph1_metorigin <- TRUE
run_graph2_domainmasst <- identical(Sys.getenv("METORIGIN_INCLUDE_DOMAINMASST", "false"), "true")
run_graph3_integrated <- run_graph2_domainmasst

stop_if_a_graph_fails <- TRUE

graph_scripts <- c(
  Graph1_MetOrigin = "scripts/plot_recorded_associations.R",
  Graph2_domainMASST = "scripts/Graph2_domainMASST_sankey_and_upset.R",
  Graph3_integrated = "scripts/Graph3_integrated_origin_evidence_and_plots.R"
)


# ==================================================================================================
# 2. Helpers
# ==================================================================================================

script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0) return(NA_character_)
  normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
}

set_project_root <- function() {
  this_script <- script_path()
  if (!is.na(this_script)) {
    setwd(normalizePath(file.path(dirname(this_script), ".."), winslash = "/", mustWork = TRUE))
  }
  invisible(getwd())
}

run_graph_script <- function(label, path) {
  cat("\n", strrep("=", 96), "\n", sep = "")
  cat("Running ", label, "\n", sep = "")
  cat("Script: ", path, "\n", sep = "")
  cat(strrep("=", 96), "\n", sep = "")

  if (!file.exists(path)) {
    stop("Missing graph script: ", path)
  }

  started <- Sys.time()
  result <- tryCatch(
    {
      suppressPackageStartupMessages(
        source(path, local = new.env(parent = globalenv()), echo = FALSE)
      )
      TRUE
    },
    error = function(e) {
      message("\nFAILED: ", label)
      message(conditionMessage(e))
      FALSE
    }
  )

  elapsed <- round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1)

  if (isTRUE(result)) {
    cat("\nFinished ", label, " in ", elapsed, " seconds.\n", sep = "")
  } else if (isTRUE(stop_if_a_graph_fails)) {
    stop("Stopping because ", label, " failed. Check its inputs and dependencies.")
  }

  invisible(result)
}


# ==================================================================================================
# 3. Run
# ==================================================================================================

project_root <- set_project_root()

cat("Project root: ", project_root, "\n", sep = "")
cat("Graph 1 exports PNG and PDF from complete recorded-source evidence.\n")

run_flags <- c(
  Graph1_MetOrigin = run_graph1_metorigin,
  Graph2_domainMASST = run_graph2_domainmasst,
  Graph3_integrated = run_graph3_integrated
)

results <- logical(0)

for (label in names(graph_scripts)) {
  if (isTRUE(run_flags[[label]])) {
    results[[label]] <- run_graph_script(label, graph_scripts[[label]])
  } else {
    cat("\nSkipping ", label, "\n", sep = "")
    results[[label]] <- NA
  }
}

cat("\n", strrep("=", 96), "\n", sep = "")
cat("Compact graph workflow complete.\n")
print(results)
cat(strrep("=", 96), "\n", sep = "")
