#===============================================================================
# Plots from simulation results
#===============================================================================

# Packages
library(tidyverse)
library(patchwork)

# Set working directory to this script's folder
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  script_path <- rstudioapi::getSourceEditorContext()$path
  if (nzchar(script_path)) setwd(dirname(script_path))
}

# Helpers
helpers_path <- file.path("R", "graph_helpers.R")
if (!file.exists(helpers_path)) {
  stop("Can't find ", helpers_path, " in ", getwd(),
       "\nSave this script and run it from its own folder.")
}
source(helpers_path)

#----NoTE VS TE-----------------------------------------------------------------

# load data (run.simulation() output)
list_results()  
file_name <- "atbt_bm.RData"                 # change this
results <- load_simulation(file_name, "batch")
model_te <- results$settings$model           # used in the saved file names

p_note_N      <- plot_note_vs_te("N",   results)          # superpopulation size
p_note_phi    <- plot_note_vs_te("phi", results)          # survival (absolute bias)
p_note_phi_rb <- plot_note_vs_te_relbias("phi", results)  # survival (percentage bias)

# save plot (change if needed)
save_plot(p_note_N,
          paste0(model_te, "_note_vs_te_n"),
          width = 6,
          height = 4)


# period-level abundance N_k (percentage bias)
# abundance_type: "alive" (numbers alive) or "available" (available for capture).
# (Note that numbers alive and available for capture are the same for NoTE. 
# Only periods both models estimate are shown.)
abundance_type <- "available"                   # change this
p_note_Nk <- plot_note_vs_te_abundance(abundance_type, results)
p_note_Nk

# save plot (change if needed)
save_plot(p_note_Nk,
          paste0(model_te, "_note_vs_te_nk_", abundance_type),
          width = 10,
          height = 4)

#----INDIVIDUAL IDENTIFICATION VS BATCH-MARKING---------------------------------

# load data (run.simulation() output for the same model)
list_results()                              # shows available files
file_name_id <- "acbt_ind.RData"            # individual identification: change this
file_name_bm <- "acbt_bm.RData"             # batch-marking: change this
results_id <- load_simulation(file_name_id, "individual")
results_bm <- load_simulation(file_name_bm, "batch")
model_idbm <- results_id$settings$model     # used in the saved file names

# make plots
p_idbm_N     <- plot_id_vs_bm("N",     results_id, results_bm)
p_idbm_phi   <- plot_id_vs_bm("phi",   results_id, results_bm)
p_idbm_alpha <- plot_id_vs_bm("alpha", results_id, results_bm)
if (model_idbm != "cr") {                                         # cr has no beta
  p_idbm_beta <- plot_id_vs_bm("beta", results_id, results_bm)
}

# save plot (change if needed)
save_plot(p_idbm_N,
          paste0(model_idbm, "_id_vs_bm_n"),
          width = 6,
          height = 4)

