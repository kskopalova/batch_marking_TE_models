# =========================================================
# R/graph_helpers.R
# Helper functions for make_graphs.R
# Assumes the working directory is the folder containing make_graphs.R
# (make_graphs.R sets this before sourcing this file).
# =========================================================

# ---- Results folder ----
results_dir <- "results"

if (!dir.exists(results_dir)) {
  stop("No 'results' folder found in ", getwd(),
       "\nThe working directory must be the folder containing make_graphs.R.")
}

# ---- List available results files ----
list_results <- function(dir = results_dir) {
  files <- list.files(dir, pattern = "\\.(rds|rdata|rda|csv)$", ignore.case = TRUE)
  if (length(files) == 0) message("No .rds, .RData, .rda or .csv files found in ", dir)
  files
}

# ---- Load an .RData / .rda file ----
# An .RData file can hold several objects. If it holds one, that object is
# returned; if it holds several, they are returned as a named list.
# Hidden objects that R saves with a workspace (e.g. .Random.seed) are ignored.
load_rdata <- function(file_path) {
  env <- new.env()
  obj_names <- load(file_path, envir = env)
  obj_names <- obj_names[!startsWith(obj_names, ".")]
  if (length(obj_names) == 0) stop("No objects found in ", file_path)
  if (length(obj_names) == 1) {
    env[[obj_names]]
  } else {
    message("File contains several objects: ", paste(obj_names, collapse = ", "),
            "\nReturning them as a list, e.g. results_file$", obj_names[1])
    mget(obj_names, envir = env)
  }
}

# ---- Load a results file by name ----
load_results <- function(file_name, dir = results_dir) {
  if (missing(file_name) || is.null(file_name) || !nzchar(file_name)) {
    stop("Please set `file_name` to a file in the '", dir, "' folder.")
  }
  
  file_path <- file.path(dir, file_name)
  
  if (!file.exists(file_path)) {
    stop("File not found: ", file_name,
         "\nAvailable files: ", paste(list_results(dir), collapse = ", "))
  }
  
  ext <- tolower(tools::file_ext(file_name))
  
  out <- switch(ext,
                rds   = readRDS(file_path),
                csv   = read.csv(file_path),
                rdata = ,
                rda   = load_rdata(file_path),
                stop("Unsupported file type: ", file_name,
                     " (use .rds, .RData, .rda or .csv)"))
  
  message("Loaded: ", file_path)
  out
}


# #########################################################
# PLOTTING: ID vs BM and NoTE vs TE boxplots
# #########################################################
# Works on the output of run.simulation() (emigration_sims/R/run_simulation.R),
# saved as results/<model>_ind.RData (individual identification) and
# results/<model>_bm.RData (batch-marking). Each file holds one object with
#   $settings              model, data.type, true values, ...
#   $emigration$results    one row per simulation: iteration, objective,
#                          convergence, <name>, <name>.Lo, <name>.Up
#   $emigration$summary    quantity, true.value, ..., coverage, ci.width
#   $noTE                  the same for the NoTE model (batch-marking, fit.noTE = TRUE)
#   $true.abundance        true numbers alive / available of every data set
# Requires ggplot2 and tidyr (loaded in make_graphs.R).

# ---- Load a run.simulation() result and check it ----
# file_name: e.g. "acbt_ind.RData" or "acbt_bm.RData" in the results folder
# data_type: "individual" or "batch" (what the file is expected to contain)
load_simulation <- function(file_name, data_type, dir = results_dir) {
  data_type <- match.arg(data_type, c("individual", "batch"))
  sim <- load_results(file_name, dir)
  
  # A saved workspace can hold several objects: pick the run.simulation() result(s)
  is_sim <- function(x) is.list(x) && !is.null(x$settings) && "emigration" %in% names(x)
  if (!is_sim(sim) && is.list(sim)) {
    found <- Filter(is_sim, sim)
    if (length(found) > 1) {
      found <- Filter(function(x) x$settings$data.type == data_type, found)
    }
    if (length(found) == 1) {
      message("Using object '", names(found), "' from ", file_name)
      sim <- found[[1]]
    } else if (length(found) > 1) {
      stop(file_name, " contains several run.simulation() results (",
           paste(names(found), collapse = ", "), "). Save each one in its own file.")
    }
  }
  if (!is_sim(sim)) {
    stop(file_name, " does not look like run.simulation() output (no $settings).")
  }
  if (is.null(sim$emigration)) {
    stop(file_name, " has no emigration-model results (was it run with fit.TE = FALSE?).")
  }
  if (sim$settings$data.type != data_type) {
    stop(file_name, " contains ", sim$settings$data.type, " data, not ", data_type, ".")
  }
  sim
}

# ---- Converged simulations and summary table of one fitted model ----
# part: "emigration" (TE model) or "noTE"
sim_tables <- function(sim, part = "emigration") {
  run <- sim[[part]]
  if (is.null(run) || is.null(run$results)) {
    stop("No ", part, " results in this file",
         if (part == "noTE") " (was it run with fit.noTE = TRUE and batch-marking data?)", ".")
  }
  res <- run$results
  if (is.null(run$summary) || sum(res$convergence == 0) == 0) {
    stop("No converged simulations for the ", part, " model in this file.")
  }
  list(full    = res[res$convergence == 0, , drop = FALSE],
       summary = run$summary,
       model   = sim$settings$model)
}

# ---- Which quantities to plot ----
# "N", "alphaprod", "alpha" (constant model) etc. match exactly;
# "phi", "p", "alpha", "beta", "gamma", "Navail", "Nalive" match phi1, phi2, ...
select_params <- function(param, available) {
  if (param %in% available) return(param)
  out <- grep(paste0("^", param, "[0-9]+$"), available, value = TRUE)
  if (length(out) == 0) {
    stop("'", param, "' not found for this model. Available: ",
         paste(available, collapse = ", "))
  }
  out[order(as.integer(sub("^[A-Za-z]+", "", out)))]
}

# ---- Axis labels: phi1 -> phi[1], alpha -> alpha, p12 -> p[12] ----
param_labels <- function(params) {
  labs <- parse(text = sub("^([A-Za-z]+)([0-9]+)$", "\\1[\\2]", params))
  names(labs) <- params
  labs
}

# ---- Legend label of an emigration model: atbt -> alpha[t] beta[t] ----
model_label <- function(model) {
  switch(model,
         cr   = quote(alpha[t]),
         atbt = quote(alpha[t] * beta[t]),
         atbc = quote(alpha[t] * beta[c]),
         acbt = quote(alpha[c] * beta[t]),
         acbc = quote(alpha[c] * beta[c]),
         model)
}

# ---- Automatic y-limits from the boxplot whiskers, means and true values ----
auto_ylim <- function(values, groups, extra = NULL, pad = 0.25) {
  ends <- unlist(lapply(split(values, groups), function(v) {
    v <- v[is.finite(v)]
    if (length(v) == 0) return(NULL)
    c(range(boxplot.stats(v)$stats), mean(v))
  }))
  ends <- c(ends, extra)   # true values, so the reference line stays in the box area
  lo <- min(ends, na.rm = TRUE); hi <- max(ends, na.rm = TRUE)
  span <- hi - lo
  if (span == 0) span <- max(abs(hi) * 0.1, 0.01)
  c(lo - pad * span, hi + pad * span)
}

# =========================================================
# SHARED PLOT: pairs of boxplots with coverage and CI width
# =========================================================
# a, b         : output of sim_tables() for the first and second group
# group_labels : legend labels (text or plotmath, e.g. list("NoTE", quote(alpha[t])))
# fill_values  : box fills of the two groups
# line_colour  : colour of the dashed true-value line
# y_label_fun, ci_label_fun : optional formatting of the y-axis and CI width labels
#                (NULL = defaults: 2 decimals for probabilities, plain for counts)
# n_style      : layout of the single-N plot. "blank" = no x labels, y title "N"
#                (ID vs BM); "tick" = one x tick labelled N, no y title (NoTE vs TE)
pair_boxplot <- function(param, a, b, group_labels, fill_values, line_colour,
                         n_style = c("blank", "tick"),
                         true_value = NULL, ylim = NULL, breaks = NULL, label_y = NULL,
                         ci_digits = 2, y_title = NULL,
                         y_label_fun = NULL, ci_label_fun = NULL) {
  
  n_style <- match.arg(n_style)
  params <- select_params(param, intersect(a$summary$quantity, b$summary$quantity))
  k <- length(params)
  lev <- c("g1", "g2")   # internal group codes; legend shows group_labels
  
  # ---- Data ----
  df <- rbind(a$full[, params, drop = FALSE], b$full[, params, drop = FALSE])
  df$mod <- factor(rep(lev, c(nrow(a$full), nrow(b$full))), levels = lev)
  dat <- tidyr::pivot_longer(df, cols = tidyselect::all_of(params), names_to = "param")
  dat$param <- factor(dat$param, levels = params)
  
  is_count <- all(grepl("^N", params))   # N, Navail, Nalive, N1..: no 2-decimal formatting
  
  # ---- Coverage / CI width, interleaved a_1, b_1, a_2, b_2, ... ----
  s_a <- a$summary[match(params, a$summary$quantity), ]
  s_b <- b$summary[match(params, b$summary$quantity), ]
  # coverage is a % of converged simulations, so it can be non-integer: 1 decimal
  coverage_vals <- as.character(round(as.vector(rbind(s_a$coverage, s_b$coverage)), 1))
  mean_ci_vals  <- as.vector(rbind(s_a$ci.width, s_b$ci.width))
  if (!is.null(ci_digits)) mean_ci_vals <- round(mean_ci_vals, ci_digits)   # NULL = not rounded
  ci_text <- if (!is.null(ci_label_fun)) ci_label_fun(mean_ci_vals)
  else if (is_count) mean_ci_vals else sprintf("%.2f", mean_ci_vals)
  
  # ---- True values (from the summary unless supplied) ----
  if (is.null(true_value)) true_value <- s_a$true.value
  if (length(true_value) > 1 && isTRUE(all.equal(min(true_value), max(true_value)))) {
    true_value <- true_value[1]   # all the same -> one dashed line
  }
  
  # ---- Layout ----
  single_N <- identical(params, "N")
  if (is.null(ylim)) {
    ylim <- auto_ylim(dat$value, interaction(dat$mod, dat$param), extra = true_value)
  }
  span <- diff(ylim)
  y_pos <- c(cov_title = ylim[2],
             cov       = ylim[2] - 0.10 * span,
             ci_title  = ylim[1] + 0.12 * span,
             ci        = ylim[1] + 0.03 * span)
  if (!is.null(label_y)) y_pos[names(label_y)] <- label_y
  
  if (single_N) {
    x_lab <- c(1, 2)
    x_mid <- 1.5
  } else {
    x_lab <- rep(seq_len(k), each = 2) + c(-0.185, 0.185)
    x_mid <- (k + 1) / 2
  }
  
  df_labels <- data.frame(label_cov  = coverage_vals,
                          label_mean = ci_text,
                          x          = x_lab,
                          y_cov      = y_pos[["cov"]],
                          y_mean     = y_pos[["ci"]])
  
  dodge <- position_dodge(width = 0.75)
  
  # ---- Plot ----
  if (single_N) {
    dat$mod_x <- as.numeric(dat$mod)
    p <- ggplot(dat, aes(x = mod_x, y = value, fill = mod, group = mod))
  } else {
    p <- ggplot(dat, aes(x = param, y = value, fill = mod))
  }
  
  p <- p +
    geom_boxplot(fatten = NULL, outlier.shape = NA, position = dodge, width = 0.75) +
    stat_summary(fun = mean,
                 geom = "errorbar",
                 aes(ymax = after_stat(y), ymin = after_stat(y)),
                 width = 0.75,
                 linetype = "solid",
                 colour = "black",
                 position = dodge) +
    scale_fill_manual(labels = group_labels,
                      values = fill_values,
                      name = NULL)
  
  if (single_N) {
    p <- p +
      scale_x_continuous(breaks = if (n_style == "tick") 1.5 else NULL,
                         labels = if (n_style == "tick") "N" else NULL,
                         limits = c(0.5, 2.5), expand = c(0, 0)) +
      scale_y_continuous(limits = ylim,
                         breaks = if (is.null(breaks)) waiver() else breaks,
                         labels = if (is.null(y_label_fun)) waiver() else y_label_fun,
                         oob = scales::oob_keep) +
      labs(x = NULL, y = if (!is.null(y_title)) y_title else if (n_style == "blank") "N" else NULL,
           title = NULL)
  } else {
    if (is.null(breaks)) {
      breaks <- pretty(ylim)
      breaks <- breaks[breaks >= ylim[1] & breaks <= ylim[2]]
    }
    p <- p +
      scale_x_discrete(labels = param_labels(params)) +
      scale_y_continuous(breaks = breaks,
                         labels = if (!is.null(y_label_fun)) y_label_fun
                         else if (is_count) waiver() else function(x) sprintf("%.2f", x),
                         limits = ylim,
                         oob = scales::oob_keep) +
      labs(y = y_title, x = NULL, title = NULL)
  }
  
  # ---- True value(s): one line, or one short line per parameter ----
  if (length(true_value) == 1) {
    p <- p + geom_hline(yintercept = true_value, linetype = "dashed",
                        colour = line_colour, linewidth = 1.1)
  } else if (length(true_value) == k) {
    p <- p + geom_segment(data = data.frame(x = seq_len(k) - 0.45,
                                            xend = seq_len(k) + 0.45,
                                            y = true_value),
                          aes(x = x, xend = xend, y = y, yend = y),
                          inherit.aes = FALSE, linetype = "dashed",
                          colour = line_colour, linewidth = 1.1)
  } else {
    stop("true_value must have length 1 or ", k, " (one per parameter).")
  }
  
  # ---- Coverage and CI width labels ----
  p <- p +
    annotate("text", label = "% Coverage:", size = 4.5,
             x = x_mid, y = y_pos[["cov_title"]]) +
    geom_text(data = df_labels,
              aes(x = x, y = y_cov, label = label_cov),
              inherit.aes = FALSE, size = 3.8) +
    annotate("text", label = "Mean CI width:", size = 4.5,
             x = x_mid, y = y_pos[["ci_title"]]) +
    geom_text(data = df_labels,
              aes(x = x, y = y_mean, label = label_mean),
              inherit.aes = FALSE, size = 3.8)
  
  # ---- Theme ----
  base <- theme(panel.border = element_blank(),
                axis.text.y = element_text(size = 12),
                legend.position = "bottom",
                legend.text = element_text(size = 12))
  if (single_N && n_style == "blank") {
    p + base + theme(axis.text.x = element_blank(),
                     axis.ticks = element_blank())
  } else if (single_N) {
    p + base + theme(axis.text.x = element_text(size = 12),
                     axis.ticks.x = element_line())
  } else {
    p + base + theme(axis.text.x = element_text(size = 12),
                     axis.title = element_text(size = 12))
  }
}

# =========================================================
# ID vs BM (same model, two runs)
# =========================================================
# param      : "N", "p", "phi", "alpha", "beta", "gamma", "alphaprod",
#              "betaprod", "Navail", "Nalive", or one exact name ("phi2")
# sim_id     : run.simulation() output for individual identification
# sim_bm     : run.simulation() output for batch-marking (same model)
# true_value : NULL = taken from the summary tables (true.value column);
#              or one value / one value per parameter to override
# ylim, breaks, label_y : optional manual overrides. Values outside ylim
#              are still used for the boxes and means (only the view is cropped).
#              label_y is a named vector with cov_title, cov, ci_title, ci (y positions).
# Only simulations that converged are plotted (as in the summary tables).
plot_id_vs_bm <- function(param, sim_id, sim_bm,
                          true_value = NULL,
                          ylim = NULL, breaks = NULL, label_y = NULL,
                          group_labels = c("ID", "BM")) {
  id <- sim_tables(sim_id)
  bm <- sim_tables(sim_bm)
  if (id$model != bm$model) {
    stop("The two runs are for different models: ", id$model, " (ID) and ", bm$model, " (BM).")
  }
  pair_boxplot(param, id, bm,
               group_labels = group_labels,
               fill_values  = c("white", "grey"),
               line_colour  = "red",
               n_style      = "blank",
               true_value = true_value, ylim = ylim, breaks = breaks, label_y = label_y)
}

# =========================================================
# NoTE vs TE (same simulated data sets, one batch-marking run)
# =========================================================
# param   : "N" (superpopulation size), "phi" (survival) or "p" (capture).
#           Only parameters that mean the same in both models are allowed
#           (NoTE's beta are entry probabilities, TE's beta return probabilities).
#           For phi, only the survival probabilities both models estimate are
#           plotted (e.g. phi1-phi4 for atbt and cr, which have no phi5).
# results : run.simulation() output with fit.TE = TRUE and fit.noTE = TRUE
# te_label: legend label of the TE model; default from the model, e.g. alpha[t]*beta[t]
plot_note_vs_te <- function(param, results,
                            true_value = NULL,
                            ylim = NULL, breaks = NULL, label_y = NULL,
                            te_label = NULL) {
  allowed <- c("N", "phi", "p")
  if (!sub("[0-9]+$", "", param) %in% allowed) {
    stop("NoTE vs TE plots are for N, phi or p (parameters with the same meaning in both models).")
  }
  te   <- sim_tables(results, "emigration")
  note <- sim_tables(results, "noTE")
  if (is.null(te_label)) te_label <- model_label(te$model)
  pair_boxplot(param, note, te,
               group_labels = list("NoTE", te_label),
               fill_values  = c("grey", "#117733"),
               line_colour  = "darkgrey",
               n_style      = "tick",
               true_value = true_value, ylim = ylim, breaks = breaks, label_y = label_y)
}

# =========================================================
# NoTE vs TE: period-level abundances N_k, relative bias
# =========================================================
# Relative bias, coverage and CI width of N_1..N_K for one fitted model,
# computed data set by data set against the true numbers of that data set
# (results$true.abundance), using only converged simulations, as in
# summarise.sims() of the simulation code.
#   part : "emigration" (TE) or "noTE"
#   type : "alive" or "available"
# Columns used:
#   NoTE                : Nalive<k> (NoTE has no emigration: alive = available)
#   TE, type = alive    : Navail1 for k = 1 (N_1 is the same), Nalive<k> for k >= 2
#   TE, type = available: Navail<k>
abundance_tables <- function(results, part, type) {
  run <- sim_tables(results, part)
  res <- run$full
  truth_mat <- if (type == "alive") results$true.abundance$alive else results$true.abundance$avail
  if (is.null(truth_mat)) stop("No true abundances ($true.abundance) in this file.")
  K <- nrow(truth_mat)
  
  col_name <- function(k) {
    if (part == "noTE") return(paste0("Nalive", k))
    if (type == "available" || k == 1) paste0("Navail", k) else paste0("Nalive", k)
  }
  
  periods <- Filter(function(k) col_name(k) %in% names(res), seq_len(K))
  full <- data.frame(row.names = seq_len(nrow(res)))
  summ <- data.frame(quantity = character(0), true.value = numeric(0),
                     coverage = numeric(0), ci.width = numeric(0))
  
  for (k in periods) {
    q   <- col_name(k)
    est <- res[[q]]; lo <- res[[paste0(q, ".Lo")]]; up <- res[[paste0(q, ".Up")]]
    tv  <- truth_mat[k, res$iteration]          # true N_k of each data set
    na.rows <- is.na(est) | is.na(lo) | is.na(up)
    n.ok    <- sum(!na.rows)
    nm <- paste0("N", k)
    full[[nm]] <- (est - tv) / tv * 100           # relative bias (%) per data set
    summ[nrow(summ) + 1, ] <- list(
      nm, 0,
      if (n.ok == 0) NA_real_ else
        sum(tv[!na.rows] >= lo[!na.rows] & tv[!na.rows] <= up[!na.rows]) / n.ok * 100,
      mean(up - lo, na.rm = TRUE))
  }
  list(full = full, summary = summ, model = run$model)
}

# type     : "alive" (numbers alive) or "available" (numbers available for capture).
#            For NoTE both use its Nalive<k>; N_1 is the same for both TE quantities.
#            Only the periods both models estimate are plotted (TE reports fewer
#            Nalive periods than NoTE).
# results  : run.simulation() output with fit.TE = TRUE and fit.noTE = TRUE
# Boxes show the relative bias (%) of every converged simulation, (estimate - true) /
# true x 100, with true = the true N_k of that data set; the dashed line is at 0.
# Coverage is against the same true values; CI width is in numbers of individuals,
# shown with two decimals (e.g. 431.60).
plot_note_vs_te_abundance <- function(type = c("alive", "available"), results,
                                      ylim = NULL, breaks = NULL, label_y = NULL,
                                      te_label = NULL,
                                      y_title = "Relative bias (%)") {
  type <- match.arg(type)
  note <- abundance_tables(results, "noTE", type)
  te   <- abundance_tables(results, "emigration", type)
  if (is.null(te_label)) te_label <- model_label(te$model)
  p <- pair_boxplot("N", note, te,
                    group_labels = list("NoTE", te_label),
                    fill_values  = c("grey", "#117733"),
                    line_colour  = "darkgrey",
                    n_style      = "tick",
                    true_value = 0, ylim = ylim, breaks = breaks, label_y = label_y,
                    ci_digits = NULL, y_title = y_title,
                    ci_label_fun = function(x) sprintf("%.2f", x))
  p
}

# =========================================================
# NoTE vs TE: parameters (N, phi, p), relative bias
# =========================================================
# Relative bias (%) of every converged simulation, (estimate - true) / true x 100,
# with the true value from the summary table (true.value). Coverage and CI width
# are taken from the summary tables (CI width on the parameter's own scale).
relbias_tables <- function(results, part, param) {
  run <- sim_tables(results, part)
  out <- data.frame(row.names = seq_len(nrow(run$full)))
  summ <- run$summary
  for (q in param) {
    tv <- summ$true.value[summ$quantity == q]
    out[[q]] <- (run$full[[q]] - tv) / tv * 100
  }
  summ <- summ[summ$quantity %in% param, c("quantity", "coverage", "ci.width")]
  summ$true.value <- 0
  list(full = out, summary = summ, model = run$model)
}

# param   : "phi" (survival), "N" (superpopulation size) or "p" (capture)
# results : run.simulation() output with fit.TE = TRUE and fit.noTE = TRUE
# Only the parameters both models estimate are plotted (e.g. phi1-phi4 for atbt and cr).
plot_note_vs_te_relbias <- function(param, results,
                                    ylim = NULL, breaks = NULL, label_y = NULL,
                                    te_label = NULL, y_title = "Relative bias (%)") {
  if (!sub("[0-9]+$", "", param) %in% c("N", "phi", "p")) {
    stop("NoTE vs TE plots are for N, phi or p (parameters with the same meaning in both models).")
  }
  note_all <- sim_tables(results, "noTE")
  te_all   <- sim_tables(results, "emigration")
  params <- select_params(param, intersect(note_all$summary$quantity, te_all$summary$quantity))
  note <- relbias_tables(results, "noTE", params)
  te   <- relbias_tables(results, "emigration", params)
  if (is.null(te_label)) te_label <- model_label(te$model)
  pair_boxplot(param, note, te,
               group_labels = list("NoTE", te_label),
               fill_values  = c("grey", "#117733"),
               line_colour  = "darkgrey",
               n_style      = "tick",
               true_value = 0, ylim = ylim, breaks = breaks, label_y = label_y,
               y_title = y_title, y_label_fun = waiver(),
               ci_label_fun = function(x) sprintf("%.2f", x))
}

# ---- Save a plot as png, pdf and eps ----
# Saved in results/plots/ by default (created if it does not exist).
save_plot <- function(plot, name, width, height, dpi = 400,
                      formats = c("png", "pdf", "eps"),
                      dir = file.path(results_dir, "plots")) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  for (f in formats) {
    ggsave(filename = file.path(dir, paste0(name, ".", f)),
           plot = plot, device = f,
           width = width, height = height, dpi = dpi)
  }
  message("Saved: ", file.path(dir, name), ".{", paste(formats, collapse = ","), "}")
  invisible(file.path(dir, paste0(name, ".", formats)))
}