# score_program.R — NPxIBD program scoring (locked manuscript formula, no new statistics)
#
# Two methods, both taken verbatim from the manuscript's own code:
#   "zmean"  : per-gene z-score across samples, score(set) = mean z(set genes present).
#              Fig.10 protein scores = NPxUC_UP zmean - NPxUC_DOWN zmean (same for CD),
#              exactly as in analysis/08_proteomics/01_pxd050335_cd_mucosa.R.
#   "ssgsea" : GSVA::gsva(ssgseaParam(...)) arm scores, the manuscript bulk method
#              (analysis/01_validation/01_validation_partB.R, GSVA v2 param-object API).
#
# Input: expr = genes x samples matrix or data.frame. Row names (or a first column named
# SYMBOL/symbol/gene) must be gene symbols. Input should be log-scale (log2 counts for
# transcriptomes, log2 LFQ for proteomes); the function does NOT transform scale.
# Gene sets: a named list of character vectors, e.g. list(NPxUC_UP=..., NPxUC_DOWN=...,
# NPxCD_UP=..., NPxCD_DOWN=...). Program sets ship with the repository as CSVs.
#
# Output: data.frame, rows = gene sets, columns = samples (arm-level scores).
# A combined program score for a UP/DOWN pair = score(UP) - score(DOWN); see
# combine_program_scores(). Coverage per set (hit genes / set size) is reported; sets
# with <50% coverage trigger an explicit warning, sets with 0 hits stop with an error.
#
# Duplicate symbols in expr are collapsed by keeping the row with maximum variance
# across samples (documented choice; manuscript matrices used here already have unique
# symbols, so this path does not trigger on manuscript data).
#
# Demo: running this file directly with SCORE_PROGRAM_DEMO=1 executes a small
# synthetic example (not manuscript data) at the end of the file.

score_program <- function(expr, sets, method = c("zmean", "ssgsea")) {
  method <- match.arg(method)
  if (is.data.frame(expr)) {
    cn <- colnames(expr)
    if (!is.null(cn) && tolower(cn[1]) %in% c("symbol", "gene", "symbols")) {
      rn <- expr[[1]]
      expr <- expr[, -1, drop = FALSE]
      rownames(expr) <- rn
    }
    expr <- as.matrix(expr)
  }
  if (is.null(rownames(expr))) stop("score_program: expr must have gene symbols as row names")
  storage.mode(expr) <- "numeric"
  # collapse duplicate symbols: keep max-variance row
  if (anyDuplicated(rownames(expr))) {
    v <- apply(expr, 1, function(x) stats::var(x, na.rm = TRUE))
    v[is.na(v)] <- -Inf
    ord <- order(rownames(expr), -v)
    expr <- expr[ord, , drop = FALSE]
    expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]
    message("score_program: duplicate symbols collapsed to max-variance row")
  }
  # coverage report
  for (nm in names(sets)) {
    hit <- sum(sets[[nm]] %in% rownames(expr))
    message(sprintf("score_program: set %s coverage %d/%d (%.1f%%)", nm, hit, length(sets[[nm]]), 100 * hit / length(sets[[nm]])))
    if (hit == 0) stop(sprintf("score_program: set %s has 0 genes in expr", nm))
    if (hit / length(sets[[nm]]) < 0.5) warning(sprintf("score_program: set %s coverage below 50%%", nm))
  }
  if (method == "zmean") {
    z <- t(scale(t(expr)))
    out <- do.call(rbind, lapply(sets, function(g) {
      g <- intersect(g, rownames(z))
      colMeans(z[g, , drop = FALSE], na.rm = TRUE)
    }))
    rownames(out) <- names(sets)
    return(as.data.frame(out))
  }
  # ssgsea (GSVA v2 param-object API)
  if (!requireNamespace("GSVA", quietly = TRUE)) stop("score_program: GSVA package required for method='ssgsea'")
  sc <- GSVA::gsva(GSVA::ssgseaParam(expr, sets))
  as.data.frame(sc)
}

# Combined program score for one UP/DOWN arm pair (Fig.10 definition).
# Returns a named numeric vector (names = samples). Note: arm_scores rows subset to a
# data.frame, so as.numeric() coercion is required (a bare row subtraction stays a data.frame).
combine_program_scores <- function(arm_scores, up_name, down_name) {
  v <- as.numeric(arm_scores[up_name, ]) - as.numeric(arm_scores[down_name, ])
  names(v) <- colnames(arm_scores)
  v
}

# Manuscript program sets (read from the shared discovery files; direction = dir_NP).
load_manuscript_sets <- function(pilot_dir = "/home/hatch/workspace/sci-np-ibd/pilot") {
  uc <- utils::read.csv(file.path(pilot_dir, "shared_NPxUC.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  cd <- utils::read.csv(file.path(pilot_dir, "shared_NPxCD.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  list(
    NPxUC_UP = uc$SYMBOL[uc$dir_NP == 1],
    NPxUC_DOWN = uc$SYMBOL[uc$dir_NP == -1],
    NPxCD_UP = cd$SYMBOL[cd$dir_NP == 1],
    NPxCD_DOWN = cd$SYMBOL[cd$dir_NP == -1]
  )
}

# ---- Demo (synthetic example; runs ONLY when SCORE_PROGRAM_DEMO=1) ----
# Not manuscript data: a tiny synthetic matrix (20 genes x 4 samples) showing
# the zmean call shape and combine_program_scores(). Run it with:
#   SCORE_PROGRAM_DEMO=1 Rscript score_program.R
# Sourcing this file (manuscript code, self-tests) never triggers the demo.
if (identical(Sys.getenv("SCORE_PROGRAM_DEMO"), "1")) {
  set.seed(1)
  demo_genes <- paste0("G", sprintf("%02d", 1:20))
  demo_expr <- matrix(stats::rnorm(20 * 4, mean = 8, sd = 1), nrow = 20,
                      dimnames = list(demo_genes, paste0("S", 1:4)))
  demo_sets <- list(DEMO_UP = demo_genes[1:10], DEMO_DOWN = demo_genes[11:20])
  demo_sc <- score_program(demo_expr, demo_sets, method = "zmean")
  print(demo_sc)
  print(combine_program_scores(demo_sc, "DEMO_UP", "DEMO_DOWN"))
  cat("score_program demo OK (synthetic data, not manuscript data)\n")
}
