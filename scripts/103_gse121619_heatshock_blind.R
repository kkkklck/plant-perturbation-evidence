suppressPackageStartupMessages({
  library(Biobase)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1 || !args[1] %in% c("freeze", "target")) {
  stop("Usage: Rscript scripts/103_gse121619_heatshock_blind.R freeze|target")
}
mode <- args[1]

# Minimal compatibility registration for a legacy Monocle 2 CellDataSet.
# We only need the inherited ExpressionSet assayData/exprs slot.
if (!isClass("CellDataSet")) {
  setClass(
    "CellDataSet",
    contains = "ExpressionSet",
    slots = c(
      reducedDimS = "ANY",
      reducedDimW = "ANY",
      reducedDimA = "ANY",
      reducedDimK = "ANY",
      minSpanningTree = "ANY",
      cellPairwiseDistances = "ANY",
      expressionFamily = "ANY",
      lowerDetectionLimit = "ANY",
      dispFitInfo = "ANY",
      dim_reduce_type = "ANY",
      auxOrderingData = "ANY",
      auxClusteringData = "ANY"
    )
  )
}

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
base <- file.path(root, "data", "external_heatshock", "GSE121619")
inp <- file.path(base, "GSE121619_Control_Heatshock_cds.rds")
meta_path <- file.path(base, "GSE121619_Control_Heatshock_pData.tsv.gz")
out <- file.path(base, "plantvc_v1_heatshock_blind")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(out, paste0(mode, ".log"))
zz <- file(log_path, open = "wt")
sink(zz, type = "output")
sink(zz, type = "message")
on.exit({
  while (sink.number(type = "message") > 2) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  try(close(zz), silent = TRUE)
}, add = TRUE)

cat("GSE121619_HEATSHOCK_", toupper(mode), "_START\n", sep = "")
cat("R=", R.version.string, "\n", sep = "")
cat("input_bytes=", file.info(inp)$size, "\n", sep = "")
flush.console()

meta <- read.delim(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(meta) == 2085)
stopifnot(all(c("cell", "total_umis", "treatment_id", "louvain_component") %in% names(meta)))
meta$louvain_component <- as.integer(meta$louvain_component)
meta$treatment_id <- as.character(meta$treatment_id)

TARGETS <- c(2L, 3L, 10L)
DONORS <- c(5L, 6L, 8L)

cat("metadata_cells=", nrow(meta), "\n", sep = "")
cat("targets=", paste(TARGETS, collapse = ","), "\n", sep = "")
cat("donors=", paste(DONORS, collapse = ","), "\n", sep = "")
print(with(meta, table(louvain_component, treatment_id)))
flush.console()

obj <- readRDS(inp)
cat("READ_RDS_OK class=", paste(class(obj), collapse = ","), "\n", sep = "")

# The GEO object is a 2018 legacy Monocle object. Depending on R/package version,
# readRDS() may restore it either as a CellDataSet/ExpressionSet or as a nested list.
# Find the actual expression matrix by CELL NAMES, never by positional guessing.
is_matrix_like <- function(x) {
  is.matrix(x) || inherits(x, "Matrix")
}

matrix_matches_cells <- function(x, cells) {
  if (!is_matrix_like(x) || length(dim(x)) != 2L) return(FALSE)
  cnx <- colnames(x)
  rnx <- rownames(x)
  col_ok <- !is.null(cnx) && length(cnx) == length(cells) && setequal(cnx, cells)
  row_ok <- !is.null(rnx) && length(rnx) == length(cells) && setequal(rnx, cells)
  col_ok || row_ok
}

extract_expr_matrix <- function(x, cells, path = "root", depth = 0L, max_depth = 8L) {
  if (matrix_matches_cells(x, cells)) {
    if (!is.null(colnames(x)) && setequal(colnames(x), cells)) {
      return(list(E = x, path = path, transposed = FALSE))
    }
    if (!is.null(rownames(x)) && setequal(rownames(x), cells)) {
      return(list(E = t(x), path = path, transposed = TRUE))
    }
  }

  if (depth >= max_depth) return(NULL)

  # Search list elements. Prefer names that commonly hold ExpressionSet assays.
  if (is.list(x)) {
    nm <- names(x)
    if (is.null(nm)) nm <- rep("", length(x))
    priority <- c("exprs", "assayData", "assay_data", "counts", "data")
    ord <- order(!(nm %in% priority), seq_along(x))
    for (ii in ord) {
      child_name <- if (nzchar(nm[ii])) nm[ii] else paste0("[[", ii, "]]" )
      ans <- try(extract_expr_matrix(
        x[[ii]], cells,
        path = paste0(path, "$", child_name),
        depth = depth + 1L,
        max_depth = max_depth
      ), silent = TRUE)
      if (!inherits(ans, "try-error") && !is.null(ans)) return(ans)
    }
  }

  # Some serialized ExpressionSet assayData containers are environments.
  if (is.environment(x)) {
    keys <- ls(x, all.names = TRUE)
    priority <- c("exprs", "counts", "data")
    keys <- c(intersect(priority, keys), setdiff(keys, priority))
    for (key in keys) {
      child <- try(get(key, envir = x, inherits = FALSE), silent = TRUE)
      if (inherits(child, "try-error")) next
      ans <- try(extract_expr_matrix(
        child, cells,
        path = paste0(path, "$", key),
        depth = depth + 1L,
        max_depth = max_depth
      ), silent = TRUE)
      if (!inherits(ans, "try-error") && !is.null(ans)) return(ans)
    }
  }

  # Legacy Monocle S4 objects remember package="monocle" in their class
  # metadata. Calling inherits()/slotNames()/exprs() on them can therefore
  # trigger a namespace load even though we only need the raw assay slot.
  # Traverse raw attributes instead; S4 slots are stored there as ordinary
  # named attributes in the serialized object. We still accept a matrix only
  # when its cell names exactly match the external 2085-cell metadata.
  at <- try(attributes(x), silent = TRUE)
  if (!inherits(at, "try-error") && !is.null(at)) {
    an <- names(at)
    if (!is.null(an)) {
      priority <- c("assayData", "exprs", "counts", "data")
      keep <- setdiff(an, c("class", "names", "dim", "dimnames", "row.names"))
      keep <- c(intersect(priority, keep), setdiff(keep, priority))
      for (key in keep) {
        child <- at[[key]]
        ans <- try(extract_expr_matrix(
          child, cells,
          path = paste0(path, "@", key),
          depth = depth + 1L,
          max_depth = max_depth
        ), silent = TRUE)
        if (!inherits(ans, "try-error") && !is.null(ans)) return(ans)
      }
    }
  }

  NULL
}

found <- extract_expr_matrix(obj, meta$cell)
if (is.null(found)) {
  if (is.list(obj)) {
    cat("top_level_names=", paste(names(obj), collapse = ","), "\n", sep = "")
    cat("top_level_classes=", paste(vapply(obj, function(z) paste(class(z), collapse = "/"), character(1)), collapse = ","), "\n", sep = "")
  }
  stop("Could not find an expression matrix whose cell names exactly match all 2085 metadata cells")
}

E <- found$E
cat("EXPR_MATRIX_FOUND path=", found$path,
    " transposed=", found$transposed,
    " genes=", nrow(E),
    " cells=", ncol(E), "\n", sep = "")

cn <- colnames(E)
if (is.null(cn)) stop("Recovered expression matrix has no cell names")
if (length(cn) != nrow(meta)) stop("Recovered expression matrix cell count differs from metadata")
if (!setequal(cn, meta$cell)) {
  missing_rds <- setdiff(meta$cell, cn)
  missing_meta <- setdiff(cn, meta$cell)
  stop("Cell-name mismatch. metadata_not_rds=", length(missing_rds), "; rds_not_metadata=", length(missing_meta))
}
meta <- meta[match(cn, meta$cell), , drop = FALSE]
stopifnot(identical(meta$cell, cn))
cat("CELL_MAPPING_OK\n")

if (!(is.numeric(E) || inherits(E, "Matrix"))) stop("Recovered exprs is not numeric/matrix-like")

gene_ids <- rownames(E)
if (is.null(gene_ids)) gene_ids <- as.character(seq_len(nrow(E)))

sum_cols <- function(M, idx) {
  if (length(idx) == 0L) stop("Empty cell subset")
  as.numeric(rowSums(M[, idx, drop = FALSE]))
}

write_profile <- function(df, filename) {
  con <- gzfile(file.path(out, filename), open = "wt")
  write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE)
  close(con)
}

if (mode == "freeze") {
  allowed <- meta$treatment_id == "Control" |
    (meta$treatment_id == "HeatShock" & meta$louvain_component %in% DONORS)
  forbidden <- meta$treatment_id == "HeatShock" & meta$louvain_component %in% TARGETS

  cat("freeze_allowed_cells=", sum(allowed), "\n", sep = "")
  cat("freeze_forbidden_target_HS_cells=", sum(forbidden), "\n", sep = "")
  cat("freeze_target_HS_cells_used=", sum(allowed & forbidden), "\n", sep = "")
  if (sum(allowed & forbidden) != 0L) stop("Target HeatShock leakage detected")

  # Check raw-count semantics using allowed cells only.
  check_idx <- which(allowed)
  lib_from_expr <- colSums(E[, check_idx, drop = FALSE])
  lib_from_meta <- meta$total_umis[check_idx]
  relerr <- abs(lib_from_expr - lib_from_meta) / pmax(1, abs(lib_from_meta))
  cat("allowed_library_size_median_relative_error=", median(relerr), "\n", sep = "")
  cat("allowed_library_size_max_relative_error=", max(relerr), "\n", sep = "")
  if (median(relerr) > 1e-6) stop("exprs does not look like raw UMI counts under metadata total_umis")

  # Technical universe: all Control cells + HeatShock donor cells, exactly as frozen.
  total_allowed <- sum_cols(E, which(allowed))
  technical <- total_allowed >= 20

  response_cols <- list()
  donor_meta <- list()
  for (d in DONORS) {
    ci <- which(meta$treatment_id == "Control" & meta$louvain_component == d)
    hi <- which(meta$treatment_id == "HeatShock" & meta$louvain_component == d)
    cc <- sum_cols(E, ci)
    hc <- sum_cols(E, hi)
    c_lib <- sum(cc)
    h_lib <- sum(hc)
    r <- log1p(hc / h_lib * 1e6) - log1p(cc / c_lib * 1e6)
    response_cols[[paste0("donor_", d)]] <- r
    donor_meta[[as.character(d)]] <- c(control_cells = length(ci), heatshock_cells = length(hi), control_umi = c_lib, heatshock_umi = h_lib)
  }

  donor_mat <- do.call(cbind, response_cols)
  shared <- rowMeans(donor_mat)

  # Pooled-donor baseline/reference.
  donor_control_idx <- which(meta$treatment_id == "Control" & meta$louvain_component %in% DONORS)
  donor_hs_idx <- which(meta$treatment_id == "HeatShock" & meta$louvain_component %in% DONORS)
  pc <- sum_cols(E, donor_control_idx)
  ph <- sum_cols(E, donor_hs_idx)
  pooled <- log1p(ph / sum(ph) * 1e6) - log1p(pc / sum(pc) * 1e6)

  pred <- data.frame(
    gene = gene_ids,
    technical = technical,
    shared_prediction = shared,
    pooled_donor_prediction = pooled,
    donor_5 = donor_mat[, "donor_5"],
    donor_6 = donor_mat[, "donor_6"],
    donor_8 = donor_mat[, "donor_8"],
    stringsAsFactors = FALSE
  )
  pred$shared_abs_rank <- rank(-abs(pred$shared_prediction), ties.method = "min")
  pred$pooled_abs_rank <- rank(-abs(pred$pooled_donor_prediction), ties.method = "min")
  pred <- pred[order(pred$shared_abs_rank, pred$gene), ]
  write_profile(pred, "frozen_prediction.tsv.gz")
  writeLines(pred$gene[pred$technical], file.path(out, "technical_genes.txt"), useBytes = TRUE)

  used_cells <- data.frame(
    cell = meta$cell[allowed],
    treatment = meta$treatment_id[allowed],
    cluster = meta$louvain_component[allowed],
    stringsAsFactors = FALSE
  )
  write.table(used_cells, file.path(out, "freeze_cells_used.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  forbidden_cells <- data.frame(
    cell = meta$cell[forbidden],
    treatment = meta$treatment_id[forbidden],
    cluster = meta$louvain_component[forbidden],
    stringsAsFactors = FALSE
  )
  write.table(forbidden_cells, file.path(out, "target_heatshock_cells_withheld.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  cat("genes_all=", nrow(pred), "\n", sep = "")
  cat("genes_technical=", sum(technical), "\n", sep = "")
  for (d in names(donor_meta)) {
    x <- donor_meta[[d]]
    cat("donor_", d, " control_cells=", x[["control_cells"]], " heatshock_cells=", x[["heatshock_cells"]], "\n", sep = "")
  }
  cat("TARGET_HEATSHOCK_VALUES_NOT_USED_FOR_ANY_PREDICTION\n")
  cat("GSE121619_HEATSHOCK_FREEZE_SUCCESS\n")

} else {
  tech_file <- file.path(out, "technical_genes.txt")
  pred_file <- file.path(out, "frozen_prediction.tsv.gz")
  if (!file.exists(tech_file) || !file.exists(pred_file)) stop("Freeze artifacts missing")
  technical_genes <- readLines(tech_file, warn = FALSE)

  target_rows <- list()
  for (d in TARGETS) {
    ci <- which(meta$treatment_id == "Control" & meta$louvain_component == d)
    hi <- which(meta$treatment_id == "HeatShock" & meta$louvain_component == d)
    cc <- sum_cols(E, ci)
    hc <- sum_cols(E, hi)
    r <- log1p(hc / sum(hc) * 1e6) - log1p(cc / sum(cc) * 1e6)
    target_rows[[as.character(d)]] <- data.frame(
      gene = gene_ids,
      target_cluster = d,
      target_response = r,
      target_abs_rank = rank(-abs(r), ties.method = "min"),
      control_cells = length(ci),
      heatshock_cells = length(hi),
      stringsAsFactors = FALSE
    )
  }
  target <- do.call(rbind, target_rows)
  target$technical <- target$gene %in% technical_genes
  write_profile(target, "target_profiles.tsv.gz")
  cat("targets_scored=", paste(TARGETS, collapse = ","), "\n", sep = "")
  cat("target_rows=", nrow(target), "\n", sep = "")
  cat("GSE121619_HEATSHOCK_TARGET_EXPORT_SUCCESS\n")
}
