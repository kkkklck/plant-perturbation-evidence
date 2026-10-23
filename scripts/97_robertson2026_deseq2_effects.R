suppressPackageStartupMessages({
  library(DESeq2)
  library(Matrix)
})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
base <- file.path(root, "data", "external_rice", "Robertson2026")
pbdir <- file.path(base, "pseudobulk")
outdir <- file.path(base, "deseq2_effects_v0")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(outdir, "97_deseq2.log")
zz <- file(logf, open = "wt")
sink(zz, type = "output")

cat("ROBERTSON2026_DESEQ2_START\n")
cat("R=", R.version.string, "\n", sep = "")
cat("DESeq2=", as.character(packageVersion("DESeq2")), "\n", sep = "")
cat("Matrix=", as.character(packageVersion("Matrix")), "\n", sep = "")

counts <- readRDS(file.path(pbdir, "pseudobulk_CellAnnotation_raw_counts.rds"))
meta <- read.delim(file.path(pbdir, "pseudobulk_CellAnnotation_meta.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
design_tab <- read.delim(file.path(base, "sample_design.tsv"), check.names = FALSE, stringsAsFactors = FALSE)

if (!inherits(counts, "Matrix")) counts <- as(counts, "dgCMatrix")
if (!all(meta$pseudobulk_id %in% colnames(counts))) stop("pseudobulk metadata IDs missing from raw-count matrix")
if (nrow(counts) != 25985) stop("unexpected gene count: ", nrow(counts))

condition_map <- c(
  WellWatered = "Well.Watered",
  MildDrought = "Mild.Drought",
  ModerateDrought = "Moderate.Drought",
  MildSalinity = "Mild.Salinity",
  ModerateSalinity = "Moderate.Salinity"
)
design_tab$condition_obj <- unname(condition_map[design_tab$condition])
if (anyNA(design_tab$condition_obj)) stop("condition mapping failed")

control <- design_tab[design_tab$condition_obj == "Well.Watered", , drop = FALSE]
control_by_batch <- setNames(control$identity, as.character(control$batch))
stress_conditions <- c("Mild.Drought", "Moderate.Drought", "Mild.Salinity", "Moderate.Salinity")
states <- sort(unique(meta$state))
if (length(states) != 14) stop("expected 14 CellAnnotation states, got ", length(states))

pb_id <- function(state, sample) {
  z <- meta[meta$state == state & meta$sample == sample, , drop = FALSE]
  if (nrow(z) != 1) stop(state, "/", sample, ": expected one pseudobulk, got ", nrow(z))
  z$pseudobulk_id[[1]]
}

all_out <- vector("list", length(states) * length(stress_conditions))
model_audit <- vector("list", length(all_out))
k <- 0L

for (state in states) {
  cat("\nSTATE ", state, "\n", sep = "")
  for (cond in stress_conditions) {
    k <- k + 1L
    sd <- design_tab[design_tab$condition_obj == cond, , drop = FALSE]
    sd <- sd[order(sd$replicate), , drop = FALSE]
    if (nrow(sd) != 2) stop(cond, ": expected two stress replicates")

    s1 <- sd$identity[[1]]; s2 <- sd$identity[[2]]
    b1 <- as.character(sd$batch[[1]]); b2 <- as.character(sd$batch[[2]])
    c1 <- unname(control_by_batch[b1]); c2 <- unname(control_by_batch[b2])
    if (is.na(c1) || is.na(c2)) stop(cond, ": missing matched control batch")

    samples <- c(c1, s1, c2, s2)
    ids <- vapply(samples, function(s) pb_id(state, s), character(1))
    X <- as.matrix(counts[, ids, drop = FALSE])
    storage.mode(X) <- "integer"
    colnames(X) <- samples

    cd <- data.frame(
      row.names = samples,
      batch = factor(c(b1, b1, b2, b2), levels = sort(unique(c(b1, b2)))),
      condition = factor(c("control", "stress", "control", "stress"), levels = c("control", "stress"))
    )

    # All-zero genes cannot contribute to an effect and are omitted from model fitting.
    keep <- rowSums(X) > 0
    Xk <- X[keep, , drop = FALSE]
    dds <- DESeqDataSetFromMatrix(countData = Xk, colData = cd, design = ~ batch + condition)
    dds <- DESeq(dds, quiet = TRUE, minReplicatesForReplace = Inf)
    rn <- resultsNames(dds)
    cond_name <- grep("^condition_.*_vs_control$", rn, value = TRUE)
    batch_name <- grep("^batch_", rn, value = TRUE)
    if (length(cond_name) != 1 || length(batch_name) != 1) {
      stop(state, "/", cond, ": unexpected resultsNames: ", paste(rn, collapse = ","))
    }

    res_cond <- results(dds, name = cond_name, independentFiltering = FALSE, cooksCutoff = FALSE)
    res_batch <- results(dds, name = batch_name, independentFiltering = FALSE, cooksCutoff = FALSE)
    shr_cond <- lfcShrink(dds, coef = cond_name, type = "normal", quiet = TRUE)
    shr_batch <- lfcShrink(dds, coef = batch_name, type = "normal", quiet = TRUE)

    genes <- rownames(X)
    out <- data.frame(
      gene = genes,
      state = state,
      condition = cond,
      baseMean = NA_real_,
      condition_lfc_raw = 0,
      condition_lfc_shrunk = 0,
      condition_lfcSE = NA_real_,
      condition_pvalue = NA_real_,
      condition_padj = NA_real_,
      batch_lfc_raw = 0,
      batch_lfc_shrunk = 0,
      batch_lfcSE = NA_real_,
      stringsAsFactors = FALSE
    )
    idx <- which(keep)
    out$baseMean[idx] <- res_cond$baseMean
    out$condition_lfc_raw[idx] <- res_cond$log2FoldChange
    out$condition_lfc_shrunk[idx] <- shr_cond$log2FoldChange
    out$condition_lfcSE[idx] <- res_cond$lfcSE
    out$condition_pvalue[idx] <- res_cond$pvalue
    out$condition_padj[idx] <- res_cond$padj
    out$batch_lfc_raw[idx] <- res_batch$log2FoldChange
    out$batch_lfc_shrunk[idx] <- shr_batch$log2FoldChange
    out$batch_lfcSE[idx] <- res_batch$lfcSE

    all_out[[k]] <- out
    model_audit[[k]] <- data.frame(
      state = state,
      condition = cond,
      control_samples = paste(c1, c2, sep = ","),
      stress_samples = paste(s1, s2, sep = ","),
      batches = paste(b1, b2, sep = ","),
      n_modeled_genes = sum(keep),
      condition_coef = cond_name,
      batch_coef = batch_name,
      stringsAsFactors = FALSE
    )
    cat(cond, ": modeled genes=", sum(keep), " coef=", cond_name, " batch=", batch_name, "\n", sep = "")
    flush.console()
    rm(X, Xk, dds, res_cond, res_batch, shr_cond, shr_batch, out)
    gc(verbose = FALSE)
  }
}

res <- do.call(rbind, all_out)
audit <- do.call(rbind, model_audit)

# Preserve explicit zeros for genes that were all-zero in one 4-sample model; any
# model-estimation NA for a nonzero gene remains NA and is handled transparently downstream.
con <- gzfile(file.path(outdir, "CellAnnotation_DESeq2_effects.tsv.gz"), open = "wt")
write.table(res, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
close(con)
write.table(audit, file.path(outdir, "model_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nmodels=", nrow(audit), " rows=", nrow(res), "\n", sep = "")
cat("ROBERTSON2026_DESEQ2_SUCCESS\n")
sink()
close(zz)
