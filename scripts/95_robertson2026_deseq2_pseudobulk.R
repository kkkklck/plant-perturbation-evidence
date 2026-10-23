suppressPackageStartupMessages({
  library(DESeq2)
  library(Matrix)
})

root <- normalizePath(".", winslash="/", mustWork=TRUE)
base <- file.path(root, "data", "external_rice", "Robertson2026")
pbdir <- file.path(base, "pseudobulk_v0")
outdir <- file.path(base, "deseq2_v0")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

logf <- file.path(outdir, "deseq2.log")
zz <- file(logf, open="wt")
sink(zz, type="output")
nsink <- TRUE
on.exit({ if (nsink) sink(); close(zz) }, add=TRUE)

cat("ROBERTSON2026_DESEQ2_START\n")
cat("R=", R.version.string, "\n", sep="")
cat("DESeq2=", as.character(packageVersion("DESeq2")), "\n", sep="")

counts <- readRDS(file.path(pbdir, "pseudobulk_raw_counts_sparse.rds"))
meta <- read.delim(file.path(pbdir, "group_metadata.tsv"), check.names=FALSE, stringsAsFactors=FALSE)
design <- read.delim(file.path(base, "sample_design.tsv"), check.names=FALSE, stringsAsFactors=FALSE)
if (!identical(colnames(counts), meta$group)) stop("pseudobulk count columns do not match group metadata")

stress_map <- data.frame(
  label=c("Mild.Drought","Moderate.Drought","Mild.Salinity","Moderate.Salinity"),
  design_condition=c("MildDrought","ModerateDrought","MildSalinity","ModerateSalinity"),
  stringsAsFactors=FALSE
)
annotations <- sort(unique(meta$cell_annotation))

control_design <- design[design$stress_type == "control", , drop=FALSE]
control_by_batch <- setNames(control_design$identity, as.character(control_design$batch))

sig_rows <- list()
summary_rows <- list()
idx_sig <- 1L
idx_sum <- 1L

run_one <- function(annotation, stress_label, design_condition, mode) {
  stress_samples <- design$identity[design$condition == design_condition]
  if (length(stress_samples) != 2) stop("Expected 2 stress samples for ", design_condition)

  if (mode == "author_like") {
    ctrl_samples <- design$identity[design$stress_type == "control"]
    samples <- c(ctrl_samples, stress_samples)
    condition <- factor(c(rep("control", length(ctrl_samples)), rep("stress", length(stress_samples))), levels=c("control","stress"))
    coldata <- data.frame(row.names=samples, condition=condition)
    formula_use <- ~ condition
  } else if (mode == "batch_matched") {
    srows <- design[match(stress_samples, design$identity), , drop=FALSE]
    ctrl_samples <- unname(control_by_batch[as.character(srows$batch)])
    if (anyNA(ctrl_samples)) stop("Missing same-batch control")
    samples <- as.vector(rbind(ctrl_samples, stress_samples))
    batch <- factor(rep(as.character(srows$batch), each=2))
    condition <- factor(rep(c("control","stress"), times=nrow(srows)), levels=c("control","stress"))
    coldata <- data.frame(row.names=samples, batch=batch, condition=condition)
    formula_use <- ~ batch + condition
  } else stop("Unknown mode")

  groups <- paste(samples, annotation, sep="||")
  if (!all(groups %in% colnames(counts))) stop("Missing pseudobulk group for ", annotation, " / ", stress_label)
  mat <- as.matrix(counts[, groups, drop=FALSE])
  storage.mode(mat) <- "integer"
  colnames(mat) <- samples

  dds <- DESeqDataSetFromMatrix(countData=mat, colData=coldata, design=formula_use)
  dds <- DESeq(dds, quiet=TRUE)
  res <- results(dds, contrast=c("condition","stress","control"), independentFiltering=TRUE)
  rdf <- as.data.frame(res)
  rdf$gene <- rownames(rdf)
  rdf <- rdf[, c("gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")]
  rdf$annotation <- annotation
  rdf$condition <- stress_label
  rdf$mode <- mode
  rdf$strong <- !is.na(rdf$padj) & rdf$padj < 0.05 & abs(rdf$log2FoldChange) > 0.5

  n_tested <- sum(!is.na(rdf$padj))
  n_fdr <- sum(!is.na(rdf$padj) & rdf$padj < 0.05)
  n_strong <- sum(rdf$strong)
  n_up <- sum(rdf$strong & rdf$log2FoldChange > 0)
  n_down <- sum(rdf$strong & rdf$log2FoldChange < 0)

  list(
    sig=rdf[!is.na(rdf$padj) & rdf$padj < 0.05, , drop=FALSE],
    summary=data.frame(
      mode=mode, condition=stress_label, cell_annotation=annotation,
      n_samples=length(samples), n_tested=n_tested, n_fdr05=n_fdr,
      n_fdr05_abs_lfc_gt_0_5=n_strong, n_up=n_up, n_down=n_down,
      stringsAsFactors=FALSE
    )
  )
}

for (mode in c("author_like","batch_matched")) {
  cat("MODE ", mode, "\n", sep="")
  for (i in seq_len(nrow(stress_map))) {
    cond <- stress_map$label[i]
    dcond <- stress_map$design_condition[i]
    cat(" CONDITION ", cond, "\n", sep="")
    for (ann in annotations) {
      cat("  ", ann, " ... ", sep="")
      flush.console()
      ans <- run_one(ann, cond, dcond, mode)
      sig_rows[[idx_sig]] <- ans$sig; idx_sig <- idx_sig + 1L
      summary_rows[[idx_sum]] <- ans$summary; idx_sum <- idx_sum + 1L
      cat("FDR05=", ans$summary$n_fdr05, " strong=", ans$summary$n_fdr05_abs_lfc_gt_0_5, "\n", sep="")
    }
  }
}

sig <- do.call(rbind, sig_rows)
summ <- do.call(rbind, summary_rows)
write.table(summ, file.path(outdir, "deseq2_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
con <- gzfile(file.path(outdir, "deseq2_fdr05_genes.tsv.gz"), open="wt")
write.table(sig, con, sep="\t", quote=FALSE, row.names=FALSE)
close(con)

# Per-condition union/intersection descriptors, useful for downstream descriptive structure.
panel_rows <- list(); ip <- 1L
for (mode in unique(sig$mode)) {
  for (cond in unique(sig$condition)) {
    sub <- sig[sig$mode == mode & sig$condition == cond & sig$strong, , drop=FALSE]
    if (!nrow(sub)) next
    tab <- sort(table(sub$gene), decreasing=TRUE)
    panel_rows[[ip]] <- data.frame(
      mode=mode, condition=cond, gene=names(tab),
      n_annotations_strong=as.integer(tab),
      stringsAsFactors=FALSE
    )
    ip <- ip + 1L
  }
}
panels <- do.call(rbind, panel_rows)
write.table(panels, file.path(outdir, "condition_deg_union_membership.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

cat("\nSUMMARY_STRONG_DEG_COUNTS\n")
print(aggregate(n_fdr05_abs_lfc_gt_0_5 ~ mode + condition, data=summ, FUN=function(x) c(median=median(x), mean=mean(x), min=min(x), max=max(x))))
cat("ROBERTSON2026_DESEQ2_SUCCESS\n")

nsink <- FALSE
sink()
close(zz)
