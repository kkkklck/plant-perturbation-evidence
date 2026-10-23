suppressPackageStartupMessages({
  library(Matrix)
})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
base <- file.path(root, "data", "external_lateral_root", "GSE158761")
inp <- file.path(base, "GSE158761_all_cells_labeled.rds")
meta_path <- file.path(base, "plantvc_v1_probe", "metadata.tsv.gz")
out <- file.path(base, "plantvc_v1_blind")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(out, "freeze.log")
zz <- file(log_path, open = "wt")
sink(zz, type = "output")
sink(zz, type = "message")
on.exit({
  while (sink.number(type = "message") > 2) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  try(close(zz), silent = TRUE)
}, add = TRUE)

ELIGIBLE <- c(
  "Ambiguous Stele Cells",
  "Columella/ Root Cap",
  "Epidermis",
  "Lateral Root Endodermis",
  "Mature Pericycle",
  "Phloem",
  "Phloem Pole Pericycle",
  "Xylem",
  "Xylem Pole Pericycle"
)
SLUG <- c(
  "Ambiguous Stele Cells" = "ambiguous_stele",
  "Columella/ Root Cap" = "columella_root_cap",
  "Epidermis" = "epidermis",
  "Lateral Root Endodermis" = "lateral_root_endodermis",
  "Mature Pericycle" = "mature_pericycle",
  "Phloem" = "phloem",
  "Phloem Pole Pericycle" = "phloem_pole_pericycle",
  "Xylem" = "xylem",
  "Xylem Pole Pericycle" = "xylem_pole_pericycle"
)
PRE_SAMPLES <- c(1L, 2L, 5L, 6L)
TARGET_SAMPLES <- c(3L, 7L)

cat("GSE158761_FREEZE_START\n")
cat("input_bytes=", file.info(inp)$size, "\n", sep = "")
meta <- read.delim(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
meta$sample_number <- as.integer(meta$sample_number)
stopifnot(all(c("cell","sample_number","experiment","celltype") %in% names(meta)))

# Re-derive metadata-only eligibility and require exact agreement with frozen list.
tab <- table(meta$celltype, meta$sample_number)
needed <- c("1","2","3","5","6","7")
if (!all(needed %in% colnames(tab))) stop("Expected six sample numbers are missing")
derived <- rownames(tab)[apply(tab[, needed, drop = FALSE] >= 20, 1, all)]
if (!setequal(derived, ELIGIBLE)) {
  stop("Frozen eligible cell types disagree with metadata-derived >=20 rule")
}
cat("eligible_celltypes=", length(ELIGIBLE), "\n", sep = "")
cat(paste(ELIGIBLE, collapse = " | "), "\n")

# Deserialize the author object. Target expression exists in the file, but all
# statistics below immediately subset to pretarget sample numbers 1,2,5,6.
o <- readRDS(inp)
a <- attr(o, "assays", exact = TRUE)
e <- attr(a, ".xData", exact = TRUE)
items <- as.list.environment(e, all.names = TRUE)
d <- items[["data"]]
ld <- attr(d, "listData", exact = TRUE)
E <- ld[["counts"]]
if (!inherits(E, "Matrix")) stop("Could not recover sparse counts matrix")
cat("counts_dim=", nrow(E), "x", ncol(E), "\n", sep = "")

cn <- colnames(E)
if (is.null(cn) || !setequal(cn, meta$cell)) stop("Count-matrix cells do not exactly match metadata cells")
meta <- meta[match(cn, meta$cell), , drop = FALSE]
stopifnot(identical(meta$cell, cn))
genes <- rownames(E)
if (is.null(genes)) stop("Counts matrix has no gene identifiers")
cat("CELL_MAPPING_OK\n")

pre_mask <- meta$sample_number %in% PRE_SAMPLES
target_mask <- meta$sample_number %in% TARGET_SAMPLES
if (any(pre_mask & target_mask)) stop("Pretarget/target sample overlap")
Epre <- E[, pre_mask, drop = FALSE]
mpre <- meta[pre_mask, , drop = FALSE]
# Drop the full object/counts before any pseudobulk statistics.
rm(E, o, a, e, items, d, ld)
gc()

used <- mpre[mpre$celltype %in% ELIGIBLE, c("cell","sample_number","experiment","celltype")]
withheld <- meta[target_mask & meta$celltype %in% ELIGIBLE, c("cell","sample_number","experiment","celltype")]
write.table(used, file.path(out, "freeze_cells_used.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(withheld, file.path(out, "target_cells_withheld.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("freeze_cells_used=", nrow(used), "\n", sep = "")
cat("target_20h_cells_withheld=", nrow(withheld), "\n", sep = "")
cat("target_20h_cells_used=0\n")

sum_group <- function(ct, sample_no) {
  idx <- which(mpre$celltype == ct & mpre$sample_number == sample_no)
  if (length(idx) < 20) stop("Unexpected low pretarget cell count for ", ct, " sample ", sample_no)
  as.numeric(Matrix::rowSums(Epre[, idx, drop = FALSE]))
}

resp <- list()
context_rows <- list()
for (ct in ELIGIBLE) {
  c1 <- sum_group(ct, 1L)
  c2 <- sum_group(ct, 2L)
  c5 <- sum_group(ct, 5L)
  c6 <- sum_group(ct, 6L)
  C <- cbind(c1, c2, c5, c6)
  libs <- colSums(C)
  if (any(libs <= 0)) stop("Zero pseudobulk library for ", ct)
  CPM <- sweep(C, 2, libs, "/") * 1e6
  L <- log2(CPM + 1)
  rA <- L[,2] - L[,1]
  rB <- L[,4] - L[,3]
  r8 <- (rA + rB) / 2
  technical <- rowSums(CPM >= 1) >= 2
  resp[[ct]] <- list(r8 = r8, rA = rA, rB = rB, technical = technical)
  context_rows[[ct]] <- data.frame(
    celltype = ct,
    slug = unname(SLUG[ct]),
    n_sample1 = sum(mpre$celltype == ct & mpre$sample_number == 1L),
    n_sample2 = sum(mpre$celltype == ct & mpre$sample_number == 2L),
    n_sample5 = sum(mpre$celltype == ct & mpre$sample_number == 5L),
    n_sample6 = sum(mpre$celltype == ct & mpre$sample_number == 6L),
    n_technical_genes = sum(technical),
    stringsAsFactors = FALSE
  )
}

R8 <- do.call(cbind, lapply(ELIGIBLE, function(ct) resp[[ct]]$r8))
colnames(R8) <- ELIGIBLE

for (ct in ELIGIBLE) {
  others <- setdiff(ELIGIBLE, ct)
  broad <- rowMeans(R8[, others, drop = FALSE])
  z <- resp[[ct]]
  ev <- data.frame(
    gene = genes,
    celltype = ct,
    technical_universe = z$technical,
    r8 = z$r8,
    r8_repA = z$rA,
    r8_repB = z$rB,
    broad_shared_r8 = broad,
    stringsAsFactors = FALSE
  )
  con <- gzfile(file.path(out, paste0("evidence_", unname(SLUG[ct]), ".tsv.gz")), open = "wt")
  write.table(ev, con, sep = "\t", quote = FALSE, row.names = FALSE)
  close(con)
}

contexts <- do.call(rbind, context_rows)
write.table(contexts, file.path(out, "contexts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("evidence_files_written=", length(ELIGIBLE), "\n", sep = "")
cat("TARGET_20H_EXPRESSION_NOT_USED_FOR_ANY_STATISTIC\n")
cat("GSE158761_FREEZE_R_SUCCESS\n")
