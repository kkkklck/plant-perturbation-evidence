suppressPackageStartupMessages({
  library(Matrix)
})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
base <- file.path(root, "data", "external_lateral_root", "GSE158761")
inp <- file.path(base, "GSE158761_all_cells_labeled.rds")
meta_path <- file.path(base, "plantvc_v1_probe", "metadata.tsv.gz")
out <- file.path(base, "plantvc_v1_blind")

log_path <- file.path(out, "target.log")
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

manifest_path <- file.path(out, "FREEZE_MANIFEST.json")
if (!file.exists(manifest_path)) stop("Freeze manifest missing; target cannot be opened")
cat("GSE158761_TARGET_START\n")
meta <- read.delim(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
meta$sample_number <- as.integer(meta$sample_number)

o <- readRDS(inp)
a <- attr(o, "assays", exact = TRUE)
e <- attr(a, ".xData", exact = TRUE)
items <- as.list.environment(e, all.names = TRUE)
d <- items[["data"]]
ld <- attr(d, "listData", exact = TRUE)
E <- ld[["counts"]]
if (!inherits(E, "Matrix")) stop("Could not recover sparse counts matrix")
cn <- colnames(E)
if (is.null(cn) || !setequal(cn, meta$cell)) stop("Count-matrix cells do not exactly match metadata cells")
meta <- meta[match(cn, meta$cell), , drop = FALSE]
genes <- rownames(E)
if (is.null(genes)) stop("Counts matrix has no gene identifiers")
cat("CELL_MAPPING_OK\n")

# Only baseline and target samples are needed for target scoring.
keep <- meta$sample_number %in% c(1L, 3L, 5L, 7L)
Et <- E[, keep, drop = FALSE]
mt <- meta[keep, , drop = FALSE]
rm(E, o, a, e, items, d, ld)
gc()

sum_group <- function(ct, sample_no) {
  idx <- which(mt$celltype == ct & mt$sample_number == sample_no)
  if (length(idx) < 20) stop("Unexpected low target/baseline cell count for ", ct, " sample ", sample_no)
  as.numeric(Matrix::rowSums(Et[, idx, drop = FALSE]))
}

for (ct in ELIGIBLE) {
  c1 <- sum_group(ct, 1L)
  c3 <- sum_group(ct, 3L)
  c5 <- sum_group(ct, 5L)
  c7 <- sum_group(ct, 7L)
  C <- cbind(c1, c3, c5, c7)
  libs <- colSums(C)
  if (any(libs <= 0)) stop("Zero pseudobulk library for ", ct)
  CPM <- sweep(C, 2, libs, "/") * 1e6
  L <- log2(CPM + 1)
  rA <- L[,2] - L[,1]
  rB <- L[,4] - L[,3]
  r20 <- (rA + rB) / 2
  z <- data.frame(
    gene = genes,
    celltype = ct,
    r20 = r20,
    r20_repA = rA,
    r20_repB = rB,
    stringsAsFactors = FALSE
  )
  con <- gzfile(file.path(out, paste0("target_", unname(SLUG[ct]), ".tsv.gz")), open = "wt")
  write.table(z, con, sep = "\t", quote = FALSE, row.names = FALSE)
  close(con)
}

cat("target_contexts_written=", length(ELIGIBLE), "\n", sep = "")
cat("GSE158761_TARGET_EXPORT_SUCCESS\n")
