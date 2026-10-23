suppressPackageStartupMessages({
  library(S4Vectors)
})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
base <- file.path(root, "data", "external_lateral_root", "GSE158761")
inp <- file.path(base, "GSE158761_all_cells_labeled.rds")
out <- file.path(base, "plantvc_v1_probe")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
log_path <- file.path(out, "probe.log")
zz <- file(log_path, open = "wt")
sink(zz, type = "output")
sink(zz, type = "message")
on.exit({
  while (sink.number(type = "message") > 2) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  try(close(zz), silent = TRUE)
}, add = TRUE)

cat("GSE158761_METADATA_PROBE_START\n")
cat("R=", R.version.string, "\n", sep = "")
cat("input_bytes=", file.info(inp)$size, "\n", sep = "")
flush.console()

o <- readRDS(inp)
cat("class=", paste(class(o), collapse = ","), "\n", sep = "")

# Avoid calling monocle3 methods. Access the serialized S4 colData slot directly.
cd <- attr(o, "colData", exact = TRUE)
if (is.null(cd)) stop("No raw colData attribute found")
ld <- attr(cd, "listData", exact = TRUE)
if (is.null(ld) || !is.list(ld)) stop("Could not access DFrame listData")
cn <- names(ld)
if (is.null(cn)) stop("colData listData has no names")

# Convert only ordinary atomic/factor columns; keep nested exotic columns out of the TSV.
keep <- vapply(ld, function(x) {
  is.atomic(x) || is.factor(x)
}, logical(1))
if (!any(keep)) stop("No atomic metadata columns found")
meta <- as.data.frame(ld[keep], stringsAsFactors = FALSE, check.names = FALSE)

# Recover row names / cell ids without invoking monocle3 methods.
rn <- attr(cd, "rownames", exact = TRUE)
if (is.null(rn)) rn <- rownames(meta)
if (is.null(rn) || length(rn) != nrow(meta)) rn <- as.character(seq_len(nrow(meta)))
meta <- cbind(cell = as.character(rn), meta, stringsAsFactors = FALSE)

cat("metadata_rows=", nrow(meta), "\n", sep = "")
cat("metadata_columns=", paste(names(meta), collapse = "|"), "\n", sep = "")

for (n in names(meta)) {
  x <- meta[[n]]
  u <- length(unique(x))
  if (u <= 50) {
    cat("\n-- ", n, " (n_unique=", u, ") --\n", sep = "")
    print(sort(table(x, useNA = "ifany"), decreasing = TRUE))
  }
}

con <- gzfile(file.path(out, "metadata.tsv.gz"), open = "wt")
write.table(meta, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)

cat("\nEXPRESSION_ASSAY_NOT_ACCESSED\n")
cat("GSE158761_METADATA_PROBE_SUCCESS\n")
