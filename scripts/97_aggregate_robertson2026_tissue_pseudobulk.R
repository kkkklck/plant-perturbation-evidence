suppressPackageStartupMessages(library(Matrix))

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
base <- file.path(root, "data", "external_rice", "Robertson2026", "pseudobulk_v0")
outdir <- file.path(root, "data", "external_rice", "Robertson2026", "tissue_pseudobulk_v0")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

meta <- read.delim(file.path(base, "group_metadata.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
pb <- readRDS(file.path(base, "pseudobulk_raw_counts_sparse.rds"))
if (!identical(colnames(pb), meta$group)) stop("pseudobulk columns do not match group metadata")

keys <- paste(meta$sample, meta$tissue_system, sep = "||")
levels <- sort(unique(keys))
f <- factor(keys, levels = levels)
Z <- sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1,
                  dims = c(length(f), length(levels)),
                  dimnames = list(meta$group, levels))
agg <- pb %*% Z
parts <- strsplit(levels, "||", fixed = TRUE)
meta2 <- data.frame(
  group = levels,
  sample = vapply(parts, `[[`, character(1), 1),
  tissue_system = vapply(parts, `[[`, character(1), 2),
  stringsAsFactors = FALSE
)
sc <- unique(meta[, c("sample", "condition")])
meta2$condition <- sc$condition[match(meta2$sample, sc$sample)]
meta2$n_cells <- as.integer(tapply(meta$n_cells, f, sum)[levels])
meta2$total_umi <- as.numeric(colSums(agg))

write.table(meta2, file.path(outdir, "group_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(agg, file.path(outdir, "pseudobulk_raw_counts_sparse.rds"), compress = TRUE)

lib <- colSums(agg)
if (any(lib <= 0)) stop("zero-library tissue pseudobulk group")
cpm <- sweep(as.matrix(agg), 2, lib / 1e6, FUN = "/")
logcpm <- log1p(cpm)
out <- data.frame(gene = rownames(logcpm), logcpm, check.names = FALSE)
con <- gzfile(file.path(outdir, "pseudobulk_log1p_cpm.tsv.gz"), "wt")
write.table(out, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)

cat("ROBERTSON2026_TISSUE_PSEUDOBULK_SUCCESS\n")
cat("groups=", ncol(agg), " tissues=", length(unique(meta2$tissue_system)), " samples=", length(unique(meta2$sample)), "\n", sep = "")
cat("cell_min=", min(meta2$n_cells), " median=", median(meta2$n_cells), " max=", max(meta2$n_cells), "\n", sep = "")
