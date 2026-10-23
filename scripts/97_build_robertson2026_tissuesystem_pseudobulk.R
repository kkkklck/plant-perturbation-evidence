suppressPackageStartupMessages({
  library(Matrix)
})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
base <- file.path(root, "data", "external_rice", "Robertson2026", "pseudobulk_v0")
outdir <- file.path(root, "data", "external_rice", "Robertson2026", "tissuesystem_pseudobulk_v0")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

pb <- readRDS(file.path(base, "pseudobulk_raw_counts_sparse.rds"))
meta <- read.delim(file.path(base, "group_metadata.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
if (!identical(colnames(pb), meta$group)) stop("pseudobulk columns do not match group metadata")

key <- paste(meta$sample, meta$tissue_system, sep = "||")
groups <- sort(unique(key))
fac <- factor(key, levels = groups)
Z <- sparseMatrix(i = seq_along(fac), j = as.integer(fac), x = 1,
                  dims = c(length(fac), length(groups)),
                  dimnames = list(meta$group, groups))
ts_counts <- pb %*% Z
parts <- strsplit(groups, "||", fixed = TRUE)
ts_meta <- data.frame(
  group = groups,
  sample = vapply(parts, `[[`, character(1), 1),
  tissue_system = vapply(parts, `[[`, character(1), 2),
  stringsAsFactors = FALSE
)
sample_condition <- unique(meta[c("sample", "condition")])
ts_meta$condition <- sample_condition$condition[match(ts_meta$sample, sample_condition$sample)]
cell_tab <- tapply(meta$n_cells, key, sum)
ts_meta$n_cells <- as.integer(cell_tab[groups])
ts_meta$total_umi <- as.numeric(colSums(ts_counts))

saveRDS(ts_counts, file.path(outdir, "tissuesystem_raw_counts_sparse.rds"), compress = TRUE)
write.table(ts_meta, file.path(outdir, "group_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(gene = rownames(ts_counts)), file.path(outdir, "genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

lib <- colSums(ts_counts)
cpm <- sweep(as.matrix(ts_counts), 2, lib / 1e6, FUN = "/")
logcpm <- log1p(cpm)
outmat <- data.frame(gene = rownames(logcpm), logcpm, check.names = FALSE)
con <- gzfile(file.path(outdir, "tissuesystem_log1p_cpm.tsv.gz"), open = "wt")
write.table(outmat, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)

cat("ROBERTSON2026_TISSUESYSTEM_PSEUDOBULK_SUCCESS\n")
cat("groups=", nrow(ts_meta), " genes=", nrow(ts_counts), " min_cells=", min(ts_meta$n_cells), " median_cells=", median(ts_meta$n_cells), "\n", sep = "")
print(xtabs(n_cells ~ tissue_system + sample, data = ts_meta))
