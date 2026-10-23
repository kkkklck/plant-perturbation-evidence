suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
inp <- file.path(root, "data", "external_rice", "Robertson2026", "GSE313726_OsLeafStressIntegrated.rds")
outdir <- file.path(root, "data", "external_rice", "Robertson2026", "pseudobulk_v0")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(outdir, "build.log")
zz <- file(logf, open = "wt")
sink(zz, type = "output")
nsink <- TRUE
on.exit({
  if (nsink) sink()
  close(zz)
}, add = TRUE)

cat("ROBERTSON2026_PSEUDOBULK_START\n")
cat("R=", R.version.string, "\n", sep = "")
cat("Seurat=", as.character(packageVersion("Seurat")), "\n", sep = "")
cat("input_bytes=", file.info(inp)$size, "\n", sep = "")
flush.console()

obj <- readRDS(inp)
cat("READ_OK cells=", ncol(obj), " features=", nrow(obj), " object_GiB=", format(as.numeric(object.size(obj))/1024^3, digits=6), "\n", sep = "")

if (!identical(DefaultAssay(obj), "RNA")) DefaultAssay(obj) <- "RNA"
ly <- Layers(obj[["RNA"]])
cat("RNA_layers=", paste(ly, collapse=","), "\n", sep="")
if (!"counts" %in% ly) stop("RNA counts layer missing")

md <- obj[[]]
req <- c("Identity", "Condition", "CellAnnotation", "TissueSystem")
miss <- setdiff(req, colnames(md))
if (length(miss)) stop("Missing metadata columns: ", paste(miss, collapse=","))

cat("n_samples=", length(unique(md$Identity)), "\n", sep="")
cat("n_conditions=", length(unique(md$Condition)), "\n", sep="")
cat("n_cell_annotations=", length(unique(md$CellAnnotation)), "\n", sep="")
cat("n_tissue_systems=", length(unique(md$TissueSystem)), "\n", sep="")

sample_condition <- unique(data.frame(sample=as.character(md$Identity), condition=as.character(md$Condition), stringsAsFactors=FALSE))
if (anyDuplicated(sample_condition$sample)) stop("A sample maps to multiple conditions")
annotation_tissue <- unique(data.frame(cell_annotation=as.character(md$CellAnnotation), tissue_system=as.character(md$TissueSystem), stringsAsFactors=FALSE))
if (anyDuplicated(annotation_tissue$cell_annotation)) stop("A CellAnnotation maps to multiple TissueSystem values")

group_key <- paste(as.character(md$Identity), as.character(md$CellAnnotation), sep = "||")
groups <- sort(unique(group_key))
group_fac <- factor(group_key, levels = groups)

Z <- sparseMatrix(i=seq_along(group_fac), j=as.integer(group_fac), x=1,
                  dims=c(length(group_fac), length(groups)), dimnames=list(rownames(md), groups))

counts <- LayerData(obj, assay="RNA", layer="counts")
if (!identical(colnames(counts), rownames(md))) stop("RNA count columns do not match metadata row order")
cat("counts_class=", paste(class(counts), collapse=","), " nnzero=", nnzero(counts), "\n", sep="")

pb_counts <- counts %*% Z
cat("PSEUDOBULK_COUNTS_OK dims=", nrow(pb_counts), "x", ncol(pb_counts), " nnzero=", nnzero(pb_counts), "\n", sep="")

parts <- strsplit(groups, "||", fixed=TRUE)
group_meta <- data.frame(
  group=groups,
  sample=vapply(parts, `[[`, character(1), 1),
  cell_annotation=vapply(parts, `[[`, character(1), 2),
  stringsAsFactors=FALSE
)
group_meta$condition <- sample_condition$condition[match(group_meta$sample, sample_condition$sample)]
group_meta$tissue_system <- annotation_tissue$tissue_system[match(group_meta$cell_annotation, annotation_tissue$cell_annotation)]
if (anyNA(group_meta$condition) || anyNA(group_meta$tissue_system)) stop("Failed to map condition or tissue system to groups")

cell_tab <- table(group_fac)
group_meta$n_cells <- as.integer(cell_tab[groups])
group_meta$total_umi <- as.numeric(colSums(pb_counts))
group_meta$n_detected_genes <- as.integer(colSums(pb_counts > 0))

write.table(group_meta, file.path(outdir, "group_metadata.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(sample_condition, file.path(outdir, "sample_condition.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(annotation_tissue, file.path(outdir, "annotation_tissue.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(data.frame(gene=rownames(pb_counts)), file.path(outdir, "genes.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
saveRDS(pb_counts, file.path(outdir, "pseudobulk_raw_counts_sparse.rds"), compress=TRUE)

lib <- colSums(pb_counts)
if (any(lib <= 0)) stop("At least one pseudobulk group has zero total UMI")
cpm <- sweep(as.matrix(pb_counts), 2, lib / 1e6, FUN="/")
logcpm <- log1p(cpm)
rm(cpm)
gc()

outmat <- data.frame(gene=rownames(logcpm), logcpm, check.names=FALSE)
con <- gzfile(file.path(outdir, "pseudobulk_log1p_cpm.tsv.gz"), open="wt")
write.table(outmat, con, sep="\t", quote=FALSE, row.names=FALSE)
close(con)
rm(outmat, logcpm, counts, Z, obj)
gc()

coverage <- xtabs(n_cells ~ cell_annotation + sample, data=group_meta)
write.table(as.data.frame.matrix(coverage), file.path(outdir, "cell_counts_by_annotation_sample.tsv"), sep="\t", quote=FALSE, col.names=NA)

cat("GROUPS=", nrow(group_meta), "\n", sep="")
cat("GROUP_CELL_MIN=", min(group_meta$n_cells), "\n", sep="")
cat("GROUP_CELL_MEDIAN=", median(group_meta$n_cells), "\n", sep="")
cat("GROUP_CELL_MAX=", max(group_meta$n_cells), "\n", sep="")
cat("OUTPUT_RAW_RDS=", file.path(outdir, "pseudobulk_raw_counts_sparse.rds"), "\n", sep="")
cat("OUTPUT_LOGCPM=", file.path(outdir, "pseudobulk_log1p_cpm.tsv.gz"), "\n", sep="")
cat("ROBERTSON2026_PSEUDOBULK_SUCCESS\n")

nsink <- FALSE
sink()
close(zz)
