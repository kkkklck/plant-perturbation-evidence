suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
inp <- file.path(root, "data", "external_rice", "Robertson2026", "GSE313726_OsLeafStressIntegrated.rds")
outdir <- file.path(root, "data", "external_rice", "Robertson2026", "probe")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(outdir, "probe.log")
zz <- file(logf, open = "wt")
sink(zz, type = "output")
nsink <- TRUE
on.exit({ if (nsink) sink(); close(zz) }, add = TRUE)

cat("ROBERTSON2026_PROBE_START\n")
cat("R=", R.version.string, "\n", sep = "")
cat("Seurat=", as.character(packageVersion("Seurat")), "\n", sep = "")
cat("SeuratObject=", as.character(packageVersion("SeuratObject")), "\n", sep = "")
cat("input=", inp, "\n", sep = "")
cat("input_bytes=", file.info(inp)$size, "\n", sep = "")
flush.console()

obj <- readRDS(inp)
cat("READ_OK\n")
cat("class=", paste(class(obj), collapse = ","), "\n", sep = "")
cat("cells=", ncol(obj), "\n", sep = "")
cat("features_default_assay=", nrow(obj), "\n", sep = "")
cat("object_GiB=", format(as.numeric(object.size(obj))/1024^3, digits = 6), "\n", sep = "")
cat("default_assay=", DefaultAssay(obj), "\n", sep = "")
cat("assays=", paste(Assays(obj), collapse = ","), "\n", sep = "")
cat("reductions=", paste(Reductions(obj), collapse = ","), "\n", sep = "")
cat("graphs=", paste(names(obj@graphs), collapse = ","), "\n", sep = "")

for (a in Assays(obj)) {
  aa <- obj[[a]]
  cat("ASSAY ", a, " class=", paste(class(aa), collapse = ","), " features=", nrow(aa), " cells=", ncol(aa), "\n", sep = "")
  ly <- tryCatch(Layers(aa), error = function(e) character(0))
  cat("LAYERS ", a, " = ", paste(ly, collapse = ","), "\n", sep = "")
}

md <- obj[[]]
md$.__cell__ <- rownames(md)
write.table(md, gzfile(file.path(outdir, "metadata.tsv.gz"), "wt"), sep = "\t", quote = FALSE, row.names = FALSE)

summ <- lapply(names(md), function(nm) {
  x <- md[[nm]]
  ux <- unique(as.character(x))
  data.frame(
    column = nm,
    class = paste(class(x), collapse = ","),
    n_unique = length(ux),
    example_values = paste(head(ux, 12), collapse = " | "),
    stringsAsFactors = FALSE
  )
})
summ <- do.call(rbind, summ)
write.table(summ, file.path(outdir, "metadata_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("META_COLUMNS=", paste(names(md), collapse = ","), "\n", sep = "")
print(summ, row.names = FALSE)

candidates <- intersect(c("orig.ident","sample","Sample","Condition","condition","Treatment","treatment","Stress","stress","TissueSystem","CellAnnotation","celltype","cell_type"), names(md))
for (nm in candidates) {
  tab <- sort(table(md[[nm]], useNA = "ifany"), decreasing = TRUE)
  cat("\nMETA_TABLE ", nm, "\n", sep = "")
  print(tab)
}

cat("ROBERTSON2026_PROBE_SUCCESS\n")
nsink <- FALSE
sink()
close(zz)
