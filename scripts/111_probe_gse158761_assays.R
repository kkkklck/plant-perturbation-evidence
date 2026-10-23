suppressPackageStartupMessages({
  library(S4Vectors)
})

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
base <- file.path(root, "data", "external_lateral_root", "GSE158761")
inp <- file.path(base, "GSE158761_all_cells_labeled.rds")
out <- file.path(base, "plantvc_v1_probe")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
log_path <- file.path(out, "assay_probe.log")
zz <- file(log_path, open = "wt")
sink(zz, type = "output")
sink(zz, type = "message")
on.exit({
  while (sink.number(type = "message") > 2) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  try(close(zz), silent = TRUE)
}, add = TRUE)

cat("GSE158761_ASSAY_STRUCTURE_PROBE_START\n")
o <- readRDS(inp)
cat("class=", paste(class(o), collapse = ","), "\n", sep = "")
a <- attr(o, "assays", exact = TRUE)
if (is.null(a)) stop("No assays attribute")
cat("assays_class=", paste(class(a), collapse = "/"), "\n", sep = "")
e <- attr(a, ".xData", exact = TRUE)
if (!is.environment(e)) stop("assays .xData is not an environment")
items <- as.list.environment(e, all.names = TRUE)
cat("assay_env_names=", paste(names(items), collapse = "|"), "\n", sep = "")
for (nm in names(items)) {
  x <- items[[nm]]
  cat("ITEM ", nm, " class=", paste(class(x), collapse = "/"), " length=", length(x), sep = "")
  if (!is.null(dim(x))) cat(" dim=", paste(dim(x), collapse = "x"), sep = "")
  if (is.list(x)) cat(" names=", paste(names(x), collapse = "|"), sep = "")
  if (is.environment(x)) cat(" env_n=", length(as.list.environment(x, all.names = TRUE)), sep = "")
  cat("\n")
  at <- attributes(x)
  if (!is.null(at)) {
    cat("  attrs=", paste(names(at), collapse = "|"), "\n", sep = "")
  }
}

if ("data" %in% names(items)) {
  d <- items[["data"]]
  ld <- attr(d, "listData", exact = TRUE)
  if (is.list(ld)) {
    cat("data_listData_names=", paste(names(ld), collapse = "|"), "\n", sep = "")
    for (nm in names(ld)) {
      x <- ld[[nm]]
      cat("ASSAY ", nm, " class=", paste(class(x), collapse = "/"), " length=", length(x), sep = "")
      if (!is.null(dim(x))) cat(" dim=", paste(dim(x), collapse = "x"), sep = "")
      cat("\n")
    }
  }
}
cat("NO_EXPRESSION_VALUES_PRINTED\n")
cat("GSE158761_ASSAY_STRUCTURE_PROBE_SUCCESS\n")
