# Recorded environments

- `python_version.txt`: archived Python version.
- `recorded_pip_freeze.txt`: historical full environment capture; includes unrelated tools. No upgrades were performed for export.
- `requirements-analysis.txt`: direct analysis/figure packages with versions extracted from that capture.
- `recorded_R_sessionInfo.txt`: archived R 4.5.3 base session; personal installation path and corrupt STDERR removed. It does not lock loaded analysis packages.
- `import_inventory.json`: Python import inventory of the selected analysis programs.

Robertson build logs independently record Seurat 5.5.1. R preprocessing also uses Matrix and SeuratObject. Complete historical R package locks are unavailable in the recorded session and must not be inferred from the current machine. Do not silently substitute newer packages for missing versions. The export smoke checks use the existing Python environment and do not establish raw-data R reproducibility.
