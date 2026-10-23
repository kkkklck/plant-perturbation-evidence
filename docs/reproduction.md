# Reproduction guide

## Environment
Use the recorded interpreter/package versions under `environments/`; do not upgrade packages merely to run this repository. `requirements-analysis.txt` is the scoped Python dependency list derived from the recorded environment. `recorded_pip_freeze.txt` is an unfiltered historical environment capture and contains unrelated installed tools; it is provenance, not a recommended universal installation command. The recorded R session is a base-only session and is not a complete R package lockfile. Robertson preprocessing logs record R 4.5.3 and Seurat 5.5.1; legacy RDS extraction uses different compatibility code. Full R environment restoration remains a documented limitation. Required R libraries are listed in `environments/R_package_inventory.txt`; this is an inventory, not a version lock.

## Tier 1: archived derived inputs
From the repository root:
```bash
python scripts/utilities/reproduce.py verify
python scripts/utilities/reproduce.py tables
python scripts/utilities/reproduce.py figures
python scripts/utilities/reproduce.py extended-data
```
The first command needs only the Python standard library. Table/figure commands copy allowlisted inputs and historical rule bytes into ignored runtime paths. Existing different files are rejected rather than overwritten. Figure tables are rebuilt and compared byte-for-byte. The existing final layouts are Figure 1 v15, Figures 2/3 v11 and Figure 4 v16. Earlier plotting modules are inherited code dependencies; they are not separate results. The two latest reviewer figures are rebuilt by scripts 166 and 167.

Figure commands write ignored working exports beside inherited modules and copy the selected final versions into `outputs/main_figures/`. Full-resolution figures and PDFs are generated locally and are excluded from Git. The root namespace for old analysis scripts is intentionally preserved to keep imports and historical script hashes valid.

## Tier 2: recompute from public study inputs
Acquire the exact named inputs in `archive/manifests/input_requirements.tsv`. The generic `scripts/preprocessing/download_geo_file.py` accepts one explicit NCBI HTTPS input URL and writes only under ignored `data/`, without overwriting files. The original Robertson download wrapper is also preserved. For the preprocessing semantic audit (script 163), reacquire the upstream GSE289716 author `scripts.zip` and the original preprocessing logs; the author script is not redistributed here. Keep source-only and target-only directories separate. Read each script's input constants and command-line mode before execution. `docs/workflow_inventory.tsv` lists categories and dependencies.

| Analysis | Ordered stages |
|---|---|
| Rice | 93 raw count-layer pseudobulk → 97 tissue-system aggregation → 98 tissue-system structure → 99 cross-replicate candidates → 100 candidate significance → 158 matched nulls → 164 expression-filter sensitivity |
| Salt bulk | 101 freeze then score for 6 h → 102 freeze then score for 12 h |
| Lateral root | 110/111 metadata and assay probes → 112 source freeze → 113 freeze/rank audit → 114 target export → 115 score |
| Heat held-out state | 103 R `freeze` → 104 Python `freeze` → 103 R `target` → 104 Python `score` |
| Cross-study heat | 118 freeze → 119 score → 148/150 controls → 157 matched nulls → 165 generic-stress sensitivity |
| GSE46205 | Historical 131–135 are retained for provenance; use 160 correction → 161 corrected statistics → 162 corrected suite for submission values |
| Sulfur | 138 freeze with 12 h absent → acquire/reveal 12 h → 139 score |
| Immune | 140 source firewall → 141 freeze → 142 target firewall → 143 score → 147 bootstrap/permutations/matched nulls |
| Publication synthesis | 149/150/152/155/156/157/158 controls as applicable → 160/161/162 correction → 159 main source tables |

Some scripts execute immediately on import and some freeze commands refuse to run with target files present. Do not import them merely to inspect functions or replay freeze in a workspace containing revealed targets. Follow the scientific split, not just file order. Baseline suite 155 also includes the optional PAT extension and requires its original inputs; the archived corrected suite can be used directly for Tier 1.

## Optional extensions
Original extension programs reside in `scripts/exploratory/`. Materialization copies their exact bytes to historical ignored paths under `scripts/`, which preserves frozen code hashes. After acquiring source/target inputs in a clean workspace, run historical scripts 116/117 for PAT or 144/145/153 for SHR. These analyses do not establish the primary conclusions.

## Boundaries of the export validation
Raw data and upstream expression objects are deliberately absent. No full raw-to-result rerun is claimed for this export. R package restoration and upstream availability are not certified by a table/figure smoke test. The authoritative snapshots are in `source_data/`; regenerated results go to ignored `data/` and `outputs/` paths. Semantic numerical checks and byte-level integrity checks are separate.
