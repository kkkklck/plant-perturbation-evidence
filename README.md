# Plant perturbation evidence structure

## Overview
This repository contains analyses supporting a task-aware evaluation framework for plant perturbation response prediction. The work asks which observed evidence remains transferable when different components of the target experiment are missing.

## Main question
How does the missing structure of a prediction task determine which source evidence can transfer, and which output resolution that evidence supports?

## Conceptual framework
Prediction target → evidence eligibility → output resolution.
Future-time, missing-state and cross-study tasks expose different information. Source capability further constrains whether full response vectors, gene candidates or candidate directions are supported. The contribution is an evaluation framework and its empirical analyses; the repository does not introduce a new universal predictor or deep-learning architecture.

## Main analyses
1. **Candidate reproducibility (Figure 1):** Robertson rice stress responses; biological-replicate vector stability, Top100 recovery, source-matched nulls and expression-filter sensitivity.
2. **Future-time prediction (Figure 2):** salt, lateral-root and sulfur experiments; recent same-state history, temporal baselines, matched nulls and exact context-level sign-flip tests.
3. **Missing-state and cross-study transfer (Figure 3):** held-out heat-shock states and independent-study transfer; source-matched nulls, alternative donors and wrong-perturbation controls.
4. **Source capability and output resolution (Figure 4):** source-only stratification of five immune stimuli, twenty nested target states, hierarchical bootstrap and exact stimulus-label permutation. Five stimuli are the independent units; geometry and response amplitude covary.

## Datasets
Primary accessions: GSE313726, GSE153103, GSE158761, GSE46205, GSE30098/GSE30166, GSE121619, GSE290176 and GSE289716. See [dataset manifest](docs/dataset_manifest.md) for assays, source/target roles and preprocessing. GSE253757 and GSE130061 are exploratory extensions.

## Repository structure
```text
archive/       Export provenance, input requirements and package checksums
config/        Frozen rules, artifact mapping and configuration examples
docs/          Scientific design, reproduction guide and provenance
environments/  Recorded versions and scoped requirements
scripts/       Original analysis namespace; figures, utilities and extensions
source_data/   Curated derived results, frozen predictions and plot inputs
```
Analysis filenames and historical result-directory labels are retained where code or frozen hashes depend on them. Historical labels are not the formal method name. The inherited figure modules are dependencies of the final renderers, not additional analyses.

## Reproduction
Run from the repository root with the recorded Python environment. No raw data are needed for these steps:
```bash
python scripts/utilities/reproduce.py verify
python scripts/utilities/reproduce.py tables
python scripts/utilities/reproduce.py figures
python scripts/utilities/reproduce.py extended-data
```
The table command verifies all four rebuilt main-figure tables against archived bytes. The figure command exports the existing Figure 1–4 layouts under `outputs/main_figures/`. The Extended Data command rebuilds the two latest reviewer-sensitivity figures. Runtime inputs are copied into ignored directories; archived source tables and frozen predictions are never overwritten.

Raw-to-result reanalysis is a separate tier: acquire the original inputs, reproduce preprocessing, and run the listed analysis stages. See [reproduction guide](docs/reproduction.md) and [workflow inventory](docs/workflow_inventory.tsv). Existing freeze commands must run in a fresh workspace with targets withheld: they are not safe commands to replay over archived freezes. Raw-input reanalysis has not been exhaustively rerun as part of this repository export.

## Freeze–reveal design
Source-only information, target sealing, frozen predictions and scoring are distinguished from later publication controls. The split is an analytical information-access restriction, not fully prospective clinical-style validation. Target gene identifiers or metadata were visible in some designs while numeric target expression was sealed. See [freeze–reveal design](docs/freeze_reveal_design.md).

## Important provenance note
### GSE46205 `corrected_post_reveal`
GEO VALUE is already **log2 GC-RMA signal**. Historical analyses applied an extra log2 transformation. Submission results were recomputed using deposited VALUE directly; treatment-minus-control responses therefore remain on the native log2 scale. The corrected results are supporting evidence labeled `corrected_post_reveal` and must not be called untouched prospective validation. Original erroneous artifacts and their hashes are retained separately. See [provenance notes](docs/provenance_notes.md).

## Exploratory extensions
PAT combination (GSE253757) and the SHR mechanistic prior (GSE130061; independent prior from GSE21338) are separated under `scripts/exploratory/` and `source_data/extensions/`.
These analyses are exploratory extensions and are not used to establish the primary task-structure conclusions.
The historical combined baseline suite includes PAT rows; those rows do not enter the primary Figure 1–4 task-structure claims.

## Data availability
Acquire source experiments from their public repositories using the accessions and input filenames in the manifests. Raw GEO matrices, expression objects, count matrices, manuscript documents and PDFs are excluded. Included tables contain study-derived statistical summaries, rankings, frozen predictions, null summaries and small visualization coordinates/metadata. No individual exported file exceeds 50 MB.

## Author and correspondence
**Congkai Liao** — College of Life Sciences, Capital Normal University, Beijing, China.  
ORCID: [0009-0004-6692-5919](https://orcid.org/0009-0004-6692-5919)  
Email: [congkailiao@gmail.com](mailto:congkailiao@gmail.com)

Associated manuscript (preprint in preparation): *Task structure determines transferable evidence for plant perturbation prediction*. This is a single-author research project.

## Citation
For reproducibility, cite the archived **v1.0.0-preprint** source release: [Zenodo DOI: 10.5281/zenodo.23973340](https://doi.org/10.5281/zenodo.23973340). The corresponding GitHub tag is [`v1.0.0-preprint`](https://github.com/kkkklck/plant-perturbation-evidence/releases/tag/v1.0.0-preprint). For current code metadata, see [CITATION.cff](CITATION.cff). The associated manuscript is being prepared for bioRxiv; its own preprint DOI has **not** been assigned yet. Metadata edits to `main` after the release do not change the Zenodo archived code snapshot.

## License
Repository-owned code and documentation are provided under the MIT License. Study-derived tables and visualization metadata retain the applicable upstream terms; this license does not relicense upstream data or dependencies. See [third-party notices](THIRD_PARTY_NOTICES.md).
