# Provenance notes

## GSE46205: corrected_post_reveal
The deposited GEO VALUE field is **log2 GC-RMA signal**. Historical scripts 131–135 applied an additional log2 transformation before calculating responses. That was an implementation error, not a scientific preprocessing choice.

Submission-level calculations use VALUE directly on the deposited scale; response is the treatment mean minus the matched control mean. Scripts 160, 161 and 162 independently produce the correction audit, corrected publication statistics and the corrected baseline suite. Script 159 reads corrected GSE46205 sources for the Figure 2 table.

Historical results remain under `source_data/figure2/results/prospective_router/GSE46205_salt/`. Corrected results remain under the separate `source_data/figure2/results/publication_controls_v1/gse46205_*` paths; the corrected combined suite is under `source_data/extended_data/results/publication_baseline_suite_v2_corrected/`. These results may support the conclusion but must be labeled `corrected_post_reveal`, not untouched prospective validation. No original erroneous artifact was deleted or overwritten.

## Export transformations
- Original analysis scripts retain their original bytes and numeric filenames. Exploratory scripts are stored separately and materialized to historical locations when needed.
- Seventeen scientific rule documents are preserved under `config/frozen_parameters/`; runtime copies keep original filenames to satisfy hash-bound code. They are scientific design records, not chat history.
- Nine lateral-root ranking manifests expose project-relative input paths instead of personal absolute paths. Original and exported hashes differ and are both recorded.
- Figure modules retain the existing inheritance chain, with repository-relative input/tool paths. The small local alignment gate is independently implemented for the explicit plot-rectangle groups in these figures; it is not a vendored skill script.
- R environment metadata have their personal installation path removed and corrupt shell STDERR excluded. Package and runtime versions are retained.

## Scientific limits
Source geometry and response amplitude covary in the five-stimulus analysis. A sample of five stimuli does not identify an independent geometric effect. The exact one-sided label-permutation result is coarse (ten assignments of a three-versus-two split); twenty nested target states must not be treated as twenty independent stimuli.

Post-reveal controls and later reviewer sensitivities are clearly distinct from source-only frozen predictions. Historical naming is preserved only where needed to reproduce artifacts and is not a claim that the paper introduces a virtual-cell model.
