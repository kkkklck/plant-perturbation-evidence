# Dataset manifest

The [machine-readable manifest](dataset_manifest.tsv) records study, organism, assay, task, source/target roles and the operational preprocessing definitions. The [input requirements](../archive/manifests/input_requirements.tsv) list acquisition locations and local runtime filenames.

Main studies: GSE313726, GSE153103, GSE158761, GSE46205, GSE30098/GSE30166, GSE121619, GSE290176, GSE289716. The sulfur analysis operationally uses six GSE30098 samples; the related/historical GSE30166 label is not an extra validation cohort. Independent sample-to-series metadata should be checked when reacquiring inputs.

Exploratory: GSE253757, GSE130061, with SHR direct-target priors from GSE21338. The optional low-pH 48 h sensitivity uses Dinneny-series source samples documented by scripts 136/137; it is not a new primary cohort.

The preprocessing semantic audit is exported in `source_data/extended_data/results/publication_controls_v1/main_dataset_preprocessing_semantic_audit/`. It records seven PASS datasets and the corrected GSE46205 exception. Its checks reflect the archived project audit, not a fresh raw-data reanalysis during export.

Obtain upstream raw/processed inputs from GEO or the authors. Do not place them in version control. GEO accessions are stable identifiers, but a supplementary file's original name, checksum and scale must be checked before rerunning an analysis. Avoid substituting a transformed matrix for the named count layer or processed VALUE table.
