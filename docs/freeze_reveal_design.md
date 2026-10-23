# Freeze–reveal design

## Information stages
1. **Source-only construction:** use permitted times, donor states or independent priors. Define the gene/probe universe, scoring rules, candidate counts and comparisons.
2. **Sealed target:** keep designated target expression values or response tables inaccessible to source fitting. This is a recorded analytical firewall; some original study objects contain both source and target values, so the code's indexing restrictions matter.
3. **Frozen predictions:** preserve rankings/predictions, parameter/rule files, source hashes and script hashes before scoring.
4. **Reveal:** open the target expression or author-response file after the freeze.
5. **Scoring:** evaluate fixed predictions on the original technical universe; retain complete-vector and candidate/direction endpoints separately.
6. **Publication controls:** later matched nulls, stronger baselines, hierarchical/context-level inference and reviewer sensitivities are post-reveal analyses.
7. **Correction:** a discovered implementation error triggers an independent correction record and a downgraded provenance label, never replacement of the historical freeze.

## Task-specific boundaries
| Task | Source information | Target information withheld | Qualification |
|---|---|---|---|
| Salt 6 h / 12 h | Earlier expression times and fixed rules | Respective later expression columns | Sequential freezes; not two independent studies |
| Lateral-root 20 h | 0/8 h counts and eligible cell metadata | 20 h numerical expression for scoring | Original processed object and metadata exist locally; code enforces access restrictions |
| Heat-shock missing state | Donor clusters 5, 6, 8; permitted technical/control information | Heat-shock cells in target clusters 2, 3, 10 | Controls/technical eligibility do not constitute sealed heat-shock outcomes |
| Independent heat transfer | Frozen GSE121619 program | GSE290176 expression values | Geneid fields and target file hash were inspected at freeze; numeric target fields were not parsed |
| Five immune stimuli | Ten named source states and source-only routing diagnostics | Four target response states per stimulus | Stimuli are the independent units; targets are nested |
| GSE46205 | Original historical task splits and rules | Original scoring targets before their historical reveal | Numerical correction is `corrected_post_reveal` |

## Archive and checksum scope
Original snapshots remain in the research archive. Exported numerical tables and frozen predictions are byte-identical except nine lateral-root ranking manifests whose absolute `input` path field was made relative. Their numerical fields and recorded historical hashes remain unchanged; both original and exported file hashes are recorded. A redacted export manifest must not be presented as the untouched original JSON byte stream.

`archive/manifests/export_manifest.json` records source-relative identity, original SHA-256, exported SHA-256 and whether a transformation occurred. `config/runtime_artifacts.json` maps curated artifacts to ignored runtime locations. `archive/checksums/SHA256SUMS` verifies the distributable package, excluding itself. These are export-time integrity records, not evidence of an independently timestamped prospective protocol.

Materialization copies original scientific rule bytes and exploratory script bytes into their historical runtime locations so historical script/rule hash checks continue to work. Freeze operations generate new artifacts and must be run in a clean reanalysis copy, never over preserved snapshots.
