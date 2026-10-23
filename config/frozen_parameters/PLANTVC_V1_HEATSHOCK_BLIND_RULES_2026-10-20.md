# PlantVC v1 — Arabidopsis 45-min heat-shock blind rules

Date frozen: 2026-10-20
Dataset: GSE121619

## Goal
Test the v1 acute-response hypothesis on a genuinely new single-cell perturbation dataset:

> For a short acute perturbation, a cross-state shared response estimated from other cell states can prioritize the strongest response genes in completely held-out target states.

This benchmark tests candidate prioritization and response direction, not exact full-transcriptome reconstruction and not causal mechanism.

## Metadata-only selection performed before reading expression values
Author metadata contains 11 Louvain states and two conditions: Control and HeatShock.

Target states are selected using cell counts only:
- HeatShock >= 100 cells
- Control >= 40 cells

This deterministically yields held-out targets:
- cluster 2
- cluster 3
- cluster 10

Donor states are selected using cell counts only:
- not a target state
- Control >= 20 cells
- HeatShock >= 20 cells

This yields donors:
- cluster 5
- cluster 6
- cluster 8

HeatShock expression from target clusters 2, 3 and 10 MUST NOT be read during freeze.

## Prediction
For each donor state d:

response_d(g) = log1p(CPM_HeatShock,d(g)) - log1p(CPM_Control,d(g))

Shared acute-response prediction:

prediction(g) = mean(response_5(g), response_6(g), response_8(g))

The same shared response vector is used to rank candidates for all three held-out target states.

Direction = sign(prediction).
Candidate strength = abs(prediction).

## Technical gene universe
Defined before target reveal using only:
- all Control cells; and
- HeatShock cells from donor states 5, 6 and 8.

A gene is technical-evaluable if it has at least 20 total raw counts across these allowed cells.
Results should also be reported for the all-gene universe as a sensitivity view.

## Evaluation after target reveal
For each held-out target state 2, 3 and 10:
- full-vector Spearman correlation: secondary descriptive metric;
- overall sign agreement: secondary;
- top-k candidate overlap for k = 100, 250, 500, 1000;
- enrichment over random expectation;
- exact hypergeometric p-value;
- sign agreement among overlapping top-k genes;
- median percentile of true target top-k genes in the frozen predicted ranking.

Headline interpretation focuses on candidate-ranking utility rather than requiring accurate prediction of every gene.

## Baselines / references
- pooled-donor response: combine donor cells before computing one HeatShock-Control response;
- each individual donor response (5, 6, 8), reported separately as diagnostic references;
- random expectation for top-k overlap.

Do not choose the best donor after seeing target outcomes and call it deployable.

## No biological hard gate
No arbitrary rho/enrichment threshold is defined as project-level PASS/FAIL.
The result is interpreted continuously and claim-by-claim.

## Integrity
Before unblinding, save:
- rule-file SHA256;
- source-object SHA256;
- technical gene universe;
- shared prediction vector and ranking;
- explicit list of expression cells read during freeze;
- explicit confirmation that HeatShock target cells from clusters 2/3/10 were excluded.
