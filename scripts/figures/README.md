# Figure rendering

The established workflow is Python/Matplotlib. Main outputs select Figure 1 v15, Figures 2/3 v11 and Figure 4 v16. Legacy modules 5/6/7/8/9/13/14 are source dependencies inherited by those final renderers. Do not remove them solely because their version number is older. Source paths and the alignment-helper import were made portable; numerical selections and drawing logic were retained.

The source inputs include sixteen rice contexts, 126,139 deposited rice nuclei coordinates, 2,085 heat-study display cells and 16,335 cross-study gene scores. Embeddings are display context, not evidence used to fit a new predictor. Two later ED sensitivity renderers are included. Earlier ED panels outside these two sensitivities are not claimed to be completely reconstructed by this export.

Use `python scripts/utilities/reproduce.py figures` or `extended-data`. Rendered PDFs/PNGs/SVGs are local outputs and are excluded from the initial repository. Figure QA is separate from statistical inference. The local alignment gate checks explicit comparable axes; heterogeneous inset roles retain declared exemptions. External audit tooling used during export is not a runtime dependency.
