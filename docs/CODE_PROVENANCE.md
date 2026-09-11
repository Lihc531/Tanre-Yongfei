# Code provenance and final-version decisions

The public repository was curated from the uploaded `code.zip`, the recovered final manuscript package, and the final supplementary workbook.

## Final RWR parameter

The author confirmed on 2026-09-10 that the final RWR restart probability was **0.5**. An older unsaved/development copy used `ALPHA = 0.7`. Because that file is inconsistent with the final numerical results and manuscript, it is excluded from the public repository rather than silently presented as an alternative final analysis.

The `alpha=0.5` implementation was independently checked against the supplied derived inputs. A single deterministic diffusion run reproduced the final observed proximity values to numerical precision (approximately 1e-16 absolute difference in diffusion scores relative to the recovered result table).

## Canonical transcriptomic script

The recovered standalone early GSE103119 script contained stale development references (including an unrelated `GSE56766` filename). The later integrated Figure 2 script contains the coherent final two-cohort processing logic for GSE103119 and GSE196399 and is therefore used as the canonical transcriptomic workflow in this repository. The stale standalone script is not distributed as final code.

## Final Figure 6

The `code.zip` contained an older cross-syndrome `figure6_V2.R`, whereas the submitted manuscript uses a final robustness/ablation Figure 6. The old cross-syndrome plotting script is excluded. This repository provides a clean `scripts/06_figures/figure6_from_robustness.py` that reconstructs the **submitted robustness Figure 6** directly from the final robustness metrics. It does not alter or recompute those metrics.

## Path refactoring

Historical absolute paths such as `E:\\PythonProject\\...` and `/mnt/data/...` were removed from active public scripts. Paths now resolve relative to the repository root or are supplied through command-line arguments. Statistical definitions and final parameter values were preserved.
