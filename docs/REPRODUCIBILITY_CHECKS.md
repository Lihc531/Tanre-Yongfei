# Reproducibility checks

The curated repository includes `scripts/verification/verify_downstream.py` and a top-level runner.

The following final quantities are checked automatically after a full run:

| Quantity | Expected value |
|---|---:|
| PPI edges (STRING 700) | 236,837 |
| PPI nodes | 16,194 |
| CAP genes in PPI | 10,245 |
| core targets total | 337 |
| core targets in PPI | 332 |
| core targets in CAP–PPI | 203 |
| CAP-up genes | 2,796 |
| CAP-down genes | 2,398 |
| observed CAP-up proximity | 5.27224330582954e-05 |
| observed CAP-down proximity | 7.37877154350432e-06 |
| top-200 CAP-up overlap | 73 |
| top-200 CAP-up OR | 1.54561463878711 |
| top-200 CAP-up Fisher P | 0.00254183083187808 |
| positive-D Spearman rho | 0.149238974678315 |
| PCI (positive alignment) | about 571.5 |

The verification uses tolerances for floating-point values rather than requiring byte-identical output.

Permutation P-values are stochastic but deterministic for a fixed implementation and seed. The final scripts use 1,000 permutations and seed 42 (CAP-down proximity uses seed 43, matching the recovered main script). Because the sparse implementation is mathematically equivalent to the recovered row-normalized implementation, deterministic diffusion scores were validated against the recovered result table before packaging.
