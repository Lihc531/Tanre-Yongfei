# Reproducibility repository: Tanre Yongfei–CAP systems pharmacology study

This repository contains the analysis code and redistributable derived inputs for the manuscript:

> **Preferential alignment of Tanre Yongfei-associated herb-target signals with activated inflammatory transcriptional modules in community-acquired pneumonia: a systems pharmacology analysis**

The repository was curated from the authors' recovered analysis scripts and the final supplementary workbook. It is intended for peer review and public archiving (e.g., GitHub + Zenodo). The final RWR parameter choice (`alpha = 0.5`) was confirmed by the author during repository curation.

## Important parameter clarification

The **final random-walk-with-restart (RWR) restart probability is `alpha = 0.5`**. The author confirmed that some older working copies were not saved after parameter tuning. An obsolete script containing `ALPHA = 0.7` is therefore **not included** in this public repository. The included main-network and robustness code uses `alpha = 0.5`, matching the final manuscript, final Figure 4/5 analysis, and the numerical values in the final robustness table.

## What is reproducible from this repository

The downstream network analysis is self-contained using the included derived inputs:

1. CAP meta-disease vector (`data/derived/CAP_D_meta.tsv`)
2. STRING-derived gene-gene edge lists at confidence thresholds 400/700/900
3. Final core and broader Tanre Yongfei target-gene sets
4. RWR diffusion, module-level localization, Mann–Whitney tests, Fisher enrichment, permutation proximity
5. Diffusion–disease correlations
6. Perturbation concordance index (PCI) and gene-level composite scores
7. Raw-vs-diffusion, STRING-threshold, and target-space robustness analyses
8. A clean reconstruction script for the final Figure 6 from the final robustness metrics

The GAT and transcriptomic preprocessing code is also included, but complete end-to-end rerunning of those upstream stages requires external/restricted source files described in `data/restricted_inputs/README.md` and `docs/DATA_AVAILABILITY.md`.

## One-command downstream reproduction

From the repository root:

```bash
python scripts/run_downstream_reproducibility.py
```

This runs the final `alpha=0.5` network analysis, correlations, PCI, all three robustness analyses, verifies the generated metrics against the final reference tables, and writes outputs to `results/generated/`.

For a faster check that skips the full 1,000-permutation robustness suite:

```bash
python scripts/run_downstream_reproducibility.py --smoke
```

## Expected headline checks

A successful full downstream run should recover, within floating-point tolerance:

- CAP genes in the 700-threshold PPI analysis space: **10,245**
- Core targets: **337 total; 332 in PPI; 203 in the CAP–PPI space**
- CAP-up genes: **2,796**; CAP-down genes: **2,398**
- CAP-up observed proximity: **5.2722433e-05**
- CAP-down observed proximity: **7.3787715e-06**
- Top-200 CAP-up enrichment: **73 genes; OR 1.5456; P 0.0025418**
- Diffusion vs positive CAP component: **Spearman rho 0.149239**
- Positive-alignment PCI: **~571.5**

The full 1,000-permutation run should also reproduce the empirical P-value pattern reported in the manuscript (CAP-up significant; CAP-down non-significant).

## Repository structure

```text
scripts/
  01_gat/                 recovered multi-label GAT implementation
  02_transcriptomics/     canonical two-cohort CAP preprocessing / Figure 2 workflow
  03_ppi/                 STRING protein-to-gene edge conversion
  04_network/             main RWR, correlations, PCI
  05_robustness/          raw-vs-diffusion / PPI-threshold / target-space sensitivity
  06_figures/             reproducible Figure 6 plotting from final metrics
  optional_acquisition/   source-acquisition helpers; not needed for downstream rerun
  verification/           numerical checks against final reference outputs

data/
  derived/                redistributable derived inputs used by downstream analyses
  restricted_inputs/      schema and placement instructions for upstream source files
results/
  reference/              final manuscript/supplementary values used for validation
  generated/              outputs created by rerunning the repository
```

## Environment

For the self-contained downstream analysis, install `environment/requirements-core.txt`. GAT-specific dependencies are separated in `environment/requirements-gat.txt`, and optional acquisition/enrichment dependencies are listed in `environment/requirements-optional.txt`. R dependencies for the transcriptomic workflow are listed in `environment/R_PACKAGES.md`.

## Reproducibility status and provenance

See:

- `docs/CODE_PROVENANCE.md`
- `docs/DATA_AVAILABILITY.md`
- `docs/REPRODUCIBILITY_CHECKS.md`
- `docs/EXCLUDED_LEGACY_FILES.md`

These documents explicitly distinguish final/validated code from recovered but obsolete development copies.

## Citation

Citation metadata are provided in `CITATION.cff`. Add the final article DOI and Zenodo DOI after publication/archive creation.

## License

Code in this repository is released under the MIT License. Third-party source data remain subject to their original licenses and terms of use; see `docs/DATA_AVAILABILITY.md`.

## Validation status

A clean downstream smoke run was executed successfully during curation and passed the automated deterministic checks. A GitHub Actions workflow (`.github/workflows/smoke-test.yml`) is included so that the same smoke test can run automatically after upload. The complete 1,000-permutation robustness suite is available through the default runner and is expected to require substantially more compute time than the smoke test.

Before public release, please read `docs/RELEASE_CHECKLIST.md`.
