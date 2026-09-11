# Public-release checklist

Before making the repository public or minting a Zenodo DOI:

1. Create a private GitHub repository and upload this folder without modification.
2. Confirm the GitHub Actions smoke test passes.
3. Run the full downstream analysis locally with `python scripts/run_downstream_reproducibility.py`.
4. Do not add the excluded historical `ALPHA = 0.7` development script; the final analysis parameter is `alpha = 0.5`.
5. Do not upload CNKI/WanFang source files or other third-party files unless redistribution rights have been confirmed.
6. Replace placeholder repository/article identifiers in `CITATION.cff` and `SUBMISSION_CODE_AVAILABILITY.txt` after GitHub/Zenodo URLs and DOIs exist.
7. Confirm author list/order and software license with all co-authors before public release.
8. Create a GitHub release (suggested tag: `v1.0.0`) and archive that release with Zenodo.
9. Put the Zenodo DOI (preferably the concept DOI) into the manuscript Code Availability statement.

## Validation performed during repository curation

The self-contained downstream smoke run was executed successfully after locking the final RWR restart probability to `alpha = 0.5`. It reproduced the deterministic headline quantities checked by `scripts/verification/verify_downstream.py`, including the network size, target counts, module proximity values, top-200 enrichment, positive-D Spearman correlation, and PCI.

The full 1,000-permutation robustness suite is intentionally provided for local/archival execution; it is substantially more computationally intensive than the smoke test.
